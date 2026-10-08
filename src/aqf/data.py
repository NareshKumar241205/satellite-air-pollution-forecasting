"""Config, manifest and raster loading.

The S5P stack is tiny (438 x 3 x 13 x 13 float32, ~0.9 MB), so it is built once into
`data/processed/s5p_cube.npz` and loaded whole. Sentinel-2 (558 x 558 x 12 per window) is too big
to hold at once on this machine, so it is always streamed one file at a time.
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import yaml

ROOT = Path(__file__).resolve().parents[2]
log = logging.getLogger(__name__)


def load_config(path: Path | str = ROOT / "configs/data.yaml") -> dict:
    with open(path) as f:
        cfg = yaml.safe_load(f)
    for k, v in cfg["paths"].items():
        cfg["paths"][k] = ROOT / v
    return cfg


def load_manifest(cfg: dict) -> pd.DataFrame:
    m = pd.read_csv(cfg["paths"]["manifest"], parse_dates=["start_date", "end_date"])
    m = m.sort_values("start_date").reset_index(drop=True)
    m["year"] = m.start_date.dt.year
    m["month"] = m.start_date.dt.month
    m["has_s2"] = m.s2_status.eq("ok")
    split = {y: s for s in ("train", "val", "test") for y in cfg["split"][f"{s}_years"]}
    m["split"] = m.year.map(split)
    return m


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def read_raster(path: Path, bands: list[int] | None = None) -> np.ndarray:
    """Read bands (1-based) as float32 with every non-finite value (the -inf nodata) set to NaN."""
    with rasterio.open(path) as src:
        a = src.read(bands).astype(np.float32)
    a[~np.isfinite(a)] = np.nan
    return a


@dataclass
class S5PCube:
    values: np.ndarray       # (T, G, H, W) mol/m^2, NaN = nodata
    dates: pd.DatetimeIndex  # window start dates, length T
    gases: list[str]
    lat: np.ndarray          # (H,) pixel-centre latitude, north to south
    lon: np.ndarray          # (W,) pixel-centre longitude

    @property
    def footprint(self) -> np.ndarray:
        """(H, W) bool: pixels that are valid in at least one window (the area of interest)."""
        return np.isfinite(self.values).any(axis=(0, 1))


def build_s5p_cube(cfg: dict, force: bool = False) -> S5PCube:
    out = cfg["paths"]["processed"] / "s5p_cube.npz"
    if out.exists() and not force:
        return load_s5p_cube(cfg)
    m = load_manifest(cfg)
    s5p_dir = cfg["paths"]["s5p_dir"]
    with rasterio.open(s5p_dir / m.s5p_file.iloc[0]) as src:
        ref_transform, ref_shape, ref_crs = src.transform, src.shape, src.crs
    stack = []
    for f in m.s5p_file:
        with rasterio.open(s5p_dir / f) as src:
            # Every window must share one grid, or the per-pixel analysis is meaningless.
            assert src.transform == ref_transform and src.shape == ref_shape and src.crs == ref_crs, f
        stack.append(read_raster(s5p_dir / f))
    values = np.stack(stack)
    H, W = ref_shape
    lon = ref_transform.c + ref_transform.a * (np.arange(W) + 0.5)
    lat = ref_transform.f + ref_transform.e * (np.arange(H) + 0.5)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, values=values, dates=m.start_date.values.astype("datetime64[D]"),
                        gases=np.array(cfg["gases"]), lat=lat, lon=lon)
    meta = {"n_windows": len(m), "shape": list(values.shape), "crs": str(ref_crs),
            "transform": list(ref_transform)[:6], "manifest_sha256": sha256(cfg["paths"]["manifest"])}
    out.with_suffix(".json").write_text(json.dumps(meta, indent=2))
    log.info("built %s %s", out, values.shape)
    return load_s5p_cube(cfg)


def load_s5p_cube(cfg: dict) -> S5PCube:
    z = np.load(cfg["paths"]["processed"] / "s5p_cube.npz")
    return S5PCube(z["values"], pd.DatetimeIndex(z["dates"]), list(z["gases"]), z["lat"], z["lon"])


def nearest_pixel(cube: S5PCube, lat: float, lon: float) -> tuple[int, int]:
    return int(np.abs(cube.lat - lat).argmin()), int(np.abs(cube.lon - lon).argmin())
