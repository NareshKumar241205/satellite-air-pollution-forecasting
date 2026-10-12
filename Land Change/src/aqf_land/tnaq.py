"""Tamil Nadu (TNAQ) export: file listing, city boxes on the fixed 0.01 deg grid, and 4 x 4 block aggregation.

Every TNAQ file shares one grid (asserted), so a city box is a fixed rasterio window. The window is read once
per file for all cities and cached in `data/processed/tnaq_<product>_cities.npz` at sub-pixel resolution;
blocking to the 0.04 deg analysis grid happens afterwards, so the block rule can change without re-reading.
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import yaml
from rasterio.windows import Window

from aqf.data import ROOT

log = logging.getLogger(__name__)


def load_land_config(path: Path | str = ROOT / "configs/land.yaml") -> dict:
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cfg["tnaq"]["dir"] = (ROOT / cfg["tnaq"]["dir"]).resolve()
    cfg["tnaq"]["extra_dirs"] = [(ROOT / d).resolve() for d in cfg["tnaq"].get("extra_dirs", [])]
    cfg["outputs"] = ROOT / cfg["outputs"]
    cfg["processed"] = ROOT / "data/processed"
    return cfg


def list_files(cfg: dict, product: str) -> pd.DataFrame:
    """One row per window file: date (window start), year, month, path. Sorted by date. The main export wins
    over an extra directory for the same window."""
    sub = cfg["tnaq"]["products"][product]
    found = {}
    for d in [cfg["tnaq"]["dir"], *cfg["tnaq"].get("extra_dirs", [])]:
        for p in sorted((d / sub).glob("*.tif")):
            m = re.search(r"_(\d{8})\.tif$", p.name)
            if m:
                found.setdefault(pd.Timestamp(m.group(1)), p)
    df = pd.DataFrame(sorted(found.items()), columns=["date", "path"])
    df["year"], df["month"] = df.date.dt.year, df.date.dt.month
    return df


@dataclass
class Box:
    city: str
    row0: int
    col0: int
    rows: int      # sub-pixels, a multiple of the block size
    cols: int

    @property
    def window(self) -> Window:
        return Window(self.col0, self.row0, self.cols, self.rows)


def city_boxes(cfg: dict, transform, shape: tuple[int, int]) -> list[Box]:
    """Box of `box_pixels` blocks centred on each city, snapped to the export grid and clipped to whole blocks."""
    b, n = cfg["grid"]["block"], cfg["grid"]["box_pixels"]
    half = n * b / 2
    boxes = []
    for city, (lat, lon) in cfg["cities"].items():
        col0 = int(round((lon - transform.c) / transform.a - half))
        row0 = int(round((lat - transform.f) / transform.e - half))
        r0, c0 = max(row0, 0), max(col0, 0)
        r1, c1 = min(row0 + n * b, shape[0]), min(col0 + n * b, shape[1])
        rows, cols = (r1 - r0) // b * b, (c1 - c0) // b * b
        if rows < b or cols < b:
            log.warning("%s is outside the export grid, skipped", city)
            continue
        boxes.append(Box(city, r0, c0, rows, cols))
    return boxes


def block_lat_lon(box: Box, transform, block: int) -> tuple[np.ndarray, np.ndarray]:
    """Block-centre latitude (north to south) and longitude."""
    lat = transform.f + transform.e * (box.row0 + block * (np.arange(box.rows // block) + 0.5))
    lon = transform.c + transform.a * (box.col0 + block * (np.arange(box.cols // block) + 0.5))
    return lat, lon


def box_bounds(box: Box, transform) -> tuple[float, float, float, float]:
    """(south, west, north, east) in degrees."""
    north = transform.f + transform.e * box.row0
    west = transform.c + transform.a * box.col0
    return north + transform.e * box.rows, west, north, west + transform.a * box.cols


def grid_boxes(cfg: dict) -> tuple[list[Box], object]:
    """City boxes and the export transform, from the first NO2 file (no full read)."""
    with rasterio.open(list_files(cfg, "NO2").path.iloc[0]) as src:
        return city_boxes(cfg, src.transform, src.shape), src.transform


def read_cities(cfg: dict, product: str, force: bool = False) -> dict:
    """{'dates', 'transform', 'boxes', city: (T, bands, rows, cols) float32 with NaN = nodata}. Cached."""
    out = cfg["processed"] / f"tnaq_{product.lower()}_cities.npz"
    if out.exists() and not force:
        return _load(out)
    files = list_files(cfg, product)
    if product == "S2":
        files = files[files.year >= 2019]              # 2018 is TOA, not interchangeable with SR
    with rasterio.open(files.path.iloc[0]) as src:
        ref = (src.transform, src.shape, src.crs)
        bands = [1] if product != "S2" else list(range(1, src.count + 1))
    boxes = city_boxes(cfg, ref[0], ref[1])
    stacks = {b.city: [] for b in boxes}
    nodata = cfg["tnaq"]["nodata"]
    for p in files.path:
        with rasterio.open(p) as src:
            assert (src.transform, src.shape, src.crs) == ref, p   # one grid, or per-pixel work is meaningless
            for b in boxes:
                a = src.read(bands, window=b.window).astype(np.float32)
                a[(a == nodata) | ~np.isfinite(a)] = np.nan
                stacks[b.city].append(a)
    out.parent.mkdir(parents=True, exist_ok=True)
    arrays = {f"city_{c}": np.stack(v) for c, v in stacks.items()}
    np.savez_compressed(out, dates=files.date.values.astype("datetime64[D]"),
                        transform=np.array(list(ref[0])[:6]), **arrays)
    out.with_suffix(".json").write_text(json.dumps({
        "product": product, "n_files": len(files), "first": str(files.date.min().date()),
        "last": str(files.date.max().date()), "source": str(files.path.iloc[0].parent),
        "boxes": [b.__dict__ for b in boxes]}, indent=2))
    log.info("cached %s: %d files, %d cities", out.name, len(files), len(boxes))
    return _load(out)


def _load(path: Path) -> dict:
    z = np.load(path)
    t = rasterio.Affine(*z["transform"])
    meta = json.loads(path.with_suffix(".json").read_text())
    d = {"dates": pd.DatetimeIndex(z["dates"]), "transform": t, "boxes": [Box(**b) for b in meta["boxes"]]}
    d.update({k[5:]: z[k] for k in z.files if k.startswith("city_")})
    return d


def block_mean(a: np.ndarray, block: int, min_valid_frac: float) -> np.ndarray:
    """Mean over block x block sub-pixels of the last two axes, NaN where fewer than min_valid_frac are valid."""
    *lead, H, W = a.shape
    r = a.reshape(*lead, H // block, block, W // block, block)
    ok = np.isfinite(r)
    n = ok.sum((-3, -1))
    with np.errstate(invalid="ignore", divide="ignore"):
        m = np.where(ok, r, 0).sum((-3, -1)) / n
    return np.where(n >= min_valid_frac * block * block, m, np.nan)
