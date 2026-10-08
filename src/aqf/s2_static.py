"""Static Sentinel-2 land-surface features aggregated onto the S5P grid.

The analysis showed that Sentinel-2 land cover is static over 2019–2024 and is missing in about half of the monsoon
windows, so it enters forecasting only as fixed per-pixel context. The features are fitted on **train years
only**, so they can be used as model inputs without leaking validation or test data.

Steps:
1. Stream every train-year S2 file and accumulate the per-100 m-pixel mean NDVI, NDBI and B12 over cloud-free observations.
2. Classify a 100 m pixel as water when its mean NDVI < 0. The sea and the large tanks/lakes are water.
3. For each S5P pixel, report the land fraction and the land-only mean of each index.
NDMI isn't used because it is exactly -NDBI in this dataset.
"""

from __future__ import annotations

import json
import logging

import numpy as np
import rasterio

from aqf.data import load_manifest, read_raster

log = logging.getLogger(__name__)
FEATURES = ["land_frac", "ndvi_land", "ndbi_land", "b12_land", "builtup_frac"]


def build(cfg: dict, s5p_lat: np.ndarray, s5p_lon: np.ndarray, force: bool = False) -> dict[str, np.ndarray]:
    out = cfg["paths"]["processed"] / "s2_static_s5pgrid.npz"
    if out.exists() and not force:
        z = np.load(out)
        return {k: z[k] for k in FEATURES}
    m = load_manifest(cfg)
    m = m[m.has_s2 & (m.split == "train")]
    names = cfg["s2_bands"]
    bands = [names.index(b) + 1 for b in ("NDVI", "NDBI", "B12")]
    s2_dir = cfg["paths"]["s2_dir"]
    with rasterio.open(s2_dir / m.s2_file.iloc[0]) as src:
        H, W, t = src.height, src.width, src.transform
    sums = np.zeros((3, H, W))
    n = np.zeros((H, W))
    for f in m.s2_file:
        a = read_raster(s2_dir / f, bands)
        ok = np.isfinite(a).all(0)
        sums += np.where(ok, a, 0)
        n += ok
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = sums / n                                  # (3, H, W); NaN where never observed
    observed = n > 0
    land = observed & (mean[0] >= 0)
    builtup = land & (mean[1] > 0)                      # NDBI > 0: impervious/built-up surface

    # Map each 100 m pixel centre to the S5P pixel containing it
    lon = t.c + t.a * (np.arange(W) + 0.5)
    lat = t.f + t.e * (np.arange(H) + 0.5)
    dlat, dlon = s5p_lat[1] - s5p_lat[0], s5p_lon[1] - s5p_lon[0]
    rows = np.floor((lat - (s5p_lat[0] - dlat / 2)) / dlat).astype(int)
    cols = np.floor((lon - (s5p_lon[0] - dlon / 2)) / dlon).astype(int)
    R, C = np.meshgrid(rows, cols, indexing="ij")
    Hs, Ws = len(s5p_lat), len(s5p_lon)
    inside = (R >= 0) & (R < Hs) & (C >= 0) & (C < Ws) & observed
    idx = (R * Ws + C)[inside]

    def agg(values: np.ndarray, weight: np.ndarray) -> np.ndarray:
        num = np.bincount(idx, (values * weight)[inside], Hs * Ws)
        den = np.bincount(idx, weight[inside].astype(float), Hs * Ws)
        with np.errstate(invalid="ignore", divide="ignore"):
            return (num / den).reshape(Hs, Ws)

    ones = np.ones((H, W))
    feats = {
        "land_frac": agg(land.astype(float), ones),
        "ndvi_land": agg(np.nan_to_num(mean[0]), land),
        "ndbi_land": agg(np.nan_to_num(mean[1]), land),
        "b12_land": agg(np.nan_to_num(mean[2]), land),
        "builtup_frac": agg(builtup.astype(float), land),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, **feats)
    out.with_suffix(".json").write_text(json.dumps({
        "years": sorted(m.year.unique().tolist()), "n_images": len(m), "water_rule": "mean NDVI < 0",
        "builtup_rule": "land and mean NDBI > 0", "features": FEATURES}, indent=2))
    log.info("built %s from %d train images", out, len(m))
    return feats
