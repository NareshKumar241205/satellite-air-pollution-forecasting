"""Added G1 sites, selected by Dynamic World built-up growth only (Land Change/SITE_SELECTION.md).

candidates(): fixed 48 x 48 sub-pixel tiles of the export grid that don't overlap an existing site and are mostly
inside the export footprint (availability, not NO2 values). select(): growth score from the Dynamic World tables,
threshold 1.5 pp, written to configs/sites_added.yaml.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import rasterio
import yaml

from aqf.data import ROOT
from aqf_land.tnaq import Box, city_boxes, list_files

TILE = 48
THRESHOLD_PP = 1.5
ADDED = ROOT / "configs/sites_added.yaml"


def footprint(cfg: dict) -> tuple[np.ndarray, object]:
    """(H, W) bool: valid in at least one 2019 NO2 window (the Tamil Nadu clip of the export)."""
    files = list_files(cfg, "NO2")
    fp, t = None, None
    for p in files[files.year == 2019].path:
        with rasterio.open(p) as src:
            a = src.read(1)
            t = src.transform
        ok = np.isfinite(a) & (a != cfg["tnaq"]["nodata"])
        fp = ok if fp is None else fp | ok
    return fp, t


def _overlaps(a: Box, b: Box) -> bool:
    return not (a.row0 + a.rows <= b.row0 or b.row0 + b.rows <= a.row0 or
                a.col0 + a.cols <= b.col0 or b.col0 + b.cols <= a.col0)


def candidates(cfg: dict) -> dict[str, list[float]]:
    fp, t = footprint(cfg)
    H, W = fp.shape
    b = cfg["grid"]["block"]
    existing = city_boxes(cfg, t, (H, W))
    out = {}
    for r0 in range(0, H - TILE + 1, TILE):
        for c0 in range(0, W - TILE + 1, TILE):
            tile = Box("", r0, c0, TILE, TILE)
            if any(_overlaps(tile, e) for e in existing):
                continue
            blocks = fp[r0:r0 + TILE, c0:c0 + TILE].reshape(TILE // b, b, TILE // b, b).any((1, 3))
            if blocks.mean() < 0.5:
                continue
            lat = t.f + t.e * (r0 + TILE / 2)
            lon = t.c + t.a * (c0 + TILE / 2)
            out[f"T{r0 // TILE:02d}_{c0 // TILE:02d}"] = [round(float(lat), 4), round(float(lon), 4)]
    return out


def growth_scores(dw_dir, names: list[str]) -> pd.Series:
    a = pd.read_csv(dw_dir / "landfrac_2019.csv")
    z = pd.read_csv(dw_dir / "landfrac_2023.csv")
    m = a.merge(z, on=["city", "row", "col"], suffixes=("_19", "_23"))
    land = (m.valid_frac_19 >= 0.5) & (m.water_19 < 0.5) & (m.valid_frac_23 >= 0.5)
    m["dL"] = (m.builtup_23 - m.builtup_19) * 100
    return m[land & m.city.isin(names)].groupby("city").dL.median().reindex(names)


def select(cand: dict[str, list[float]], dw_dir) -> tuple[dict, pd.DataFrame]:
    s = growth_scores(dw_dir, list(cand))
    table = pd.DataFrame({"lat": [v[0] for v in cand.values()], "lon": [v[1] for v in cand.values()],
                          "growth_median_pp": s.values}, index=list(cand))
    table["selected"] = table.growth_median_pp >= THRESHOLD_PP
    chosen = {k: cand[k] for k in table.index[table.selected]}
    ADDED.write_text("# Written by aqf_land.sites (rule: Land Change/SITE_SELECTION.md). Do not edit by hand.\n"
                     + yaml.safe_dump(chosen, sort_keys=True), encoding="utf-8")
    return chosen, table
