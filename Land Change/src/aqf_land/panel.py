"""C2 panel: EXPLORATORY land -> NO2 estimates, run although G1 failed at 0.5 %/pp (decision of 12 Oct 2026).

Every number here must be read against the G1 detection limit: the pooled long-difference MDE is 0.73 % NO2 per
pp of built-up gain (80 % power, 5 % two-sided). An estimate smaller than that is not evidence of no effect.
NO2 is the raw same-season (Jan-Apr) mean: there is no ERA5 deweathering yet. Tuning years 2019-2023 only; 2024
and 2025 are never read here.

Model (pixel = one 0.04 deg block, t = year):
    y_pt = a_p + g_{site,t} + b * L_pt + th * WL_pt + e_pt
y = 100 * log(NO2 seasonal mean), so b is in % NO2 per pp. L is the Dynamic World built-up fraction (pp), and WL
is the mean L of the 4 rook neighbours in the same site (a spatial lag of the land variable, SLX). Pixel fixed
effects remove everything static; site-year effects remove weather, the lockdown and processor steps common to a
site. The fixed effects are removed by alternating projections.

Inference: spatial-block bootstrap. Blocks of 4 x 4 pixels (~17 km, about the scale where the neighbour
correlation of the residual fades) are resampled with replacement within each site, keeping each block's whole
time series.
Placebo: the land trajectories of whole 4 x 4 blocks are permuted to other blocks of the same site, and b is
re-estimated. That gives the null distribution of b when land and NO2 are unrelated in space.
Lead/lag: L_{t+1} (a lead: future land must not predict NO2 now) and L_{t-1} enter together with L_t.
Long difference: dy = 100 dlog NO2 (2019 -> 2023) on dL with site fixed effects, the G1 design.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

BLOCK = 4
EXPLORATORY = "EXPLORATORY (G1 failed at 0.5 %/pp; MDE 0.73 %/pp); raw Jan-Apr NO2, no ERA5 deweathering"


def build_panel(npz: dict, sites: list[str], years: list[int], k_built: int,
                y_log: dict | None = None) -> pd.DataFrame:
    """y_log: optional {(site, year): (Hb, Wb)} already in 100*log units (e.g. deweathered) instead of the npz NO2."""
    rows = []
    for s in sites:
        for y in years:
            land = npz[f"{s}_{y}"][k_built] * 100
            no2 = npz[f"{s}_{y}_no2"]
            H, W = land.shape
            pad = np.pad(land, 1, constant_values=np.nan)
            nb = np.nanmean(np.stack([pad[:-2, 1:-1], pad[2:, 1:-1], pad[1:-1, :-2], pad[1:-1, 2:]]), 0)
            r, c = np.indices((H, W))
            with np.errstate(invalid="ignore", divide="ignore"):
                yv = y_log[(s, y)] if y_log is not None else np.where(no2 > 1e-6, 100 * np.log(no2), np.nan)
            rows.append(pd.DataFrame({"site": s, "year": y, "row": r.ravel(), "col": c.ravel(), "L": land.ravel(),
                                      "WL": nb.ravel(), "y": yv.ravel()}))
    p = pd.concat(rows, ignore_index=True)
    p["pixel"] = p.site + ":" + p.row.astype(str) + ":" + p.col.astype(str)
    p["blk"] = p.site + ":" + (p.row // BLOCK).astype(str) + ":" + (p.col // BLOCK).astype(str)
    p = p.sort_values(["pixel", "year"])
    g = p.groupby("pixel").L
    p["L_lead"], p["L_lag"] = g.shift(-1), g.shift(1)
    return p.reset_index(drop=True)


def demean(df: pd.DataFrame, cols: list[str], fe: tuple[str, str], tol: float = 1e-8, it: int = 200) -> np.ndarray:
    """Two-way within transformation by alternating projections (unbalanced panels)."""
    X = df[cols].to_numpy(float).copy()
    a, b = df[fe[0]].to_numpy(), df[fe[1]].to_numpy()
    ia, ib = pd.factorize(a)[0], pd.factorize(b)[0]
    na, nb = ia.max() + 1, ib.max() + 1
    ca, cb = np.bincount(ia, minlength=na), np.bincount(ib, minlength=nb)
    for _ in range(it):
        prev = X.copy()
        for idx, cnt, n in ((ia, ca, na), (ib, cb, nb)):
            means = np.stack([np.bincount(idx, X[:, j], n) for j in range(X.shape[1])], 1) / cnt[:, None]
            X -= means[idx]
        if np.abs(X - prev).max() < tol:
            break
    return X


def fit(df: pd.DataFrame, xcols: list[str]) -> np.ndarray:
    d = df.dropna(subset=["y", *xcols])
    d = d[d.groupby("pixel").y.transform("size") >= 2]
    Z = demean(d.assign(sy=d.site + ":" + d.year.astype(str)), ["y", *xcols], ("pixel", "sy"))
    beta, *_ = np.linalg.lstsq(Z[:, 1:], Z[:, 0], rcond=None)
    return beta


def block_bootstrap(df: pd.DataFrame, xcols: list[str], n: int, rng: np.random.Generator) -> np.ndarray:
    groups = {s: g for s, g in df.groupby("site")}
    blocks = {s: {b: gb for b, gb in g.groupby("blk")} for s, g in groups.items()}
    out = []
    for i in range(n):
        parts = []
        for s, bl in blocks.items():
            keys = list(bl)
            for j, k in enumerate(rng.choice(len(keys), len(keys), replace=True)):
                part = bl[keys[k]].copy()
                part["pixel"] = part.pixel + f"#{j}"        # a repeated block is a new set of pixels
                parts.append(part)
        out.append(fit(pd.concat(parts, ignore_index=True), xcols))
    return np.array(out)


def placebo(df: pd.DataFrame, n: int, rng: np.random.Generator) -> np.ndarray:
    """b with land trajectories of 4x4 blocks moved to other blocks of the same site (NO2 stays in place)."""
    out = []
    land_cols = ["L", "WL"]
    for i in range(n):
        parts = []
        for s, g in df.groupby("site"):
            bl = sorted(g.blk.unique())
            perm = dict(zip(bl, rng.permutation(bl)))
            src = g.set_index(["blk", "row", "col", "year"])[land_cols]
            moved = g.copy()
            # same position inside the block, same year, land from the permuted block (NaN if that block is smaller)
            key = pd.MultiIndex.from_arrays([moved.blk.map(perm), moved.row % BLOCK, moved.col % BLOCK, moved.year])
            src2 = src.copy()
            src2.index = pd.MultiIndex.from_arrays([src.index.get_level_values(0),
                                                    src.index.get_level_values(1) % BLOCK,
                                                    src.index.get_level_values(2) % BLOCK,
                                                    src.index.get_level_values(3)])
            src2 = src2[~src2.index.duplicated()]
            moved[land_cols] = src2.reindex(key).to_numpy()
            parts.append(moved)
        out.append(fit(pd.concat(parts, ignore_index=True), land_cols)[0])
    return np.array(out)


def long_difference_fit(df: pd.DataFrame, y0: int, y1: int) -> tuple[float, pd.DataFrame]:
    a = df[df.year == y0].set_index("pixel")
    b = df[df.year == y1].set_index("pixel")
    d = pd.DataFrame({"site": a.site, "blk": a.blk, "dy": b.y - a.y, "dL": b.L - a.L, "dWL": b.WL - a.WL}).dropna()
    for c in ("dy", "dL", "dWL"):
        d[c] = d[c] - d.groupby("site")[c].transform("mean")
    beta, *_ = np.linalg.lstsq(d[["dL", "dWL"]].to_numpy(), d.dy.to_numpy(), rcond=None)
    return float(beta[0]), d


def long_difference_bootstrap(d: pd.DataFrame, n: int, rng: np.random.Generator) -> np.ndarray:
    blocks = {s: [g for _, g in gs.groupby("blk")] for s, gs in d.groupby("site")}
    out = []
    for _ in range(n):
        parts = []
        for s, bl in blocks.items():
            parts += [bl[k] for k in rng.choice(len(bl), len(bl), replace=True)]
        x = pd.concat(parts)
        for c in ("dy", "dL", "dWL"):
            x[c] = x[c] - x.groupby("site")[c].transform("mean")
        out.append(np.linalg.lstsq(x[["dL", "dWL"]].to_numpy(), x.dy.to_numpy(), rcond=None)[0][0])
    return np.array(out)
