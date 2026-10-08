"""F2 feature matrix: one row per (origin t, lead h, footprint pixel), for one gas.

Everything is in standardised anomaly units: z = anomaly / sd, where the anomaly is relative to the
train-fitted climatology and sd is the train-year anomaly standard deviation of the gas. Missing values
(zeros, unusable windows) are filled with 0, which is the climatology, and flagged.

Origin features use windows t, t-1, ... only, never t+1 or later. `tests/test_forecast.py` checks this.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.ndimage import convolve

from aqf_forecast.problem import Fitted, Problem


def _trailing_nanmean(x: np.ndarray, n: int) -> np.ndarray:
    """Mean of the finite values among x[t-n+1 .. t] along axis 0 (NaN if none). x: (T, ...)."""
    fin = np.isfinite(x)
    cs = np.cumsum(np.where(fin, x, 0.0), 0)
    cc = np.cumsum(fin, 0).astype(float)
    s, c = cs.copy(), cc.copy()
    s[n:] -= cs[:-n]
    c[n:] -= cc[:-n]
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(c > 0, s / c, np.nan)


@dataclass
class FeatureSet:
    origin: np.ndarray        # (T, P, n_origin) per origin and footprint pixel
    target: np.ndarray        # (T_ext, P, n_target) per target window and pixel
    static: np.ndarray        # (P, n_static)
    names: list[str]          # origin + target + static + ["lead"]
    z: np.ndarray             # (T, P) standardised observed anomaly (the target), NaN = missing
    sd: float                 # anomaly sd (mol/m^2) used for standardising
    pix: tuple[np.ndarray, np.ndarray]   # footprint pixel indices (rows, cols)


def build(p: Problem, f: Fitted, g: int, s2: dict[str, np.ndarray] | None, fit_mask: np.ndarray,
          sd: float | None = None, other_gases: bool | None = None) -> FeatureSet:
    """sd: anomaly scale; by default the train-window anomaly sd of the gas (pass it to reuse a fitted scale).
    other_gases: add the other gases' anomalies (default from `f2.other_gases` in the config)."""
    fc = p.fcfg["f2"]
    fp = p.cube.footprint
    rows, cols = np.nonzero(fp)
    T = p.T
    a = f.anom[:, g]                                                   # (T, Hp, Wp)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        sd = float(np.nanstd(a[fit_mask][:, fp])) if sd is None else sd
        z = a / sd
        area = np.nanmean(z[:, fp], 1)                                 # (T,) area-mean anomaly
    names, feats = [], []

    def add(name, arr):                                                # arr: (T, Hp, Wp) or (T,) broadcast
        names.append(name)
        feats.append(arr[:, rows, cols] if arr.ndim == 3 else np.repeat(arr[:, None], len(rows), 1))

    for lag in range(fc["lags"]):
        zl = np.full_like(z, np.nan)
        zl[lag:] = z[: T - lag] if lag else z
        add(f"z_lag{lag}", np.nan_to_num(zl))
        if lag < 3:
            add(f"missing_lag{lag}", (~np.isfinite(zl)).astype(float))
    k = np.ones((1, 3, 3))
    k[0, 1, 1] = 0
    fin = np.isfinite(z)
    with np.errstate(invalid="ignore", divide="ignore"):
        nb = convolve(np.where(fin, z, 0), k, mode="constant") / convolve(fin.astype(float), k, mode="constant")
    for lag in range(fc["neighbour_lags"]):
        nl = np.full_like(nb, np.nan)
        nl[lag:] = nb[: T - lag] if lag else nb
        add(f"neigh_lag{lag}", np.nan_to_num(nl))
    add("area_lag0", np.nan_to_num(area))
    for n in fc["trailing"]:
        add(f"area_mean{n}", np.nan_to_num(_trailing_nanmean(area, n)))
        add(f"pixel_mean{n}", np.nan_to_num(_trailing_nanmean(z, n)))
    add("n_valid_last6", _trailing_nanmean(np.isfinite(z).astype(float), 6) * 6)
    if fc.get("other_gases") if other_gases is None else other_gases:
        for g2, gas2 in enumerate(p.gases):
            if g2 != g:
                z2 = f.anom[:, g2] / float(np.nanstd(f.anom[:, g2][fit_mask][:, fp]))
                add(f"{gas2}_lag0", np.nan_to_num(z2))
                add(f"{gas2}_mean6", np.nan_to_num(_trailing_nanmean(z2, 6)))
    origin = np.stack(feats, -1).astype(np.float32)                    # (T, P, n_origin)

    # Target-window features: standardised climatology level and season
    clim = f.clim_ext[:, g][:, rows, cols]                             # (T_ext, P)
    c_fit = clim[: T][fit_mask]
    clim_z = (clim - np.nanmean(c_fit)) / np.nanstd(c_fit)
    tdates = p.target_dates(f.clim_ext.shape[0] - T)
    w = 2 * np.pi * np.asarray(tdates.dayofyear, float) / 365.25
    target = np.stack([clim_z, np.repeat(np.sin(w)[:, None], len(rows), 1),
                       np.repeat(np.cos(w)[:, None], len(rows), 1)], -1).astype(np.float32)
    t_names = ["target_clim_z", "target_doy_sin", "target_doy_cos"]

    s_names = ["pixel_row", "pixel_col"]
    stat = [rows.astype(float), cols.astype(float)]
    if s2 is not None:
        for k_, v in s2.items():
            s_names.append(f"s2_{k_}")
            stat.append(v[rows, cols])                                 # NaN over the sea (LightGBM handles NaN)
    static = np.stack(stat, -1).astype(np.float32)
    return FeatureSet(origin, target, static, names + t_names + s_names + ["lead"], z[:, rows, cols], sd, (rows, cols))


def rows_for(fs: FeatureSet, origins: np.ndarray, leads: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Stack the feature rows for every (origin, lead) pair over all footprint pixels.

    Returns X (n_pairs * P, n_features), y (standardised target anomaly, NaN = missing or beyond the record),
    and pair_index (n_pairs * P,), the position of each row's pair in the input arrays.
    """
    P = fs.origin.shape[1]
    tgt = origins + leads
    X = np.concatenate([fs.origin[origins], fs.target[tgt],
                        np.broadcast_to(fs.static, (len(origins), P, fs.static.shape[1])),
                        np.broadcast_to(leads[:, None, None].astype(np.float32), (len(origins), P, 1))], -1)
    y = np.full((len(origins), P), np.nan, np.float32)
    inside = tgt < fs.z.shape[0]
    y[inside] = fs.z[tgt[inside]]
    pair = np.repeat(np.arange(len(origins)), P)
    return X.reshape(-1, X.shape[-1]), y.ravel(), pair


def pairs(p: Problem, leads: list[int], target_mask: np.ndarray, min_origin: int,
          exclude: tuple[str, str] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """All (origin, lead) pairs with origin >= min_origin whose target window is inside the record and in target_mask.

    exclude=(start, end) drops pairs whose origin or target window starts in that date range.
    """
    o, l = [], []
    for h in leads:
        t = np.arange(min_origin, p.T - h)
        keep = target_mask[t + h]
        if exclude is not None:
            lo, hi = pd.Timestamp(exclude[0]), pd.Timestamp(exclude[1])
            in_ex = (p.dates >= lo) & (p.dates <= hi)
            keep &= ~in_ex[t] & ~in_ex[t + h]
        o.append(t[keep])
        l.append(np.full(keep.sum(), h))
    return np.concatenate(o), np.concatenate(l)
