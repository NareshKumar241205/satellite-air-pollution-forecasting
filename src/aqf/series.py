"""Per-window quality flags and area-averaged (domain) time series for each gas."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from aqf.data import S5PCube


def clean_values(cube: S5PCube, cfg: dict) -> np.ndarray:
    """Cube values with exact zeros set to NaN for the gases in `quality.zeros_as_missing`."""
    v = cube.values.copy()
    for g, gas in enumerate(cube.gases):
        if gas in cfg["quality"]["zeros_as_missing"]:
            v[:, g][v[:, g] == 0] = np.nan
    return v


def analysis_values(cube: S5PCube, cfg: dict) -> np.ndarray:
    """Clean values with every unusable (window, gas) map blanked: the input for all pixel-level analyses."""
    wt = window_table(cube, cfg)
    usable = wt.pivot(index="date", columns="gas", values="usable")[cube.gases].values  # (T, G)
    return np.where(usable[:, :, None, None], clean_values(cube, cfg), np.nan)


def window_table(cube: S5PCube, cfg: dict) -> pd.DataFrame:
    """One row per (window, gas): coverage, zero share, usability flag and area statistics (mol/m^2).

    `coverage` is the share of area pixels holding a valid value after zeros are treated as missing.
    Area statistics use only those pixels, and only usable windows. Outlier windows are flagged with a
    robust z-score but never removed.
    """
    fp = cube.footprint
    n_fp = fp.sum()
    clean = clean_values(cube, cfg)
    rows = []
    for g, gas in enumerate(cube.gases):
        raw = cube.values[:, g][:, fp]                      # (T, n_fp)
        x = clean[:, g][:, fp]
        n_finite = np.isfinite(raw).sum(1)
        zeros = (raw == 0).sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            zero_frac = np.where(n_finite > 0, zeros / n_finite, np.nan)
        all_zero = (n_finite > 0) & (zeros == n_finite)
        n_valid = np.isfinite(x).sum(1)
        usable = (n_valid / n_fp >= cfg["quality"]["min_coverage"]) & ~all_zero
        xm = np.where(usable[:, None], x, np.nan)
        with warnings.catch_warnings():  # all-NaN rows (unusable windows) give NaN, which is intended
            warnings.simplefilter("ignore", RuntimeWarning)
            mean, median = np.nanmean(xm, 1), np.nanmedian(xm, 1)
            p90 = np.nanpercentile(xm, 90, axis=1)
        med = np.nanmedian(mean)
        mad = 1.4826 * np.nanmedian(np.abs(mean - med))
        rz = (mean - med) / mad if mad > 0 else np.zeros_like(mean)
        rows.append(pd.DataFrame({
            "date": cube.dates, "gas": gas, "coverage": n_valid / n_fp, "zero_frac": zero_frac,
            "all_zero": all_zero, "usable": usable, "mean": mean, "median": median, "p90": p90,
            "outlier": np.abs(rz) > cfg["quality"]["outlier_robust_z"],
        }))
    df = pd.concat(rows, ignore_index=True)
    df["year"], df["month"] = df.date.dt.year, df.date.dt.month
    return df


def monthly(wt: pd.DataFrame, col: str = "mean") -> pd.DataFrame:
    """Monthly mean of usable window values: columns gas, year, month, value, n_windows."""
    u = wt[wt.usable]
    out = u.groupby(["gas", "year", "month"])[col].agg(["mean", "size"]).reset_index()
    return out.rename(columns={"mean": "value", "size": "n_windows"})


def season_of(months: pd.Series | np.ndarray, cfg: dict) -> np.ndarray:
    lut = {m: s for s, ms in cfg["seasons"].items() for m in ms}
    return np.array([lut[int(m)] for m in np.asarray(months)])
