"""Deweathered same-season NO2 per pixel and year, from ERA5 scalars per site and window (tuning years only).

Per pixel: y_w = 100 * log NO2 of each Jan-Apr window w (block mean). Regress it on standardised ERA5 site scalars
(10 m wind u, v and speed, boundary-layer height, 2 m temperature, precipitation) plus year dummies:
    y_w = c_year + X_w b + e_w
The year dummies mean b is identified only from weather variation within a season, so slow land-driven change
can't leak into the weather coefficients. The deweathered seasonal value is mean_w(y_w - (X_w - Xbar) b) per
year, where Xbar is the pixel's mean over the tuning windows.

Note: ERA5 scalars are site-wide, so in a panel with site x year effects the deweathering only changes results through
pixel-specific weather responses (e.g. plumes moving with wind direction), not through the site mean.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

MET = ["u10", "v10", "ws10", "blh", "t2m", "tp"]


def deweather_site(blocks: np.ndarray, dates: pd.DatetimeIndex, met: pd.DataFrame, years: list[int],
                   months: list[int], min_windows: int) -> dict[int, np.ndarray]:
    """blocks: (T, Hb, Wb) NO2 block means; met: rows indexed by window start. Returns {year: (Hb, Wb)} of
    deweathered seasonal 100*log NO2 (NaN where too few windows)."""
    sel = dates.year.isin(years) & dates.month.isin(months)
    d = dates[sel]
    X = met.reindex(d)[MET].to_numpy(float)
    ok_met = np.isfinite(X).all(1)
    Xs = (X - np.nanmean(X[ok_met], 0)) / np.nanstd(X[ok_met], 0)
    yrs = np.array(d.year)
    Y = np.where(blocks[sel] > 1e-6, 100 * np.log(np.where(blocks[sel] > 1e-6, blocks[sel], 1)), np.nan)
    T, H, W = Y.shape
    out = {y: np.full((H, W), np.nan) for y in years}
    D = np.stack([(yrs == y).astype(float) for y in years], 1)
    for i in range(H):
        for j in range(W):
            y = Y[:, i, j]
            ok = np.isfinite(y) & ok_met
            if ok.sum() < len(MET) + len(years) + 5:
                continue
            A = np.column_stack([D[ok], Xs[ok]])
            coef, *_ = np.linalg.lstsq(A, y[ok], rcond=None)
            b = coef[len(years):]
            adj = y - np.where(ok_met, np.nan_to_num(Xs - Xs[ok].mean(0)) @ b, np.nan)
            for k, yr in enumerate(years):
                m = ok & (yrs == yr)
                if m.sum() >= min_windows:
                    out[yr][i, j] = adj[m].mean()
    return out
