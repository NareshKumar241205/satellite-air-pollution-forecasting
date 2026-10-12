"""ERA5 reduced to one scalar per city, 5-day window and variable (the meteorology covariates of ADR-001).

The NetCDF reading lives in scripts/cds_era5.py (it needs xarray); this module is the pure reduction so it can be
tested without the downloads. Windows use the TNAQ grid: 5-day steps anchored at 2018-01-01, continuous across
years. Hours are the ones requested (around the 08 UTC S5P overpass), so a window mean is an overpass mean.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

ANCHOR = pd.Timestamp("2018-01-01")
SUM_VARS = {"tp"}                      # accumulated precipitation: summed, everything else averaged


def window_start(times: pd.DatetimeIndex) -> pd.DatetimeIndex:
    days = (times.normalize() - ANCHOR).days
    return ANCHOR + pd.to_timedelta(days // 5 * 5, unit="D")


def box_mask(lat: np.ndarray, lon: np.ndarray, box: tuple[float, float, float, float], pad: float) -> np.ndarray:
    """(lat, lon) bool: grid cells whose centre lies in [south, west, north, east] widened by `pad` degrees."""
    s, w, n, e = box
    la = (lat >= s - pad) & (lat <= n + pad)
    lo = (lon >= w - pad) & (lon <= e + pad)
    return la[:, None] & lo[None, :]


def window_scalars(times: pd.DatetimeIndex, fields: dict[str, np.ndarray], lat: np.ndarray, lon: np.ndarray,
                   box: tuple[float, float, float, float], pad: float) -> pd.DataFrame:
    """fields: {short_name: (T, lat, lon)}. Returns one row per window with the box-mean of each variable, plus
    the wind speed (mean of hourly speeds) when u10 and v10 are present."""
    m = box_mask(lat, lon, box, pad)
    if not m.any():
        raise ValueError(f"no ERA5 cell inside {box} + {pad} deg")
    area = {k: np.nanmean(v[:, m], 1) for k, v in fields.items()}
    if "u10" in area and "v10" in area:
        area["ws10"] = np.nanmean(np.hypot(fields["u10"][:, m], fields["v10"][:, m]), 1)
    df = pd.DataFrame(area, index=times)
    df["window"] = window_start(times)
    agg = {k: ("sum" if k in SUM_VARS else "mean") for k in area}
    out = df.groupby("window").agg(agg)
    out["n_hours"] = df.groupby("window").size()
    return out.reset_index()
