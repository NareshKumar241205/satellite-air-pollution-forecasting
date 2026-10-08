"""The forecasting problem: cleaned values, splits, fitted climatology and anomalies, (origin, lead) pairs.

Index convention: window t is the forecast origin, the data available is t and earlier, and the target is
window t + h for lead h = 1..H. A pair belongs to the split of its *target* window.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from aqf.anomaly import Climatology, fit_climatology
from aqf.data import S5PCube, build_s5p_cube, load_config
from aqf.series import analysis_values

STEP_DAYS = 5


@dataclass
class Problem:
    values: np.ndarray            # (T, G, H, W) mol/m^2, NaN = missing (zeros, unusable windows, off-footprint)
    dates: pd.DatetimeIndex       # (T,) window start dates
    split: np.ndarray             # (T,) 'train' | 'val' | 'test'
    gases: list[str]
    cube: S5PCube
    cfg: dict
    fcfg: dict

    @property
    def T(self) -> int:
        return len(self.dates)

    def target_dates(self, n_extra: int = 0) -> pd.DatetimeIndex:
        """The window grid extended n_extra steps past the last window (for future targets)."""
        return pd.date_range(self.dates[0], periods=self.T + n_extra, freq=f"{STEP_DAYS}D")

    def windows(self, *splits: str) -> np.ndarray:
        return np.isin(self.split, splits)


def load_problem() -> Problem:
    cfg = load_config()
    import yaml

    from aqf.data import ROOT

    fcfg = yaml.safe_load(open(ROOT / "configs/forecast.yaml"))
    fcfg["outputs"] = ROOT / fcfg["outputs"]
    cube = build_s5p_cube(cfg)
    year_split = {y: s for s in ("train", "val", "test") for y in cfg["split"][f"{s}_years"]}
    split = np.array([year_split[d.year] for d in cube.dates])
    return Problem(analysis_values(cube, cfg), cube.dates, split, cube.gases, cube, cfg, fcfg)


@dataclass
class Fitted:
    """Everything fitted on a subset of windows: climatology, anomalies and damped-persistence coefficients."""
    clim: Climatology
    clim_ext: np.ndarray          # (T + H_extra, G, H, W) climatology on the extended window grid
    anom: np.ndarray              # (T, G, H, W) observed anomaly, NaN where missing
    phi: np.ndarray               # (G, max_lag + 1) anomaly regression coefficient at each lag (phi[:, 0] = 1)
    fit_mask: np.ndarray          # (T,) windows used for fitting


def fit(p: Problem, fit_mask: np.ndarray, n_extra: int | None = None) -> Fitted:
    fc = p.fcfg
    H = fc["horizon"]
    lb = fc["persistence"]["max_lookback"]
    n_extra = H if n_extra is None else n_extra
    vals = np.where(fit_mask[:, None, None, None], p.values, np.nan)
    clim = fit_climatology(vals, p.dates, fc["climatology"]["harmonics"], fc["climatology"]["ridge"],
                           fc["climatology"]["min_obs"])
    clim_ext = clim(p.target_dates(n_extra))
    anom = p.values - clim_ext[: p.T]
    # Damped persistence: least-squares slope through the origin of a(t+L) on a(t), pooled over pixels and
    # over pairs whose origin and target are both fitting windows.
    max_lag = H + lb
    phi = np.ones((len(p.gases), max_lag + 1))
    a = np.where(fit_mask[:, None, None, None], anom, np.nan)
    for L in range(1, max_lag + 1):
        x, y = a[:-L], a[L:]
        ok = np.isfinite(x) & np.isfinite(y)
        num = np.where(ok, x * y, 0).sum(axis=(0, 2, 3))
        den = np.where(ok, x * x, 0).sum(axis=(0, 2, 3))
        phi[:, L] = num / den
    return Fitted(clim, clim_ext, anom, phi, fit_mask)


def last_valid(x: np.ndarray, lookback: int) -> tuple[np.ndarray, np.ndarray]:
    """For every origin t, the most recent finite value of x among t, t-1, ..., t-lookback and its age.

    x: (T, ...). Returns (value, age), both (T, ...); value NaN and age -1 where nothing is valid.
    """
    val = np.full(x.shape, np.nan)
    age = np.full(x.shape, -1, int)
    for k in range(lookback, -1, -1):              # smaller ages overwrite larger ones
        shifted = np.full(x.shape, np.nan)
        shifted[k:] = x[: x.shape[0] - k] if k else x
        ok = np.isfinite(shifted)
        val[ok], age[ok] = shifted[ok], k
    return val, age
