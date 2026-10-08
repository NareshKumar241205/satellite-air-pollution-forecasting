"""Trend statistics: seasonal Mann-Kendall test, seasonal Sen slope, Benjamini-Hochberg FDR.

The seasonal Mann-Kendall test (Hirsch et al., 1982) compares each month only with the same month in
other years. That makes it robust to the strong seasonal cycle in the gas columns. It doesn't need
normally distributed data and it tolerates missing values. With 6 years of data its power is low,
so a non-significant result means "no detectable trend", not "no trend".
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


def _mk_s_var(x: np.ndarray) -> tuple[float, float]:
    """Mann-Kendall S and its variance (with the ties correction) for one series. NaNs are dropped."""
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 2:
        return 0.0, 0.0
    s = float(np.sign(x[None, :] - x[:, None])[np.triu_indices(n, 1)].sum())
    _, counts = np.unique(x, return_counts=True)
    ties = (counts * (counts - 1) * (2 * counts + 5)).sum()
    var = (n * (n - 1) * (2 * n + 5) - ties) / 18.0
    return s, float(var)


@dataclass
class TrendResult:
    slope: float       # seasonal Sen slope, data units per year
    z: float
    p: float           # two-sided p-value
    n_pairs: int       # number of within-season pairs behind the slope

    @property
    def direction(self) -> str:
        return "increasing" if self.slope > 0 else "decreasing" if self.slope < 0 else "flat"


def seasonal_trend(values: np.ndarray, years: np.ndarray, seasons: np.ndarray) -> TrendResult:
    """Seasonal Mann-Kendall test plus the seasonal Sen slope.

    values, years, seasons: 1-D arrays of equal length, e.g. monthly means with year and month labels.
    The slope is the median of all within-season pairwise slopes (value change per year).
    """
    values, years, seasons = map(np.asarray, (values, years, seasons))
    S = V = 0.0
    slopes = []
    for m in np.unique(seasons):
        sel = (seasons == m) & np.isfinite(values)
        order = np.argsort(years[sel])
        x, t = values[sel][order], years[sel][order].astype(float)
        s, v = _mk_s_var(x)
        S, V = S + s, V + v
        if len(x) >= 2:
            i, j = np.triu_indices(len(x), 1)
            dt = t[j] - t[i]
            ok = dt != 0
            slopes.append((x[j] - x[i])[ok] / dt[ok])
    slopes = np.concatenate(slopes) if slopes else np.array([])
    if V <= 0:
        return TrendResult(np.nan, np.nan, np.nan, len(slopes))
    z = (S - np.sign(S)) / np.sqrt(V)  # continuity correction
    p = 2 * stats.norm.sf(abs(z))
    return TrendResult(float(np.median(slopes)) if len(slopes) else np.nan, float(z), float(p), len(slopes))


def fdr_bh(p: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Benjamini-Hochberg: a boolean array (same shape as p) of discoveries at FDR alpha. NaNs are never discoveries."""
    p = np.asarray(p, float)
    flat = p.ravel()
    ok = np.isfinite(flat)
    out = np.zeros(flat.shape, bool)
    q = flat[ok]
    if q.size:
        order = np.argsort(q)
        thresh = alpha * np.arange(1, q.size + 1) / q.size
        passed = q[order] <= thresh
        k = np.nonzero(passed)[0].max() + 1 if passed.any() else 0
        sig = np.zeros(q.size, bool)
        sig[order[:k]] = True
        out[np.nonzero(ok)[0]] = sig
    return out.reshape(p.shape)
