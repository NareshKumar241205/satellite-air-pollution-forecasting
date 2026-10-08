"""Verification metrics on (origin, lead) pairs: errors, skill vs climatology, anomaly correlation,
interval coverage, CRPS from residual samples, and moving-block bootstrap CIs over origins.

Only observed targets count. Missing values (zeros, unusable windows, off-footprint pixels) never enter a metric.
"""

from __future__ import annotations

import numpy as np


def target_obs(values: np.ndarray, H: int) -> np.ndarray:
    """(T, H, G, Hp, Wp) observed value of window t+h (NaN beyond the record)."""
    T = values.shape[0]
    out = np.full((T, H) + values.shape[1:], np.nan)
    for h in range(1, H + 1):
        out[: T - h, h - 1] = values[h:]
    return out


def crps_from_sorted(err: np.ndarray, r_sorted: np.ndarray) -> np.ndarray:
    """Empirical CRPS of the predictive distribution pred + r (r = residual samples) at observations y.

    err = y - pred (any shape). CRPS = E|r - err| - 0.5 E|r - r'|, computed exactly from the sorted residuals.
    """
    n = len(r_sorted)
    csum = np.concatenate([[0.0], np.cumsum(r_sorted)])
    k = np.searchsorted(r_sorted, err)                     # residuals below err
    below = k * err - csum[k]
    above = (csum[n] - csum[k]) - (n - k) * err
    e_abs = (below + above) / n
    spread = 2.0 / n**2 * np.sum((2 * np.arange(1, n + 1) - n - 1) * r_sorted)
    return e_abs - 0.5 * spread


def per_origin_sse(pred: np.ndarray, obs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """pred, obs: (n_origin, ...) -> SSE and count per origin."""
    ok = np.isfinite(pred) & np.isfinite(obs)
    d = np.where(ok, pred - obs, 0.0)
    ax = tuple(range(1, pred.ndim))
    return (d * d).sum(ax), ok.sum(ax)


def block_bootstrap_skill(sse_m: np.ndarray, sse_ref: np.ndarray, block: int, n: int, seed: int = 0,
                          ci: float = 0.95) -> tuple[float, float]:
    """CI for skill = 1 - ΣSSE_m/ΣSSE_ref, resampling contiguous blocks of origins (moving block bootstrap)."""
    m = len(sse_m)
    if m == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    b = min(block, m)
    starts = np.arange(m - b + 1)
    n_blocks = int(np.ceil(m / b))
    draws = rng.choice(starts, size=(n, n_blocks))
    idx = (draws[:, :, None] + np.arange(b)[None, None, :]).reshape(n, -1)[:, :m]
    sk = 1 - sse_m[idx].sum(1) / sse_ref[idx].sum(1)
    lo, hi = np.quantile(sk, [(1 - ci) / 2, 1 - (1 - ci) / 2])
    return float(lo), float(hi)
