"""F1 baselines. Each returns predicted values in mol/m^2 for every origin and lead: shape (T, H, G, Hp, Wp).

- climatology: the per-pixel seasonal expectation for the target date
- persistence: the last valid observed value (up to `max_lookback` windows old), else climatology
- damped: climatology + phi(lead + age) x the last valid anomaly. Anomaly persistence decays toward climatology
  at the rate measured on the fitting years.
"""

from __future__ import annotations

import numpy as np

from aqf_forecast.problem import Fitted, Problem, last_valid

MODELS = ["climatology", "persistence", "damped"]


def _target_clim(f: Fitted, T: int, H: int) -> np.ndarray:
    idx = np.arange(T)[:, None] + np.arange(1, H + 1)[None, :]       # (T, H) target index on the extended grid
    return f.clim_ext[idx]                                            # (T, H, G, Hp, Wp)


def predict(name: str, p: Problem, f: Fitted, H: int | None = None) -> np.ndarray:
    H = H or p.fcfg["horizon"]
    lb = p.fcfg["persistence"]["max_lookback"]
    clim_t = _target_clim(f, p.T, H)
    if name == "climatology":
        return clim_t
    if name == "persistence":
        v, _ = last_valid(p.values, lb)                               # (T, G, Hp, Wp)
        pred = np.broadcast_to(v[:, None], clim_t.shape)
        return np.where(np.isfinite(pred), pred, clim_t)
    if name == "damped":
        a, age = last_valid(f.anom, lb)
        G = a.shape[1]
        lead = np.arange(1, H + 1)[None, :, None, None, None]
        lag = np.where(age[:, None] >= 0, lead + age[:, None], 0)     # (T, H, G, Hp, Wp)
        g_idx = np.arange(G)[None, None, :, None, None]
        coef = f.phi[g_idx, lag]
        anom_hat = np.where(age[:, None] >= 0, coef * a[:, None], 0.0)
        return clim_t + anom_hat
    raise ValueError(name)
