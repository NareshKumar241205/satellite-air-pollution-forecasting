"""Per-pixel harmonic climatology and anomalies, the common target space for every forecasting model.

The climatology is fitted only on the windows passed in (normally the train years) by least squares on
`1, cos(2πk·doy/365.25), sin(2πk·doy/365.25)` for k = 1..K, separately for each gas and pixel, using only
valid values (zeros already treated as missing).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


def harmonics(dates: pd.DatetimeIndex, k: int) -> np.ndarray:
    """(T, 1 + 2k) design matrix of annual harmonics on day of year."""
    doy = np.asarray(dates.dayofyear, float)
    cols = [np.ones_like(doy)]
    for j in range(1, k + 1):
        w = 2 * np.pi * j * doy / 365.25
        cols += [np.cos(w), np.sin(w)]
    return np.stack(cols, 1)


@dataclass
class Climatology:
    coef: np.ndarray        # (1 + 2k, G, H, W)
    k: int

    def __call__(self, dates: pd.DatetimeIndex) -> np.ndarray:
        """Expected value for each date: (T, G, H, W). NaN where the pixel had too few observations."""
        return np.einsum("tp,pghw->tghw", harmonics(dates, self.k), self.coef)


def fit_climatology(values: np.ndarray, dates: pd.DatetimeIndex, k: int = 3, ridge: float = 1e-3,
                    min_obs: int = 20) -> Climatology:
    """values: (T, G, H, W) with NaN = missing. `ridge` is relative to the mean diagonal of XᵀX."""
    X = harmonics(dates, k)
    P = X.shape[1]
    T, G, H, W = values.shape
    coef = np.full((P, G, H, W), np.nan)
    for g in range(G):
        for i in range(H):
            for j in range(W):
                y = values[:, g, i, j]
                ok = np.isfinite(y)
                if ok.sum() < min_obs:
                    continue
                Xo = X[ok]
                A = Xo.T @ Xo
                A += ridge * np.trace(A) / P * np.eye(P)
                coef[:, g, i, j] = np.linalg.solve(A, Xo.T @ y[ok])
    return Climatology(coef, k)
