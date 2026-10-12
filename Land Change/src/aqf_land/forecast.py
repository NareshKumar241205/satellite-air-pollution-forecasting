"""NO2 forecasting on the TNAQ sites (Land Change/DECISION_RULES.md, sections C-E).

Unit: one 0.04 deg block x one 5-day window, all 18 sites pooled. Values in umol/m^2. Anomaly = value minus the
per-pixel harmonic climatology (K = 3) fitted on the training years only (train + validation for the final refit).
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from aqf.anomaly import fit_climatology

SEASONS = {"Winter": [1, 2], "Pre-monsoon": [3, 4, 5], "SW monsoon": [6, 7, 8, 9], "NE monsoon": [10, 11, 12]}
MET = ["u10", "v10", "ws10", "blh", "t2m", "tp"]
LGB_PARAMS = dict(objective="regression", learning_rate=0.05, num_leaves=31, min_data_in_leaf=200,
                  feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1, lambda_l2=1.0, verbose=-1, seed=0)


class Data:
    """Stacked pixels: Y (T, P) umol/m^2 with NaN = missing; meta per pixel; dates per window."""

    def __init__(self, cfg: dict, no2: dict, land: dict[str, dict[int, np.ndarray]], era5: pd.DataFrame,
                 start: str = "2019-01-01"):
        from aqf_land.tnaq import block_mean

        g = cfg["grid"]
        scale = cfg["tnaq"]["gas_units"][cfg.get("gas", "NO2")][0]
        keep = no2["dates"] >= pd.Timestamp(start)
        self.dates = no2["dates"][keep]
        cols, meta = [], []
        for box in no2["boxes"]:
            b = block_mean(no2[box.city][keep][:, 0], g["block"], g["min_valid_frac"]) * scale   # (T, Hb, Wb)
            T, H, W = b.shape
            cols.append(b.reshape(T, H * W))
            r, c = np.indices((H, W))
            meta.append(pd.DataFrame({"site": box.city, "row": r.ravel(), "col": c.ravel()}))
        self.Y = np.concatenate(cols, 1)
        self.meta = pd.concat(meta, ignore_index=True)
        self.site_idx = pd.factorize(self.meta.site)[0]
        self.sites = list(dict.fromkeys(self.meta.site))
        self.nb = self._neighbours()
        self.land = land
        self.era5 = era5

    def _neighbours(self) -> list[np.ndarray]:
        key = {(s, r, c): i for i, (s, r, c) in enumerate(self.meta[["site", "row", "col"]].itertuples(index=False))}
        out = []
        for s, r, c in self.meta[["site", "row", "col"]].itertuples(index=False):
            out.append(np.array([key[(s, r + dr, c + dc)] for dr in (-1, 0, 1) for dc in (-1, 0, 1)
                                 if (dr or dc) and (s, r + dr, c + dc) in key]))
        return out

    def years(self) -> np.ndarray:
        return np.asarray(self.dates.year)


def climatology(d: Data, fit_mask: np.ndarray):
    """Per-pixel harmonic climatology fitted on the windows in fit_mask. Returns (T, P) expected values."""
    vals = d.Y[fit_mask][:, None, None, :]
    clim = fit_climatology(vals, d.dates[fit_mask], k=3, ridge=1e-3, min_obs=20)
    return clim(d.dates)[:, 0, 0, :]


def last_valid(a: np.ndarray, lookback: int = 3) -> np.ndarray:
    """(T, P) last valid anomaly within `lookback` windows back (inclusive of t), else 0."""
    out = np.zeros_like(a)
    have = np.zeros(a.shape, bool)
    for k in range(lookback + 1):
        sh = np.full_like(a, np.nan)
        sh[k:] = a[:len(a) - k] if k else a
        take = ~have & np.isfinite(sh)
        out[take] = sh[take]
        have |= take
    return out


def damping(a: np.ndarray, alast: np.ndarray, fit_rows: np.ndarray, hmax: int = 30) -> np.ndarray:
    """phi_h for h = 1..hmax: least squares through the origin of a(t+h) on a_last(t), origins in fit_rows."""
    phi = np.zeros(hmax + 1)
    T = len(a)
    for h in range(1, hmax + 1):
        o = np.flatnonzero(fit_rows[: T - h])
        x, y = alast[o], a[o + h]
        ok = np.isfinite(y) & (x != 0)
        phi[h] = float((x[ok] * y[ok]).sum() / (x[ok] ** 2).sum()) if ok.any() else 0.0
    return phi


def features(d: Data, a: np.ndarray, land_year_of, use_era5: bool = True, use_land: bool = True,
             shuffle_land: np.random.Generator | None = None) -> tuple[np.ndarray, list[str]]:
    """(T, P, F) features known at each origin t (for the target t+1)."""
    T, P = a.shape
    f, names = [], []
    a0 = np.nan_to_num(a)
    for k in range(6):
        sh = np.full_like(a, np.nan)
        sh[k:] = a[:T - k] if k else a
        f += [np.nan_to_num(sh), np.isnan(sh).astype(float)]
        names += [f"lag{k}", f"miss{k}"]
    warnings.simplefilter("ignore", RuntimeWarning)              # all-missing neighbourhoods are expected
    nbm = np.stack([np.nanmean(a[:, n], 1) if len(n) else np.full(T, np.nan) for n in d.nb], 1)
    sm = np.stack([np.nanmean(a[:, d.site_idx == s], 1) for s in range(len(d.sites))], 1)[:, d.site_idx]
    for name, x in (("nb", nbm), ("site", sm)):
        for k in range(3):
            sh = np.full_like(x, np.nan)
            sh[k:] = x[:T - k] if k else x
            f.append(np.nan_to_num(sh))
            names.append(f"{name}_lag{k}")
    for w in (6, 12, 18):
        cs = pd.DataFrame(a).rolling(w, min_periods=1).mean().to_numpy()
        f.append(np.nan_to_num(cs))
        names.append(f"trail{w}")
    doy = np.asarray(d.dates.dayofyear, float)
    for k in (1, 2):
        for fn in (np.cos, np.sin):
            f.append(np.repeat(fn(2 * np.pi * k * doy / 365.25)[:, None], P, 1))
            names.append(f"{fn.__name__}{k}")
    if use_era5:
        e = d.era5
        for v in MET:
            piv = e.pivot(index="window", columns="city", values=v).reindex(index=d.dates, columns=d.sites)
            x = piv.to_numpy()[:, d.site_idx]
            for k in (0, 1):
                sh = np.full_like(x, np.nan)
                sh[k:] = x[:T - k] if k else x
                f.append(np.nan_to_num(sh))
                names.append(f"{v}_t-{k}")
    if use_land:
        L = np.full((T, P), np.nan)
        for t in range(T):
            y = land_year_of(d.dates[t])
            col = np.concatenate([d.land[s][y].ravel() if y in d.land[s] else
                                  np.full((d.meta.site == s).sum(), np.nan) for s in d.sites])
            L[t] = col
        if shuffle_land is not None:
            for s in range(len(d.sites)):
                idx = np.flatnonzero(d.site_idx == s)
                L[:, idx] = L[:, shuffle_land.permutation(idx)]
        f.append(np.nan_to_num(L, nan=np.nanmean(L)))
        names.append("builtup")
    return np.stack(f, 2).astype(np.float32), names


def land_year(date: pd.Timestamp) -> int:
    """Most recent complete Jan-Apr season strictly before the origin window."""
    return date.year - 1 if date.month <= 4 else date.year


def skill(y: np.ndarray, pred: np.ndarray, clim: np.ndarray) -> float:
    ok = np.isfinite(y) & np.isfinite(pred) & np.isfinite(clim)
    return float(1 - ((y[ok] - pred[ok]) ** 2).mean() / ((y[ok] - clim[ok]) ** 2).mean())


def scores(y: np.ndarray, pred: np.ndarray, clim: np.ndarray, dates: pd.DatetimeIndex, site_idx: np.ndarray,
           rng: np.random.Generator, n_boot: int = 500, block: int = 6) -> dict:
    """y, pred, clim: (T, P) on the scored origins' targets. Block bootstrap over target windows."""
    ok = np.isfinite(y) & np.isfinite(pred) & np.isfinite(clim)
    e_m, e_c = np.where(ok, (y - pred) ** 2, 0), np.where(ok, (y - clim) ** 2, 0)
    n = ok.sum(1)
    sm, sc = e_m.sum(1), e_c.sum(1)
    T = len(y)
    nb = int(np.ceil(T / block))
    boot = []
    for _ in range(n_boot):
        starts = rng.integers(0, max(T - block, 1), nb)
        idx = np.concatenate([np.arange(s, min(s + block, T)) for s in starts])
        boot.append(1 - sm[idx].sum() / sc[idx].sum())
    res = {"skill": float(1 - sm.sum() / sc.sum()), "skill_lo": float(np.percentile(boot, 2.5)),
           "skill_hi": float(np.percentile(boot, 97.5)),
           "rmse": float(np.sqrt(sm.sum() / n.sum())), "mae": float(np.abs(y - pred)[ok].mean()),
           "rmse_clim": float(np.sqrt(sc.sum() / n.sum())), "n": int(n.sum())}
    pcs = []
    for t in range(T):
        for s in np.unique(site_idx):
            m = ok[t] & (site_idx == s)
            if m.sum() >= 10:
                a, b = (y - clim)[t, m], (pred - clim)[t, m]
                if a.std() > 0 and b.std() > 0:
                    pcs.append(np.corrcoef(a, b)[0, 1])
    res["pattern_corr"] = float(np.mean(pcs)) if pcs else np.nan
    for name, months in SEASONS.items():
        mm = dates.month.isin(months)
        res[f"skill_{name}"] = float(1 - sm[mm].sum() / sc[mm].sum()) if sc[mm].sum() > 0 else np.nan
    return res
