"""Extensions X2, X3, X5 (Land Change/EXTENSIONS_PLAN.md). Validation years only; 2025 is never read or re-scored."""

from __future__ import annotations

import json
import logging

import numpy as np
import pandas as pd

from aqf.data import ROOT
from aqf_land import forecast as fc
from aqf_land import forecast_run as fr

log = logging.getLogger("aqf_land.ext")
HMAX = fr.HMAX


def _setup(cfg: dict):
    d = fr.load(cfg)
    years = d.years()
    tr = fr._in(years, fr.TRAIN)
    clim = fc.climatology(d, tr)
    a = d.Y - clim
    return d, years, tr, clim, a


# ------------------------------------------------------------------ X2: CNN validation check
def x2_cnn(cfg: dict) -> dict:
    d, years, tr, clim, a = _setup(cfg)
    cnn_pred, note, epochs = fr.cnn_fit(d, a, years, fr.TRAIN, fr.VAL)
    gas = cfg.get("gas", "NO2")
    out = cfg["outputs"] / "X2_cnn"
    out.mkdir(parents=True, exist_ok=True)
    if cnn_pred is None:
        res = {"gas": gas, "status": note}
    else:
        tv, pv = fr.h1_rows(a, years, fr.VAL)
        mask = np.zeros(d.Y.shape, bool)
        mask[tv + 1, pv] = True
        pr = np.where(mask, cnn_pred(fr.VAL) + clim, np.nan)
        y = np.where(mask, d.Y, np.nan)
        s = fc.scores(y, pr, clim, d.dates, d.site_idx, np.random.default_rng(0))
        ladder = pd.read_csv(cfg["outputs"] / "C3_forecast" / "val_h1_ladder.csv").set_index("model")
        best = ladder.skill.idxmax()
        res = {"gas": gas, "status": note, "cnn_skill": s["skill"], "cnn_ci": [s["skill_lo"], s["skill_hi"]],
               "ridge_skill": float(ladder.loc["ridge", "skill"]), "chosen_model": best,
               "chosen_skill": float(ladder.loc[best, "skill"]),
               "cnn_beats_chosen": bool(s["skill"] > ladder.loc[best, "skill"]),
               "reading": "validation only; the recorded 2025 score stands. If the CNN wins, test deferred to 2026."}
    (out / "cnn_validation.json").write_text(json.dumps(res, indent=2))
    log.info("X2 %s: %s", gas, res)
    return res


# ------------------------------------------------------------------ X3: interval bands
def x3_bands(cfg: dict) -> pd.DataFrame:
    """Per lead: coverage of the 10-90 % band on 2024 targets, bands fitted on 2023 targets only.
    M0 global residual quantiles, M1 quantiles of residual / pixel training sd. Persistence-based forecast, g = 1."""
    d, years, tr, clim, a = _setup(cfg)
    alast = fc.last_valid(a)
    phi = fc.damping(a, alast, tr)
    sd_p = np.nanstd(np.where(tr[:, None], a, np.nan), 0)
    sd_p = np.where(np.isfinite(sd_p) & (sd_p > 0), sd_p, np.nanmean(sd_p))
    rows = []
    T = len(a)
    for h in range(1, HMAX + 1):
        ok = np.zeros(a.shape, bool)
        ok[:T - h] = np.isfinite(a[h:])
        t, p = np.nonzero(ok)
        tgt_year = years[t + h]
        A = phi[h] * alast[t, p]
        r = a[t + h, p] - A
        fit, ev = tgt_year == 2023, tgt_year == 2024
        q0 = np.percentile(r[fit], [10, 90])
        z = r / sd_p[p]
        q1 = np.percentile(z[fit], [10, 90])
        c0 = float(((r[ev] >= q0[0]) & (r[ev] <= q0[1])).mean())
        c1 = float(((z[ev] >= q1[0]) & (z[ev] <= q1[1])).mean())
        w0 = float(q0[1] - q0[0])
        w1 = float((sd_p[p][ev] * (q1[1] - q1[0])).mean())
        rows.append({"lead": h, "cov_M0": c0, "cov_M1": c1, "width_M0": w0, "width_M1": w1})
    r = pd.DataFrame(rows)
    gas = cfg.get("gas", "NO2")
    out = cfg["outputs"] / "X3_bands"
    out.mkdir(parents=True, exist_ok=True)
    r.to_csv(out / "bands_2023fit_2024eval.csv", index=False)
    l1 = r.iloc[0]
    log.info("X3 %s lead1 coverage M0 %.3f M1 %.3f; mean over leads M0 %.3f M1 %.3f", gas, l1.cov_M0, l1.cov_M1,
             r.cov_M0.mean(), r.cov_M1.mean())
    return r


# ------------------------------------------------------------------ X5: climate indices
def _monthly(path, skiprows=1, ncol=13) -> pd.Series:
    rows = {}
    for line in open(path):
        parts = line.split()
        if len(parts) == ncol:
            try:
                y = int(parts[0])
                vals = [float(v) for v in parts[1:]]
            except ValueError:
                continue
            for m, v in enumerate(vals, 1):
                if abs(v) < 90:                                      # -99.99 / -9999 mark missing
                    rows[pd.Timestamp(y, m, 1)] = v
    return pd.Series(rows).sort_index()


def index_at(dates: pd.DatetimeIndex, series: pd.Series, lag_months: int = 1) -> np.ndarray:
    """Index value of the month `lag_months` before each origin window's month (past only)."""
    keys = [(pd.Timestamp(d.year, d.month, 1) - pd.DateOffset(months=lag_months)) for d in dates]
    return np.array([series.get(k, np.nan) for k in keys])


def x5_climate(cfg: dict, alpha: float = 1.0, n_boot: int = 500) -> dict:
    d, years, tr, clim, a = _setup(cfg)
    oni = _monthly(ROOT / "data/climate/oni.data")
    dmi = _monthly(ROOT / "data/climate/dmi.data")
    X = np.column_stack([index_at(d.dates, oni), index_at(d.dates, dmi)])
    alast = fc.last_valid(a)
    phi = fc.damping(a, alast, tr)
    T, P = a.shape
    Xz = (X - np.nanmean(X[tr], 0)) / np.nanstd(X[tr], 0)
    site_mean = np.stack([np.nanmean(a[:, d.site_idx == s], 1) for s in range(len(d.sites))], 1)
    SM = {h: np.zeros(T) for h in range(1, HMAX + 1)}
    SC = {h: np.zeros(T) for h in range(1, HMAX + 1)}
    SD = {h: np.zeros(T) for h in range(1, HMAX + 1)}
    for h in range(2, HMAX + 1):
        coef = []
        for s in range(len(d.sites)):
            o = np.flatnonzero(tr[:T - h] & np.isfinite(site_mean[h:, s]) & np.isfinite(Xz[:T - h]).all(1))
            A_ = Xz[o]
            y_ = site_mean[o + h, s]
            w = np.linalg.solve(A_.T @ A_ + alpha * np.eye(2), A_.T @ (y_ - y_.mean())) if len(o) > 5 else np.zeros(2)
            coef.append(w)
        coef = np.array(coef)                                                 # (S, 2)
        pred_site = np.nan_to_num(Xz) @ coef.T                                # (T, S)
        for t in range(T - h):
            if years[t + h] not in (2023, 2024) or not np.isfinite(Xz[t]).all():
                continue
            y = a[t + h]
            ok = np.isfinite(y)
            pm = pred_site[t][d.site_idx]
            pd_ = phi[h] * alast[t]
            SM[h][t] = ((y - pm)[ok] ** 2).sum()
            SD[h][t] = ((y - pd_)[ok] ** 2).sum()
            SC[h][t] = (y[ok] ** 2).sum()
    leads = range(2, HMAX + 1)

    def mean_skill(sm, idx):
        return float(np.mean([1 - sm[h][idx].sum() / SC[h][idx].sum() for h in leads]))

    all_t = np.arange(T)
    rng = np.random.default_rng(0)
    boot = {"climate": [], "damped": []}
    block, nb = 6, int(np.ceil(T / 6))
    for _ in range(n_boot):
        starts = rng.integers(0, T - block, nb)
        idx = np.concatenate([np.arange(s, s + block) for s in starts])
        boot["climate"].append(mean_skill(SM, idx))
        boot["damped"].append(mean_skill(SD, idx))
    res = {"gas": cfg.get("gas", "NO2"), "alpha": alpha,
           "mean_skill_leads_2_30_climate_indices": mean_skill(SM, all_t),
           "ci": [float(np.percentile(boot["climate"], 2.5)), float(np.percentile(boot["climate"], 97.5))],
           "mean_skill_leads_2_30_damped_persistence": mean_skill(SD, all_t),
           "ci_damped": [float(np.percentile(boot["damped"], 2.5)), float(np.percentile(boot["damped"], 97.5))],
           "by_lead_skill": {h: float(1 - SM[h].sum() / SC[h].sum()) for h in (2, 6, 12, 18, 24, 30)}}
    res["adopt_for_2026"] = bool(res["ci"][0] > 0)
    out = cfg["outputs"] / "X5_climate"
    out.mkdir(parents=True, exist_ok=True)
    (out / "climate_indices_validation.json").write_text(json.dumps(res, indent=2))
    log.info("X5 %s: %s", res["gas"], {k: v for k, v in res.items() if k != "by_lead_skill"})
    return res
