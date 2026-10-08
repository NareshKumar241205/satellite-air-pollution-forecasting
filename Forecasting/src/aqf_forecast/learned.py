"""F2: learned models. LightGBM (mean + 10/90 % quantiles) and a Ridge reference, one per gas, pooled over pixels.

- Task A: a dedicated lead-1 model.
- Task B: one direct multi-horizon model with the lead as a feature (leads 1..30).
- Targets are standardised anomalies against the train climatology; missing targets are excluded from training.
- Fitting uses train-year targets. LightGBM early stopping and the Ridge alpha use validation 2023, so validation
  scores are optimistic; the single 2024 test (F5) decides.
- Interval and CRPS residuals come from leave-one-train-year-out refits with the chosen rounds/alpha, the same
  out-of-sample method as the baselines. Quantile LightGBM bands are scored separately.
"""

from __future__ import annotations

import json
import logging
import time

import lightgbm as lgb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from aqf import s2_static
from aqf_forecast import evaluate
from aqf_forecast.baselines import MODELS as BASELINES
from aqf_forecast.baselines import predict as predict_baseline
from aqf_forecast.features import build, pairs, rows_for
from aqf_forecast.metrics import target_obs
from aqf_forecast.plotting import GAS_COLORS, MODEL_STYLE, save
from aqf_forecast.problem import Problem, fit

log = logging.getLogger(__name__)
MODEL_STYLE.update({
    "lgbm": {"color": "#27ae60", "ls": "-", "lw": 1.8, "label": "LightGBM"},
    "ridge": {"color": "#8e44ad", "ls": "-", "lw": 1.4, "label": "Ridge"},
    "lgbm_no2020q2": {"color": "#16a085", "ls": ":", "lw": 1.6, "label": "LightGBM, 2020 Q2 excluded"},
    "ridge_no2020q2": {"color": "#c39bd3", "ls": ":", "lw": 1.6, "label": "Ridge, 2020 Q2 excluded"},
})


def _lgb_params(cfg: dict, objective: str = "regression", alpha: float | None = None, seed: int = 0) -> dict:
    keys = ("learning_rate", "num_leaves", "min_data_in_leaf", "feature_fraction", "bagging_fraction",
            "bagging_freq", "lambda_l2")
    prm = {k: cfg[k] for k in keys} | {"objective": objective, "verbosity": -1, "seed": seed,
                                       "deterministic": True, "num_threads": 0}
    if alpha is not None:
        prm["alpha"] = alpha
    return prm


def _fit_lgb(X, y, Xv, yv, cfg, objective="regression", alpha=None, rounds=None, seed=0):
    prm = _lgb_params(cfg, objective, alpha, seed)
    dtr = lgb.Dataset(X, y, free_raw_data=True)
    if rounds is not None:
        return lgb.train(prm, dtr, num_boost_round=rounds)
    dva = lgb.Dataset(Xv, yv, reference=dtr)
    return lgb.train(prm, dtr, num_boost_round=cfg["max_rounds"], valid_sets=[dva],
                     callbacks=[lgb.early_stopping(cfg["early_stopping"], verbose=False)])


class _Ridge:
    """Ridge on standardised features. NaNs (S2 over the sea) are filled with the train mean, plus an indicator."""

    def __init__(self, alpha: float):
        self.alpha = alpha

    def _prep(self, X):
        nan = np.isnan(X)
        Xf = np.where(nan, self.fill, X)
        return np.concatenate([(Xf - self.mu) / self.sd, nan[:, self.nan_cols].astype(np.float32)], 1)

    def fit(self, X, y):
        self.fill = np.nanmean(X, 0)
        self.nan_cols = np.nonzero(np.isnan(X).any(0))[0]
        Xf = np.where(np.isnan(X), self.fill, X)
        self.mu, self.sd = Xf.mean(0), Xf.std(0) + 1e-9
        Z = self._prep(X).astype(np.float64)
        self.b0 = y.mean()
        A = Z.T @ Z + self.alpha * np.eye(Z.shape[1])
        self.w = np.linalg.solve(A, Z.T @ (y - self.b0))
        return self

    def predict(self, X):
        return self._prep(X) @ self.w + self.b0


def _scatter(pred_rows, pair, origins, leads, fs, shape, clim_ext, g):
    """Put row predictions (standardised anomaly) back into a (T, H, G, Hp, Wp) value array for gas g.

    pair: the (origin, lead) pair index of each row, as returned by `rows_for`.
    """
    out = np.full(shape, np.nan)
    o, l = origins[pair], leads[pair]
    rows, cols = fs.pix
    P = len(rows)
    pix = np.tile(np.arange(P), len(origins))
    vals = clim_ext[o + l, g, rows[pix], cols[pix]] + pred_rows * fs.sd
    out[o, l - 1, g, rows[pix], cols[pix]] = vals
    return out


def run(p: Problem) -> dict:
    fc = p.fcfg
    cfg2 = fc["f2"]
    H = fc["horizon"]
    out = fc["outputs"] / "F2_learned"
    out.mkdir(parents=True, exist_ok=True)
    train = p.windows("train")
    val = p.windows("val")
    years = np.array([d.year for d in p.dates])
    f = fit(p, train)
    s2 = s2_static.build(p.cfg, p.cube.lat, p.cube.lon)
    shape = (p.T, H, len(p.gases)) + p.values.shape[2:]
    lg = cfg2["lightgbm"]

    base = {m: predict_baseline(m, p, f, H) for m in BASELINES}
    log.info("baseline residuals (leave one train year out)")
    resid = evaluate.loyo_residuals(p, BASELINES, H)
    learned = {k: np.full(shape, np.nan) for k in ("lgbm", "ridge", "lgbm_no2020q2", "ridge_no2020q2", "lgbmA",
                                                           "ridgeA")}
    bands = {k: (np.full(shape, np.nan), np.full(shape, np.nan)) for k in ("lgbm", "lgbmA")}
    info, importance = [], []

    for g, gas in enumerate(p.gases):
        fs = build(p, f, g, s2, train)
        for task, leads_list in (("A", [1]), ("B", list(range(1, H + 1)))):
            t0 = time.time()
            o_tr, l_tr = pairs(p, leads_list, train, cfg2["min_origin"])
            o_va, l_va = pairs(p, leads_list, val, cfg2["min_origin"])
            X, y, pr_tr = rows_for(fs, o_tr, l_tr)
            Xv, yv, pr_va = rows_for(fs, o_va, l_va)
            ok, okv = np.isfinite(y), np.isfinite(yv)
            # mean model + quantile bands
            m = _fit_lgb(X[ok], y[ok], Xv[okv], yv[okv], lg, seed=fc["seed"])
            rounds = m.best_iteration or lg["max_rounds"]
            q = [_fit_lgb(X[ok], y[ok], Xv[okv], yv[okv], lg, "quantile", a, seed=fc["seed"])
                 for a in fc["intervals"]["quantiles"]]
            # ridge reference, alpha chosen on validation
            best = None
            for alpha in cfg2["ridge_alphas"]:
                r = _Ridge(alpha).fit(X[ok], y[ok])
                mse = np.mean((r.predict(Xv[okv]) - yv[okv]) ** 2)
                if best is None or mse < best[0]:
                    best = (mse, alpha, r)
            _, alpha, rdg = best
            key_l, key_r = ("lgbm", "ridge") if task == "B" else ("lgbmA", "ridgeA")
            arrs = {key_l: m.predict(Xv), key_r: rdg.predict(Xv)}
            for k, v in arrs.items():
                filled = _scatter(v, pr_va, o_va, l_va, fs, shape, f.clim_ext, g)
                learned[k] = np.where(np.isfinite(filled), filled, learned[k])
            for bi, qm in enumerate(q):
                filled = _scatter(qm.predict(Xv), pr_va, o_va, l_va, fs, shape, f.clim_ext, g)
                bands[key_l][bi][...] = np.where(np.isfinite(filled), filled, bands[key_l][bi])
            # out-of-sample residuals: leave one train year out with the chosen rounds / alpha
            target_year = years[o_tr + l_tr]
            row_year = np.repeat(target_year, len(fs.pix[0]))
            row_lead = np.repeat(l_tr, len(fs.pix[0]))
            for k in (key_l, key_r):
                for h in leads_list:
                    resid.setdefault((k, g, h - 1), [])
            for yr in p.cfg["split"]["train_years"]:
                tr_rows = ok & (row_year != yr)
                te_rows = ok & (row_year == yr)
                mf = _fit_lgb(X[tr_rows], y[tr_rows], None, None, lg, rounds=rounds, seed=fc["seed"])
                rf = _Ridge(alpha).fit(X[tr_rows], y[tr_rows])
                for k, model in ((key_l, mf), (key_r, rf)):
                    e = (y[te_rows] - model.predict(X[te_rows])) * fs.sd
                    for h in leads_list:
                        resid[(k, g, h - 1)].append(e[row_lead[te_rows] == h])
            for k in (key_l, key_r):
                for h in leads_list:
                    resid[(k, g, h - 1)] = np.concatenate(resid[(k, g, h - 1)])
            imp = m.feature_importance("gain")
            importance += [{"gas": gas, "task": task, "feature": n, "gain_share": v / imp.sum()}
                           for n, v in zip(fs.names, imp)]
            info.append({"gas": gas, "task": task, "train_rows": int(ok.sum()), "val_rows": int(okv.sum()),
                         "lgbm_rounds": int(rounds), "ridge_alpha": alpha, "seconds": round(time.time() - t0, 1)})
            log.info("%s task %s: %d train rows, %d rounds, ridge alpha %g, %.0fs", gas, task, ok.sum(), rounds,
                     alpha, time.time() - t0)
            # 2020 Q2 check for the long-persistence signal (task B only)
            if task == "B":
                o_x, l_x = pairs(p, leads_list, train, cfg2["min_origin"], tuple(cfg2["exclude_2020q2"]))
                Xx, yx, _ = rows_for(fs, o_x, l_x)
                okx = np.isfinite(yx)
                mx = _fit_lgb(Xx[okx], yx[okx], Xv[okv], yv[okv], lg, seed=fc["seed"])
                filled = _scatter(mx.predict(Xv), pr_va, o_va, l_va, fs, shape, f.clim_ext, g)
                learned["lgbm_no2020q2"] = np.where(np.isfinite(filled), filled, learned["lgbm_no2020q2"])
                rx = _Ridge(alpha).fit(Xx[okx], yx[okx])
                filled = _scatter(rx.predict(Xv), pr_va, o_va, l_va, fs, shape, f.clim_ext, g)
                learned["ridge_no2020q2"] = np.where(np.isfinite(filled), filled, learned["ridge_no2020q2"])
                for h in leads_list:   # intervals for this check reuse the full models' residuals (not used for claims)
                    resid[("lgbm_no2020q2", g, h - 1)] = resid[("lgbm", g, h - 1)]
                    resid[("ridge_no2020q2", g, h - 1)] = resid[("ridge", g, h - 1)]

    pd.DataFrame(info).to_csv(out / "training_info.csv", index=False)
    imp = pd.DataFrame(importance)
    imp.to_csv(out / "feature_importance.csv", index=False)

    # Task A: lead-1 slice of baselines + dedicated lead-1 models
    predsA = {m: base[m][:, :1] for m in BASELINES} | {"lgbm": learned["lgbmA"][:, :1], "ridge": learned["ridgeA"][:, :1]}
    residA = {(m, g, 0): resid[(m if m in BASELINES else m + "A", g, 0)] for m in predsA for g in range(len(p.gases))}
    bandsA = {"lgbm": (bands["lgbmA"][0][:, :1], bands["lgbmA"][1][:, :1])}
    metA = evaluate.score(p, predsA, residA, "val", bandsA)
    metA.to_csv(out / "metrics_taskA_val.csv", index=False)
    a, _ = evaluate.summarise(metA, out, "val", prefix="A_")
    # Task B
    predsB = base | {k: learned[k] for k in ("lgbm", "ridge", "lgbm_no2020q2", "ridge_no2020q2")}
    metB = evaluate.score(p, predsB, resid, "val", {"lgbm": bands["lgbm"]})
    metB.to_csv(out / "metrics_taskB_val.csv", index=False)
    _, b = evaluate.summarise(metB, out, "val", prefix="B_")

    seas = _season_breakdown(p, predsA, predsB)
    seas.to_csv(out / "season_breakdown_val.csv", index=False)
    q2 = _q2_check(p, f, metB)
    q2.to_csv(out / "check_2020q2.csv", index=False)
    log.info("2020 Q2 check\n%s", q2.round(4).to_string(index=False))
    _figures(p, metB, imp, out)
    (out / "run_info.json").write_text(json.dumps({"lightgbm": lg, "ridge_alphas": cfg2["ridge_alphas"],
                                                   "features": build(p, f, 0, s2, train).names}, indent=2))
    return {"A": metA, "B": metB, "season": seas, "q2": q2, "importance": imp}


def _season_breakdown(p, predsA, predsB) -> pd.DataFrame:
    """Skill vs climatology per target season: task A (lead 1) and task B averaged over leads 1-30 (SSE pooled)."""
    from aqf.series import season_of

    rows = []
    for task, preds in (("A", predsA), ("B", predsB)):
        H = preds["climatology"].shape[1]
        obs = target_obs(p.values, H)
        tsplit = evaluate._target_split(p, H)
        tdates = p.target_dates(H)
        tseason = np.empty((p.T, H), dtype=object)
        for h in range(H):
            tseason[:, h] = season_of(tdates[np.arange(p.T) + h + 1].month, p.cfg)
        for g, gas in enumerate(p.gases):
            for season in p.cfg["seasons"]:
                sel = (tsplit == "val") & (tseason == season)
                o = obs[:, :, g][sel]
                c = preds["climatology"][:, :, g][sel]
                okc = np.isfinite(o) & np.isfinite(c)
                sse_c = np.sum((c - o)[okc] ** 2)
                for m, pm in preds.items():
                    if m in ("climatology", "persistence") or m.endswith("no2020q2"):
                        continue
                    pr = pm[:, :, g][sel]
                    ok = okc & np.isfinite(pr)
                    rows.append({"task": task, "gas": gas, "season": season, "model": m,
                                 "skill_vs_clim": 1 - np.sum((pr - o)[ok] ** 2) / np.sum((c - o)[ok] ** 2),
                                 "n": int(ok.sum())})
    return pd.DataFrame(rows)


def _q2_check(p, f, metB) -> pd.DataFrame:
    """Does the NO2 long-lead signal survive without 2020 Q2? Area-ACF with/without Q2, and LightGBM skill at long leads."""
    import warnings

    lo, hi = (pd.Timestamp(x) for x in p.fcfg["f2"]["exclude_2020q2"])
    train = p.windows("train")
    q2 = (p.dates >= lo) & (p.dates <= hi)
    rows = []
    for g, gas in enumerate(p.gases):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            area = np.nanmean(np.where(train[:, None, None], f.anom[:, g], np.nan).reshape(p.T, -1), 1)
        for label, mask in (("all train", train), ("train excl. 2020 Q2", train & ~q2)):
            r = []
            for L in (6, 12, 18):
                both = mask[:-L] & mask[L:]
                x, y = area[:-L][both], area[L:][both]
                ok = np.isfinite(x) & np.isfinite(y)
                r.append(np.corrcoef(x[ok], y[ok])[0, 1])
            d = metB[(metB.gas == gas)]
            sk = {m: d[(d.model == m) & (d.lead >= 7)].skill_vs_clim.mean()
                  for m in ("lgbm", "lgbm_no2020q2", "ridge", "ridge_no2020q2", "damped")}
            suffix = "" if label == "all train" else "_no2020q2"
            rows.append({"gas": gas, "data": label, "area_acf_30d": r[0], "area_acf_60d": r[1], "area_acf_90d": r[2],
                         "lgbm_mean_skill_leads7_30": sk["lgbm" + suffix],
                         "ridge_mean_skill_leads7_30": sk["ridge" + suffix],
                         "damped_mean_skill_leads7_30": sk["damped"]})
    return pd.DataFrame(rows)


def _figures(p, metB, imp, out):
    G = len(p.gases)
    fig, axes = plt.subplots(2, G, figsize=(14, 7), sharex=True)
    for g, gas in enumerate(p.gases):
        for m in ("climatology", "damped", "ridge", "lgbm", "ridge_no2020q2"):
            d = metB[(metB.gas == gas) & (metB.model == m)]
            st = MODEL_STYLE[m]
            axes[0, g].plot(d.lead_days, d.skill_vs_clim, **st)
            if m in ("ridge", "damped"):
                axes[0, g].fill_between(d.lead_days, d.skill_ci_lo, d.skill_ci_hi, color=st["color"], alpha=0.15)
            if not m.endswith("no2020q2"):
                axes[1, g].plot(d.lead_days, d.crps * p.cfg["display"][gas]["scale"], **st)
        axes[0, g].axhline(0, color="k", lw=0.8)
        lo = metB[metB.gas == gas].query("model in ['damped','lgbm','ridge']").skill_vs_clim.min()
        axes[0, g].set_ylim(max(-0.15, lo - 0.02), None)
        axes[0, g].set_title(f"{gas}: MSE skill vs climatology (95 % CI)")
        axes[1, g].set_title(f"{gas}: CRPS ({p.cfg['display'][gas]['unit']})")
        axes[1, g].set_xlabel("lead (days)")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("F2 learned models vs baselines by lead, validation 2023")
    fig.tight_layout()
    save(fig, out / "fig_f2_by_lead_val.png")

    fig, axes = plt.subplots(1, G, figsize=(15, 5))
    for g, gas in enumerate(p.gases):
        d = imp[(imp.gas == gas) & (imp.task == "B")].sort_values("gain_share").tail(12)
        axes[g].barh(d.feature, d.gain_share, color=GAS_COLORS[gas])
        axes[g].set_title(f"{gas}: LightGBM task B, top features (gain share)")
    fig.tight_layout()
    save(fig, out / "fig_feature_importance.png")
