"""F4: ablations, rolling-origin robustness and a seasonal gate for the models chosen in F2.

Chosen models (Forecasting/REPORT.md, F2 decision): NO2 Ridge for tasks A and B, CO LightGBM for task A.
Their hyperparameters are **frozen** from the F2 run (`results/F2_learned/training_info.csv`: LightGBM rounds,
Ridge alpha), so nothing here is tuned on the year it is scored on.

1. Ablations on validation 2023. Each variant changes one thing, and its skill is compared with the full model
   (block-bootstrap CI on 1 - SSE_variant / SSE_full).
2. Rolling origin: for Y = 2021, 2022, 2023, refit climatology, scales and models on the years before Y and score Y.
   This tests whether a choice made on one validation year (2023, the highest-NO2 year) holds in other years.
3. Seasonal gate: use climatology instead of the model for targets in the NE monsoon (Oct-Dec), where F2 was weak.
"""

from __future__ import annotations

import dataclasses
import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from aqf import s2_static
from aqf.series import season_of
from aqf_forecast.baselines import predict as predict_baseline
from aqf_forecast.features import build, pairs, rows_for
from aqf_forecast.learned import _fit_lgb, _Ridge, _scatter
from aqf_forecast.metrics import block_bootstrap_skill, per_origin_sse, target_obs
from aqf_forecast.plotting import save
from aqf_forecast.problem import Problem, fit

log = logging.getLogger(__name__)

# (gas, task, model kind)
CHOSEN = [("NO2", "A", "ridge"), ("NO2", "B", "ridge"), ("CO", "A", "lgbm")]
VARIANTS = {
    "full": {},
    "no_s2": {"drop": ["s2_"]},
    "s2_shuffled": {"shuffle_s2": True},
    "no_neighbours": {"drop": ["neigh_"]},
    "no_missing_flags": {"drop": ["missing_", "n_valid_"]},
    "no_long_means": {"drop": ["area_mean", "pixel_mean"]},
    "multi_gas": {"other_gases": True},
    "no_2020q2": {"exclude_q2": True},
}
LEAD_GROUPS = {"A": {"lead 1 (5 d)": [1]},
               "B": {"lead 1 (5 d)": [1], "leads 2-6 (10-30 d)": list(range(2, 7)),
                     "leads 7-30 (35-150 d)": list(range(7, 31))}}


def _frozen(p: Problem) -> pd.DataFrame:
    path = p.fcfg["outputs"] / "F2_learned" / "training_info.csv"
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing: run `aqf-forecast learned` (F2) first")
    return pd.read_csv(path).set_index(["gas", "task"])


def fit_predict(p: Problem, gas: str, task: str, kind: str, hp: dict, train_mask: np.ndarray, eval_mask: np.ndarray,
                opts: dict | None = None, seed: int = 0) -> np.ndarray:
    """Fit one chosen model on targets in train_mask (climatology and scales fitted on those windows too) and
    return predicted values (T, H, G, Hp, Wp) for every pair whose target lies in eval_mask (NaN elsewhere)."""
    opts = opts or {}
    fc = p.fcfg
    cfg2 = fc["f2"]
    g = p.gases.index(gas)
    leads = [1] if task == "A" else list(range(1, fc["horizon"] + 1))
    H = len(leads) if task == "A" else fc["horizon"]
    f = fit(p, train_mask)
    s2 = s2_static.build(p.cfg, p.cube.lat, p.cube.lon)
    fs = build(p, f, g, s2, train_mask, other_gases=opts.get("other_gases"))
    if opts.get("shuffle_s2"):
        cols = [i for i, n in enumerate(fs.names[-fs.static.shape[1] - 1:-1]) if n.startswith("s2_")]
        perm = np.random.default_rng(seed).permutation(fs.static.shape[0])
        fs.static[:, cols] = fs.static[perm][:, cols]
    keep = np.array([not any(n.startswith(d) for d in opts.get("drop", [])) for n in fs.names])
    exclude = tuple(cfg2["exclude_2020q2"]) if opts.get("exclude_q2") else None
    o_tr, l_tr = pairs(p, leads, train_mask, cfg2["min_origin"], exclude)
    o_ev, l_ev = pairs(p, leads, eval_mask, cfg2["min_origin"])
    X, y, _ = rows_for(fs, o_tr, l_tr)
    Xe, _, pr_ev = rows_for(fs, o_ev, l_ev)
    ok = np.isfinite(y)
    X, Xe = X[:, keep], Xe[:, keep]
    if kind == "ridge":
        model = _Ridge(hp["ridge_alpha"]).fit(X[ok], y[ok])
    else:
        model = _fit_lgb(X[ok], y[ok], None, None, cfg2["lightgbm"], rounds=int(hp["lgbm_rounds"]), seed=seed)
    shape = (p.T, H, len(p.gases)) + p.values.shape[2:]
    return _scatter(model.predict(Xe), pr_ev, o_ev, l_ev, fs, shape, f.clim_ext, g)


def _skills(p: Problem, preds: dict, g: int, eval_split: np.ndarray, leads: list[int], ref: str) -> dict:
    """Skill of every model vs `ref`, pooled over `leads`, with block-bootstrap CIs over origins."""
    H = next(iter(preds.values())).shape[1]
    obs = target_obs(p.values, H)
    tsplit = np.zeros((p.T, H), bool)
    for h in range(1, H + 1):
        tsplit[: p.T - h, h - 1] = eval_split[h:]
    origins = np.nonzero(tsplit[:, np.array(leads) - 1].any(1))[0]     # every origin with a scored target
    sse = {}
    for m, pm in preds.items():
        tot = np.zeros(len(origins))
        for h in leads:
            # all pairs at this lead whose target is in the evaluation years, aligned on the origin index
            o_h = np.nonzero(tsplit[:, h - 1])[0]
            s, _ = per_origin_sse(pm[o_h, h - 1, g], obs[o_h, h - 1, g])
            tot += pd.Series(s, index=o_h).reindex(origins, fill_value=0.0).to_numpy()
        sse[m] = tot
    bs = (p.fcfg["bootstrap"]["block"], p.fcfg["bootstrap"]["n"], p.fcfg["seed"])
    out = {}
    for m in preds:
        if m == ref:
            continue
        lo, hi = block_bootstrap_skill(sse[m], sse[ref], *bs)
        out[m] = (float(1 - sse[m].sum() / sse[ref].sum()), lo, hi)
    return out


def run(p: Problem) -> dict:
    out = p.fcfg["outputs"] / "F4_ablations"
    out.mkdir(parents=True, exist_ok=True)
    hp_all = _frozen(p)
    years = np.array([d.year for d in p.dates])
    train, val = p.windows("train"), p.windows("val")

    # 1. ablations on validation 2023
    rows = []
    for gas, task, kind in CHOSEN:
        g = p.gases.index(gas)
        hp = hp_all.loc[(gas, task)].to_dict()
        preds = {}
        for name, opts in VARIANTS.items():
            preds[name] = fit_predict(p, gas, task, kind, hp, train, val, opts, p.fcfg["seed"])
            log.info("ablation %s %s %s: %s", gas, task, kind, name)
        H = preds["full"].shape[1]
        f = fit(p, train)
        clim = predict_baseline("climatology", p, f, H)
        for grp, leads in LEAD_GROUPS[task].items():
            vs_full = _skills(p, preds, g, val, leads, "full")
            vs_clim = _skills(p, preds | {"climatology": clim}, g, val, leads, "climatology")
            for name in VARIANTS:
                r = {"gas": gas, "task": task, "model": kind, "leads": grp, "variant": name,
                     "skill_vs_clim": vs_clim[name][0], "clim_ci_lo": vs_clim[name][1], "clim_ci_hi": vs_clim[name][2]}
                if name != "full":
                    r |= {"skill_vs_full": vs_full[name][0], "full_ci_lo": vs_full[name][1], "full_ci_hi": vs_full[name][2]}
                rows.append(r)
    abl = pd.DataFrame(rows)
    abl.to_csv(out / "ablations_val.csv", index=False)
    log.info("ablations\n%s", abl.round(4).to_string(index=False))

    # 2. rolling origin + 3. seasonal gate, for the full chosen models
    rows, gate_rows = [], []
    for Y in (2021, 2022, 2023):
        tr = years < Y
        ev = years == Y
        pY = dataclasses.replace(p, split=np.where(tr, "train", np.where(ev, "val", "other")))
        fY = fit(pY, tr)
        for gas, task, kind in CHOSEN:
            g = p.gases.index(gas)
            hp = hp_all.loc[(gas, task)].to_dict()
            model = fit_predict(pY, gas, task, kind, hp, tr, ev, None, p.fcfg["seed"])
            H = model.shape[1]
            clim = predict_baseline("climatology", pY, fY, H)
            damped = predict_baseline("damped", pY, fY, H)
            tdates = pY.target_dates(H)
            ne = np.zeros((p.T, H), bool)
            for h in range(1, H + 1):
                ne[:, h - 1] = season_of(tdates[np.arange(p.T) + h].month, p.cfg) == "NE monsoon"
            gated = np.where(ne[:, :, None, None, None], clim, model)
            preds = {"model": model, "model_gated": gated, "damped": damped, "climatology": clim}
            for grp, leads in LEAD_GROUPS[task].items():
                vs_c = _skills(pY, preds, g, ev, leads, "climatology")
                vs_d = _skills(pY, {k: v for k, v in preds.items() if k != "climatology"}, g, ev, leads, "damped")
                rows.append({"year": Y, "train_years": f"2019-{Y - 1}", "gas": gas, "task": task, "model": kind,
                             "leads": grp,
                             "model_skill": vs_c["model"][0], "model_ci_lo": vs_c["model"][1], "model_ci_hi": vs_c["model"][2],
                             "damped_skill": vs_c["damped"][0], "damped_ci_lo": vs_c["damped"][1],
                             "model_vs_damped": vs_d["model"][0], "vs_damped_ci_lo": vs_d["model"][1],
                             "vs_damped_ci_hi": vs_d["model"][2]})
                gate_rows.append({"year": Y, "gas": gas, "task": task, "leads": grp,
                                  "model_skill": vs_c["model"][0], "gated_skill": vs_c["model_gated"][0],
                                  "gate_gain": vs_c["model_gated"][0] - vs_c["model"][0]})
            log.info("rolling %d %s %s done", Y, gas, task)
    roll = pd.DataFrame(rows)
    roll.to_csv(out / "rolling_origin.csv", index=False)
    gate = pd.DataFrame(gate_rows)
    gate.to_csv(out / "seasonal_gate.csv", index=False)
    log.info("rolling origin\n%s", roll.round(4).to_string(index=False))
    log.info("seasonal gate\n%s", gate.round(4).to_string(index=False))
    _figures(abl, roll, out)
    return {"ablations": abl, "rolling": roll, "gate": gate}


def _figures(abl, roll, out):
    sel = abl[abl.variant != "full"]
    groups = list(sel.groupby(["gas", "task", "leads"], sort=False))
    fig, axes = plt.subplots(1, len(groups), figsize=(4.2 * len(groups), 4), sharey=True)
    for ax, ((gas, task, leads), d) in zip(np.atleast_1d(axes), groups):
        y = np.arange(len(d))
        err = np.vstack([d.skill_vs_full - d.full_ci_lo, d.full_ci_hi - d.skill_vs_full])
        ax.barh(y, d.skill_vs_full * 100, xerr=err * 100, color=np.where(d.full_ci_hi < 0, "#c0392b", "#95a5a6"))
        ax.set_yticks(y, d.variant)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_title(f"{gas} task {task}, {leads}", fontsize=9)
        ax.set_xlabel("skill vs full model (%); red = significantly worse")
    fig.suptitle("F4 ablations on validation 2023 (one change at a time, frozen hyperparameters)")
    fig.tight_layout()
    save(fig, out / "fig_ablations.png")

    groups = list(roll.groupby(["gas", "task", "leads"], sort=False))
    fig, axes = plt.subplots(1, len(groups), figsize=(4.2 * len(groups), 3.6), sharey=False)
    for ax, ((gas, task, leads), d) in zip(np.atleast_1d(axes), groups):
        x = np.arange(len(d))
        ax.bar(x - 0.18, d.model_skill * 100, 0.36, label="chosen model",
               yerr=np.vstack([d.model_skill - d.model_ci_lo, d.model_ci_hi - d.model_skill]) * 100, color="#8e44ad")
        ax.bar(x + 0.18, d.damped_skill * 100, 0.36, label="damped persistence", color="#2471a3")
        ax.set_xticks(x, d.year)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(f"{gas} task {task}, {leads}", fontsize=9)
        ax.set_ylabel("skill vs climatology (%)")
    np.atleast_1d(axes)[0].legend(fontsize=7)
    fig.suptitle("Rolling origin: train on the years before Y, score Y")
    fig.tight_layout()
    save(fig, out / "fig_rolling_origin.png")
