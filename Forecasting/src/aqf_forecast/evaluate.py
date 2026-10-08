"""F1: evaluate the baselines for task A (lead 1) and task B (leads 1..30) on the validation year.

Fitting uses train years only. The 10/90 % intervals and CRPS samples come from out-of-sample residuals:
leave-one-train-year-out refits (fit on three train years, predict the fourth), so the interval width is not
estimated on data the model was fitted to. The test year is evaluated only with `split="test"` (F5, once).
"""

from __future__ import annotations

import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from aqf_forecast.baselines import MODELS, predict
from aqf_forecast.metrics import block_bootstrap_skill, crps_from_sorted, per_origin_sse, target_obs
from aqf_forecast.plotting import GAS_COLORS, MODEL_STYLE, save, scaled, unit
from aqf_forecast.problem import Problem, fit

log = logging.getLogger(__name__)


def _target_split(p: Problem, H: int) -> np.ndarray:
    ts = np.full((p.T, H), "", dtype=object)
    for h in range(1, H + 1):
        ts[: p.T - h, h - 1] = p.split[h:]
    return ts


def loyo_residuals(p: Problem, models: list[str], H: int) -> dict:
    """Out-of-sample residuals (obs - pred, mol/m^2) per (model, gas, lead) from leave-one-train-year-out fits."""
    obs = target_obs(p.values, H)
    years = np.array([d.year for d in p.dates])
    train_years = p.cfg["split"]["train_years"]
    tyear = np.full((p.T, H), -1)
    for h in range(1, H + 1):
        tyear[: p.T - h, h - 1] = years[h:]
    res = {(m, g, h): [] for m in models for g in range(len(p.gases)) for h in range(H)}
    for y in train_years:
        f = fit(p, p.windows("train") & (years != y))
        for m in models:
            e = obs - predict(m, p, f, H)                                # (T, H, G, Hp, Wp)
            for h in range(H):
                sel = tyear[:, h] == y
                for g in range(len(p.gases)):
                    v = e[sel, h, g].ravel()
                    res[(m, g, h)].append(v[np.isfinite(v)])
    return {k: np.concatenate(v) for k, v in res.items()}


def run(p: Problem, split: str = "val") -> pd.DataFrame:
    fc = p.fcfg
    H = fc["horizon"]
    out = fc["outputs"] / ("F1_baselines" if split == "val" else "F5_test_baselines")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(fc["seed"])
    q_lo, q_hi = fc["intervals"]["quantiles"]

    f = fit(p, p.windows("train"))
    preds = {m: predict(m, p, f, H) for m in MODELS}
    obs = target_obs(p.values, H)
    clim_t = preds["climatology"]
    sel_split = _target_split(p, H) == split
    log.info("computing out-of-sample residuals (leave one train year out)")
    resid = loyo_residuals(p, MODELS, H)

    rows, boot = [], []
    for g, gas in enumerate(p.gases):
        for h in range(H):
            orig = np.nonzero(sel_split[:, h])[0]
            o = obs[orig, h, g]
            c = clim_t[orig, h, g]
            sse_c, _ = per_origin_sse(c, o)
            for m in MODELS:
                pr = preds[m][orig, h, g]
                ok = np.isfinite(pr) & np.isfinite(o)
                err = (o - pr)[ok]
                pa, oa = (pr - c)[ok], (o - c)[ok]                       # predicted / observed anomaly
                r = resid[(m, g, h)]
                lo, hi = np.quantile(r, [q_lo, q_hi])
                rs = np.sort(rng.choice(r, size=min(fc["intervals"]["crps_samples"], len(r)), replace=False))
                # area mean over the pixels observed in each target map
                okm = np.isfinite(o)
                with warnings.catch_warnings():               # target maps with no valid pixel give NaN
                    warnings.simplefilter("ignore", RuntimeWarning)
                    area_o = np.nanmean(np.where(okm, o, np.nan).reshape(len(orig), -1), 1)
                    area_p = np.nanmean(np.where(okm, pr, np.nan).reshape(len(orig), -1), 1)
                    area_c = np.nanmean(np.where(okm, c, np.nan).reshape(len(orig), -1), 1)
                aok = np.isfinite(area_o)
                sse_m, _ = per_origin_sse(pr, o)
                rows.append({
                    "split": split, "gas": gas, "model": m, "lead": h + 1, "lead_days": 5 * (h + 1),
                    "n_origins": len(orig), "n_pixels": int(ok.sum()),
                    "rmse": float(np.sqrt(np.mean(err**2))), "mae": float(np.mean(np.abs(err))),
                    "bias": float(np.mean(-err)),
                    "skill_vs_clim": float(1 - sse_m.sum() / sse_c.sum()),
                    "acc": float(np.corrcoef(pa, oa)[0, 1]) if m != "climatology" else np.nan,
                    "coverage_80": float(np.mean((err >= lo) & (err <= hi))),
                    "interval_width": float(hi - lo),
                    "crps": float(np.mean(crps_from_sorted(err, rs))),
                    "area_rmse": float(np.sqrt(np.mean((area_p - area_o)[aok] ** 2))),
                    "area_skill_vs_clim": float(1 - np.sum((area_p - area_o)[aok] ** 2) / np.sum((area_c - area_o)[aok] ** 2)),
                })
                if m != "climatology":
                    b_lo, b_hi = block_bootstrap_skill(sse_m, sse_c, fc["bootstrap"]["block"], fc["bootstrap"]["n"],
                                                       fc["seed"])
                    boot.append({"gas": gas, "model": m, "lead": h + 1, "skill_ci_lo": b_lo, "skill_ci_hi": b_hi})
    met = pd.DataFrame(rows).merge(pd.DataFrame(boot), on=["gas", "model", "lead"], how="left")
    met.to_csv(out / f"metrics_{split}.csv", index=False)

    # Task A table (lead 1) and a compact task B summary
    a = met[met.lead == 1][["gas", "model", "rmse", "mae", "skill_vs_clim", "skill_ci_lo", "skill_ci_hi", "acc",
                            "coverage_80", "crps", "area_skill_vs_clim"]]
    a.to_csv(out / f"taskA_lead1_{split}.csv", index=False)
    b = met.groupby(["gas", "model"]).agg(mean_skill=("skill_vs_clim", "mean"), mean_crps=("crps", "mean"),
                                          mean_coverage=("coverage_80", "mean")).reset_index()
    def useful_leads(r) -> int:
        """Number of consecutive leads from lead 1 whose 95 % skill CI lies above 0 (isolated later hits are noise)."""
        d = met[(met.gas == r.gas) & (met.model == r.model)].sort_values("lead")
        above = (d.skill_ci_lo > 0).to_numpy()
        return int(np.argmin(above)) if (~above).any() else len(above)

    b["useful_leads"] = b.apply(useful_leads, axis=1)
    b["useful_days"] = 5 * b.useful_leads
    b.to_csv(out / f"taskB_summary_{split}.csv", index=False)
    log.info("task A (lead 1, %s)\n%s", split, a.round(3).to_string(index=False))
    log.info("task B summary (%s)\n%s", split, b.round(3).to_string(index=False))

    _figures(p, met, out, split)
    _example_maps(p, preds, obs, sel_split, out, split)
    return met


def _figures(p, met, out, split):
    G = len(p.gases)
    fig, axes = plt.subplots(3, G, figsize=(14, 9), sharex=True)
    for g, gas in enumerate(p.gases):
        for m in MODELS:
            d = met[(met.gas == gas) & (met.model == m)]
            st = MODEL_STYLE[m]
            axes[0, g].plot(d.lead_days, d.skill_vs_clim, **st)
            if m != "climatology":
                axes[0, g].fill_between(d.lead_days, d.skill_ci_lo, d.skill_ci_hi, color=st["color"], alpha=0.15)
            axes[1, g].plot(d.lead_days, d.crps * p.cfg["display"][gas]["scale"], **st)
            axes[2, g].plot(d.lead_days, d.coverage_80, **st)
        axes[0, g].axhline(0, color="k", lw=0.8)
        axes[0, g].set_ylim(max(-0.6, axes[0, g].get_ylim()[0]), None)
        axes[0, g].set_title(f"{gas}: MSE skill vs climatology (95 % CI)")
        axes[1, g].set_title(f"{gas}: CRPS ({unit(p.cfg, gas)}, lower is better)")
        axes[2, g].axhline(0.8, color="k", lw=0.8, ls="--")
        axes[2, g].set_ylim(0.5, 1)
        axes[2, g].set_title(f"{gas}: coverage of the 80 % interval")
        axes[2, g].set_xlabel("lead (days)")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f"Baselines by lead, {split} year (task A = 5 days, task B = 5–150 days)")
    fig.tight_layout()
    save(fig, out / f"fig_baselines_by_lead_{split}.png")


def _example_maps(p, preds, obs, sel_split, out, split):
    """One origin in the split: observed vs damped-persistence forecast at leads 1, 6, 18, 30."""
    leads = [1, 6, 18, 30]
    cand = np.nonzero(sel_split[:, max(leads) - 1])[0]
    if not len(cand):
        return
    t = int(cand[len(cand) // 2])
    for g, gas in enumerate(p.gases):
        fig, axes = plt.subplots(2, len(leads), figsize=(14, 6.4))
        fields = [obs[t, h - 1, g] for h in leads] + [preds["damped"][t, h - 1, g] for h in leads]
        vmin, vmax = np.nanpercentile(np.concatenate([x.ravel() for x in fields]), [2, 98])
        for k, h in enumerate(leads):
            for r, (name, src) in enumerate((("observed", obs), ("damped persistence", preds["damped"]))):
                ax = axes[r, k]
                im = ax.imshow(scaled(p.cfg, gas, np.where(p.cube.footprint, src[t, h - 1, g], np.nan))[:12, 1:],
                               cmap="YlOrRd", vmin=scaled(p.cfg, gas, vmin), vmax=scaled(p.cfg, gas, vmax))
                ax.set_title(f"{name}, +{5 * h} d ({(p.target_dates(30)[t + h]).date()})", fontsize=8)
                ax.set_xticks([])
                ax.set_yticks([])
                ax.grid(False)
        fig.colorbar(im, ax=axes, fraction=0.02).set_label(unit(p.cfg, gas))
        fig.suptitle(f"{gas}: forecast issued {p.dates[t].date()} (grey = missing)")
        save(fig, out / f"fig_example_maps_{gas}_{split}.png")
