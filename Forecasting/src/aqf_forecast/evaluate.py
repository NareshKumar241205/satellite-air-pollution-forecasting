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
    """F1: fit the baselines on train, score them on `split`, write tables and figures."""
    H = p.fcfg["horizon"]
    out = p.fcfg["outputs"] / ("F1_baselines" if split == "val" else "F5_test_baselines")
    out.mkdir(parents=True, exist_ok=True)
    f = fit(p, p.windows("train"))
    preds = {m: predict(m, p, f, H) for m in MODELS}
    log.info("computing out-of-sample residuals (leave one train year out)")
    resid = loyo_residuals(p, MODELS, H)
    met = score(p, preds, resid, split)
    met.to_csv(out / f"metrics_{split}.csv", index=False)
    summarise(met, out, split)
    _figures(p, met, out, split)
    _example_maps(p, preds, target_obs(p.values, H), _target_split(p, H) == split, out, split)
    return met


def score(p: Problem, preds: dict, resid: dict, split: str, bands: dict | None = None) -> pd.DataFrame:
    """Score every model in `preds` (values (T, H, G, Hp, Wp), mol/m^2) on target windows of `split`.

    resid[(model, g, lead_index)]: out-of-sample residuals giving the 10/90 % interval and CRPS samples.
    bands[model] = (lo, hi) arrays like preds: model-native quantile bands (e.g. quantile LightGBM), scored too.
    Skill is relative to climatology. When "damped" is present, each other model also gets its skill relative
    to damped persistence, with a block-bootstrap CI.
    """
    fc = p.fcfg
    H = next(iter(preds.values())).shape[1]
    rng = np.random.default_rng(fc["seed"])
    q_lo, q_hi = fc["intervals"]["quantiles"]
    obs = target_obs(p.values, H)
    clim_t = preds["climatology"]
    sel_split = _target_split(p, H) == split
    bs = (fc["bootstrap"]["block"], fc["bootstrap"]["n"], fc["seed"])

    rows = []
    for g, gas in enumerate(p.gases):
        for h in range(H):
            orig = np.nonzero(sel_split[:, h])[0]
            o = obs[orig, h, g]
            c = clim_t[orig, h, g]
            sse_c, _ = per_origin_sse(c, o)
            sse_d = per_origin_sse(preds["damped"][orig, h, g], o)[0] if "damped" in preds else None
            for m, pm in preds.items():
                pr = pm[orig, h, g]
                ok = np.isfinite(pr) & np.isfinite(o)
                if not ok.any():
                    continue
                err = (o - pr)[ok]
                pa, oa = (pr - c)[ok], (o - c)[ok]                       # predicted / observed anomaly
                r = resid[(m, g, h)]
                lo, hi = np.quantile(r, [q_lo, q_hi])
                rs = np.sort(rng.choice(r, size=min(fc["intervals"]["crps_samples"], len(r)), replace=False))
                okm = np.isfinite(o)
                with warnings.catch_warnings():                           # target maps with no valid pixel give NaN
                    warnings.simplefilter("ignore", RuntimeWarning)
                    area_o = np.nanmean(np.where(okm, o, np.nan).reshape(len(orig), -1), 1)
                    area_p = np.nanmean(np.where(okm, pr, np.nan).reshape(len(orig), -1), 1)
                    area_c = np.nanmean(np.where(okm, c, np.nan).reshape(len(orig), -1), 1)
                aok = np.isfinite(area_o)
                sse_m, _ = per_origin_sse(pr, o)
                row = {
                    "split": split, "gas": gas, "model": m, "lead": h + 1, "lead_days": 5 * (h + 1),
                    "n_origins": len(orig), "n_pixels": int(ok.sum()),
                    "rmse": float(np.sqrt(np.mean(err**2))), "mae": float(np.mean(np.abs(err))),
                    "bias": float(np.mean(-err)),
                    "skill_vs_clim": float(1 - sse_m.sum() / sse_c.sum()),
                    "acc": float(np.corrcoef(pa, oa)[0, 1]) if m != "climatology" and np.std(pa) > 0 else np.nan,
                    "coverage_80": float(np.mean((err >= lo) & (err <= hi))),
                    "interval_width": float(hi - lo),
                    "crps": float(np.mean(crps_from_sorted(err, rs))),
                    "area_rmse": float(np.sqrt(np.mean((area_p - area_o)[aok] ** 2))),
                    "area_skill_vs_clim": float(1 - np.sum((area_p - area_o)[aok] ** 2)
                                                / np.sum((area_c - area_o)[aok] ** 2)),
                }
                if m != "climatology":
                    row["skill_ci_lo"], row["skill_ci_hi"] = block_bootstrap_skill(sse_m, sse_c, *bs)
                if sse_d is not None and m not in ("climatology", "damped"):
                    row["skill_vs_damped"] = float(1 - sse_m.sum() / sse_d.sum())
                    row["vs_damped_ci_lo"], row["vs_damped_ci_hi"] = block_bootstrap_skill(sse_m, sse_d, *bs)
                if bands and m in bands:
                    blo, bhi = bands[m][0][orig, h, g][ok], bands[m][1][orig, h, g][ok]
                    oo = o[ok]
                    row["band_coverage_80"] = float(np.mean((oo >= blo) & (oo <= bhi)))
                    row["band_width"] = float(np.mean(bhi - blo))
                rows.append(row)
    return pd.DataFrame(rows)


def useful_leads(met: pd.DataFrame, gas: str, model: str, col: str = "skill_ci_lo") -> int:
    """Number of consecutive leads from lead 1 whose 95 % CI lower bound lies above 0 (later isolated hits are noise)."""
    d = met[(met.gas == gas) & (met.model == model)].sort_values("lead")
    if col not in d or d[col].isna().all():
        return 0
    above = (d[col] > 0).to_numpy()
    return int(np.argmin(above)) if (~above).any() else len(above)


def summarise(met: pd.DataFrame, out, split: str, prefix: str = "") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Task A table (lead 1) and task B summary (all leads)."""
    cols = ["gas", "model", "rmse", "mae", "skill_vs_clim", "skill_ci_lo", "skill_ci_hi", "acc", "coverage_80",
            "crps", "area_skill_vs_clim"] + [c for c in ("skill_vs_damped", "vs_damped_ci_lo", "vs_damped_ci_hi",
                                                          "band_coverage_80") if c in met]
    a = met[met.lead == 1][cols]
    a.to_csv(out / f"{prefix}taskA_lead1_{split}.csv", index=False)
    agg = {"mean_skill": ("skill_vs_clim", "mean"), "mean_crps": ("crps", "mean"), "mean_coverage": ("coverage_80", "mean")}
    if "band_coverage_80" in met:
        agg["mean_band_coverage"] = ("band_coverage_80", "mean")
    b = met.groupby(["gas", "model"]).agg(**agg).reset_index()
    b["useful_leads"] = [useful_leads(met, r.gas, r.model) for r in b.itertuples()]
    b["useful_days"] = 5 * b.useful_leads
    if "vs_damped_ci_lo" in met:
        b["leads_beating_damped"] = [useful_leads(met, r.gas, r.model, "vs_damped_ci_lo") for r in b.itertuples()]
    b.to_csv(out / f"{prefix}taskB_summary_{split}.csv", index=False)
    log.info("task A (lead 1, %s)\n%s", split, a.round(3).to_string(index=False))
    log.info("task B summary (%s)\n%s", split, b.round(3).to_string(index=False))
    return a, b


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
        cmap = plt.get_cmap("YlOrRd").copy()
        cmap.set_bad("#bdbdbd")                                  # missing pixels in grey
        fields = [obs[t, h - 1, g] for h in leads] + [preds["damped"][t, h - 1, g] for h in leads]
        vmin, vmax = np.nanpercentile(np.concatenate([x.ravel() for x in fields]), [2, 98])
        for k, h in enumerate(leads):
            for r, (name, src) in enumerate((("observed", obs), ("damped persistence", preds["damped"]))):
                ax = axes[r, k]
                im = ax.imshow(scaled(p.cfg, gas, np.where(p.cube.footprint, src[t, h - 1, g], np.nan))[:12, 1:],
                               cmap=cmap, vmin=scaled(p.cfg, gas, vmin), vmax=scaled(p.cfg, gas, vmax))
                ax.set_title(f"{name}, +{5 * h} d ({(p.target_dates(30)[t + h]).date()})", fontsize=8)
                ax.set_xticks([])
                ax.set_yticks([])
                ax.grid(False)
        fig.colorbar(im, ax=axes, fraction=0.02).set_label(unit(p.cfg, gas))
        fig.suptitle(f"{gas}: forecast issued {p.dates[t].date()} (grey = missing)")
        save(fig, out / f"fig_example_maps_{gas}_{split}.png")
