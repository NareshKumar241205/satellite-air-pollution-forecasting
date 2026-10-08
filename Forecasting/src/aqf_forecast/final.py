"""F5: the single 2024 test of the frozen choice, then the operational Jan–May 2025 outlook.

1. Test: every component is fitted on 2019–2023 (climatology too) with the frozen F2 hyperparameters and scored
   once on 2024. The frozen composite (`final.frozen` in configs/forecast.yaml) is scored together with the
   alternatives, which are listed for the record only. A lock file stops a silent second scoring.
2. Outlook: everything is refitted on 2019–2024 (climatology recomputed on those years). The forecast is issued
   from the last window (25 Dec 2024) for 30 steps, 30 Dec 2024 to 29 May 2025.
   - 80 % bands come from leave-one-year-out residual quantiles per gas and lead: pixel-level for the maps,
     area-level for the area-mean outlook.
   - Every output names the model that produced each lead.
"""

from __future__ import annotations

import dataclasses
import json
import logging
import warnings
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

from aqf.data import nearest_pixel
from aqf_forecast import evaluate
from aqf_forecast.ablations import _frozen, _skills, fit_predict
from aqf_forecast.baselines import predict as predict_baseline
from aqf_forecast.metrics import target_obs
from aqf_forecast.plotting import GAS_COLORS, save, scaled, unit
from aqf_forecast.problem import STEP_DAYS, Problem, fit

log = logging.getLogger(__name__)
LEARNED = {"ridgeA": ("A", "ridge"), "ridgeB": ("B", "ridge"), "lgbmA": ("A", "lgbm")}
GROUPS = {"5 d (lead 1)": [1], "10-30 d (leads 2-6)": list(range(2, 7)), "35-150 d (leads 7-30)": list(range(7, 31)),
          "all 5-150 d": list(range(1, 31))}


def lead_models(p: Problem, gas: str) -> list[str]:
    """Component name for each lead 1..H under the frozen choice."""
    H = p.fcfg["horizon"]
    out = [None] * H
    for a, b, comp in p.fcfg["final"]["frozen"][gas]:
        for h in range(a, b + 1):
            out[h - 1] = comp
    assert all(out), f"frozen spec for {gas} leaves leads uncovered"
    return out


def component(p: Problem, gas: str, comp: str, train_mask, eval_mask, hp_all, eval_origins=None) -> np.ndarray:
    """Values (T, H, G, Hp, Wp) of one component fitted on train_mask; only gas `gas` is filled for learned models."""
    H = p.fcfg["horizon"]
    if comp in ("climatology", "damped", "persistence"):
        return predict_baseline(comp, p, fit(p, train_mask), H)
    task, kind = LEARNED[comp]
    hp = hp_all.loc[(gas, task)].to_dict()
    pr = fit_predict(p, gas, task, kind, hp, train_mask, eval_mask, None, p.fcfg["seed"], eval_origins)
    if pr.shape[1] < H:                                   # task-A models only exist at lead 1
        pad = np.full((pr.shape[0], H - pr.shape[1]) + pr.shape[2:], np.nan)
        pr = np.concatenate([pr, pad], 1)
    return pr


def compose(p: Problem, comps: dict[tuple[str, str], np.ndarray]) -> np.ndarray:
    """Build the frozen composite: for each gas and lead, take the component the frozen choice names."""
    H = p.fcfg["horizon"]
    out = np.full((p.T, H, len(p.gases)) + p.values.shape[2:], np.nan)
    for g, gas in enumerate(p.gases):
        for h, comp in enumerate(lead_models(p, gas)):
            out[:, h, g] = comps[(gas, comp)][:, h, g]
    return out


def needed(p: Problem, with_alternatives: bool) -> list[tuple[str, str]]:
    need = {(gas, c) for gas in p.gases for c in lead_models(p, gas)}
    if with_alternatives:
        need |= {(gas, c) for gas, cs in p.fcfg["final"]["alternatives"].items() for c in cs}
    need |= {(gas, c) for gas in p.gases for c in ("climatology", "damped")}
    return sorted(need)


def loyo(p: Problem, years_in: list[int], hp_all, with_alternatives: bool) -> dict:
    """Out-of-sample predictions for every needed component: each year in years_in is predicted by a fit on the others."""
    years = np.array([d.year for d in p.dates])
    H = p.fcfg["horizon"]
    tyear = np.full((p.T, H), -1)
    for h in range(1, H + 1):
        tyear[: p.T - h, h - 1] = years[h:]
    out = {}
    for gas, comp in needed(p, with_alternatives):
        acc = None
        for y in years_in:
            tr = np.isin(years, [v for v in years_in if v != y])
            pr = component(p, gas, comp, tr, years == y, hp_all)
            keep = (tyear == y)[:, :, None, None, None]
            acc = np.where(keep, pr, np.nan) if acc is None else np.where(keep, pr, acc)
        out[(gas, comp)] = acc
        log.info("LOYO %s %s over %s", gas, comp, years_in)
    return out


def residuals(p: Problem, preds: dict, keys: dict[str, tuple[str, str] | None]) -> tuple[dict, dict]:
    """Pixel and area-mean residuals (obs - pred) per (model label, gas index, lead index)."""
    H = p.fcfg["horizon"]
    obs = target_obs(p.values, H)
    pix, area = {}, {}
    for label, arr in preds.items():
        for g in range(len(p.gases)):
            e = obs[:, :, g] - arr[:, :, g]                      # (T, H, Hp, Wp)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                ok = np.isfinite(obs[:, :, g]) & np.isfinite(arr[:, :, g])
                ao = np.nanmean(np.where(ok, obs[:, :, g], np.nan).reshape(p.T, H, -1), 2)
                ap = np.nanmean(np.where(ok, arr[:, :, g], np.nan).reshape(p.T, H, -1), 2)
            for h in range(H):
                v = e[:, h].ravel()
                pix[(label, g, h)] = v[np.isfinite(v)]
                a = (ao - ap)[:, h]
                area[(label, g, h)] = a[np.isfinite(a)]
    return pix, area


# --------------------------------------------------------------------------------------------------- test
def test(p: Problem, rescore: bool = False) -> pd.DataFrame:
    fc = p.fcfg
    out = fc["outputs"] / "F5_final"
    out.mkdir(parents=True, exist_ok=True)
    lock = out / "TEST_SCORED.json"
    if lock.exists() and not rescore:
        raise RuntimeError(f"2024 was already scored ({json.loads(lock.read_text())['scored_at']}). "
                           "The test is scored once. Use --rescore-test only to reproduce the identical numbers.")
    Y = fc["final"]["test_year"]
    years = np.array([d.year for d in p.dates])
    hp_all = _frozen(p)
    pT = dataclasses.replace(p, split=np.where(years < Y, "train", np.where(years == Y, "test", "")))
    tr, te = years < Y, years == Y

    comps = {k: component(pT, k[0], k[1], tr, te, hp_all) for k in needed(pT, True)}
    preds = {"climatology": comps[(pT.gases[0], "climatology")], "damped": comps[(pT.gases[0], "damped")],
             "frozen": compose(pT, comps)}
    for gas, cs in fc["final"]["alternatives"].items():
        for c in cs:
            arr = preds.setdefault(c, np.full_like(preds["frozen"], np.nan))
            g = pT.gases.index(gas)
            arr[:, :, g] = comps[(gas, c)][:, :, g]
    # interval residuals: leave one year out over the fitting years, same components
    oos = loyo(pT, sorted(set(years[tr])), hp_all, True)
    oos_preds = {"climatology": oos[(pT.gases[0], "climatology")], "damped": oos[(pT.gases[0], "damped")],
                 "frozen": compose(pT, oos)}
    for gas, cs in fc["final"]["alternatives"].items():
        for c in cs:
            arr = oos_preds.setdefault(c, np.full_like(oos_preds["frozen"], np.nan))
            g = pT.gases.index(gas)
            arr[:, :, g] = oos[(gas, c)][:, :, g]
    resid, _ = residuals(pT, oos_preds, {})
    for k in list(preds):
        for g in range(len(pT.gases)):
            for h in range(fc["horizon"]):
                resid.setdefault((k, g, h), np.array([0.0]))

    met = evaluate.score(pT, preds, resid, "test")
    met = met[met.n_pixels > 0]
    met.to_csv(out / "test_metrics_by_lead.csv", index=False)
    rows = []
    for g, gas in enumerate(pT.gases):
        models = ["frozen", "damped"] + fc["final"]["alternatives"].get(gas, [])
        for grp, leads in GROUPS.items():
            # a model is summarised over a lead group only if it forecasts every lead in it (task-A models: lead 1)
            has = [m for m in models if all(np.isfinite(preds[m][:, l - 1, g]).any() for l in leads)]
            sk = _skills(pT, {m: preds[m] for m in has} | {"climatology": preds["climatology"]}, g, te, leads,
                         "climatology")
            cov = met[(met.gas == gas) & met.lead.isin(leads)].groupby("model").coverage_80.mean()
            for m in has:
                if m in sk and np.isfinite(sk[m][0]):
                    rows.append({"gas": gas, "leads": grp, "model": m, "skill_vs_clim": sk[m][0],
                                 "ci_lo": sk[m][1], "ci_hi": sk[m][2], "coverage_80": cov.get(m, np.nan),
                                 "frozen_components": ", ".join(sorted(set(lead_models(pT, gas)[l - 1] for l in leads)))
                                 if m == "frozen" else m})
    summ = pd.DataFrame(rows)
    summ.to_csv(out / "test_summary.csv", index=False)
    lock.write_text(json.dumps({"scored_at": datetime.now().isoformat(timespec="seconds"), "test_year": Y,
                                "frozen": fc["final"]["frozen"], "alternatives": fc["final"]["alternatives"]},
                               indent=2))
    log.info("2024 TEST (scored once)\n%s", summ.round(4).to_string(index=False))
    _test_figure(pT, met, out)
    return summ


def _test_figure(p, met, out):
    fig, axes = plt.subplots(1, len(p.gases), figsize=(15, 3.8))
    sty = {"frozen": ("#c0392b", "-", 2.2, "frozen choice"), "damped": ("#2471a3", "-", 1.2, "damped persistence"),
           "ridgeB": ("#8e44ad", ":", 1.4, "Ridge (task B), alternative"),
           "ridgeA": ("#8e44ad", "o", 1.4, "Ridge (task A), alternative"),
           "lgbmA": ("#27ae60", "o", 1.4, "LightGBM (task A), alternative")}
    for ax, gas in zip(axes, p.gases):
        for m, (c, ls, lw, lab) in sty.items():
            d = met[(met.gas == gas) & (met.model == m)]
            if d.empty:
                continue
            if ls == "o":
                ax.plot(d.lead_days, d.skill_vs_clim, "o", color=c, label=lab)
            else:
                ax.plot(d.lead_days, d.skill_vs_clim, ls, color=c, lw=lw, label=lab)
                if m == "frozen":
                    ax.fill_between(d.lead_days, d.skill_ci_lo, d.skill_ci_hi, color=c, alpha=0.15)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(f"{gas}: 2024 test, skill vs climatology (95 % CI)", fontsize=9)
        ax.set_xlabel("lead (days)")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    save(fig, out / "fig_test_2024_by_lead.png")


# ------------------------------------------------------------------------------------------------ outlook
def outlook(p: Problem) -> pd.DataFrame:
    fc = p.fcfg
    H = fc["horizon"]
    out = fc["outputs"] / "F5_final" / "outlook_2025"
    out.mkdir(parents=True, exist_ok=True)
    years = np.array([d.year for d in p.dates])
    hp_all = _frozen(p)
    all_w = np.ones(p.T, bool)
    origin = np.array([p.T - 1])
    comps = {k: component(p, k[0], k[1], all_w, all_w, hp_all, eval_origins=origin) for k in needed(p, False)}
    fc_vals = compose(p, comps)[p.T - 1]                       # (H, G, Hp, Wp)
    oos = loyo(p, sorted(set(years)), hp_all, False)
    pix_r, area_r = residuals(p, {"frozen": compose(p, oos)}, {})
    q_lo, q_hi = fc["intervals"]["quantiles"]
    tdates = p.target_dates(H)[p.T:]
    fp = p.cube.footprint
    lo = np.full_like(fc_vals, np.nan)
    hi = np.full_like(fc_vals, np.nan)
    rows = []
    for g, gas in enumerate(p.gases):
        models = lead_models(p, gas)
        for h in range(H):
            r_lo, r_hi = np.quantile(pix_r[("frozen", g, h)], [q_lo, q_hi])
            lo[h, g] = np.maximum(fc_vals[h, g] + r_lo, 0)        # columns are non-negative
            hi[h, g] = fc_vals[h, g] + r_hi
            a_lo, a_hi = np.quantile(area_r[("frozen", g, h)], [q_lo, q_hi])
            mean = float(np.nanmean(fc_vals[h, g][fp]))
            rows.append({"gas": gas, "lead": h + 1, "lead_days": STEP_DAYS * (h + 1),
                         "window_start": tdates[h].date(), "window_end": (tdates[h] + pd.Timedelta(days=STEP_DAYS)).date(),
                         "model": models[h], "area_mean_mol_m2": mean,
                         "area_p10_mol_m2": max(mean + a_lo, 0.0), "area_p90_mol_m2": mean + a_hi,
                         "area_mean_display": scaled(p.cfg, gas, mean),
                         "area_p10_display": scaled(p.cfg, gas, max(mean + a_lo, 0.0)),
                         "area_p90_display": scaled(p.cfg, gas, mean + a_hi), "display_unit": unit(p.cfg, gas)})
        _geotiff(p, out / f"outlook_{gas}_2025.tif", gas, fc_vals[:, g], lo[:, g], hi[:, g], models, tdates)
    tab = pd.DataFrame(rows)
    tab.to_csv(out / "outlook_area_mean.csv", index=False)

    hot = []
    for gas, hub in fc["final"]["hotspots"]:
        g = p.gases.index(gas)
        i, j = nearest_pixel(p.cube, *p.cfg["hubs"][hub])
        for h in range(H):
            hot.append({"gas": gas, "hub": hub, "pixel_row": i, "pixel_col": j, "lead_days": STEP_DAYS * (h + 1),
                        "window_start": tdates[h].date(), "model": lead_models(p, gas)[h],
                        "mean_display": scaled(p.cfg, gas, fc_vals[h, g, i, j]),
                        "p10_display": scaled(p.cfg, gas, lo[h, g, i, j]),
                        "p90_display": scaled(p.cfg, gas, hi[h, g, i, j]), "display_unit": unit(p.cfg, gas)})
    hot = pd.DataFrame(hot)
    hot.to_csv(out / "outlook_hotspots.csv", index=False)
    np.savez_compressed(out / "outlook_maps.npz", mean=fc_vals, p10=lo, p90=hi, gases=np.array(p.gases),
                        window_start=np.array([str(d.date()) for d in tdates]))
    _fan_charts(p, tab, hot, out)
    _maps(p, fc_vals, out, tdates)
    (out / "README.md").write_text(_readme(p, tdates))
    log.info("outlook\n%s", tab[tab.lead.isin([1, 2, 6, 7, 18, 30])][
        ["gas", "lead_days", "window_start", "model", "area_mean_display", "area_p10_display", "area_p90_display",
         "display_unit"]].round(2).to_string(index=False))
    return tab


def _geotiff(p, path, gas, mean, lo, hi, models, tdates):
    """One GeoTIFF per gas: bands 1-30 mean, 31-60 10 % quantile, 61-90 90 % quantile (mol/m^2, NaN = no data)."""
    with rasterio.open(p.cfg["paths"]["s5p_dir"] / sorted(p.cfg["paths"]["s5p_dir"].glob("s5p_*.tif"))[0].name) as src:
        prof = src.profile
    H = len(models)
    prof.update(count=3 * H, dtype="float32", nodata=np.nan, compress="deflate")
    fp = p.cube.footprint
    with rasterio.open(path, "w", **prof) as dst:
        for k, (name, arr) in enumerate((("mean", mean), ("p10", lo), ("p90", hi))):
            for h in range(H):
                b = k * H + h + 1
                dst.write(np.where(fp, arr[h], np.nan).astype("float32"), b)
                dst.set_band_description(b, f"{gas} {name} +{STEP_DAYS * (h + 1)}d {tdates[h].date()} [{models[h]}]")
                dst.update_tags(b, statistic=name, lead_days=STEP_DAYS * (h + 1), window_start=str(tdates[h].date()),
                                model=models[h])
        dst.update_tags(gas=gas, unit="mol/m^2", issued_from_window="2024-12-25",
                        layout="bands 1-30 mean, 31-60 p10, 61-90 p90", interval="80% (leave-one-year-out residuals)")


def _fan_charts(p, tab, hot, out):
    colors = {"damped": "#2471a3", "ridgeB": "#8e44ad", "climatology": "#7f8c8d"}
    years = np.array([d.year for d in p.dates])
    from aqf.series import window_table

    wt = window_table(p.cube, p.cfg)
    fig, axes = plt.subplots(len(p.gases), 1, figsize=(11, 9), sharex=True)
    for ax, gas in zip(axes, p.gases):
        h = wt[(wt.gas == gas) & wt.usable & (wt.date >= "2023-07-01")]
        ax.plot(h.date, scaled(p.cfg, gas, h["mean"]), "-", color="k", lw=0.9, label="observed area mean")
        d = tab[tab.gas == gas]
        x = pd.to_datetime(d.window_start)
        ax.fill_between(x, d.area_p10_display, d.area_p90_display, color=GAS_COLORS[gas], alpha=0.2, label="80 % band")
        for m, seg in d.groupby("model", sort=False):
            ax.plot(pd.to_datetime(seg.window_start), seg.area_mean_display, "o-", ms=3, color=colors.get(m, "k"),
                    label=f"forecast: {m}")
        ax.axvline(pd.Timestamp("2024-12-30"), color="k", ls="--", lw=0.8)
        ax.set_ylabel(f"{gas} ({unit(p.cfg, gas)})")
        ax.legend(fontsize=7, ncol=3, loc="upper left")
    axes[0].set_title("Jan–May 2025 outlook: area mean with 80 % band (issued from the 25 Dec 2024 window)")
    save(fig, out / "fig_outlook_fan_area.png")

    fig, axes = plt.subplots(len(hot.groupby(["gas", "hub"])), 1, figsize=(11, 6), sharex=True)
    for ax, ((gas, hub), d) in zip(np.atleast_1d(axes), hot.groupby(["gas", "hub"], sort=False)):
        i, j = int(d.pixel_row.iloc[0]), int(d.pixel_col.iloc[0])
        g = p.gases.index(gas)
        sel = p.dates >= "2023-07-01"
        ax.plot(p.dates[sel], scaled(p.cfg, gas, p.values[sel, g, i, j]), "-", color="k", lw=0.8, label="observed")
        x = pd.to_datetime(d.window_start)
        ax.fill_between(x, d.p10_display, d.p90_display, color=GAS_COLORS[gas], alpha=0.2, label="80 % band")
        ax.plot(x, d.mean_display, "o-", ms=3, color=GAS_COLORS[gas], label="forecast")
        ax.axvline(pd.Timestamp("2024-12-30"), color="k", ls="--", lw=0.8)
        ax.set_title(f"{gas} at the {hub} pixel (models by lead: {', '.join(dict.fromkeys(d.model))})", fontsize=9)
        ax.set_ylabel(unit(p.cfg, gas))
        ax.legend(fontsize=7)
    save(fig, out / "fig_outlook_hotspots.png")


def _maps(p, fc_vals, out, tdates):
    steps = [1, 6, 12, 18, 24, 30]
    fig, axes = plt.subplots(len(p.gases), len(steps), figsize=(17, 8))
    for g, gas in enumerate(p.gases):
        f = np.where(p.cube.footprint, fc_vals[:, g], np.nan)
        vmin, vmax = np.nanpercentile(f[[s - 1 for s in steps]], [2, 98])
        models = lead_models(p, gas)
        for k, s in enumerate(steps):
            ax = axes[g, k]
            im = ax.imshow(scaled(p.cfg, gas, f[s - 1])[:12, 1:], cmap="YlOrRd", vmin=scaled(p.cfg, gas, vmin),
                           vmax=scaled(p.cfg, gas, vmax))
            ax.set_title(f"{gas} {tdates[s - 1].date()} (+{5 * s} d)\n[{models[s - 1]}]", fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
        fig.colorbar(im, ax=axes[g], fraction=0.015).set_label(unit(p.cfg, gas))
    fig.suptitle("Jan–May 2025 outlook: forecast mean maps (12×12 grid, north up)")
    save(fig, out / "fig_outlook_maps.png")


def _readme(p, tdates) -> str:
    lines = ["# Jan–May 2025 outlook (F5)", "",
             f"Issued from the last observed window (25 Dec 2024). 30 steps of 5 days: {tdates[0].date()} to "
             f"{tdates[-1].date()}. Everything was refitted on 2019–2024, with the climatology recomputed on those years.",
             "", "Model used per lead (frozen before the 2024 test, never changed afterwards):", ""]
    for gas in p.gases:
        segs = [f"leads {a}-{b} ({5 * a}-{5 * b} d): {c}" for a, b, c in p.fcfg["final"]["frozen"][gas]]
        lines.append(f"- **{gas}**: " + "; ".join(segs))
    lines += ["", "Files:",
              "- `outlook_<GAS>_2025.tif`: bands 1–30 are the mean, 31–60 the 10 % quantile, 61–90 the 90 % quantile, in mol/m².",
              "  Each band's description and tags give the lead, the window start and the model used.",
              "- `outlook_area_mean.csv`: area-mean forecast with an 80 % band, one row per gas and lead.",
              "- `outlook_hotspots.csv`: forecast at the Manali (NO2) and Ennore (SO2) hub pixels.",
              "- `fig_outlook_fan_area.png`, `fig_outlook_hotspots.png`, `fig_outlook_maps.png`.",
              "", "Interpretation: beyond about 30 days the forecast is the seasonal climatology (Forecasting/REPORT.md, F4).",
              "The bands are 80 % intervals from leave-one-year-out residuals. Pixel bands are clipped at 0."]
    return "\n".join(lines) + "\n"
