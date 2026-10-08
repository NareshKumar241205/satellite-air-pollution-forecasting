"""Step 3: temporal behaviour and trends of the area-averaged gas columns.

All series are area means over the area of interest, computed only from usable windows (see
`series.window_table`). This step is descriptive and covers all six years. Nothing here is fitted for
later modelling, so the 2024 test year doesn't leak into a model.
"""

from __future__ import annotations

import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.tsa.seasonal import STL

from aqf.data import S5PCube
from aqf.plotting import GAS_COLORS, gas_map, save, scaled, unit
from aqf.series import analysis_values, monthly, season_of, window_table
from aqf.stats import fdr_bh, seasonal_trend

log = logging.getLogger(__name__)
MONTHS = [pd.Timestamp(2000, k, 1).strftime("%b") for k in range(1, 13)]
WINDOWS_PER_YEAR = 73


def run(cube: S5PCube, cfg: dict) -> dict:
    out = cfg["paths"]["outputs"] / "03_trends"
    out.mkdir(parents=True, exist_ok=True)
    wt = window_table(cube, cfg)
    mon = monthly(wt)
    mon.to_csv(out / "monthly_domain_mean.csv", index=False)
    gases = cube.gases

    _time_series(wt, cfg, gases, out)
    _seasonal_cycle(mon, cfg, gases, out)
    annual = _annual(wt, cfg, gases, out)
    _stl(wt, cfg, gases, out)
    tests = _trend_tests(mon, cfg, gases, out)
    _pixel_trends(cube, wt, cfg, out)
    lock = _lockdown(wt, cfg, gases, out)
    _gas_correlation(mon, gases, out)
    return {"annual": annual, "trend_tests": tests, "lockdown": lock}


def _time_series(wt, cfg, gases, out):
    fig, axes = plt.subplots(len(gases), 1, figsize=(11, 7.5), sharex=True)
    for ax, gas in zip(axes, gases):
        d = wt[(wt.gas == gas)].set_index("date")
        y = scaled(cfg, gas, d["mean"].where(d.usable))
        ax.plot(y.index, y, lw=0.6, color=GAS_COLORS[gas], alpha=0.45, label="5-day area mean")
        roll = y.rolling(12, center=True, min_periods=6).mean()  # ~2 months
        ax.plot(roll.index, roll, lw=1.8, color=GAS_COLORS[gas], label="2-month running mean")
        ax.axhline(y.mean(), color="k", lw=0.6, ls=":")
        for yr in cfg["split"]["val_years"] + cfg["split"]["test_years"]:
            ax.axvspan(pd.Timestamp(yr, 1, 1), pd.Timestamp(yr, 12, 31), color="#d5d8dc", alpha=0.35, lw=0)
        ax.set_ylabel(f"{gas} ({unit(cfg, gas)})")
        ax.legend(loc="upper left", fontsize=7, ncol=2)
    axes[0].set_title("Area-mean column density over Chennai, 2019–2024 (shaded: validation 2023 and test 2024 years)")
    save(fig, out / "fig06_time_series.png")


def _seasonal_cycle(mon, cfg, gases, out):
    fig, axes = plt.subplots(1, len(gases), figsize=(14, 3.6))
    years = sorted(mon.year.unique())
    cmap = plt.get_cmap("viridis", len(years))
    for ax, gas in zip(axes, gases):
        d = mon[mon.gas == gas]
        data = [scaled(cfg, gas, d[d.month == k].value.values) for k in range(1, 13)]
        ax.boxplot(data, tick_labels=MONTHS, widths=0.6, showfliers=False,
                   medianprops={"color": GAS_COLORS[gas]})
        for i, yr in enumerate(years):
            dy = d[d.year == yr]
            ax.plot(dy.month, scaled(cfg, gas, dy.value), "-", lw=0.8, alpha=0.8, color=cmap(i), label=str(yr))
        ax.set_title(f"{gas} seasonal cycle (monthly area means)")
        ax.set_ylabel(unit(cfg, gas))
        ax.tick_params(axis="x", labelsize=7)
    axes[-1].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    save(fig, out / "fig07_seasonal_cycle.png")

    clim = mon.groupby(["gas", "month"]).value.agg(["mean", "std", "min", "max"]).reset_index()
    clim.to_csv(out / "monthly_climatology.csv", index=False)
    mon["season"] = season_of(mon.month, cfg)
    mon.groupby(["gas", "season"]).value.mean().unstack(0).to_csv(out / "season_means.csv")

    fig, axes = plt.subplots(1, len(gases), figsize=(14, 3.2))
    for ax, gas in zip(axes, gases):
        piv = mon[mon.gas == gas].pivot(index="year", columns="month", values="value")
        im = ax.imshow(scaled(cfg, gas, piv.values), cmap="YlOrRd", aspect="auto")
        ax.set_xticks(range(12), [m[0] for m in MONTHS])
        ax.set_yticks(range(len(piv)), piv.index)
        ax.grid(False)
        ax.set_title(f"{gas} monthly mean by year")
        plt.colorbar(im, ax=ax, fraction=0.04).set_label(unit(cfg, gas), fontsize=7)
    fig.tight_layout()
    save(fig, out / "fig08_year_month_heatmap.png")


def _annual(wt, cfg, gases, out):
    u = wt[wt.usable]
    ann = u.groupby(["gas", "year"])["mean"].agg(["mean", "std", "size"]).reset_index()
    base = ann[ann.year == ann.year.min()].set_index("gas")["mean"]
    ann["pct_vs_first_year"] = 100 * (ann["mean"] / ann.gas.map(base) - 1)
    ann["pct_vs_prev_year"] = 100 * ann.groupby("gas")["mean"].pct_change()
    ann.to_csv(out / "annual_means.csv", index=False)

    fig, axes = plt.subplots(1, len(gases), figsize=(13, 3.2))
    for ax, gas in zip(axes, gases):
        d = ann[ann.gas == gas]
        ax.bar(d.year, scaled(cfg, gas, d["mean"]), color=GAS_COLORS[gas], alpha=0.85,
               yerr=scaled(cfg, gas, d["std"] / np.sqrt(d["size"])), capsize=3)
        for x, y, p in zip(d.year, scaled(cfg, gas, d["mean"]), d.pct_vs_first_year):
            ax.annotate(f"{p:+.1f}%", (x, y), xytext=(0, 4), textcoords="offset points", ha="center", fontsize=7)
        lo, hi = scaled(cfg, gas, d["mean"]).min(), scaled(cfg, gas, d["mean"]).max()
        ax.set_ylim(lo - 0.6 * (hi - lo), hi + 0.5 * (hi - lo))
        ax.set_title(f"{gas} annual mean (± s.e.; % vs {int(d.year.min())})")
        ax.set_ylabel(unit(cfg, gas))
    fig.tight_layout()
    save(fig, out / "fig09_annual_means.png")
    return ann


def _stl(wt, cfg, gases, out):
    """STL on the 5-day series: trend + seasonal (period one year = 73 windows) + residual.

    Unusable windows are linearly interpolated for STL only, which needs an unbroken series.
    """
    fig, axes = plt.subplots(4, len(gases), figsize=(14, 8), sharex=True)
    rows = []
    for j, gas in enumerate(gases):
        d = wt[wt.gas == gas].set_index("date")
        y = scaled(cfg, gas, d["mean"].where(d.usable)).interpolate("time", limit_direction="both")
        res = STL(y, period=WINDOWS_PER_YEAR, robust=True).fit()
        for i, (name, comp) in enumerate([("observed", y), ("trend", res.trend), ("seasonal", res.seasonal),
                                          ("residual", res.resid)]):
            ax = axes[i, j]
            ax.plot(comp.index, comp, lw=0.8 if name != "trend" else 1.6, color=GAS_COLORS[gas])
            if i == 0:
                ax.set_title(f"{gas} ({unit(cfg, gas)})")
            if j == 0:
                ax.set_ylabel(name)
        var = np.var(y)
        rows.append({"gas": gas, "interpolated_windows": int((~d.usable).sum()),
                     "seasonal_strength": max(0.0, 1 - np.var(res.resid) / np.var(res.seasonal + res.resid)),
                     "trend_strength": max(0.0, 1 - np.var(res.resid) / np.var(res.trend + res.resid)),
                     "share_var_seasonal": np.var(res.seasonal) / var, "share_var_trend": np.var(res.trend) / var,
                     "share_var_residual": np.var(res.resid) / var})
    fig.suptitle("STL decomposition of the 5-day area-mean series (period = 1 year)")
    fig.tight_layout()
    save(fig, out / "fig10_stl.png")
    pd.DataFrame(rows).round(3).to_csv(out / "stl_strength.csv", index=False)


def _trend_tests(mon, cfg, gases, out):
    """Seasonal Mann-Kendall + Sen slope on monthly area means: all years, and without 2020 (COVID)."""
    rows = []
    for gas in gases:
        d = mon[mon.gas == gas]
        mean = d.value.mean()
        for label, sel in (("2019-2024", d), ("2019-2024 excl. 2020", d[d.year != 2020]),
                           ("train years 2019-2022", d[d.year.isin(cfg["split"]["train_years"])])):
            r = seasonal_trend(sel.value.values, sel.year.values, sel.month.values)
            rows.append({"gas": gas, "period": label, "n_months": len(sel),
                         f"sen_slope_per_year": scaled(cfg, gas, r.slope), "unit": unit(cfg, gas) + "/yr",
                         "pct_per_year": 100 * r.slope / mean, "mk_z": r.z, "p_value": r.p,
                         "significant_5pct": r.p < 0.05, "direction": r.direction})
    t = pd.DataFrame(rows)
    t.round(4).to_csv(out / "trend_tests.csv", index=False)
    log.info("trend tests\n%s", t[["gas", "period", "pct_per_year", "p_value"]].round(4).to_string())
    return t


def _pixel_trends(cube, wt, cfg, out):
    """Per-pixel seasonal Sen slope (% of the pixel mean per year) on monthly means; dots = FDR-significant."""
    vals = analysis_values(cube, cfg)
    ym = pd.PeriodIndex(cube.dates, freq="M")
    keys = sorted(set(ym))
    with warnings.catch_warnings():  # pixels with no valid value in a month give NaN, which is intended
        warnings.simplefilter("ignore", RuntimeWarning)
        mon = np.stack([np.nanmean(vals[ym == k], axis=0) for k in keys])  # (M, G, H, W)
    years = np.array([k.year for k in keys])
    months = np.array([k.month for k in keys])
    G, H, W = mon.shape[1:]
    fp = cube.footprint
    fig, axes = plt.subplots(1, G, figsize=(14, 4))
    rows = []
    for g, gas in enumerate(cube.gases):
        pct = np.full((H, W), np.nan)
        p = np.full((H, W), np.nan)
        for i, j in zip(*np.nonzero(fp)):
            x = mon[:, g, i, j]
            if np.isfinite(x).sum() < 24:
                continue
            r = seasonal_trend(x, years, months)
            pct[i, j] = 100 * r.slope / np.nanmean(x)
            p[i, j] = r.p
        sig = fdr_bh(p)
        lim = np.nanpercentile(np.abs(pct), 98)
        gas_map(axes[g], cube, pct, cfg, f"{gas}: trend, % of pixel mean per year", cmap="RdBu_r",
                vmin=-lim, vmax=lim, label="%/yr")
        ii, jj = np.nonzero(sig)
        axes[g].plot(cube.lon[jj], cube.lat[ii], "k.", ms=3)
        rows.append({"gas": gas, "pixels": int(np.isfinite(p).sum()), "sig_fdr05": int(sig.sum()),
                     "sig_increasing": int((sig & (pct > 0)).sum()), "sig_decreasing": int((sig & (pct < 0)).sum()),
                     "median_pct_per_year": float(np.nanmedian(pct))})
    fig.suptitle("Per-pixel seasonal Sen slope 2019–2024 (black dots: significant after Benjamini–Hochberg FDR 5%)")
    fig.tight_layout()
    save(fig, out / "fig11_pixel_trend_maps.png")
    pd.DataFrame(rows).round(3).to_csv(out / "pixel_trend_summary.csv", index=False)


def _lockdown(wt, cfg, gases, out):
    """COVID-19 lockdown natural experiment: same calendar window in every year."""
    s, e = pd.Timestamp(cfg["lockdown"]["start"]), pd.Timestamp(cfg["lockdown"]["end"])
    u = wt[wt.usable].copy()
    doy = u.date.dt.dayofyear
    u = u[(doy >= s.dayofyear) & (doy <= e.dayofyear)]
    tab = u.groupby(["gas", "year"])["mean"].mean().unstack("year")
    others = [y for y in tab.columns if y != s.year]
    tab["mean_other_years"] = tab[others].mean(1)
    tab["pct_2020_vs_others"] = 100 * (tab[s.year] / tab["mean_other_years"] - 1)
    tab.to_csv(out / "lockdown_comparison.csv")

    fig, axes = plt.subplots(1, len(gases), figsize=(13, 3.2))
    for ax, gas in zip(axes, gases):
        r = tab.loc[gas]
        yrs = [c for c in tab.columns if isinstance(c, (int, np.integer))]
        ax.bar(yrs, scaled(cfg, gas, r[yrs].astype(float)),
               color=[GAS_COLORS[gas] if y == s.year else "#bdc3c7" for y in yrs])
        ax.set_title(f"{gas}, {s:%d %b}–{e:%d %b}: 2020 is {r.pct_2020_vs_others:+.1f}% vs other years", fontsize=8)
        ax.set_ylabel(unit(cfg, gas))
        vals = scaled(cfg, gas, r[yrs].astype(float))
        ax.set_ylim(vals.min() * 0.8, vals.max() * 1.08)
    fig.tight_layout()
    save(fig, out / "fig12_lockdown.png")
    return tab


def _gas_correlation(mon, gases, out):
    """Correlation between gases of the monthly anomalies (seasonal cycle removed)."""
    m = mon.copy()
    m["anom"] = m.value - m.groupby(["gas", "month"]).value.transform("mean")
    piv = m.pivot_table(index=["year", "month"], columns="gas", values="anom")[gases]
    raw = m.pivot_table(index=["year", "month"], columns="gas", values="value")[gases]
    pd.concat({"raw_monthly": raw.corr(method="spearman"), "anomaly": piv.corr(method="spearman")}).round(3) \
        .to_csv(out / "gas_correlation_spearman.csv")
