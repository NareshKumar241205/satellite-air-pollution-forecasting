"""Step 2: Sentinel-5P data quality.

How much of each gas map is valid, how often values are exactly zero (retrieval noise floored at
zero), which windows are unusable, and where in the area the problems sit.
"""

from __future__ import annotations

import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from aqf.data import S5PCube
from aqf.plotting import GAS_COLORS, gas_map, save, scaled, unit
from aqf.series import window_table

log = logging.getLogger(__name__)


def run(cube: S5PCube, cfg: dict) -> pd.DataFrame:
    out = cfg["paths"]["outputs"] / "02_quality"
    out.mkdir(parents=True, exist_ok=True)
    wt = window_table(cube, cfg)
    wt.to_csv(out / "window_quality.csv", index=False)

    summ = wt.groupby("gas", sort=False).agg(
        windows=("usable", "size"), usable=("usable", "sum"), all_zero_maps=("all_zero", "sum"),
        low_coverage=("coverage", lambda c: int((c < cfg["quality"]["min_coverage"]).sum())),
        mean_coverage=("coverage", "mean"), mean_zero_frac=("zero_frac", "mean"), outlier_windows=("outlier", "sum"))
    summ["usable_pct"] = 100 * summ.usable / summ.windows
    fp = cube.footprint
    summ["grid"] = f"{cube.values.shape[2]}x{cube.values.shape[3]} ({fp.sum()} pixels in the area of interest)"
    summ.round(4).to_csv(out / "quality_summary.csv")
    by_year = wt.groupby(["gas", "year"]).agg(usable_pct=("usable", "mean"), zero_frac=("zero_frac", "mean"))
    (by_year * [100, 1]).round(3).unstack(0).to_csv(out / "quality_by_year.csv")
    log.info("quality summary\n%s", summ.round(3).to_string())

    # Timeline: coverage and zero share per window
    fig, axes = plt.subplots(len(cube.gases), 1, figsize=(11, 6), sharex=True)
    for ax, gas in zip(axes, cube.gases):
        d = wt[wt.gas == gas]
        ax.plot(d.date, d.coverage, lw=0.8, color="#7f8c8d", label="valid share")
        ax.plot(d.date, d.zero_frac, lw=0.8, color=GAS_COLORS[gas], label="zero share (of valid)")
        bad = d[~d.usable]
        ax.scatter(bad.date, np.full(len(bad), 1.05), marker="|", color="k", s=40, label="unusable window")
        ax.set_ylim(-0.02, 1.12)
        ax.set_ylabel(gas)
        ax.legend(loc="upper left", fontsize=7, ncol=3)
    axes[0].set_title("Sentinel-5P map quality per 5-day window")
    save(fig, out / "fig02_quality_timeline.png")

    # Zero share by calendar month (is the zero floor seasonal?)
    fig, ax = plt.subplots(figsize=(7, 3))
    for gas in cube.gases:
        d = wt[wt.gas == gas].groupby("month").zero_frac.mean()
        ax.plot(d.index, d.values, "o-", color=GAS_COLORS[gas], label=gas)
    ax.set_xticks(range(1, 13), [pd.Timestamp(2000, k, 1).strftime("%b") for k in range(1, 13)])
    ax.set_ylabel("mean share of zero pixels")
    ax.set_title("Zero-floored pixels by month")
    ax.legend()
    save(fig, out / "fig03_zero_share_by_month.png")

    # Where are the zeros? Per-pixel share of windows that are exactly zero
    fig, axes = plt.subplots(1, len(cube.gases), figsize=(13, 3.8))
    for g, (ax, gas) in enumerate(zip(axes, cube.gases)):
        x = cube.values[:, g]
        with np.errstate(invalid="ignore"):
            zf = (x == 0).sum(0) / np.isfinite(x).sum(0)
        gas_map(ax, cube, zf, cfg, f"{gas}: share of windows at exactly 0", cmap="magma", vmin=0, vmax=1)
    fig.tight_layout()
    save(fig, out / "fig04_zero_share_maps.png")

    # Value distributions of all valid pixels (log x; zeros shown separately in the title)
    fig, axes = plt.subplots(1, len(cube.gases), figsize=(13, 3))
    for g, (ax, gas) in enumerate(zip(axes, cube.gases)):
        v = cube.values[:, g][:, fp].ravel()
        v = v[np.isfinite(v)]
        pos = scaled(cfg, gas, v[v > 0])
        ax.hist(pos, bins=np.logspace(np.log10(pos.min()), np.log10(pos.max()), 60), color=GAS_COLORS[gas])
        ax.set_xscale("log")
        ax.set_xlabel(unit(cfg, gas))
        ax.set_title(f"{gas}: {np.mean(v == 0):.1%} of valid pixels are 0, {np.mean(v < 0):.1%} negative", fontsize=8)
    fig.tight_layout()
    save(fig, out / "fig05_value_histograms.png")
    return wt
