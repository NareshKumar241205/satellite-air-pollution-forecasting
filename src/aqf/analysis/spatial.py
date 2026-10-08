"""Step 4: spatial patterns. Long-term and seasonal mean maps, persistent hotspots, and the industrial hubs."""

from __future__ import annotations

import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from aqf.data import S5PCube, nearest_pixel
from aqf.plotting import gas_map, save, scaled, unit
from aqf.series import analysis_values, season_of

log = logging.getLogger(__name__)


def run(cube: S5PCube, cfg: dict) -> pd.DataFrame:
    out = cfg["paths"]["outputs"] / "04_spatial"
    out.mkdir(parents=True, exist_ok=True)
    v = analysis_values(cube, cfg)
    gases = cube.gases
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = np.nanmean(v, 0)                      # (G, H, W)
        cv = np.nanstd(v, 0) / mean

    fig, axes = plt.subplots(2, len(gases), figsize=(14, 7.6))
    for g, gas in enumerate(gases):
        gas_map(axes[0, g], cube, scaled(cfg, gas, mean[g]), cfg, f"{gas}: mean 2019–2024", cmap="YlOrRd",
                label=unit(cfg, gas))
        gas_map(axes[1, g], cube, cv[g], cfg, f"{gas}: coefficient of variation", cmap="PuBu", label="std / mean")
    fig.tight_layout()
    save(fig, out / "fig13_mean_maps.png")
    np.save(out / "mean_maps.npy", mean)

    # Seasonal mean maps (IMD seasons), shared colour scale per gas
    seasons = list(cfg["seasons"])
    s_of = season_of(cube.dates.month, cfg)
    fig, axes = plt.subplots(len(gases), len(seasons), figsize=(15, 3.4 * len(gases)))
    for g, gas in enumerate(gases):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            maps = [np.nanmean(v[s_of == s, g], 0) for s in seasons]
        lo, hi = np.nanpercentile(np.stack(maps), [2, 98])
        for k, s in enumerate(seasons):
            gas_map(axes[g, k], cube, scaled(cfg, gas, maps[k]), cfg, f"{gas} — {s}", cmap="YlOrRd",
                    vmin=scaled(cfg, gas, lo), vmax=scaled(cfg, gas, hi), hubs=(k == 0), label=unit(cfg, gas))
    fig.tight_layout()
    save(fig, out / "fig14_seasonal_maps.png")

    # Hotspot persistence: share of usable windows in which a pixel is in that window's top 10%
    fig, axes = plt.subplots(1, len(gases), figsize=(14, 4))
    for g, gas in enumerate(gases):
        x = v[:, g]
        ok = np.isfinite(x).any(axis=(1, 2))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            thr = np.nanpercentile(x[ok], 90, axis=(1, 2))
        hot = (x[ok] >= thr[:, None, None]).sum(0) / ok.sum()
        gas_map(axes[g], cube, hot, cfg, f"{gas}: share of windows in top-10% pixels", cmap="inferno", vmin=0,
                vmax=max(0.5, float(np.nanmax(hot))))
    fig.tight_layout()
    save(fig, out / "fig15_hotspot_persistence.png")

    # Industrial hubs: long-term enhancement over the area median, and monthly series at each hub
    rows = []
    fig, axes = plt.subplots(len(gases), 1, figsize=(11, 7.5), sharex=True)
    ym = pd.PeriodIndex(cube.dates, freq="M")
    for g, gas in enumerate(gases):
        area_med = np.nanmedian(mean[g][cube.footprint])
        with warnings.catch_warnings():  # unusable windows are all-NaN by design
            warnings.simplefilter("ignore", RuntimeWarning)
            dom = pd.Series(np.nanmean(v[:, g].reshape(len(v), -1), 1), index=cube.dates)
        axes[g].plot(dom.groupby(ym).mean().index.to_timestamp(), scaled(cfg, gas, dom.groupby(ym).mean()),
                     "k--", lw=1.2, label="area mean")
        for name, (la, lo) in cfg["hubs"].items():
            i, j = nearest_pixel(cube, la, lo)
            s = pd.Series(v[:, g, i, j], index=cube.dates)
            ms = s.groupby(ym).mean()
            axes[g].plot(ms.index.to_timestamp(), scaled(cfg, gas, ms), lw=0.9, label=name)
            rows.append({"gas": gas, "hub": name, "pixel_row": i, "pixel_col": j, "in_aoi": bool(cube.footprint[i, j]),
                         "mean": scaled(cfg, gas, mean[g, i, j]), "unit": unit(cfg, gas),
                         "ratio_to_area_median": mean[g, i, j] / area_med,
                         "pixel_rank_of_" + str(int(cube.footprint.sum())): int((mean[g][cube.footprint] > mean[g, i, j]).sum() + 1)})
        axes[g].set_ylabel(f"{gas} ({unit(cfg, gas)})")
    axes[0].legend(fontsize=7, ncol=3)
    axes[0].set_title("Monthly mean at the pixel containing each industrial hub")
    save(fig, out / "fig16_hub_series.png")
    hubs = pd.DataFrame(rows)
    hubs.round(3).to_csv(out / "hub_enhancement.csv", index=False)
    log.info("hubs\n%s", hubs.round(2).to_string())
    return hubs
