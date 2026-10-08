"""Step 5: Sentinel-2 summary, streamed one file at a time (the full stack would not fit in RAM).

Covers cloud-free coverage per image, area-median spectral indices over time (land-surface change),
and long-term index maps for context against the pollution hotspots.
"""

from __future__ import annotations

import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

from aqf.data import load_manifest, read_raster
from aqf.plotting import save
from aqf.series import season_of
from aqf.stats import seasonal_trend

log = logging.getLogger(__name__)
INDICES = ["NDVI", "NDBI", "NDMI"]


def run(cfg: dict) -> pd.DataFrame:
    out = cfg["paths"]["outputs"] / "05_sentinel2"
    out.mkdir(parents=True, exist_ok=True)
    m = load_manifest(cfg)
    m = m[m.has_s2].reset_index(drop=True)
    names = cfg["s2_bands"]
    idx = [names.index(b) + 1 for b in INDICES]
    refl = [names.index(b) + 1 for b in ("B4", "B8", "B11", "B12")]

    with rasterio.open(cfg["paths"]["s2_dir"] / m.s2_file.iloc[0]) as src:
        H, W = src.shape
        b = src.bounds
    sums = np.zeros((len(INDICES), H, W))
    counts = np.zeros((H, W))
    ever_valid = np.zeros((H, W), bool)
    rows = []
    for _, r in m.iterrows():
        a = read_raster(cfg["paths"]["s2_dir"] / r.s2_file, idx + refl)  # (7, H, W)
        valid = np.isfinite(a).all(0)
        ever_valid |= valid
        row = {"date": r.start_date, "valid_share": valid.mean()}
        for k, name in enumerate(INDICES + ["B4", "B8", "B11", "B12"]):
            row[name] = float(np.median(a[k][valid])) if valid.any() else np.nan
        rows.append(row)
        if r.split == "train":  # long-term maps from train years only, so they can be reused as model inputs later
            sums += np.where(valid, a[: len(INDICES)], 0)
            counts += valid
    ts = pd.DataFrame(rows)
    aoi_share = ever_valid.mean()
    ts["valid_share_of_aoi"] = ts.valid_share / aoi_share
    ts.to_csv(out / "s2_window_stats.csv", index=False)

    fig, axes = plt.subplots(4, 1, figsize=(11, 8), sharex=True)
    axes[0].bar(ts.date, ts.valid_share_of_aoi, width=5, color="#2e86c1")
    axes[0].set_ylabel("cloud-free share")
    axes[0].set_title(f"Sentinel-2 per available window ({len(ts)} of {len(load_manifest(cfg))}); "
                      f"{1 - aoi_share:.1%} of the raster is never valid (outside the area)")
    for ax, name, c in zip(axes[1:], INDICES, ["#27ae60", "#7f8c8d", "#2980b9"]):
        ok = ts.valid_share_of_aoi >= 0.5
        ax.plot(ts.date[ok], ts[name][ok], ".", ms=3, color=c, alpha=0.6)
        ax.plot(ts.date[ok], ts[name][ok].rolling(9, center=True, min_periods=4).median(), color=c, lw=1.6)
        ax.set_ylabel(f"area median {name}")
    save(fig, out / "fig17_s2_time_series.png")

    # Trends in land surface (monthly medians, seasonal MK); only windows with >=50% cloud-free area
    ok = ts[ts.valid_share_of_aoi >= 0.5].copy()
    ok["year"], ok["month"] = ok.date.dt.year, ok.date.dt.month
    mon = ok.groupby(["year", "month"])[INDICES].median().reset_index()
    trows = []
    for name in INDICES:
        r = seasonal_trend(mon[name].values, mon.year.values, mon.month.values)
        trows.append({"index": name, "n_months": int(mon[name].notna().sum()), "sen_slope_per_year": r.slope,
                      "mk_z": r.z, "p_value": r.p})
    pd.DataFrame(trows).round(5).to_csv(out / "s2_index_trends.csv", index=False)
    ok["season"] = season_of(ok.month, cfg)
    ok.groupby("season")[INDICES + ["valid_share_of_aoi"]].median().round(4).to_csv(out / "s2_season_medians.csv")

    with np.errstate(invalid="ignore", divide="ignore"):
        maps = sums / counts
    extent = [b.left, b.right, b.bottom, b.top]
    fig, axes = plt.subplots(1, 4, figsize=(17, 4))
    for ax, k, name, cmap, lim in [(axes[0], 0, "NDVI", "RdYlGn", (-0.2, 0.8)), (axes[1], 1, "NDBI", "RdGy_r", (-0.5, 0.3)),
                                   (axes[2], 2, "NDMI", "BrBG", (-0.3, 0.5))]:
        im = ax.imshow(np.where(ever_valid, maps[k], np.nan), extent=extent, cmap=cmap, vmin=lim[0], vmax=lim[1])
        plt.colorbar(im, ax=ax, fraction=0.046)
        ax.set_title(f"Mean {name}, train years 2019–2022")
    im = axes[3].imshow(np.where(ever_valid, counts, np.nan), extent=extent, cmap="viridis")
    plt.colorbar(im, ax=axes[3], fraction=0.046)
    axes[3].set_title("Cloud-free observations per pixel (train years)")
    for ax in axes:
        ax.grid(False)
        for name, (la, lo) in cfg["hubs"].items():
            ax.plot(lo, la, "w^", ms=5, mec="k", mew=0.6)
    fig.tight_layout()
    save(fig, out / "fig18_s2_index_maps.png")
    log.info("s2: %d images, median cloud-free share %.2f", len(ts), ts.valid_share_of_aoi.median())
    return ts
