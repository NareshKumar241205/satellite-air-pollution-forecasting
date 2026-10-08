"""Step 1: dataset inventory and integrity.

Checks the manifest against the files on disk, checks that every raster of a sensor shares one grid,
records SHA-256 provenance, and plots data availability over time.
"""

from __future__ import annotations

import json
import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

from aqf.data import load_manifest, sha256
from aqf_analysis.plotting import save

log = logging.getLogger(__name__)


def run(cfg: dict, hash_files: bool = True) -> dict:
    out = cfg["paths"]["outputs"] / "01_inventory"
    m = load_manifest(cfg)
    report: dict = {"n_windows": len(m), "first": str(m.start_date.min().date()),
                    "last": str(m.start_date.max().date()),
                    "step_days": m.start_date.diff().dt.days.dropna().value_counts().to_dict(),
                    "duplicate_window_ids": int(m.window_id.duplicated().sum()),
                    "windows_per_year": m.groupby("year").size().to_dict(),
                    "status": m.status.value_counts().to_dict()}

    rows = []
    for sensor, col, d in (("s5p", "s5p_file", cfg["paths"]["s5p_dir"]), ("s2", "s2_file", cfg["paths"]["s2_dir"])):
        listed = set(m[col] if sensor == "s5p" else m.loc[m.has_s2, col])
        on_disk = {p.name for p in d.glob("*.tif")}
        report[f"{sensor}_missing_files"] = sorted(listed - on_disk)
        report[f"{sensor}_unlisted_files"] = sorted(on_disk - set(m[col]))
        grids = set()
        for name in sorted(listed & on_disk):
            p = d / name
            with rasterio.open(p) as src:
                grids.add((src.count, src.shape, str(src.crs), tuple(round(v, 12) for v in list(src.transform)[:6])))
                rows.append({"sensor": sensor, "file": name, "bytes": p.stat().st_size, "bands": src.count,
                             "height": src.height, "width": src.width, "crs": str(src.crs), "nodata": src.nodata,
                             "sha256": sha256(p) if hash_files else ""})
        report[f"{sensor}_distinct_grids"] = len(grids)
        report[f"{sensor}_grid"] = [list(map(str, g)) for g in grids]
    prov = pd.DataFrame(rows)
    out.mkdir(parents=True, exist_ok=True)
    prov.to_csv(out / "provenance.csv", index=False)
    report["dataset_bytes"] = int(prov.bytes.sum())

    # S2 availability: complete monsoon dropout is a key constraint for using S2 at all
    avail = m.pivot_table(index="year", columns="month", values="has_s2", aggfunc="mean")
    avail.round(3).to_csv(out / "s2_availability_year_month.csv")
    report["s2_availability_by_month"] = m.groupby("month").has_s2.mean().round(3).to_dict()

    fig, axes = plt.subplots(2, 1, figsize=(11, 4.2), gridspec_kw={"height_ratios": [1, 2.2]})
    ax = axes[0]
    ax.bar(m.start_date, np.ones(len(m)), width=5, color=np.where(m.has_s2, "#2e86c1", "#e5e7e9"))
    ax.set_yticks([])
    ax.set_title("Sentinel-2 availability per 5-day window (blue = image, grey = cloud gap); S5P is present in every window")
    ax.grid(False)
    ax = axes[1]
    im = ax.imshow(avail.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(12), [pd.Timestamp(2000, k, 1).strftime("%b") for k in avail.columns])
    ax.set_yticks(range(len(avail)), avail.index)
    for (i, j), v in np.ndenumerate(avail.values):
        ax.text(j, i, f"{v:.0%}", ha="center", va="center", fontsize=7, color="w" if v > 0.6 else "k")
    ax.grid(False)
    ax.set_title("Share of windows with a Sentinel-2 image, by year and month")
    plt.colorbar(im, ax=ax, fraction=0.02)
    fig.tight_layout()
    save(fig, out / "fig01_availability.png")

    (out / "inventory.json").write_text(json.dumps(report, indent=2, default=str))
    log.info("inventory: %s", {k: report[k] for k in ("n_windows", "s5p_distinct_grids", "s2_distinct_grids")})
    return report
