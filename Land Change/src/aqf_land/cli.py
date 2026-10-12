"""Command line: `uv run aqf-land <step>`. Results go to Land Change/results/<step>/."""

from __future__ import annotations

import argparse
import json
import logging

import numpy as np
import pandas as pd

from aqf.data import ROOT

log = logging.getLogger("aqf_land")


def cache(cfg: dict, force: bool) -> tuple[dict, dict]:
    from aqf_land.tnaq import read_cities

    return read_cities(cfg, "NO2", force), read_cities(cfg, "S2", force)


def provenance(cfg: dict, pairs: pd.DataFrame | None = None) -> None:
    from aqf_land import provenance as pv
    from aqf_land.tnaq import list_files

    out = cfg["outputs"] / "P1_provenance"
    out.mkdir(parents=True, exist_ok=True)
    files = list_files(cfg, "NO2")
    meta = pv.file_metadata(list(files.path))
    man = sorted((cfg["tnaq"]["dir"] / "TNAQ_MANIFEST").glob("manifest_s5p_no2_*.csv"))
    cols = list(pd.read_csv(man[0], nrows=1).columns) if man else []
    acq = pv.acquisition_info((ROOT / cfg["tnaq"]["acquisition_script"]).resolve())
    sv = pv.season_versions(files, cfg["season_months"])
    sv.to_csv(out / "season_windows_by_version.csv", index=False)
    lines = [
        "# TROPOMI NO2 processor versions in the TNAQ export", "",
        f"- Files checked: {len(files)} ({files.date.min().date()} .. {files.date.max().date()}).",
        f"- Collection: `{acq.get('collection', '?')}`, band `{acq.get('band', '?')}` (from `{acq['script']}`).",
        f"- Distinct metadata sets in the GeoTIFFs: {len(meta)}. Dataset tags: `{meta.dataset_tags.iloc[0]}`; "
        f"bands: `{meta.band_names.iloc[0]}`.",
        f"- Manifest columns: `{', '.join(cols)}`.",
        "",
        "**Result: the processor version is not recorded anywhere in the local export.** Neither the GeoTIFFs nor",
        "the manifests carry it. The exact per-window versions need `scripts/ee_s5p_versions.py` (Earth Engine).",
        "",
        "The collection is the operational OFFL stream, not a consistent reprocessing. By the processor timeline",
        "(`provenance.py`, months verified with Earth Engine below) the Jan-Apr season windows fall into:", "",
        sv.pivot(index="year", columns="version", values="windows").fillna(0).astype(int).to_markdown(), "",
        "Every Jan-Apr season of 2019-2024 uses a different processor version. 2019-2021 are v1.x and 2022-2024",
        "are v2.3.1 or later, so the **2021 -> 2022 pair",
        "straddles the v1 -> v2 change** (v2.2, 1 Jul 2021: new cloud and surface-albedo inputs that raise",
        "tropospheric NO2, more over polluted pixels). City-year fixed effects absorb a uniform step but not a",
        "level-dependent one, so G1 is also reported without that pair, and the main analysis must either use a",
        "consistently reprocessed product (S5P-PAL / RPRO v2.4) or model the step as level-dependent.",
    ]
    ee_csv = out / "ee_processor_versions.csv"
    if ee_csv.exists():
        v = pd.read_csv(ee_csv)
        lines += ["", "## Exact versions from Earth Engine (`scripts/ee_s5p_versions.py`)", "",
                  v.pivot_table(index="month", columns="value", values="images", aggfunc="sum").fillna(0)
                  .astype(int).to_markdown()]
    else:
        lines += ["", "_Exact versions: not yet queried (`scripts/ee_s5p_versions.py` needs an Earth Engine account)._"]
    if pairs is not None:
        lines += ["", "Mean same-season dNO2 (%, log change) per year pair and city, all available years. A retrieval",
                  "step would show up in every city in the same pair (2020 and 2021 also carry the lockdown and",
                  "its rebound):", "", pv.common_step(pairs).to_markdown()]
    (out / "PROVENANCE.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    log.info("wrote %s", out / "PROVENANCE.md")


def layers(cfg: dict, no2: dict, s2: dict) -> pd.DataFrame:
    from aqf_land.layers import CLASSES, land_fractions, land_from_csv, no2_seasonal, pair_table
    from aqf_land.tnaq import block_lat_lon

    out = cfg["outputs"] / "P1_layers"
    out.mkdir(parents=True, exist_ok=True)
    years = cfg["years"]["tuning"] + cfg["years"]["held_out_record"]   # all local years; G1 filters to tuning
    tables, fr, sources = [], {}, {}
    ee_dir = ROOT / "Land Change/data/ee_landfrac"
    ee_src = ee_dir / "source.json"
    ee_label = json.loads(ee_src.read_text())["label"] if ee_src.exists() else "Earth Engine table"
    for box in no2["boxes"]:
        c = box.city
        lat, lon = block_lat_lon(box, no2["transform"], cfg["grid"]["block"])
        land = land_from_csv(ee_dir, c, (len(lat), len(lon)), years, cfg["grid"]["min_valid_frac"])
        sources[c] = ee_label if land is not None else "provisional 1 km TNAQ S2"
        if land is None:
            land = land_fractions(s2[c], s2["dates"], years, cfg)
        n = no2_seasonal(no2[c], no2["dates"], years, cfg)
        fr.update({f"{c}_{y}": land[y] for y in years})
        fr.update({f"{c}_{y}_no2": n[y] for y in years})
        tables.append(pair_table(c, lat, lon, land, n, cfg["dL_class"]))
    pairs = pd.concat(tables, ignore_index=True)
    pairs["split"] = np.where(pairs.y1.isin(cfg["years"]["tuning"]), "tuning", "held_out_record")
    pairs.to_csv(out / "pixel_year_pairs.csv", index=False)
    np.savez_compressed(out / "land_fractions_no2_season.npz", classes=np.array(CLASSES), **fr)
    (out / "layers.json").write_text(json.dumps({
        "season_months": cfg["season_months"], "years": years, "grid_deg": 0.01 * cfg["grid"]["block"],
        "dL_class": cfg["dL_class"], "classes": cfg["classes"], "land_source": sources,
        "note": "The 1 km proxy moves a block fraction in 6.25 pp steps and tracks greenness (see G1)."}, indent=2))
    log.info("pairs: %d rows, %d with both dL and dNO2", len(pairs), pairs[["dL", "dNO2"]].notna().all(1).sum())
    return pairs


def g1(cfg: dict, pairs: pd.DataFrame) -> pd.DataFrame:
    from aqf_land.g1 import run_g1

    out = cfg["outputs"] / "G1_power"
    out.mkdir(parents=True, exist_ok=True)
    res = run_g1(pairs, cfg)
    res.to_csv(out / "g1_power.csv", index=False)
    _g1_figures(cfg, pairs, out)
    keep = ["n", "dL_sd_within", "dL_noise_sd", "dL_signal_sd", "corr_dL_dveg", "dNO2_resid_sd", "n_eff",
            "mde_pct_per_pp", "sd_required_pp", "passes"]
    pooled = res[res.city == "ALL"].set_index("pairs")[keep]
    passes = bool(pooled.passes.any())
    for name, r in res[res.city == "ALL_weighted"].set_index("pairs").iterrows():
        log.info("G1 %-18s city-weighted MDE %.2f %%/pp (sensitivity, not the gate)", name, r.mde_pct_per_pp)
    (out / "G1_RESULT.json").write_text(json.dumps({
        "passes": passes, "effect_pct_per_pp": cfg["g1"]["effect_pct_per_pp"],
        "tuning_years": cfg["years"]["tuning"], "cities": list(cfg["cities"]),
        "dL_source": json.loads((cfg["outputs"] / "P1_layers" / "layers.json").read_text())["land_source"],
        "pooled": json.loads(pooled.replace([np.inf], None).to_json(orient="index")),
        "sensitivity_city_weighted_mde": res[res.city == "ALL_weighted"].set_index("pairs").mde_pct_per_pp
        .replace([np.inf], None).to_dict(),
        "note": "passes uses the pooled MDE (the gate as defined); the city-weighted MDE is a sensitivity only"},
        indent=2))
    for name, r in pooled.iterrows():
        log.info("G1 %-18s signal sd %.2f pp (noise %.1f), MDE %s %%/pp, needs true sd >= %.2f pp -> %s", name,
                 r.dL_signal_sd, r.dL_noise_sd, f"{r.mde_pct_per_pp:.2f}", r.sd_required_pp,
                 "PASS" if r.passes else "FAIL")
    return res


def _g1_figures(cfg: dict, pairs: pd.DataFrame, out) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tune = cfg["years"]["tuning"]
    p = pairs[pairs.y0.isin(tune) & pairs.y1.isin(tune)].dropna(subset=["dL"])
    cities = list(cfg["cities"])
    nr = int(np.ceil(len(cities) / 3))
    fig, axes = plt.subplots(nr, 3, figsize=(11, 3 * nr), sharex=True)
    bins = np.arange(-40, 41, 6.25) - 3.125
    thr = cfg["g1"]["dL_meaningful_pp"]
    for ax, c in zip(axes.flat, cities):
        d = p[p.city == c].dL
        ax.hist(d, bins=bins, color="#4C78A8")
        for x in (-thr, thr):
            ax.axvline(x, color="0.4", lw=0.8, ls="--")
        ax.set_title(f"{c} (n={len(d)}, |dL|>{thr:g}: {(d.abs() > thr).mean():.0%})", fontsize=9)
        ax.set_yscale("log")
    for ax in axes[-1]:
        ax.set_xlabel("dL built-up (pp per year)")
    fig.suptitle(f"G1: year-to-year built-up change per S5P pixel, tuning pairs {tune[0]}-{tune[-1]} "
                 "(provisional 1 km S2)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "fig_g1_dL_hist.png", dpi=130)
    plt.close(fig)

    # cumulative change over the tuning period, per city map
    fig, axes = plt.subplots(nr, 3, figsize=(11, 3.3 * nr))
    y0, y1 = tune[0], tune[-1]
    for ax, c in zip(axes.flat, cities):
        g = pairs[pairs.city == c]
        cum = g[(g.y0 >= y0) & (g.y1 <= y1)].groupby(["row", "col"]).dL.sum(min_count=len(tune) - 1)
        grid = cum.unstack("col").to_numpy()
        im = ax.imshow(grid, cmap="RdBu_r", vmin=-30, vmax=30)
        ax.set_title(f"{c}: built-up {y0}->{y1} (pp)", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(im, ax=axes, shrink=0.6)
    fig.savefig(out / "fig_g1_cumulative_dL_maps.png", dpi=130)
    plt.close(fig)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="aqf-land", description=__doc__)
    ap.add_argument("--force", action="store_true", help="re-read the TNAQ GeoTIFFs instead of the cache")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("cache", help="read the city boxes from every TNAQ NO2 and S2 file into data/processed/")
    sub.add_parser("provenance", help="TROPOMI processor versions: what the local export records")
    sub.add_parser("layers", help="yearly land fractions, same-season NO2, per-pixel year-pair table")
    sub.add_parser("g1", help="gate G1: dL distribution and power over the 12x12 pixels (tuning years only)")
    sub.add_parser("all", help="cache, layers, provenance, g1")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s: %(message)s", datefmt="%H:%M:%S")
    logging.getLogger("rasterio").setLevel(logging.ERROR)

    from aqf_land.tnaq import load_land_config

    cfg = load_land_config()
    if args.cmd == "provenance":
        p = cfg["outputs"] / "P1_layers" / "pixel_year_pairs.csv"
        provenance(cfg, pd.read_csv(p) if p.exists() else None)
        return
    no2, s2 = cache(cfg, args.force)
    if args.cmd == "cache":
        return
    p = cfg["outputs"] / "P1_layers" / "pixel_year_pairs.csv"
    pairs = layers(cfg, no2, s2) if args.cmd in ("layers", "all") or not p.exists() else pd.read_csv(p)
    if args.cmd == "all":
        provenance(cfg, pairs)
    if args.cmd in ("g1", "all"):
        g1(cfg, pairs)
