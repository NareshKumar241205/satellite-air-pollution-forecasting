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

    return read_cities(cfg, cfg.get("gas", "NO2"), force), read_cities(cfg, "S2", force)


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
    from aqf_land.layers import ALT_LAYERS, CLASSES, land_fractions, land_from_csv, land_from_layer, no2_seasonal, pair_table
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
        layer = cfg.get("land_layer", "dw_label")
        if layer in ALT_LAYERS:
            land = land_from_layer(ROOT, layer, c, (len(lat), len(lon)), years)
            sources[c] = f"alternative layer {layer} {ALT_LAYERS[layer]}"
        else:
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
        "gate": {"design": "pooled long difference", **{f"passes_at_{e:g}_pct_per_pp": bool(
            pooled.loc[f"long_{cfg['years']['tuning'][0]}_{cfg['years']['tuning'][-1]}", "mde_pct_per_pp"] <= e)
            for e in cfg["g1"].get("effects_reported", [cfg["g1"]["effect_pct_per_pp"]])}},
        "note": "passes uses the pooled MDE (the gate as defined); the city-weighted MDE is a sensitivity only"},
        indent=2))
    for name, r in pooled.iterrows():
        log.info("G1 %-18s signal sd %.2f pp (noise %.1f), MDE %s %%/pp, needs true sd >= %.2f pp -> %s", name,
                 r.dL_signal_sd, r.dL_noise_sd, f"{r.mde_pct_per_pp:.2f}", r.sd_required_pp,
                 "PASS" if r.passes else "FAIL")
    return res


def _deweathered(cfg: dict, tune: list[int]) -> dict | None:
    """{(site, year): (Hb, Wb)} deweathered 100*log NO2, or None without ERA5 for every tuning year."""
    from aqf_land.deweather import deweather_site
    from aqf_land.tnaq import block_mean, read_cities

    f = ROOT / "Land Change/results/P1_era5/era5_city_windows.csv"      # shared by every layer
    if not f.exists():
        return None
    era = pd.read_csv(f, parse_dates=["window"])
    if not set(tune) <= set(era.window.dt.year):
        return None
    no2 = read_cities(cfg, cfg.get("gas", "NO2"))
    g = cfg["grid"]
    out = {}
    for s in cfg["cities"]:
        met = era[era.city == s].set_index("window")
        blocks = block_mean(no2[s][:, 0], g["block"], g["min_valid_frac"])
        for y, v in deweather_site(blocks, no2["dates"], met, tune, cfg["season_months"], cfg["min_windows"]).items():
            out[(s, y)] = v
    return out


def panel(cfg: dict, n_boot: int = 500, n_placebo: int = 500, exclude_outcome_years: tuple = ()) -> None:
    from aqf_land import panel as pn
    from aqf_land.layers import CLASSES

    excl = tuple(sorted(exclude_outcome_years))
    tag = "" if not excl else "_POSTHOC_D7_without_" + "_".join(map(str, excl))
    out = cfg["outputs"] / f"C2_panel_EXPLORATORY{tag}"
    out.mkdir(parents=True, exist_ok=True)
    z = np.load(cfg["outputs"] / "P1_layers" / "land_fractions_no2_season.npz")
    tune = cfg["years"]["tuning"]                                         # 2024 and 2025 are never read
    k = CLASSES.index(cfg["dL_class"])
    variants = {"raw Jan-Apr NO2 (log of seasonal mean)": None}
    dw = _deweathered(cfg, tune)
    if dw is not None:
        variants["ERA5-deweathered Jan-Apr NO2 (mean of window logs)"] = dw
    rng = np.random.default_rng(cfg.get("seed", 0))
    rows, notes = [], []

    def add(variant, model, term, est, draws):
        lo, hi = np.percentile(draws, [2.5, 97.5])
        rows.append({"outcome": variant, "model": model, "term": term, "estimate": est,
                     "boot_se": float(np.std(draws, ddof=1)), "ci95_lo": lo, "ci95_hi": hi, "n_draws": len(draws)})

    ll_years = [y for y in tune[1:-1] if y not in excl]          # outcomes that have both a lead and a lag
    ll_label = f"lead/lag ({ll_years[0]}-{ll_years[-1]} outcomes)" if len(ll_years) > 1 else f"lead/lag ({ll_years[0]})"
    for variant, ylog in variants.items():
        p = pn.build_panel(z, list(cfg["cities"]), tune, k, ylog)
        if excl:                                   # D7: drop those outcomes; their land stays as leads/lags
            p.loc[p.year.isin(excl), "y"] = np.nan
        b = pn.fit(p, ["L", "WL"])
        bs = pn.block_bootstrap(p, ["L", "WL"], n_boot, rng)
        add(variant, "panel two-way FE + SLX", "L (own built-up)", b[0], bs[:, 0])
        add(variant, "panel two-way FE + SLX", "WL (neighbour built-up)", b[1], bs[:, 1])
        bl = pn.fit(p, ["L", "L_lead", "L_lag", "WL"])
        bls = pn.block_bootstrap(p, ["L", "L_lead", "L_lag", "WL"], n_boot, rng)
        for j, term in enumerate(["L", "L_lead (placebo: must be ~0)", "L_lag", "WL"]):
            add(variant, ll_label, term, bl[j], bls[:, j])
        b_ld, d = pn.long_difference_fit(p, tune[0], tune[-1])
        add(variant, f"long difference {tune[0]}->{tune[-1]} + site FE (G1 design)", "dL", b_ld,
            pn.long_difference_bootstrap(d, n_boot, rng))
        plac = pn.placebo(p, n_placebo, rng)
        p_perm = float((np.abs(plac) >= abs(b[0])).mean())
        notes.append(f"- {variant}: spatial placebo ({n_placebo} draws) median b = {np.median(plac):.3f}, 95 % of |b| "
                     f"below {np.percentile(np.abs(plac), 95):.3f}; permutation p for the panel b = {p_perm:.3f}.")
        log.info("C2 (exploratory, %s): b = %.3f %%/pp [%.3f, %.3f], long-diff b = %.3f, lead = %.3f, placebo p = %.3f",
                 variant, b[0], *np.percentile(bs[:, 0], [2.5, 97.5]), b_ld, bl[1], p_perm)
    res = pd.DataFrame(rows)
    res.to_csv(out / "estimates.csv", index=False)
    mde = json.loads((cfg["outputs"] / "G1_power" / "G1_RESULT.json").read_text())["pooled"]
    mde_ld = mde[f"long_{tune[0]}_{tune[-1]}"]["mde_pct_per_pp"]
    head = pn.EXPLORATORY if dw is None else pn.EXPLORATORY.replace(
        "raw Jan-Apr NO2, no ERA5 deweathering", "raw and ERA5-deweathered Jan-Apr NO2")
    if excl:
        head = f"POST-HOC DIAGNOSTIC D7, outcome years {list(excl)} removed (DIAGNOSTICS_PLAN.md addendum). " + head
    lines = [f"# C2 panel: {head}", "",
             f"G1 detection limit (pooled long difference, 80 % power): **{mde_ld:.2f} % {cfg.get("gas", "NO2")} per pp**. "
             "Estimates below that size can't be told apart from zero by this design; a null here is not evidence of "
             "no effect.",
             f"Sites: {len(cfg['cities'])}, years {tune[0]}-{tune[-1]} (Jan-Apr). Spatial-block bootstrap: {n_boot} "
             "draws of 4x4-pixel blocks within site.", "",
             res.round(3).to_markdown(index=False), "", *notes]
    (out / "C2_PANEL_EXPLORATORY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def diagnostics(cfg: dict) -> None:
    from aqf_land import diag
    from aqf_land import panel as pn
    from aqf_land.layers import CLASSES

    out = cfg["outputs"] / "D_lead_diagnostics"
    out.mkdir(parents=True, exist_ok=True)
    z = np.load(cfg["outputs"] / "P1_layers" / "land_fractions_no2_season.npz")
    tune = cfg["years"]["tuning"]
    p = pn.build_panel(z, list(cfg["cities"]), tune, CLASSES.index(cfg["dL_class"]))
    parts = ["# Lead-test diagnostics (plan: `Land Change/DIAGNOSTICS_PLAN.md`)", ""]
    cov_dir = ROOT / "data/ee_covariates"
    if any(cov_dir.glob("covariates_*.csv")):
        cov = pd.concat([pd.read_csv(f) for f in sorted(cov_dir.glob("covariates_*.csv"))])
        d1 = diag.d1_covariates(p, cov)
        d1.to_csv(out / "D1_covariates.csv", index=False)
        piv = d1.pivot(index="covariate", columns="with", values="r_within_site_year").round(3)
        parts += ["## D1: scene/atmosphere covariates (within site-year r)", "", piv.to_markdown(), ""]
    d2 = diag.d2_level_dependence(p)
    d3 = diag.d3_pixel_trends(p)
    d4 = diag.d4_leave_one_out(p)
    d5 = diag.d5_drop_edges(p)
    d6 = diag.d6_spatial(p)
    for name, obj in (("D2", d2), ("D4", d4)):
        obj.to_csv(out / f"{name}.csv", index=False)
    (out / "D3_D5_D6.json").write_text(json.dumps({"D3": d3, "D5": d5, "D6": d6}, indent=2))
    parts += ["## D2: level dependence of the yearly label change (within site)", "", d2.round(3).to_markdown(index=False),
              "", "## D3: pixel-specific linear trends", "", f"`{json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in d3.items()})}`",
              "", "## D4: leave one out (lead coefficient)", "", d4.round(3).to_markdown(index=False), "",
              "## D5: drop the outer ring of blocks", "", f"`{json.dumps({k: round(v, 3) for k, v in d5.items()})}`", "",
              "## D6: residual spatial structure", "", f"`{json.dumps({k: round(v, 3) for k, v in d6.items()})}`"]
    (out / "DIAGNOSTICS.md").write_text("\n".join(parts) + "\n", encoding="utf-8")
    log.info("diagnostics -> %s", out)


def b_decision(cfg: dict) -> dict:
    """Apply Land Change/DECISION_RULES.md section B mechanically to the per-layer G1 and C2 results.
    'Clean' is read on the primary (raw NO2) panel outcome."""
    rows, base = [], cfg["outputs"]
    tune = cfg["years"]["tuning"]
    ld = f"long_{tune[0]}_{tune[-1]}"
    for layer in ("dw_label", "dw_prob_trail3", "obt"):
        root = base if layer == "dw_label" else base / "layers" / layer
        g1f, c2f = root / "G1_power" / "G1_RESULT.json", root / "C2_panel_EXPLORATORY" / "estimates.csv"
        if not (g1f.exists() and c2f.exists()):
            rows.append({"layer": layer, "available": False})
            continue
        g = json.loads(g1f.read_text())
        e = pd.read_csv(c2f)
        e = e[e.outcome.str.startswith("raw")]
        lead = e[e.term.str.startswith("L_lead")].iloc[0]
        wl = e[(e.model == "panel two-way FE + SLX") & e.term.str.startswith("WL")].iloc[0]
        mde = g["pooled"][ld]["mde_pct_per_pp"]
        clean = bool(lead.ci95_lo <= 0 <= lead.ci95_hi and wl.ci95_lo <= 0 <= wl.ci95_hi)
        rows.append({"layer": layer, "available": True, "mde": mde, "g1_pass_0.5": bool(mde is not None and mde <= 0.5),
                     "lead": lead.estimate, "lead_lo": lead.ci95_lo, "lead_hi": lead.ci95_hi, "WL": wl.estimate,
                     "WL_lo": wl.ci95_lo, "WL_hi": wl.ci95_hi, "clean": clean})
    r = pd.DataFrame(rows)
    steady = r[r.layer.isin(["dw_prob_trail3", "obt"]) & r.available.astype(bool)]
    steady = steady[steady.clean.astype(bool)]
    chosen = None if steady.empty else steady.sort_values("mde").layer.iloc[0]
    passes = bool(chosen is not None and r.set_index("layer").loc[chosen, "g1_pass_0.5"])
    mde_ref = r.set_index("layer").loc[chosen or "dw_label", "mde"]
    dec = {"steadier_layer": chosen, "rule_passes": passes, "B": "confirmatory C2" if passes else 0,
           "detection_limit_pct_per_pp": mde_ref, "layers": json.loads(r.to_json(orient="records")),
           "rule": "Land Change/DECISION_RULES.md section B"}
    (base / "B_decision.json").write_text(json.dumps(dec, indent=2))
    log.info("B decision: steadier layer %s, rule passes %s -> B = %s", chosen, passes, dec["B"])
    return dec


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
    ap.add_argument("--gas", default="NO2", choices=["NO2", "CO", "SO2"],
                    help="gas for layers/g1/panel/forecast; non-NO2 results go to Land Change/results/gas/<gas>/")
    ap.add_argument("--layer", default="dw_label", help="built-up layer: dw_label (default), dw_prob_trail3 or obt; "
                    "alternatives write to Land Change/results/layers/<layer>/")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("cache", help="read the city boxes from every TNAQ NO2 and S2 file into data/processed/")
    sub.add_parser("provenance", help="TROPOMI processor versions: what the local export records")
    sub.add_parser("layers", help="yearly land fractions, same-season NO2, per-pixel year-pair table")
    sub.add_parser("g1", help="gate G1: dL distribution and power over the 12x12 pixels (tuning years only)")
    sub.add_parser("figures", help="report figures from the saved result tables")
    sub.add_parser("bdecision", help="apply decision rule B to the per-layer G1 and C2 results")
    sf = sub.add_parser("forecast", help="DECISION_RULES C-E: validation ladder, fusion, then the one-time 2025 test")
    sf.add_argument("--rescore-test", action="store_true", help="only to reproduce identical numbers; never to tune")
    sf.add_argument("--validation-only", action="store_true", help="debug run of the validation stage; 2025 untouched")
    sub.add_parser("diag", help="lead-test diagnostics D1-D6 (Land Change/DIAGNOSTICS_PLAN.md)")
    spn = sub.add_parser("panel", help="C2 panel, EXPLORATORY: two-way FE + SLX, block bootstrap, placebo, lead/lag")
    spn.add_argument("--exclude-outcome-years", type=int, nargs="*", default=[],
                     help="post-hoc diagnostic D7 only (e.g. 2020); results go to a separate, labelled folder")
    sub.add_parser("all", help="cache, layers, provenance, g1")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s: %(message)s", datefmt="%H:%M:%S")
    logging.getLogger("rasterio").setLevel(logging.ERROR)

    from aqf_land.tnaq import load_land_config

    cfg = load_land_config()
    if args.gas != "NO2":
        cfg["gas"] = args.gas
        cfg["outputs"] = cfg["outputs"] / "gas" / args.gas
    if args.layer != "dw_label":
        cfg["land_layer"] = args.layer
        cfg["outputs"] = cfg["outputs"] / "layers" / args.layer
    if args.cmd == "diag":
        diagnostics(cfg)
        return
    if args.cmd == "figures":
        from aqf_land import figures

        figures.run(cfg)
        return
    if args.cmd == "bdecision":
        b_decision(cfg)
        return
    if args.cmd == "forecast":
        from aqf_land import forecast_run

        forecast_run.run(cfg, rescore=args.rescore_test, val_only=args.validation_only)
        return
    if args.cmd == "panel":
        panel(cfg, exclude_outcome_years=tuple(args.exclude_outcome_years))
        return
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
