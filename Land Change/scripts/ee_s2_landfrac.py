"""Land fractions per S5P block and year from 20 m Sentinel-2 SR: the real dL for ADR-001 (needs Earth Engine).

For every city box in configs/land.yaml and every year, this script:
1. Builds a same-season median composite of COPERNICUS/S2_SR_HARMONIZED (Jan-Apr by default). Clouds are masked
   with Cloud Score+ cs_cdf >= 0.60, as in AirQ/Scripts/Data_acquisition.js.
2. Classifies each 20 m pixel with the land.yaml rules: water MNDWI > 0, vegetation NDVI >= 0.3, built-up NDBI > 0,
   other. The order is the same as aqf_land.layers.classify.
3. Takes the class fractions and the valid fraction of each 0.04 deg block (the same blocks as the local analysis).
It exports one CSV per year to Google Drive, in the folder given by --folder.

--classifier dynamicworld instead takes the same-season mode of the GOOGLE/DYNAMICWORLD/V1 label (10 m, from
the same S2 scenes) and maps it to the four classes: water; trees, grass, flooded vegetation, crops and shrub ->
vegetation; built -> builtup; bare and snow -> other. The 1 km proxy showed that an index rule on one season
follows rainfall (the 2019 drought reads as built-up), so this is the more robust choice for dL.

Run it from the repository root (Earth Engine account needed; authenticate once):
    uv run --with earthengine-api earthengine authenticate
    uv run --with earthengine-api python "Land Change/scripts/ee_s2_landfrac.py" --project <gcp-project> --years 2019 2020 2021 2022 2023 2024 2025
Then download the CSVs into `Land Change/data/ee_landfrac/`. After that, `uv run aqf-land all` uses them
instead of the 1 km proxy.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

import ee

from aqf_land.tnaq import box_bounds, grid_boxes, load_land_config


def block_features(cfg: dict) -> ee.FeatureCollection:
    boxes, t = grid_boxes(cfg)
    b = cfg["grid"]["block"]
    step = b * t.a
    feats = []
    for box in boxes:
        s, w, n, e = box_bounds(box, t)
        for r in range(box.rows // b):
            for c in range(box.cols // b):
                north, west = n - r * step, w + c * step
                geom = ee.Geometry.Rectangle([west, north - step, west + step, north], "EPSG:4326", False)
                feats.append(ee.Feature(geom, {"city": box.city, "row": r, "col": c}))
    return ee.FeatureCollection(feats)


def composite(year: int, months: list[int], region: ee.Geometry) -> ee.Image:
    start = ee.Date.fromYMD(year, min(months), 1)
    end = ee.Date.fromYMD(year, max(months), 1).advance(1, "month")
    cs = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED")
    col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterBounds(region).filterDate(start, end)
           .linkCollection(cs, ["cs_cdf"])
           .map(lambda img: img.select(["B3", "B4", "B8", "B11"], ["green", "red", "nir", "swir1"])
                .multiply(0.0001).updateMask(img.select("cs_cdf").gte(0.60))))
    return col.median()


def classes(img: ee.Image) -> ee.Image:
    mndwi = img.normalizedDifference(["green", "swir1"])
    ndvi = img.normalizedDifference(["nir", "red"])
    ndbi = img.normalizedDifference(["swir1", "nir"])
    water = mndwi.gt(0)
    veg = ndvi.gte(0.3).And(water.Not())
    built = ndbi.gt(0).And(water.Not()).And(veg.Not())
    other = water.Not().And(veg.Not()).And(built.Not())
    valid = img.select("swir1").mask().gt(0)
    one = ee.Image.constant(1)
    # class bands are masked where the composite is missing, so their block means are fractions of valid pixels
    return (ee.Image.cat([water, veg, built, other]).rename(["water", "vegetation", "builtup", "other"])
            .updateMask(valid).addBands(one.updateMask(valid).unmask(0).rename("valid_frac")))


def classes_dynamicworld(year: int, months: list[int], region: ee.Geometry) -> ee.Image:
    start = ee.Date.fromYMD(year, min(months), 1)
    end = ee.Date.fromYMD(year, max(months), 1).advance(1, "month")
    label = (ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1").filterBounds(region).filterDate(start, end)
             .select("label").reduce(ee.Reducer.mode()))
    k = label.remap([0, 1, 2, 3, 4, 5, 6, 7, 8], [0, 1, 1, 1, 1, 1, 2, 3, 3])
    valid = label.mask().gt(0)
    bands = [k.eq(i) for i in range(4)]
    return (ee.Image.cat(bands).rename(["water", "vegetation", "builtup", "other"])
            .addBands(ee.Image.constant(1).updateMask(valid).unmask(0).rename("valid_frac")))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True, help="Google Cloud project registered for Earth Engine")
    ap.add_argument("--years", type=int, nargs="+", required=True)
    ap.add_argument("--folder", default="aqf_land_ee", help="Google Drive folder name (flat)")
    ap.add_argument("--scale", type=float, default=20.0, help="metres; B11 is native 20 m")
    ap.add_argument("--classifier", choices=["indices", "dynamicworld"], default="indices")
    ap.add_argument("--sites", type=Path, help="YAML {name: [lat, lon]} to use instead of the configured sites")
    ap.add_argument("--out-dir", type=Path, help="for --direct (default Land Change/data/ee_landfrac)")
    ap.add_argument("--direct", action="store_true",
                    help="fetch each table now (getInfo) into Land Change/data/ee_landfrac/ instead of a Drive export")
    args = ap.parse_args()
    ee.Initialize(project=args.project)
    cfg = load_land_config()
    if args.sites:
        cfg["cities"] = yaml.safe_load(args.sites.read_text())
    blocks = block_features(cfg)
    region = blocks.geometry().bounds()
    for y in args.years:
        img = (classes_dynamicworld(y, cfg["season_months"], region) if args.classifier == "dynamicworld"
               else classes(composite(y, cfg["season_months"], region)))
        table = img.reduceRegions(collection=blocks, reducer=ee.Reducer.mean(), scale=args.scale, tileScale=4)
        table = table.map(lambda f: f.set("year", y).setGeometry(None))
        if args.direct:
            import pandas as pd

            from aqf.data import ROOT

            out = (args.out_dir or ROOT / "Land Change/data/ee_landfrac") / f"landfrac_{y}.csv"
            out.parent.mkdir(parents=True, exist_ok=True)
            rows = [f["properties"] for f in table.getInfo()["features"]]
            cols = ["city", "row", "col", "year", "valid_frac", "water", "vegetation", "builtup", "other"]
            pd.DataFrame(rows).reindex(columns=cols).to_csv(out, index=False)
            label = ("Earth Engine Dynamic World V1, same-season mode label, 10 m" if args.classifier == "dynamicworld"
                     else f"Earth Engine S2_SR_HARMONIZED index rules, {args.scale:g} m")
            (out.parent / "source.json").write_text(json.dumps({"label": label, "classifier": args.classifier}))
            print(f"{y}: {len(rows)} blocks -> {out}")
            continue
        task = ee.batch.Export.table.toDrive(
            collection=table, description=f"aqf_land_{args.classifier}_{y}", folder=args.folder,
            fileNamePrefix=f"landfrac_{y}",
            fileFormat="CSV", selectors=["city", "row", "col", "year", "valid_frac", "water", "vegetation",
                                         "builtup", "other"])
        task.start()
        print(f"{y}: started task {task.id} -> Drive/{args.folder}/landfrac_{y}.csv")


if __name__ == "__main__":
    main()
