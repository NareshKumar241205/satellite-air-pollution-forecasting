"""Per block and year (Jan-Apr) Earth Engine layers for Land Change/DIAGNOSTICS_PLAN.md.

--what dwprob      mean Dynamic World `built` probability (input to the trailing 3-year median, layer S1)
--what obt         Open Buildings Temporal mean building_presence and building_height (layer S2; annual, 2016-2023)
--what covariates  D1 scene and atmosphere covariates: Dynamic World image count, Sentinel-2 cloud probability
                   (1 - Cloud Score+ cs_cdf), MAIAC AOD 0.55 um, and S2 visible brightness (median of mean B2-B4)

Each year is written to data/ee_<what>/<what>_<year>.csv with the columns city, row, col, year, <variables>.

    uv run --with earthengine-api python "Land Change/scripts/ee_block_layers.py" --project <gcp> --what dwprob --years 2017 ... 2024
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import ee
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from ee_s2_landfrac import block_features  # noqa: E402

from aqf.data import ROOT  # noqa: E402
from aqf_land.tnaq import load_land_config  # noqa: E402


def season(y: int, months: list[int]) -> tuple[ee.Date, ee.Date]:
    return ee.Date.fromYMD(y, min(months), 1), ee.Date.fromYMD(y, max(months), 1).advance(1, "month")


def dwprob(y, months, region):
    s, e = season(y, months)
    return ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1").filterBounds(region).filterDate(s, e).select("built").mean() \
        .rename("built_prob"), 20


def obt(y, months, region):
    col = ee.ImageCollection("GOOGLE/Research/open-buildings-temporal/v1").filterBounds(region) \
        .filter(ee.Filter.calendarRange(y, y, "year"))
    return col.select(["building_presence", "building_height"]).mosaic(), 20


def covariates(y, months, region):
    s, e = season(y, months)
    dw = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1").filterBounds(region).filterDate(s, e)
    n_dw = dw.select("label").count().rename("dw_images")
    cs = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED").filterBounds(region).filterDate(s, e)
    cloud = cs.select("cs_cdf").map(lambda i: ee.Image(1).subtract(i)).mean().rename("cloud_prob")
    aod = (ee.ImageCollection("MODIS/061/MCD19A2_GRANULES").filterBounds(region).filterDate(s, e)
           .select("Optical_Depth_055").mean().multiply(0.001).rename("aod055"))
    s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterBounds(region).filterDate(s, e)
          .linkCollection(cs, ["cs_cdf"])
          .map(lambda i: i.select(["B2", "B3", "B4"]).reduce(ee.Reducer.mean()).multiply(0.0001)
               .updateMask(i.select("cs_cdf").gte(0.6)).rename("brightness")))
    return ee.Image.cat([n_dw, cloud, aod.unmask(), s2.median()]), 100


WHAT = {"dwprob": dwprob, "obt": obt, "covariates": covariates}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--what", choices=list(WHAT), required=True)
    ap.add_argument("--years", type=int, nargs="+", required=True)
    args = ap.parse_args()
    ee.Initialize(project=args.project)
    cfg = load_land_config()
    blocks = block_features(cfg)
    region = blocks.geometry().bounds()
    out_dir = ROOT / f"data/ee_{args.what}"
    out_dir.mkdir(parents=True, exist_ok=True)
    for y in args.years:
        out = out_dir / f"{args.what}_{y}.csv"
        if out.exists():
            print(y, "exists")
            continue
        img, scale = WHAT[args.what](y, cfg["season_months"], region)
        table = img.reduceRegions(collection=blocks, reducer=ee.Reducer.mean(), scale=scale, tileScale=8) \
            .map(lambda f: f.set("year", y).setGeometry(None))
        for attempt in range(5):
            try:
                rows = [f["properties"] for f in table.getInfo()["features"]]
                break
            except Exception as err:                       # EEException or HttpError (computation timed out)
                if attempt == 4:
                    raise
                print("retry", y, err)
                time.sleep(15 * 2 ** attempt)
        pd.DataFrame(rows).to_csv(out, index=False)
        print(y, len(rows), "blocks ->", out)


if __name__ == "__main__":
    main()
