"""5-day S5P NO2 composites for new windows (the 2025 test year), made exactly like the TNAQ export (needs Earth Engine).

Same method as AirQ/Scripts/Data_acquisition.js: COPERNICUS/S5P/OFFL/L3_NO2 tropospheric column, values below
-0.0005 mol/m^2 masked, masked mean over the orbits of each 5-day window plus the orbit count, clipped to Tamil
Nadu (FAO GAUL level 1), nodata -9999. Windows are 5-day steps anchored at 2018-01-01.

Files are downloaded straight onto the local export grid (0.01 deg, 414 x 547, origin 76.22 E / 13.54 N) into
data/tnaq_extra/TNAQ_S5P_NO2/, never into Dataset/. `--check YYYYMMDD` re-downloads one existing window and
compares it with the local file, to confirm the method matches.

    uv run --with earthengine-api --with requests python "Land Change/scripts/ee_s5p_no2_windows.py" --project <gcp> --year 2025
"""

from __future__ import annotations

import argparse
import io

import ee
import numpy as np
import pandas as pd
import rasterio
import requests

from aqf.data import ROOT
from aqf_land.tnaq import list_files, load_land_config

TRANSFORM = [0.01, 0, 76.22, 0, -0.01, 13.54]
WIDTH, HEIGHT = 414, 547
FLOOR = -0.0005
ANCHOR = pd.Timestamp("2018-01-01")


def composite(day: pd.Timestamp, aoi: ee.Geometry) -> ee.Image:
    d0 = ee.Date(day.strftime("%Y-%m-%d"))
    col = (ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2").filterBounds(aoi).filterDate(d0, d0.advance(5, "day"))
           .select(["tropospheric_NO2_column_number_density"], ["value"])
           .map(lambda img: img.updateMask(img.gte(FLOOR))))
    img = ee.Image(ee.Algorithms.If(col.size().gt(0),
                                    col.mean().rename("value").addBands(col.count().rename("n_obs")),
                                    ee.Image([0, 0]).rename(["value", "n_obs"]).selfMask()))
    return img.toFloat().clip(aoi).unmask(-9999)


def download(img: ee.Image, path) -> None:
    url = img.getDownloadURL({"format": "GEO_TIFF", "crs": "EPSG:4326", "crs_transform": TRANSFORM,
                              "dimensions": f"{WIDTH}x{HEIGHT}"})
    r = requests.get(url, timeout=300)
    r.raise_for_status()
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(io.BytesIO(r.content)) as src:
        a = src.read()
        prof = src.profile
        assert (src.width, src.height) == (WIDTH, HEIGHT) and np.allclose(list(src.transform)[:6], TRANSFORM), path
    prof.update(nodata=-9999.0, compress="deflate")
    with rasterio.open(path, "w", **prof) as dst:
        dst.write(a)
        dst.descriptions = ("value", "n_obs")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--year", type=int)
    ap.add_argument("--check", help="YYYYMMDD of an existing local window to reproduce and compare")
    args = ap.parse_args()
    ee.Initialize(project=args.project)
    cfg = load_land_config()
    aoi = (ee.FeatureCollection("FAO/GAUL_SIMPLIFIED_500m/2015/level1")
           .filter(ee.Filter.And(ee.Filter.eq("ADM0_NAME", "India"), ee.Filter.eq("ADM1_NAME", "Tamil Nadu")))
           .geometry())
    out_dir = ROOT / "data/tnaq_extra/TNAQ_S5P_NO2"
    if args.check:
        day = pd.Timestamp(args.check)
        tmp = ROOT / "data/tnaq_check" / f"s5p_no2_{args.check}.tif"
        download(composite(day, aoi), tmp)
        local = list_files(cfg, "NO2").set_index("date").path[day]
        with rasterio.open(local) as a, rasterio.open(tmp) as b:
            x, y = a.read(1), b.read(1)
        both = (x != -9999) & (y != -9999) & np.isfinite(x) & np.isfinite(y)
        print(f"check {args.check}: valid local {np.mean(x != -9999):.3f}, new {np.mean(y != -9999):.3f}, "
              f"same mask {np.mean((x == -9999) == (y == -9999)):.4f}, "
              f"max |diff| {np.abs(x[both] - y[both]).max():.3g} mol/m2, "
              f"median rel diff {np.median(np.abs(x[both] - y[both]) / np.abs(x[both]).clip(1e-6)):.2e}")
        return
    starts = pd.date_range(ANCHOR, f"{args.year}-12-31", freq="5D")
    starts = starts[starts.year == args.year]
    for day in starts:
        path = out_dir / f"s5p_no2_{day:%Y%m%d}.tif"
        if path.exists():
            continue
        download(composite(day, aoi), path)
        print("saved", path.name)


if __name__ == "__main__":
    main()
