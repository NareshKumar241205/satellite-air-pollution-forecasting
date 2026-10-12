"""Exact TROPOMI NO2 processor versions behind the TNAQ export (needs Earth Engine).

The local GeoTIFFs and manifests don't record the version (Land Change/results/P1_provenance/PROVENANCE.md).
This script reads the image properties of COPERNICUS/S5P/OFFL/L3_NO2 over the TNAQ area and counts images per
month and version. It first prints the property names of one image, so the version keys are found rather than
guessed. The output is Land Change/results/P1_provenance/ee_processor_versions.csv.

    uv run --with earthengine-api python "Land Change/scripts/ee_s5p_versions.py" --project <gcp-project>
"""

from __future__ import annotations

import argparse

import ee
import pandas as pd

from aqf_land.tnaq import load_land_config

BBOX = [76.22, 8.07, 80.36, 13.54]        # TNAQ export extent (W, S, E, N)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--start", default="2018-07-01")
    ap.add_argument("--end", default="2026-01-01")
    args = ap.parse_args()
    ee.Initialize(project=args.project)
    cfg = load_land_config()
    col = (ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2").filterBounds(ee.Geometry.Rectangle(BBOX))
           .filterDate(args.start, args.end))
    names = ee.Image(col.first()).propertyNames().getInfo()
    print("image properties:", sorted(names))
    keys = [k for k in names if "VERSION" in k.upper() or "PROCESSOR" in k.upper() or k.upper() == "PRODUCT_ID"]
    rows = []
    for month in pd.date_range(args.start, args.end, freq="MS")[:-1]:
        sub = col.filterDate(month.strftime("%Y-%m-%d"), (month + pd.offsets.MonthBegin()).strftime("%Y-%m-%d"))
        for k in keys:
            if k.upper() == "PRODUCT_ID":
                continue
            hist = sub.aggregate_histogram(k).getInfo()
            rows += [{"month": month.strftime("%Y-%m"), "property": k, "value": v, "images": n} for v, n in hist.items()]
        print(month.strftime("%Y-%m"), {r["value"]: r["images"] for r in rows if r["month"] == month.strftime("%Y-%m")})
    out = cfg["outputs"] / "P1_provenance" / "ee_processor_versions.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print("wrote", out)


if __name__ == "__main__":
    main()
