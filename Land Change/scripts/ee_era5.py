"""ERA5 scalars per site and 5-day window from Earth Engine's ECMWF/ERA5/HOURLY (the same ERA5 single levels as CDS).

Equivalent to cds_era5.py (download + reduce) without a CDS account:
- the hours are 06-09 UTC (the S5P overpass is ~08 UTC);
- the box is each site box widened by 0.25 deg, so it covers several 0.25 deg ERA5 cells;
- 5-day windows are anchored at 2018-01-01;
- values are box means of the hourly fields averaged over the window. Precipitation is summed, and wind speed is the mean of the hourly speeds.

Output columns match cds_era5.py: Land Change/results/P1_era5/era5_city_windows.csv.

    uv run --with earthengine-api python "Land Change/scripts/ee_era5.py" --project <gcp> --years 2019 2020 2021 2022 2023 2024 2025
"""

from __future__ import annotations

import argparse
import time

import ee
import pandas as pd

from aqf_land.tnaq import box_bounds, grid_boxes, load_land_config

BANDS = {"u_component_of_wind_10m": "u10", "v_component_of_wind_10m": "v10", "boundary_layer_height": "blh",
         "temperature_2m": "t2m", "surface_pressure": "sp", "total_precipitation": "tp",
         "total_cloud_cover": "tcc"}
ANCHOR = pd.Timestamp("2018-01-01")


def window_image(start: ee.Date) -> ee.Image:
    col = (ee.ImageCollection("ECMWF/ERA5/HOURLY").filterDate(start, start.advance(5, "day"))
           .filter(ee.Filter.calendarRange(6, 9, "hour")).select(list(BANDS), list(BANDS.values())))
    ws = col.map(lambda i: i.expression("hypot(u, v)", {"u": i.select("u10"), "v": i.select("v10")}).rename("ws10"))
    mean = col.select(["u10", "v10", "blh", "t2m", "sp", "tcc"]).mean()
    return (mean.addBands(col.select("tp").sum()).addBands(ws.mean())
            .addBands(ee.Image.constant(col.size()).rename("n_hours")).set("window", start.format("YYYY-MM-dd")))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--years", type=int, nargs="+", required=True)
    ap.add_argument("--pad", type=float, default=0.25)
    args = ap.parse_args()
    ee.Initialize(project=args.project)
    cfg = load_land_config()
    boxes, t = grid_boxes(cfg)
    feats = []
    for b in boxes:
        s, w, n, e = box_bounds(b, t)
        p = args.pad
        feats.append(ee.Feature(ee.Geometry.Rectangle([w - p, s - p, e + p, n + p], "EPSG:4326", False),
                                {"city": b.city}))
    sites = ee.FeatureCollection(feats)
    out = cfg["outputs"] / "P1_era5" / "era5_city_windows.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    frames = [pd.read_csv(out)] if out.exists() else []
    done = set(pd.to_datetime(frames[0].window).dt.year) if frames else set()
    for y in args.years:
        if y in done:
            print(y, "already in", out.name)
            continue
        starts = pd.date_range(ANCHOR, f"{y}-12-31", freq="5D")
        starts = starts[starts.year == y]
        rows = []
        for chunk in range(0, len(starts), 4):           # small requests: EE limits concurrent aggregations
            fc = ee.FeatureCollection([
                window_image(ee.Date(d.strftime("%Y-%m-%d"))).reduceRegions(sites, ee.Reducer.mean(), 27830)
                .map(lambda f, d=d: f.set("window", d.strftime("%Y-%m-%d")).setGeometry(None))
                for d in starts[chunk:chunk + 4]]).flatten()
            for attempt in range(5):
                try:
                    rows += [f["properties"] for f in fc.getInfo()["features"]]
                    break
                except ee.EEException as err:
                    if attempt == 4:
                        raise
                    print("retry", chunk, err)
                    time.sleep(10 * 2 ** attempt)
        df = pd.DataFrame(rows)
        frames.append(df[["city", "window", "u10", "v10", "blh", "t2m", "sp", "tp", "tcc", "ws10", "n_hours"]])
        pd.concat(frames, ignore_index=True).sort_values(["city", "window"]).to_csv(out, index=False)
        print(y, len(df), "rows ->", out)


if __name__ == "__main__":
    main()
