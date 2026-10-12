"""ERA5 for ADR-001: one request per year over the padded TNAQ box, then scalars per city and 5-day window.

download: one CDS request per year (reanalysis-era5-single-levels, hours around the 08 UTC S5P overpass,
          the TNAQ extent padded by `era5.pad_deg`), saved as data/era5/era5_<year>.nc. Existing files are skipped.
reduce:   box-mean scalars per city and window (the city box widened by one ERA5 cell, 0.25 deg, so each box
          covers several cells) -> Land Change/results/P1_era5/era5_city_windows.csv.

Needs a CDS account and ~/.cdsapirc (https://cds.climate.copernicus.eu/how-to-api), with the ERA5 licence accepted:
    uv run --with cdsapi python "Land Change/scripts/cds_era5.py" download
    uv run --with xarray --with netcdf4 python "Land Change/scripts/cds_era5.py" reduce
"""

from __future__ import annotations

import argparse
import zipfile
from pathlib import Path

import pandas as pd

from aqf.data import ROOT
from aqf_land.era5 import window_scalars
from aqf_land.tnaq import box_bounds, grid_boxes, load_land_config

BBOX = [76.22, 8.07, 80.36, 13.54]        # TNAQ export extent (W, S, E, N)
SHORT = {"10m_u_component_of_wind": "u10", "10m_v_component_of_wind": "v10", "boundary_layer_height": "blh",
         "2m_temperature": "t2m", "surface_pressure": "sp", "total_precipitation": "tp", "total_cloud_cover": "tcc"}


def download(cfg: dict, years: list[int], out_dir: Path) -> None:
    import cdsapi

    e = cfg["era5"]
    p = e["pad_deg"]
    area = [BBOX[3] + p, BBOX[0] - p, BBOX[1] - p, BBOX[2] + p]          # N, W, S, E
    c = cdsapi.Client()
    out_dir.mkdir(parents=True, exist_ok=True)
    for y in years:
        target = out_dir / f"era5_{y}.nc"
        if target.exists():
            print("exists", target)
            continue
        c.retrieve("reanalysis-era5-single-levels", {
            "product_type": ["reanalysis"], "variable": e["variables"], "year": [str(y)],
            "month": [f"{m:02d}" for m in range(1, 13)], "day": [f"{d:02d}" for d in range(1, 32)],
            "time": e["hours_utc"], "area": area, "data_format": "netcdf", "download_format": "unarchived",
        }, str(target))
        print("downloaded", target)


def open_year(path: Path):
    import xarray as xr

    if zipfile.is_zipfile(path):                 # the new CDS sometimes zips one file per variable type
        d = path.with_suffix("")
        with zipfile.ZipFile(path) as z:
            z.extractall(d)
        return xr.merge([xr.open_dataset(f) for f in sorted(d.glob("*.nc"))], compat="override")
    return xr.open_dataset(path)


def reduce(cfg: dict, in_dir: Path) -> None:
    boxes, t = grid_boxes(cfg)
    rows = []
    for path in sorted(in_dir.glob("era5_*.nc")):
        ds = open_year(path)
        tname = "valid_time" if "valid_time" in ds.dims else "time"
        times = pd.DatetimeIndex(ds[tname].values)
        lat, lon = ds["latitude"].values, ds["longitude"].values
        fields = {SHORT.get(v, v): ds[v].values for v in ds.data_vars if ds[v].ndim == 3}
        for box in boxes:
            df = window_scalars(times, fields, lat, lon, box_bounds(box, t), pad=0.25)
            df.insert(0, "city", box.city)
            rows.append(df)
        print("reduced", path.name)
    out = cfg["outputs"] / "P1_era5" / "era5_city_windows.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(rows, ignore_index=True).to_csv(out, index=False)
    print("wrote", out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["download", "reduce"])
    ap.add_argument("--years", type=int, nargs="+", help="default: era5.years in configs/land.yaml")
    ap.add_argument("--dir", type=Path, default=ROOT / "data/era5")
    args = ap.parse_args()
    cfg = load_land_config()
    if args.step == "download":
        download(cfg, args.years or cfg["era5"]["years"], args.dir)
    else:
        reduce(cfg, args.dir)


if __name__ == "__main__":
    main()
