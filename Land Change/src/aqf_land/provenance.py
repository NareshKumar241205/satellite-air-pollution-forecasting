"""Which TROPOMI NO2 processor versions are in the S5P export?

What can be checked locally:
- The GeoTIFFs: every file's dataset and band tags (Earth Engine writes no product metadata into exports).
- The manifests: their columns.
- The acquisition script: the Earth Engine collection and band.

The collection COPERNICUS/S5P/OFFL/L3_NO2 is built from the operational OFFL L2 stream as it was processed at
the time, not from one consistent reprocessing, so 2019-2024 spans several processor versions. The per-image
version is an Earth Engine image property, so the exact count per window needs scripts/ee_s5p_versions.py.
Until that runs, the approximate public timeline below is used and flagged as such.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

# Approximate OFFL NO2 processor timeline (start dates) from the S5P product readme / release notes.
# Only the v1 -> v2 boundary matters for this design; the rest is context. VERIFY with ee_s5p_versions.py.
TIMELINE = [
    ("2018-04-30", "1.2.x"),
    ("2019-03-20", "1.3.x"),
    ("2020-12-02", "1.4.0"),
    ("2021-07-01", "2.2.0"),   # FRESCO-wide clouds + new surface albedo: tropospheric columns up, most over polluted areas
    ("2021-11-14", "2.3.1"),
    ("2022-07-17", "2.4.0"),
    ("2023-03-12", "2.5.0 and later (approximate)"),
]
V2_START = pd.Timestamp("2021-07-01")


def version_at(dates: pd.DatetimeIndex) -> np.ndarray:
    starts = pd.DatetimeIndex([s for s, _ in TIMELINE])
    idx = starts.searchsorted(dates, side="right") - 1
    return np.array([TIMELINE[i][1] if i >= 0 else "none" for i in idx])


def file_metadata(paths: list[Path]) -> pd.DataFrame:
    """Distinct tag sets over all files (normally one row: Earth Engine writes only AREA_OR_POINT)."""
    seen = {}
    for p in paths:
        with rasterio.open(p) as src:
            key = (repr(sorted(src.tags().items())), repr(src.descriptions),
                   repr([sorted(src.tags(i).items()) for i in range(1, src.count + 1)]))
        seen.setdefault(key, []).append(p.name)
    return pd.DataFrame([{"dataset_tags": k[0], "band_names": k[1], "band_tags": k[2], "n_files": len(v),
                          "example": v[0]} for k, v in seen.items()])


def acquisition_info(script: Path) -> dict:
    if not script.exists():
        return {"script": f"{script} (not found)"}
    txt = script.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"NO2:\s*\{\s*collection:\s*'([^']+)',\s*band:\s*'([^']+)'", txt)
    return {"script": str(script), "collection": m.group(1) if m else "?", "band": m.group(2) if m else "?"}


def season_versions(files: pd.DataFrame, months: list[int]) -> pd.DataFrame:
    """Season windows per year and (approximate) processor version."""
    s = files[files.month.isin(months)].copy()
    s["version"] = version_at(pd.DatetimeIndex(s.date))
    return s.groupby(["year", "version"]).size().rename("windows").reset_index()


def common_step(pairs: pd.DataFrame) -> pd.DataFrame:
    """Mean dNO2 (%) per year pair in each city: a retrieval step would show up in every city at once."""
    t = pairs.groupby(["y0", "y1", "city"]).dNO2.mean().unstack("city")
    t["mean_all_cities"] = t.mean(1)
    return t.round(1)
