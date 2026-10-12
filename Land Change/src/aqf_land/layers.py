"""The two clocks on the 0.04 deg block grid: yearly land fractions (L) and same-season NO2 means, plus their changes.

Land: a same-season (Jan-Apr) median composite per sub-pixel and year, classified by the rules in
configs/land.yaml, then the class fractions of each block. dL is the change in percentage points between
consecutive years.

NO2: each window is block-averaged first, then averaged over the same season of each year. dNO2 is the log
change between consecutive years, reported in percent (100 * dlog).

The local S2 export is 0.01 deg (~1 km), so one sub-pixel changing class moves a block fraction by 6.25 pp.
That is a provisional proxy; the real dL comes from 10-20 m S2_SR_HARMONIZED composites
(scripts/ee_s2_landfrac.py).
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

CLASSES = ["water", "vegetation", "builtup", "other"]
NO2_FLOOR = 1e-6          # mol/m^2; seasonal means below this are not logged (none occur in practice)


def seasonal_composite(stack: np.ndarray, dates: pd.DatetimeIndex, year: int, months: list[int],
                       min_windows: int) -> np.ndarray:
    """Median over the season windows of one year: (bands, H, W), NaN where fewer than min_windows are valid."""
    sel = (dates.year == year) & dates.month.isin(months)
    s = stack[sel]
    if not len(s):
        return np.full(stack.shape[1:], np.nan, np.float32)
    n = np.isfinite(s).all(1).sum(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)       # all-NaN slices are expected (sea, clouds)
        med = np.nanmedian(s, 0)
    med[:, n < min_windows] = np.nan
    return med


def classify(comp: np.ndarray, band_names: list[str]) -> np.ndarray:
    """(H, W) class index into CLASSES, -1 where the composite is missing. First match wins, as in land.yaml."""
    b = dict(zip(band_names, comp))
    with np.errstate(invalid="ignore", divide="ignore"):
        mndwi = (b["green"] - b["swir1"]) / (b["green"] + b["swir1"])
        ndvi = (b["nir"] - b["red"]) / (b["nir"] + b["red"])
        ndbi = (b["swir1"] - b["nir"]) / (b["swir1"] + b["nir"])
    out = np.full(mndwi.shape, 3, np.int8)
    out[ndbi > 0] = 2
    out[ndvi >= 0.3] = 1
    out[mndwi > 0] = 0
    out[~np.isfinite(comp).all(0)] = -1
    return out


def class_fractions(cls: np.ndarray, block: int, min_valid_frac: float) -> np.ndarray:
    """(len(CLASSES), Hb, Wb) fraction of valid sub-pixels in each class, NaN where too few are valid."""
    from aqf_land.tnaq import block_mean

    onehot = np.stack([(cls == k).astype(np.float32) for k in range(len(CLASSES))])
    onehot[:, cls < 0] = np.nan
    return block_mean(onehot, block, min_valid_frac)


def land_fractions(s2: np.ndarray, dates: pd.DatetimeIndex, years: list[int], cfg: dict) -> dict[int, np.ndarray]:
    """{year: (C, Hb, Wb)} for one city."""
    g = cfg["grid"]
    out = {}
    for y in years:
        comp = seasonal_composite(s2, dates, y, cfg["season_months"], cfg["min_windows"])
        out[y] = class_fractions(classify(comp, cfg["tnaq"]["s2_bands"]), g["block"], g["min_valid_frac"])
    return out


def no2_seasonal(no2: np.ndarray, dates: pd.DatetimeIndex, years: list[int], cfg: dict) -> dict[int, np.ndarray]:
    """{year: (Hb, Wb)} same-season mean of block-averaged NO2 (mol/m^2) for one city. no2: (T, 1, H, W)."""
    from aqf_land.tnaq import block_mean

    g = cfg["grid"]
    blocks = block_mean(no2[:, 0], g["block"], g["min_valid_frac"])      # (T, Hb, Wb)
    out = {}
    for y in years:
        s = blocks[(dates.year == y) & dates.month.isin(cfg["season_months"])]
        n = np.isfinite(s).sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            m = np.nansum(s, 0) / n
        out[y] = np.where(n >= cfg["min_windows"], m, np.nan)
    return out


def pair_table(city: str, lat: np.ndarray, lon: np.ndarray, land: dict[int, np.ndarray],
               no2: dict[int, np.ndarray], dl_class: str) -> pd.DataFrame:
    """One row per block and consecutive-year pair: dL (pp, all classes) and dNO2 (% log change)."""
    k = CLASSES.index(dl_class)
    years = sorted(set(land) & set(no2))
    rows = []
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    for y0, y1 in zip(years[:-1], years[1:]):
        d = {"city": city, "y0": y0, "y1": y1, "row": np.indices(LA.shape)[0].ravel(),
             "col": np.indices(LA.shape)[1].ravel(), "lat": LA.ravel(), "lon": LO.ravel(),
             "L0": land[y0][k].ravel() * 100, "dL": (land[y1][k] - land[y0][k]).ravel() * 100}
        for j, c in enumerate(CLASSES):
            d[f"d_{c}"] = (land[y1][j] - land[y0][j]).ravel() * 100
        a, b = no2[y0].ravel(), no2[y1].ravel()
        ok = (a > NO2_FLOOR) & (b > NO2_FLOOR)
        with np.errstate(invalid="ignore", divide="ignore"):
            d["no2_0"] = a
            d["dNO2"] = np.where(ok, 100 * np.log(b / a), np.nan)
        rows.append(pd.DataFrame(d))
    return pd.concat(rows, ignore_index=True)


def land_from_csv(csv_dir, city: str, shape: tuple[int, int], years: list[int],
                  min_valid_frac: float) -> dict[int, np.ndarray] | None:
    """{year: (C, Hb, Wb)} from the Earth Engine 20 m tables (scripts/ee_s2_landfrac.py), or None if a year is
    missing. Blocks with valid_frac below min_valid_frac are NaN."""
    files = {int(p.stem.split("_")[-1]): p for p in csv_dir.glob("landfrac_*.csv")}
    if not set(years) <= set(files):
        return None
    out = {}
    for y in years:
        d = pd.read_csv(files[y])
        d = d[d.city == city]
        a = np.full((len(CLASSES), *shape), np.nan)
        ok = d.valid_frac >= min_valid_frac
        for k, c in enumerate(CLASSES):
            a[k, d.row[ok], d.col[ok]] = d[c][ok]
        out[y] = a
    return out


ALT_LAYERS = {
    # name: (data dir, file prefix, column, trailing window in years); see Land Change/DIAGNOSTICS_PLAN.md
    "dw_prob_trail3": ("data/ee_dwprob", "dwprob", "mean", 3),          # EE names a single-band mean "mean"
    "obt": ("data/ee_obt", "obt", "building_presence", 1),
    # X4: night-time lights as an emission proxy. v = 0.1 ln(1 + radiance), so L = 100 v = 10 ln(1 + radiance)
    # and a 1-unit change is about a 10 % change in lights.
    "viirs_log": ("data/ee_viirs", "viirs", "mean", 1),
}
ALT_TRANSFORM = {"viirs_log": lambda v: 0.1 * np.log1p(np.maximum(v, 0))}


def land_from_layer(root, layer: str, city: str, shape: tuple[int, int], years: list[int]) -> dict[int, np.ndarray]:
    """{year: (C, Hb, Wb)} with only the built-up channel filled from an alternative layer. A trailing window takes
    the median of years y-w+1 .. y (never later years, so L(t) can't contain L(t+1))."""
    d, prefix, col, w = ALT_LAYERS[layer]
    cache: dict[int, np.ndarray] = {}
    tf = ALT_TRANSFORM.get(layer)

    def year_grid(y: int) -> np.ndarray:
        if y not in cache:
            g = np.full(shape, np.nan)
            f = root / d / f"{prefix}_{y}.csv"
            if f.exists():
                t = pd.read_csv(f)
                t = t[(t.city == city) & t[col].notna()]
                g[t.row, t.col] = tf(t[col].to_numpy()) if tf else t[col]
            cache[y] = g
        return cache[y]

    out = {}
    for y in years:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            stack = np.stack([year_grid(k) for k in range(y - w + 1, y + 1)])
            v = np.nanmedian(stack, 0) if w > 1 else stack[0]
        if w > 1:
            v[np.isfinite(stack).sum(0) < w] = np.nan                    # need the full trailing window
        a = np.full((len(CLASSES), *shape), np.nan)
        a[CLASSES.index("builtup")] = v
        out[y] = a
    return out
