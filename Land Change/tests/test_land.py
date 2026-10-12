import numpy as np
import pandas as pd
import pytest
from rasterio import Affine

from aqf_land.g1 import power_row
from aqf_land.layers import CLASSES, class_fractions, classify, no2_seasonal, pair_table, seasonal_composite
from aqf_land.provenance import version_at
from aqf_land.tnaq import block_mean, city_boxes, load_land_config

cfg = load_land_config()
BANDS = cfg["tnaq"]["s2_bands"]


def test_block_mean_needs_half_the_subpixels():
    a = np.ones((8, 8))
    a[:4, :4][:3] = np.nan                       # 12 of 16 missing -> NaN
    a[:4, 4:][:2] = np.nan                       # 8 of 16 missing -> still a value
    b = block_mean(a, 4, 0.5)
    assert np.isnan(b[0, 0]) and b[0, 1] == 1 and b[1, 1] == 1


def test_city_boxes_are_snapped_and_clipped_to_whole_blocks():
    t = Affine(0.01, 0, 76.22, 0, -0.01, 13.54)  # the TNAQ export grid
    boxes = {b.city: b for b in city_boxes(cfg, t, (547, 414))}
    ch = boxes["Chennai"]
    assert ch.rows == 48 and ch.cols == 40       # export ends at 80.36 E: 12 x 10 blocks
    assert abs(t.c + t.a * ch.col0 - 79.96) < 1e-9 and abs(t.f + t.e * ch.row0 - 13.30) < 1e-9
    assert all(b.rows % 4 == 0 and b.cols % 4 == 0 for b in boxes.values())
    assert boxes["Madurai"].rows == boxes["Madurai"].cols == 48


def _pixel(**r):
    return np.array([r[b] for b in BANDS], np.float32)[:, None, None]


def test_classify_rules_and_order():
    water = _pixel(blue=.05, green=.08, red=.05, nir=.03, swir1=.02, swir2=.01)
    veg = _pixel(blue=.03, green=.06, red=.04, nir=.35, swir1=.18, swir2=.09)
    built = _pixel(blue=.10, green=.12, red=.14, nir=.18, swir1=.24, swir2=.22)
    other = _pixel(blue=.08, green=.10, red=.12, nir=.25, swir1=.22, swir2=.18)
    comp = np.concatenate([water, veg, built, other], 2)
    comp[0, 0, 3] = np.nan
    cls = classify(comp, BANDS)[0]
    assert list(cls) == [0, 1, 2, -1]
    assert CLASSES[2] == "builtup"


def test_class_fractions_sum_to_one_where_valid():
    cls = np.zeros((8, 8), np.int8)
    cls[:4, :4] = 2
    cls[4:, 4:][:2] = -1                         # 8 missing of 16: still valid, fractions over valid only
    cls[4:6, :4] = 1
    f = class_fractions(cls, 4, 0.5)
    assert np.allclose(np.nansum(f, 0), 1)
    assert f[2, 0, 0] == 1 and f[1, 1, 0] == 0.5 and f[0, 1, 1] == 1


def test_seasonal_composite_uses_only_the_season_of_that_year():
    dates = pd.date_range("2020-12-01", "2021-06-30", freq="5D")
    stack = np.zeros((len(dates), 2, 1, 1), np.float32)
    stack[:, :, 0, 0] = np.where(dates.month.isin([1, 2, 3, 4]) & (dates.year == 2021), 1.0, 9.0)[:, None]
    comp = seasonal_composite(stack, dates, 2021, [1, 2, 3, 4], 4)
    assert comp[:, 0, 0].tolist() == [1.0, 1.0]
    s = stack.copy()
    s[np.flatnonzero((dates.year == 2021) & (dates.month < 5))[3:]] = np.nan
    assert np.isnan(seasonal_composite(s, dates, 2021, [1, 2, 3, 4], 4)).all()  # only 3 valid windows


def test_pair_table_dno2_is_log_change_and_dl_in_pp():
    dates = pd.date_range("2019-01-01", "2020-12-31", freq="5D")
    no2 = np.full((len(dates), 1, 4, 4), 2e-5, np.float32)
    no2[dates.year == 2020] = 3e-5
    n = no2_seasonal(no2, dates, [2019, 2020], cfg)
    land = {2019: np.zeros((4, 1, 1)), 2020: np.zeros((4, 1, 1))}
    land[2019][2], land[2020][2] = 0.25, 0.375
    t = pair_table("X", np.array([10.0]), np.array([78.0]), land, n, "builtup")
    assert len(t) == 1 and t.dL.iloc[0] == pytest.approx(12.5)
    assert t.dNO2.iloc[0] == pytest.approx(100 * np.log(1.5), rel=1e-5)


def test_g1_mde_matches_the_formula_and_detects_a_planted_effect():
    rng = np.random.default_rng(0)
    rows = []
    for city in "AB":
        for y0 in (2019, 2020, 2021):
            r, c = np.indices((12, 12))
            dl = rng.gamma(2.0, 3.0, 144)        # real built-up change: growth only
            rows.append(pd.DataFrame({"city": city, "y0": y0, "y1": y0 + 1, "row": r.ravel(), "col": c.ravel(),
                                      "dL": dl, "dNO2": 2.0 * dl + rng.normal(0, 10, 144) + y0 % 7}))
    df = pd.concat(rows, ignore_index=True)
    res = power_row(df, cfg)
    assert res["n"] == 864 and abs(res["rho_neighbour"]) < 0.15
    from scipy.stats import norm

    z = norm.ppf(0.975) + norm.ppf(0.8)
    assert res["dL_noise_sd"] == 0 and res["dL_signal_sd"] == pytest.approx(res["dL_sd_within"])
    assert res["mde_pct_per_pp"] == pytest.approx(z * res["dNO2_resid_sd"] / (res["dL_signal_sd"] * np.sqrt(res["n_eff"])))
    assert res["mde_pct_per_pp"] < 0.5 < 2.0 and res["passes"]


def test_g1_symmetric_dl_counts_as_noise_not_signal():
    rng = np.random.default_rng(1)
    r, c = np.indices((12, 12))
    df = pd.DataFrame({"city": "A", "y0": 2019, "y1": 2020, "row": r.ravel(), "col": c.ravel(),
                       "dL": rng.choice([-12.5, -6.25, 0, 6.25, 12.5], 144), "dNO2": rng.normal(0, 8, 144)})
    res = power_row(df, cfg)
    assert res["dL_signal_sd"] < 0.5 * res["dL_sd_within"] and not res["passes"]


def test_processor_timeline_v1_v2_boundary():
    v = version_at(pd.DatetimeIndex(["2021-04-26", "2021-07-01", "2022-02-01"]))
    assert v[0].startswith("1.") and v[1] == "2.2.0" and v[2] == "2.3.1"


def test_era5_window_scalars_mean_sum_and_wind_speed():
    from aqf_land.era5 import window_scalars, window_start

    times = pd.date_range("2019-01-01", periods=10 * 4, freq="6h")   # 10 days, 4 hours a day
    lat, lon = np.array([14.0, 13.25, 13.0, 12.75]), np.array([79.75, 80.0, 80.25, 80.5])
    shape = (len(times), 4, 4)
    f = {"u10": np.full(shape, 3.0), "v10": np.full(shape, 4.0), "tp": np.full(shape, 0.001)}
    f["u10"][:, 0, 0] = 100.0                                          # outside the box: must be ignored
    df = window_scalars(times, f, lat, lon, (12.8, 79.96, 13.3, 80.44), pad=0.25)
    assert window_start(pd.DatetimeIndex(["2019-01-01", "2019-01-05 23:00", "2019-01-06"])).day.tolist() == [1, 1, 6]
    assert len(df) == 2 and (df.n_hours == 20).all()
    assert np.allclose(df.u10, 3) and np.allclose(df.ws10, 5) and np.allclose(df.tp, 0.02)


def test_land_from_csv_needs_every_year_and_valid_blocks(tmp_path):
    from aqf_land.layers import land_from_csv

    for y in (2019, 2020):
        pd.DataFrame({"city": "X", "row": [0, 1], "col": [0, 1], "year": y, "valid_frac": [0.9, 0.2],
                      "water": 0.0, "vegetation": 0.5, "builtup": [0.3 + 0.1 * (y - 2019), 0.2], "other": 0.2}
                     ).to_csv(tmp_path / f"landfrac_{y}.csv", index=False)
    assert land_from_csv(tmp_path, "X", (2, 2), [2019, 2020, 2021], 0.5) is None
    land = land_from_csv(tmp_path, "X", (2, 2), [2019, 2020], 0.5)
    assert land[2020][2, 0, 0] == pytest.approx(0.4) and np.isnan(land[2020][2, 1, 1]) and np.isnan(land[2019][0, 0, 1])
