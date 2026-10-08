import numpy as np

from aqf.stats import fdr_bh, seasonal_trend


def _monthly(trend_per_year, n_years=6, noise=0.0, seed=0):
    rng = np.random.default_rng(seed)
    years = np.repeat(np.arange(2019, 2019 + n_years), 12)
    months = np.tile(np.arange(1, 13), n_years)
    seasonal = 10 * np.sin(2 * np.pi * months / 12)  # strong seasonal cycle the test must ignore
    x = 100 + seasonal + trend_per_year * (years - 2019) + noise * rng.standard_normal(len(years))
    return x, years, months


def test_sen_slope_recovers_exact_linear_trend():
    x, y, m = _monthly(2.5)
    r = seasonal_trend(x, y, m)
    assert np.isclose(r.slope, 2.5)
    assert r.p < 1e-6 and r.direction == "increasing"


def test_no_trend_in_pure_noise_is_not_significant():
    x, y, m = _monthly(0.0, noise=3.0, seed=1)
    assert seasonal_trend(x, y, m).p > 0.05


def test_missing_values_are_tolerated():
    x, y, m = _monthly(-1.0, noise=0.1)
    x[::5] = np.nan
    r = seasonal_trend(x, y, m)
    assert r.slope < 0 and r.p < 0.001


def test_fdr_bh_matches_hand_computation():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, np.nan])
    # thresholds 0.05*k/6: .0083 .0167 .025 .033 .0417 .05 -> only the first two pass
    assert fdr_bh(p).tolist() == [True, True, False, False, False, False, False]
