import numpy as np
import pandas as pd
import pytest

from aqf.anomaly import fit_climatology, harmonics
from aqf_forecast.metrics import block_bootstrap_skill, crps_from_sorted, target_obs
from aqf_forecast.problem import last_valid


def test_climatology_recovers_a_known_seasonal_cycle_with_gaps():
    dates = pd.date_range("2019-01-01", periods=292, freq="5D")
    X = harmonics(dates, 2)
    true = np.array([5.0, 1.0, -0.5, 0.3, 0.2])
    y = (X @ true)[:, None, None, None].repeat(2, 2)
    y[::7] = np.nan                                    # missing windows must not matter
    clim = fit_climatology(y, dates, k=2, ridge=0.0)
    assert np.allclose(clim.coef[:, 0, 0, 0], true, atol=1e-8)
    assert np.allclose(clim(dates)[:, 0, 0, 0], X @ true)


def test_last_valid_uses_only_the_past():
    x = np.array([1.0, np.nan, np.nan, 4.0, np.nan, np.nan, np.nan, np.nan])
    v, age = last_valid(x, lookback=2)
    assert v[0] == 1 and age[0] == 0
    assert v[2] == 1 and age[2] == 2                   # 2 windows old, inside the lookback
    assert v[3] == 4 and age[3] == 0                   # never takes a future value
    assert np.isnan(v[6]) and age[6] == -1             # 3 windows old: outside the lookback


def test_target_obs_alignment():
    vals = np.arange(10.0)[:, None]
    t = target_obs(vals, 3)
    assert t[0, 0, 0] == 1 and t[0, 2, 0] == 3 and np.isnan(t[8, 1, 0])


def test_crps_matches_brute_force():
    rng = np.random.default_rng(0)
    r = np.sort(rng.normal(size=300))
    err = rng.normal(size=50) * 2
    brute = np.abs(r[None] - err[:, None]).mean(1) - 0.5 * np.abs(r[:, None] - r[None]).mean()
    assert np.allclose(crps_from_sorted(err, r), brute)


def test_bootstrap_ci_contains_point_skill():
    rng = np.random.default_rng(1)
    ref = rng.uniform(1, 2, 200)
    m = ref * 0.8
    lo, hi = block_bootstrap_skill(m, ref, block=6, n=500)
    assert lo <= 0.2 <= hi + 1e-12


def test_no_lookahead_in_baselines():
    """Changing data after the origin must not change any baseline forecast issued at that origin."""
    from aqf_forecast.baselines import MODELS, predict
    from aqf_forecast.problem import fit, load_problem

    try:
        p = load_problem()
    except FileNotFoundError:
        pytest.skip("raw dataset not present")
    t0 = 100
    f = fit(p, p.windows("train"))
    before = {m: predict(m, p, f, 5)[t0] for m in MODELS}
    p.values = p.values.copy()
    p.values[t0 + 1:] += 1.0                            # corrupt everything after the origin
    f.anom = p.values - f.clim_ext[: p.T]               # same fitted parameters, corrupted observations
    for m in MODELS:
        after = predict(m, p, f, 5)[t0]
        assert np.allclose(before[m], after, equal_nan=True), m


def test_trailing_mean_uses_only_the_past():
    from aqf_forecast.features import _trailing_nanmean

    x = np.array([1.0, np.nan, 3.0, 5.0, 100.0])
    m = _trailing_nanmean(x, 3)
    assert np.isclose(m[2], 2.0) and np.isclose(m[3], 4.0) and np.isclose(m[1], 1.0)


def test_no_lookahead_in_f2_features():
    """Every feature of an (origin t0, lead) row must be unchanged when all windows after t0 are altered."""
    from aqf_forecast.features import build, rows_for
    from aqf_forecast.problem import fit, load_problem

    try:
        p = load_problem()
    except FileNotFoundError:
        pytest.skip("raw dataset not present")
    t0 = 150
    train = p.windows("train")
    f = fit(p, train)
    fs = build(p, f, 0, None, train)
    X0, _, _ = rows_for(fs, np.array([t0, t0]), np.array([1, 30]))
    f.anom = f.anom.copy()
    f.anom[t0 + 1:] = np.random.default_rng(0).normal(size=f.anom[t0 + 1:].shape)
    fs2 = build(p, f, 0, None, train, sd=fs.sd)           # same fitted scale, corrupted future
    X1, _, _ = rows_for(fs2, np.array([t0, t0]), np.array([1, 30]))
    assert np.array_equal(X0, X1, equal_nan=True)
