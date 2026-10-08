import numpy as np
import pytest

from aqf.data import load_config, load_manifest
from aqf.series import window_table

cfg = load_config()
pytestmark = pytest.mark.skipif(not cfg["paths"]["manifest"].exists(), reason="raw dataset not present")


def test_manifest_is_a_contiguous_5_day_grid():
    m = load_manifest(cfg)
    assert len(m) == 438 and m.window_id.is_unique
    assert (m.start_date.diff().dropna().dt.days == 5).all()
    assert set(m.split) == {"train", "val", "test"}


def test_s5p_cube_grid_and_footprint():
    from aqf.data import build_s5p_cube

    cube = build_s5p_cube(cfg)
    assert cube.values.shape == (438, 3, 13, 13)
    fp = cube.footprint
    assert not fp[12].any() and not fp[:, 0].any()  # always-nodata edge row and column
    assert fp[:12, 1:].all()                          # everything else is inside the area: 12 x 12
    assert np.all(np.diff(cube.lat) < 0) and np.all(np.diff(cube.lon) > 0)


def test_window_table_flags():
    from aqf.data import build_s5p_cube

    wt = window_table(build_s5p_cube(cfg), cfg)
    assert len(wt) == 438 * 3
    assert not wt.loc[wt.all_zero, "usable"].any()
    assert wt.loc[wt.usable, "mean"].notna().all()


def test_zeros_are_missing_in_analysis_values():
    from aqf.data import build_s5p_cube
    from aqf.series import analysis_values

    cube = build_s5p_cube(cfg)
    v = analysis_values(cube, cfg)
    assert (cube.values == 0).any()          # the raw data does contain zeros...
    assert not (v == 0).any()                # ...but no analysis ever sees them as values
    assert np.isnan(v[:, :, 12]).all() and np.isnan(v[:, :, :, 0]).all()
