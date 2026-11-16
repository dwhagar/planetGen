# tests/test_moist_limit.py

"""
GEN.1579: the moist limit of an air (`physics/moist_limit.py`) replaces the one 375 K the steam limit of GEN.1329 held every body under.

The rule is tested on a small table of its own (what it does with a curve), and the real table is checked against what the sources say of a saturated air.
"""

import numpy as np
import pytest

from planetgen.physics import greenhouse_table, moist_limit


def _grid(curves):
    """A table of two non-condensing pressures by two carbon dioxide pressures at one gravity; `curves` gives the four curves in the order (n0 c0, n0 c1, n1 c0, n1 c1)."""
    temperatures = [300.0, 320.0, 340.0, 360.0, 380.0]
    data = {"temperature_k": temperatures, "non_condensing_pa": [0.0, 1.0e5], "carbon_dioxide_pa": [0.0, 1.0e5], "gravity_m_s2": [9.81],
            "outgoing_w_m2": [[[curves[0], curves[1]], [curves[2], curves[3]]]]}
    return moist_limit.Table(data)


FLAT_EARLY = [100.0, 200.0, 280.0, 281.0, 281.0]            # flattens at 340 K
FLAT_LATE = [100.0, 150.0, 200.0, 250.0, 300.0]           # still rising at the top of the grid
FLAT_FROM_START = [200.0, 200.0, 200.0, 200.0, 200.0]


def test_the_ceiling_is_where_the_outgoing_flux_stops_rising(monkeypatch):
    grid = _grid([FLAT_EARLY] * 4)
    monkeypatch.setattr(moist_limit, "table", lambda: grid)
    limit = moist_limit.limit_k({"N2": 5.0e4}, 9.81)
    assert 330.0 < limit < 360.0


def test_a_curve_still_rising_at_the_top_of_the_grid_has_the_top_as_its_ceiling(monkeypatch):
    grid = _grid([FLAT_LATE] * 4)
    monkeypatch.setattr(moist_limit, "table", lambda: grid)
    assert moist_limit.limit_k({"N2": 5.0e4}, 9.81) == 380.0


def test_a_curve_flat_from_the_start_has_the_bottom_of_the_grid_as_its_ceiling(monkeypatch):
    grid = _grid([FLAT_FROM_START] * 4)
    monkeypatch.setattr(moist_limit, "table", lambda: grid)
    assert moist_limit.limit_k({"N2": 5.0e4}, 9.81) == 300.0


def test_an_air_between_the_grids_airs_has_a_ceiling_between_theirs(monkeypatch):
    grid = _grid([FLAT_EARLY, FLAT_EARLY, FLAT_LATE, FLAT_LATE])
    monkeypatch.setattr(moist_limit, "table", lambda: grid)
    thin = moist_limit.limit_k({}, 9.81)
    thick = moist_limit.limit_k({"N2": 1.0e5}, 9.81)
    middle = moist_limit.limit_k({"N2": 3.0e4}, 9.81)
    assert thin < thick and thin < middle <= thick


def test_only_carbon_dioxide_and_water_are_told_apart_from_the_other_gases(monkeypatch):
    grid = _grid([FLAT_EARLY, FLAT_EARLY, FLAT_LATE, FLAT_LATE])
    monkeypatch.setattr(moist_limit, "table", lambda: grid)
    assert moist_limit.limit_k({"N2": 5.0e4, "O2": 5.0e4}, 9.81) == pytest.approx(moist_limit.limit_k({"N2": 1.0e5}, 9.81))
    assert moist_limit.limit_k({"N2": 1.0e5, "H2O": 2.0e5}, 9.81) == pytest.approx(moist_limit.limit_k({"N2": 1.0e5}, 9.81))


def test_a_hole_in_a_curve_is_read_from_its_neighbours():
    curve = np.array([100.0, np.nan, 300.0])
    assert moist_limit._filled(np.array([300.0, 320.0, 340.0]), curve)[1] == pytest.approx(200.0)


def test_the_gravity_must_be_positive():
    with pytest.raises(ValueError):
        moist_limit.limit_k({"N2": 1.0e5}, 0.0)


def test_the_table_file_was_built_from_the_catalogs_data():
    """`table()` refuses a file built from other k-tables; the file in the repository must be the catalog's."""
    import json
    with open(moist_limit.DATA_FILE, encoding="utf-8") as handle:
        assert json.load(handle)["data_digest"] == greenhouse_table.data_digest()


def test_a_thin_nitrogen_air_flattens_where_the_saturated_models_curve_does():
    """The curve of the saturated clear-sky air rises about 1.4 W/m2 per K near 290 K and is flat by 330 K (the published limit, 282 to 296 W/m2 [R1607], is the flat of it)."""
    assert 300.0 < moist_limit.limit_k({"N2": 1.0e5, "CO2": 40.0}, 9.81) < 350.0


def test_the_ceiling_does_not_fall_when_the_carbon_dioxide_rises():
    previous = 0.0
    for carbon_dioxide in (0.0, 1.0e4, 1.0e5, 1.0e6, 1.0e7):
        limit = moist_limit.limit_k({"N2": 1.0e5, "CO2": carbon_dioxide}, 9.81)
        assert limit >= previous - 1.0
        previous = limit


def test_a_thick_carbon_dioxide_air_is_allowed_above_400_k():
    """GEN.1568: the air of such a body held a stable moist surface above the old 375 K [R4773]."""
    assert moist_limit.limit_k({"CO2": 1.0e7}, 9.81) > 400.0
    assert moist_limit.limit_k({"N2": 1.0e5}, 9.81) < 400.0
