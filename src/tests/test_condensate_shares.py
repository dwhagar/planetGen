# tests/test_condensate_shares.py

"""GEN.950: the condensates of Lodders 2003 Tables 10 and 11 by disk temperature and the water share of the solids (`physics/condensate_shares.py`)."""

import pytest

from planetgen.physics import condensate_shares as shares
from planetgen.physics import constants as physical_constants

SUN_W = physical_constants.SOLAR_LUMINOSITY


def test_nothing_but_rock_is_whole_above_the_water_line():
    assert shares.condensed_masses(300.0) == {"rock": 0.489}
    assert shares.water_share_of_solids(300.0) == 0.0
    assert shares.ice_over_rock(300.0) == 0.0


def test_the_water_share_steps_down_at_each_ice_line_as_table_11_prints():
    """Lodders 2003 Table 11: rock 0.489, water 0.571, ammonia 0.097, methane 0.330 (water over rock 1.17, all ices over rock 2.04)."""
    assert shares.water_share_of_solids(150.0) == pytest.approx(0.571 / (0.489 + 0.571))
    assert shares.water_share_of_solids(100.0) == pytest.approx(0.571 / (0.489 + 0.571 + 0.097))
    assert shares.water_share_of_solids(30.0) == pytest.approx(0.571 / 1.487)
    assert shares.ice_over_rock(150.0) == pytest.approx(1.17, abs=0.005)
    assert shares.ice_over_rock(30.0) == pytest.approx(0.998 / 0.489, abs=0.005)


def test_a_line_belongs_to_the_colder_side():
    assert "H2O" not in shares.condensed_masses(182.0)
    assert "H2O" in shares.condensed_masses(181.9)
    assert "NH3" in shares.condensed_masses(130.0) and "NH3" not in shares.condensed_masses(131.0)


def test_the_boost_of_the_solids_at_the_snow_line_is_the_water_over_rock_plus_one():
    assert shares.solid_boost(150.0) == pytest.approx(1.0 + 0.571 / 0.489)
    assert shares.solid_boost(150.0) < 30.0 / 7.0


def test_oxygen_after_the_rock_makes_the_water():
    """22.75 percent of the Sun's oxygen is in rock; a star with that much oxygen per rock atom or less has no water ice."""
    assert shares.condensed_masses(150.0, oxygen=shares.ROCK_OXYGEN_SHARE)["H2O"] == 0.0
    assert shares.condensed_masses(150.0, oxygen=0.0)["H2O"] == 0.0
    assert shares.condensed_masses(150.0, oxygen=2.0)["H2O"] == pytest.approx(0.571 * (2.0 - 0.2275) / 0.7725)
    assert shares.water_share_of_solids(150.0, oxygen=1.5) > shares.water_share_of_solids(150.0) > shares.water_share_of_solids(150.0, oxygen=0.6)


def test_ammonia_and_methane_follow_their_own_element():
    assert shares.condensed_masses(30.0, carbon=2.0)["CH4"] == pytest.approx(0.660)
    assert shares.condensed_masses(100.0, nitrogen=0.5)["NH3"] == pytest.approx(0.0485)


def test_a_negative_ratio_or_a_cold_zero_is_refused():
    with pytest.raises(ValueError):
        shares.condensed_masses(100.0, oxygen=-0.1)
    with pytest.raises(ValueError):
        shares.condensed_masses(0.0)


def test_the_blackbody_temperature_falls_as_the_square_root_of_the_distance():
    assert shares.disk_temperature_k(1.0, SUN_W) == pytest.approx(278.6, abs=0.5)
    assert shares.disk_temperature_k(4.0, SUN_W) == pytest.approx(shares.disk_temperature_k(1.0, SUN_W) / 2.0)


def test_the_water_line_of_the_sun_is_where_a_blackbody_is_182_k():
    line = shares.ice_line_au(SUN_W)
    assert line == pytest.approx(2.34, abs=0.02)
    assert shares.disk_temperature_k(line, SUN_W) == pytest.approx(182.0)
    assert shares.ice_line_au(SUN_W, "CH4") > shares.ice_line_au(SUN_W, "NH3") > line
    assert shares.ice_line_au(4.0 * SUN_W) == pytest.approx(2.0 * line)


def test_the_share_by_distance_is_zero_inside_the_line_and_never_rises_outward():
    line = shares.ice_line_au(SUN_W)

    def share_at(distance_au):
        return shares.water_share_of_solids(shares.disk_temperature_k(distance_au, SUN_W))

    assert share_at(0.9 * line) == 0.0
    outward = [share_at(line * 1.05 ** step) for step in range(1, 200)]
    assert outward[0] == pytest.approx(0.539, abs=0.001)
    assert all(a >= b for a, b in zip(outward, outward[1:]))
    assert outward[-1] == pytest.approx(0.571 / 1.487, abs=0.001)


def test_an_unknown_ice_is_refused():
    with pytest.raises(ValueError):
        shares.ice_line_au(SUN_W, "CO2")
