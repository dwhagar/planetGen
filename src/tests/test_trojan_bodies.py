# tests/test_trojan_bodies.py

"""GEN.1151 pure parts: which whole bodies may sit at L4 and L5 of a host and a secondary, and the formation gate."""

import pytest

from planetgen.physics import trojan_bodies as tb
from planetgen.physics import trojans

SUN_EARTH, JUPITER_EARTH = 332946.0, 317.8


def test_a_sun_and_a_jupiter_hold_a_giant_and_a_rocky_body():
    room = tb.room_earth(SUN_EARTH, JUPITER_EARTH)
    assert room == pytest.approx(0.0385209 / (1 - 0.0385209) * SUN_EARTH - JUPITER_EARTH, rel=1e-3)
    assert tb.fits(SUN_EARTH, JUPITER_EARTH, 35.0) and tb.fits(SUN_EARTH, JUPITER_EARTH, 1.0)


def test_a_dwarf_star_limits_the_trojan_to_a_few_jupiters():
    host = 0.1 * SUN_EARTH
    assert tb.room_earth(host, JUPITER_EARTH) / JUPITER_EARTH < 4.0
    assert not tb.fits(host, 4.15 * JUPITER_EARTH, 35.0)


def test_two_stars_hold_nothing_and_a_brown_dwarf_holds_a_giant():
    assert tb.room_earth(SUN_EARTH, 0.5 * SUN_EARTH) == 0.0
    assert tb.fits(SUN_EARTH, 0.03 * SUN_EARTH, 35.0)


def test_a_planet_holds_a_moon_sized_trojan_only_when_big_enough():
    assert tb.fits(JUPITER_EARTH, 1.0, 1.0)
    assert not tb.fits(3.0, 0.1, 1.0)
    assert tb.room_earth(JUPITER_EARTH, 1.0) == pytest.approx(11.7, abs=0.3)


def test_the_cases():
    assert tb.case("star", "planet", 1.0, 0.001) == "star_giant"
    assert tb.case("star", "star", 1.0, 0.03) == "star_substellar"
    assert tb.case("star", "star", 1.0, 0.5) == "star_star"
    assert tb.case("planet", "moon", 1.0, 0.001) == "planet_moon"
    assert tb.case("moon", "moon", 1.0, 0.001) == "moon_moon"
    with pytest.raises(trojans.TrojanError):
        tb.case("moon", "planet", 1.0, 0.1)


def test_fits_follows_the_room_and_the_orbit_limits():
    assert tb.fits(SUN_EARTH, JUPITER_EARTH, 1.0, 0.02, 3.0)
    assert not tb.fits(SUN_EARTH, JUPITER_EARTH, 1.0, 0.6, 3.0)
    assert not tb.fits(SUN_EARTH, 0.5 * SUN_EARTH, 1.0)


def test_the_disk_life_is_a_ramp_from_the_foot_to_the_full_growth_time():
    low, high = tb.DISK_LIFETIME_RAMP_MYR
    assert tb.disk_life_factor(low - 0.1) == 0.0 and tb.disk_life_factor(high + 5.0) == 1.0
    assert tb.disk_life_factor((low + high) / 2.0) == pytest.approx(0.5)


def test_the_zone_holds_the_solids_of_two_tadpoles_and_the_ceiling_is_its_share_less_the_room():
    # A Jupiter at 5.2 au of the Sun in a nebula of 10 g/cm^2 of solids: about 4 Earth masses in the zone
    solids = tb.zone_solids_earth(10.0, 5.2, JUPITER_EARTH, SUN_EARTH)
    assert 2.0 < solids < 8.0
    assert tb.zone_solids_earth(20.0, 5.2, JUPITER_EARTH, SUN_EARTH) == pytest.approx(2.0 * solids)
    full = tb.ceiling_earth(solids, 3.0, SUN_EARTH, JUPITER_EARTH)
    assert full == pytest.approx(tb.FORMATION_EFFICIENCY * solids)
    assert tb.ceiling_earth(solids, 1.65, SUN_EARTH, JUPITER_EARTH) == pytest.approx(0.5 * full)
    assert tb.ceiling_earth(solids, 3.0, SUN_EARTH, 0.5 * SUN_EARTH) == 0.0                # a pair over the limit holds nothing


def test_the_mass_is_drawn_between_the_lightest_planet_and_the_ceiling_or_not_at_all():
    assert tb.draw_mass_earth(0.05, 0.1, 0.5) is None
    assert tb.draw_mass_earth(2.0, 0.1, 0.0) == pytest.approx(0.1) and tb.draw_mass_earth(2.0, 0.1, 1.0) == pytest.approx(2.0)
    assert tb.draw_mass_earth(10.0, 0.1, 0.5) == pytest.approx(1.0)                           # log-uniform: the middle is the geometric mean
