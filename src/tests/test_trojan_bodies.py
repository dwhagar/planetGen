# tests/test_trojan_bodies.py

"""GEN.1151 pure parts: which whole bodies may sit at L4 and L5 of a host and a secondary, and the formation gate."""

import pytest

from planetgen.physics import trojan_bodies as tb
from planetgen.physics import trojans

SUN_EARTH, JUPITER_EARTH = 332946.0, 317.8


def test_a_sun_and_a_jupiter_admit_class_j_and_class_m():
    result = tb.classes_allowed(SUN_EARTH, JUPITER_EARTH)
    assert result["class_j"] and result["class_m"]
    assert result["room_earth"] == pytest.approx(0.0385209 / (1 - 0.0385209) * SUN_EARTH - JUPITER_EARTH, rel=1e-3)


def test_a_dwarf_star_limits_class_j_to_a_few_jupiters():
    host = 0.1 * SUN_EARTH
    assert tb.classes_allowed(host, JUPITER_EARTH)["room_earth"] / JUPITER_EARTH < 4.0
    assert tb.classes_allowed(host, 4.15 * JUPITER_EARTH)["class_j"] is False


def test_two_stars_admit_nothing_and_a_brown_dwarf_admits_class_m():
    assert tb.classes_allowed(SUN_EARTH, 0.5 * SUN_EARTH) == {"room_earth": 0.0, "class_j": False, "class_m": False}
    brown = 0.03 * SUN_EARTH
    result = tb.classes_allowed(SUN_EARTH, brown)
    assert result["class_m"] and result["class_j"]


def test_a_planet_holds_a_moon_and_a_class_m_sized_trojan_only_when_big_enough():
    assert tb.classes_allowed(JUPITER_EARTH, 1.0)["class_m"]
    assert not tb.classes_allowed(3.0, 0.1)["class_m"]
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


def test_the_formation_gate_needs_solids_and_a_long_enough_disc():
    assert tb.can_form(1.0, solids_in_zone_earth=12.0, disc_lifetime_myr=3.0)
    assert not tb.can_form(1.0, solids_in_zone_earth=5.0, disc_lifetime_myr=3.0)
    assert not tb.can_form(1.0, solids_in_zone_earth=12.0, disc_lifetime_myr=0.5)
    assert not tb.can_form(0.0, 12.0, 3.0)
