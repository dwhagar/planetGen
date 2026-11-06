# tests/test_irregular_moons.py

"""
GEN.1083: the irregular moons of a planet (`generation/irregular_moons.py`): the second tier of moons above the regular ones, caught from the outer disk and the planet's feeding zone in the capture windows,
and the captured-origin weight of the regular moons (`moon_composition.draw_origin`).
"""

import types

import pytest

from planetgen.generation import irregular_moons
from planetgen.generation.planet import Planet
from planetgen.physics import capture_windows, irregular_zone, moon_composition
from planetgen.util import draw

SUN_KG = 1.989e30
JUPITER_KG = 1.898e27
EARTH_KG = 5.972e24


def _planet(distance_au, mass_kg, radius_km, envelope=0.9, moons=(), disk_myr=3.0):
    return types.SimpleNamespace(distance=distance_au, mass=mass_kg, radius=radius_km, envelope_mass_fraction=envelope, moons=list(moons), moon_disk_lifetime_myr=disk_myr)


def _source(reservoir=30.0, encounters=130.0, age_gy=4.6):
    return irregular_moons.Source(reservoir, 30.0, 50.0, encounters, 3.828e26, 5.0e-5, SUN_KG, age_gy)


def test_a_planet_with_nothing_to_catch_from_has_none():
    assert irregular_moons.for_planet(_planet(5.2, JUPITER_KG, 69911.0), _source(reservoir=0.0)) is irregular_moons.NO_IRREGULARS


def test_a_jupiter_after_an_instability_keeps_irregulars_that_stay_in_its_zone():
    planet = _planet(5.2, JUPITER_KG, 69911.0)
    caught = irregular_moons.for_planet(planet, _source())
    assert len(caught.moons) >= 3
    zone = irregular_zone.zone(5.2, JUPITER_KG, SUN_KG)
    masses = [moon.mass_kg for moon in caught.moons]
    assert masses == sorted(masses, reverse=True) and [moon.rank for moon in caught.moons] == list(range(len(caught.moons)))
    for moon in caught.moons:
        assert irregular_zone.verdict(zone, moon.orbit.semi_major_axis_au, moon.orbit.eccentricity, moon.orbit.inclination_deg) is None
        assert moon.window in ("gas", "encounter", "exchange") and moon.radius_km > 0.0
    assert caught.remainder_kg > 0.0 and caught.above_floor >= len(caught.moons)


def test_the_same_planet_gets_the_same_irregulars():
    planet = _planet(5.2, JUPITER_KG, 69911.0)
    assert irregular_moons.for_planet(planet, _source()) == irregular_moons.for_planet(planet, _source())


def test_an_earth_at_one_au_catches_nothing():
    caught = irregular_moons.for_planet(_planet(1.0, EARTH_KG, 6371.0, envelope=0.0), _source())
    assert not caught.moons


def test_a_hot_jupiter_has_no_zone_for_them():
    caught = irregular_moons.for_planet(_planet(0.04, JUPITER_KG, 80000.0), _source())
    assert not caught.moons


def test_a_system_with_no_instability_has_fewer_than_one_that_had_one():
    planet = _planet(5.2, JUPITER_KG, 69911.0)
    quiet = irregular_moons.for_planet(planet, _source(encounters=0.0))
    loud = irregular_moons.for_planet(planet, _source())
    assert sum(moon.mass_kg for moon in quiet.moons) <= sum(moon.mass_kg for moon in loud.moons)


def test_the_source_is_the_outer_disks_of_the_bodies_and_the_windows_of_the_system():
    disk = types.SimpleNamespace(body_type='a', origin="outer-disk", primordial_mass_earth=12.0, lower_limit=30.0, upper_limit=45.0)
    other = types.SimpleNamespace(body_type='a', origin="inner", primordial_mass_earth=99.0, lower_limit=2.0, upper_limit=3.0)
    planet = types.SimpleNamespace(body_type='g')
    star = types.SimpleNamespace(luminosity=3.828e26, mass=SUN_KG, age=4.6)
    source = irregular_moons.source_of([planet, other, disk], star, capture_windows.Windows(1, 130.0), 5.0e-5)
    assert (source.reservoir_mass_earth, source.inner_au, source.outer_au, source.encounters) == (12.0, 30.0, 45.0, 130.0)
    assert irregular_moons.source_of([planet], star, capture_windows.NO_WINDOWS, 5.0e-5).reservoir_mass_earth == 0.0


def test_a_heavy_moon_is_rarely_a_captured_one_and_a_light_one_often():
    with draw.bound(7):
        heavy = sum(moon_composition.draw_origin(moon_composition.GIANT, 1.0e22) == "captured" for _ in range(400))
        light = sum(moon_composition.draw_origin(moon_composition.GIANT, 1.0e15) == "captured" for _ in range(400))
    assert heavy < 20
    assert 150 < light < 250
