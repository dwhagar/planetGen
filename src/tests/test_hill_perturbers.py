"""
GEN.1449, slice 1: a planet's Hill sphere against every body that pulls on it (physics/hill_perturbers.py), and the moons it keeps (physics/planets.py `moon_orbit_bounds_km`,
`moon_stable_limit_km`, `drop_unstable_moons`; physics/irregular_zone.py; physics/moon_history.py).
"""

import math
import types

import pytest

from planetgen.physics import constants, hill_perturbers, orbits, planets, tide_sources


class _Body:
    def __init__(self, mass, radius, distance=None, eccentricity=0.0, star=None, moons=(), name="b"):
        self.mass, self.radius, self.distance, self.eccentricity, self.star, self.moons, self.name = mass, radius, distance, eccentricity, star, list(moons), name
        self.is_moon = False
        self.scale_height = 8.0
        self.hill_radius = None


def _planet_with_moon(fraction):
    """The Earth at 1 AU round the Sun with one moon at `fraction` of its Hill radius."""
    star = _Body(constants.SOLAR_MASS_TO_KG, 6.9e5, name="star")
    earth = _Body(constants.EARTH_MASS_TO_KG, 6371.0, 1.0, 0.0, star, name="earth")
    earth.primary_mass_kg = star.mass
    planets.update_hill_sphere(earth)
    moon = _Body(7.35e22, 1737.0, fraction * earth.hill_radius / constants.AU_TO_M * 1000.0, 0.0, star, name="moon")
    moon.is_moon = True
    earth.moons = [moon]
    return star, earth, moon


def test_a_planet_nothing_but_its_host_pulls_on_keeps_its_whole_hill_sphere():
    star, earth, moon = _planet_with_moon(0.3)
    assert hill_perturbers.tidal_ratio(earth) == 0.0
    assert hill_perturbers.hill_scale(earth) == 1.0
    assert hill_perturbers.tidal_hill_radius_km(earth) == earth.hill_radius
    earth._tide_siblings = (earth,)       # itself is not a neighbour
    assert hill_perturbers.hill_scale(earth) == 1.0


def test_a_moon_or_a_body_with_no_orbit_has_no_ratio():
    star, earth, moon = _planet_with_moon(0.3)
    assert hill_perturbers.tidal_ratio(moon) == 0.0
    assert hill_perturbers.tidal_ratio(_Body(1.0e24, 6000.0, None, 0.0, star)) == 0.0
    assert hill_perturbers.tidal_hill_radius_km(_Body(1.0e24, 6000.0)) is None


def test_the_ratio_of_one_sibling_is_twice_its_mass_over_the_host_field_by_the_mean_of_one_over_d_cubed():
    star, earth, moon = _planet_with_moon(0.3)
    neighbour = _Body(10.0 * constants.EARTH_MASS_TO_KG, 12000.0, 1.2, 0.0, star, name="neighbour")
    earth._tide_siblings = (neighbour,)
    g3 = tide_sources.orbit_mean_inverse_cube(1.0, 1.2, 0.0) / constants.AU_TO_M ** 3
    host_field = 3.0 * (star.mass + earth.mass) / (constants.AU_TO_M) ** 3
    expected = 2.0 * neighbour.mass * g3 / host_field
    assert hill_perturbers.tidal_ratio(earth) == pytest.approx(expected, rel=1e-9)
    assert hill_perturbers.hill_scale(earth) == pytest.approx((1.0 + expected) ** (-1.0 / 3.0), rel=1e-12)
    assert 0.0 < expected < 1.0e-2     # a neighbour at 0.2 AU of a few Earths is a small field against the Sun's


def test_a_massive_sibling_that_passes_close_shrinks_the_sphere_and_the_moons_limit():
    star, earth, moon = _planet_with_moon(0.3)
    alone = planets.moon_stable_limit_km(earth)
    earth._tide_siblings = (_Body(3000.0 * constants.EARTH_MASS_TO_KG, 70000.0, 1.15, 0.2, star, name="giant"),)
    assert hill_perturbers.hill_scale(earth) < 0.8
    assert planets.moon_stable_limit_km(earth) < 0.8 * alone
    assert planets.moon_orbit_bounds_km(earth)[1] < 0.8 * planets.moon_orbit_bounds_km(_planet_with_moon(0.3)[1])[1]


def test_a_far_sibling_leaves_the_moons_limit_within_a_hair_of_the_hosts_alone():
    star, earth, moon = _planet_with_moon(0.3)
    alone = planets.moon_stable_limit_km(earth)
    earth._tide_siblings = (_Body(317.8 * constants.EARTH_MASS_TO_KG, 69911.0, 5.2, 0.048, star, name="jupiter"),)
    assert 0.99 * alone < planets.moon_stable_limit_km(earth) < alone


def test_a_far_moon_is_lost_to_a_close_massive_sibling_and_a_near_one_stays():
    star, earth, far = _planet_with_moon(0.45)
    near = _Body(7.0e21, 1200.0, 0.1 * earth.hill_radius / constants.AU_TO_M * 1000.0, 0.0, star, name="near")
    near.is_moon = True
    earth.moons = [near, far]
    assert planets.drop_unstable_moons(earth) == 0       # alone, both are held
    earth.moons = [near, far]
    earth._tide_siblings = (_Body(3000.0 * constants.EARTH_MASS_TO_KG, 70000.0, 1.15, 0.2, star, name="giant"),)
    assert planets.drop_unstable_moons(earth) == 1
    assert earth.moons == [near]


def test_a_companion_star_shrinks_the_sphere_by_its_tide():
    star, earth, moon = _planet_with_moon(0.3)
    companion = types.SimpleNamespace(star=_Body(0.8 * constants.SOLAR_MASS_TO_KG, 5.0e5, name="b"), axis_au=8.0, eccentricity=0.3, member_no=1)
    star._companions = (companion,)
    assert 0.5 < hill_perturbers.hill_scale(earth) < 1.0
    star._companions = ()
    assert hill_perturbers.hill_scale(earth) == 1.0


def test_the_irregular_zone_scales_with_the_share_of_the_hill_radius_kept():
    from planetgen.physics import irregular_zone
    star = constants.SOLAR_MASS_TO_KG
    full = irregular_zone.zone(5.2, 1.9e27, star, 0.01)
    assert irregular_zone.zone(5.2, 1.9e27, star, 0.01, 1.0) == full
    half = irregular_zone.zone(5.2, 1.9e27, star, 0.01, 0.5)
    assert half.hill_radius_au == pytest.approx(0.5 * full.hill_radius_au)
    assert half.prograde_limit_au == pytest.approx(0.5 * full.prograde_limit_au)
    assert orbits.calculate_hill_sphere(5.2 * constants.AU_TO_M, 1.9e27, star) / constants.AU_TO_M == pytest.approx(full.hill_radius_au)
