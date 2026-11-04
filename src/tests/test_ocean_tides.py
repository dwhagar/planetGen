# tests/test_ocean_tides.py

"""
GEN.1034: the open-ocean equilibrium tide (`planetgen.physics.ocean_tides`): the Earth's published equilibrium range of the Moon and the Sun, spring and neap, a world with no
moon, a world locked to its star, the eccentric extremes, the loading factor and its refusal, and an ice-covered ocean.
"""

import math

import pytest

from planetgen.physics import constants as pc
from planetgen.physics import ocean_tides
from planetgen.physics.ocean_tides import Liquid, OceanTideError, Raiser

EARTH_MASS = 5.972e24
EARTH_RADIUS = 6.371e6
EARTH_SPIN = 2.0 * math.pi / 86164.0
MOON = Raiser.of("Moon", 7.346e22, 3.844e8, EARTH_MASS, 0.0549)
SUN = Raiser.of("Sun", 1.989e30, pc.AU_TO_M, EARTH_MASS, 0.0167)


def _earth(raisers=(MOON, SUN), spin=EARTH_SPIN, **kwargs):
    kwargs.setdefault("k2", 0.30)
    kwargs.setdefault("h2", 0.61)
    return ocean_tides.ocean_tide(EARTH_MASS, EARTH_RADIUS, spin, raisers, **kwargs)


def test_the_equilibrium_tide_is_the_published_half_metre_of_the_moon_and_quarter_metre_of_the_sun():
    assert ocean_tides.equilibrium_range_m(7.346e22, 3.844e8, EARTH_MASS, EARTH_RADIUS) == pytest.approx(0.54, abs=0.02)
    assert ocean_tides.equilibrium_range_m(1.989e30, pc.AU_TO_M, EARTH_MASS, EARTH_RADIUS) == pytest.approx(0.25, abs=0.02)
    tide = _earth(k2=0.0, h2=0.0, liquid=Liquid("weightless", 0.0))
    assert tide.spring_range_m == pytest.approx(0.79, abs=0.03)
    assert tide.neap_range_m == pytest.approx(0.29, abs=0.03)


def test_the_floor_and_the_weight_of_the_sea_move_the_range_by_their_two_factors():
    tide = _earth()
    assert tide.love_factor == pytest.approx(0.69)
    assert tide.loading_factor == pytest.approx(1.0 / (1.0 - 0.6 * 1000.0 / 5514.0), rel=1e-3)
    assert tide.spring_range_m == pytest.approx(_earth(k2=0.0, h2=0.0, liquid=Liquid("weightless", 0.0)).spring_range_m * 0.69 * tide.loading_factor)
    assert ocean_tides.love_factor(0.3) == pytest.approx(0.8)


def test_an_ocean_world_with_no_moon_shows_only_the_stars_tide():
    tide = _earth(raisers=(SUN,))
    assert [t.name for t in tide.raisers] == ["Sun"]
    assert tide.spring_range_m == tide.neap_range_m > 0.0


def test_a_world_locked_to_its_star_has_no_tide_from_it():
    locked = _earth(raisers=(SUN,), spin=SUN.orbital_rate_rad_s)
    assert not locked.raisers[0].oscillates and locked.spring_range_m == locked.neap_range_m == 0.0
    mixed = _earth(spin=SUN.orbital_rate_rad_s)
    assert mixed.spring_range_m == pytest.approx(mixed.raisers[0].range_m) and mixed.neap_range_m == mixed.spring_range_m
    backwards = Raiser.of("x", 1.0, pc.AU_TO_M, EARTH_MASS, retrograde=True)
    assert not ocean_tides.oscillates(backwards.orbital_rate_rad_s, backwards)
    assert ocean_tides.oscillates(-backwards.orbital_rate_rad_s, backwards)


def test_a_planet_that_turns_with_its_moon_sees_no_tide_from_it():
    charon_like = Raiser("Charon", 1.0e21, 1.9e7, 0.0, 1.14e-5)
    assert not ocean_tides.oscillates(1.14e-5 * 1.0005, charon_like)
    assert ocean_tides.oscillates(1.14e-5 * 1.01, charon_like)


def test_an_eccentric_orbit_varies_the_range_as_the_cube_of_the_distance():
    tide = _earth(raisers=(MOON,))
    moon = tide.raisers[0]
    assert moon.range_closest_m == pytest.approx(moon.range_m / (1.0 - 0.0549) ** 3)
    assert moon.range_farthest_m == pytest.approx(moon.range_m / (1.0 + 0.0549) ** 3)
    assert tide.spring_range_max_m == pytest.approx(moon.range_closest_m)


def test_the_neap_range_is_the_largest_less_the_rest_and_never_below_zero():
    assert ocean_tides.ranges([0.54, 0.25]) == pytest.approx((0.79, 0.29))
    assert ocean_tides.ranges([0.5, 0.4, 0.3]) == pytest.approx((1.2, 0.0))
    assert ocean_tides.ranges([0.5]) == (0.5, 0.5)
    assert ocean_tides.ranges([]) == (0.0, 0.0)
    tide = _earth()
    assert 0.0 < tide.neap_range_min_m < tide.neap_range_m < tide.spring_range_m < tide.spring_range_max_m


def test_a_liquid_too_heavy_for_a_thin_layer_is_refused_and_a_light_one_is_not():
    with pytest.raises(OceanTideError, match="not a thin layer"):
        _earth(liquid=Liquid("mercury", 13534.0))
    ethane = _earth(liquid=Liquid("ethane", 540.0))
    assert ethane.loading_factor < _earth().loading_factor
    assert ocean_tides.loading_factor(0.0, 5514.0) == 1.0
    # the same tide in height without the loading: the equilibrium height does not depend on the liquid
    assert ethane.raisers[0].equilibrium_range_m == _earth().raisers[0].equilibrium_range_m


def test_an_ice_covered_ocean_has_no_free_surface_tide():
    tide = _earth(ice_covered=True)
    assert tide.spring_range_m == tide.neap_range_m == 0.0 and "ice" in tide.reason
    assert _earth().reason is None
