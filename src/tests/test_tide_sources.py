"""
GEN.1446: every body of a system raises a tide on every other (physics/tide_sources.py): the means of `1 / d^3`, the list of sources of a planet and a moon, the extra raisers of the
ocean tide, and the effect ledger entries.
"""

import math

import pytest

from planetgen.physics import constants, ocean_tides, tide_sources


def test_the_mean_of_inverse_cube_of_two_circles_is_the_quadrature():
    for a, b in ((1.0, 1.5), (1.0, 1.05), (2.0, 1.0), (1.0, 5.0)):
        n = 20000
        mean = sum((a * a + b * b - 2 * a * b * math.cos(2 * math.pi * (k + 0.5) / n)) ** -1.5 for k in range(n)) / n
        assert tide_sources.circles_mean_inverse_cube(a, b) == pytest.approx(mean, rel=1e-6)


def test_a_far_source_is_one_over_its_distance_cubed_and_the_kepler_mean_has_the_eccentricity_factor():
    assert tide_sources.circles_mean_inverse_cube(1.0, 100.0) == pytest.approx(1.0e-6, rel=1e-3)
    assert tide_sources.kepler_mean_inverse_cube(2.0, 0.0) == pytest.approx(1.0 / 8.0)
    assert tide_sources.kepler_mean_inverse_cube(2.0, 0.6) == pytest.approx(1.0 / (8.0 * 0.64 ** 1.5))
    assert tide_sources.orbit_mean_inverse_cube(1.0, 100.0, 0.5) > tide_sources.orbit_mean_inverse_cube(1.0, 100.0, 0.0)


class _Body:
    def __init__(self, mass, radius, distance=None, eccentricity=0.0, star=None, moons=(), name="b"):
        self.mass, self.radius, self.distance, self.eccentricity, self.star, self.moons, self.name = mass, radius, distance, eccentricity, star, list(moons), name


def _system():
    star = _Body(constants.SOLAR_MASS_TO_KG, 6.9e5, name="star")
    earth = _Body(constants.EARTH_MASS_TO_KG, 6371.0, 1.0, 0.0167, star, name="earth")
    moon = _Body(7.35e22, 1737.0, 0.00257, 0.05, star, name="moon")
    jupiter = _Body(317.8 * constants.EARTH_MASS_TO_KG, 69911.0, 5.2, 0.048, star, name="jupiter")
    earth.moons = [moon]
    return star, earth, moon, jupiter


def test_a_planet_is_raised_by_its_host_a_sibling_and_its_own_moon_and_the_sibling_is_tiny():
    star, earth, moon, jupiter = _system()
    sources = tide_sources.sources_of(earth, None, [earth, jupiter])
    by_name = {s.name: s for s in sources}
    assert set(by_name) == {"moon", "star", "jupiter"}
    accelerations = {name: tide_sources.acceleration_m_s2(earth, s) for name, s in by_name.items()}
    assert accelerations["moon"] > accelerations["star"] > accelerations["jupiter"] > 0.0
    assert accelerations["jupiter"] / accelerations["star"] < 1e-4


def test_only_the_sources_the_ocean_tide_did_not_have_become_raisers():
    star, earth, moon, jupiter = _system()
    sources = tide_sources.sources_of(earth, None, [earth, jupiter])
    extra = tide_sources.raisers_of(sources, tide_sources.counted_by_ocean_tides(earth, None))
    assert [r.name for r in extra] == ["jupiter"]
    assert isinstance(extra[0], ocean_tides.Raiser) and extra[0].distance_m > 3.0e11 and extra[0].eccentricity == 0.0
    assert extra[0].orbital_rate_rad_s > 0.0


def test_a_moon_is_raised_by_its_planet_the_host_the_siblings_and_the_other_moons():
    star, earth, moon, jupiter = _system()
    second = _Body(1.0e22, 1000.0, 0.0008, 0.01, star, name="inner moon")
    earth.moons.append(second)
    sources = tide_sources.sources_of(moon, earth, [earth, jupiter])
    names = [s.name for s in sources]
    assert names[0] == "earth" and {"star", "jupiter", "inner moon"} <= set(names)
    extra = tide_sources.raisers_of(sources, tide_sources.counted_by_ocean_tides(moon, earth))
    assert {r.name for r in extra} == {"jupiter", "inner moon"}


def test_a_lone_body_has_no_extra_raiser_and_the_ocean_tide_is_what_it_was():
    star, earth, moon, jupiter = _system()
    earth.moons = []
    sources = tide_sources.sources_of(earth, None, [earth])
    assert [s.name for s in sources] == ["star"]
    assert tide_sources.raisers_of(sources, tide_sources.counted_by_ocean_tides(earth, None)) == []


def test_extra_raisers_add_to_the_ranges_of_the_ocean_tide():
    base = ocean_tides.Raiser.of("star", constants.SOLAR_MASS_TO_KG, 1.496e11, constants.EARTH_MASS_TO_KG)
    extra = ocean_tides.Raiser("jupiter", 317.8 * constants.EARTH_MASS_TO_KG, 6.3e11, 0.0, 1.0e-8)
    args = (constants.EARTH_MASS_TO_KG, 6.371e6, 7.292e-5, [base], 0.3)
    without = ocean_tides.ocean_tide(*args)
    with_jupiter = ocean_tides.ocean_tide(constants.EARTH_MASS_TO_KG, 6.371e6, 7.292e-5, [base, extra], 0.3)
    assert with_jupiter.spring_range_m > without.spring_range_m
    assert with_jupiter.spring_range_m / without.spring_range_m - 1.0 < 1e-3
