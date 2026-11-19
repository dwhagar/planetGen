"""GEN.1802: a planet is built where the spacing pass would put it, so the zone its dials are drawn for is the zone it ends in."""

import types

import pytest

from planetgen.generation import validation
from planetgen.generation.system import StarSystem

MSUN_KG = 1.98847e30
EARTH_KG = 5.972e24


def _planet(distance, mass_earth, eccentricity, star, body_type="t"):
    return types.SimpleNamespace(body_type=body_type, distance=distance, eccentricity=eccentricity, mass=mass_earth * EARTH_KG, star=star)


def _star():
    return types.SimpleNamespace(mass=MSUN_KG)


def test_a_planet_with_nothing_inside_it_is_born_where_the_plan_put_it():
    seed = types.SimpleNamespace(mass_kg=5 * EARTH_KG, eccentricity=0.1)
    assert StarSystem._distance_clear_of([], seed, _star()) == 0.0


def test_the_birth_distance_is_the_one_the_spacing_pass_asks_for():
    star = _star()
    inner = [_planet(1.0, 300.0, 0.2, star), _planet(2.1, 20.0, 0.05, star)]
    for mass_earth, eccentricity in ((3.0, 0.0), (3.0, 0.3), (200.0, 0.1)):
        seed = types.SimpleNamespace(mass_kg=mass_earth * EARTH_KG, eccentricity=eccentricity)
        born = StarSystem._distance_clear_of(inner, seed, star)
        asked = validation.min_distance_after_au(_planet(born, mass_earth, eccentricity, star), inner)
        assert born == pytest.approx(asked)
        assert born > inner[-1].distance


def test_a_belt_inside_is_left_to_the_pass_that_runs_once_the_belts_are_in():
    star = _star()
    belt = types.SimpleNamespace(body_type="a", distance=1.5, eccentricity=0.0, mass=1e18, star=star, lower_limit=1.4, upper_limit=1.6)
    seed = types.SimpleNamespace(mass_kg=5 * EARTH_KG, eccentricity=0.1)
    assert StarSystem._distance_clear_of([belt], seed, star) == 0.0
