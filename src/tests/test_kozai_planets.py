# tests/test_kozai_planets.py

"""GEN.1447 slice 0: the Kozai-Lidov screen of a planet by a companion star or an outer sibling (physics/kozai_planets.py), pure and unwired."""

import math

import pytest

from planetgen.physics import constants, kozai_planets as kp
from planetgen.physics.constants import EARTH_MASS_TO_KG

SUN_KG = constants.SOLAR_MASS_TO_KG
SUN_RADIUS_KM = 695700.0
JUPITER_KG = 1.898e27
AGE = 5.0e9


def planet(axis=1.0, mass=EARTH_MASS_TO_KG, radius=6371.0, e=0.02, inc=0.0, node=0.0, key="p"):
    return kp.PlanetInput(key, mass, radius, axis, e, math.radians(inc), math.radians(node), )


def companion(axis=300.0, mass=0.5, e=0.3, inc=60.0, node=0.0, key="c"):
    return kp.Perturber(key, "companion", mass, axis, e, math.radians(inc), math.radians(node))


def one(planets, companions=()):
    return kp.screen_planets(SUN_KG, SUN_RADIUS_KM, AGE, planets, companions)


def test_a_lone_planet_with_no_companion_is_quiet():
    (result,) = one([planet()])
    assert result.outcome is None and result.perturber is None and result.cycling_perturbers == () and result.ledger_rates == {}


def test_a_coplanar_companion_cannot_pump_it():
    (result,) = one([planet(inc=0.0)], [companion(inc=10.0)])
    assert result.outcome is not None and not result.outcome.cycles and result.outcome.cause == "window"
    assert result.ledger_rates["c"] == pytest.approx(1.0 / result.outcome.timescale_yr)


def test_a_distant_companion_at_sixty_degrees_pumps_a_giant_of_a_few_au():
    (result,) = one([planet(axis=3.0, mass=JUPITER_KG, radius=71492.0)], [companion(axis=200.0, inc=60.0)])
    assert result.outcome.cycles and result.perturber == "c" and result.cycling_perturbers == ("c",)
    assert result.outcome.eccentricity_max > 0.5
    assert 0.0 <= result.outcome.eccentricity_now <= result.outcome.eccentricity_max + 1e-9


def test_a_near_polar_companion_swallows_a_planet_whose_periapsis_falls_inside_the_star():
    (result,) = one([planet(axis=3.0, mass=JUPITER_KG, radius=71492.0, e=0.05)], [companion(axis=150.0, mass=0.8, e=0.1, inc=88.0)])
    assert result.outcome.cycles and result.outcome.lost


def test_a_companion_that_is_too_close_is_left_to_the_stability_cut():
    (result,) = one([planet(axis=1.0)], [companion(axis=2.0)])
    assert result.cycling_perturbers == () and result.ledger_rates == {}


def test_an_outer_giant_sibling_is_a_perturber_but_its_own_pull_quenches_a_close_pair():
    planets = [planet(axis=1.0, key="in"), planet(axis=6.0, mass=JUPITER_KG, radius=71492.0, inc=60.0, key="out")]
    inner, outer = one(planets)
    assert "out" in inner.ledger_rates and outer.ledger_rates == {}
    assert not inner.outcome.cycles                 # the giant's secular pull (and relativity) turns the inner orbit faster than the Kozai rate


def test_a_lighter_outer_sibling_is_not_a_perturber():
    planets = [planet(axis=1.0, mass=JUPITER_KG, key="in"), planet(axis=6.0, mass=EARTH_MASS_TO_KG, inc=60.0, key="out")]
    inner, _outer = one(planets)
    assert inner.ledger_rates == {}


def test_the_strongest_cycling_perturber_is_chosen():
    near, far = companion(axis=150.0, key="near"), companion(axis=600.0, key="far")
    (result,) = one([planet(axis=3.0, mass=JUPITER_KG, radius=71492.0)], [far, near])
    assert result.cycling_perturbers[0] == "near" and result.perturber == "near"
    assert result.ledger_rates["near"] > result.ledger_rates["far"]


def test_the_outcome_is_repeatable_and_does_not_depend_on_other_planets_order():
    planets = [planet(axis=3.0, mass=JUPITER_KG, radius=71492.0, key="a"), planet(axis=40.0, mass=EARTH_MASS_TO_KG, key="b")]
    first = one(planets, [companion()])
    assert first == one(planets, [companion()])


def test_bad_inputs_are_refused():
    with pytest.raises(ValueError):
        kp.screen_planets(0.0, SUN_RADIUS_KM, AGE, [planet()])
    assert kp.screen_planets(SUN_KG, SUN_RADIUS_KM, AGE, []) == []
