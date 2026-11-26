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


# GEN.1645: the screen wired into generation


def test_the_plane_of_a_normal_is_the_plane_it_came_from():
    from planetgen.physics import spin
    for inclination, node in ((0.3, 1.0), (2.0, 4.0), (0.0, 0.0)):
        got_i, got_node = kp.plane_of_normal(spin.orbit_normal(math.degrees(inclination), math.degrees(node)))
        assert got_i == pytest.approx(inclination)
        if inclination > 0.0:
            assert got_node == pytest.approx(node)


def test_a_tilted_normal_is_at_the_asked_angle_from_the_perturbers_and_on_the_planets_side():
    from planetgen.physics import spin
    planet_normal = spin.orbit_normal(20.0, 40.0)
    perturber_normal = spin.orbit_normal(70.0, 100.0)
    tilted = kp.tilted_normal(planet_normal, perturber_normal, math.radians(85.0))
    assert math.sqrt(sum(c * c for c in tilted)) == pytest.approx(1.0)
    assert sum(a * b for a, b in zip(tilted, perturber_normal)) == pytest.approx(math.cos(math.radians(85.0)))
    # on the planet's side: the three normals share a plane (their triple product is nil) and the tilted one leans toward the planet's
    a, b, c = planet_normal, perturber_normal, tilted
    triple = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])
    assert triple == pytest.approx(0.0, abs=1e-9)
    assert sum(x * y for x, y in zip(tilted, planet_normal)) > 0.0


def test_a_tilt_of_a_plane_on_the_perturbers_is_left_alone():
    normal = (0.0, 0.0, 1.0)
    assert kp.tilted_normal(normal, normal, 1.0) == normal


def test_a_planet_a_companion_pumps_keeps_the_cycle_and_the_eccentricity_it_left_it_with():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    cfg = SystemConfig()
    cfg.STARS = 2
    for seed in range(26, 70):      # the first seed whose two-star system has a planet a companion pumps (the draws move with every change to the generator, so no seed is fixed)
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        pumped = [p for p in system.all_planets if p.body_type != "a" and p.build.kozai is not None]
        if pumped:
            break
    assert system.kozai_planet_outcomes["screened"] > 0
    assert pumped and system.kozai_planet_outcomes["cycling"] >= len(pumped)
    for body in pumped:
        outcome = body.build.kozai
        assert outcome.cycles and not outcome.lost
        assert body.build.orbit_birth == pytest.approx(outcome.eccentricity_now)
        assert body.eccentricity <= outcome.eccentricity_now + 1.0e-9      # the tide only wears it down


def test_a_planet_the_pump_leaves_on_a_very_eccentric_orbit_is_held_to_the_bound_of_every_orbit(monkeypatch):
    """A cycle that carries a planet past `orbit_elements.MAX_ECCENTRICITY` (a class C planet at 1.06 AU left at 0.96 by the pump) was left there; every orbit of a planet is held to the bound the checks hold it to."""
    import dataclasses
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import kozai_planets, orbit_elements, planets
    from planetgen.util import draw
    cfg = SystemConfig()
    cfg.STARS = 2
    with draw.bound(26):
        system = StarSystem(system_config=cfg)
    star = system.star
    bodies = [p for p in system.primary.orbiters if p.body_type in ("t", "g") and p.mass and p.distance and p.radius and p.period]
    real = kozai_planets.screen_planets

    def forced(*args, **kwargs):
        results = real(*args, **kwargs)
        return [dataclasses.replace(r, outcome=dataclasses.replace(r.outcome, cycles=True, lost=False, eccentricity_now=0.97)) if r.outcome is not None and r.perturber is not None else r for r in results]
    monkeypatch.setattr(kozai_planets, "screen_planets", forced)
    planets.kozai_pump_planets(star, bodies)
    pumped = [p for p in bodies if p.build.kozai is not None and p.build.kozai.eccentricity_now is not None and p.build.orbit_birth is not None]
    assert pumped
    for body in pumped:
        cap = orbit_elements.max_planet_eccentricity(body.distance, star, body)
        assert body.eccentricity <= cap + 1.0e-12 and body.eccentricity <= orbit_elements.MAX_ECCENTRICITY
        assert body.build.orbit_birth <= cap + 1.0e-12 and body.build.kozai.eccentricity_now <= cap + 1.0e-12
