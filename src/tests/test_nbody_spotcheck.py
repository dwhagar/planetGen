"""GEN.1468: the N-body check of the stability screens (`planetgen/physics/nbody_check.py`). The integrator is checked on systems whose answer is known; the generated systems of a few
seeds are then integrated and none the spacing screen passed may be broken by the integration."""

import pytest

from planetgen.physics import nbody_check as nb


def planet(earths, a, e=0.02, inc=0.5):
    return nb.Planet(earths * 3.003e-6, 4.26e-5 * earths ** 0.3, a, e, inc, 0.0)


def second_at(inner, gap, earths):
    """A planet of `earths` Earth masses `gap` mutual Hill radii outside `inner` (one sun's mass of star)."""
    a = inner.a_au
    while True:
        a += 0.001
        other = planet(earths, a)
        if nb.hill_separation(inner, other, 1.0) >= gap:
            return other


def test_a_wide_pair_stays_and_a_pair_inside_the_hill_spacing_has_an_encounter():
    inner = planet(10, 1.0)
    wide = nb.integrate(1.0, [inner, second_at(inner, 12.0, 10)], orbits=5000, seed=1)
    assert wide.verdict == "stable" and wide.energy_change < nb.ENERGY_TOLERANCE and wide.worst_axis_change < 0.01
    close = nb.integrate(1.0, [inner, second_at(inner, 2.0, 10)], orbits=5000, seed=1)
    assert close.verdict == "encounter" and close.bodies == (0, 1) and close.min_hill_separation < nb.ENCOUNTER_HILL


def test_a_planet_on_a_bound_orbit_is_given_back_its_elements_and_a_lone_planet_does_not_drift():
    one = planet(1, 1.0, e=0.1, inc=3.0)
    pos, vel = nb.state_vector(one, 1.0, 40.0, 100.0)
    mu = nb.G * (1.0 + one.mass_sol)
    energy = 0.5 * (vel ** 2).sum() - mu / (pos ** 2).sum() ** 0.5
    assert energy == pytest.approx(-mu / (2.0 * one.a_au), rel=1e-12)
    result = nb.integrate(1.0, [one, planet(1, 3.0)], orbits=2000, seed=3)
    assert result.verdict == "stable" and result.worst_axis_change < 1e-3


def test_a_system_the_spacing_screen_passed_is_not_one_the_integration_breaks():
    from planetgen.generation import validation
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from tests.fuzz_support import deterministic_entropy
    broken, checked = [], 0
    for seed in (1, 2, 3, 4):
        with deterministic_entropy(seed):
            system = StarSystem(system_config=SystemConfig())
        bodies = sorted((p for p in validation.slot_planets(system.all_planets) if p.body_type != "a" and getattr(p, "distance", None)), key=lambda p: p.distance)       # not the belt members promoted to planets (GEN.1510): the spacing screen exempts them, as it does their belt
        if len(system.stars) != 1 or len(bodies) < 2 or validation.check_orbit_spacing(system.all_planets):
            continue
        sun = nb.SOLAR_MASS_KG
        planets = [nb.Planet(p.mass / sun, p.radius * nb.AU_PER_KM, p.distance, p.eccentricity or 0.0, p.orbital_inclination_deg or 0.0, p.orbital_ascending_node_deg or 0.0) for p in bodies]
        outcome = nb.integrate(system.stars[0].mass / sun, planets, orbits=3000, seed=seed)
        checked += 1
        if outcome.verdict in ("collision", "ejection", "encounter"):
            broken.append((seed, outcome.verdict, outcome.bodies))
    assert checked >= 2
    assert not broken, f"the screens passed these systems and the integration broke them: {broken}"
