"""GEN.1493 slice 1: the air reads the system's late accretion (`generation/late_accretion.py`, `planets.remake_air`)."""

import types

import pytest

from planetgen.generation import late_accretion
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import impact_history, planets as planet_physics
from planetgen.util import draw
from tests.fuzz_support import deterministic_entropy

EARTH_KG = 5.9722e24


def _earth_like(seed):
    """A rocky planet of a G star at 1 AU on its own, with the late accretion the first air drew."""
    seen = []
    real = impact_history.late_accreted_fraction

    def record(mass_earth, constants):
        fraction = real(mass_earth, constants)
        seen.append(fraction * 100.0)
        return fraction

    impact_history.late_accreted_fraction = record
    try:
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        with deterministic_entropy(seed):
            star = Star(cfg)
            planet = Planet(cfg, star, star.habitable_zone, 1.0, mass=EARTH_KG, moon_count=0)
    finally:
        impact_history.late_accreted_fraction = real
    return planet, seen[-1]


def _air(planet):
    return (planet.atmospheric_pressure, planet.surface_temperature, dict(planet.atmosphere) if isinstance(planet.atmosphere, dict) else planet.atmosphere)


def test_making_the_air_again_with_the_late_accretion_it_drew_gives_the_planet_it_was():
    """The air solves twice (a first pass, then the one that is kept), each with its own draw of the late accretion; the replay holds both to the one percent, so a
    planet that drew nearly the same twice comes back as it was and one that drew differently comes back within a hair of it."""
    checked = 0
    for seed in range(1, 6):
        planet, percent = _earth_like(seed)
        if planet.body_type != "t" or getattr(planet, "_air_replay", None) is None:
            continue
        before = _air(planet)
        assert planet_physics.remake_air(planet, percent)
        after = _air(planet)
        assert after[0] == pytest.approx(before[0], rel=0.01)
        assert after[1] == pytest.approx(before[1], rel=0.01)
        checked += 1
    assert checked >= 3


def test_making_the_air_again_twice_with_one_percent_is_the_same_air():
    planet, percent = _earth_like(4)
    state = planet._air_replay
    assert planet_physics.remake_air(planet, percent)
    first = _air(planet)
    planet._air_replay = state
    assert planet_physics.remake_air(planet, percent)
    assert _air(planet) == first


def test_a_different_late_accretion_changes_the_air_and_not_the_run_stream():
    planet, percent = _earth_like(2)
    before = _air(planet)
    draw.set_run_seed(77)
    state = draw.getstate()
    assert planet_physics.remake_air(planet, percent * 3.0)
    assert draw.getstate() == state                      # no draw of the run's stream moved
    assert _air(planet) != before


def test_a_moon_or_a_planet_with_no_replay_is_left_alone():
    planet, percent = _earth_like(2)
    planet._air_replay = None
    assert not planet_physics.remake_air(planet, percent)
    moon = types.SimpleNamespace(is_moon=True, _air_replay=(None, None, None))
    assert not planet_physics.remake_air(moon, percent)


def test_the_floor_and_the_ceiling_hold():
    planet = types.SimpleNamespace(mass=EARTH_KG, radius=6371.0, distance=1.0)
    alone = [types.SimpleNamespace(body_type="t", is_moon=False, mass=EARTH_KG, radius=6371.0, distance=1.0)]
    mean = impact_history.late_accreted_percent_mean(1.0)
    assert late_accretion.percent_for(planet, alone, 1.0) == pytest.approx(late_accretion.FLOOR_SHARE * mean)        # no belt: the floor
    huge = types.SimpleNamespace(body_type="a", is_moon=False, primordial_mass_earth=1.0e4, lower_limit=2.0, upper_limit=3.0)
    assert late_accretion.percent_for(planet, [*alone, huge], 1.0) == impact_history.LATE_ACCRETED_PERCENT_MAX      # a vast belt: the ceiling


def test_a_giant_beside_a_belt_leaves_the_planet_inside_less_than_the_same_system_without_it():
    planet = types.SimpleNamespace(mass=EARTH_KG, radius=6371.0, distance=1.0)
    inner = types.SimpleNamespace(body_type="t", is_moon=False, mass=EARTH_KG, radius=6371.0, distance=1.0)
    belt = types.SimpleNamespace(body_type="a", is_moon=False, primordial_mass_earth=10.0, lower_limit=2.1, upper_limit=3.3)
    giant = types.SimpleNamespace(body_type="g", is_moon=False, mass=1.9e27, radius=69911.0, distance=5.2)
    assert late_accretion.percent_for(planet, [inner, belt, giant], 1.0) < 0.3 * late_accretion.percent_for(planet, [inner, belt], 1.0)


def test_a_whole_system_has_its_rocky_planets_set_by_the_belts_and_leaves_the_run_stream_to_the_rest():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = True
    cfg.MOONS = False
    found = 0
    for seed in range(1, 12):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=cfg)
        for planet in system.all_planets:
            if planet.body_type == "t" and not planet.is_moon and getattr(planet, "_late_accreted_percent", None) is not None:
                found += 1
                assert 0.0 < planet._late_accreted_percent <= impact_history.LATE_ACCRETED_PERCENT_MAX
                assert planet._air_replay is None
    assert found > 5
