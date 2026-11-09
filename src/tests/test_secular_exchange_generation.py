"""GEN.1450: the secular exchange in the generator (`StarSystem._secular_exchange`, `planets.secular_exchange`)."""

import math

import pytest

from planetgen.generation import validation
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import planets as planetPhysics
from planetgen.util import draw


def _system(seed):
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    with draw.bound(seed):
        return StarSystem(system_config=cfg)


@pytest.fixture(scope="module")
def systems():
    return [_system(seed) for seed in range(1, 11)]


def test_every_star_reports_an_outcome_the_exchange_knows(systems):
    seen = set()
    for system in systems:
        assert system.secular_outcomes
        seen.update(system.secular_outcomes.values())
    assert seen <= set(planetPhysics.SECULAR_OUTCOMES)
    assert seen - {"unchanged"}, "ten full systems should hold a multi-planet star the exchange acts on"


def test_a_star_with_one_planet_is_left_as_it_was(systems):
    for system in systems:
        for star, planets in validation.star_lists(system):
            if len([p for p in planets if p.body_type != 'a']) < 2:
                assert system.secular_outcomes[star.name] == "unchanged"


def test_the_spacing_and_the_moons_stay_valid_after_the_exchange(systems):
    for system in systems:
        assert validation.check_star_system(system) == []


def test_eccentricities_and_inclinations_stay_physical(systems):
    for system in systems:
        for body in system.all_planets:
            if body.body_type == 'a':
                continue
            assert 0.0 <= body.eccentricity < 1.0
            assert math.isfinite(body.orbital_inclination_deg) and 0.0 <= body.orbital_inclination_deg <= 180.0
            orbit = body.orbit_from_vector()
            assert orbit["eccentricity"] == pytest.approx(body.eccentricity, abs=1e-6)
            assert math.degrees(orbit["inclination"]) == pytest.approx(body.orbital_inclination_deg, abs=1e-6)


def test_the_same_seed_gives_the_same_orbits():
    first, second = _system(4), _system(4)
    assert [(p.distance, p.eccentricity, p.orbital_inclination_deg) for p in first.planets if p.body_type != 'a'] == \
           [(p.distance, p.eccentricity, p.orbital_inclination_deg) for p in second.planets if p.body_type != 'a']


def test_a_belt_just_short_of_its_keep_out_is_pushed_and_one_on_the_limit_is_not():
    """`space_orbits` and the exchange's guard share one rule (`validation.belt_needs_pushing`); a belt placed on the limit, as generation does, is not pushed by a rounding error."""
    from types import SimpleNamespace
    last = SimpleNamespace(body_type="t", distance=1.0, eccentricity=0.0, min_orbit_distance=0.05)
    assert not validation.belt_needs_pushing(SimpleNamespace(body_type="a", distance=1.05), last)
    assert not validation.belt_needs_pushing(SimpleNamespace(body_type="a", distance=1.05 - 1e-12), last)
    assert validation.belt_needs_pushing(SimpleNamespace(body_type="a", distance=1.05 - 1e-6), last)
    last.eccentricity = 0.01        # a raised eccentricity moves the apoapsis out, and the belt is then short of its keep-out
    assert validation.belt_needs_pushing(SimpleNamespace(body_type="a", distance=1.05), last)


def test_a_rocky_planet_the_exchange_moved_has_its_air_made_on_its_final_orbit(systems):
    """The tide wears an eccentricity down after the exchange, so the air has to be made after both: the stored temperature follows from the stored albedo, greenhouse and the light at the
    final orbit (the check of test_surface_air, over full systems where the exchange acts)."""
    from planetgen.physics import climate_draws, orbit_elements, radiation_field
    checked = 0
    moved = 0
    for system in systems:
        for planet in system.planets:
            if getattr(planet, "body_type", None) != "t" or planet.is_moon or (planet.water_mass_fraction or 0) <= 1.0e-4:
                continue
            flux = radiation_field.starlight_w_m2(planet.star, orbit_elements.flux_distance_au(float(planet.distance), planet.eccentricity)) * (planet.light_factor or 1.0)
            assert planet.surface_temperature == pytest.approx(climate_draws.surface_temperature_k(planet.bond_albedo, planet.greenhouse_factor, flux), rel=1e-6), planet.name
            checked += 1
            moved += system.secular_outcomes.get(planet.star.name) in ("exchanged", "held by the spacing")
    assert checked > 0 and moved > 0


def test_the_stale_air_mark_is_cleared_by_the_final_age_refresh(systems):
    for system in systems:
        for planet in system.planets:
            for body in (*([planet] if planet.body_type == "t" else []), *getattr(planet, "moons", ())):      # (a giant is not made again by the refresh)
                if getattr(body, "_age_drawn", None) is not None:
                    assert not planetPhysics.age_moved(body), (body.name, body._age_drawn)
