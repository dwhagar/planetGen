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


def _planets_with_axes(systems):
    return [p for s in systems for p in s.planets if hasattr(p, "moons") and getattr(p, "spin_axis_x", None) is not None]


def _obliquity(planet):
    from planetgen.physics import spin
    return spin.obliquity_between((planet.spin_axis_x, planet.spin_axis_y, planet.spin_axis_z),
                                  spin.orbit_normal(planet.orbital_inclination_deg, planet.orbital_ascending_node_deg))


def test_every_planet_keeps_a_spin_axis_that_matches_its_stored_tilt_after_the_exchange(systems):
    """GEN.1504: the exchange moves a planet's orbital plane, and its spin axis goes with it, so the obliquity is the stored tilt."""
    checked = 0
    for planet in _planets_with_axes(systems):
        assert _obliquity(planet) == pytest.approx(planet.axial_tilt_deg, abs=1e-6)
        checked += 1
    assert checked > 10


def test_a_free_spinner_keeps_its_obliquity_when_its_plane_moves_and_a_locked_one_keeps_its_share_of_the_inclination(systems):
    from planetgen.physics import spin
    free = locked = 0
    for planet in _planets_with_axes(systems):
        old_i, old_n, tilt = planet.orbital_inclination_deg, planet.orbital_ascending_node_deg, planet.axial_tilt_deg
        state = spin.spin_state(planet.rotation_period_hours, planet.period)[0]
        planet.orbital_inclination_deg, planet.orbital_ascending_node_deg = old_i + 0.7, (old_n + 40.0) % 360.0
        planetPhysics.carry_spin_with_orbit(planet, old_i, old_n)
        if state in ("locked", "resonant") and min(old_i, 180.0 - old_i) > 1e-9:
            share = min(planet.orbital_inclination_deg, 180.0 - planet.orbital_inclination_deg) / min(old_i, 180.0 - old_i)
            assert planet.axial_tilt_deg == pytest.approx(tilt * share, rel=1e-9)
            locked += 1
        else:
            assert planet.axial_tilt_deg == pytest.approx(tilt, abs=1e-9)
            free += 1
        assert _obliquity(planet) == pytest.approx(planet.axial_tilt_deg, abs=1e-6)
    assert free > 5


def test_carry_axis_keeps_the_angle_round_the_normal_and_a_plane_turned_over_flips_the_axis():
    from planetgen.physics import spin
    old, new = spin.orbit_normal(5.0, 30.0), spin.orbit_normal(12.0, 200.0)
    axis = spin.tilted_axis(old, 23.0, 1.3)
    turned = spin.carry_axis(axis, old, new)
    assert spin.obliquity_between(turned, new) == pytest.approx(23.0, abs=1e-9)
    back = spin.carry_axis(turned, new, old)
    assert back == pytest.approx(axis, abs=1e-9)
    assert spin.obliquity_between(spin.carry_axis(axis, old, new, 4.0), new) == pytest.approx(4.0, abs=1e-9)
    flipped = spin.carry_axis(axis, old, tuple(-c for c in old))
    assert spin.obliquity_between(flipped, tuple(-c for c in old)) == pytest.approx(23.0, abs=1e-9)


def test_a_belt_pulls_on_the_planets_of_its_star_and_a_system_without_belt_pull_differs(monkeypatch):
    """GEN.1464: the same system made with the belts' rings left out of the exchange leaves the planets on other orbits; the belts are not in the exchange's result."""
    with_belt = _system(1)
    belts = [p for p in with_belt.planets if p.body_type == 'a']
    assert belts and belts[0].mass_kg > 0.0, "seed 1 makes a belt with a mass"
    monkeypatch.setattr(planetPhysics, "belt_rings_of", lambda planets, axes: ([], []))
    without = _system(1)
    planets_with = [p for p in with_belt.planets if p.body_type != 'a']
    planets_without = [p for p in without.planets if p.body_type != 'a']
    assert len(planets_with) == len(planets_without) > 1
    assert any(abs(a.eccentricity - b.eccentricity) > 1e-6 for a, b in zip(planets_with, planets_without)), "the belt's pull moves at least one planet's orbit"


def test_the_rings_of_a_belt_keep_clear_of_the_planets_and_carry_the_belts_mass():
    for seed in range(2, 20):           # the first system with a belt that has a mass and a ring clear of the planets (the seed that does moves with the generator)
        system = _system(seed)
        belts = [p for p in system.planets if p.body_type == 'a' and p.mass_kg > 0.0]
        bodies = [p for p in system.planets if p.body_type != 'a']
        axes, masses = planetPhysics.belt_rings_of(system.planets, [p.distance for p in bodies])
        if axes:
            break
    assert axes and sum(masses) <= sum(b.mass_kg for b in belts) * (1.0 + 1e-9)
    assert all(min(a, p.distance) / max(a, p.distance) <= 0.9 for a in axes for p in bodies)
    assert planetPhysics.belt_rings_of([p for p in system.planets if p.body_type != 'a'], [1.0]) == ([], [])
