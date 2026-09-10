# tests/test_orbital_motion.py

"""
Regression tests for the orbital-motion feature (docs/TODO.md's
"Introducing realistic orbital paths and speeds..." entry; CHANGELOG.md
[5.11.0]/[5.11.1]): `Planet.orbital_inclination_deg`/
`orbital_ascending_node_deg`/`orbital_phase_deg`/`rotation_period_hours`,
the primary-mass-aware `period` formula (Kepler's third law using the
actual body a planet/moon orbits, not always the star), and
`planetPhysics.generate_orbital_motion_properties`/
`_tidal_locking_timescale_seconds` (real tidal-despinning physics decides
whether a moon is tidally locked, not a flat probability -- see [5.11.1]).

DB-backed tests for the persistence/migration/update-script side of this
feature (`insert_planet`/`insert_moon`'s new columns,
`advance_orbital_phases`, `migrate_database`'s v8->v9 step) live in
`test_db_persistence.py` instead, alongside every other `_db.py` test.
"""

import math
import types

import pytest

from stellarObjects import physical_constants as pc
from stellarObjects.config import SystemConfig
from stellarObjects.planetPhysics import _tidal_locking_timescale_seconds, generate_orbital_motion_properties
from stellarObjects.starData import Star
from stellarObjects.systemData import StarSystem

N_SYSTEMS = 20


def _generate_systems_with_moons(n=N_SYSTEMS):
    """Generates `n` independent systems biased toward having moons,
    retrying (bounded) any generation that comes up moonless."""
    systems = []
    for _ in range(n * 3):
        if len(systems) >= n:
            break
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.MOONS = True
        cfg.MAX_PLANETS = True
        system = StarSystem(system_config=cfg)
        planets = [obj for obj in system.planets if obj.body_type != "a"]
        if any(p.moons for p in planets):
            systems.append(system)
    if not systems:
        pytest.fail("could not generate any system with at least one moon")
    return systems


def _all_planets_and_moons(systems):
    planets = []
    moons = []
    for system in systems:
        for obj in system.planets:
            if obj.body_type == "a":
                continue
            planets.append(obj)
            moons.extend(obj.moons)
    return planets, moons


@pytest.fixture(scope="module")
def systems_with_moons():
    return _generate_systems_with_moons()


@pytest.fixture(scope="module")
def bodies(systems_with_moons):
    planets, moons = _all_planets_and_moons(systems_with_moons)
    assert planets and moons
    return planets, moons


def test_orbital_orientation_and_phase_are_in_range(bodies):
    planets, moons = bodies
    for body in planets + moons:
        assert 0 <= body.orbital_ascending_node_deg < 360
        assert 0 <= body.orbital_phase_deg < 360

    for planet in planets:
        assert 0 <= planet.orbital_inclination_deg <= pc.PLANET_ORBITAL_INCLINATION_MAX_DEG
    for moon in moons:
        assert 0 <= moon.orbital_inclination_deg <= pc.MOON_ORBITAL_INCLINATION_MAX_DEG


def _is_tidally_locked(moon):
    expected = moon.period * pc.SECONDS_PER_YEAR / 3600
    return moon.rotation_period_hours == pytest.approx(expected, rel=1e-9)


def test_rotation_period_within_body_type_range_or_tidally_locked(bodies):
    planets, moons = bodies
    for planet in planets:
        min_h, max_h = pc.ROTATION_PERIOD_RANGE_HOURS[planet.body_type]
        assert min_h <= planet.rotation_period_hours <= max_h

    for moon in moons:
        min_h, max_h = pc.ROTATION_PERIOD_RANGE_HOURS[moon.body_type]
        assert _is_tidally_locked(moon) or (min_h <= moon.rotation_period_hours <= max_h)


def test_moon_tidal_lock_is_not_universal(bodies):
    """
    Real tidal physics (`_tidal_locking_timescale_seconds`) decides
    locking, not a flat probability -- across many independently
    generated moons (spanning a wide range of mass/radius/distance/
    primary-mass/system-age combinations), that should produce both
    locked and unlocked outcomes, not push every moon to one extreme.
    """
    _, moons = bodies
    assert len(moons) >= 30, "sample too small for a meaningful check"
    locked = sum(1 for m in moons if _is_tidally_locked(m))
    assert 0 < locked < len(moons)


def test_tidal_locking_timescale_matches_real_earth_moon_order_of_magnitude():
    """
    Real Earth-Moon system: real estimates for how long tidal forces took
    to lock the Moon are commonly cited as tens of millions of years (the
    Moon has long since been observed in its locked state). Generously
    bounded -- this formula is an order-of-magnitude estimate (fixed
    Q/k2), not a precision calculation.
    """
    moon = types.SimpleNamespace(distance=384400 / pc.AU_TO_KM, mass=7.342e22, radius=1737.4)
    earth_mass_kg = 5.972e24
    t_lock_years = _tidal_locking_timescale_seconds(moon, earth_mass_kg, initial_rotation_period_hours=10.0) / pc.SECONDS_PER_YEAR
    assert 1e6 < t_lock_years < 5e8


def test_tidal_locking_timescale_increases_with_distance_and_decreases_with_primary_mass():
    close = types.SimpleNamespace(distance=0.001, mass=1e22, radius=1500)
    far = types.SimpleNamespace(distance=0.01, mass=1e22, radius=1500)
    t_close = _tidal_locking_timescale_seconds(close, 6e24, 10.0)
    t_far = _tidal_locking_timescale_seconds(far, 6e24, 10.0)
    assert t_far > t_close

    moon = types.SimpleNamespace(distance=0.01, mass=1e22, radius=1500)
    t_light_primary = _tidal_locking_timescale_seconds(moon, 1e23, 10.0)
    t_heavy_primary = _tidal_locking_timescale_seconds(moon, 1e26, 10.0)
    assert t_heavy_primary < t_light_primary


def test_generate_orbital_motion_properties_locks_a_close_moon_around_an_old_system():
    """
    Small/close moon, massive primary, old system -- the locking
    timescale at either end of the candidate-rotation-period range (10h
    or 1400h) comes out to a couple hundred years at most (verified
    directly), utterly dwarfed by a 5 Gy system age, so this must lock
    regardless of which candidate period the random draw lands on.
    """
    star = types.SimpleNamespace(age=5.0)
    moon = types.SimpleNamespace(
        is_moon=True, body_type="t", distance=0.0005, mass=5e20, radius=800,
        period=0.001, star=star,
    )
    generate_orbital_motion_properties(moon, primary_mass_kg=6e24)
    expected_locked_hours = moon.period * pc.SECONDS_PER_YEAR / 3600
    assert moon.rotation_period_hours == pytest.approx(expected_locked_hours)


def test_generate_orbital_motion_properties_does_not_lock_a_far_moon_around_a_young_system():
    """
    Large-orbit moon, light primary, young system -- the locking
    timescale at either end of the candidate-rotation-period range comes
    out to 1e14+ years (verified directly), so this must NOT lock
    regardless of which candidate period the random draw lands on.
    """
    star = types.SimpleNamespace(age=1.0)
    moon = types.SimpleNamespace(
        is_moon=True, body_type="t", distance=0.05, mass=5e22, radius=2000,
        period=5.0, star=star,
    )
    generate_orbital_motion_properties(moon, primary_mass_kg=6e23)
    min_h, max_h = pc.ROTATION_PERIOD_RANGE_HOURS["t"]
    assert min_h <= moon.rotation_period_hours <= max_h
    expected_locked_hours = moon.period * pc.SECONDS_PER_YEAR / 3600
    assert moon.rotation_period_hours != pytest.approx(expected_locked_hours)


def test_generate_orbital_motion_properties_never_evaluates_locking_for_a_planet():
    """A planet (is_moon=False) should never take the tidal-locking
    branch at all -- rotation_period_hours always comes from the plain
    body_type range, regardless of star.age."""
    star = types.SimpleNamespace(age=100.0)
    planet = types.SimpleNamespace(
        is_moon=False, body_type="t", distance=0.0005, mass=5e20, radius=800,
        period=0.001, star=star,
    )
    generate_orbital_motion_properties(planet, primary_mass_kg=6e30)
    min_h, max_h = pc.ROTATION_PERIOD_RANGE_HOURS["t"]
    assert min_h <= planet.rotation_period_hours <= max_h


def test_period_uses_the_actual_primary_not_always_the_star(bodies):
    """
    Kepler's third law: T(years) = sqrt(a(AU)^3 / M_primary(Msun)). A
    planet's primary is its host star; a moon's primary is its parent
    planet -- NOT the star, even though `moon.star` is set to the same
    (grand-owning) star object for other purposes (life chemistry, zone
    comparisons). Before this fix, every body used `self.star.mass`
    unconditionally, which is wrong for a moon by many orders of magnitude
    (a planet's mass is a tiny fraction of its star's).
    """
    planets, moons = bodies
    for planet in planets:
        primary_mass_sol = planet.star.mass / pc.SOLAR_MASS_TO_KG
        expected_period = math.sqrt(planet.distance ** 3 / primary_mass_sol)
        assert planet.period == pytest.approx(expected_period, rel=1e-9)

    for planet in planets:
        for moon in planet.moons:
            primary_mass_sol = planet.mass / pc.SOLAR_MASS_TO_KG
            expected_period = math.sqrt(moon.distance ** 3 / primary_mass_sol)
            assert moon.period == pytest.approx(expected_period, rel=1e-9)
            # The (wrong) star-mass-based formula would give a wildly
            # different answer -- confirms this isn't accidentally passing
            # because the two happen to coincide.
            star_mass_sol = moon.star.mass / pc.SOLAR_MASS_TO_KG
            wrong_period = math.sqrt(moon.distance ** 3 / star_mass_sol)
            if abs(star_mass_sol - primary_mass_sol) / primary_mass_sol > 0.01:
                assert moon.period != pytest.approx(wrong_period, rel=1e-6)


def test_serializable_fields_include_orbital_motion_attributes():
    from stellarObjects.planetData import Planet
    for field in (
        "orbital_inclination_deg", "orbital_ascending_node_deg",
        "orbital_phase_deg", "rotation_period_hours",
    ):
        assert field in Planet.SERIALIZABLE_FIELDS
