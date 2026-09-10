# tests/test_orbital_motion.py

"""
Regression tests for the orbital-motion feature (docs/TODO.md's
"Introducing realistic orbital paths and speeds..." entry; CHANGELOG.md
[5.11.0]): `Planet.orbital_inclination_deg`/`orbital_ascending_node_deg`/
`orbital_phase_deg`/`rotation_period_hours`, the primary-mass-aware
`period` formula (Kepler's third law using the actual body a planet/moon
orbits, not always the star), and
`planetPhysics.generate_orbital_motion_properties`.

DB-backed tests for the persistence/migration/update-script side of this
feature (`insert_planet`/`insert_moon`'s new columns,
`advance_orbital_phases`, `migrate_database`'s v8->v9 step) live in
`test_db_persistence.py` instead, alongside every other `_db.py` test.
"""

import math

import pytest

from stellarObjects import physical_constants as pc
from stellarObjects.config import SystemConfig
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


def test_moon_tidal_lock_rate_is_close_to_configured_probability(bodies):
    """
    Statistical check, not exact -- across many independently generated
    moons, the tidally-locked fraction should land near
    MOON_TIDAL_LOCK_PROBABILITY (0.75), well outside plausible binomial
    noise for this sample size (used as a sanity check that the tidal-lock
    branch is wired up at all, not a precise calibration test).
    """
    _, moons = bodies
    assert len(moons) >= 30, "sample too small for a meaningful rate check"
    rate = sum(1 for m in moons if _is_tidally_locked(m)) / len(moons)
    assert 0.55 < rate < 0.9


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
