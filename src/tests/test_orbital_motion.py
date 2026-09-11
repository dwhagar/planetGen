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
Also covers the later position/speed follow-up: `Planet.position_x/y/z`/
`orbital_speed_kms`, derived from the orbital elements above via
`utils.orbital_position_au`/`circular_orbital_speed_kms` and
`planetPhysics.update_orbital_position` (called both at generation and by
`StarSystem.validate_system` whenever it corrects a planet's `distance`);
and the floating-point update-guard follow-up after that:
`Planet.min_update_interval_years` (`utils.minimum_update_interval_years`,
`period_years * math.ulp(360.0) / 360` -- the shortest `elapsed_years`
worth calling `_db.advance_orbital_phases` for before the phase delta
added would be too small to change the stored value at all).

DB-backed tests for the persistence/migration/update-script side of this
feature (`insert_planet`/`insert_moon`'s new columns,
`advance_orbital_phases`, `migrate_database`'s v8->v9/v10->v11/v11->v12
steps) live in `test_db_persistence.py` instead, alongside every other
`_db.py` test.
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
        "orbital_phase_deg", "position_x", "position_y", "position_z",
        "orbital_speed_kms", "min_update_interval_years", "rotation_period_hours",
    ):
        assert field in Planet.SERIALIZABLE_FIELDS


# ---------------------------------------------------------------------------
# Position/speed (CHANGELOG's "Add planet/moon position" entry): pure-math
# helpers (`utils.orbital_position_au`/`circular_orbital_speed_kms`), plus
# their wiring into generation via `generate_orbital_motion_properties`'s
# call to `planetPhysics.update_orbital_position`.
# ---------------------------------------------------------------------------

from stellarObjects.utils import circular_orbital_speed_kms, orbital_position_au


def test_orbital_position_au_lands_on_the_orbit_radius_sphere():
    """
    Regardless of inclination/ascending-node/phase, the derived position
    must sit at exactly `distance_au` from the origin (the orbital
    anchor) -- a circular orbit by definition.
    """
    for inclination in (0, 15, 45, 89.9):
        for node in (0, 90, 180, 270):
            for phase in (0, 45, 90, 180, 270, 359):
                x, y, z = orbital_position_au(10.0, inclination, node, phase)
                radius = math.sqrt(x ** 2 + y ** 2 + z ** 2)
                assert radius == pytest.approx(10.0, rel=1e-9)


def test_orbital_position_au_zero_inclination_is_flat_with_degenerate_node():
    """
    An uninclined orbit (i=0) is flat in the primary's own reference plane
    (z=0 always). Node and phase become degenerate at i=0 -- with no
    inclined plane for the ascending node to describe the crossing of,
    only their sum matters (the standard, expected behavior at this
    degenerate case, not a bug): `x = r*cos(node+phase)`,
    `y = r*sin(node+phase)`.
    """
    for node in (0, 90, 180, 275):
        x, y, z = orbital_position_au(5.0, 0.0, node, 30.0)
        assert z == pytest.approx(0.0, abs=1e-9)
        combined = math.radians(node + 30.0)
        assert x == pytest.approx(5.0 * math.cos(combined))
        assert y == pytest.approx(5.0 * math.sin(combined))


def test_orbital_position_au_zero_phase_lands_on_the_x_axis():
    """
    A simple, hand-checkable case: node=0, phase=0 puts the body on the
    primary's own +X axis regardless of inclination (u=0 means cos(u)=1,
    sin(u)=0, so the sin(u) terms all vanish).
    """
    for inclination in (0, 30, 60, 90):
        x, y, z = orbital_position_au(7.0, inclination, 0.0, 0.0)
        assert (x, y, z) == pytest.approx((7.0, 0.0, 0.0), abs=1e-9)


def test_circular_orbital_speed_kms_matches_earth_orbit_order_of_magnitude():
    """Earth: ~1 AU, ~1 year period, real orbital speed ~29.8 km/s."""
    speed = circular_orbital_speed_kms(1.0, 1.0)
    assert 29.0 <= speed <= 30.5


def test_circular_orbital_speed_kms_matches_manual_formula():
    speed = circular_orbital_speed_kms(2.5, 3.2)
    expected = (2 * math.pi * 2.5 * pc.AU_TO_KM) / (3.2 * pc.SECONDS_PER_YEAR)
    assert speed == pytest.approx(expected, rel=1e-9)


def test_generated_bodies_have_position_on_the_distance_sphere_and_matching_speed(bodies):
    planets, moons = bodies
    for body in planets + moons:
        radius = math.sqrt(body.position_x ** 2 + body.position_y ** 2 + body.position_z ** 2)
        assert radius == pytest.approx(body.distance, rel=1e-9)

        expected_speed = circular_orbital_speed_kms(body.distance, body.period)
        assert body.orbital_speed_kms == pytest.approx(expected_speed, rel=1e-9)


def test_moon_position_is_relative_to_its_planet_not_the_star(bodies):
    """
    A moon's `distance`/orbital elements (and so its derived position) are
    relative to its parent planet, not the star -- a moon's `position_x/y/z`
    magnitude should match its own (planet-relative) `distance`.
    """
    planets, moons = bodies
    for planet in planets:
        for moon in planet.moons:
            moon_radius = math.sqrt(moon.position_x ** 2 + moon.position_y ** 2 + moon.position_z ** 2)
            assert moon_radius == pytest.approx(moon.distance, rel=1e-9)


# ---------------------------------------------------------------------------
# Floating-point update guard (CHANGELOG's "Floating-point update guard"
# entry): `utils.minimum_update_interval_years`, plus its wiring into
# generation via `planetPhysics.update_orbital_position`.
# ---------------------------------------------------------------------------

from stellarObjects.utils import minimum_update_interval_years


def test_minimum_update_interval_years_matches_manual_formula():
    interval = minimum_update_interval_years(10.0)
    expected = 10.0 * math.ulp(360.0) / 360
    assert interval == pytest.approx(expected, rel=1e-9)


def test_minimum_update_interval_years_scales_linearly_with_period():
    """A pure function of period -- doubling the period should exactly
    double the guard interval (both sides of the formula are linear in
    period_years)."""
    base = minimum_update_interval_years(5.0)
    doubled = minimum_update_interval_years(10.0)
    assert doubled == pytest.approx(2 * base, rel=1e-9)


def test_minimum_update_interval_years_is_far_below_any_realistic_cadence():
    """
    For a 1-year period, the guard floor should be many orders of
    magnitude below `updateOrbits.py`'s own "once a month or so" cadence
    (a month is roughly 0.083 years) -- confirming this guard exists for
    correctness against pathological callers, not because real usage ever
    comes close to it.
    """
    interval = minimum_update_interval_years(1.0)
    assert interval < 1e-10
    assert interval > 0


def test_advancing_phase_by_less_than_the_guard_interval_is_a_true_float_noop():
    """
    Directly confirms the floating-point claim the guard is built on: at a
    delta smaller than the guard interval, `MOD(phase + delta, 360)` (the
    same expression `advance_orbital_phases` uses) really does round back
    to the exact original phase in Python's own `float`, the same IEEE 754
    double `DOUBLE` uses.
    """
    period_years = 3.0
    phase_deg = 271.3384217  # arbitrary, not near a clean boundary
    interval = minimum_update_interval_years(period_years)

    # A comfortable margin below the guard interval (not exactly at it,
    # to avoid a rounding-to-even tie-break right at the ULP boundary).
    delta_deg = (interval * 0.1 / period_years) * 360
    new_phase = math.fmod(phase_deg + delta_deg, 360)
    assert new_phase == phase_deg


def test_generated_bodies_have_finite_positive_min_update_interval(bodies):
    planets, moons = bodies
    for body in planets + moons:
        assert math.isfinite(body.min_update_interval_years)
        assert body.min_update_interval_years > 0

        expected = minimum_update_interval_years(body.period)
        assert body.min_update_interval_years == pytest.approx(expected, rel=1e-9)
