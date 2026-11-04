# tests/test_day_night.py

"""
GEN.1097: every object has a day-and-night state and counter. The pure functions (`physics/day_night.py`) name the state from the spin, the orbit and the
polar circle; the Solar System bodies below are test targets, never inputs. The state, the hours of night and the solar days in a year are stored on every
planet and moon, kept through a save and a load, and migration 0193 gives a v192 galaxy the columns, empty.
"""

import pytest

from planetgen.db import store
from planetgen.describe import sheet
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import daylight, day_night
from planetgen.util import draw
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

HOURS_PER_YEAR = 8766.0
FIELDS = ("day_night_state", "night_hours", "days_per_year")


def _state(rotation_hours, period_years, tilt_deg=23.4, **kwargs):
    return day_night.state(rotation_hours, period_years, daylight.polar_circle_deg(tilt_deg, 0.0), **kwargs)


@pytest.mark.parametrize("name, rotation_hours, period_years, tilt_deg, state", [
    ("Earth", 23.934, 1.0, 23.44, day_night.CYCLE),
    ("Mars", 24.623, 1.881, 25.19, day_night.CYCLE),
    ("Venus", 5832.5, 0.6152, 177.4, day_night.CYCLE),               # slow and backwards, still a day and a night
    ("Mercury", 1407.6, 0.2408, 0.03, day_night.RESONANT),           # 3:2, the star returns every two years
    ("a locked M-dwarf planet", 12.0 * 24.0, 12.0 * 24.0 / HOURS_PER_YEAR, 0.0, day_night.LOCKED),
    ("Io round Jupiter's year", 42.46, 11.86, 3.1, day_night.CYCLE),     # turns once in 1.77 d round Jupiter, so the Sun is up and down about every 42 h
    ("a belt body", 7.0, 4.6, 40.0, day_night.CYCLE),
    ("Uranus", 17.24, 84.0, 97.77, day_night.POLAR),
    ("a body with a year-long day", 1.0e6, 1.0, 10.0, day_night.NON_ROTATING),
])
def test_the_state_word_of_the_known_bodies(name, rotation_hours, period_years, tilt_deg, state):
    assert _state(rotation_hours, period_years, tilt_deg) == state, name


def test_a_star_has_no_night_and_a_giant_has_no_solid_surface_whatever_its_spin():
    assert day_night.state(609.0, None, None, star=True) == day_night.NO_NIGHT        # the Sun's 25 d
    assert day_night.state(9.9, 11.86, 66.0, solid=False) == day_night.NO_SURFACE      # Jupiter
    assert day_night.state(None, None, None, solid=False) == day_night.NO_SURFACE


def test_a_planet_whose_spin_is_not_stored_has_no_state():
    assert day_night.state(None, 1.0, 66.0) is None
    assert day_night.state(24.0, None, 66.0) is None


def test_the_counter_is_the_year_over_the_solar_day():
    assert day_night.days_per_year(HOURS_PER_YEAR, daylight.solar_day_hours(23.934, HOURS_PER_YEAR)) == pytest.approx(365.25, abs=0.5)
    assert day_night.days_per_year(HOURS_PER_YEAR * 0.2408, daylight.solar_day_hours(1407.6, HOURS_PER_YEAR * 0.2408)) == pytest.approx(0.5, abs=0.01)
    assert day_night.days_per_year(HOURS_PER_YEAR, None) is None
    assert day_night.days_per_year(HOURS_PER_YEAR, 0.0) is None
    assert day_night.night_hours(24.0, 12.0) == 12.0
    assert day_night.night_hours(None, 12.0) is None


def _bodies(system):
    return [body for planet in system.planets if hasattr(planet, "moons") for body in [planet, *planet.moons]]


def _system():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(200):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        bodies = _bodies(system)
        if any(body.body_type == "g" for body in bodies) and any(b.day_night_state == day_night.CYCLE for b in bodies):
            return system, cfg
    pytest.fail("no system with a giant and a day and night cycle in 200 seeds")


def test_every_generated_planet_and_moon_has_a_state_that_agrees_with_its_spin():
    system, _cfg = _system()
    for body in _bodies(system):
        assert body.day_night_state in day_night.STATES, body.day_night_state
        if body.body_type == "g":
            assert body.day_night_state == day_night.NO_SURFACE
        if body.day_night_state == day_night.LOCKED:
            assert body.solar_day_hours is None and body.days_per_year is None and body.night_hours is None
        if body.day_night_state in (day_night.CYCLE, day_night.POLAR) and body.solar_day_hours:
            assert body.days_per_year > 0.0
            assert 0.0 <= body.night_hours <= body.solar_day_hours


def _stored(bodies):
    return [{name: getattr(body, name, None) for name in FIELDS} for body in bodies]


def test_a_saved_system_keeps_the_state_and_the_counter(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        assert _stored(_bodies(store.load_star_system(conn, system_id))) == _stored(_bodies(system))
    finally:
        conn.rollback()
        conn.close()


def test_the_sheet_lists_the_state_the_hours_of_night_and_the_days_in_a_year():
    for kind in ("planet", "moon"):
        keys = {spec.key: spec for spec in sheet.specs_for(kind)}
        assert keys["day_night_state"].label == "State of its day and night" and keys["day_night_state"].quantity == sheet.TEXT
        assert keys["night_hours"].label == "Hours of night in a solar day, equator"
        assert keys["days_per_year"].label == "Solar days in one year"
        assert sheet._value(keys["day_night_state"], "planets", "tidally locked") == ["Tidally locked"]


def test_a_v192_galaxy_gets_the_columns_empty(mysql_config):
    load_old_schema(mysql_config, 192)
    assert store.schema_status(mysql_config)[0] == 192
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert not conn.execute("SHOW COLUMNS FROM " + table + " LIKE 'day_night_state'").fetchall()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            for column in FIELDS:
                assert conn.execute("SHOW COLUMNS FROM " + table + " LIKE ?", (column,)).fetchall(), (table, column)
    finally:
        conn.close()
