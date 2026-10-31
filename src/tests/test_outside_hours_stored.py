# tests/test_outside_hours_stored.py

"""
GEN.768 part 2: the hours of one local day a person can be outside, and the limit that sets them, are stored on every scored rocky planet
and moon. A system saved and loaded keeps them, a dose refresh gives the same ones from the stored columns alone, a gas giant has none, and
migration 0187 gives a v186 galaxy the two columns, empty.
"""

import pytest

from planetgen.db import store
from planetgen.describe import sheet
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import outside_hours
from planetgen.util import draw
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

FIELDS = outside_hours.FIELDS


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
        if any(body.outside_hours_day is not None for body in bodies) and any(body.body_type == "g" for body in bodies):
            return system, cfg
    pytest.fail("no system with a scored world and a gas giant in 200 seeds")


def _stored(bodies):
    return [{name: getattr(body, name, None) for name in FIELDS} for body in bodies]


def test_a_scored_world_has_hours_in_a_day_and_the_limit_that_sets_them():
    system, _cfg = _system()
    scored = [body for body in _bodies(system) if body.outside_hours_day is not None]
    assert scored
    for body in scored:
        assert body.outside_hours_limit in outside_hours.LIMITS
        assert 0.0 <= body.outside_hours_day <= 24.0
    for giant in (body for body in _bodies(system) if body.body_type == "g"):
        assert (giant.outside_hours_day, giant.outside_hours_limit) == (None, None)


def test_a_body_that_needs_a_vehicle_has_no_hours_and_says_so():
    system, _cfg = _system()
    for body in _bodies(system):
        if body.outside_hours_limit == "vehicle":
            assert body.outside_hours_day == 0.0


def test_a_saved_system_keeps_the_hours_and_a_dose_refresh_gives_the_same(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        before = _stored(_bodies(system))
        assert _stored(_bodies(store.load_star_system(conn, system_id))) == before
        with conn:
            store.refresh_habitability(conn, system_id)
        assert _stored(_bodies(store.load_star_system(conn, system_id))) == before
    finally:
        conn.rollback()
        conn.close()


def test_every_limit_has_words_for_a_page_and_the_sheet_lists_the_rows():
    assert set(outside_hours.LIMIT_WORDS) == set(outside_hours.LIMITS)
    for kind in ("planet", "moon"):
        keys = {spec.key: spec for spec in sheet.specs_for(kind)}
        assert keys["outside_hours_day"].label == "Hours a person can be outside in a day"
        assert keys["outside_hours_limit"].quantity == sheet.LIMIT
        assert keys["outside_hours"].label == "Hours a person can work outside per year (dose limit)"
        assert sheet._value(keys["outside_hours_limit"], "planets", "uv") == [outside_hours.LIMIT_WORDS["uv"]]
        assert sheet._value(keys["outside_hours_limit"], "planets", None) == [sheet.NOT_STORED]
    assert not sheet.applies("planet", "J", "outside_hours_day") and sheet.applies("planet", "M", "outside_hours_day")


def test_a_v186_galaxy_gets_the_columns_empty(mysql_config):
    load_old_schema(mysql_config, 186)
    assert store.schema_status(mysql_config)[0] == 186
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert not conn.execute("SHOW COLUMNS FROM " + table + " LIKE 'outside_hours_day'").fetchall()
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
