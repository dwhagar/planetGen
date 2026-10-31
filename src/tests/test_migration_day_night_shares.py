# tests/test_migration_day_night_shares.py

"""
GEN.767: migration 0170 adds the day and night shares and the hours of daylight (`physics/daylight.py`) to `planets` and `moons`. A body
from before v170 keeps its rows with NULL there; a saved system keeps what it was made with.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import daylight
from planetgen.util import draw
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _types(conn, table):
    return {row["COLUMN_NAME"]: row["DATA_TYPE"] for row in conn.execute(
        "SELECT COLUMN_NAME, DATA_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", (table,)).fetchall()}


def test_a_v169_galaxy_gets_the_columns_empty(mysql_config):
    load_old_schema(mysql_config, 169)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert not set(daylight.FIELDS) & set(_types(conn, table))
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            types = _types(conn, table)
            assert all(types[name] == ("varchar" if name == "day_night_state" else "double") for name in daylight.FIELDS), table  # GEN.1097: the state is a word
    finally:
        conn.close()


def test_a_saved_system_keeps_the_shares_of_every_body(mysql_config):
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(getattr(planet, "moons", None) for planet in system.planets):
            break
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        made = {planet.name: planet for planet in system.planets if planet.body_type != "a"}
        checked = 0
        for planet in loaded.planets:
            if planet.body_type == "a":
                continue
            for name in daylight.FIELDS:
                assert getattr(planet, name) == pytest.approx(getattr(made[planet.name], name)), (planet.name, name)
            for moon, original in zip(planet.moons, made[planet.name].moons):
                assert moon.parent_shadow_share == pytest.approx(original.parent_shadow_share)
                checked += 1
        assert checked > 0
    finally:
        conn.rollback()
        conn.close()
