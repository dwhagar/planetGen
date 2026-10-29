# tests/test_migration_lid_regime.py

"""
GEN.721: migration 0169 adds `lid_regime` to `planets` and `moons`. A body from before v169 keeps its rows with NULL there (read as a
stagnant lid); a saved system keeps every rocky body's regime, mobile or stagnant.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _types(conn, table):
    return {row["COLUMN_NAME"]: row["DATA_TYPE"] for row in conn.execute(
        "SELECT COLUMN_NAME, DATA_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", (table,)).fetchall()}


def test_a_v168_galaxy_gets_the_column_empty(mysql_config):
    load_old_schema(mysql_config, 168)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert "lid_regime" not in _types(conn, table)
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert _types(conn, table)["lid_regime"] == "varchar", table
    finally:
        conn.close()


def test_a_saved_system_keeps_every_rocky_bodys_regime(mysql_config):
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        rocky = [p for p in system.planets if p.body_type == "t"]
        if rocky:
            break
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        made = {p.name: p for p in rocky}
        for planet in loaded.planets:
            if planet.body_type == "t" and planet.name in made:
                assert planet.lid_regime == made[planet.name].lid_regime in ("mobile", "stagnant")
            elif planet.body_type == "g":
                assert planet.lid_regime is None
    finally:
        conn.rollback()
        conn.close()
