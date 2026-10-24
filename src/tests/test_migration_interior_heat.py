# tests/test_migration_interior_heat.py

"""
GEN.699: migration 0153 adds the mantle's two draws and the internal heat flux to `planets` and `moons`. A body from before v153 keeps
its rows with NULL there (and the Urey ratio it was drawn with); a saved system keeps its bodies' heat.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

COLUMNS = ("mantle_eta0_pa_s", "mantle_activation_kj_mol", "internal_heat_flux_w_m2")


def _columns(conn, table):
    return {row["COLUMN_NAME"]: row["DATA_TYPE"] for row in conn.execute(
        "SELECT COLUMN_NAME, DATA_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", (table,)).fetchall()}


def test_a_v152_galaxy_gets_the_columns_empty(mysql_config):
    load_old_schema(mysql_config, 152)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert not set(COLUMNS) & set(_columns(conn, table))
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            found = _columns(conn, table)
            assert all(found[column] == "double" for column in COLUMNS), table
    finally:
        conn.close()


def test_a_saved_system_keeps_every_bodys_heat(mysql_config):
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    system = None
    for seed in range(40):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(p.body_type == "t" for p in system.planets):
            break
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        made = {p.name: p for p in system.planets if p.body_type == "t"}
        assert made
        for planet in loaded.planets:
            if planet.body_type == "t" and planet.name in made:
                before = made[planet.name]
                for name in (*COLUMNS, "urey_ratio"):
                    assert getattr(planet, name) == pytest.approx(getattr(before, name))
                assert planet.internal_heat_flux_w_m2 > 0.0 and 0.0 < planet.urey_ratio < 1.5
    finally:
        conn.rollback()
        conn.close()
