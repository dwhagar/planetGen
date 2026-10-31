# tests/test_migration_moon_tidal_flux.py

"""
A moon keeps the tide heat flow its ocean's life counts: a saved moon carries `tidal_flux_w_m2`, and migration 0152 adds the
column to a v151 galaxy with NULL for every moon it already holds.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests import db_schema_support

pytestmark = pytest.mark.fresh_database


def _system_with_a_tide():
    for seed in range(120):
        cfg = SystemConfig()
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(getattr(m, "tidal_flux_w_m2", None) for p in system.planets if p.body_type != "a" for m in p.moons):
            return system, cfg
    pytest.skip("no icy moon under ice in 120 seeds")


def test_a_saved_moon_keeps_its_tide(mysql_config):
    system, cfg = _system_with_a_tide()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        stored = [row["tidal_flux_w_m2"] for row in conn.execute(
            "SELECT tidal_flux_w_m2 FROM moons WHERE star_system_id = ?", (system_id,)).fetchall()]
        assert any(value for value in stored)
        assert all(value is None or value >= 0 for value in stored)   # GEN.778: a circular orbit stores 0.0
        loaded = store.load_star_system(conn, system_id)
        original = sorted(m.tidal_flux_w_m2 or 0.0 for p in system.planets if p.body_type != "a" for m in p.moons)
        again = sorted(m.tidal_flux_w_m2 or 0.0 for p in loaded.planets if p.body_type != "a" for m in p.moons)
        assert again == pytest.approx(original)
    finally:
        conn.rollback()
        conn.close()


def test_a_v151_galaxy_gets_the_column_and_every_moon_keeps_its_row(mysql_config):
    system, cfg = _system_with_a_tide()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        count = conn.execute("SELECT COUNT(*) AS n FROM moons WHERE star_system_id = ?", (system_id,)).fetchone()["n"]
    finally:
        conn.close()
    db_schema_support.downgrade_to(mysql_config, 151)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = conn.execute("SELECT tidal_flux_w_m2 FROM moons WHERE star_system_id = ?", (system_id,)).fetchall()
        assert len(rows) == count and all(row["tidal_flux_w_m2"] is None for row in rows)
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION  # a second run changes nothing
