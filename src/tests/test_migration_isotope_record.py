# tests/test_migration_isotope_record.py

"""
GEN.1377, GEN.1420: the isotope record as child table rows. A system saved and loaded keeps every body's record; migration 0235 gives a v234 galaxy the three
tables (empty: nothing is backfilled), and a second run changes nothing.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import isotope_record
from planetgen.util import draw
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

TABLES = ("planet_isotopes", "moon_isotopes", "rogue_planet_isotopes")


def _system():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(getattr(p, "moons", None) for p in system.planets):
            return system, cfg
    pytest.fail("no system with a moon in 60 seeds")


def _records(system):
    return [getattr(body, isotope_record.ATTRIBUTE) for planet in system.planets if hasattr(planet, "moons") for body in [planet, *planet.moons]]


def test_a_saved_system_keeps_every_bodys_record(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before, after = _records(system), _records(loaded)
        assert any(before) and after == before
        stored = sum(conn.execute(f"SELECT COUNT(*) AS n FROM {kind}_isotopes i JOIN {kind}s b ON b.id = i.{kind}_id WHERE b.star_system_id = ?",
                                  (system_id,)).fetchone()["n"] for kind in ("planet", "moon"))
        assert stored == sum(len(record or ()) for record in before)
    finally:
        conn.rollback()
        conn.close()


def test_a_v234_galaxy_gets_the_three_tables_empty_and_a_second_run_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 234)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        existing = {next(iter(row.values())) for row in conn.execute("SHOW TABLES").fetchall()}
        assert not set(TABLES) & existing
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in TABLES:
            columns = [row["Field"] for row in conn.execute(f"SHOW COLUMNS FROM {table}").fetchall()]
            assert columns[1:] == [table.replace("_isotopes", "") + "_id", "element", "isotope", "fraction", "initial_ratio", "origin"]
            assert conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"] == 0
    finally:
        conn.close()
