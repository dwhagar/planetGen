# tests/test_migration_bulk_composition.py

"""
DB.69, GEN.428: the bulk composition columns. A system saved and loaded keeps every body's composition; migration 0144 gives a galaxy made at
v143 the columns, empty, and a body without them reads back as having no composition.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import bulk_composition as bc
from planetgen.util import draw
from tests import db_schema_support

pytestmark = pytest.mark.fresh_database


def _system():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(p.moons for p in system.planets) and any(p.body_type == "g" for p in system.planets):
            return system, cfg
    pytest.fail("no system with a giant and a moon in 60 seeds")


def _compositions(system):
    return [bc.bulk_fractions(body) for planet in system.planets if hasattr(planet, "moons") for body in [planet, *planet.moons]]


def test_a_saved_system_keeps_every_bodys_composition(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before, after = _compositions(system), _compositions(loaded)
        assert before and all(item is not None for item in before)
        assert after == before
        row = conn.execute("SELECT carbon_rock, composition_origin FROM moons WHERE star_system_id = ? LIMIT 1", (system_id,)).fetchone()
        assert row["carbon_rock"] in (0, 1) and row["composition_origin"] in bc.ORIGINS
    finally:
        conn.rollback()
        conn.close()


def test_a_v143_galaxy_gets_the_columns_empty(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
    finally:
        conn.close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        db_schema_support.drop_bulk_composition(conn)
        conn.execute("UPDATE alembic_version SET version_num = '0143'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (143)")
        conn.commit()
        assert store._stored_version(conn, "schema_migrations") == 143
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = conn.execute("SELECT metal_mass_fraction, composition_origin, water_mass_fraction FROM planets WHERE star_system_id = ?",
                            (system_id,)).fetchall()
        assert rows and all(row["metal_mass_fraction"] is None and row["composition_origin"] is None for row in rows)
        loaded = store.load_star_system(conn, system_id)
        assert all(item is None for item in _compositions(loaded))
    finally:
        conn.close()
