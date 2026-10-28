# tests/test_migration_surface_composition.py

"""
GEN.441: the surface composition columns. A system saved and loaded keeps every rocky body's ground; migration 0159 gives a galaxy made at
v151 the columns, empty, and a body without them reads back as having no ground.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import surface_composition as sc
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


def _grounds(system):
    return [{field: getattr(body, field) for field in sc.COLUMN_FIELDS}
            for planet in system.planets if hasattr(planet, "moons") for body in [planet, *planet.moons]]


def test_a_saved_system_keeps_every_bodys_ground(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before, after = _grounds(system), _grounds(loaded)
        assert any(item["ground_type"] for item in before) and any(item["ground_type"] is None for item in before)
        assert after == before
        for item in before:
            if item["ground_type"]:
                assert sum(item[field] for field in sc.FRACTION_FIELDS) == pytest.approx(1.0)
    finally:
        conn.rollback()
        conn.close()


def test_a_v151_galaxy_gets_the_columns_empty(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
    finally:
        conn.close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        db_schema_support.drop_surface_composition(conn)
        conn.execute("UPDATE alembic_version SET version_num = '0151'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (151)")
        conn.commit()
        assert store._stored_version(conn, "schema_migrations") == 151
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = conn.execute("SELECT ground_type, ground_basaltic_fraction FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()
        assert rows and all(row["ground_type"] is None and row["ground_basaltic_fraction"] is None for row in rows)
        loaded = store.load_star_system(conn, system_id)
        assert all(item["ground_type"] is None for item in _grounds(loaded))
    finally:
        conn.close()
