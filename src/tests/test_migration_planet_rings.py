# tests/test_migration_planet_rings.py

"""
GEN.732: the ring system of a planet. A system saved and loaded keeps each planet's ring; migration 0160 gives a galaxy made at v159 the
columns, empty, and a body without them reads back as having no ring.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import rings
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
        if any(getattr(p, "ring_kind", None) for p in system.planets):
            return system, cfg
    pytest.fail("no system with a ringed planet in 60 seeds")


def _rings_of(system):
    return [tuple(getattr(body, name, None) for name in rings.BODY_FIELDS) for body in system.planets if body.body_type != "a"]


def test_a_saved_system_keeps_every_planets_ring(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before, after = _rings_of(system), _rings_of(loaded)
        assert any(row[0] is not None for row in before)
        assert after == before
    finally:
        conn.rollback()
        conn.close()


def test_a_ring_kind_the_check_does_not_know_is_refused(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        with pytest.raises(Exception):
            conn.execute("UPDATE planets SET ring_kind = 'granite' WHERE star_system_id = ? LIMIT 1", (system_id,))
        conn.rollback()
        with pytest.raises(Exception):
            conn.execute("UPDATE planets SET ring_origin = 'granite' WHERE star_system_id = ? LIMIT 1", (system_id,))
    finally:
        conn.rollback()
        conn.close()


def test_a_v159_galaxy_gets_the_columns_empty(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
    finally:
        conn.close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        db_schema_support.drop_rings(conn)
        conn.execute("UPDATE alembic_version SET version_num = '0159'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (159)")
        conn.commit()
        assert store._stored_version(conn, "schema_migrations") == 159
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = conn.execute("SELECT ring_kind, ring_mass_kg FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()
        assert rows and all(row["ring_kind"] is None and row["ring_mass_kg"] is None for row in rows)
        loaded = store.load_star_system(conn, system_id)
        assert all(row == (None,) * len(rings.BODY_FIELDS) for row in _rings_of(loaded))
    finally:
        conn.close()
