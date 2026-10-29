# tests/test_migration_atmosphere_colour.py

"""
DB.96: the atmosphere colour of a planet or moon. A system saved and loaded keeps each body's colour; migration 0162 gives a galaxy made at
v161 the column, empty, and a body without it reads back as having none.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import body_colour
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
        if any(getattr(p, "atmosphere_colour", None) for p in system.planets):
            return system, cfg
    pytest.fail("no system with a planet with an atmosphere colour in 60 seeds")


def _colours_of(system):
    return [tuple(getattr(body, name, None) for name in body_colour.BODY_FIELDS) for body in system.planets if body.body_type != "a"]


def test_a_saved_system_keeps_every_planets_atmosphere_colour(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        before, after = _colours_of(system), _colours_of(loaded)
        assert any(row[0] is not None for row in before)
        assert after == before
    finally:
        conn.rollback()
        conn.close()


def test_a_v161_galaxy_gets_the_columns_empty(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
    finally:
        conn.close()
    db_schema_support.downgrade_to(mysql_config, 161)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = conn.execute("SELECT atmosphere_colour FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()
        assert rows and all(row["atmosphere_colour"] is None for row in rows)
        loaded = store.load_star_system(conn, system_id)
        assert all(row == (None,) * len(body_colour.BODY_FIELDS) for row in _colours_of(loaded))
    finally:
        conn.close()
