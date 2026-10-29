# tests/test_migration_rogue_rings.py

"""
GEN.836: the ring system of a rogue planet. A saved rogue planet keeps its ring; migration 0171 gives a v170 galaxy the columns, empty.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.physics import rings
from planetgen.util import draw
from tests import db_schema_support

pytestmark = pytest.mark.fresh_database


def _ringed_rogue():
    for seed in range(80):
        with draw.bound(seed):
            planet = RoguePlanet(SystemConfig(), mass_bin="saturn")
        if planet.ring_kind is not None:
            return planet
    pytest.fail("no ringed rogue planet in 80 seeds")


def _downgrade(config):
    db_schema_support.downgrade_to(config, 170)


def test_a_saved_rogue_planet_keeps_its_ring(mysql_config):
    planet = _ringed_rogue()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = store.insert_rogue_planet(conn, planet)
        row = conn.execute("SELECT * FROM rogue_planets WHERE id = ?", (row_id,)).fetchone()
        for name in rings.BODY_FIELDS:
            assert row[name] == getattr(planet, name)
    finally:
        conn.rollback()
        conn.close()


def test_a_v170_galaxy_gets_the_columns_empty(mysql_config):
    planet = _ringed_rogue()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = store.insert_rogue_planet(conn, planet)
    finally:
        conn.close()
    _downgrade(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT * FROM rogue_planets WHERE id = ?", (row_id,)).fetchone()
        assert all(row[name] is None for name in rings.BODY_FIELDS)
    finally:
        conn.close()
