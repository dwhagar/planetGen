# tests/test_migration_rogue_habitability.py

"""
A saved free-floating brown dwarf keeps no habitability score (every field `None`), and migration 0158 leaves a v157 galaxy's brown dwarf row
without one (GEN.1111: the rogue planets that carried a score are now built by the planet chain and stored by the planet columns).
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.physics import rogue_habitability
from planetgen.util import draw
from tests import db_schema_support

pytestmark = pytest.mark.fresh_database


def _brown_dwarf():
    with draw.bound(3):
        return RoguePlanet(SystemConfig())


def test_a_saved_brown_dwarf_keeps_no_score(mysql_config):
    planet = _brown_dwarf()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = store.insert_rogue_planet(conn, planet)
        row = conn.execute("SELECT * FROM rogue_planets WHERE id = ?", (row_id,)).fetchone()
        assert all(row[name] is None for name in rogue_habitability.ROGUE_FIELDS)
    finally:
        conn.rollback()
        conn.close()


def test_a_v157_galaxy_migrates_with_its_brown_dwarf_unscored(mysql_config):
    planet = _brown_dwarf()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = store.insert_rogue_planet(conn, planet)
    finally:
        conn.close()
    db_schema_support.downgrade_to(mysql_config, 157)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT * FROM rogue_planets WHERE id = ?", (row_id,)).fetchone()
        assert all(row[name] is None for name in rogue_habitability.ROGUE_FIELDS)
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION  # a second run changes nothing
