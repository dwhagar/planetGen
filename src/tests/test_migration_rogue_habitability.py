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
