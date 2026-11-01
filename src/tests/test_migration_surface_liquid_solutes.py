# tests/test_migration_surface_liquid_solutes.py

"""
GEN.1215: migration 0197 gives a v196 galaxy the acid strength, the acid place and the ammonia columns of the surface liquid, empty, on both bodies tables.
"""

import pytest

from planetgen.db import store
from planetgen.physics import surface_liquid_world
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

NEW = ("surface_liquid_acid_wt", "surface_liquid_acid_layer", "surface_liquid_ammonia_wt")


def test_a_v196_galaxy_gets_the_columns_empty_after_the_depth(mysql_config):
    load_old_schema(mysql_config, 196)
    assert store.schema_status(mysql_config)[0] == 196
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            for column in NEW:
                assert not conn.execute("SHOW COLUMNS FROM " + table + " LIKE ?", (column,)).fetchall()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            names = [row["Field"] for row in conn.execute("SHOW COLUMNS FROM " + table).fetchall()]
            start = names.index("surface_liquid_depth_km")
            assert names[start + 1:start + 4] == list(NEW), table
            assert set(NEW) <= set(surface_liquid_world.COLUMN_FIELDS)
    finally:
        conn.close()
