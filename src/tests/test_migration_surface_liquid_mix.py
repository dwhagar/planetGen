# tests/test_migration_surface_liquid_solutes.py

"""
GEN.1218b: migration 0200 gives a v199 galaxy the makeup, salt and second-liquid columns of the surface liquid, empty, on both bodies tables.
"""

import pytest

from planetgen.db import store
from planetgen.physics import surface_liquid_world
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

NEW = ("surface_liquid_methane_mol", "surface_liquid_ethane_mol", "surface_liquid_nitrogen_mol", "surface_liquid_salt_wt", "surface_liquid_2", "surface_liquid_2_fraction")


def test_a_v199_galaxy_gets_the_columns_empty_after_the_ammonia(mysql_config):
    load_old_schema(mysql_config, 199)
    assert store.schema_status(mysql_config)[0] == 199
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
            start = names.index("surface_liquid_ammonia_wt")
            assert names[start + 1:start + 7] == list(NEW), table
            assert set(NEW) <= set(surface_liquid_world.COLUMN_FIELDS)
    finally:
        conn.close()
