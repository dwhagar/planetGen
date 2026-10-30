# tests/test_migration_ice_species.py

"""
GEN.834: migration 0179 gives a v178 galaxy the eight ice species columns, empty, on both bodies tables.
"""

import pytest

from planetgen.db import store
from planetgen.physics import ice_species
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def test_a_v178_galaxy_gets_the_columns_empty(mysql_config):
    load_old_schema(mysql_config, 178)
    assert store.schema_status(mysql_config)[0] == 178
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            assert not conn.execute("SHOW COLUMNS FROM " + table + " LIKE 'ice_co2_mass_fraction'").fetchall()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            for column in ice_species.BODY_FIELDS:
                assert conn.execute("SHOW COLUMNS FROM " + table + " LIKE ?", (column,)).fetchall(), (table, column)
    finally:
        conn.close()
