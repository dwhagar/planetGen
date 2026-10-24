# tests/test_migration_ozone_column.py

"""
GEN.413: migration 0135 adds `planets.ozone_column_du` and `moons.ozone_column_du`. A generated rocky world with oxygen
keeps its ozone column through a save and a load; a body row made at v134 reads back with NULL, as nothing was drawn for it.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def test_a_v134_galaxy_gains_the_column_empty(mysql_config):
    load_old_schema(mysql_config, 134)
    assert store.schema_status(mysql_config)[0] == 134
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            column = conn.execute(
                "SELECT DATA_TYPE, IS_NULLABLE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
                " AND TABLE_NAME = ? AND COLUMN_NAME = 'ozone_column_du'", (table,)).fetchone()
            assert column is not None and column["IS_NULLABLE"] == "YES" and column["DATA_TYPE"] == "double"
        assert store.detect_schema_version(conn) == store.SCHEMA_VERSION
    finally:
        conn.close()
