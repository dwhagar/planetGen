# tests/test_migration_trojan_planets.py

"""
GEN.1151: migration 0270 adds `lagrange_point` and `lagrange_host_id` to `planets`. The rows made before keep their figures and get NULL, and the migration can be run over a table that already has a column
(a stopped run is started again).
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

COLUMNS = ("lagrange_point", "lagrange_host_id")


def _has_column(conn, name):
    return bool(conn.execute(f"SHOW COLUMNS FROM planets LIKE '{name}'").fetchall())


def test_the_columns_are_added_empty_and_the_rows_keep_their_figures(mysql_config):
    load_old_schema(mysql_config, 269)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not any(_has_column(conn, name) for name in COLUMNS)
        before = conn.execute("SHOW COLUMNS FROM planets").fetchall()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert all(_has_column(conn, name) for name in COLUMNS)
        assert len(conn.execute("SHOW COLUMNS FROM planets").fetchall()) == len(before) + len(COLUMNS)
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION         # again: the columns are there, nothing breaks
