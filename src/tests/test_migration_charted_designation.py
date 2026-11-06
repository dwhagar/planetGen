# tests/test_migration_charted_designation.py

"""GEN.1257: migration 0221 lets `designations.reason` be `charted`; a database from before it keeps its rows and refuses any other new reason."""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _insert(conn, reason):
    conn.execute("INSERT INTO designations (designation_kind, designation_name, object_kind, object_row_id, retired_epoch, reason)"
                 " VALUES ('name', 'Banhunm Finkap', 'star_systems', 1, 1, ?)", (reason,))


def test_a_v220_database_takes_the_charted_reason_once_migrated(mysql_config):
    load_old_schema(mysql_config, 220)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        with pytest.raises(Exception):
            _insert(conn, "charted")
        conn.rollback()
        _insert(conn, "rename")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM designations WHERE reason = 'rename'").fetchone()["n"] == 1
        _insert(conn, "charted")
        conn.commit()
        with pytest.raises(Exception):
            _insert(conn, "invented")
        conn.rollback()
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
