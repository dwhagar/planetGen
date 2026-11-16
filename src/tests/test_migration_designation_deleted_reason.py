# tests/test_migration_designation_deleted_reason.py

"""GEN.1141: migration 0265 lets `designations.reason` be `deleted`; a database from before it keeps its rows and refuses any other new reason."""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _insert(conn, reason):
    conn.execute("INSERT INTO designations (designation_kind, designation_name, object_kind, object_row_id, retired_epoch, reason)"
                 " VALUES ('name', 'Gone Name', 'star_systems', 1, 1, ?)", (reason,))


def test_a_v264_database_takes_the_deleted_reason_once_migrated(mysql_config):
    load_old_schema(mysql_config, 264)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        with pytest.raises(Exception):
            _insert(conn, "deleted")
        conn.rollback()
        _insert(conn, "rename")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _insert(conn, "deleted")
        conn.commit()
        assert {row["reason"] for row in conn.execute("SELECT reason FROM designations").fetchall()} == {"rename", "deleted"}
        with pytest.raises(Exception):
            _insert(conn, "invented")
        conn.rollback()
    finally:
        conn.close()
