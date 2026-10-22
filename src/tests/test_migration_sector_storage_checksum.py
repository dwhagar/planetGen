# tests/test_migration_sector_storage_checksum.py

"""
DB.44: migration 0115 adds `sectors.storage_checksum` (BINARY(32), NULL) and `sectors.parity_dirty` (TINYINT NOT NULL,
default 1) in one `ALGORITHM=COPY` rebuild. A sector saved before the revision has no checksum and counts as not yet
protected, a migration run again leaves the columns as they are, and the columns hold what is written to them.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _column(conn, name):
    return conn.execute(
        "SELECT column_type AS t, is_nullable AS n, column_default AS d FROM information_schema.COLUMNS"
        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'sectors' AND COLUMN_NAME = ?", (name,)).fetchone()


def test_a_v114_sector_has_no_checksum_and_is_not_protected_after_the_migration(mysql_config):
    load_old_schema(mysql_config, 114)
    assert store.schema_status(mysql_config)[0] == 114
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _column(conn, "parity_dirty") is None
        conn.execute("INSERT INTO sectors (name, edge_mpc) VALUES ('old', 4.0)")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT storage_checksum, parity_dirty FROM sectors WHERE name = 'old'").fetchone()
        assert row["storage_checksum"] is None
        assert row["parity_dirty"] == 1
        assert _column(conn, "storage_checksum")["t"] == "binary(32)"
        assert _column(conn, "storage_checksum")["n"] == "YES"
        assert _column(conn, "parity_dirty")["n"] == "NO"
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_columns_as_they_are(mysql_config):
    load_old_schema(mysql_config, 114)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        digest = bytes(range(32))
        conn.execute("INSERT INTO sectors (name, edge_mpc, storage_checksum, parity_dirty) VALUES ('kept', 4.0, ?, 0)",
                     (digest,))
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0114')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT storage_checksum, parity_dirty FROM sectors WHERE name = 'kept'").fetchone()
        assert bytes(row["storage_checksum"]) == bytes(range(32))
        assert row["parity_dirty"] == 0
    finally:
        conn.close()
