# tests/test_migration_object_locks.py

"""
DB.37: migration 0123 makes `object_locks` and touches nothing else. A migration run again leaves the table and its rows as
they are, and an object holds at most one lock per aspect (`*` being the whole object).
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

UID = b"\x80" * 12


def _exists(conn):
    return conn.execute("SELECT COUNT(*) AS n FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
                        " AND TABLE_NAME = 'object_locks'").fetchone()["n"]


def _lock(conn, uid=UID, **more):
    columns = {"object_uid": uid, "locked_by_name": "Boss", **more}
    conn.execute(f"INSERT INTO object_locks ({', '.join(columns)}) VALUES ({', '.join(['%s'] * len(columns))})",
                 tuple(columns.values()))


def test_a_v122_galaxy_gets_an_empty_table_and_keeps_its_rows(mysql_config):
    load_old_schema(mysql_config, 122)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _exists(conn) == 0
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _exists(conn) == 1
        assert conn.execute("SELECT COUNT(*) AS n FROM object_locks").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_table_and_its_rows(mysql_config):
    load_old_schema(mysql_config, 122)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _lock(conn, reason="the capital")
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0122')")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT aspect, reason FROM object_locks").fetchall() == [{"aspect": "*", "reason": "the capital"}]
    finally:
        conn.close()


def test_an_object_has_one_lock_per_aspect(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _lock(conn)
        _lock(conn, aspect="class")
        _lock(conn, aspect="orbit")
        _lock(conn, uid=b"\x81" * 12, aspect="class")
        with pytest.raises(Exception, match="Duplicate entry"):
            _lock(conn, aspect="class")
        with pytest.raises(Exception, match="Duplicate entry"):
            _lock(conn)
        assert conn.execute("SELECT COUNT(*) AS n FROM object_locks WHERE object_uid = %s", (UID,)).fetchone()["n"] == 3
        row = conn.execute("SELECT locked_by_id, locked_at IS NOT NULL AS stamped FROM object_locks WHERE aspect = 'orbit'").fetchone()
        assert row == {"locked_by_id": None, "stamped": 1}
    finally:
        conn.close()
