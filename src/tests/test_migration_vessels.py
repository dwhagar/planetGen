# tests/test_migration_vessels.py

"""
DB.63: migration 0129 makes `vessels` and touches nothing else. A migration run again leaves the table and its rows as they are, a
vessel is a ship or a station and was generated or made by a user, and an object ID names one vessel.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _tables(conn):
    return {row["t"] for row in conn.execute(
        "SELECT TABLE_NAME AS t FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()").fetchall()}


def _vessel(conn, uid=b"\x01" * 12, kind="ship", origin="generated"):
    return conn.execute("INSERT INTO vessels (uid, kind, name, origin) VALUES (?, ?, 'Test', ?)", (uid, kind, origin)).lastrowid


def test_a_v128_galaxy_gets_an_empty_table(mysql_config):
    load_old_schema(mysql_config, 128)
    assert store.schema_status(mysql_config)[0] == 128
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "vessels" in _tables(conn)
        assert conn.execute("SELECT COUNT(*) AS n FROM vessels").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_table_and_its_rows(mysql_config):
    load_old_schema(mysql_config, 128)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _vessel(conn)
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0128')")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT kind, owner_id, velocity_x_kms AS vx FROM vessels").fetchall() == [
            {"kind": "ship", "owner_id": None, "vx": 0.0}]
    finally:
        conn.close()


def test_a_vessel_has_a_known_kind_and_origin_and_a_unique_object_id(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _vessel(conn, b"\x01" * 12, "ship", "generated")
        _vessel(conn, b"\x02" * 12, "station", "user")
        with pytest.raises(Exception, match="chk_vessels_kind"):
            _vessel(conn, b"\x03" * 12, kind="probe")
        with pytest.raises(Exception, match="chk_vessels_origin"):
            _vessel(conn, b"\x04" * 12, origin="imported")
        with pytest.raises(Exception, match="uq_vessels_uid"):
            _vessel(conn, b"\x01" * 12)
    finally:
        conn.close()
