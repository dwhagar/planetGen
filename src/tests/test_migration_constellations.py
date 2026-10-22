# tests/test_migration_constellations.py

"""
DB.64: migration 0128 makes `constellations` and touches nothing else. A migration run again leaves the table and its rows as they
are, and a system keeps one row for each constellation number.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

UID = b"\x01" * 12


def _tables(conn):
    return {row["t"] for row in conn.execute(
        "SELECT TABLE_NAME AS t FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()").fetchall()}


def _constellation(conn, index=0, system=UID):
    return conn.execute(
        "INSERT INTO constellations (system_uid, constellation_index, members, edges, seed, version_key)"
        " VALUES (?, ?, '[1, 2]', '[[1, 2]]', 7, 'AAAAAAAAAAAAAAAAAAAAAA')", (system, index)).lastrowid


def test_a_v127_galaxy_gets_an_empty_table(mysql_config):
    load_old_schema(mysql_config, 127)
    assert store.schema_status(mysql_config)[0] == 127
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "constellations" in _tables(conn)
        assert conn.execute("SELECT COUNT(*) AS n FROM constellations").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_table_and_its_rows(mysql_config):
    load_old_schema(mysql_config, 127)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _constellation(conn)
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0127')")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT constellation_index AS i FROM constellations").fetchall() == [{"i": 0}]
    finally:
        conn.close()


def test_a_system_has_one_row_for_each_constellation_number(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _constellation(conn, 0)
        _constellation(conn, 1)
        _constellation(conn, 0, system=b"\x02" * 12)
        with pytest.raises(Exception, match="uq_constellations_system_index"):
            _constellation(conn, 0)
    finally:
        conn.close()
