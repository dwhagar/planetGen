# tests/test_migration_sector_leaf_digest.py

"""
DB.49: migration 0120 adds `sectors.leaf_digest` and makes `sector_tree_nodes`, and touches no row. A migration run again
leaves the column, the table and their rows as they are, and the table holds one band node per shell and band and one shell
node per shell.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

SECTOR = ("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES ('Aldebaran', 4000, 3, 0, 1)")


def _digest_column(conn):
    return conn.execute("SELECT COUNT(*) AS n FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
                        " AND TABLE_NAME = 'sectors' AND COLUMN_NAME = 'leaf_digest'").fetchone()["n"]


def _node(conn, *values):
    conn.execute("INSERT INTO sector_tree_nodes (level, shell_index, band_index) VALUES (%s, %s, %s)", values)


def test_a_v119_galaxy_gets_the_column_and_an_empty_table_and_keeps_its_sectors(mysql_config):
    load_old_schema(mysql_config, 119)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute(SECTOR)
        conn.commit()
        assert _digest_column(conn) == 0
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert store.schema_status(mysql_config)[0] >= 120
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _digest_column(conn) == 1
        assert conn.execute("SELECT name, leaf_digest FROM sectors").fetchall() == [{"name": "Aldebaran", "leaf_digest": None}]
        assert conn.execute("SELECT COUNT(*) AS n FROM sector_tree_nodes").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_column_the_table_and_their_rows(mysql_config):
    load_old_schema(mysql_config, 119)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute(SECTOR)
        conn.execute("UPDATE sectors SET leaf_digest = UNHEX(REPEAT('ab', 32))")
        _node(conn, 1, 3, 0)
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0119')")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT HEX(leaf_digest) AS d FROM sectors").fetchone()["d"] == "AB" * 32
        assert conn.execute("SELECT COUNT(*) AS n FROM sector_tree_nodes").fetchone()["n"] == 1
    finally:
        conn.close()


def test_a_node_is_one_row_per_level_shell_and_band_and_starts_dirty_with_no_digest(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _node(conn, 1, 3, -2)
        _node(conn, 1, 3, 0)
        _node(conn, 2, 3, 0)
        with pytest.raises(Exception, match="Duplicate entry"):
            _node(conn, 1, 3, 0)
        row = conn.execute("SELECT digest, child_count, dirty FROM sector_tree_nodes WHERE level = 2").fetchone()
        assert row == {"digest": None, "child_count": 0, "dirty": 1}
    finally:
        conn.close()


def test_a_shell_node_has_no_band_and_a_level_is_one_or_two(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        with pytest.raises(Exception, match="chk_sector_tree_nodes_shell_band"):
            _node(conn, 2, 3, 4)
        with pytest.raises(Exception, match="chk_sector_tree_nodes_level"):
            _node(conn, 3, 3, 0)
        with pytest.raises(Exception, match="chk_sector_tree_nodes_dirty"):
            conn.execute("INSERT INTO sector_tree_nodes (level, shell_index, dirty) VALUES (2, 1, 2)")
    finally:
        conn.close()


def test_the_digest_holds_thirty_two_bytes(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("INSERT INTO sector_tree_nodes (level, shell_index, band_index, digest) VALUES (1, 1, 1, %s)", (b"\x01" * 32,))
        assert conn.execute("SELECT digest FROM sector_tree_nodes").fetchone()["digest"] == b"\x01" * 32
        with pytest.raises(Exception, match="Data too long"):
            conn.execute("INSERT INTO sector_tree_nodes (level, shell_index, band_index, digest) VALUES (1, 1, 2, %s)", (b"\x01" * 33,))
    finally:
        conn.close()
