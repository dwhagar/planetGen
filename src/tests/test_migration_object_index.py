# tests/test_migration_object_index.py

"""
DB.79: migration 0113 makes `object_index` and touches nothing else. A migration run again leaves the table as it is, and the
table holds one row per kind and object ID, finds a name with its FULLTEXT index and reads a page after a row from a keyset
key.
"""

import importlib

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

migration = importlib.import_module("planetgen.db.migrations.versions.0113_object_index")


def _insert(conn, kind, uid, name, **more):
    columns = {"kind": kind, "uid": uid, "source_id": 1, "name": name, **more}
    conn.execute(f"INSERT INTO object_index ({', '.join(columns)}) VALUES ({', '.join(['%s'] * len(columns))})",
                 tuple(columns.values()))


def test_a_v112_galaxy_gets_an_empty_index_and_keeps_its_rows(mysql_config):
    load_old_schema(mysql_config, 112)
    assert store.schema_status(mysql_config)[0] == 112
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
                            " AND TABLE_NAME = 'object_index'").fetchone()["n"] == 0
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert store.schema_status(mysql_config)[0] >= 113
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM object_index").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_table_and_its_rows(mysql_config):
    load_old_schema(mysql_config, 112)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _insert(conn, "planet", b"\x80" * 12, "Kepler")
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0112')")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT name FROM object_index").fetchall() == [{"name": "Kepler"}]
    finally:
        conn.close()


def test_an_object_is_one_row_per_kind_and_object_id(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _insert(conn, "planet", b"\x80" * 12, "Kepler b")
        _insert(conn, "moon", b"\x80" * 12, "Kepler b I")
        with pytest.raises(Exception, match="Duplicate entry"):
            _insert(conn, "planet", b"\x80" * 12, "Kepler b again")
        with pytest.raises(Exception, match="chk_object_index_placed"):
            _insert(conn, "planet", b"\x81" * 12, "Odd", placed=2)
    finally:
        conn.close()


def test_the_name_index_finds_a_word_and_a_keyset_key_reads_the_next_page(mysql_config):
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for number in range(5):
            _insert(conn, "planet", bytes([0x80 + number]) * 12, f"Aldebaran {chr(97 + number)}", mass_kg=number)
        _insert(conn, "planet", b"\x90" * 12, "Sirius b")
        conn.commit()
        found = conn.execute("SELECT name FROM object_index WHERE MATCH(name) AGAINST (%s IN BOOLEAN MODE)"
                             " ORDER BY name", ("Aldebaran",)).fetchall()
        assert [row["name"] for row in found] == [f"Aldebaran {chr(97 + number)}" for number in range(5)]
        page = conn.execute("SELECT name FROM object_index WHERE kind = 'planet' AND (name, uid) > (%s, %s)"
                            " ORDER BY name, uid LIMIT 2", ("Aldebaran b", bytes([0x81]) * 12)).fetchall()
        assert [row["name"] for row in page] == ["Aldebaran c", "Aldebaran d"]
    finally:
        conn.close()


def test_the_revision_is_the_head_after_this_one():
    assert migration.revision == "0113" and migration.down_revision == "0112"
