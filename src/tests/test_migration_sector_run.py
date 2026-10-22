# tests/test_migration_sector_run.py

"""
DB.16: migration 0121 gives `sectors` a `generator_epoch` and a `run_id` in place of `version_key`, `planetgen_version`,
`python_version` and `platform`. Each distinct set of the old values becomes one `generation_runs` row that its sectors point
at, a sector with no recorded version keeps an empty run, and a migration stopped half way starts again where it stopped.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

OLD = ("version_key", "planetgen_version", "python_version", "platform")
KEY_A = "0008000000030F030D1000"
KEY_B = "0007007F000160030C0300"


def _columns(conn):
    return {row["c"] for row in conn.execute(
        "SELECT COLUMN_NAME AS c FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'sectors'").fetchall()}


def _old_sectors(conn):
    """Four sectors: two from one release, one from another, one with no version recorded."""
    rows = [("a1", KEY_A, "8.0.3", "3.13.5", "Linux x86_64"), ("a2", KEY_A, "8.0.3", "3.13.5", "Linux x86_64"),
            ("b1", KEY_B, "7.127.352", "3.12.3", "Linux x86_64"), ("none", None, None, None, None)]
    for name, key, release, python, platform in rows:
        conn.execute("INSERT INTO sectors (name, edge_mpc, version_key, planetgen_version, python_version, platform)"
                     " VALUES (?, 4.0, ?, ?, ?, ?)", (name, key, release, python, platform))
    conn.commit()


def _by_name(conn):
    return {row["name"]: row for row in conn.execute(
        "SELECT s.name, s.generator_epoch, s.run_id, r.command, r.version_key, r.planetgen_version, r.python_version,"
        " r.platform, r.outcome FROM sectors s LEFT JOIN generation_runs r ON r.id = s.run_id").fetchall()}


def test_the_sectors_point_at_a_run_for_each_set_of_old_values(mysql_config):
    load_old_schema(mysql_config, 120)
    assert store.schema_status(mysql_config)[0] == 120
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _old_sectors(conn)
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        columns = _columns(conn)
        assert {"generator_epoch", "run_id"} <= columns and not columns & set(OLD)
        rows = _by_name(conn)
        assert rows["a1"]["run_id"] == rows["a2"]["run_id"] != rows["b1"]["run_id"]
        assert (rows["a1"]["version_key"], rows["a1"]["planetgen_version"], rows["a1"]["python_version"],
                rows["a1"]["platform"]) == (KEY_A, "8.0.3", "3.13.5", "Linux x86_64")
        assert rows["b1"]["version_key"] == KEY_B and rows["b1"]["planetgen_version"] == "7.127.352"
        assert rows["a1"]["command"] == "migrated" and rows["a1"]["outcome"] == "ok"
        assert rows["none"]["run_id"] is None
        assert {row["generator_epoch"] for row in rows.values()} == {None}
        assert conn.execute("SELECT COUNT(*) AS n FROM generation_runs").fetchone()["n"] == 2
        assert not conn.execute("SELECT 1 FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = DATABASE()"
                                " AND TABLE_NAME = 'sectors' AND INDEX_NAME = 'idx_sectors_version_key'").fetchall()
    finally:
        conn.close()


def test_a_migration_stopped_after_the_new_columns_goes_on_from_there(mysql_config):
    load_old_schema(mysql_config, 120)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _old_sectors(conn)
        conn.execute("ALTER TABLE sectors ADD COLUMN generator_epoch SMALLINT UNSIGNED NULL,"
                     " ADD COLUMN run_id BIGINT UNSIGNED NULL, ALGORITHM=COPY")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = _by_name(conn)
        assert rows["a1"]["version_key"] == KEY_A and rows["b1"]["version_key"] == KEY_B
        assert not _columns(conn) & set(OLD)
        assert conn.execute("SELECT COUNT(*) AS n FROM generation_runs").fetchone()["n"] == 2
    finally:
        conn.close()


def test_a_migration_run_again_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 120)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _old_sectors(conn)
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        before = _by_name(conn)
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0120')")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _by_name(conn) == before
        assert conn.execute("SELECT COUNT(*) AS n FROM generation_runs").fetchone()["n"] == 2
    finally:
        conn.close()
