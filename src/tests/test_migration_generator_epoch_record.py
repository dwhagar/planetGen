# tests/test_migration_generator_epoch_record.py

"""
OPS.58: migration 0131 adds `generator_epoch`, `fp_spec`, `battery_digest` and `battery_cases` to `galaxy_shape` and
`generator_epoch` to `generation_runs`. A galaxy planned before it keeps its rows with NULL in the new columns (the
epoch it was planned under was never recorded) and the version check says nothing of them; a migration run again
changes nothing.
"""

import pytest

from planetgen.db import store
from planetgen.galaxy import version_check
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

NEW = {"galaxy_shape": ("generator_epoch", "fp_spec", "battery_digest", "battery_cases"),
       "generation_runs": ("generator_epoch",)}


def _columns(conn, table):
    return {row["c"]: row["t"] for row in conn.execute(
        "SELECT COLUMN_NAME AS c, COLUMN_TYPE AS t FROM information_schema.COLUMNS"
        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", (table,)).fetchall()}


def test_a_v130_galaxy_gains_the_columns_and_keeps_its_rows(mysql_config):
    load_old_schema(mysql_config, 130)
    assert store.schema_status(mysql_config)[0] == 130
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table, names in NEW.items():
            assert not set(names) & set(_columns(conn, table))
        conn.execute(
            "INSERT INTO generation_runs (command, version_key, planetgen_version, python_version, platform)"
            " VALUES ('plan', ?, '8.0.1', '3.12.3', 'Linux x86_64')", ("0" * 22,))
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _columns(conn, "galaxy_shape")["generator_epoch"] == "smallint unsigned"
        assert _columns(conn, "galaxy_shape")["fp_spec"] == "smallint unsigned"
        assert _columns(conn, "galaxy_shape")["battery_digest"] == "char(64)"
        assert _columns(conn, "galaxy_shape")["battery_cases"] == "json"
        assert _columns(conn, "generation_runs")["generator_epoch"] == "smallint unsigned"
        assert conn.execute("SELECT command, generator_epoch FROM generation_runs").fetchall() == [
            {"command": "plan", "generator_epoch": None}]
        assert version_check.galaxy_check(conn, compute_battery=True) == ([], [])  # no galaxy, nothing recorded
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_columns_and_their_values(mysql_config):
    load_old_schema(mysql_config, 130)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        before = {table: _columns(conn, table) for table in NEW}
        conn.execute(
            "INSERT INTO generation_runs (command, version_key, planetgen_version, python_version, platform, generator_epoch)"
            " VALUES ('plan', ?, '8.0.1', '3.12.3', 'Linux x86_64', 55)", ("0" * 22,))
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0130')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert {table: _columns(conn, table) for table in NEW} == before
        assert conn.execute("SELECT generator_epoch FROM generation_runs").fetchall() == [{"generator_epoch": 55}]
    finally:
        conn.close()
