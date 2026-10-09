# tests/test_db_columns_and_indexes.py

"""
DB.13: every stored value is in its own column, not a JSON block, and the
columns the searches filter and sort on are indexed.

- No galaxy or control table keeps a JSON or serialized blob: the text
  columns that remain are listed here, each prose a person reads (or a
  display summary of rows kept in child tables), and a new one has to be
  added to the list on purpose.
- `SEARCHED_COLUMNS` lists what `queryDb.search`, `list_*` and the API's
  filters and sorts use; each needs an index that starts with it.
- A run's command line and a job's command line are one row per argument.
"""

import json
import uuid

import pymysql
import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema, scratch_database

PROSE_TEXT_COLUMNS = {
    # Galaxy schema: generated prose and the summaries shown on pages.
    ("facilities", "description"), ("planets", "description"), ("planets", "atmosphere"),
    ("planets", "composition"), ("planets", "flavor_text"), ("moons", "description"), ("moons", "atmosphere"),
    ("moons", "composition"), ("moons", "flavor_text"),
    ("planet_evolutionary_paragraphs", "paragraph"), ("moon_evolutionary_paragraphs", "paragraph"),
    ("star_systems", "location"), ("star_systems", "system_flavor_text"),
    ("asteroid_belts", "composition_summary"), ("asteroid_fields", "composition_summary"),
    ("comets", "composition_summary"), ("interstellar_comets", "composition_summary"),
    ("rogue_planets", "composition"), ("nebulae", "composition"), ("nebulae", "formation_cause"),
    # Control schema: an error message and the audit log's note.
    ("work_tasks", "error"), ("admin_audit_log", "detail"),
}
BLOB_TYPES = ("tinytext", "text", "mediumtext", "longtext", "json", "tinyblob", "blob", "mediumblob", "longblob")

SEARCHED_COLUMNS = {
    "sectors": ("name", "galactic_radius_pc"),
    "star_systems": ("name", "sector_id", "quadrant"),
    "stars": ("name", "star_type", "yerkes_class", "radius_km"),
    "planets": ("name", "planet_class", "body_type", "life_chemical", "radius_km"),
    "moons": ("name", "planet_class", "body_type", "life_chemical", "radius_km"),
    "asteroid_belts": ("density",),
}


def _text_columns(config):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        return {(row["t"], row["c"]) for row in conn.execute(
            "SELECT table_name AS t, column_name AS c FROM information_schema.columns"
            " WHERE table_schema = DATABASE() AND data_type IN ({})".format(", ".join("?" * len(BLOB_TYPES))),
            BLOB_TYPES).fetchall()}
    finally:
        conn.close()


@pytest.fixture
def control_config(mysql_config):
    name = f"planetgen_test_ctl_{uuid.uuid4().hex[:12]}"
    admin = pymysql.connect(host=mysql_config.host, port=mysql_config.port, user=mysql_config.user,
                            password=mysql_config.password)
    admin.cursor().execute(f"CREATE DATABASE `{name}`")
    config = store.MySQLConfig(host=mysql_config.host, port=mysql_config.port, user=mysql_config.user,
                               password=mysql_config.password, database=name)
    store.get_control_connection(config, ensure_schema=True).close()
    try:
        yield config
    finally:
        store.close_pool(config)
        admin.cursor().execute(f"DROP DATABASE `{name}`")
        admin.close()


def test_no_galaxy_table_keeps_a_json_or_serialized_blob(mysql_config):
    store.get_connection(mysql_config).close()
    found = _text_columns(mysql_config)
    unexpected = sorted(found - PROSE_TEXT_COLUMNS)
    assert not unexpected, f"text columns not listed as prose (JSON blocks belong in columns or child tables): {unexpected}"


def test_no_control_table_keeps_a_json_or_serialized_blob(control_config):
    found = _text_columns(control_config)
    unexpected = sorted(found - PROSE_TEXT_COLUMNS)
    assert not unexpected, f"text columns not listed as prose: {unexpected}"


def test_the_prose_list_has_no_stale_entries(mysql_config, control_config):
    store.get_connection(mysql_config).close()
    assert PROSE_TEXT_COLUMNS == _text_columns(mysql_config) | _text_columns(control_config)


def test_every_searched_column_has_an_index_that_starts_with_it(mysql_config):
    store.get_connection(mysql_config).close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        first = {(row["t"], row["c"]) for row in conn.execute(
            "SELECT table_name AS t, column_name AS c FROM information_schema.statistics"
            " WHERE table_schema = DATABASE() AND seq_in_index = 1").fetchall()}
    finally:
        conn.close()
    missing = [(table, column) for table, columns in SEARCHED_COLUMNS.items() for column in columns
               if (table, column) not in first]
    assert not missing, f"searched columns without an index: {missing}"


def test_a_generation_runs_command_line_is_one_row_per_argument(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        run_id = store.start_generation_run(conn, "sector", ["--ring", "3", "--name", "a b"])
        rows = conn.execute("SELECT position, value FROM generation_run_arguments WHERE run_id = ?"
                            " ORDER BY position", (run_id,)).fetchall()
        assert [(row["position"], row["value"]) for row in rows] == [(0, "--ring"), (1, "3"), (2, "--name"),
                                                                     (3, "a b")]
        conn.execute("DELETE FROM generation_runs WHERE id = ?", (run_id,))
        conn.commit()
        assert conn.execute("SELECT COUNT(*) AS n FROM generation_run_arguments").fetchone()["n"] == 0
    finally:
        conn.close()


def test_the_v63_migration_moves_the_json_command_lines_into_rows(mysql_config):
    load_old_schema(mysql_config, 61)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute(
            "INSERT INTO generation_runs (command, arguments, version_key, planetgen_version, python_version,"
            " platform) VALUES ('sector', ?, 'k', 'v', 'p', 'x')", (json.dumps(["--ring", "3", "x" * 2000]),))
        conn.execute(
            "INSERT INTO generation_runs (command, arguments, version_key, planetgen_version, python_version,"
            " platform) VALUES ('plan', 'not json', 'k', 'v', 'p', 'x')")
        conn.commit()
    finally:
        conn.close()

    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION

    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = conn.execute("SELECT r.command, a.position, a.value FROM generation_runs r"
                            " JOIN generation_run_arguments a ON a.run_id = r.id ORDER BY a.position").fetchall()
        assert [(row["command"], row["position"], row["value"][:3]) for row in rows] == [
            ("sector", 0, "--r"), ("sector", 1, "3"), ("sector", 2, "xxx")]
        assert len(rows[2]["value"]) == store.GENERATION_ARGUMENT_LENGTH
        columns = {row["c"] for row in conn.execute(
            "SELECT column_name AS c FROM information_schema.columns WHERE table_schema = DATABASE()"
            " AND table_name = 'generation_runs'").fetchall()}
        assert "arguments" not in columns
    finally:
        conn.close()


def test_control_schema_moves_job_command_lines_into_rows_and_drops_unread_results(control_config):
    conn = store.get_control_connection(control_config)
    try:
        conn.execute("ALTER TABLE work_jobs ADD COLUMN argv TEXT NULL")
        conn.execute("ALTER TABLE work_tasks ADD COLUMN result TEXT NULL")
        conn.execute(
            "INSERT INTO work_jobs (id, title, holder, state, workers, created_at, heartbeat_at, argv)"
            " VALUES ('j1', 't', 'h', 'done', 0, NOW(6), NOW(6), ?)", (json.dumps(["galaxy", "--fast"]),))
        conn.execute(
            "INSERT INTO work_jobs (id, title, holder, state, workers, created_at, heartbeat_at, argv)"
            " VALUES ('j2', 't', 'h', 'done', 0, NOW(6), NOW(6), '{broken')")
        conn.commit()
    finally:
        conn.close()
    store.close_pool(control_config)

    conn = store.get_control_connection(control_config, ensure_schema=True)
    try:
        rows = conn.execute("SELECT job_id, position, value FROM work_job_args ORDER BY job_id, position").fetchall()
        assert [(row["job_id"], row["position"], row["value"]) for row in rows] == [
            ("j1", 0, "galaxy"), ("j1", 1, "--fast")]
        columns = {(row["t"], row["c"]) for row in conn.execute(
            "SELECT table_name AS t, column_name AS c FROM information_schema.columns"
            " WHERE table_schema = DATABASE()").fetchall()}
        assert ("work_jobs", "argv") not in columns and ("work_tasks", "result") not in columns
    finally:
        conn.close()
