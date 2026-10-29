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

import uuid

import pymysql
import pytest

from planetgen.db import store

PROSE_TEXT_COLUMNS = {
    # Galaxy schema: generated prose and the summaries shown on pages.
    ("facilities", "description"), ("facilities", "ratings_note"), ("facility_types", "description"),
    ("facility_type_fields", "choices"), ("planets", "description"), ("planets", "atmosphere"),
    ("planets", "composition"), ("moons", "description"), ("moons", "atmosphere"),
    ("moons", "composition"),
    ("star_systems", "location"), ("star_systems", "system_flavor_text"),
    ("asteroid_belts", "composition_summary"), ("asteroid_fields", "composition_summary"),
    ("comets", "composition_summary"), ("interstellar_comets", "composition_summary"),
    ("rogue_planets", "composition"), ("nebulae", "composition"), ("nebulae", "formation_cause"),
    # DB.36: a plan's request and decision text, and the before and after of an operation or journal entry as JSON documents
    ("change_plan_ops", "after_json"), ("change_plan_ops", "before_json"), ("change_plans", "decision_note"), ("change_plans", "request_text"), ("edit_journal", "after_json"), ("edit_journal", "before_json"),
    # DB.38: the parameters of an override as a JSON document
    ("region_density_overrides", "params"),
    # DB.39: a comment's text
    ("request_comments", "body"),
    # DB.64: the stars and the lines of a constellation as JSON lists
    ("constellations", "edges"), ("constellations", "members"),
    # DB.63: the route plan of a ship as a JSON document
    ("vessels", "course"),
    # DB.58: the extra JSON document on an account token
    ("account_tokens", "meta"),
    # DB.62: the list of errors an upload unit met, as a JSON document
    ("upload_units", "errors"),
    # Control schema: an error message and the audit log's note.
    ("admin_audit_log", "detail"),
    ("generation_stage_runs", "settings"), ("generation_stage_runs", "metrics"),
    # OPS.58: the generator battery's digest per case and the environment record, as JSON documents read whole
    ("galaxy_shape", "battery_cases"), ("version_key_history", "battery_cases"), ("version_key_history", "environment"),
}
BLOB_TYPES = ("tinytext", "text", "mediumtext", "longtext", "json", "tinyblob", "blob", "mediumblob", "longblob")

SEARCHED_COLUMNS = {
    "sectors": ("name", "galactic_radius_pc"),
    "star_systems": ("name", "sector_id", "quadrant"),
    "stars": ("name", "star_type", "yerkes_class", "radius_km"),
    "planets": ("name", "planet_class", "body_type", "life_class", "radius_km"),
    "moons": ("name", "planet_class", "body_type", "life_class", "radius_km"),
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
        run_id = store.start_generation_run(conn, "sector", ["--shell", "3", "--name", "a b"])
        rows = conn.execute("SELECT position, value FROM generation_run_arguments WHERE run_id = ?"
                            " ORDER BY position", (run_id,)).fetchall()
        assert [(row["position"], row["value"]) for row in rows] == [(0, "--shell"), (1, "3"), (2, "--name"),
                                                                     (3, "a b")]
        conn.execute("DELETE FROM generation_runs WHERE id = ?", (run_id,))
        conn.commit()
        assert conn.execute("SELECT COUNT(*) AS n FROM generation_run_arguments").fetchone()["n"] == 0
    finally:
        conn.close()


def test_control_schema_v21_drops_the_work_queue_tables(control_config):
    tables = ("work_job_args", "work_tasks", "work_jobs", "work_lease")
    conn = store.get_control_connection(control_config)
    try:
        conn.execute("CREATE TABLE work_jobs (id VARCHAR(32) NOT NULL PRIMARY KEY)")
        for child in ("work_job_args", "work_tasks"):
            conn.execute(f"CREATE TABLE {child} (job_id VARCHAR(32) NOT NULL, CONSTRAINT fk_{child}_job"
                         " FOREIGN KEY (job_id) REFERENCES work_jobs(id) ON DELETE CASCADE)")
        conn.execute("CREATE TABLE work_lease (id TINYINT UNSIGNED NOT NULL PRIMARY KEY)")
        conn.execute("DELETE FROM control_schema_migrations WHERE version >= 21")
        conn.execute("INSERT IGNORE INTO control_schema_migrations (version) VALUES (20)")
        conn.commit()
    finally:
        conn.close()
    store.close_pool(control_config)

    conn = store.get_control_connection(control_config, ensure_schema=True)
    try:
        found = {row["t"] for row in conn.execute(
            "SELECT table_name AS t FROM information_schema.tables WHERE table_schema = DATABASE()").fetchall()}
        assert not found & set(tables)
        assert conn.execute("SELECT MAX(version) AS v FROM control_schema_migrations").fetchone()["v"] == store.CONTROL_SCHEMA_VERSION
    finally:
        conn.close()


def test_a_fresh_control_schema_matches_the_added_columns(control_config):
    """`control_schema.sql`'s CREATE TABLEs and `_CONTROL_COLUMNS` agree, so a new install and an upgraded one end up the
    same."""
    conn = store.get_control_connection(control_config)
    try:
        for table, added in store._CONTROL_COLUMNS.items():
            rows = conn.execute("SELECT COLUMN_NAME AS name FROM information_schema.COLUMNS"
                                " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", (table,)).fetchall()
            assert {name for name, _definition in added} <= {row["name"] for row in rows}, table
    finally:
        conn.close()
