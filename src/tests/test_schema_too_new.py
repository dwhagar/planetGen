# tests/test_schema_too_new.py

"""
TEST.10: a galaxy or control database whose schema version is above this
code's (`_db.SCHEMA_VERSION`, `_db.CONTROL_SCHEMA_VERSION`) -- one a newer
planetGen has already migrated -- is refused with a clear message instead
of being migrated or written to:

- `migrateDb.py` (and so update.sh/install.sh) exits 1 naming both
  versions, for the galaxy database and for the control schema, and
  `--status` fails the same way rather than reporting "0 pending".
- A read-write connection (`get_connection`, what the generator uses)
  raises `SchemaTooNewError` before any of this code's older DDL runs.
- `/api/health` says the database is newer than the code, rather than
  telling the operator to run migrateDb.py.
"""

import sys
import uuid

import pymysql
import pytest

import migrateDb
from api.app import create_app
from api.config import Config
from stellarObjects import _db
from tests.bughunt_support import mysql_argv
from tests.conftest import _test_server_kwargs


def _run_main(module, argv, monkeypatch):
    monkeypatch.setattr(sys, "argv", [module.__name__ + ".py", *argv])
    try:
        module.main()
    except SystemExit as exc:
        return exc.code or 0
    return 0


def _bump_galaxy(config, version):
    """Lays the current schema down, then records `version` as if a newer
    planetGen had migrated it."""
    conn = _db.get_connection(config)
    try:
        conn.execute("INSERT INTO schema_migrations (version) VALUES (?)", (version,))
        conn.commit()
    finally:
        conn.close()
    _db.close_pool(config)  # forget that this process already ensured it


def _views(config):
    conn = _db.get_connection(config, ensure_schema=False)
    try:
        return conn.execute(
            "SELECT view_definition AS d FROM information_schema.views WHERE table_schema = ? "
            "AND table_name = 'sector_objects'", (config.database,)).fetchone()["d"]
    finally:
        conn.close()


def test_migrate_refuses_a_newer_galaxy_database(mysql_config, monkeypatch, capsys):
    newer = _db.SCHEMA_VERSION + 1
    _bump_galaxy(mysql_config, newer)
    assert _run_main(migrateDb, mysql_argv(mysql_config), monkeypatch) == 1
    err = capsys.readouterr().err
    assert err.startswith("error:")
    assert f"v{newer}" in err and f"v{_db.SCHEMA_VERSION}" in err
    assert "newer than this code" in err and "Update planetGen" in err
    assert "Traceback" not in err


def test_migrate_status_refuses_a_newer_galaxy_database(mysql_config, monkeypatch, capsys):
    _bump_galaxy(mysql_config, _db.SCHEMA_VERSION + 3)
    assert _run_main(migrateDb, mysql_argv(mysql_config) + ["--status"], monkeypatch) == 1
    captured = capsys.readouterr()
    assert captured.out == ""  # update.sh never reads "0 pending" for it
    assert "newer than this code" in captured.err


def test_a_write_connection_refuses_a_newer_database_before_any_ddl(mysql_config):
    _bump_galaxy(mysql_config, _db.SCHEMA_VERSION + 1)
    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        # Stand in for a newer schema's view; this code's schema.sql would
        # put its own CREATE OR REPLACE VIEW back.
        conn.execute("CREATE OR REPLACE VIEW sector_objects AS SELECT 1 AS newer_view")
        conn.commit()
    finally:
        conn.close()
    with pytest.raises(_db.SchemaTooNewError) as caught:
        _db.get_connection(mysql_config)
    assert caught.value.version == _db.SCHEMA_VERSION + 1
    assert caught.value.expected == _db.SCHEMA_VERSION
    assert mysql_config.database in str(caught.value)
    assert "newer_view" in _views(mysql_config)
    # Still refused on the next try: a refusal is never cached as "ensured".
    with pytest.raises(_db.SchemaTooNewError):
        _db.get_connection(mysql_config)


def test_migrate_database_and_schema_status_refuse_and_change_nothing(mysql_config):
    newer = _db.SCHEMA_VERSION + 2
    _bump_galaxy(mysql_config, newer)
    with pytest.raises(_db.SchemaTooNewError):
        _db.schema_status(mysql_config)
    with pytest.raises(_db.SchemaTooNewError):
        _db.migrate_database(mysql_config)
    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT MAX(version) AS v FROM schema_migrations").fetchone()["v"] == newer
    finally:
        conn.close()


def test_the_current_and_older_versions_are_still_accepted(mysql_config):
    _db.get_connection(mysql_config).close()
    assert _db.schema_status(mysql_config) == (_db.SCHEMA_VERSION, 0)
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION


@pytest.fixture
def control_schema(monkeypatch):
    """A control schema of this test's own (the session's is shared by
    every test in this process, so it must never be made "newer")."""
    name = f"planetgen_test_ctl_{uuid.uuid4().hex[:12]}"
    monkeypatch.setenv(_db.CONTROL_DB_ENV_VAR, name)
    yield name
    conn = pymysql.connect(**_test_server_kwargs())
    try:
        with conn.cursor() as cur:
            cur.execute(f"DROP DATABASE IF EXISTS `{name}`")
        conn.commit()
    finally:
        conn.close()


@pytest.fixture
def client(mysql_config):
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.config["TESTING"] = True
    return app.test_client()


def test_migrate_refuses_a_newer_control_schema(mysql_config, control_schema, monkeypatch, capsys):
    control = _db.control_mysql_config(mysql_config)
    assert control.database == control_schema
    assert _run_main(migrateDb, mysql_argv(mysql_config), monkeypatch) == 0
    capsys.readouterr()
    newer = _db.CONTROL_SCHEMA_VERSION + 1
    conn = _db.get_control_connection(control)
    try:
        conn.execute("INSERT INTO control_schema_migrations (version) VALUES (?)", (newer,))
        conn.commit()
    finally:
        conn.close()
    assert _run_main(migrateDb, mysql_argv(mysql_config), monkeypatch) == 1
    err = capsys.readouterr().err
    assert f"control schema v{newer}" in err and f"v{_db.CONTROL_SCHEMA_VERSION}" in err
    assert "newer than this code" in err
    with pytest.raises(_db.SchemaTooNewError):
        _db.get_control_connection(control, ensure_schema=True)


def test_health_says_the_database_is_newer_than_the_code(mysql_config, client):
    _bump_galaxy(mysql_config, _db.SCHEMA_VERSION + 1)
    body = client.get("/api/health").get_json()
    assert body["schema_current"] is False
    assert body["schema_version"] == _db.SCHEMA_VERSION + 1
    assert "newer than this code" in body["detail"]
    assert "run migrateDb.py" not in body["detail"]
