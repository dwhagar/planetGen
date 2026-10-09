# tests/test_update_database_step.py

"""
TEST.62: the database step update.sh and install.sh share
(`migrate_or_reset_db` and `offer_population_pass` in
scripts/deploy-common.sh), run for real in bash against a live database,
the way a scheduled update runs it (no terminal to ask on):

- a new, empty database is brought to the current schema;
- a database that's up to date is left alone, and nothing is asked;
- one that needs migrating is migrated, keeping its data;
- one newer than the code is refused, untouched, and one older than the
  Alembic baseline too;
- an unreachable server stops the step with an error;
- a migration that fails (an account without CREATE/ALTER) stops the
  step with an error and leaves the version where it was.

A stopped step stops update.sh and install.sh (`set -e`); the CI job
`linux-update` runs the whole update.sh on top of this.
"""

import os
import subprocess
import sys
import uuid

import pymysql
import pytest

from planetgen.db import store
from tests.bughunt_support import mysql_argv, run_cli
from tests.conftest import _test_server_kwargs
from tests.db_schema_support import migrations_with_probes

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

pytestmark = pytest.mark.skipif(sys.platform.startswith("win"), reason="update.sh is Linux and macOS only")


@pytest.fixture
def control_db(monkeypatch):
    """A control schema of the test's own, dropped afterwards."""
    name = f"planetgen_test_ctl_{uuid.uuid4().hex[:12]}"
    yield name
    _admin(f"DROP DATABASE IF EXISTS `{name}`")


def _admin(*statements):
    conn = pymysql.connect(**_test_server_kwargs())
    try:
        with conn.cursor() as cur:
            for statement in statements:
                cur.execute(statement)
        conn.commit()
    finally:
        conn.close()


def _run_step(config, control_db, **overrides):
    """Runs the shared database step in bash with stdin not a terminal,
    pointed at `config` through the same PLANETGEN_MYSQL_* variables
    update.sh's environment would carry."""
    env = {k: v for k, v in os.environ.items() if not k.startswith("PLANETGEN_")}
    env.update(
        PLANETGEN_MYSQL_HOST=config.host, PLANETGEN_MYSQL_PORT=str(config.port),
        PLANETGEN_MYSQL_USER=config.user, PLANETGEN_MYSQL_PASSWORD=config.password,
        PLANETGEN_MYSQL_DATABASE=config.database, PLANETGEN_CONTROL_DATABASE=control_db,
    )
    env.update(overrides)
    script = ('set -euo pipefail; SCRIPT_DIR="$1"; PYTHON="$2"; '
              'source "$SCRIPT_DIR/scripts/deploy-common.sh"; migrate_or_reset_db; offer_population_pass')
    return subprocess.run(["bash", "-c", script, "update-step", REPO_DIR, sys.executable], env=env,
                          stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=300)


def _version(config):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        return conn.execute("SELECT MAX(version) AS v FROM schema_migrations").fetchone()["v"]
    finally:
        conn.close()


def _systems(config):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        return conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"]
    finally:
        conn.close()


def _extra_revision(tmp_path):
    """`PLANETGEN_MIGRATIONS_DIR` for the step's own Python: the real
    revisions plus a test revision 0062, so a database at v61 has one
    migration step to run."""
    return {"PLANETGEN_MIGRATIONS_DIR": migrations_with_probes(tmp_path, up_to=62)}


def _make_older_than_the_baseline(config):
    store.get_connection(config).close()
    conn = store.get_connection(config, ensure_schema=False)
    try:
        conn.execute("ALTER TABLE star_systems DROP COLUMN velocity_x_kms")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (55)")
        conn.commit()
    finally:
        conn.close()
    store.close_pool(config)


def test_a_new_empty_database_is_brought_current(mysql_config, control_db):
    result = _run_step(mysql_config, control_db)
    assert result.returncode == 0, result.stderr
    assert f"Database is at schema v{store.SCHEMA_VERSION} (current)." in result.stdout
    assert "Control schema (admin logins) is up to date." in result.stdout
    assert "Skipping the population pass" in result.stdout
    assert _version(mysql_config) == store.SCHEMA_VERSION


def test_an_up_to_date_database_is_left_alone_and_nothing_is_asked(mysql_config, control_db):
    assert _run_step(mysql_config, control_db).returncode == 0
    result = _run_step(mysql_config, control_db)
    assert result.returncode == 0, result.stderr
    assert "migration step(s)" not in result.stdout
    assert "Delete all galaxy data" not in result.stdout
    assert f"schema v{store.SCHEMA_VERSION} (current)" in result.stdout
    assert "password:" not in result.stdout  # the first admin password only shows once


def test_a_database_that_needs_migrating_is_migrated_and_keeps_its_data(mysql_config, control_db, tmp_path):
    run_cli("system", mysql_argv(mysql_config))
    systems = _systems(mysql_config)
    assert systems > 0
    result = _run_step(mysql_config, control_db, **_extra_revision(tmp_path))
    assert result.returncode == 0, result.stderr
    assert (f"is at schema v{store.SCHEMA_VERSION}; this version needs v{store.SCHEMA_VERSION + 1} "
            f"(1 migration step(s))") in result.stdout
    assert "No terminal to ask on: keeping the data and migrating it." in result.stdout
    assert _version(mysql_config) == store.SCHEMA_VERSION + 1
    assert _systems(mysql_config) == systems


def test_a_database_older_than_the_baseline_is_refused_untouched(mysql_config, control_db):
    _make_older_than_the_baseline(mysql_config)
    result = _run_step(mysql_config, control_db)
    assert result.returncode != 0
    assert "older than v61" in result.stderr
    assert "Traceback" not in result.stderr
    assert _version(mysql_config) == 55


def test_a_database_newer_than_the_code_is_refused_untouched(mysql_config, control_db):
    store.get_connection(mysql_config).close()
    newer = store.SCHEMA_VERSION + 1
    _admin(f"INSERT INTO `{mysql_config.database}`.schema_migrations (version) VALUES ({newer})")
    result = _run_step(mysql_config, control_db)
    assert result.returncode != 0
    assert "newer than this code" in result.stderr
    assert "Traceback" not in result.stderr
    assert _version(mysql_config) == newer


def test_an_unreachable_server_stops_the_step(mysql_config, control_db):
    result = _run_step(mysql_config, control_db, PLANETGEN_MYSQL_PORT="1")
    assert result.returncode != 0
    assert result.stderr.startswith("error:")
    assert "Traceback" not in result.stderr


def test_a_failed_migration_stops_the_step_and_keeps_the_version(mysql_config, control_db, tmp_path):
    store.get_connection(mysql_config).close()
    user = f"pgro_{uuid.uuid4().hex[:10]}"
    _admin(f"CREATE USER '{user}'@'%' IDENTIFIED BY 'read-only-pw'",
           f"GRANT SELECT ON `{mysql_config.database}`.* TO '{user}'@'%'")
    try:
        result = _run_step(mysql_config, control_db, PLANETGEN_MYSQL_USER=user,
                           PLANETGEN_MYSQL_PASSWORD="read-only-pw", **_extra_revision(tmp_path))
    finally:
        _admin(f"DROP USER IF EXISTS '{user}'@'%'")
    assert f"is at schema v{store.SCHEMA_VERSION}" in result.stdout
    assert result.returncode != 0
    assert "error:" in result.stderr and "denied" in result.stderr
    assert _version(mysql_config) == store.SCHEMA_VERSION
