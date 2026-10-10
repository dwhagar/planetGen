# tests/test_db_schema_version_detect.py

"""
DB.4: a database whose `schema_migrations` table was emptied (or lost)
is no longer taken for current. Its version is read from its tables
(`store.detect_schema_version`), and `migrate_database` runs the revisions it
is missing. One older than the Alembic baseline (v61) is refused.
"""

import pytest

from planetgen.db import alembic_runner, store
from tests.db_schema_support import load_old_schema, old_schema_versions


def _forget_version(config, drop=False):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        conn.execute("DROP TABLE schema_migrations" if drop else "DELETE FROM schema_migrations")
        conn.commit()
    finally:
        conn.close()


def _detected(config):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        return store.detect_schema_version(conn)
    finally:
        conn.close()


def _make_older_than_the_baseline(config):
    """A current database with the baseline's marker column, and every later marker, gone."""
    store.get_connection(config).close()
    conn = store.get_connection(config, ensure_schema=False)
    try:
        conn.execute("ALTER TABLE planets DROP COLUMN equipment_tier")  # v76
        conn.execute("ALTER TABLE planets DROP COLUMN surface_dose_msv_yr")  # v75
        conn.execute("ALTER TABLE planets DROP COLUMN ocean_class")  # v74
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN bright_star_mass_limit_sol")  # v73
        conn.execute("ALTER TABLE stars DROP COLUMN l_xuv_w")  # v72
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN phenomenon_min_mass_solar")  # v71
        conn.execute("ALTER TABLE planets DROP COLUMN mantle_redox")  # v70
        conn.execute("ALTER TABLE stars DROP COLUMN axial_tilt_deg")  # v68
        conn.execute("ALTER TABLE sectors DROP COLUMN version_key")  # v67
        conn.execute("ALTER TABLE planets DROP COLUMN next_update_due")
        conn.execute("DROP TABLE phenomenon_scatter")
        conn.execute("DROP TABLE generation_run_arguments")
        conn.execute("DROP TABLE sector_path_knots")  # the later revisions' markers too
        conn.execute("DROP TABLE sector_paths")
        conn.execute("ALTER TABLE facilities DROP COLUMN velocity_x_kms")
        conn.execute("ALTER TABLE star_systems DROP COLUMN velocity_x_kms")
        conn.commit()
    finally:
        conn.close()


@pytest.mark.parametrize("version", old_schema_versions())
def test_every_checked_in_schema_is_recognized(mysql_config, version):
    load_old_schema(mysql_config, version)
    _forget_version(mysql_config)
    assert _detected(mysql_config) == version


def test_a_current_database_is_recognized(mysql_config):
    store.get_connection(mysql_config).close()
    _forget_version(mysql_config)
    assert _detected(mysql_config) == store.SCHEMA_VERSION


def test_a_new_database_is_current(mysql_config):
    assert _detected(mysql_config) == store.SCHEMA_VERSION


@pytest.mark.parametrize("drop", [False, True], ids=["emptied", "dropped"])
def test_a_baseline_database_without_its_version_is_migrated(mysql_config, drop):
    load_old_schema(mysql_config, alembic_runner.BASELINE_VERSION)
    _forget_version(mysql_config, drop=drop)

    assert store.schema_status(mysql_config) == (61, store.SCHEMA_VERSION - 61)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION

    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        versions = sorted(row["version"] for row in conn.execute("SELECT version FROM schema_migrations").fetchall())
    finally:
        conn.close()
    assert versions == list(range(61, store.SCHEMA_VERSION + 1))


def test_first_connection_records_the_detected_version(mysql_config):
    """A plain first connection (not a migration) records what the tables
    show rather than `SCHEMA_VERSION`."""
    load_old_schema(mysql_config, alembic_runner.BASELINE_VERSION)
    _forget_version(mysql_config)
    store._schema_ensured.discard(mysql_config._key())
    store.get_connection(mysql_config).close()
    assert store.schema_status(mysql_config) == (61, store.SCHEMA_VERSION - 61)


def test_a_database_older_than_the_baseline_is_refused_untouched(mysql_config):
    _make_older_than_the_baseline(mysql_config)
    _forget_version(mysql_config)
    assert _detected(mysql_config) == alembic_runner.BASELINE_VERSION - 1
    with pytest.raises(store.SchemaTooOldError, match="older than v61"):
        store.migrate_database(mysql_config)
    with pytest.raises(store.SchemaTooOldError):
        store.schema_status(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM schema_migrations").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_recorded_old_version_is_refused(mysql_config):
    store.get_connection(mysql_config).close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (55)")
        conn.commit()
    finally:
        conn.close()
    store._schema_ensured.discard(mysql_config._key())
    with pytest.raises(store.SchemaTooOldError, match="v55"):
        store.migrate_database(mysql_config)
