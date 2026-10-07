# tests/test_db_schema_version_detect.py

"""
DB.4: a database whose `schema_migrations` table was emptied (or lost)
is no longer taken for current. Its version is read from its tables
(`store.detect_schema_version`), and `migrate_database` runs the steps it is
missing. Checked against every released schema in `fixtures/old_schemas/`.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema, old_schema_versions

# Versions whose tables look the same as the one before (their step only
# changed rows): v30 reads as v29 (its step is safe to repeat), v35 as
# v35 even when it is v34 (repeating its step would delete sectors).
_READ_AS = {30: 29, 34: 35}


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


@pytest.mark.parametrize("version", old_schema_versions())
def test_every_released_schema_is_recognized(mysql_config, version):
    load_old_schema(mysql_config, version)
    _forget_version(mysql_config)
    assert _detected(mysql_config) == _READ_AS.get(version, version)


def test_a_current_database_is_recognized(mysql_config):
    store.get_connection(mysql_config).close()
    _forget_version(mysql_config)
    assert _detected(mysql_config) == store.SCHEMA_VERSION


def test_a_new_database_is_current(mysql_config):
    assert _detected(mysql_config) == store.SCHEMA_VERSION


@pytest.mark.parametrize("drop", [False, True], ids=["emptied", "dropped"])
def test_an_old_database_without_its_version_is_migrated(mysql_config, drop):
    load_old_schema(mysql_config, 44)
    _forget_version(mysql_config, drop=drop)

    assert store.schema_status(mysql_config) == (44, store.SCHEMA_VERSION - 44)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION

    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        versions = sorted(row["version"] for row in conn.execute("SELECT version FROM schema_migrations").fetchall())
        # v45's table and v50's columns, which only the migration adds.
        has_id_blocks = conn.execute(
            "SELECT COUNT(*) AS n FROM information_schema.tables"
            " WHERE table_schema = DATABASE() AND table_name = 'id_blocks'").fetchone()["n"]
    finally:
        conn.close()
    assert versions == list(range(44, store.SCHEMA_VERSION + 1))
    assert has_id_blocks == 1
    assert store._has_column(store.get_connection(mysql_config, ensure_schema=False), "system_configs", "comets")


def test_first_connection_records_the_detected_version(mysql_config):
    """A plain first connection (not a migration) records what the tables
    show rather than `SCHEMA_VERSION`, so update.sh still sees steps
    pending."""
    load_old_schema(mysql_config, 47)
    _forget_version(mysql_config)
    store._schema_ensured.discard(mysql_config._key())
    store.get_connection(mysql_config).close()
    assert store.schema_status(mysql_config) == (47, store.SCHEMA_VERSION - 47)
