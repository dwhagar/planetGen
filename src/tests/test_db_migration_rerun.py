# tests/test_db_migration_rerun.py

"""
TEST.9: a finished migration run again changes nothing, and a missing or
emptied `schema_migrations` on a current database is taken as current.
"""

from planetgen.db import alembic_runner, store
from tests.db_schema_support import (
    load_old_schema,
    schema_differences,
    schema_snapshot,
    scratch_database,
)

_fresh_snapshots = {}


def fresh_snapshot(beside):
    server = (beside.host, beside.port)
    if server not in _fresh_snapshots:
        with scratch_database(beside) as fresh:
            store.get_connection(fresh).close()
            _fresh_snapshots[server] = schema_snapshot(fresh)
    return _fresh_snapshots[server]


def assert_like_new(config):
    differences = schema_differences(schema_snapshot(config), fresh_snapshot(config))
    assert not differences, "\n".join(differences)
    assert store.schema_status(config) == (store.SCHEMA_VERSION, 0)


def test_rerunning_a_finished_migration_changes_nothing(mysql_config):
    load_old_schema(mysql_config, alembic_runner.BASELINE_VERSION)
    store.migrate_database(mysql_config)
    before = schema_snapshot(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert schema_snapshot(mysql_config) == before


def test_gapped_schema_migrations_counts_its_highest_version(mysql_config):
    """Rows 8 and the baseline only: the database is the baseline's version."""
    load_old_schema(mysql_config, alembic_runner.BASELINE_VERSION)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (8), (?)", (alembic_runner.BASELINE_VERSION,))
        conn.commit()
    finally:
        conn.close()
    assert store.schema_status(mysql_config) == (
        alembic_runner.BASELINE_VERSION, store.SCHEMA_VERSION - alembic_runner.BASELINE_VERSION)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert_like_new(mysql_config)
