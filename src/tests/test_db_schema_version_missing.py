# tests/test_db_schema_version_missing.py

"""
DB.97: a database that has the galaxy tables but no recorded schema version
(`schema_migrations` emptied or dropped) is refused with a clear message, not
guessed at from the shape of its tables. A new, empty database is created at
the current version, and a recorded version older than the Alembic baseline
(v61) is still refused.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

# A fresh database with no schema, which these tests build themselves (TEST.139 stage 3).
pytestmark = pytest.mark.fresh_database


def _forget_version(config, drop=False):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        conn.execute("DROP TABLE schema_migrations" if drop else "DELETE FROM schema_migrations")
        conn.commit()
    finally:
        conn.close()


def _recorded_rows(config):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        found = conn.execute("SELECT COUNT(*) AS n FROM information_schema.tables"
                             " WHERE table_schema = DATABASE() AND table_name = 'schema_migrations'").fetchone()["n"]
        return conn.execute("SELECT COUNT(*) AS n FROM schema_migrations").fetchone()["n"] if found else None
    finally:
        conn.close()


def test_a_new_database_is_created_at_the_current_version(mysql_config):
    assert store.schema_status(mysql_config) == (store.SCHEMA_VERSION, 0)
    store.get_connection(mysql_config).close()
    assert store.schema_status(mysql_config) == (store.SCHEMA_VERSION, 0)


@pytest.mark.parametrize("drop", [False, True], ids=["emptied", "dropped"])
def test_a_database_without_its_version_is_refused_untouched(mysql_config, drop):
    store.get_connection(mysql_config).close()
    _forget_version(mysql_config, drop=drop)
    before = _recorded_rows(mysql_config)
    store._schema_ensured.discard(mysql_config._key())

    with pytest.raises(store.SchemaVersionMissingError, match="no schema version recorded"):
        store.get_connection(mysql_config)
    with pytest.raises(store.SchemaVersionMissingError):
        store.schema_status(mysql_config)
    with pytest.raises(store.SchemaVersionMissingError):
        store.migrate_database(mysql_config)
    assert _recorded_rows(mysql_config) == before


def test_a_recorded_old_version_is_refused(mysql_config):
    load_old_schema(mysql_config, 61)
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
