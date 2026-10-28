# tests/test_db_schema_current.py

"""
PERF.174: a connection to a database that already holds the current schema version does not apply `schema.sql` again
(its `CREATE TABLE IF NOT EXISTS` statements and views ran once per work unit, a new process each); a new database and
one behind the code still get it.
"""

import pytest

from planetgen.db import store

# A fresh database with no schema, which these tests build themselves (TEST.139 stage 3).
pytestmark = pytest.mark.fresh_database


@pytest.fixture
def applied(monkeypatch):
    """The databases `store._apply_schema` was run against, in order."""
    calls = []
    real = store._apply_schema

    def counting(conn):
        calls.append(conn.execute("SELECT DATABASE() AS db").fetchone()["db"])
        return real(conn)

    monkeypatch.setattr(store, "_apply_schema", counting)
    return calls


def test_a_new_database_gets_the_schema_once_and_a_current_one_is_left_alone(mysql_config, applied):
    store._schema_ensured.discard(mysql_config._key())
    store.get_connection(mysql_config).close()
    assert applied == [mysql_config.database]
    assert store.schema_status(mysql_config) == (store.SCHEMA_VERSION, 0)

    for _ in range(3):                                          # a new process each: nothing remembered
        store._schema_ensured.discard(mysql_config._key())
        store.get_connection(mysql_config).close()
    assert applied == [mysql_config.database]


def test_only_the_code_s_own_version_counts_as_current(mysql_config):
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not store._schema_is_current(conn)               # nothing built yet
        store._ensure_schema(conn)
        assert store._schema_is_current(conn)
        conn.execute("UPDATE schema_migrations SET version = ?", (store.SCHEMA_VERSION - 1,))
        conn.commit()
        assert not store._schema_is_current(conn)
    finally:
        conn.close()
