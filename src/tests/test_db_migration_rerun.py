# tests/test_db_migration_rerun.py

"""
TEST.9: a migration that stops partway is finished by running it again.
MySQL commits each ALTER/CREATE/DROP as it runs, so a step that fails
after its first change leaves that change behind; the next
`migrate_database` must still end at a database exactly like a new one.
Also: every step applied twice in a row, and an empty or gapped
`schema_migrations`.
"""

import pytest

from stellarObjects import _db
from tests.db_schema_support import (
    load_old_schema,
    old_schema_versions,
    schema_differences,
    schema_snapshot,
    scratch_database,
)

_fresh_snapshots = {}

_CHANGES = ("ALTER", "CREATE", "DROP", "UPDATE", "DELETE", "INSERT")


class InjectedCrash(Exception):
    pass


def fresh_snapshot(beside):
    server = (beside.host, beside.port)
    if server not in _fresh_snapshots:
        with scratch_database(beside) as fresh:
            _db.get_connection(fresh).close()
            _fresh_snapshots[server] = schema_snapshot(fresh)
    return _fresh_snapshots[server]


def assert_like_new(config):
    differences = schema_differences(schema_snapshot(config), fresh_snapshot(config))
    assert not differences, "\n".join(differences)
    assert _db.schema_status(config) == (_db.SCHEMA_VERSION, 0)


def _crashing(step, when):
    """`step`, stopped with `InjectedCrash` either right after its first
    change ran, or just before it records its version."""
    def run(conn):
        real = conn.execute

        def execute(sql, params=()):
            statement = sql.lstrip().upper()
            if when == "before_version" and statement.startswith("INSERT INTO SCHEMA_MIGRATIONS"):
                raise InjectedCrash(sql)
            result = real(sql, params)
            if when == "after_first_change" and statement.startswith(_CHANGES):
                raise InjectedCrash(sql)
            return result

        conn.execute = execute
        step(conn)
    return run


def _start_for(target):
    """The newest checked-in schema older than `target`."""
    return max(version for version in old_schema_versions() if version < target)


STEP_TARGETS = [target for target, _step in _db._migration_steps()]


@pytest.mark.parametrize("when", ["after_first_change", "before_version"])
@pytest.mark.parametrize("target", STEP_TARGETS)
def test_a_step_that_crashed_is_finished_by_the_next_run(mysql_config, monkeypatch, target, when):
    load_old_schema(mysql_config, _start_for(target))
    real_steps = _db._migration_steps()
    crashing = [(version, _crashing(step, when) if version == target else step) for version, step in real_steps]
    monkeypatch.setattr(_db, "_migration_steps", lambda: crashing)
    with pytest.raises(InjectedCrash):
        _db.migrate_database(mysql_config)
    assert _db.schema_status(mysql_config)[0] < target

    monkeypatch.setattr(_db, "_migration_steps", lambda: real_steps)
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION
    assert_like_new(mysql_config)


def test_every_step_applied_twice(mysql_config, monkeypatch):
    load_old_schema(mysql_config, 8)
    real_steps = _db._migration_steps()

    def twice(step):
        def run(conn):
            step(conn)
            step(conn)
        return run

    monkeypatch.setattr(_db, "_migration_steps", lambda: [(version, twice(step)) for version, step in real_steps])
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION
    assert_like_new(mysql_config)


def test_rerunning_a_finished_migration_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 20)
    _db.migrate_database(mysql_config)
    before = schema_snapshot(mysql_config)
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION
    assert schema_snapshot(mysql_config) == before


def test_empty_schema_migrations_on_a_current_database_is_taken_as_current(mysql_config):
    _db.get_connection(mysql_config).close()
    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM schema_migrations")
        conn.commit()
    finally:
        conn.close()
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION
    assert_like_new(mysql_config)


def test_gapped_schema_migrations_counts_its_highest_version(mysql_config):
    """Rows 8 and 20 only (9-19 lost): the database is v20, and the steps
    from 21 on run."""
    load_old_schema(mysql_config, 20)
    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (8), (20)")
        conn.commit()
    finally:
        conn.close()
    assert _db.schema_status(mysql_config) == (20, _db.SCHEMA_VERSION - 20)
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION
    assert_like_new(mysql_config)
