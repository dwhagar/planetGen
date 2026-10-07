# tests/test_db_schema_race.py

"""
DB.5: several first connections reaching an empty database at once. Each
used to run `_ensure_schema`, and all but one failed inserting the same
baseline `schema_migrations` row (IntegrityError 1062). Now one creates
the schema under the database's schema lock and the others wait, then
find it current.
"""

import threading

import pytest

from planetgen.db import store

_CONNECTIONS = 6


def _run_together(count, work):
    """Runs `work(i)` in `count` threads released at the same moment;
    returns the exceptions raised."""
    barrier = threading.Barrier(count)
    errors = []

    def run(index):
        try:
            barrier.wait()
            work(index)
        except BaseException as exc:  # noqa: BLE001 - reported by the test
            errors.append(exc)

    threads = [threading.Thread(target=run, args=(i,)) for i in range(count)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=120)
    assert not any(thread.is_alive() for thread in threads), "a connection never finished"
    return errors


def _versions(config):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        return [row["version"] for row in conn.execute("SELECT version FROM schema_migrations").fetchall()]
    finally:
        conn.close()


@pytest.mark.parametrize("with_migration", [False, True], ids=["connections", "connections-and-migrate"])
def test_first_connections_at_once_create_the_schema_once(mysql_config, with_migration):
    store._schema_ensured.discard(mysql_config._key())

    def work(index):
        if with_migration and index == 0:
            store.migrate_database(mysql_config)
        else:
            conn = store.get_connection(mysql_config)
            try:
                conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()
            finally:
                conn.close()

    errors = _run_together(_CONNECTIONS, work)
    assert not errors, errors
    assert _versions(mysql_config) == [store.SCHEMA_VERSION]


def test_first_control_connections_at_once(mysql_config):
    control = store.MySQLConfig(host=mysql_config.host, port=mysql_config.port, user=mysql_config.user,
                              password=mysql_config.password, database=mysql_config.database)

    def work(_index):
        store.get_control_connection(control, ensure_schema=True).close()

    errors = _run_together(_CONNECTIONS, work)
    assert not errors, errors
    conn = store.get_connection(control, ensure_schema=False)
    try:
        rows = conn.execute("SELECT version FROM control_schema_migrations").fetchall()
    finally:
        conn.close()
    assert [row["version"] for row in rows] == [store.CONTROL_SCHEMA_VERSION]


def test_a_first_connection_waits_for_the_schema_lock(mysql_config):
    """Whether threads above collide depends on timing; this checks the
    lock itself: a first connection waits while another holds it."""
    holder = store.get_connection(mysql_config, ensure_schema=False)
    finished = threading.Event()
    try:
        with store._schema_lock(holder):
            store._schema_ensured.discard(mysql_config._key())
            thread = threading.Thread(target=lambda: (store.get_connection(mysql_config).close(), finished.set()))
            thread.start()
            assert not finished.wait(1.0), "the first connection did not wait for the lock"
        thread.join(timeout=60)
        assert finished.is_set()
    finally:
        holder.close()
    assert _versions(mysql_config) == [store.SCHEMA_VERSION]
