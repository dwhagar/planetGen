# tests/test_drop_database_holders.py

"""
A test database is dropped without waiting for connections that still hold it (TEST.193).

In the gate, the module fixture of the browser tests (`site_db`) ran `DROP DATABASE` and waited 853 s: another
connection (state Sleep, an open transaction on one of the database's tables) held a metadata lock. The drop
helper of `tests/dbpool.py` ends such connections first and never waits more than a few seconds.
"""

import time
import uuid

import pymysql
import pytest

from tests import dbpool
from tests.conftest import _test_server_kwargs

pytestmark = pytest.mark.db


@pytest.fixture
def server(_mysql_server_available):
    kwargs = _test_server_kwargs()
    name = f"planetgen_test_{uuid.uuid4().hex[:16]}"
    admin = pymysql.connect(**kwargs)
    with admin.cursor() as cur:
        cur.execute(f"CREATE DATABASE `{name}`")
        cur.execute(f"CREATE TABLE `{name}`.t (id INT PRIMARY KEY)")
        cur.execute(f"INSERT INTO `{name}`.t VALUES (1)")
    admin.commit()
    opened = []
    try:
        yield kwargs, name, opened, admin
    finally:
        for conn in opened:
            try:
                conn.close()
            except pymysql.MySQLError:
                pass
        with admin.cursor() as cur:
            cur.execute("SET SESSION lock_wait_timeout = 10")
            cur.execute(f"DROP DATABASE IF EXISTS `{name}`")
        admin.close()


def _exists(admin, name):
    with admin.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM information_schema.SCHEMATA WHERE SCHEMA_NAME = %s", (name,))
        return cur.fetchone()[0] == 1


def test_a_connection_with_an_open_transaction_on_a_table_no_longer_stalls_the_drop(server):
    kwargs, name, opened, admin = server
    holder = pymysql.connect(**kwargs, database=name)  # a pooled connection: its current database is the test's
    opened.append(holder)
    with holder.cursor() as cur:
        cur.execute("SELECT * FROM t")  # a transaction that was never ended: holds the table's metadata lock
        holder_id = holder.thread_id()
    started = time.monotonic()
    with pytest.warns(UserWarning, match="TEST.193"):
        ended = dbpool.drop_database(kwargs, name, wait_seconds=5)
    assert time.monotonic() - started < 5
    assert [row[0] for row in ended] == [holder_id]
    assert not _exists(admin, name)


def test_an_idle_connection_whose_current_database_it_is_is_ended_too(server):
    kwargs, name, opened, admin = server
    holder = pymysql.connect(**kwargs, database=name)
    opened.append(holder)
    holder_id = holder.thread_id()
    with pytest.warns(UserWarning, match="TEST.193"):
        ended = dbpool.drop_database(kwargs, name)
    assert [row[0] for row in ended] == [holder_id]
    assert not _exists(admin, name)


def test_a_connection_that_has_nothing_to_do_with_the_database_is_left_alone(server):
    kwargs, name, opened, admin = server
    bystander = pymysql.connect(**kwargs)
    opened.append(bystander)
    dbpool.drop_database(kwargs, name)
    with bystander.cursor() as cur:
        cur.execute("SELECT 1")
        assert cur.fetchone()[0] == 1


def test_dropping_a_database_nobody_holds_ends_nothing_and_says_nothing(server, recwarn):
    kwargs, name, _opened, admin = server
    assert dbpool.drop_database(kwargs, name) == []
    assert not [w for w in recwarn if "TEST.193" in str(w.message)]
    assert not _exists(admin, name)


def test_a_pooled_connection_that_was_never_returned_is_named_by_the_code_that_took_it(server):
    from planetgen.db import store

    kwargs, name, _opened, admin = server
    config = store.MySQLConfig(database=name, **kwargs)
    leaked = store.get_connection(config, ensure_schema=False)   # never closed: the shape of the leak in the gate
    try:
        leaked.execute("SELECT * FROM t")
        with pytest.warns(UserWarning, match="taken in thread") as caught:
            ended = dbpool.drop_database(kwargs, name)
        assert len(ended) == 1
        assert "test_a_pooled_connection_that_was_never_returned" in str(caught[0].message)
    finally:
        store.close_pool(config)


def test_a_connection_that_reads_a_table_for_a_moment_is_not_ended(server):
    """TEST.262: `list_databases` of another worker counts the rows of every test database; it holds the table's
    metadata lock until it commits, and ending it cut that worker's query off ("Lost connection to MySQL server")."""
    import threading

    kwargs, name, opened, admin = server
    reader = pymysql.connect(**kwargs)   # no current database: it reads the table by its qualified name
    opened.append(reader)
    with reader.cursor() as cur:
        cur.execute(f"SELECT COUNT(*) FROM `{name}`.t")   # a transaction: holds the table's lock until it ends
    timer = threading.Timer(0.5, reader.commit)
    timer.start()
    try:
        ended = dbpool.drop_database(kwargs, name)
    finally:
        timer.join()
    assert ended == []
    with reader.cursor() as cur:
        cur.execute("SELECT 1")
        assert cur.fetchone()[0] == 1
    assert not _exists(admin, name)


def test_a_connection_that_keeps_reading_a_table_is_ended_after_the_gentle_wait(server):
    kwargs, name, opened, admin = server
    reader = pymysql.connect(**kwargs)
    opened.append(reader)
    with reader.cursor() as cur:
        cur.execute(f"SELECT COUNT(*) FROM `{name}`.t")   # never committed
        reader_id = reader.thread_id()
    started = time.monotonic()
    with pytest.warns(UserWarning, match="TEST.193"):
        ended = dbpool.drop_database(kwargs, name)
    assert [row[0] for row in ended] == [reader_id]
    assert dbpool.GENTLE_WAIT_SECONDS - 1 <= time.monotonic() - started < dbpool.DROP_WAIT_SECONDS + 3
    assert not _exists(admin, name)
