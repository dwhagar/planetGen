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
