# tests/test_dbpool.py

"""The worker's reused database (src/tests/dbpool.py, TEST.139 stage 3): a reset leaves it the way a new one is,
and a test that changes the schema does not hand its changes to the next one."""

import pymysql
import pytest

from planetgen.db import store
from tests import dbpool
from tests.conftest import _test_server_kwargs
from tests.test_api import _place_sector


@pytest.fixture
def worker(_mysql_server_available):
    database = dbpool.WorkerDatabase(_test_server_kwargs(), deep_check=True)
    yield database
    database.close()


def _count(config, table):
    conn = store.get_connection(config)
    try:
        return conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]
    finally:
        conn.close()


def test_a_reset_database_is_the_same_as_a_new_one(worker):
    """The deep check compares every table's definition, row and auto-increment counter with a new database."""
    config = worker.lease()
    _place_sector(config, "One", (5.0, 5.0, 5.0))
    _place_sector(config, "Two", (50.0, 5.0, 5.0))
    assert _count(config, "star_systems") == 2
    assert worker.release() == ""  # raises AssertionError when anything differs from a new database
    config = worker.lease()
    assert _count(config, "star_systems") == 0 and _count(config, "sectors") == 0
    assert _count(config, "schema_migrations") == 1
    _place_sector(config, "Again", (5.0, 5.0, 5.0))
    conn = store.get_connection(config)
    try:
        assert conn.execute("SELECT MIN(id) AS first FROM star_systems").fetchone()["first"] == 1
    finally:
        conn.close()
    assert worker.release() == ""
    assert worker.resets == 2 and worker.dropped == 0


def test_an_advanced_counter_is_set_back_even_when_the_rows_are_gone(worker):
    config = worker.lease()
    _place_sector(config, "One", (5.0, 5.0, 5.0))
    conn = store.get_connection(config)
    try:
        conn.execute("DELETE FROM star_systems")  # the rows go, the counter stays
        conn.commit()
    finally:
        conn.close()
    assert worker.release() == ""


def test_a_test_that_changes_the_schema_does_not_leave_it_to_the_next(worker):
    config = worker.lease()
    name = config.database
    admin = pymysql.connect(**_test_server_kwargs(), database=name)
    try:
        with admin.cursor() as cur:
            cur.execute("ALTER TABLE schema_migrations ADD COLUMN leftover INT")
        admin.commit()
    finally:
        admin.close()
    assert worker.release() == "the schema changed"
    assert worker.dropped == 1
    again = worker.lease()
    assert again.database != name
    conn = store.get_connection(again)
    try:
        columns = [row["Field"] for row in conn.execute("SHOW COLUMNS FROM schema_migrations").fetchall()]
    finally:
        conn.close()
    assert "leftover" not in columns


def test_an_added_table_is_noticed(worker):
    config = worker.lease()
    admin = pymysql.connect(**_test_server_kwargs(), database=config.database)
    try:
        with admin.cursor() as cur:
            cur.execute("CREATE TABLE stray (id INT)")
        admin.commit()
    finally:
        admin.close()
    assert worker.release() == "the schema changed"
