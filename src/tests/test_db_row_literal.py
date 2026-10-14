# tests/test_db_row_literal.py

"""
PERF.82: a batched INSERT's rows are written as SQL text by `store._row_literal`, not by the driver's per-value
formatting. The text has to be exactly what pymysql would have made, for every kind of value a row holds.
Against a real server's connection (`conftest.py`'s `mysql_config`); skipped when none is reachable.
"""

import random

import pymysql
import pytest

from planetgen.db import store

VALUES = [0, 1, -5, 12345678901234567890, 2.5, 0.1, 1e-9, 1e22, -0.0, float(3), None, True, False,
          "plain", "it's", 'say "hi"', "back\\slash", "line\nbreak", "nul\0byte", "cr\rx", "\x1a", "é☃", "", b"ab"]


@pytest.fixture
def driver(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        yield conn._conn
    finally:
        conn.close()


def _pymysql_text(driver, shape, row):
    return driver.cursor().mogrify(shape.replace("?", "%s"), row)


@pytest.mark.parametrize("shape", ["(?)", "(?, ?)", "(?, ?, ?, ?, ?)", "(?, 0, ?, NOW(), ?)"])
def test_a_row_is_written_as_pymysql_writes_it(driver, shape):
    marks = shape.count("?")
    pieces = None if store._ALL_MARKS.fullmatch(shape) else shape.split("?")
    rng = random.Random(82)
    for _ in range(200):
        row = tuple(rng.choice(VALUES) for _ in range(marks))
        assert store._row_literal(driver.escape, pieces, row) == _pymysql_text(driver, shape, row), row


def test_every_value_kind_is_covered_alone(driver):
    for value in VALUES:
        assert store._row_literal(driver.escape, None, (value,)) == _pymysql_text(driver, "(?)", (value,)), value


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_a_float_mysql_cannot_hold_is_refused_as_the_driver_refuses_it(driver, bad):
    with pytest.raises(pymysql.err.ProgrammingError):
        store._row_literal(driver.escape, None, (bad,))
    with pytest.raises(pymysql.err.ProgrammingError):
        _pymysql_text(driver, "(?)", (bad,))


def test_a_batched_insert_stores_the_values_it_was_given(mysql_config):
    """Names with quotes, backslashes and non-ASCII letters, and floats, come back as they went in."""
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("CREATE TABLE literal_check (id INT PRIMARY KEY, name VARCHAR(100), x DOUBLE, n BIGINT NULL)")
        rows = [(1, "O'Brien \\ \"q\"", 0.1, None), (2, "Zoë ☃", 1e22, 5), (3, "line\nbreak", -2.5e-9, -7)]
        head = "INSERT INTO literal_check (id, name, x, n) VALUES "
        cur = conn._conn.cursor()
        cur.execute(head + ",".join(store._row_literal(conn._conn.escape, None, row) for row in rows))
        conn.commit()
        stored = conn.execute("SELECT id, name, x, n FROM literal_check ORDER BY id").fetchall()
        assert [tuple(r.values()) for r in stored] == rows
    finally:
        conn.close()
