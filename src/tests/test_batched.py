"""
OPS.185: a long statement of a clear runs a batch of rows at a time (`db/batched.py`), commits each batch, tells the running unit it moved and says the rows done
so far; the clears of the scatters use it, so a unit that empties millions of rows never goes quiet for the length of one statement.
"""

import pytest

from planetgen.db import batched, store
from planetgen.queue import unitbeat

pytestmark = pytest.mark.db


@pytest.fixture
def moved(monkeypatch):
    calls = []
    monkeypatch.setattr(unitbeat, "moving", lambda: calls.append(1))
    return calls


def _scratch(conn, rows):
    conn.execute("CREATE TABLE IF NOT EXISTS batched_probe (id INT PRIMARY KEY, level DOUBLE NOT NULL)")
    conn.execute("DELETE FROM batched_probe")
    for number in range(1, rows + 1):
        conn.execute("INSERT INTO batched_probe (id, level) VALUES (?, ?)", (number, float(number)))
    conn.commit()


def test_a_delete_runs_in_batches_until_nothing_is_left(mysql_config, moved):
    conn = store.get_connection(mysql_config)
    try:
        _scratch(conn, 10)
        removed = batched.repeat(conn, "DELETE FROM batched_probe WHERE level > ?", (2.0,), label="Probe", batch=3)
        left = conn.execute("SELECT COUNT(*) AS n FROM batched_probe").fetchone()["n"]
    finally:
        conn.execute("DROP TABLE IF EXISTS batched_probe")
        conn.close()
    assert removed == 8 and left == 2
    assert len(moved) == 3, "8 rows in batches of 3 is three batches that changed rows"


def test_an_update_that_takes_a_row_out_of_its_own_where_runs_in_batches(mysql_config, moved):
    conn = store.get_connection(mysql_config)
    try:
        _scratch(conn, 7)
        changed = batched.repeat(conn, "UPDATE batched_probe SET level = -1 WHERE level > 0", label="Probe", batch=2)
        levels = {row["level"] for row in conn.execute("SELECT level FROM batched_probe").fetchall()}
    finally:
        conn.execute("DROP TABLE IF EXISTS batched_probe")
        conn.close()
    assert changed == 7 and levels == {-1.0} and len(moved) == 4


def test_a_statement_with_nothing_to_change_moves_nothing(mysql_config, moved):
    conn = store.get_connection(mysql_config)
    try:
        _scratch(conn, 0)
        assert batched.repeat(conn, "DELETE FROM batched_probe", label="Probe") == 0
    finally:
        conn.execute("DROP TABLE IF EXISTS batched_probe")
        conn.close()
    assert moved == []


def test_each_batch_is_said_by_the_ticker(mysql_config, moved):
    said = []

    class Pace:
        def tick(self, done, detail=None):
            said.append(done)

    conn = store.get_connection(mysql_config)
    try:
        _scratch(conn, 5)
        batched.repeat(conn, "DELETE FROM batched_probe", label="Probe", batch=2, pace=Pace())
    finally:
        conn.execute("DROP TABLE IF EXISTS batched_probe")
        conn.close()
    assert said == [2, 4, 5]


def test_the_indexes_of_the_phenomenon_scatter_come_off_one_statement_each(mysql_config, moved):
    conn = store.get_connection(mysql_config)
    try:
        store.add_phenomenon_scatter_indexes(conn)
        wanted = len(store.PHENOMENON_SCATTER_INDEXES)
        store.drop_phenomenon_scatter_indexes(conn)
        left = store.missing_phenomenon_scatter_indexes(conn)
    finally:
        conn.close()
    assert len(moved) == wanted and len(left) == wanted
