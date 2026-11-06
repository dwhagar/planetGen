"""The scatter units' batch writes start over after a MySQL deadlock (error 1213) instead of stopping the stage (a deadlock killed unit 11 of the massive-star stage)."""
import pymysql
import pytest

from planetgen.db import store
from planetgen.generation import scatter_companions


class _Conn:
    """Records the statements and commits of a connection; no database needed."""

    def __init__(self):
        self.log = []

    def commit(self):
        self.log.append("commit")

    def rollback(self):
        self.log.append("rollback")

    def execute(self, sql, *args):
        self.log.append(sql)


def _deadlock(code=1213):
    return pymysql.err.OperationalError(code, "Deadlock found when trying to get lock; try restarting transaction")


@pytest.fixture(autouse=True)
def _no_sleep(monkeypatch):
    monkeypatch.setattr(scatter_companions.time, "sleep", lambda seconds: None)


def test_a_batch_that_meets_a_deadlock_is_rolled_back_and_written_again():
    conn, calls = _Conn(), []

    def write(c, rows, epoch):
        calls.append(list(rows))
        if len(calls) == 1:
            raise _deadlock()
        return scatter_companions.Written(len(rows), 1)

    written = scatter_companions.write_committed(conn, write, [1, 2, 3], 0.0)
    assert written == (3, 1)
    assert calls == [[1, 2, 3], [1, 2, 3]]
    assert "rollback" in conn.log and conn.log[-1] == "commit"
    assert any("READ COMMITTED" in entry for entry in conn.log if entry not in ("commit", "rollback"))


def test_a_lock_wait_timeout_is_retried_too():
    conn, calls = _Conn(), []

    def write(c, rows, epoch):
        calls.append(1)
        if len(calls) < 3:
            raise _deadlock(1205)
        return scatter_companions.Written(1, 0)

    assert scatter_companions.write_committed(conn, write, [1], 0.0).rows == 1


def test_a_deadlock_that_persists_gives_up_after_the_attempt_limit():
    calls = []

    def write(c, rows, epoch):
        calls.append(1)
        raise _deadlock()

    with pytest.raises(pymysql.err.OperationalError):
        scatter_companions.write_committed(_Conn(), write, [1], 0.0)
    assert len(calls) == store.SECTOR_SAVE_ATTEMPTS


def test_any_other_database_error_is_not_retried():
    calls = []

    def write(c, rows, epoch):
        calls.append(1)
        raise pymysql.err.OperationalError(1146, "Table doesn't exist")

    with pytest.raises(pymysql.err.OperationalError):
        scatter_companions.write_committed(_Conn(), write, [1], 0.0)
    assert len(calls) == 1
