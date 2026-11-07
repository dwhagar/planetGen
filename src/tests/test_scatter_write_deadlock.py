"""DB.113: a unit that meets a MySQL deadlock (1213) or lock wait timeout (1205) tries again with dynamic waits, and after the eighth failed try asks for the step to itself; all logged."""
import pymysql
import pytest

from planetgen.db import contention
from planetgen.generation import scatter_companions
from planetgen.util import log


class _Conn:
    def __init__(self):
        self.log = []
        self._config = None

    def commit(self):
        self.log.append("commit")

    def rollback(self):
        self.log.append("rollback")

    def execute(self, sql, *args):
        self.log.append(sql)


class _Gate:
    """A gate with no database: closed while `closed` says so, and a grant that is given or refused."""

    def __init__(self, closed=(), grants=True):
        self.closed, self.grants, self.events = list(closed), grants, []

    def is_closed(self):
        return self.closed.pop(0) if self.closed else False

    def grant(self, need):
        self.events.append(("grant", str(need)))
        return self.grants

    def release(self):
        self.events.append(("release",))


def _error(code=1213):
    return pymysql.err.OperationalError(code, "Deadlock found when trying to get lock; try restarting transaction")


@pytest.fixture(autouse=True)
def _quiet(monkeypatch):
    sleeps, lines = [], []
    monkeypatch.setattr(contention.time, "sleep", sleeps.append)
    monkeypatch.setattr(log, "normal", lines.append)
    monkeypatch.setattr(contention, "_pressure", 0.0)
    return sleeps, lines


def _failing(times, code=1213):
    calls = []

    def attempt(conn):
        calls.append(1)
        if len(calls) <= times:
            raise _error(code)
        return "done"
    return attempt, calls


def test_a_deadlock_is_retried_with_a_logged_wait(_quiet):
    sleeps, lines = _quiet
    attempt, calls = _failing(2)
    gate = _Gate()
    assert contention.run_with_retries(_Conn(), attempt, "g", "Unit 1", gate_object=gate) == "done"
    assert len(calls) == 3 and len(sleeps) == 2 and gate.events == []
    assert sum("error 1213" in line and "trying again in" in line for line in lines) == 2


def test_a_lock_wait_timeout_is_retried_too():
    attempt, calls = _failing(1, code=1205)
    assert contention.run_with_retries(_Conn(), attempt, "g", "Unit 1", gate_object=_Gate()) == "done"


def test_the_waits_grow_but_never_pass_the_cap(monkeypatch):
    monkeypatch.setattr(contention, "_pressure", 1.0)
    wait = 0.0
    for number in range(1, 30):
        wait = contention.next_wait(number, wait)
        assert contention.BASE_WAIT_S <= wait <= contention.MAX_WAIT_S


def test_trouble_in_this_process_lengthens_the_waits(monkeypatch):
    monkeypatch.setattr(contention._rng, "uniform", lambda low, high: high)
    calm = contention.next_wait(1, 0.5)
    monkeypatch.setattr(contention, "_pressure", 1.0)
    assert contention.next_wait(1, 0.5) > calm


def test_after_the_eighth_failed_try_the_unit_is_granted_the_step_to_itself(_quiet):
    _sleeps, lines = _quiet
    attempt, calls = _failing(8)
    gate = _Gate()
    assert contention.run_with_retries(_Conn(), attempt, "scatter-ids", "Unit 1", gate_object=gate) == "done"
    assert len(calls) == 9
    assert [event[0] for event in gate.events] == ["grant", "release"]
    assert "needs exclusive use" in gate.events[0][1]
    text = "\n".join(lines)
    assert "asking for it" in text and "granted exclusive use" in text and "released exclusive use" in text


def test_the_grant_is_released_when_the_granted_try_fails_too():
    attempt, calls = _failing(99)
    gate = _Gate()
    with pytest.raises(pymysql.err.OperationalError):
        contention.run_with_retries(_Conn(), attempt, "g", "Unit 1", gate_object=gate)
    assert len(calls) == 9 and [event[0] for event in gate.events] == ["grant", "release"]


def test_a_grant_that_does_not_come_gives_up_with_a_lock_error():
    attempt, calls = _failing(99)
    with pytest.raises(pymysql.err.OperationalError) as caught:
        contention.run_with_retries(_Conn(), attempt, "g", "Unit 1", gate_object=_Gate(grants=False))
    assert caught.value.args[0] == 1205 and len(calls) == 8


def test_a_unit_holds_off_while_another_has_the_grant(_quiet):
    sleeps, lines = _quiet
    attempt, calls = _failing(0)
    contention.run_with_retries(_Conn(), attempt, "g", "Unit 1", gate_object=_Gate(closed=[True, True, False]))
    assert len(sleeps) == 2 and any("holds off" in line for line in lines) and any("goes on after" in line for line in lines)


def test_any_other_database_error_is_not_retried():
    calls = []

    def attempt(conn):
        calls.append(1)
        raise pymysql.err.OperationalError(1146, "Table doesn't exist")

    with pytest.raises(pymysql.err.OperationalError):
        contention.run_with_retries(_Conn(), attempt, "g", "Unit 1", gate_object=_Gate())
    assert len(calls) == 1


def test_a_scatter_batch_is_rolled_back_and_written_again(monkeypatch):
    monkeypatch.setattr(contention, "MySQLGate", lambda conn, gate: _Gate())
    conn, calls = _Conn(), []

    def write(c, rows, epoch):
        calls.append(list(rows))
        if len(calls) == 1:
            raise _error()
        return scatter_companions.Written(len(rows), 1)

    assert scatter_companions.write_committed(conn, write, [1, 2, 3], 0.0) == (3, 1)
    assert calls == [[1, 2, 3], [1, 2, 3]]
    assert "rollback" in conn.log and conn.log[-1] == "commit"
    assert any("READ COMMITTED" in entry for entry in conn.log if entry not in ("commit", "rollback"))
