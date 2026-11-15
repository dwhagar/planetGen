"""OPS.220: failure by kind (`planetgen/queue/failures.py`). The classifier names a failure passing, killed, idle or lasting; a passing one is run again by the engine after a doubling, spread wait inside the
run's budget (every attempt in the ledger); anything else is not run again."""

import json
import random

import pymysql
import pytest
import redis

from planetgen.db import store
from planetgen.queue import commands, failures, planrows, redisqueue, runrecord, units, wrapper
from tests.test_attempts import _control, _ledger_run
from tests.test_work_units import _handlers

pytestmark = pytest.mark.fresh_database


def _mysql(code):
    return pymysql.err.OperationalError(code, "x")


@pytest.mark.parametrize("exc, kind", [
    (_mysql(1213), failures.PASSING), (_mysql(1205), failures.PASSING), (_mysql(2006), failures.PASSING), (_mysql(2013), failures.PASSING),
    (pymysql.err.InterfaceError(0, "closed"), failures.PASSING),
    (redis.exceptions.ConnectionError("gone"), failures.PASSING), (redis.exceptions.TimeoutError("slow"), failures.PASSING),
    (_mysql(1062), failures.LASTING), (_mysql(1146), failures.LASTING), (ValueError("bug"), failures.LASTING), (KeyError("k"), failures.LASTING),
    ("JobTimeoutException", failures.KILLED), ("WorkerLostError", failures.KILLED), ("AbandonedJobError", failures.KILLED), ("UnitStuck", failures.IDLE),
])
def test_every_failure_is_one_kind(exc, kind):
    assert failures.classify(exc) == kind


def test_the_waits_double_from_two_seconds_to_five_minutes_and_are_spread():
    rng = random.Random(1)
    for retry, ceiling in [(1, 2.0), (2, 4.0), (3, 8.0), (6, 64.0), (8, 256.0), (9, 300.0), (30, 300.0)]:
        waits = [failures.wait_seconds(retry, rng) for _ in range(200)]
        assert ceiling * 0.5 <= min(waits) and max(waits) <= ceiling, retry
        assert max(waits) - min(waits) > ceiling * 0.2, "the waits are spread, not one number"


def test_the_budget_of_a_run_is_a_number_of_retries_and_seconds_of_waiting(redis_server):
    conn = redisqueue.connect()
    run = "budget-run"
    conn.delete(runrecord.failures_key(run))
    assert all(failures.charge(conn, run, 1.0) for _ in range(failures.RETRIES_PER_RUN))
    assert not failures.charge(conn, run, 1.0), "no retry left"
    assert failures.spent(conn, run) == (failures.RETRIES_PER_RUN, float(failures.RETRIES_PER_RUN))
    other = "budget-run-2"
    conn.delete(runrecord.failures_key(other))
    assert failures.charge(conn, other, failures.WAIT_PER_RUN_SECONDS - 1.0) and not failures.charge(conn, other, 2.0), "waiting is limited too"
    assert failures.charge(conn, other, 1.0)
    conn.delete(runrecord.failures_key(run), runrecord.failures_key(other))


def test_only_a_passing_failure_inside_the_budget_gets_a_wait(redis_server):
    conn = redisqueue.connect()
    run = "wait-run"
    conn.delete(runrecord.failures_key(run))
    assert failures.retry_wait(conn, run, _mysql(1213), 0) is not None
    assert failures.retry_wait(conn, run, ValueError("bug"), 0) is None
    assert failures.retry_wait(conn, run, KeyboardInterrupt(), 0) is None
    conn.delete(runrecord.failures_key(run))


# --- through the unit ---------------------------------------------------------------------------------------------------------------------------------

def _started(monkeypatch, mysql_config, handler):
    """A one-unit run (stage `test`) whose unit is `handler`, started and not run; returns `(conn, run_id, row, args)`."""
    _handlers(monkeypatch, test=handler)
    database = wrapper.database_args(mysql_config)
    run_id = runrecord.start_run(redisqueue.connect(), "Failures", [runrecord.Stage("test", "T", units=1)],
                                 steps=[commands.Step("test", "test", "T", {"database": database}).as_dict()])
    conn = redisqueue.connect()
    args = {"n": 1, "database": database}
    rows = planrows.record_cut(conn, run_id, 0, "test", [(args, 1.0)], 1.0, {})
    return conn, run_id, rows[0], args


def _flaky(errors):
    """A unit handler that raises each of `errors` in turn, then answers."""
    pending = list(errors)

    def handler(args, report):
        report(1, 1)
        if pending:
            raise pending.pop(0)
        return {"n": args["n"]}

    return handler


def _outcomes(mysql_config, row):
    control = _control(mysql_config)
    try:
        return [(item["outcome"], item["error_text"]) for item in control.execute("SELECT outcome, error_text FROM engine_attempts WHERE job_id = ? ORDER BY attempt", (row,)).fetchall()]
    finally:
        control.close()


def test_a_passing_failure_is_run_again_after_a_wait_and_every_attempt_is_in_the_ledger(monkeypatch, mysql_config):
    conn, run_id, row, args = _started(monkeypatch, mysql_config, _flaky([_mysql(1213), redis.exceptions.ConnectionError("gone")]))
    waits = []
    monkeypatch.setattr(units, "SLEEP", waits.append)
    result = units._do_unit(conn, run_id, 0, 0, "test", args, 1.0, row)
    assert "seconds" in result and conn.sismember(runrecord.done_key(run_id, 0), 0)
    assert [outcome for outcome, _text in _outcomes(mysql_config, row)] == ["failed", "failed", "done"]
    assert all(text.startswith("[passing]") for outcome, text in _outcomes(mysql_config, row) if outcome == "failed")
    assert 1.0 <= sum(waits[:2]) <= 2.0 and waits[2:] and sum(waits[2:]) >= 2.0, "the first wait is up to 2 s, the second up to 4 s, taken a second at a time"
    assert failures.spent(conn, run_id)[0] == 2


def test_a_lasting_failure_is_not_run_again(monkeypatch, mysql_config):
    conn, run_id, row, args = _started(monkeypatch, mysql_config, _flaky([ValueError("a bug")]))
    monkeypatch.setattr(units, "SLEEP", lambda seconds: pytest.fail("a lasting failure must not wait"))
    with pytest.raises(ValueError):
        units._do_unit(conn, run_id, 0, 0, "test", args, 1.0, row)
    assert [outcome for outcome, _text in _outcomes(mysql_config, row)] == ["failed"]
    assert _outcomes(mysql_config, row)[0][1].startswith("[lasting]")
    assert failures.spent(conn, run_id) == (0, 0.0)


def test_a_passing_failure_that_outlasts_the_runs_budget_ends_the_unit_and_says_so(monkeypatch, mysql_config):
    conn, run_id, row, args = _started(monkeypatch, mysql_config, _flaky([_mysql(2013)] * 99))
    monkeypatch.setattr(units, "SLEEP", lambda seconds: None)
    monkeypatch.setattr(failures, "RETRIES_PER_RUN", 3)
    with pytest.raises(pymysql.err.OperationalError):
        units._do_unit(conn, run_id, 0, 0, "test", args, 1.0, row)
    assert [outcome for outcome, _text in _outcomes(mysql_config, row)] == ["failed"] * 4, "the first run and three more"
    assert "could not clear" in failures.describe(conn, run_id, _mysql(2013)) and "3 time(s)" in failures.describe(conn, run_id, _mysql(2013))
    assert failures.describe(conn, run_id, ValueError("x")) == ""


def test_a_run_that_ends_while_a_unit_waits_to_run_again_stops_the_unit(monkeypatch, mysql_config):
    conn, run_id, row, args = _started(monkeypatch, mysql_config, _flaky([_mysql(1205)]))
    monkeypatch.setattr(units, "SLEEP", lambda seconds: runrecord.finish_run(conn, run_id, "canceled"))
    assert units._do_unit(conn, run_id, 0, 0, "test", args, 1.0, row) == {"stopped": True}
    assert not conn.sismember(runrecord.done_key(run_id, 0), 0)


def test_a_unit_that_does_not_fail_takes_nothing_from_the_budget(monkeypatch, mysql_config):
    conn, run_id, row, args = _started(monkeypatch, mysql_config, _flaky([]))
    monkeypatch.setattr(units, "SLEEP", lambda seconds: pytest.fail("no wait"))
    assert "seconds" in units._do_unit(conn, run_id, 0, 0, "test", args, 1.0, row)
    assert failures.spent(conn, run_id) == (0, 0.0)


# --- a unit run again continues from its cursor (OPS.261, OPS.219) -------------------------------------------------------------------------------------

RESUME_PIECES = 6
DONE_PIECES = []
FAIL_AT = {}


def _pieces(args, cursor):
    return list(range(args["count"]))[cursor:]


def _piece(args, piece, ctx):
    DONE_PIECES.append(piece)
    if FAIL_AT.get("piece") == piece:
        FAIL_AT.clear()
        raise _mysql(1213)
    ctx.db.execute("INSERT INTO retry_probe (piece) VALUES (?)", (piece,))
    return {"piece": piece}


def _together(args, results, ctx):
    return {"results": results}


@pytest.mark.parametrize("failed", [0, 3, RESUME_PIECES - 1])
def test_a_unit_run_again_after_a_passing_failure_goes_on_after_the_pieces_it_finished(monkeypatch, mysql_config, failed):
    """The end of the piece-0 restart: the second attempt of the unit reads its cursor and does the pieces from the one that failed, writes no row twice, and gives the result an undisturbed unit gives."""
    DONE_PIECES.clear()
    FAIL_AT.clear()
    monkeypatch.setitem(commands.REGISTRY, "test", commands.CommandSpec("test", "A test", expand=(__name__, "_pieces"), pieces=(__name__, "_pieces"), do=(__name__, "_piece"),
                                                                       gather=(__name__, "_together"), resumes=True))
    conn = redisqueue.connect()
    database = wrapper.database_args(mysql_config)
    run_id = runrecord.start_run(conn, "Failures", [runrecord.Stage("test", "T", units=1)], steps=[commands.Step("test", "test", "T", {"database": database}).as_dict()])
    args = {"count": RESUME_PIECES, "database": database}
    row = planrows.record_cut(conn, run_id, 0, "test", [(args, 1.0)], 1.0, {})[0]
    galaxy = store.get_connection(mysql_config)
    try:
        galaxy.execute("CREATE TABLE IF NOT EXISTS retry_probe (piece INT NOT NULL)")
        galaxy.execute("DELETE FROM retry_probe")
        galaxy.commit()
        FAIL_AT["piece"] = failed
        monkeypatch.setattr(units, "SLEEP", lambda seconds: None)
        result = units._do_unit(conn, run_id, 0, 0, "test", args, 1.0, row)
        assert "seconds" in result
        assert DONE_PIECES == [*range(failed + 1), *range(failed, RESUME_PIECES)], "the pieces before the failed one ran once, the failed one twice, the rest once"
        galaxy.commit()
        assert [r["piece"] for r in galaxy.execute("SELECT piece FROM retry_probe ORDER BY piece").fetchall()] == list(range(RESUME_PIECES)), "no piece wrote its data twice and none was lost"
        assert [outcome for outcome, _text in _outcomes(mysql_config, row)] == ["failed", "done"]
        stored = json.loads(conn.hget(runrecord.results_key(run_id, 0), 0))
        assert stored["results"] == [{"piece": n} for n in range(RESUME_PIECES)], "the pieces done before the failure are in the result, as in an undisturbed unit's"
    finally:
        galaxy.execute("DROP TABLE IF EXISTS retry_probe")
        galaxy.commit()
        galaxy.close()
