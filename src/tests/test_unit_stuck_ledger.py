"""OPS.234: a stuck unit of a run in the ledger is stopped through the supervisor (`controls`, `stop_unit`) and handed out again by the launcher, with no RQ job in it."""

import json
import signal
import types

import pytest

from planetgen.db import store
from planetgen.queue import beatjudge, commands, controls, ledger, planrows, redisqueue, runrecord, supervisor, unitbeat, unitstuck


@pytest.fixture
def conn(redis_server):
    connection = redisqueue.connect(redis_server)
    for key in connection.scan_iter(match=f"{runrecord.PREFIX}:*"):
        connection.delete(key)
    return connection


@pytest.fixture
def control(mysql_config):
    connection = store.get_control_connection(mysql_config, ensure_schema=True)
    connection.execute("DELETE FROM engine_runs")
    connection.commit()
    try:
        yield connection
    finally:
        connection.execute("DELETE FROM engine_runs")
        connection.commit()
        connection.close()


@pytest.fixture
def new_control(mysql_config):
    return lambda: store.get_control_connection(mysql_config, ensure_schema=True)


@pytest.fixture
def ledgered(conn, control, monkeypatch, mysql_config):
    """A live run of one stage with two units in the ledger, unit 1 running (attempt 1, a process `4001`); returns `(run id, engine run, {unit: job row id})`."""
    monkeypatch.setitem(commands.REGISTRY, "stuckkind", commands.CommandSpec("stuckkind", "A test", unit=(__name__, "_unit"), expand=(__name__, "_unit")))
    run = runrecord.start_run(conn, "T", [runrecord.Stage("stuckkind", "Stuck", units=2)])
    runrecord.begin_stage(conn, run, 0)
    key = planrows._stage_key(conn, run, 0)
    engine = ledger.add_run(control, "generate", "g1")
    ledger.add_jobs(control, engine, [{"stage": key, "command": "stuckkind", "item_key": str(n), "arguments": {"args": {"n": n}}, "claim_mb": 100} for n in range(2)])
    rows = {n: row["id"] for n, row in enumerate(control.execute("SELECT id FROM engine_jobs WHERE run_id = ? ORDER BY id", (engine,)).fetchall())}
    control.execute("UPDATE engine_jobs SET state = 'ready' WHERE run_id = ?", (engine,))
    control.execute("UPDATE engine_jobs SET state = 'running' WHERE id = ?", (rows[1],))
    ledger.begin_attempt(control, rows[1], pid=4001)
    ledger.set_run_state(control, engine, "running")
    control.commit()
    conn.hset(runrecord.run_key(run), planrows.ENGINE_RUN_FIELD, engine)
    monkeypatch.setattr(planrows, "control_connection", lambda *_a, **_k: store.get_control_connection(mysql_config, ensure_schema=True))
    return run, engine, rows


class Job:
    def __init__(self):
        self.killed = False

    def kill(self, number=signal.SIGKILL):
        self.killed = True


def _beat(conn, run, unit):
    now = unitbeat.server_time(conn)
    conn.hset(runrecord.progress_key(run, 0), unit, json.dumps({"done": 3.0, "estimate": 9.0, "beat": now, "started": now - 3000, "moved": now - 2900, "phase": "p",
                                                                "waiting_on": "", "restarts": 0, "pid": 999_999_999, "host": "elsewhere"}))


def _quiet():
    return beatjudge.Verdict(beatjudge.QUIET, reason="the count has stood at 3 for 2900 s")


def _follower(conn, run, errors=None):
    return unitstuck.Stuck(conn, run, 0, dump=lambda raw: "STACK", record_error=lambda message, **fields: (errors if errors is not None else []).append(message))


def _duty(conn, control, new_control, job):
    launcher = types.SimpleNamespace(control=control, running={})
    duty = controls.Controls(conn, launches=types.SimpleNamespace(launcher=launcher), control=new_control)
    return duty, launcher


def _row(control, row_id):
    return control.execute("SELECT state, attempt FROM engine_jobs WHERE id = ?", (row_id,)).fetchone()


def _attempt(control, row_id):
    return control.execute("SELECT outcome, error_kind, error_text, stack_file FROM engine_attempts WHERE job_id = ? ORDER BY attempt DESC LIMIT 1", (row_id,)).fetchone()


def test_a_stuck_unit_in_the_ledger_is_killed_by_the_supervisor_and_its_attempt_ends_stuck(conn, control, new_control, ledgered, tmp_path, monkeypatch):
    run, _engine, rows = ledgered
    monkeypatch.setattr("planetgen.queue.attemptlog.stack_path", lambda row, number: str(tmp_path / f"{row}-{number}.txt"))
    supervisor.RedisStore(conn).publish(4242, 1.0, 1.0, 1)
    duty, launcher = _duty(conn, control, new_control, None)
    job = Job()
    launcher.running[rows[1]] = types.SimpleNamespace(row={"id": rows[1], "run_id": _engine}, job=job)
    _beat(conn, run, 1)
    errors = []
    follower = _follower(conn, run, errors)
    follower.step([1], {1: _quiet()})
    record = follower.records()[1]
    assert record["phase"] == "stopping" and record["attempt"] == 1 and errors
    duty.turn()
    control.commit()
    assert job.killed
    attempt = _attempt(control, rows[1])
    assert (attempt["outcome"], attempt["error_kind"]) == ("stuck", "stuck") and "quiet" in attempt["error_text"]
    assert open(attempt["stack_file"], encoding="utf-8").read() == "STACK"
    assert _row(control, rows[1])["state"] == "running", "the launcher settles the end of the process and hands the unit out again"
    follower.step([], {})
    assert follower.records()[1]["phase"] == "stopping", "a stop under way is waited for"
    control.execute("UPDATE engine_jobs SET state = 'ready' WHERE id = ?", (rows[1],))     # what the launcher does when it sees the process end
    control.commit()
    _beat(conn, run, 1)
    follower.step([], {})
    assert follower.records()[1]["phase"] == "requeued"
    assert conn.hget(runrecord.progress_key(run, 0), 1) is None, "the beat of the stopped attempt is cleared"


def test_a_stuck_unit_whose_stop_never_takes_effect_fails_the_run(conn, control, ledgered):
    run, _engine, _rows = ledgered
    _beat(conn, run, 1)
    follower = _follower(conn, run)
    follower.step([1], {1: _quiet()})        # no supervisor is alive: nothing was asked
    assert follower.records()[1]["phase"] == "stopping"
    record = follower.records()[1]
    conn.hset(follower.key, 1, json.dumps({**record, "at": float(record["at"]) - unitstuck.STOP_WAIT_SECONDS - 1}))
    follower.step([], {})
    state = runrecord.read_run(conn, run)
    assert state["state"] == "failed" and "could not be stopped" in state["error"]


def test_a_dead_beat_of_a_ledgered_unit_that_is_ready_again_is_not_acted_on(conn, control, ledgered):
    run, _engine, rows = ledgered
    control.execute("UPDATE engine_jobs SET state = 'ready' WHERE id = ?", (rows[1],))
    control.commit()
    _beat(conn, run, 1)
    follower = _follower(conn, run)
    follower.step([1], {1: beatjudge.Verdict(beatjudge.DEAD, reason="the beat stopped")})
    assert follower.records() == {}, "its process is gone because the attempt ended and the launcher has it ready: it is waiting, not dead"


def test_a_stuck_request_for_a_unit_that_does_not_run_is_answered_with_no_jobs(conn, control, new_control, ledgered):
    run, _engine, _rows = ledgered
    duty, _launcher = _duty(conn, control, new_control, None)
    reply = duty._stop_unit({}, {"run_id": run, "stage": 0, "unit": 0, "reason": "r", "stack": ""}, "x")
    assert reply == {"ok": True, "jobs": 0}
