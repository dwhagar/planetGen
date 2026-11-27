"""OPS.284: a run's record is kept in the control database (`queue/runstore.py`), so a run Redis lost is put back, ended as a run whose process is gone, and carried on by `finish-run`.
Real MySQL and Redis; the command is the test's own (`tests/test_window_spawn.py`)."""

import pytest

from planetgen.db import store
from planetgen.generation import finish_run
from planetgen.queue import commands, controls, planrows, redisqueue, runrecord, runstore, supervisor, wrapper
from tests.test_plan_rows import _jobs
from tests.test_window_spawn import TOTAL, WINDOW, _cut, _rows, _units_of, _windowed   # noqa: F401 -- the fixture and helpers

pytestmark = [pytest.mark.fresh_database]

SEED, BUILD = "AB" * 16, "77F01600909B06030C0300"


def _empty_copies(mysql_config):
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        control.execute("DELETE FROM engine_run_records")
        control.commit()
    finally:
        control.close()


@pytest.fixture(autouse=True)
def _no_copies_left_behind(mysql_config):
    """A copy this file's runs leave in the control database would be put back into a later test's Redis (that is what it is for): the file starts and ends with none."""
    _empty_copies(mysql_config)
    yield
    _empty_copies(mysql_config)
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        controls.end_rows(control)          # the rows of a run whose Redis this file dropped stay ready: the inline launcher of a later test would start them and write a hash for each run
        control.commit()
    finally:
        control.close()


def _started(mysql_config):
    conn = redisqueue.connect()
    step = commands.Step("a", "windowed", "A", {"database": wrapper.database_args(mysql_config)})
    run_id = runrecord.start_run(conn, "Windows", [runrecord.Stage("a", "A", units=TOTAL)], steps=[step.as_dict()], database=mysql_config.database, galaxy_seed=SEED, build=BUILD)
    planrows.record_plan(conn, run_id, [step])
    first = _cut({"window": WINDOW})
    planrows.record_cut(conn, run_id, 0, "windowed", first.units, float(TOTAL), first.plan, total=TOTAL, later=first.later, window=WINDOW)
    planrows.release_units(conn, run_id, 0)
    runrecord.begin_stage(conn, run_id, 0)
    return conn, run_id, step


def _lose_redis(conn, *run_ids):
    """What a restart of the machine does to a Redis that keeps nothing, for the runs of this test (the other tests' runs are left where they are)."""
    for run_id in run_ids:
        conn.delete(*runrecord._keys_of(conn, run_id, every=True))


def _ledger_state(conn, run_id, mysql_config):
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        engine = int(runrecord._s(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)))
        return engine, control.execute("SELECT state FROM engine_runs WHERE id = ?", (engine,)).fetchone()["state"]
    finally:
        control.close()


def test_a_run_with_rows_in_the_ledger_has_a_durable_copy_that_follows_its_changes(mysql_config):
    conn, run_id, _step = _started(mysql_config)
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        record = runstore.stored(control, run_id)
        assert record["run"]["galaxy_seed"] == SEED and record["run"]["build"] == BUILD and record["run"]["state"] == "running"
        assert "started" in record["stages"][0], "a stage that began is in the copy"
        runrecord.set_conclusion(conn, run_id, 0, {"made": 3})
        control.commit()            # a new look, not the snapshot of the first
        assert '"made": 3' in runstore.stored(control, run_id)["stages"][0]["conclusion"]
        runrecord.hold_run(conn, run_id, "unit 3 failed")
        control.commit()
        assert runstore.stored(control, run_id)["run"]["state"] == "held"
    finally:
        control.close()


def test_a_run_without_rows_in_the_ledger_has_no_copy(mysql_config):
    conn = redisqueue.connect()
    run_id = runrecord.start_run(conn, "By hand", [runrecord.Stage("a", "A", units=1)], database=mysql_config.database)
    runrecord.begin_stage(conn, run_id, 0)
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        assert runstore.stored(control, run_id) is None
    finally:
        control.close()


def test_a_run_redis_lost_comes_back_and_a_live_one_is_ended_with_its_ledger_rows(mysql_config):
    conn, run_id, _step = _started(mysql_config)
    before = runrecord.read_run(conn, run_id)
    _lose_redis(conn, run_id)
    assert runrecord.read_run(conn, run_id) is None
    assert runstore.settle_lost(conn) == [run_id]
    after = runrecord.read_run(conn, run_id)
    assert (after["title"], after["galaxy_seed"], after["build"], after["database"]) == (before["title"], SEED, BUILD, before["database"])
    assert after["state"] == "failed" and after["killed"] and "restarted" in after["error"]
    engine, state = _ledger_state(conn, run_id, mysql_config)
    assert state == "canceled", "the ledger's run row is ended, not left running"
    rows = _jobs(mysql_config, engine)
    assert {job["state"] for job in rows} <= {"done", "canceled", "failed"}
    assert runstore.settle_lost(conn) == [], "a run that is back is not restored or ended again"


def test_a_held_run_comes_back_held_and_is_left_alone(mysql_config):
    conn, run_id, _step = _started(mysql_config)
    runrecord.hold_run(conn, run_id, "unit 4 failed: boom")
    _lose_redis(conn, run_id)
    assert runstore.settle_lost(conn) == []
    assert runrecord.read_run(conn, run_id)["state"] == "held"


def test_finish_run_carries_on_a_run_whose_redis_was_lost(mysql_config, monkeypatch):
    """The case of the lost trial: the machine restarted in the middle of a stage; `finish-run` finds the run in the control database and plans it from its first stage that did not finish."""
    monkeypatch.setattr(finish_run, "FINISHABLE", finish_run.FINISHABLE | {"windowed"})
    conn, run_id, _step = _started(mysql_config)
    _lose_redis(conn, run_id)
    found = finish_run.find(conn, SEED, BUILD, mysql_config.database)
    assert found.run["id"] == run_id and found.first == 0 and found.steps[0].params[finish_run.resumed.RESUMES_FROM] == run_id


def test_the_supervisors_duty_ends_a_lost_run_at_most_every_half_minute(mysql_config):
    conn, run_id, _step = _started(mysql_config)
    _lose_redis(conn, run_id)
    now = [1000.0]
    duty = supervisor.LostRuns(conn, clock=lambda: now[0])
    assert duty.turn() == [run_id]
    assert duty.turn() == []
    _lose_redis(conn, run_id)
    now[0] += supervisor.LOST_RUNS_SECONDS - 1
    assert duty.turn() == [] and runrecord.read_run(conn, run_id) is None, "not yet"
    now[0] += 2
    assert duty.turn() == [] and runrecord.read_run(conn, run_id)["state"] == "failed", "it is back as it ended, and ended once"


def test_a_run_the_queue_page_deleted_is_not_restored_and_a_finished_run_comes_back_only_while_redis_would_keep_it(mysql_config, monkeypatch):
    conn, run_id, _step = _started(mysql_config)
    runrecord.finish_run(conn, run_id, "failed", error="boom", killed=True)
    assert runrecord.delete_run(conn, run_id)
    assert runstore.restore_lost(conn) == [], "deleted, so not restored"
    conn, kept, _step = _started(mysql_config)
    runrecord.finish_run(conn, kept, "failed", error="boom", killed=True)
    _lose_redis(conn, kept)
    assert runstore.restore_lost(conn) == [kept], "failed a moment ago: Redis would have it"
    _lose_redis(conn, kept)
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        control.execute("UPDATE engine_run_records SET updated_at = NOW(6) - INTERVAL 31 DAY WHERE run_id = ?", (kept,))
        control.commit()
        assert runstore.restore_lost(conn, control) == [], "a finished run older than Redis keeps it is not restored"
    finally:
        control.close()


def test_a_run_that_ended_in_redis_with_its_ledger_rows_open_has_them_ended_and_a_live_run_is_left_alone(mysql_config):
    """OPS.284: a run whose process could not even write its own end (a full disk) or was found dead by `finish-run` ends in Redis, and nothing ended its ledger rows (the OPS.247 trial's run 13)."""
    conn, ended_run, _step = _started(mysql_config)
    runrecord.finish_run(conn, ended_run, "failed", error="its process was gone", killed=True)
    live_run = _started(mysql_config)[1]
    engine, state = _ledger_state(conn, ended_run, mysql_config)
    assert state in ("running", "waiting"), "the end in Redis leaves the ledger's run row open"
    assert runstore.settle_ended(conn) == [ended_run]
    assert _ledger_state(conn, ended_run, mysql_config) == (engine, "canceled")
    assert {job["state"] for job in _jobs(mysql_config, engine)} <= {"done", "canceled", "failed"}
    assert _ledger_state(conn, live_run, mysql_config)[1] in ("running", "waiting"), "a live run's rows are not touched"
    assert runstore.settle_ended(conn) == [], "ended once"
    runrecord.finish_run(conn, live_run, "canceled", error="the test is over")
    runrecord.delete_run(conn, live_run)
    runrecord.delete_run(conn, ended_run)


def _set_jobs(mysql_config, engine, state, failed_error=None):
    """Puts every job of the ledger run `engine` into `state`; the first of them is `failed` with `failed_error` kept on its attempt when that is given."""
    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        control.execute("UPDATE engine_jobs SET state = ? WHERE run_id = ?", (state, engine))
        if failed_error is not None:
            first = control.execute("SELECT id, attempt FROM engine_jobs WHERE run_id = ? ORDER BY id LIMIT 1", (engine,)).fetchone()
            control.execute("UPDATE engine_jobs SET state = 'failed', attempt = 1 WHERE id = ?", (first["id"],))
            control.execute("INSERT INTO engine_attempts (job_id, attempt, started_at, ended_at, outcome, error_kind, error_text) VALUES (?, 1, NOW(6), NOW(6), 'failed', 'ValueError', ?)",
                            (first["id"], failed_error))
        control.commit()
    finally:
        control.close()


def test_a_run_nothing_can_run_in_ends_failed_with_the_cause_after_ten_minutes(mysql_config):
    """OPS.314: Boss's Run 57 stood Running for a day at 21 of 233 jobs, 1 failed, no worker alive; the supervisor ends it and says why."""
    conn, run_id, _step = _started(mysql_config)
    engine = _ledger_state(conn, run_id, mysql_config)[0]
    _set_jobs(mysql_config, engine, "waiting", failed_error="ValueError: the air solve did not converge\nTraceback ...")
    runstore._idle_since.clear()
    assert runstore.settle_stuck(conn, now=1000.0) == [], "first seen: not yet"
    assert runstore.settle_stuck(conn, now=1000.0 + runstore.STUCK_FAILED_SECONDS - 1) == []
    assert runrecord.read_run(conn, run_id)["state"] in runrecord.LIVE_STATES
    assert runstore.settle_stuck(conn, now=1000.0 + runstore.STUCK_FAILED_SECONDS + 1) == [run_id]
    after = runrecord.read_run(conn, run_id)
    assert after["state"] == "failed"
    assert "nothing could run for 10m" in after["error"] and "the air solve did not converge" in after["error"] and "Traceback" not in after["error"]
    assert {job["state"] for job in _jobs(mysql_config, engine)} <= {"done", "canceled", "failed"}
    assert runstore.settle_stuck(conn, now=5000.0) == [], "ended once"


def test_a_run_with_a_job_running_or_ready_is_left_alone_however_long_it_stands(mysql_config):
    conn, run_id, _step = _started(mysql_config)
    engine = _ledger_state(conn, run_id, mysql_config)[0]
    runstore._idle_since.clear()
    for state in ("ready", "running"):
        _set_jobs(mysql_config, engine, state)
        assert runstore.settle_stuck(conn, now=1000.0) == []
        assert runstore.settle_stuck(conn, now=1000.0 + 10 * runstore.STUCK_WAITING_SECONDS) == []
        assert runrecord.read_run(conn, run_id)["state"] in runrecord.LIVE_STATES
    assert run_id not in runstore._idle_since, "a run that can move is not counted as idle"
    runrecord.finish_run(conn, run_id, "canceled", error="the test is over")
    runrecord.delete_run(conn, run_id)
