# tests/test_work_queue_failures.py

"""
The work queue (`planetgen/queue/work.py`) when things go wrong:

- TEST.20, failure paths: a worker process dies, `on_done` raises, a
  payload or result won't pickle, a result isn't JSON, two tasks share a
  key, the heartbeat fails, the control database drops mid-run, and the
  lease under clock skew.
- TEST.21, cancelling a run: SIGTERM, Ctrl+C and `SystemExit` during a
  parallel run end it as cancelled, free the lease and leave no
  half-written sector.

Every run ends with its job row finished, no task row left `running` and
the lease free, whatever stopped it.
"""

import os
import signal
import subprocess
import sys
import threading
import time

import pytest

from planetgen.queue import redisqueue
from planetgen.queue import work as workQueue
from planetgen.db import store as _db

from tests.test_galaxy_gen import _mysql_argv, _plan_wide_galaxy
from tests.test_work_queue import _rows, control_config  # noqa: F401 -- fixture

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _square(payload):
    return payload * payload


def _slow_square(payload):
    time.sleep(0.2)
    return payload * payload


def _die_on_two(payload):
    if payload == 2:
        os._exit(7)  # the worker process dies, as an OOM kill would
    time.sleep(0.05)
    return payload


def _a_set(payload):
    return {payload, "x"}


def _a_generator(payload):
    return (n for n in range(payload))


def _ended_cleanly(config, job_id, state):
    """The run's row ended in `state`, none of its tasks still claims to
    run, and the lease is free."""
    job = _rows(config, "SELECT state, finished_at FROM work_jobs WHERE id = ?", (job_id,))[0]
    assert job["state"] == state and job["finished_at"] is not None
    running = _rows(config, "SELECT COUNT(*) AS n FROM work_tasks WHERE job_id = ?"
                            " AND state IN ('queued', 'running')", (job_id,))[0]["n"]
    assert running == 0
    lease = _rows(config, "SELECT holder FROM work_lease")
    assert not lease or lease[0]["holder"] is None


# ---------------------------------------------------------------------------
# TEST.20 Failure paths
# ---------------------------------------------------------------------------

def test_a_dead_worker_fails_the_run_and_frees_the_lease(control_config):
    """PERF.24: RQ runs a task whose worker died once more on a fresh one;
    when that one dies too, the run fails with `WorkerDied`."""
    with pytest.raises(redisqueue.WorkerDied):
        with workQueue.WorkQueue("Worker dies", workers=2, control_config=control_config) as queue:
            for n in range(8):
                queue.submit("die", f"n{n}", _die_on_two, n)
    _ended_cleanly(control_config, queue.job_id, "failed")
    errors = _rows(control_config, "SELECT error FROM work_tasks WHERE job_id = ? AND state = 'failed'",
                   (queue.job_id,))
    assert errors and any("WorkerDied" in row["error"] for row in errors)



@pytest.mark.parametrize("attempt", range(3))
def test_a_dead_worker_never_hangs_the_run(control_config, attempt):
    """PERF.22: a dead worker fails the run, with the other worker busy
    when it dies (the case that hung the old process pool on 3.12)."""
    finished = threading.Event()
    outcome = {}

    def run():
        try:
            with workQueue.WorkQueue("Worker dies, other busy", workers=2, control_config=control_config) as queue:
                for n in range(12):
                    queue.submit("die", f"n{n}", _die_on_two, n)
        except BaseException as exc:  # noqa: BLE001 -- checked below
            outcome["raised"] = exc
        finally:
            finished.set()

    thread = threading.Thread(target=run, daemon=True)
    thread.start()
    assert finished.wait(60), "the run hung after a worker died"
    assert isinstance(outcome.get("raised"), redisqueue.WorkerDied)


def test_on_done_raising_fails_the_run(control_config):
    def broken(result, seconds, weight):
        raise RuntimeError("the progress bar broke")

    with pytest.raises(RuntimeError, match="progress bar broke"):
        with workQueue.WorkQueue("on_done raises", workers=2, control_config=control_config) as queue:
            for n in range(6):
                queue.submit("square", f"n{n}", _slow_square, n, on_done=broken)
    _ended_cleanly(control_config, queue.job_id, "failed")


def test_on_done_raising_with_one_worker(control_config):
    def broken(result, seconds, weight):
        raise RuntimeError("the progress bar broke")

    with pytest.raises(RuntimeError):
        with workQueue.WorkQueue("on_done raises here", workers=1, control_config=control_config) as queue:
            queue.submit("square", "only", _square, 3, on_done=broken)
    _ended_cleanly(control_config, queue.job_id, "failed")


def test_a_payload_that_wont_pickle_fails_its_task(control_config):
    with pytest.raises(Exception) as raised:
        with workQueue.WorkQueue("Unpicklable payload", workers=2, control_config=control_config) as queue:
            queue.submit("square", "fine", _square, 2)
            queue.submit("square", "lock", _square, threading.Lock())
    assert "pickle" in str(raised.value).lower() or "lock" in str(raised.value).lower()
    _ended_cleanly(control_config, queue.job_id, "failed")


def test_a_result_that_wont_pickle_fails_its_task(control_config):
    with pytest.raises(Exception) as raised:
        with workQueue.WorkQueue("Unpicklable result", workers=2, control_config=control_config) as queue:
            queue.submit("gen", "gen", _a_generator, 3)
    assert "pickle" in str(raised.value).lower() or "generator" in str(raised.value).lower()
    _ended_cleanly(control_config, queue.job_id, "failed")


@pytest.mark.parametrize("workers", [1, 2])
def test_a_result_that_isnt_json_is_recorded_as_text(control_config, workers):
    results = []
    with workQueue.WorkQueue("Not JSON", workers=workers, control_config=control_config) as queue:
        queue.submit("set", "s", _a_set, 5, on_done=lambda result, *_: results.append(result))
    assert results == [{5, "x"}]
    row = _rows(control_config, "SELECT state, result FROM work_tasks WHERE job_id = ?", (queue.job_id,))[0]
    assert row["state"] == "done" and row["result"]
    _ended_cleanly(control_config, queue.job_id, "done")


def test_two_tasks_with_one_key_both_run_with_one_seed(control_config):
    from tests.test_work_queue import _draw

    results = []
    with workQueue.WorkQueue("Same key", workers=2, control_config=control_config, run_seed=11) as queue:
        for _ in range(2):
            queue.submit("draw", "1,2,3", _draw, None, on_done=lambda result, *_: results.append(result))
    assert len(results) == 2 and results[0] == results[1]
    rows = _rows(control_config, "SELECT task_key, state FROM work_tasks WHERE job_id = ?", (queue.job_id,))
    assert [(row["task_key"], row["state"]) for row in rows] == [("1,2,3", "done")] * 2
    _ended_cleanly(control_config, queue.job_id, "done")


def test_a_failing_heartbeat_doesnt_stop_the_run(control_config, monkeypatch):
    monkeypatch.setattr(workQueue, "HEARTBEAT_SECONDS", 0.05)
    beats = []

    def broken_heartbeat(self, job_id, holder):
        beats.append(job_id)
        raise ConnectionError("control database went away")

    monkeypatch.setattr(workQueue._ControlStore, "heartbeat", broken_heartbeat)
    results = []
    with workQueue.WorkQueue("Heartbeat fails", workers=2, control_config=control_config) as queue:
        for n in range(6):
            queue.submit("square", f"n{n}", _slow_square, n, on_done=lambda result, *_: results.append(result))
    assert sorted(results) == [n * n for n in range(6)]
    assert beats
    _ended_cleanly(control_config, queue.job_id, "done")


@pytest.mark.parametrize("workers", [1, 2])
def test_the_control_database_dropping_mid_run_doesnt_stop_it(control_config, workers):
    results = []

    def done(result, seconds, weight):
        results.append(result)
        if len(results) == 2:
            def gone():
                raise ConnectionError("control database went away")
            queue._store._connect = gone

    with workQueue.WorkQueue("Control DB drops", workers=workers, control_config=control_config) as queue:
        for n in range(8):
            queue.submit("square", f"n{n}", _square, n, on_done=done)
    assert sorted(results) == [n * n for n in range(8)]
    # Its row can't be finished; it reads as a dead run, and the lease
    # goes stale and is taken over by the next run.
    job = _rows(control_config, "SELECT state FROM work_jobs WHERE id = ?", (queue.job_id,))[0]
    assert job["state"] == "running"
    if workers > 1:
        _age_lease(control_config, workQueue.STALE_SECONDS + 5)
        with workQueue.WorkQueue("Next run", workers=2, control_config=control_config) as nxt:
            nxt.submit("square", "only", _square, 2)
        _ended_cleanly(control_config, nxt.job_id, "done")


def _age_lease(config, seconds):
    conn = _db.get_control_connection(config)
    try:
        with conn:
            conn.execute("UPDATE work_lease SET heartbeat_at = NOW(6) - INTERVAL ? SECOND", (seconds,))
            conn.execute("UPDATE work_jobs SET heartbeat_at = NOW(6) - INTERVAL ? SECOND"
                         " WHERE state IN ('waiting', 'running', 'paused')", (seconds,))
    finally:
        conn.close()


def _hold_lease(config):
    conn = _db.get_control_connection(config)
    try:
        with conn:
            conn.execute("INSERT INTO work_lease (id, holder, job_id, heartbeat_at)"
                         " VALUES (1, 'otherhost:1:abcd', 'other-job', NOW(6))"
                         " ON DUPLICATE KEY UPDATE holder = VALUES(holder), job_id = VALUES(job_id),"
                         " heartbeat_at = VALUES(heartbeat_at)")
    finally:
        conn.close()


def test_lease_expiry_uses_the_database_clock(control_config, monkeypatch):
    """A lease is stale by the database server's clock alone: this
    machine's clock running an hour fast (or slow) neither steals a live
    lease nor keeps a dead one."""
    monkeypatch.setattr(workQueue, "WAIT_POLL_SECONDS", 0.05)
    real_time, real_monotonic = time.time, time.monotonic
    store = workQueue._ControlStore(control_config)
    for skew in (3600, -3600):
        monkeypatch.setattr(time, "time", lambda skew=skew: real_time() + skew)
        _hold_lease(control_config)
        _age_lease(control_config, workQueue.STALE_SECONDS - 10)
        taken, current = store.take_lease("skewed-job", "here:1:skew")
        assert not taken and current["holder"] == "otherhost:1:abcd"
        _age_lease(control_config, workQueue.STALE_SECONDS + 5)
        taken, _current = store.take_lease("skewed-job", "here:1:skew")
        assert taken
    monkeypatch.setattr(time, "time", real_time)
    assert time.monotonic is real_monotonic


# ---------------------------------------------------------------------------
# TEST.21 Cancelling a run
# ---------------------------------------------------------------------------

def _raise_after_first(exc_factory):
    fired = []

    def done(result, seconds, weight):
        if not fired:
            fired.append(True)
            exc_factory()

    return done


def _sigterm_self():
    os.kill(os.getpid(), signal.SIGTERM)
    time.sleep(1)  # the handler raises before this ends


@pytest.mark.skipif(not hasattr(signal, "SIGTERM") or os.name == "nt", reason="POSIX signal")
@pytest.mark.parametrize("how, raised", [
    ("sigterm", SystemExit),
    ("ctrl-c", KeyboardInterrupt),
    ("exit", SystemExit),
])
def test_a_cancelled_parallel_run_ends_cancelled_and_frees_the_lease(control_config, how, raised):
    if threading.current_thread() is not threading.main_thread():
        pytest.skip("signals reach the main thread only")
    trigger = {
        "sigterm": _sigterm_self,
        "ctrl-c": lambda: (_ for _ in ()).throw(KeyboardInterrupt()),
        "exit": lambda: (_ for _ in ()).throw(SystemExit(1)),
    }[how]
    old = signal.getsignal(signal.SIGTERM)
    with pytest.raises(raised):
        with workQueue.WorkQueue(f"Cancelled by {how}", workers=2, control_config=control_config) as queue:
            for n in range(12):
                queue.submit("square", f"n{n}", _slow_square, n, on_done=_raise_after_first(trigger))
    assert signal.getsignal(signal.SIGTERM) is old  # handler put back
    _ended_cleanly(control_config, queue.job_id, "cancelled")
    states = {row["state"] for row in _rows(control_config, "SELECT state FROM work_tasks WHERE job_id = ?",
                                             (queue.job_id,))}
    assert states <= {"done", "cancelled"}


GALAXY_RING = 3
# Enough systems per sector that a loaded machine cannot finish the whole
# run between the first sector being saved and the test's signal landing
# (TEST.101: with 6 the run sometimes ended on its own, exit status 0).
NUM_SYSTEMS = 40


def test_a_stop_signal_swallowed_by_a_database_call_still_stops_the_run(control_config, monkeypatch):
    """TEST.73: a SIGTERM's `SystemExit` can land inside pymysql, which
    turns it into a database error that the queue's bookkeeping logs and
    carries on past; the run then finished as if nothing had happened.
    The queue remembers the signal and stops at the next step."""
    with pytest.raises(SystemExit) as raised:
        with workQueue.WorkQueue("Swallowed signal", workers=2, control_config=control_config) as queue:
            real = queue._store.add_tasks
            sent = []

            def add_tasks(job_id, tasks):
                if not sent:
                    sent.append(True)
                    try:
                        os.kill(os.getpid(), signal.SIGTERM)
                        time.sleep(0.5)
                    except BaseException:  # noqa: BLE001 -- what pymysql did with it
                        pass
                return real(job_id, tasks)

            monkeypatch.setattr(queue._store, "add_tasks", add_tasks)
            for n in range(8):
                queue.submit("square", f"n{n}", _slow_square, n)
    assert sent and raised.value.code == 128 + signal.SIGTERM
    _ended_cleanly(control_config, queue.job_id, "cancelled")


def _sector_counts(config):
    conn = _db.get_connection(config)
    try:
        sectors = conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"]
        short = conn.execute(
            "SELECT s.id, COUNT(ss.id) AS n FROM sectors s LEFT JOIN star_systems ss ON ss.sector_id = s.id"
            " GROUP BY s.id HAVING n <> ?", (NUM_SYSTEMS,)).fetchall()
        orphans = conn.execute(
            "SELECT COUNT(*) AS n FROM star_systems ss LEFT JOIN sectors s ON s.id = ss.sector_id"
            " WHERE ss.sector_id IS NOT NULL AND s.id IS NULL").fetchone()["n"]
        conn.rollback()
        return sectors, short, orphans
    finally:
        conn.close()


@pytest.mark.slow
@pytest.mark.skipif(os.name == "nt", reason="POSIX process groups")
@pytest.mark.parametrize("sig, group", [
    (signal.SIGINT, True),    # Ctrl+C in a terminal reaches the run and its workers
    (signal.SIGTERM, True),   # Cancel on the Generate page stops the step's whole process group
    (signal.SIGTERM, False),  # `kill <pid>` of the run alone
])
def test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector(control_config, sig, group):
    _plan_wide_galaxy(control_config)
    env = dict(os.environ, PLANETGEN_CONTROL_DATABASE=control_config.database, PLANETGEN_WORKERS="2")
    env.pop(workQueue.PARENT_ENV_VAR, None)
    run = subprocess.Popen(
        [sys.executable, "-m", "planetgen.cli.generate", "galaxy", "--ring", str(GALAXY_RING),
         "--num-systems", str(NUM_SYSTEMS), "--workers", "2"] + _mysql_argv(control_config),
        cwd=REPO, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, start_new_session=True)
    try:
        deadline = time.time() + 120
        while _sector_counts(control_config)[0] < 1:
            assert run.poll() is None, run.stderr.read().decode(errors="replace")[-2000:]
            assert time.time() < deadline, "no sector was saved"
            time.sleep(0.05)
        if group:
            os.killpg(run.pid, sig)
        else:
            os.kill(run.pid, sig)
        run.wait(timeout=120)
    finally:
        if run.poll() is None:
            os.killpg(run.pid, signal.SIGKILL)
            run.wait()
    assert run.returncode != 0

    sectors, short, orphans = _sector_counts(control_config)
    assert sectors >= 1 and short == [] and orphans == 0
    jobs = _rows(control_config, "SELECT kind, state FROM work_jobs WHERE parent_id IS NULL OR kind = 'queue'")
    assert jobs and all(job["state"] == "cancelled" for job in jobs), jobs
    assert _rows(control_config, "SELECT COUNT(*) AS n FROM work_tasks WHERE state IN ('queued', 'running')"
                 )[0]["n"] == 0
    assert _rows(control_config, "SELECT holder FROM work_lease")[0]["holder"] is None
