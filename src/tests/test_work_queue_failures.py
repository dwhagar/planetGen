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

import json
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


def _sleep_long(payload):
    time.sleep(600)


def _die_on_two(payload):
    if payload == 2:
        os._exit(7)  # the worker process dies, as an OOM kill would
    time.sleep(0.05)
    return payload


def _a_set(payload):
    return {payload, "x"}


def _a_generator(payload):
    return (n for n in range(payload))


_BOOKKEEPING_ERRORS = []
"""list: The bookkeeping the queue could not write in the running test (`WorkQueue._book` carries on past it and only
debug-logs it), so a missing row can say why (TEST.297)."""


@pytest.fixture(autouse=True)
def _bookkeeping_errors_are_kept(monkeypatch):
    _BOOKKEEPING_ERRORS.clear()
    real = workQueue.log.debug

    def keep(message, *args, **kwargs):
        if isinstance(message, str) and message.startswith("Work queue: could not"):
            _BOOKKEEPING_ERRORS.append(message)
        return real(message, *args, **kwargs)

    monkeypatch.setattr(workQueue.log, "debug", keep)


def _ended_cleanly(config, job_id, state):
    """The run's row ended in `state`, none of its tasks still claims to
    run, and the lease is free."""
    found = _rows(config, "SELECT state, finished_at FROM work_jobs WHERE id = ?", (job_id,))
    assert found, (f"no work_jobs row for job {job_id}; it has {_rows(config, 'SELECT id, state FROM work_jobs')}; the "
                   f"bookkeeping the queue could not write: {_BOOKKEEPING_ERRORS}")
    job = found[0]
    assert job["state"] == state and job["finished_at"] is not None
    running = _rows(config, "SELECT COUNT(*) AS n FROM work_tasks WHERE job_id = ?"
                            " AND state IN ('queued', 'running')", (job_id,))[0]["n"]
    assert running == 0
    lease = _rows(config, "SELECT holder FROM work_lease")
    assert not lease or lease[0]["holder"] is None


def _run_ends(finished, config, idle_seconds=180):
    """
    Waits for `finished` (the run's thread) and says whether it came, giving up only when the run has stood still: the
    run's record (how many tasks are in each state) has not changed for `idle_seconds`. A fixed limit on the whole run
    also fails one that is merely slow, as a loaded machine makes the start of a worker process (TEST.277); a hung one
    stops changing its record, and the pytest limit still ends a run that keeps changing it for ever. The window is 180 s:
    at a load of 6 to 8 a replacement worker took longer than 60 s to start and the run was taken for hung (TEST.297).
    """
    last, since = None, time.monotonic()
    while not finished.wait(1):
        seen = [(row["state"], row["n"]) for row in
                _rows(config, "SELECT state, COUNT(*) AS n FROM work_tasks GROUP BY state ORDER BY state")]
        if seen != last:
            last, since = seen, time.monotonic()
        elif time.monotonic() - since > idle_seconds:
            return False
    return True


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
    assert _run_ends(finished, control_config), "the run hung after a worker died: its record stood still for 180 s"
    assert isinstance(outcome.get("raised"), redisqueue.WorkerDied)


def test_workers_that_die_at_start_fail_the_run_instead_of_being_started_for_ever(monkeypatch, redis_server):
    """TEST.295: with every worker exiting at once, tasks waiting and the queue never empty, the executor started a
    worker every few milliseconds (752 in 21 s) and the run never ended."""
    started = []

    def dying(names, url, name=None, python=None, burst=True):
        started.append(name)
        return [sys.executable, "-c", "raise SystemExit(3)"]

    monkeypatch.setattr(redisqueue, "worker_argv", dying)
    began = time.monotonic()
    with pytest.raises(redisqueue.WorkerDied, match=r"10 workers in a row exited with a failure.*code 3"):
        with workQueue.WorkQueue("Workers die at start", workers=2) as queue:
            for n in range(4):
                queue.submit("die", f"n{n}", _square, n)
    assert time.monotonic() - began < 60
    assert redisqueue.WORKER_CRASH_LIMIT <= len(started) <= redisqueue.WORKER_CRASH_LIMIT + 2, started


def test_a_run_whose_record_stands_still_is_taken_for_hung(control_config):
    """TEST.277: the wait of the test above gives up on a run that changes nothing, and not before."""
    never = threading.Event()
    started = time.monotonic()
    assert not _run_ends(never, control_config, idle_seconds=2)
    assert 2 <= time.monotonic() - started < 10
    done = threading.Event()
    done.set()
    assert _run_ends(done, control_config, idle_seconds=2)


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
def test_a_result_that_isnt_json_still_completes(control_config, workers):
    results = []
    with workQueue.WorkQueue("Not JSON", workers=workers, control_config=control_config) as queue:
        queue.submit("set", "s", _a_set, 5, on_done=lambda result, *_: results.append(result))
    assert results == [{5, "x"}]
    row = _rows(control_config, "SELECT state FROM work_tasks WHERE job_id = ?", (queue.job_id,))[0]
    assert row["state"] == "done"
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


_STOPPED_RUN = """
import json, sys
from planetgen.db.store import MySQLConfig
from planetgen.queue import work
from tests.test_work_queue_failures import _sleep_long

with work.WorkQueue("Stopped run", workers=2, control_config=MySQLConfig(**json.loads(sys.argv[1]))) as queue:
    print(queue.job_id, flush=True)
    for n in range(4):
        queue.submit("sleep", f"n{n}", _sleep_long, n)
    queue.drain()
"""


def _descendants(pid):
    """The pids of every process below `pid` (workers and their work horses), from /proc."""
    parents = {}
    for entry in os.listdir("/proc"):
        if entry.isdigit():
            try:
                with open(f"/proc/{entry}/stat") as handle:
                    parents[int(entry)] = int(handle.read().rsplit(")", 1)[1].split()[1])
            except (OSError, ValueError, IndexError):
                continue
    found, queue = set(), [pid]
    while queue:
        parent = queue.pop()
        for child, ppid in parents.items():
            if ppid == parent and child not in found:
                found.add(child)
                queue.append(child)
    return found


def _alive(pid):
    """Whether `pid` still runs (a zombie is dead)."""
    try:
        with open(f"/proc/{pid}/stat") as handle:
            return handle.read().rsplit(")", 1)[1].split()[0] != "Z"
    except OSError:
        return False


def test_a_sigterm_stops_the_workers_busy_with_a_long_task_and_leaves_no_process(control_config):
    """OPS.59: the run's workers were left running their task (and the run waiting for them) after a SIGTERM."""
    kwargs = {"host": control_config.host, "port": control_config.port, "user": control_config.user,
              "password": control_config.password, "database": control_config.database}
    env = dict(os.environ, PLANETGEN_CONTROL_DATABASE=control_config.database)
    env.pop(workQueue.PARENT_ENV_VAR, None)
    run = subprocess.Popen([sys.executable, "-c", _STOPPED_RUN, json.dumps(kwargs)], cwd=os.path.join(REPO, "src"),
                           env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True,
                           preexec_fn=lambda: signal.signal(signal.SIGINT, signal.SIG_DFL))
    try:
        job_id = run.stdout.readline().decode().strip()
        deadline = time.time() + 120
        while len(_rows(control_config, "SELECT id FROM work_tasks WHERE job_id = ? AND state = 'running'",
                        (job_id,))) < 2:
            assert run.poll() is None, run.stderr.read().decode(errors="replace")[-2000:]
            assert time.time() < deadline, "no task started running"
            time.sleep(0.1)
        while len(_descendants(run.pid)) < 4:   # two workers, and a work horse each
            assert time.time() < deadline, "the workers did not start their tasks"
            time.sleep(0.1)
        below = _descendants(run.pid)
        os.kill(run.pid, signal.SIGTERM)
        assert run.wait(timeout=60) == 128 + signal.SIGTERM
        for _ in range(100):
            if not any(_alive(pid) for pid in below):
                break
            time.sleep(0.1)
        assert not [pid for pid in below if _alive(pid)], "a worker or work horse outlived the run"
    finally:
        if run.poll() is None:
            os.killpg(run.pid, signal.SIGKILL)
            run.wait()
        for pid in _descendants(run.pid):
            try:
                os.kill(pid, signal.SIGKILL)
            except OSError:
                pass
    _ended_cleanly(control_config, job_id, "cancelled")


GALAXY_SHELL = 5
# Enough systems per sector that a loaded machine cannot finish the whole
# run between the first sector being saved and the test's signal landing
# (TEST.101: with 6 the run sometimes ended on its own, exit status 0).
NUM_SYSTEMS = 16


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


_LAUNCH_GENERATE = ("from planetgen.generation import run_galaxy; run_galaxy.require_datasets = lambda: None; "
                    "from planetgen.cli import generate; generate.main()")
"""str: `python -m planetgen.cli.generate` with the data gate (GEN.590) passed: a test machine downloads no datasets."""


@pytest.mark.slow
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
        [sys.executable, "-c", _LAUNCH_GENERATE, "galaxy", "--shell", str(GALAXY_SHELL),
         "--num-systems", str(NUM_SYSTEMS), "--workers", "2"] + _mysql_argv(control_config),
        cwd=REPO, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, start_new_session=True,
        # A process inherits an ignored Ctrl+C, and the run then leaves it
        # alone ("someone else handles it") and finishes with status 0.
        # Another test in this pytest worker can leave SIGINT ignored
        # (TEST.101).
        preexec_fn=lambda: signal.signal(signal.SIGINT, signal.SIG_DFL))
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
        run.wait(timeout=300)
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
    assert all(row["holder"] is None for row in _rows(control_config, "SELECT holder FROM work_lease"))
