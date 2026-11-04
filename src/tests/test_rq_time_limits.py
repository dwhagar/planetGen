# tests/test_rq_time_limits.py

"""
Every job queued to RQ has a finite time limit taken from the best estimate of how long it takes (PERF.144): RQ then
expires the entry of a job whose worker was killed, and ends a job that hangs. The sites that queue:

* a run's units (`units.enqueue_units`): the stage's expected seconds;
* the Generate page's job (`jobs._spawn`): the sum of its steps' recorded times;
* an API job (`api_jobs.submit`): the route's estimate, else a fixed limit;
* the node prefetch and rebuild (`tileprefetch._enqueue_chunks`): the size of the chunk.

`redisqueue.TimedQueue` refuses a job without a limit, and no module makes an RQ queue or job of its own.
"""

import ast
import json
import math
import os
import pathlib
import secrets
import signal
import sys
import time

import pytest
import rq

from planetgen.queue import api_jobs, redisqueue, units, workers
from planetgen.web import job_runner, jobs
from planetgen.web.lib import tileprefetch

pytestmark = pytest.mark.real_queue

SOURCE = pathlib.Path(__file__).resolve().parents[1] / "planetgen"


@pytest.fixture
def limits(monkeypatch):
    """Every `(queue name, function, time limit)` a queue created through `redisqueue.queue` is given."""
    seen = []
    real = redisqueue.TimedQueue.create_job

    def create_job(self, func, *args, timeout=None, **kwargs):
        job = real(self, func, *args, timeout=timeout, **kwargs)
        seen.append((self.name, func, job.timeout))
        return job

    monkeypatch.setattr(redisqueue.TimedQueue, "create_job", create_job)
    monkeypatch.setattr(workers, "ensure_workers", lambda *a, **k: None)
    return seen


@pytest.fixture
def own_queue(redis_server):
    """A queue of this test's own and its connection; deleted afterwards."""
    connection = redisqueue.connect(redis_server)
    name = f"test-timeout-{secrets.token_hex(6)}"
    yield redisqueue.queue(name, connection)
    redisqueue.queue(name, connection).delete(delete_jobs=True)


# --- The limit ------------------------------------------------------------------------------------------

def test_a_limit_is_a_multiple_of_the_estimate_between_a_floor_and_a_ceiling():
    assert redisqueue.job_timeout(1000) == 10 * 1000
    assert redisqueue.job_timeout(1) == redisqueue.MIN_TIMEOUT_SECONDS
    assert redisqueue.job_timeout(10 ** 9) == redisqueue.MAX_TIMEOUT_SECONDS
    assert redisqueue.job_timeout(100, ceiling=500, floor=200) == 500
    assert redisqueue.job_timeout(1, ceiling=500, floor=200) == 200


@pytest.mark.parametrize("nothing", [None, 0, -5, float("nan"), float("inf"), -float("inf")])
def test_a_job_with_no_estimate_gets_the_fixed_ceiling_never_no_limit(nothing):
    assert redisqueue.job_timeout(nothing) == redisqueue.MAX_TIMEOUT_SECONDS
    assert redisqueue.job_timeout(nothing, ceiling=900) == 900
    assert redisqueue.job_timeout(nothing, ceiling=900, unestimated=60) == 60
    assert redisqueue.job_timeout(nothing) > 0


@pytest.mark.parametrize("missing", [{}, {"timeout": None}, {"timeout": -1}, {"timeout": 0}, {"timeout": True},
                                     {"timeout": "soon"}])
def test_a_queue_refuses_a_job_without_a_finite_limit(own_queue, missing):
    with pytest.raises(ValueError, match="time limit"):
        own_queue.enqueue_call(math.sqrt, args=(4,), **missing)
    with pytest.raises(ValueError, match="time limit"):
        own_queue.enqueue(math.sqrt, 4)
    with pytest.raises(ValueError, match="time limit"):
        own_queue.enqueue(math.sqrt, 4, job_timeout=missing.get("timeout"))
    with pytest.raises(ValueError, match="time limit"):
        own_queue.enqueue_many([own_queue.prepare_data(math.sqrt, args=(4,), **missing)])
    assert own_queue.count == 0


def test_a_queue_takes_a_job_with_a_limit_in_every_way_it_queues(own_queue):
    own_queue.enqueue(math.sqrt, 4, job_timeout=30)
    own_queue.enqueue_call(math.sqrt, args=(4,), timeout=31)
    own_queue.enqueue_many([own_queue.prepare_data(math.sqrt, args=(4,), timeout=32)])
    assert sorted(job.timeout for job in own_queue.jobs) == [30, 31, 32]


def test_the_queue_every_site_asks_for_is_the_timed_one(redis_server):
    assert isinstance(redisqueue.queue(redisqueue.BULK_QUEUE, redisqueue.connect(redis_server)), redisqueue.TimedQueue)


def _calls(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            yield node


def test_no_module_makes_an_rq_queue_or_job_of_its_own():
    """A new kind of job goes through `redisqueue.queue`, whose queue refuses a job without a limit."""
    found = []
    for path in sorted(SOURCE.rglob("*.py")):
        if path.name == "redisqueue.py":
            continue
        for call in _calls(path):
            function = call.func
            name = function.attr if isinstance(function, ast.Attribute) else getattr(function, "id", "")
            owner = ast.unparse(function.value) if isinstance(function, ast.Attribute) else ""
            if (name in ("Queue", "enqueue_job", "enqueue_at", "enqueue_in", "schedule_job") and owner.endswith("rq")) \
                    or (name in ("create", "enqueue_job") and "Job" in owner) or name == "enqueue_job":
                found.append(f"{path.relative_to(SOURCE)}:{call.lineno} {ast.unparse(function)}")
    assert found == []


def test_every_call_that_queues_a_job_names_its_limit():
    """The grep-style guard next to `TimedQueue`: a call that queues without a limit fails here even if no test runs it."""
    found = []
    for path in sorted(SOURCE.rglob("*.py")):
        for call in _calls(path):
            if isinstance(call.func, ast.Attribute) and call.func.attr in ("enqueue", "enqueue_call", "prepare_data") \
                    and not any(keyword.arg in ("job_timeout", "timeout") for keyword in call.keywords):
                found.append(f"{path.relative_to(SOURCE)}:{call.lineno} {ast.unparse(call.func)}")
    assert found == []


def test_no_site_passes_a_negative_or_missing_limit():
    """The grep-style guard behind `TimedQueue`: a literal -1 or None for a job's limit is a bug at the call."""
    found = []
    for path in sorted(SOURCE.rglob("*.py")):
        for call in _calls(path):
            for keyword in call.keywords:
                if keyword.arg in ("job_timeout", "timeout") and isinstance(keyword.value, ast.Constant) \
                        and isinstance(call.func, ast.Attribute) and call.func.attr in (
                            "enqueue", "enqueue_call", "enqueue_many", "prepare_data") \
                        and (keyword.value.value is None or keyword.value.value <= 0):
                    found.append(f"{path.relative_to(SOURCE)}:{call.lineno}")
    assert found == []


# --- A job over its limit ---------------------------------------------------------------------------------

def test_a_job_over_its_limit_is_ended_and_reported_failed(own_queue):
    job = own_queue.enqueue(time.sleep, 30, job_timeout=1, result_ttl=0, failure_ttl=60)
    rq.SimpleWorker([own_queue], connection=own_queue.connection).work(burst=True)
    job.refresh()
    assert job.get_status() == "failed" and "maximum timeout value (1 seconds)" in job.latest_result().exc_string


# --- Each site ----------------------------------------------------------------------------------------------

def test_the_units_of_a_stage_get_a_multiple_of_its_expected_time_or_the_ceiling(redis_server, limits):
    connection = redisqueue.connect(redis_server)
    run_id = f"timed-{secrets.token_hex(4)}"
    try:
        units.enqueue_units(connection, run_id, 0, "scatter-mass", [({}, 1.0)], expected_seconds=2000)
        units.enqueue_units(connection, run_id, 1, "scatter-mass", [({}, 1.0)])
    finally:
        redisqueue.queue(redisqueue.BULK_QUEUE, connection).empty()
    assert [timeout for _name, _func, timeout in limits] == [20_000, redisqueue.MAX_TIMEOUT_SECONDS]


@pytest.fixture
def jobs_root(tmp_path, monkeypatch):
    root = tmp_path / "jobs"
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(root))
    monkeypatch.setattr(jobs, "_start_detached", lambda *args: None)
    return str(root)


@pytest.mark.parametrize("expected, limit", [
    (None, jobs.JOB_CEILING_SECONDS),
    (60.0, jobs.JOB_FLOOR_SECONDS),
    (2 * 3600.0, 20 * 3600),
    (10 ** 9, jobs.JOB_CEILING_SECONDS)])
def test_the_generate_pages_job_gets_a_limit_from_the_time_its_steps_are_recorded_to_take(
        jobs_root, redis_server, limits, monkeypatch, expected, limit):
    monkeypatch.setattr(job_runner, "expected_seconds", lambda job: expected)
    job_id = jobs.start_job("reset", "Timed", [{"label": "x", "argv": [sys.executable, "-c", "pass"]}], root=jobs_root)
    try:
        assert [(name, timeout) for name, _func, timeout in limits] == [(jobs.QUEUE_PREFIX + job_id, limit)]
    finally:
        jobs.release_lock(job_id)       # the job lock is one Redis key: a job left holding it makes the next case's start raise JobBusy
        redisqueue.queue(jobs.QUEUE_PREFIX + job_id, redisqueue.connect(redis_server)).delete(delete_jobs=True)


def test_a_jobs_expected_seconds_are_the_sum_of_what_its_steps_took_before(monkeypatch, tmp_path):
    steps = [{"label": "Reset", "argv": ["reset"]}, {"label": "Other", "argv": ["other"]}]
    control = object()
    both = {"Reset": 30.0, "Other": 90.0}
    assert job_runner.step_estimates(control, steps, both)[0] == [30.0, 90.0]
    assert job_runner.step_estimates(control, steps, {"Reset": 30.0})[0] == [30.0, None]

    class Connection:
        def close(self):
            pass

    from planetgen.db import store
    from planetgen.web import jobs
    monkeypatch.setattr(store, "get_control_connection", lambda config: Connection())
    monkeypatch.setattr(jobs, "configured_jobs_dir", lambda: str(tmp_path))
    job = {"steps": steps, "env": {}}
    assert job_runner.expected_seconds(job) is None             # nothing recorded yet
    for label, seconds in both.items():
        job_runner._record_step_seconds(str(tmp_path), label, seconds)
    assert job_runner.expected_seconds(job) == 120.0
    job_runner._record_step_seconds(str(tmp_path), "Reset", 50.0)
    assert job_runner.recorded_step_seconds(str(tmp_path)) == {"Reset": 40.0, "Other": 90.0}
    monkeypatch.setattr(store, "get_control_connection", lambda config: (_ for _ in ()).throw(OSError("no database")))
    assert job_runner.expected_seconds(job) is None


def test_a_job_the_queue_ends_for_its_limit_stops_its_step_and_says_why(tmp_path):
    job_dir = tmp_path / "20261021-120000-abcd"
    job_dir.mkdir()
    marker = tmp_path / "step.pid"
    step = f"import os, time; open({str(marker)!r}, 'w').write(str(os.getpid())); time.sleep(60)"
    (job_dir / "job.json").write_text(json.dumps({
        "id": job_dir.name, "kind": "reset", "title": "Slow", "admin": None, "database": None,
        "created_at": time.time(), "cwd": str(tmp_path), "env": {},
        "steps": [{"label": "Sleep", "argv": [sys.executable, "-c", step]}]}), encoding="utf-8")
    old = signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGINT)
    try:
        with rq.timeouts.UnixSignalDeathPenalty(2, rq.timeouts.JobTimeoutException):
            assert job_runner.run(str(job_dir)) == 1
    finally:
        signal.signal(signal.SIGTERM, old[0])
        signal.signal(signal.SIGINT, old[1])
    state = json.loads((job_dir / "state.json").read_text(encoding="utf-8"))
    assert state["status"] == "failed" and "past its time limit" in state["error"]
    pid = int(marker.read_text(encoding="utf-8"))
    for _ in range(100):
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            break
        time.sleep(0.1)
    else:
        os.kill(pid, signal.SIGKILL)
        pytest.fail("the step kept running after its job was ended")


def test_an_api_job_gets_the_routes_estimate_else_a_fixed_limit(redis_server, limits):
    connection = redisqueue.connect(redis_server)
    try:
        api_jobs.submit(math.pow, 2, 3)
        api_jobs.submit(math.pow, 2, 3, expected_seconds=10)
        api_jobs.submit(math.pow, 2, 3, expected_seconds=10 ** 9)
    finally:
        redisqueue.queue(api_jobs.QUEUE_NAME, connection).empty()
    assert [timeout for _name, _func, timeout in limits] == [
        api_jobs.UNESTIMATED_SECONDS, redisqueue.MIN_TIMEOUT_SECONDS, api_jobs.ESTIMATED_CEILING_SECONDS]


def test_an_api_jobs_estimate_is_read_from_the_recorded_task_seconds(monkeypatch):
    from planetgen.generation import stats
    monkeypatch.setenv(stats.STATS_ENV_VAR, "0")
    assert api_jobs.recorded_task_seconds("sector") is None
    monkeypatch.setenv(stats.STATS_ENV_VAR, "1")
    monkeypatch.setattr(stats, "GenerationStats", lambda config: type("S", (), {
        "seconds_per_task": lambda self, kind, workers: 12.5})())
    assert api_jobs.recorded_task_seconds("sector") == 12.5

    def broken(config):
        raise OSError("no control database")

    monkeypatch.setattr(stats, "GenerationStats", broken)
    assert api_jobs.recorded_task_seconds("sector") is None


def test_a_node_chunk_gets_a_limit_from_its_size(redis_server, limits):
    assert tileprefetch.chunk_timeout(1) == tileprefetch.MIN_JOB_SECONDS
    assert tileprefetch.chunk_timeout(10 ** 6) == redisqueue.MAX_TIMEOUT_SECONDS
    connection = redisqueue.connect(redis_server)
    queue = redisqueue.queue(f"test-timeout-{secrets.token_hex(6)}", connection)
    try:
        tileprefetch._enqueue_chunks(queue, "db", [f"t{n}" for n in range(tileprefetch.CHUNK + 1)], "prefetch")
    finally:
        queue.delete(delete_jobs=True)
    assert [timeout for _name, _func, timeout in limits] == [tileprefetch.chunk_timeout(tileprefetch.CHUNK),
                                                             tileprefetch.chunk_timeout(1)]
