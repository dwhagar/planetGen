# tests/test_work_queue.py

"""
Tests for the generation work queue (PERF.8, `planetgen/queue/work.py`)
and parallel sector generation (PERF.7): worker counts, per-task seeds,
the pool, failures, the control database's lease and task rows, and a
`planetgen galaxy` run on several workers producing the same links
between neighboring sectors a one-at-a-time run would.

Tests that take the `mysql_config` fixture (see `conftest.py`) are
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import math
import random
import sys
import threading
import time

import pytest

from planetgen.util import draw
from planetgen.cli import generate as generate_cli
from planetgen.galaxy.geometry import ring_sector_count
from planetgen.generation import run_common
from planetgen.queue import work as workQueue
from planetgen.db import store

from tests.test_galaxy_gen import _mysql_argv, _plan_wide_galaxy


def _square(payload):
    return payload * payload


def _draw(payload):
    return draw.random()


def _fail_on_three(payload):
    if payload == 3:
        raise ValueError("task three broke")
    time.sleep(0.05)
    return payload


def _collector():
    results = []

    def done(result, seconds, weight):
        assert seconds >= 0
        results.append((result, weight))

    return results, done


# ---------------------------------------------------------------------------
# Worker counts, priority, seeds
# ---------------------------------------------------------------------------

def test_worker_count_is_80_percent_of_the_cores(monkeypatch):
    monkeypatch.delenv(workQueue.WORKERS_ENV_VAR, raising=False)
    monkeypatch.setattr(workQueue, "cpu_count", lambda: 10)
    assert workQueue.worker_count() == 8
    assert workQueue.worker_count(0, "db.example.org") == 8
    assert workQueue.worker_count(None, "127.0.0.1") == 7
    assert workQueue.worker_count("auto", "localhost") == 7
    monkeypatch.setattr(workQueue, "cpu_count", lambda: 1)
    assert workQueue.worker_count(None, "localhost") == 1


def test_an_explicit_worker_count_wins(monkeypatch):
    monkeypatch.setenv(workQueue.WORKERS_ENV_VAR, "3")
    assert workQueue.worker_count() == 3
    assert workQueue.worker_count(5, "localhost") == 5
    monkeypatch.setenv(workQueue.WORKERS_ENV_VAR, "auto")
    monkeypatch.setattr(workQueue, "cpu_count", lambda: 5)
    assert workQueue.worker_count() == 4


def test_task_seeds_depend_only_on_the_run_and_the_key():
    assert workQueue.task_seed(7, "1,2,3") == workQueue.task_seed(7, "1,2,3")
    assert workQueue.task_seed(7, "1,2,3") != workQueue.task_seed(7, "1,2,4")
    assert workQueue.task_seed(7, "1,2,3") != workQueue.task_seed(8, "1,2,3")


def test_workers_get_a_lower_priority(monkeypatch):
    calls = []
    monkeypatch.setattr(workQueue.os, "nice", lambda step: calls.append(step), raising=False)
    monkeypatch.setattr(workQueue.os, "name", "posix")
    workQueue.lower_priority()
    assert calls == [workQueue.WORKER_NICENESS]


# ---------------------------------------------------------------------------
# The pool
# ---------------------------------------------------------------------------

def test_one_worker_runs_each_task_here_in_order():
    results, done = _collector()
    with workQueue.WorkQueue("inline", workers=1) as queue:
        assert not queue.parallel
        for n in range(5):
            queue.submit("square", n, _square, n, weight=n, on_done=done)
            assert results[-1] == (n * n, n)
    assert queue.finished == 5


def test_several_workers_run_every_task():
    results, done = _collector()
    with workQueue.WorkQueue("pool", workers=2) as queue:
        for n in range(12):
            queue.submit("square", n, _square, n, on_done=done)
    assert sorted(result for result, _weight in results) == [n * n for n in range(12)]
    assert queue.submitted == queue.finished == 12


def test_submit_each_gives_the_same_results_in_few_tasks():
    """PERF.79: many small items go out in a few chunks, each item still finishing with its own weight."""
    def run(workers):
        results, done = _collector()
        with workQueue.WorkQueue("each", workers=workers) as queue:
            queue.submit_each("square", [(f"n {n}", n, n + 1, done) for n in range(40)], _square)
        return sorted(results), queue.submitted

    inline, inline_tasks = run(1)
    pooled, pooled_tasks = run(2)
    assert pooled == inline == sorted((n * n, n + 1) for n in range(40))
    assert inline_tasks == 40
    assert pooled_tasks == 2 * workQueue.CHUNKS_PER_WORKER


def test_task_results_do_not_depend_on_the_worker_count():
    def run(workers):
        drawn = {}
        with workQueue.WorkQueue("seeds", workers=workers, run_seed=1234) as queue:
            for n in range(6):
                queue.submit("draw", f"task-{n}", _draw, n,
                             on_done=lambda result, _s, _w, n=n: drawn.__setitem__(n, result))
        return drawn

    two, three = run(2), run(3)
    assert two == three
    assert len(set(two.values())) == 6


def test_a_failed_task_stops_the_run_with_its_own_error():
    results, done = _collector()
    with pytest.raises(ValueError, match="task three broke"):
        with workQueue.WorkQueue("failing", workers=2) as queue:
            for n in range(40):
                queue.submit("maybe", n, _fail_on_three, n, on_done=done)
    finished = [result for result, _weight in results]
    assert 3 not in finished
    # Tasks queued after the failure was seen were never started.
    assert len(finished) < 39


# ---------------------------------------------------------------------------
# The control database's lease and rows
# ---------------------------------------------------------------------------

@pytest.fixture
def control_config(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.close()
    return mysql_config


def _rows(config, sql, params=()):
    conn = store.get_control_connection(config)
    try:
        return conn.execute(sql, params).fetchall()
    finally:
        conn.close()


def test_a_run_records_its_job_and_tasks_and_frees_the_lease(control_config):
    with workQueue.WorkQueue("Recorded run", workers=2, control_config=control_config) as queue:
        for n in range(5):
            queue.submit("square", f"n{n}", _square, n)
    job = _rows(control_config, "SELECT * FROM work_jobs")[0]
    assert job["id"] == queue.job_id and job["title"] == "Recorded run"
    assert (job["state"], job["workers"], job["tasks_queued"], job["tasks_done"]) == ("done", 2, 5, 5)
    tasks = _rows(control_config, "SELECT task_key, state, seconds FROM work_tasks ORDER BY id")
    assert [t["task_key"] for t in tasks] == [f"n{n}" for n in range(5)]
    assert {t["state"] for t in tasks} == {"done"}
    assert all(t["seconds"] is not None for t in tasks)
    assert _rows(control_config, "SELECT holder FROM work_lease")[0]["holder"] is None


def test_a_failed_run_is_recorded(control_config):
    with pytest.raises(ValueError):
        with workQueue.WorkQueue("Failing run", workers=2, control_config=control_config) as queue:
            for n in range(6):
                queue.submit("maybe", f"n{n}", _fail_on_three, n)
    job = _rows(control_config, "SELECT state, tasks_failed FROM work_jobs WHERE id = ?", (queue.job_id,))[0]
    assert (job["state"], job["tasks_failed"]) == ("failed", 1)
    error = _rows(control_config, "SELECT error FROM work_tasks WHERE state = 'failed'")[0]["error"]
    assert "task three broke" in error


def test_a_second_run_waits_for_the_lease(control_config, monkeypatch):
    monkeypatch.setattr(workQueue, "WAIT_POLL_SECONDS", 0.1)
    conn = store.get_control_connection(control_config)
    try:
        with conn:
            conn.execute("INSERT INTO work_lease (id, holder, job_id, heartbeat_at)"
                         " VALUES (1, 'otherhost:1:abcd', 'other-job', NOW(6))")
    finally:
        conn.close()

    def release():
        time.sleep(1.0)
        other = store.get_control_connection(control_config)
        try:
            with other:
                other.execute("UPDATE work_lease SET holder = NULL WHERE id = 1")
        finally:
            other.close()

    waits = []
    releaser = threading.Thread(target=release)
    releaser.start()
    started = time.monotonic()
    with workQueue.WorkQueue("Second run", workers=2, control_config=control_config,
                             on_wait=waits.append) as queue:
        queue.submit("square", "only", _square, 4)
    releaser.join()
    assert time.monotonic() - started >= 0.9
    assert waits and waits[0]["holder"] == "otherhost:1:abcd"


def test_a_dead_runs_lease_is_taken_over(control_config):
    conn = store.get_control_connection(control_config)
    try:
        with conn:
            conn.execute("INSERT INTO work_jobs (id, title, holder, state, workers, created_at, heartbeat_at)"
                         " VALUES ('dead-job', 'Died', 'gone:1:dead', 'running', 2,"
                         " NOW(6) - INTERVAL 5 MINUTE, NOW(6) - INTERVAL 2 MINUTE)")
            conn.execute("INSERT INTO work_tasks (job_id, kind, task_key, state, created_at)"
                         " VALUES ('dead-job', 'sector', '1,0,0', 'running', NOW(6))")
            conn.execute("INSERT INTO work_lease (id, holder, job_id, heartbeat_at)"
                         " VALUES (1, 'gone:1:dead', 'dead-job', NOW(6) - INTERVAL 2 MINUTE)")
    finally:
        conn.close()
    with workQueue.WorkQueue("Next run", workers=2, control_config=control_config) as queue:
        queue.submit("square", "only", _square, 2)
    assert _rows(control_config, "SELECT state FROM work_jobs WHERE id = 'dead-job'")[0]["state"] == "cancelled"
    assert _rows(control_config, "SELECT state FROM work_tasks WHERE job_id = 'dead-job'")[0]["state"] == "cancelled"


def test_without_the_control_tables_the_pool_still_runs(mysql_config):
    results, done = _collector()
    with workQueue.WorkQueue("No control tables", workers=2, control_config=mysql_config) as queue:
        for n in range(3):
            queue.submit("square", n, _square, n, on_done=done)
    assert sorted(result for result, _weight in results) == [0, 1, 4]


# ---------------------------------------------------------------------------
# Parallel generation
# ---------------------------------------------------------------------------

def test_nearest_search_without_shells_matches_the_shell_walk():
    rng = random.Random(5)
    systems = [(n, (rng.uniform(-6, 6), rng.uniform(-6, 6), rng.uniform(-6, 6))) for n in range(3000)]
    dense = store._SystemGrid(systems)
    sparse = store._SystemGrid(systems[:40])
    for _ in range(50):
        point = (rng.uniform(-8, 8), rng.uniform(-8, 8), rng.uniform(-8, 8))
        assert len(dense.systems) >= (2 * (math.ceil(store.NEAREST_SYSTEMS_SEARCH_PC) + 1) + 1) ** 3
        expected = sorted(
            (math.dist(point, p), n) for n, p in systems[:40] if math.dist(point, p) <= store.NEAREST_SYSTEMS_SEARCH_PC
        )[:store.NEAREST_SYSTEMS_COUNT]
        assert sparse.nearest(point) == expected
        expected_dense = sorted(
            (math.dist(point, p), n) for n, p in systems if math.dist(point, p) <= store.NEAREST_SYSTEMS_SEARCH_PC
        )[:store.NEAREST_SYSTEMS_COUNT]
        assert dense.nearest(point) == expected_dense


def _run_galaxy(argv):
    old_argv = sys.argv
    try:
        sys.argv = ["planetgen", "galaxy"] + argv
        generate_cli.main()
    finally:
        sys.argv = old_argv


def test_a_parallel_galaxy_run_links_neighbors_like_a_serial_one(control_config, monkeypatch):
    """Two workers fill one ring at once: every sector is saved once, the
    run counts them all, and every stored nearest-system list matches a
    from-scratch recompute (neighbors saved at the same moment still see
    each other)."""
    monkeypatch.setenv("PLANETGEN_CONTROL_DATABASE", control_config.database)
    _plan_wide_galaxy(control_config)
    for counter in run_common.RUN_COUNTS:
        run_common.RUN_COUNTS[counter] = 0
    _run_galaxy(["--ring", "1", "--num-systems", "6", "--workers", "2"] + _mysql_argv(control_config))

    conn = store.get_connection(control_config)
    try:
        sectors = conn.execute("SELECT ring_index, layer_index, ring_slot_index FROM sectors").fetchall()
        expected = ring_sector_count(1)
        assert len(sectors) == len({tuple(row.values()) for row in sectors}) == expected
        assert run_common.RUN_COUNTS["sectors"] == expected
        assert run_common.RUN_COUNTS["systems"] == conn.execute(
            "SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"]

        query = "SELECT object_table, object_id, neighbor_rank, neighbor_system_id FROM nearest_systems"
        stored = {tuple(row.values()) for row in conn.execute(query).fetchall()}
        ids = [row["id"] for row in conn.execute("SELECT id FROM sectors").fetchall()]
        store.refresh_nearest_systems(conn, ids)
        recomputed = {tuple(row.values()) for row in conn.execute(query).fetchall()}
        conn.rollback()
        assert stored and stored == recomputed
    finally:
        conn.close()

    job = _rows(control_config, "SELECT state, workers, tasks_done FROM work_jobs WHERE kind = 'queue'")[0]
    assert (job["state"], job["workers"], job["tasks_done"]) == ("done", 2, expected)
