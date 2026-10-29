# tests/test_work_queue.py

"""
Tests for the worker count and task seeds (PERF.8, `planetgen/queue/work.py`)
and parallel sector generation (PERF.7): a `planetgen galaxy` run on several
workers producing the same links between neighboring sectors a one-at-a-time
run would.

Tests that take the `mysql_config` fixture (see `conftest.py`) are
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import math
import random
import sys

import pytest

from planetgen.util import detmath
from planetgen.cli import generate as generate_cli
from planetgen.galaxy import spherical
from planetgen.generation import run_common
from planetgen.queue import work as workQueue
from planetgen.db import store

from tests.test_galaxy_gen import _mysql_argv, _plan_wide_galaxy


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


@pytest.fixture
def control_config(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.close()
    return mysql_config


def _run_states(redis_url, database, prefix=None):
    """The states of the runs on Redis that were made on `database` (the run records, PERF.98), newest first. `prefix` is the
    key prefix of a run made by another process (a test's own process has its private one)."""
    from unittest import mock

    from planetgen.queue import redisqueue, runrecord
    conn = redisqueue.connect(redis_url)
    with mock.patch.object(runrecord, "PREFIX", prefix or runrecord.PREFIX):
        runs = [runrecord.read_run(conn, run_id) for run_id in runrecord.recent_runs(conn)]
    return [run["state"] for run in runs if run and run["database"] == database]


def _rows(config, sql, params=()):
    conn = store.get_control_connection(config)
    try:
        return conn.execute(sql, params).fetchall()
    finally:
        conn.close()


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
            (detmath.dist(point, p), n) for n, p in systems[:40] if detmath.dist(point, p) <= store.NEAREST_SYSTEMS_SEARCH_PC
        )[:store.NEAREST_SYSTEMS_COUNT]
        assert sparse.nearest(point) == expected
        expected_dense = sorted(
            (detmath.dist(point, p), n) for n, p in systems if detmath.dist(point, p) <= store.NEAREST_SYSTEMS_SEARCH_PC
        )[:store.NEAREST_SYSTEMS_COUNT]
        assert dense.nearest(point) == expected_dense


def _run_galaxy(argv):
    old_argv = sys.argv
    try:
        sys.argv = ["planetgen", "galaxy"] + argv
        generate_cli.main()
    finally:
        sys.argv = old_argv


def test_a_parallel_galaxy_run_links_neighbors_like_a_serial_one(control_config, monkeypatch, redis_server):
    """Two workers fill one shell at once: every sector is saved once, the
    run counts them all, and every stored nearest-system list matches a
    from-scratch recompute (neighbors saved at the same moment still see
    each other)."""
    monkeypatch.setenv("PLANETGEN_CONTROL_DATABASE", control_config.database)
    _plan_wide_galaxy(control_config)
    for counter in run_common.RUN_COUNTS:
        run_common.RUN_COUNTS[counter] = 0
    _run_galaxy(["--shell", "2", "--num-systems", "6", "--workers", "2"] + _mysql_argv(control_config))

    conn = store.get_connection(control_config)
    try:
        sectors = conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors").fetchall()
        expected = spherical.cell_count(2)
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

    states = _run_states(redis_server, control_config.database)        # the runs' own records on Redis (PERF.98)
    assert states and set(states) == {"done"}
