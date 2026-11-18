# tests/test_job_time_limits.py

"""
Every job has a finite time limit taken from the best estimate of how long it takes (PERF.144). OPS.235: RQ is gone, so what is left is the limit itself (`redisqueue.job_timeout`) and the
sites that record it for the supervisor to end a job by: an API job (`api_jobs.submit`: the route's estimate, else a fixed limit) and the node prefetch (`tileprefetch._add_chunks`: the size of the chunk).
"""

import math

import pytest

from planetgen.queue import api_jobs, redisqueue
from planetgen.web.lib import tileprefetch

pytestmark = pytest.mark.real_queue


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


def test_an_api_job_gets_the_routes_estimate_else_a_fixed_limit(redis_server, monkeypatch):
    """OPS.230: the supervisor ends a light job after the limit its record holds (`lightjobs.time_limit`), the one RQ's `job_timeout` was."""
    from planetgen.queue import lightjobs
    monkeypatch.setattr(lightjobs, "ask_supervisor", lambda connection=None, url=None: None)
    connection = redisqueue.connect(redis_server)
    ids = [api_jobs.submit(math.pow, 2, 3), api_jobs.submit(math.pow, 2, 3, expected_seconds=10), api_jobs.submit(math.pow, 2, 3, expected_seconds=10 ** 9)]
    try:
        limits = [int(connection.hget(lightjobs.job_key(job_id), "limit")) for job_id in ids]
    finally:
        for job_id in ids:
            lightjobs.cancel(connection, job_id)
    assert limits == [lightjobs.UNESTIMATED_SECONDS, redisqueue.MIN_TIMEOUT_SECONDS, lightjobs.ESTIMATED_CEILING_SECONDS]


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


def test_a_node_chunk_gets_a_limit_from_its_size(redis_server):
    """OPS.230: the node prefetch's jobs are light jobs; the limit the supervisor ends one by is on its record."""
    from planetgen.queue import lightjobs
    assert tileprefetch.chunk_timeout(1) == tileprefetch.MIN_JOB_SECONDS
    assert tileprefetch.chunk_timeout(10 ** 6) == redisqueue.MAX_TIMEOUT_SECONDS
    connection = redisqueue.connect(redis_server)
    before = set(connection.zrange(lightjobs.pending_key(), 0, -1))
    try:
        assert tileprefetch._add_chunks(connection, "db", "gen", [f"t{n}" for n in range(tileprefetch.CHUNK + 1)]) == tileprefetch.CHUNK + 1
        ids = [item for item in connection.zrange(lightjobs.pending_key(), 0, -1) if item not in before]
        found = sorted((int(connection.hget(lightjobs.job_key(job_id.decode() if isinstance(job_id, bytes) else job_id), "limit")) for job_id in ids))
    finally:
        for job_id in ids:
            lightjobs.cancel(connection, job_id.decode() if isinstance(job_id, bytes) else job_id)
    assert found == sorted([tileprefetch.chunk_timeout(tileprefetch.CHUNK), tileprefetch.chunk_timeout(1)])
