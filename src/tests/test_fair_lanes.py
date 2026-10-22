"""PERF.151: the interactive worker is reserved outside the bulk count and runs at a higher priority, and an API job
expected to run long goes to a lane of its own, so it can never hold the worker kept for short jobs."""

import math
import os
import time
import uuid

import pytest

from planetgen.queue import api_jobs, redisqueue, work, workers
from planetgen.util import settings

pytestmark = pytest.mark.real_queue


@pytest.fixture
def lanes(redis_server, monkeypatch):
    """The two API lanes under names of this test's own, and the worker starts recorded instead of made."""
    token = uuid.uuid4().hex[:12]
    monkeypatch.setattr(api_jobs, "QUEUE_NAME", f"planetgen-api-test-short-{token}")
    monkeypatch.setattr(api_jobs, "LONG_QUEUE_NAME", f"planetgen-api-test-long-{token}")
    started = []
    monkeypatch.setattr(workers.redisqueue, "start_detached", lambda argv, cwd: started.append(argv))
    connection = redisqueue.connect()
    connection.delete(workers.WANTED_KEY)
    return connection, started


def _waiting(connection, name):
    return redisqueue.queue(name, connection).count


def test_a_job_expected_to_run_over_the_limit_goes_to_the_long_lane(lanes):
    connection, started = lanes
    api_jobs.submit(math.pow, 2, 3)                                          # no estimate: short
    api_jobs.submit(math.pow, 2, 4, expected_seconds=60)                     # at the limit, not over: short
    api_jobs.submit(math.pow, 2, 5, expected_seconds=60.5)                   # over it: long
    api_jobs.submit(math.pow, 2, 6, expected_seconds=3 * 86400)
    assert _waiting(connection, api_jobs.QUEUE_NAME) == 2
    assert _waiting(connection, api_jobs.LONG_QUEUE_NAME) == 2
    assert not api_jobs.is_long(None) and not api_jobs.is_long(float("nan")) and api_jobs.is_long(1e9)


def test_the_long_job_limit_is_a_setting(lanes, monkeypatch):
    section = settings.get_settings().queue
    monkeypatch.setattr(section, "long_job_seconds", 5.0)
    assert api_jobs.is_long(6) and not api_jobs.is_long(5)
    assert api_jobs.lane(6) == (api_jobs.LONG_QUEUE_NAME, section.long_workers)
    assert api_jobs.lane(5) == (api_jobs.QUEUE_NAME, section.interactive_workers)


def test_each_lane_is_capped_by_its_own_setting_and_not_by_the_bulk_count(lanes, monkeypatch):
    connection, started = lanes
    section = settings.get_settings().queue
    monkeypatch.setattr(section, "interactive_workers", 2)
    monkeypatch.setattr(section, "long_workers", 1)
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 1)       # the bulk count: 1
    for n in range(5):
        api_jobs.submit(math.pow, 2, n)                                        # five short jobs
    for n in range(5):
        api_jobs.submit(math.pow, 3, n, expected_seconds=600)                  # five long ones
    short = [argv for argv in started if api_jobs.QUEUE_NAME in argv]
    long = [argv for argv in started if api_jobs.LONG_QUEUE_NAME in argv]
    assert len(short) == 2 and len(long) == 1


def test_a_long_job_running_leaves_a_worker_for_a_short_one_even_with_one_wanted(lanes, monkeypatch):
    connection, started = lanes
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 1)
    long_queue = redisqueue.queue(api_jobs.LONG_QUEUE_NAME, connection)
    import rq
    busy = rq.Worker([long_queue], connection=connection)                      # the long lane's one worker, alive
    busy.register_birth()
    connection.hset(busy.key, "pid", os.getpid())
    try:
        api_jobs.submit(math.pow, 2, 3, expected_seconds=600)                  # waits: the long lane is full
        assert not [argv for argv in started if api_jobs.LONG_QUEUE_NAME in argv]
        api_jobs.submit(math.pow, 2, 4)                                        # short: a worker of its own starts
        assert len([argv for argv in started if api_jobs.QUEUE_NAME in argv]) == 1
    finally:
        busy.register_death()


def test_the_interactive_worker_runs_at_a_higher_priority_than_the_bulk_and_long_ones():
    assert workers.niceness([redisqueue.INTERACTIVE_QUEUE]) == workers.INTERACTIVE_NICENESS < work.WORKER_NICENESS
    assert workers.niceness([redisqueue.BULK_QUEUE]) == work.WORKER_NICENESS
    assert workers.niceness([redisqueue.LONG_QUEUE]) == work.WORKER_NICENESS


def test_the_worker_command_lowers_its_priority_by_its_lane(redis_server, monkeypatch):
    from planetgen.cli import worker
    seen = []
    monkeypatch.setattr(worker.work, "lower_priority", lambda niceness=None: seen.append(niceness))
    monkeypatch.setattr(worker, "warm", lambda queues=None: None)
    for queue in (redisqueue.INTERACTIVE_QUEUE, redisqueue.LONG_QUEUE):
        assert worker.main(["--burst", "--url", redis_server, queue]) in (0, None)
    assert seen == [workers.INTERACTIVE_NICENESS, work.WORKER_NICENESS]


def test_the_lane_queues_carry_json_and_the_update_guard_counts_them():
    assert redisqueue.serializer_for(redisqueue.LONG_QUEUE) is redisqueue.JsonSerializer
    assert redisqueue.LONG_QUEUE in redisqueue.QUEUES
    assert workers.lane_cap(redisqueue.BULK_QUEUE) is None


def test_with_one_bulk_worker_wanted_a_long_job_does_not_hold_up_a_short_one(redis_server, monkeypatch):
    """End to end, with real workers: a long job running on the long lane, a short job answered in seconds."""
    token = uuid.uuid4().hex[:12]
    monkeypatch.setattr(api_jobs, "QUEUE_NAME", f"planetgen-api-test-short-{token}")
    monkeypatch.setattr(api_jobs, "LONG_QUEUE_NAME", f"planetgen-api-test-long-{token}")
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 1)
    connection = redisqueue.connect()
    connection.delete(workers.WANTED_KEY)
    long_job = None
    try:
        long_job = api_jobs.submit(time.sleep, 25, expected_seconds=600)
        deadline = time.monotonic() + 60
        while api_jobs.status(long_job)["state"] != "running" and time.monotonic() < deadline:
            time.sleep(0.2)
        assert api_jobs.status(long_job)["state"] == "running"
        started = time.monotonic()
        short_job = api_jobs.submit(math.pow, 2, 5)
        job = api_jobs.wait(short_job, 60)
        assert job["state"] == "succeeded" and job["result"] == 32.0
        assert time.monotonic() - started < 20
        assert api_jobs.status(long_job)["state"] == "running"                  # the long job is still going
    finally:
        for name in (api_jobs.QUEUE_NAME, api_jobs.LONG_QUEUE_NAME):
            for worker in workers.live_workers(redisqueue.queue(name, connection)):
                workers.command.send_shutdown_command(connection, worker.name)
