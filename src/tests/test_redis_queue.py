# tests/test_redis_queue.py

"""
The work queue's Redis side (PERF.24, step 1): `planetgen.queue.
redisqueue` and the `planetgen.cli.worker` burst workers it starts.

Needs a Redis server: `PLANETGEN_TEST_REDIS_URL` (CI sets it); skips
without one. Each test uses its own queue names, so parallel runs share
the server safely. The work queue's own tests
(`test_work_queue*.py`) run whole runs on it.
"""

import os
import pickle
import queue as queue_module
import secrets

import pytest

from planetgen.queue import redisqueue
from planetgen.util import settings

URL = os.environ.get("PLANETGEN_TEST_REDIS_URL")


@pytest.fixture
def connection():
    if not URL:
        pytest.skip("PLANETGEN_TEST_REDIS_URL is not set")
    connection = redisqueue.connect(URL)
    try:
        connection.ping()
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"no Redis at {URL}: {exc}")
    return connection


@pytest.fixture
def names(connection):
    names = [f"test-{secrets.token_hex(6)}-{n}" for n in range(2)]
    yield names
    for name in names:
        redisqueue.queue(name, connection).delete(delete_jobs=True)


def test_redis_url_prefers_the_environment(monkeypatch):
    monkeypatch.setenv("PLANETGEN_REDIS_URL", "redis://example:1/2")
    assert redisqueue.redis_url() == "redis://example:1/2"
    monkeypatch.delenv("PLANETGEN_REDIS_URL")
    assert redisqueue.redis_url() == settings.get_settings().redis.url


def test_worker_class_forks_where_it_can():
    import rq
    assert redisqueue.worker_class() is (rq.Worker if hasattr(os, "fork") else rq.SpawnWorker)


def test_worker_argv_names_the_worker():
    argv = redisqueue.worker_argv(["q"], "redis://h:1/0", "q-1")
    assert argv[1:] == ["-m", "planetgen.cli.worker", "--burst", "--url", "redis://h:1/0", "--name", "q-1", "q"]


def test_the_executor_runs_tasks_on_its_own_burst_workers(connection, names):
    import math

    executor = redisqueue.RQExecutor(2, math.hypot, names[0], url=URL)
    futures = [executor.submit(3, n) for n in (4, 4, 4)]
    done = []
    while len(done) < len(futures):
        done += executor.wait([f for f in futures if f not in done], block=True)
    assert [future.job.return_value() for future in futures] == [5.0, 5.0, 5.0]
    assert len(executor._procs) <= 2
    executor.shutdown(futures)
    assert all(proc.poll() is not None for proc in executor._procs.values()) or not executor._procs
    assert redisqueue.queue(names[0], connection).count == 0


def test_a_task_that_wont_pickle_fails_without_being_queued(connection, names):
    executor = redisqueue.RQExecutor(1, print, names[0], url=URL)
    future = executor.submit(lambda: None)
    assert future.done()
    with pytest.raises(Exception):
        future.result()
    executor.shutdown([future])


def test_no_redis_server_is_unavailable():
    with pytest.raises(redisqueue.Unavailable):
        redisqueue.RQExecutor(1, print, "nowhere", url="redis://127.0.0.1:1/0")


def test_a_channel_carries_items_through_a_pickle(connection):
    # The scatter's progress channel rides in a task's payload to a burst
    # worker, so it has to survive pickling and still reach the same list.
    channel = redisqueue.Channel(URL, f"test-{secrets.token_hex(6)}:progress")
    try:
        pickle.loads(pickle.dumps(channel)).put((3, 10, 20.5))
        assert channel.get(timeout=1) == (3, 10, 20.5)
        with pytest.raises(queue_module.Empty):
            channel.get(timeout=0.1)
    finally:
        channel.close()
