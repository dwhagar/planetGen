"""PERF.150: the cache warm-up, the tile prefetch and the database check run on their own queue, with workers of
their own that a generation stage of thousands of units cannot shut out (`planetgen/queue/workers.py`)."""

import os
import uuid

import pytest
import rq

from planetgen.queue import redisqueue, runrecord, workers, wrapper
from planetgen.util import settings

pytestmark = pytest.mark.real_queue


@pytest.fixture
def lane(redis_server, key_prefix, monkeypatch):
    """The two queues under names of this test's own, a key prefix of its own, three workers wanted, and the worker
    starts and the shutdown commands recorded instead of made."""
    token = uuid.uuid4().hex[:12]
    monkeypatch.setattr(redisqueue, "BULK_QUEUE", f"planetgen-bulk-test-{token}")
    monkeypatch.setattr(redisqueue, "MAINTENANCE_QUEUE", f"planetgen-maintenance-test-{token}")
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    connection = redisqueue.connect()
    started, told = [], []
    monkeypatch.setattr(workers.redisqueue, "start_detached", lambda argv, cwd: started.append(argv))
    monkeypatch.setattr(workers.command, "send_shutdown_command", lambda conn, name: told.append(name))
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 3)
    monkeypatch.setattr(workers, "host_worker_processes", lambda: 0)         # the host's own processes are not the test's
    made = []

    def register(queue_name, state="idle"):
        queue = redisqueue.queue(queue_name, connection)
        worker = rq.Worker([queue], connection=connection)
        worker.register_birth()
        connection.hset(worker.key, "pid", os.getpid())
        worker.set_state(state)
        made.append(worker)
        return worker

    yield connection, started, told, register
    for worker in made:
        worker.register_death()
    for name in (redisqueue.BULK_QUEUE, redisqueue.MAINTENANCE_QUEUE):
        redisqueue.queue(name, connection).empty()
    connection.delete(workers.REFUSED_KEY)


def _queue(connection, name, jobs):
    queue = redisqueue.queue(name, connection)
    for n in range(jobs):
        queue.enqueue(os.getpid, job_timeout=600)
    return queue


def _generation_run(connection):
    return runrecord.start_run(connection, "Fill", [runrecord.Stage("mass", "Mass", units=200)])


def test_a_run_goes_on_the_queue_of_its_kind():
    assert wrapper.queue_of(runrecord.GENERATION) == redisqueue.BULK_QUEUE
    assert wrapper.queue_of(runrecord.MAINTENANCE) == redisqueue.MAINTENANCE_QUEUE      # PERF.245: a queue for each task type of the lane
    assert wrapper.queue_of(runrecord.CHECK) == redisqueue.CHECK_QUEUE
    assert redisqueue.CHECK_QUEUE in redisqueue.QUEUES
    assert redisqueue.serializer_for(redisqueue.MAINTENANCE_QUEUE) is redisqueue.JsonSerializer
    assert redisqueue.MAINTENANCE_QUEUE in redisqueue.QUEUES


def test_perf185_the_host_cap_is_sized_from_memory_and_cores_unless_set(monkeypatch):
    """PERF.185: 60 percent of memory in 300 MB workers, no more than four a core, never fewer than cores plus four."""
    monkeypatch.setattr(settings.get_settings().queue, "max_host_workers", 0)
    monkeypatch.setattr(workers.work, "cpu_count", lambda: 8)
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: int(15.6 * 2 ** 30))
    assert workers.host_worker_cap() == 31                          # a 16 GB droplet: 31 by memory, 32 by cores
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 2 * 2 ** 30)
    assert workers.host_worker_cap() == 12                          # a small host keeps the cores plus four
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 256 * 2 ** 30)
    assert workers.host_worker_cap() == 32                          # a big one is held to four a core
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: None)
    assert workers.host_worker_cap() == 32
    monkeypatch.setattr(settings.get_settings().queue, "max_host_workers", 5)
    assert workers.host_worker_cap() == 5


def test_perf185_a_work_horse_is_not_a_second_worker():
    """PERF.185: the horse a worker forks has its command line; only the worker counts."""
    listing = "\n".join([
        "1 0 /sbin/init",
        "100 1 /venv/bin/python -m planetgen.cli.worker --url redis://x planetgen-bulk",
        "101 100 /venv/bin/python -m planetgen.cli.worker --url redis://x planetgen-bulk",       # its horse
        "200 1 /venv/bin/python -m planetgen.cli.worker --url redis://x planetgen-maintenance",
        "300 1 grep planetgen",
    ])
    assert workers.count_worker_processes(listing) == 2
