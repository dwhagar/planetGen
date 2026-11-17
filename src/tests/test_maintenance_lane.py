"""PERF.150: the cache warm-up, the tile prefetch and the database check run on their own queue, with workers of
their own that a generation stage of thousands of units cannot shut out (`planetgen/queue/workers.py`)."""

import os
import uuid

import pytest
import rq

from planetgen.queue import redisqueue, runrecord, supervisor, workers, wrapper
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
    connection.delete(workers.WANTED_KEY)
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
        connection.delete(workers.STARTING_PREFIX + name)
    connection.delete(workers.REFUSED_KEY, workers.REFUSED_NOTE_KEY)


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


def test_maintenance_gets_the_whole_count_when_no_generation_run_is_live(lane):
    connection, started, told, register = lane
    _queue(connection, redisqueue.MAINTENANCE_QUEUE, 10)
    assert workers.maintenance_cap(connection) == 3
    assert workers.start_workers(redisqueue.MAINTENANCE_QUEUE, connection) == 3 and len(started) == 3 and told == []


def test_maintenance_is_held_to_its_setting_while_a_generation_run_is_live(lane, monkeypatch):
    connection, started, told, register = lane
    _queue(connection, redisqueue.MAINTENANCE_QUEUE, 10)
    _generation_run(connection)
    assert workers.maintenance_cap(connection) == 1
    monkeypatch.setattr(settings.get_settings().queue, "maintenance_workers", 2)
    assert workers.maintenance_cap(connection) == 2
    assert workers.start_workers(redisqueue.MAINTENANCE_QUEUE, connection) == 2


def test_the_bulk_lane_keeps_at_least_one_worker(lane, monkeypatch):
    connection, started, told, register = lane
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 1)
    _queue(connection, redisqueue.BULK_QUEUE, 5)
    _queue(connection, redisqueue.MAINTENANCE_QUEUE, 5)
    _generation_run(connection)
    assert workers.start_workers(redisqueue.BULK_QUEUE, connection) == 1


def test_perf181_many_callers_at_once_start_no_more_workers_than_the_cap(lane):
    """PERF.181: the web app asks for upkeep workers from every request thread. A burst of threads (a map being
    dragged during a run) each read "none starting" before any wrote its start, so each started the whole cap: on
    2026-10-24 about 130 workers of 125 MB took a 16 GB host and the kernel killed MySQL."""
    import threading
    connection, started, told, register = lane
    _queue(connection, redisqueue.MAINTENANCE_QUEUE, 10)
    _generation_run(connection)
    assert workers.maintenance_cap(connection) == 1
    gate = threading.Barrier(24)

    def ask():
        own = redisqueue.connect()
        gate.wait()
        workers.ensure_workers(redisqueue.MAINTENANCE_QUEUE, own)          # OPS.209: a caller only leaves an ask

    threads = [threading.Thread(target=ask) for _ in range(24)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert started == [], "no caller started anything"
    supervisor.WorkerStarts(connection).turn()          # the one decider acts on the asks
    assert len(started) == 1


def test_perf181_a_worker_that_is_slow_to_register_is_counted_for_longer_than_a_few_seconds(lane):
    """PERF.181: a worker on a host short of memory can take far longer than 15 s to import and register. Until it
    does it counts as starting, so asking again does not start another for the same need."""
    connection, started, told, register = lane
    _queue(connection, redisqueue.MAINTENANCE_QUEUE, 10)
    _generation_run(connection)
    assert workers.start_workers(redisqueue.MAINTENANCE_QUEUE, connection) == 1
    assert connection.ttl(workers.STARTING_PREFIX + redisqueue.MAINTENANCE_QUEUE) > 60
    assert workers.start_workers(redisqueue.MAINTENANCE_QUEUE, connection) == 0
    assert len(started) == 1


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


def test_perf185_many_callers_together_never_pass_the_host_cap(lane, monkeypatch):
    """PERF.185: workers of every queue count against one cap; a start past it is refused, counted and logged once."""
    import threading
    connection, started, told, register = lane
    monkeypatch.setattr(workers, "host_worker_cap", lambda: 4)
    monkeypatch.setattr(workers, "host_worker_processes", lambda: len(started))       # a started worker is a process
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 10)
    _queue(connection, redisqueue.BULK_QUEUE, 50)
    _queue(connection, redisqueue.MAINTENANCE_QUEUE, 50)
    lines = []
    monkeypatch.setattr(workers.log, "normal", lines.append)
    gate = threading.Barrier(16)

    def ask(name):
        own = redisqueue.connect()
        gate.wait()
        for _ in range(3):
            workers.ensure_workers(name, own)

    threads = [threading.Thread(target=ask, args=(name,)) for name in
               [redisqueue.BULK_QUEUE, redisqueue.MAINTENANCE_QUEUE] * 8]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert started == []                                      # OPS.209: the callers only left asks
    supervisor.WorkerStarts(connection).turn()
    assert len(started) <= 4                                  # never past the cap
    for _ in range(10):                                       # the asks come again at the callers' next poll
        for name in (redisqueue.BULK_QUEUE, redisqueue.MAINTENANCE_QUEUE):
            workers.ensure_workers(name, connection)
        supervisor.WorkerStarts(connection).turn()
    assert len(started) == 4        # the whole of the host cap, no places kept back (OPS.231: the short and long lanes are light jobs)
    assert int(connection.get(workers.REFUSED_KEY)) >= 1
    assert len(lines) == 1 and "cap is 4" in lines[0]                                 # said once, not at every poll
