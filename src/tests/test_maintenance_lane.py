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
    for worker in made:
        connection.delete(workers.LENT_PREFIX + worker.name)


def _queue(connection, name, jobs):
    queue = redisqueue.queue(name, connection)
    for n in range(jobs):
        queue.enqueue(os.getpid, job_timeout=600)
    return queue


def _generation_run(connection):
    return runrecord.start_run(connection, "Fill", [runrecord.Stage("mass", "Mass", units=200)])


def test_a_run_goes_on_the_queue_of_its_kind():
    assert wrapper.queue_of(runrecord.GENERATION) == redisqueue.BULK_QUEUE
    assert wrapper.queue_of(runrecord.MAINTENANCE) == wrapper.queue_of(runrecord.CHECK) == redisqueue.MAINTENANCE_QUEUE
    assert redisqueue.serializer_for(redisqueue.MAINTENANCE_QUEUE) is redisqueue.JsonSerializer
    assert redisqueue.MAINTENANCE_QUEUE in redisqueue.QUEUES


def test_maintenance_work_behind_200_bulk_units_gets_a_worker_by_asking_a_bulk_worker_to_leave(lane):
    connection, started, told, register = lane
    bulk = [register(redisqueue.BULK_QUEUE, "busy") for _ in range(3)]            # every worker of the host is on bulk
    _queue(connection, redisqueue.BULK_QUEUE, 200)
    _queue(connection, redisqueue.MAINTENANCE_QUEUE, 1)
    _generation_run(connection)
    assert workers.ensure_workers(redisqueue.MAINTENANCE_QUEUE, connection) == 1
    assert len(started) == 1 and redisqueue.MAINTENANCE_QUEUE in started[0]       # the unit does not wait for 200 others
    assert len(told) == 1 and told[0] in {worker.name for worker in bulk}          # a bulk worker finishes and leaves
    workers.ensure_workers(redisqueue.MAINTENANCE_QUEUE, connection)
    assert len(started) == 1 and len(told) == 1                                    # not again for the same need
    # the bulk lane does not take the place back while the maintenance lane needs it
    assert workers.ensure_workers(redisqueue.BULK_QUEUE, connection) == 0
    leaving = next(worker for worker in bulk if worker.name == told[0])
    leaving.register_death()
    assert workers.ensure_workers(redisqueue.BULK_QUEUE, connection) == 0         # 2 bulk workers, cap 3 - 1


def test_after_the_maintenance_lane_drains_the_bulk_pool_is_back_to_wanted(lane):
    connection, started, told, register = lane
    bulk = [register(redisqueue.BULK_QUEUE, "busy") for _ in range(2)]
    _queue(connection, redisqueue.BULK_QUEUE, 200)
    upkeep = register(redisqueue.MAINTENANCE_QUEUE, "busy")
    _generation_run(connection)
    assert workers.maintenance_share(connection) == 1
    assert workers.ensure_workers(redisqueue.BULK_QUEUE, connection) == 0         # 2 live, cap 2
    upkeep.set_state("idle")                                                       # drained: waiting out its idle time
    assert workers.maintenance_share(connection) == 0
    assert workers.ensure_workers(redisqueue.BULK_QUEUE, connection) == 1         # back to wanted (3)
    assert len(started) == 1 and redisqueue.BULK_QUEUE in started[0] and len(bulk) == 2


def test_maintenance_gets_the_whole_count_when_no_generation_run_is_live(lane):
    connection, started, told, register = lane
    _queue(connection, redisqueue.MAINTENANCE_QUEUE, 10)
    assert workers.maintenance_cap(connection) == 3
    assert workers.ensure_workers(redisqueue.MAINTENANCE_QUEUE, connection) == 3 and len(started) == 3 and told == []


def test_maintenance_is_held_to_its_setting_while_a_generation_run_is_live(lane, monkeypatch):
    connection, started, told, register = lane
    _queue(connection, redisqueue.MAINTENANCE_QUEUE, 10)
    _generation_run(connection)
    assert workers.maintenance_cap(connection) == 1
    monkeypatch.setattr(settings.get_settings().queue, "maintenance_workers", 2)
    assert workers.maintenance_cap(connection) == 2
    assert workers.ensure_workers(redisqueue.MAINTENANCE_QUEUE, connection) == 2


def test_the_bulk_lane_keeps_at_least_one_worker(lane, monkeypatch):
    connection, started, told, register = lane
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 1)
    _queue(connection, redisqueue.BULK_QUEUE, 5)
    _queue(connection, redisqueue.MAINTENANCE_QUEUE, 5)
    _generation_run(connection)
    assert workers.ensure_workers(redisqueue.BULK_QUEUE, connection) == 1


def test_shed_bulk_asks_idle_workers_first_and_never_twice(lane):
    connection, started, told, register = lane
    busy = register(redisqueue.BULK_QUEUE, "busy")
    idle = register(redisqueue.BULK_QUEUE, "idle")
    assert workers.shed_bulk(connection, 1) == 1 and told == [idle.name]
    assert workers.shed_bulk(connection, 1) == 0 and told == [idle.name]            # it is leaving: counted as gone
    assert workers.shed_bulk(connection, 0) == 1 and told == [idle.name, busy.name]
