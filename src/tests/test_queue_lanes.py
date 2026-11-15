"""PERF.245: the maintenance lane has a queue for each task type served in turn, the places kept for the short and long job lanes, and the
app a worker builds once for the jobs that need it."""

import math

import pytest
import rq

from planetgen.cli import worker as worker_cli
from planetgen.queue import redisqueue, runrecord, wrapper, workers

pytestmark = pytest.mark.real_queue


def test_each_task_type_of_the_maintenance_lane_has_a_queue_of_its_own():
    assert redisqueue.maintenance_lane() == (redisqueue.MAINTENANCE_QUEUE, redisqueue.CHECK_QUEUE, redisqueue.NODES_QUEUE)
    assert len(set(redisqueue.maintenance_lane())) == 3 and set(redisqueue.maintenance_lane()) <= set(redisqueue.QUEUES)
    assert wrapper.queue_of(runrecord.MAINTENANCE) == redisqueue.MAINTENANCE_QUEUE
    assert wrapper.queue_of(runrecord.CHECK) == redisqueue.CHECK_QUEUE
    assert wrapper.queue_of(runrecord.GENERATION) == redisqueue.BULK_QUEUE
    for name in redisqueue.maintenance_lane():
        assert redisqueue.serializer_for(name) is redisqueue.JsonSerializer      # one serializer, so one worker can serve all three


def test_a_worker_for_any_queue_of_the_lane_serves_the_whole_lane():
    for name in redisqueue.maintenance_lane():
        assert workers.lane_of(name) == redisqueue.maintenance_lane()
    assert workers.lane_of(redisqueue.BULK_QUEUE) == (redisqueue.BULK_QUEUE,)
    assert workers.lane_of(redisqueue.INTERACTIVE_QUEUE) == (redisqueue.INTERACTIVE_QUEUE,)


def test_a_worker_of_several_queues_takes_one_job_from_each_in_turn():
    plain, turns = worker_cli.worker_class(), worker_cli.worker_class(round_robin=True)
    assert issubclass(turns, rq.worker.RoundRobinWorker) and not issubclass(plain, rq.worker.RoundRobinWorker)


def test_the_short_and_long_lanes_keep_places_the_bulk_and_maintenance_lanes_may_not_take(monkeypatch):
    monkeypatch.setattr(workers, "host_worker_cap", lambda: 10)
    kept = workers.reserved_places(redisqueue.BULK_QUEUE)
    assert kept >= 2 and kept < 10
    assert all(workers.reserved_places(name) == kept for name in redisqueue.maintenance_lane())
    assert workers.reserved_places(redisqueue.INTERACTIVE_QUEUE) == 0 and workers.reserved_places(redisqueue.LONG_QUEUE) == 0
    assert workers.reserved_places("planetgen-api-test-x") == 0
    monkeypatch.setattr(workers, "host_worker_cap", lambda: 1)
    assert workers.reserved_places(redisqueue.BULK_QUEUE) == 0           # a host of one place still lets the bulk lane start


@pytest.fixture
def lane_queues(redis_server):
    connection = redisqueue.connect()
    names = (*redisqueue.maintenance_lane(), redisqueue.INTERACTIVE_QUEUE)
    for name in names:
        redisqueue.queue(name, connection).empty()
    connection.delete(workers.WANTED_KEY, *(workers.STARTING_PREFIX + name for name in names))
    yield connection
    for name in names:
        redisqueue.queue(name, connection).empty()
    connection.delete(*(workers.STARTING_PREFIX + name for name in names))


def test_a_job_on_any_queue_of_the_lane_starts_workers_that_serve_all_of_them(lane_queues, monkeypatch):
    started = []
    monkeypatch.setattr(workers.redisqueue, "start_detached", lambda argv, cwd: started.append(argv))
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 2)
    monkeypatch.setattr(workers, "host_worker_processes", lambda: 0)
    monkeypatch.setattr(workers, "live_workers", lambda queue: [])
    redisqueue.queue(redisqueue.CHECK_QUEUE, lane_queues).enqueue(math.pow, 2, 3, job_timeout=600)      # only the check queue has work
    assert workers.start_workers(redisqueue.CHECK_QUEUE, lane_queues) == 1
    argv = started[0]
    assert all(name in argv for name in redisqueue.maintenance_lane()) and redisqueue.INTERACTIVE_QUEUE not in argv
    assert int(lane_queues.get(workers.STARTING_PREFIX + redisqueue.MAINTENANCE_QUEUE)) == 1       # counted once, under the lane's name


def test_a_host_whose_places_the_warm_up_used_up_still_starts_the_short_lanes_worker(lane_queues, monkeypatch):
    started = []
    monkeypatch.setattr(workers.redisqueue, "start_detached", lambda argv, cwd: started.append(argv))
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 6)
    monkeypatch.setattr(workers, "host_worker_cap", lambda: 6)
    monkeypatch.setattr(workers, "live_workers", lambda queue: [])
    monkeypatch.setattr(workers, "host_worker_processes", lambda: 4)       # four places taken: the cap less the two kept
    for _ in range(5):
        redisqueue.queue(redisqueue.MAINTENANCE_QUEUE, lane_queues).enqueue(math.pow, 2, 3, job_timeout=600)
    assert workers.start_workers(redisqueue.MAINTENANCE_QUEUE, lane_queues) == 0 and started == []
    redisqueue.queue(redisqueue.INTERACTIVE_QUEUE, lane_queues).enqueue(math.pow, 2, 3, job_timeout=600)
    assert workers.start_workers(redisqueue.INTERACTIVE_QUEUE, lane_queues, count=1) == 1
    assert redisqueue.INTERACTIVE_QUEUE in started[0]


def test_the_workers_whose_jobs_build_the_site_app_build_it_once_before_they_fork():
    assert worker_cli.APP_QUEUES == {*redisqueue.maintenance_lane(), redisqueue.INTERACTIVE_QUEUE, redisqueue.LONG_QUEUE}
    assert redisqueue.BULK_QUEUE not in worker_cli.APP_QUEUES


def test_the_whole_map_build_leaves_cores_to_the_short_lanes(monkeypatch):
    from planetgen.web import tilebuild
    monkeypatch.setattr(tilebuild.os, "cpu_count", lambda: 8)
    assert tilebuild.workers_setting() <= 6
    monkeypatch.setattr(tilebuild.os, "cpu_count", lambda: 1)
    assert tilebuild.workers_setting() == 1
