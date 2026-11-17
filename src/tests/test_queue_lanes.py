"""PERF.245, OPS.229: the maintenance lane has a queue for each task type (a worker serves one), the places kept for the short and long job lanes, and the
app the lane's template builds once for the units that need it."""

import math

import pytest
import rq

from planetgen.cli import worker as worker_cli
from planetgen.queue import redisqueue, runrecord, wrapper, workers

pytestmark = pytest.mark.real_queue


def test_each_task_type_of_the_maintenance_lane_has_a_queue_of_its_own():
    assert redisqueue.maintenance_lane() == (redisqueue.MAINTENANCE_QUEUE, redisqueue.CHECK_QUEUE)
    assert len(set(redisqueue.maintenance_lane())) == 2 and set(redisqueue.maintenance_lane()) <= set(redisqueue.QUEUES)
    assert wrapper.queue_of(runrecord.MAINTENANCE) == redisqueue.MAINTENANCE_QUEUE
    assert wrapper.queue_of(runrecord.CHECK) == redisqueue.CHECK_QUEUE
    assert wrapper.queue_of(runrecord.GENERATION) == redisqueue.BULK_QUEUE
    for name in redisqueue.maintenance_lane():
        assert redisqueue.serializer_for(name) is redisqueue.JsonSerializer      # one serializer for the whole lane


def test_a_worker_serves_one_queue_and_the_round_robin_worker_is_gone():
    for name in (*redisqueue.maintenance_lane(), redisqueue.BULK_QUEUE):
        assert workers.lane_of(name) == (name,)
    assert not issubclass(worker_cli.worker_class(), rq.worker.RoundRobinWorker)
    assert not hasattr(worker_cli, "RoundRobinWorker") and not hasattr(worker_cli, "APP_QUEUES")


@pytest.fixture
def lane_queues(redis_server):
    connection = redisqueue.connect()
    names = redisqueue.maintenance_lane()
    for name in names:
        redisqueue.queue(name, connection).empty()
    connection.delete(workers.WANTED_KEY, *(workers.STARTING_PREFIX + name for name in names))
    yield connection
    for name in names:
        redisqueue.queue(name, connection).empty()
    connection.delete(*(workers.STARTING_PREFIX + name for name in names))


def test_a_job_on_a_queue_of_the_lane_starts_a_worker_for_that_queue_only(lane_queues, monkeypatch):
    started = []
    monkeypatch.setattr(workers.redisqueue, "start_detached", lambda argv, cwd: started.append(argv))
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 2)
    monkeypatch.setattr(workers, "host_worker_processes", lambda: 0)
    monkeypatch.setattr(workers, "live_workers", lambda queue: [])
    redisqueue.queue(redisqueue.CHECK_QUEUE, lane_queues).enqueue(math.pow, 2, 3, job_timeout=600)      # only the check queue has work
    assert workers.start_workers(redisqueue.CHECK_QUEUE, lane_queues) == 1
    argv = started[0]
    assert redisqueue.CHECK_QUEUE in argv and redisqueue.MAINTENANCE_QUEUE not in argv
    assert int(lane_queues.get(workers.STARTING_PREFIX + redisqueue.CHECK_QUEUE)) == 1


def test_a_host_whose_places_the_warm_up_used_up_starts_no_more_of_its_workers(lane_queues, monkeypatch):
    started = []
    monkeypatch.setattr(workers.redisqueue, "start_detached", lambda argv, cwd: started.append(argv))
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 6)
    monkeypatch.setattr(workers, "host_worker_cap", lambda: 6)
    monkeypatch.setattr(workers, "live_workers", lambda queue: [])
    monkeypatch.setattr(workers, "host_worker_processes", lambda: 6)       # every place the cap gives is taken
    for _ in range(5):
        redisqueue.queue(redisqueue.MAINTENANCE_QUEUE, lane_queues).enqueue(math.pow, 2, 3, job_timeout=600)
    assert workers.start_workers(redisqueue.MAINTENANCE_QUEUE, lane_queues) == 0 and started == []


def test_the_template_of_the_maintenance_lane_builds_the_site_app_once_and_no_other_does():
    from planetgen.queue import allowlist, jobhost, launcher
    assert jobhost.WARM_SETS["app"] == ("planetgen.web.warmup.warm_app",)
    assert launcher.TEMPLATES == {"maintenance": "app"} and all(key in jobhost.WARM_SETS for key in launcher.TEMPLATES.values())
    assert jobhost.WARM_SETS["plain"] == () and "planetgen.web.warmup.warm_app" not in jobhost.WARM_SETS["tables"]
    allowlist.check_function("planetgen.web.warmup.warm_app", "a template's warm set")


def test_the_whole_map_build_leaves_cores_to_the_short_lanes(monkeypatch):
    from planetgen.web import tilebuild
    monkeypatch.setattr(tilebuild.os, "cpu_count", lambda: 8)
    assert tilebuild.workers_setting() <= 6
    monkeypatch.setattr(tilebuild.os, "cpu_count", lambda: 1)
    assert tilebuild.workers_setting() == 1
