"""PERF.245, OPS.229: the maintenance lane has a queue for each task type (a worker serves one), the places kept for the short and long job lanes, and the
app the lane's template builds once for the units that need it."""


import pytest

from planetgen.queue import redisqueue, runrecord, wrapper

pytestmark = pytest.mark.real_queue


def test_each_task_type_of_the_maintenance_lane_has_a_queue_of_its_own():
    assert redisqueue.maintenance_lane() == (redisqueue.MAINTENANCE_QUEUE, redisqueue.CHECK_QUEUE)
    assert len(set(redisqueue.maintenance_lane())) == 2 and set(redisqueue.maintenance_lane()) <= set(redisqueue.QUEUES)
    assert wrapper.queue_of(runrecord.MAINTENANCE) == redisqueue.MAINTENANCE_QUEUE
    assert wrapper.queue_of(runrecord.CHECK) == redisqueue.CHECK_QUEUE
    assert wrapper.queue_of(runrecord.GENERATION) == redisqueue.BULK_QUEUE


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
