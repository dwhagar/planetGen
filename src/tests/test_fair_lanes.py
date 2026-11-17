"""PERF.151, OPS.230: an API job expected to run long goes to a lane of its own, so it can never hold the room kept for short jobs, and the short jobs run at a higher priority than the bulk and long ones."""

import math

import pytest

from planetgen.queue import api_jobs, lightjobs, memorygate, redisqueue, work
from planetgen.util import settings
from tests import inline_queue


@pytest.fixture
def lanes(redis_server, monkeypatch):
    """The supervisor asked for, and not acting: the jobs stay on the line."""
    monkeypatch.setattr(lightjobs, "ask_supervisor", lambda connection=None, url=None: None)
    return redisqueue.connect()


def _lane(connection, job_id):
    return connection.hget(lightjobs.job_key(job_id), "lane").decode()


def test_a_job_expected_to_run_over_the_limit_goes_to_the_long_lane(lanes):
    lanes_of = [_lane(lanes, api_jobs.submit(math.pow, 2, n, expected_seconds=expected)) for n, expected in ((3, None), (4, 60), (5, 60.5), (6, 3 * 86400))]
    assert lanes_of == ["interactive", "interactive", "long", "long"]
    assert not lightjobs.is_long(None) and not lightjobs.is_long(float("nan")) and lightjobs.is_long(1e9)


def test_the_long_job_limit_is_a_setting(lanes, monkeypatch):
    monkeypatch.setattr(settings.get_settings().queue, "long_job_seconds", 5.0)
    assert lightjobs.is_long(6) and not lightjobs.is_long(5)
    assert lightjobs.lane_of(6) == memorygate.LANE_LONG and lightjobs.lane_of(5) == memorygate.LANE_INTERACTIVE


def test_the_light_job_of_the_short_lane_runs_at_a_higher_priority_than_the_bulk_and_long_ones():
    assert lightjobs.PROCESS_NICE[memorygate.LANE_INTERACTIVE] < lightjobs.PROCESS_NICE[memorygate.LANE_LONG] == work.WORKER_NICENESS


def test_a_long_job_running_leaves_the_room_kept_for_a_short_one(lanes):
    """The short and long lanes are judged against the places kept for them (`dispatch.dispatch`), not against each other's claims: a long job that holds its claim does not keep a short one waiting."""
    held = []
    duty = lightjobs.LightJobs(lanes, inline_queue.InlineLightForker(), read_host=lambda: memorygate.Host(16_000, 6_000, 4, 0.5), held_mb=lambda pid: 0.0, also_held=lambda: held)
    held.append(memorygate.Holder("long-1", memorygate.Claim("api:x", 2_000, 1.0, memorygate.LANE_LONG), 1, "light", 2_000, 2_000))
    short = api_jobs.submit(math.pow, 2, 5)
    assert duty.turn() == [short]
    duty.ended()
    assert api_jobs.status(short)["result"] == 32.0
