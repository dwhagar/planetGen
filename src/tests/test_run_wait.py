"""PERF.153: a run that finds the machine taken waits in line (`planetgen/queue/runwait.py`), and an API job shows
that it waits (`wrapper.report_wait`, `api_jobs.status`)."""

import rq
import pytest

from planetgen.queue import api_jobs, redisqueue, runrecord, runwait, wrapper
from planetgen.queue.runrecord import Stage


@pytest.fixture
def conn(redis_server, key_prefix, monkeypatch):
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    return redisqueue.connect(redis_server)


def _take(conn):
    return lambda: runrecord.start_run(conn, "Waiting", [Stage("mass", "Mass")])


class Clock:
    """A clock the `sleep` of a test moves, so no test waits for real."""

    def __init__(self, on_sleep=None):
        self.now = 1000.0
        self.sleeps = 0
        self.on_sleep = on_sleep

    def time(self):
        return self.now

    def sleep(self, seconds):
        self.now += seconds
        self.sleeps += 1
        if self.on_sleep:
            self.on_sleep(self.sleeps)


def test_a_free_machine_is_taken_at_once_without_a_place_in_line(conn):
    shown = []
    run_id = runwait.start_when_free(conn, _take(conn), 3600, on_wait=lambda *a: shown.append(a))
    assert runrecord.active_run(conn) == run_id and shown == [] and runwait.line(conn) == []


def test_a_run_waits_for_the_one_that_holds_the_machine_and_then_takes_it(conn):
    holder = runrecord.start_run(conn, "Holder", [Stage("mass", "Mass")])
    shown = []
    clock = Clock(lambda n: runrecord.finish_run(conn, holder, "done") if n == 3 else None)
    run_id = runwait.start_when_free(conn, _take(conn), 3600, on_wait=lambda *a: shown.append(a),
                                     clock=clock.time, sleep=clock.sleep)
    assert runrecord.active_run(conn) == run_id != holder
    assert shown == [(holder, 1), (None, None)] and runwait.line(conn) == []


def test_the_wait_ends_after_the_bound_and_names_the_run(conn):
    holder = runrecord.start_run(conn, "Holder", [Stage("mass", "Mass")])
    clock = Clock()
    with pytest.raises(runrecord.Busy) as refused:
        runwait.start_when_free(conn, _take(conn), 120, clock=clock.time, sleep=clock.sleep)
    assert holder in str(refused.value) and "waited 2 minutes" in str(refused.value)
    assert runwait.line(conn) == [] and runrecord.active_run(conn) == holder


def test_a_wait_of_zero_does_not_wait(conn):
    runrecord.start_run(conn, "Holder", [Stage("mass", "Mass")])
    clock = Clock()
    with pytest.raises(runrecord.Busy):
        runwait.start_when_free(conn, _take(conn), 0, clock=clock.time, sleep=clock.sleep)
    assert clock.sleeps == 0


def test_the_line_is_first_come_first_served_and_shows_the_place(conn):
    holder = runrecord.start_run(conn, "Holder", [Stage("mass", "Mass")])
    conn.zadd(runwait._waiters_key(), {"ahead": 1.0})
    conn.set(runwait._alive_key("ahead"), "1", ex=60)
    shown = []

    def move(n):
        if n == 2:                                   # the run ahead gives up
            conn.zrem(runwait._waiters_key(), "ahead")
        if n == 4:
            runrecord.finish_run(conn, holder, "done")

    clock = Clock(move)
    run_id = runwait.start_when_free(conn, _take(conn), 3600, on_wait=lambda *a: shown.append(a),
                                     clock=clock.time, sleep=clock.sleep)
    assert [position for _holder, position in shown[:-1]] == [2, 1] and runrecord.active_run(conn) == run_id


def test_a_waiter_that_was_killed_drops_out_of_the_line(conn):
    conn.zadd(runwait._waiters_key(), {"dead": 1.0})              # no sign of life
    assert runwait.line(conn) == []
    run_id = runwait.start_when_free(conn, _take(conn), 60)        # not held up by it
    assert runrecord.active_run(conn) == run_id


def test_a_waiting_api_job_reports_where_it_waits(conn, monkeypatch):
    class Job:
        meta = {}
        saved = 0

        def save_meta(self):
            Job.saved += 1

    monkeypatch.setattr(rq, "get_current_job", lambda: Job())
    monkeypatch.setattr(wrapper.progress_file, "report", lambda *a, **k: None)
    wrapper.report_wait("20260101-000000-abcd", 2)
    assert Job.meta["waiting"] == {"run": "20260101-000000-abcd", "position": 2,
                                   "text": "Waiting at position 2 for run 20260101-000000-abcd"}
    wrapper.report_wait(None, None)
    assert Job.meta["waiting"] is None and Job.saved == 2


def test_the_api_status_of_a_running_job_carries_its_wait(redis_server, monkeypatch):
    import math
    connection = redisqueue.connect()
    queue = redisqueue.queue(api_jobs.QUEUE_NAME, connection)
    job_id = api_jobs.new_job_id()
    job = queue.enqueue(math.pow, 2, 3, job_id=api_jobs.QUEUE_PREFIX + job_id, job_timeout=60)
    job.set_status(rq.job.JobStatus.STARTED)
    assert api_jobs.status(job_id)["waiting"] is None
    job.meta["waiting"] = {"run": "r1", "position": 1, "text": "Waiting at position 1 for run r1"}
    job.save_meta()
    body = api_jobs.status(job_id)
    assert body["state"] == "running" and body["waiting"]["position"] == 1
    job.cancel()
