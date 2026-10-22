"""PERF.152: the Queue page's Pause is one flag in Redis that the workers honour: a paused queue starts no new unit,
a resume releases it, and the two API lanes are never held."""

import math
import subprocess
import sys
import time
import uuid

import pytest
import rq

from planetgen.queue import api_jobs, pause, redisqueue

pytestmark = pytest.mark.real_queue


@pytest.fixture
def flag(redis_server):
    """The pause flag, cleared before and after the test."""
    connection = redisqueue.connect()
    pause.resume(connection)
    yield connection
    pause.resume(connection)


def test_pause_and_resume_set_and_clear_one_flag_that_names_who_and_when(flag):
    assert pause.state(flag) == {"paused": False, "by": None, "at": None}
    before = time.time()
    assert pause.pause("boss", flag) and not pause.pause("someone else", flag)       # the first pause is kept
    state = pause.state(flag)
    assert state["paused"] and state["by"] == "boss" and before <= state["at"] <= time.time()
    assert pause.resume(flag) and not pause.resume(flag)
    assert not pause.is_paused(flag)
    assert pause.key() == "planetgen:queue:paused" or pause.key().endswith(":queue:paused")


def test_without_redis_nothing_is_paused():
    class Down:
        def get(self, key):
            raise OSError("no server")
    assert pause.state(Down()) == {"paused": False, "by": None, "at": None}


def test_the_api_lanes_and_the_generate_pages_queues_are_never_held():
    assert not pause.holds(redisqueue.INTERACTIVE_QUEUE) and not pause.holds(redisqueue.LONG_QUEUE)
    assert not pause.holds("planetgen-web-20260101-000000-abcd") and not pause.holds(api_jobs.QUEUE_NAME)
    assert pause.holds(redisqueue.BULK_QUEUE) and pause.holds("planetgen-bulk-test-abc") and pause.holds("planetgen-generation")


def _run_worker(queue_name, connection, redis_url, burst=True):
    """A worker of the worker command, as its own process (a burst worker, or one that leaves after 5 idle seconds);
    returns it once it is listening on the queue."""
    argv = redisqueue.worker_argv([queue_name], redis_url, python=sys.executable, burst=burst)
    if not burst:
        argv[argv.index(queue_name):argv.index(queue_name)] = ["--idle-seconds", "5"]
    process = subprocess.Popen(argv,
                               stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    queue = redisqueue.queue(queue_name, connection)
    deadline = time.monotonic() + 60
    while rq.Worker.count(queue=queue) == 0 and process.poll() is None and time.monotonic() < deadline:
        time.sleep(0.1)
    return process


def test_a_paused_queue_starts_no_unit_and_a_resume_releases_it(flag, redis_server):
    name = f"planetgen-bulk-test-{uuid.uuid4().hex[:12]}"
    queue = redisqueue.queue(name, flag)
    job = queue.enqueue(math.pow, 2, 5, job_timeout=60)
    pause.pause("boss", flag)
    process = _run_worker(name, flag, redis_server)
    time.sleep(3)
    assert job.get_status(refresh=True) == "queued"          # held, not started
    pause.resume(flag)
    assert process.wait(60) == 0
    assert job.get_status(refresh=True) == "finished" and job.return_value() == 32.0


def test_a_job_that_arrives_after_the_pause_is_held_not_run(flag, redis_server):
    name = f"planetgen-bulk-test-{uuid.uuid4().hex[:12]}"
    queue = redisqueue.queue(name, flag)
    process = _run_worker(name, flag, redis_server, burst=False)       # blocked asking for a job, nothing paused yet
    time.sleep(0.5)
    pause.pause("boss", flag)
    job = queue.enqueue(math.pow, 2, 6, job_timeout=60)        # the blocked worker takes it, then must hold it
    time.sleep(3)
    assert job.get_status(refresh=True) != "finished" and job.get_status(refresh=True) != "started"
    pause.resume(flag)
    assert process.wait(60) == 0
    assert job.get_status(refresh=True) == "finished"


def test_a_paused_queue_does_not_hold_the_api_lanes(flag, redis_server):
    name = f"planetgen-api-test-{uuid.uuid4().hex[:12]}"
    job = redisqueue.queue(name, flag).enqueue(math.pow, 2, 3, job_timeout=60)
    pause.pause("boss", flag)
    assert _run_worker(name, flag, redis_server).wait(60) == 0
    assert job.get_status(refresh=True) == "finished"
