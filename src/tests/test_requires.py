"""PERF.247: a job can require other jobs; the queue never hands a worker one that cannot run yet, a failed requirement fails what waits
for it with a clear line, and jobs that declare no requirement still run together."""

import subprocess

import pytest
import rq

from planetgen.queue import commands, redisqueue, requires, runrecord
from planetgen.queue.runrecord import Stage

pytestmark = pytest.mark.real_queue

JOBS = "tests.requires_jobs.mark"
ORDER = "test:requires:order"


@pytest.fixture
def lane(redis_server):
    connection = redisqueue.connect(redis_server)
    name = redisqueue.CHECK_QUEUE
    queue = redisqueue.queue(name, connection)
    queue.empty()
    requires.drop_held(connection, name)
    connection.delete(ORDER)
    yield connection, queue
    requires.drop_held(connection, name)
    queue.empty()
    connection.delete(ORDER)


def order(connection):
    return [item.decode() for item in connection.lrange(ORDER, 0, -1)]


def burst(queue, workers=1):
    """Runs `workers` burst workers on `queue` until it is empty."""
    url = queue.connection.connection_pool.connection_kwargs
    address = f"redis://{url.get('host', '127.0.0.1')}:{url.get('port', 6379)}/{url.get('db', 0)}"
    processes = [subprocess.Popen(redisqueue.worker_argv([queue.name], address, burst=True)) for _ in range(workers)]
    for process in processes:
        process.wait(timeout=120)


def queue_job(queue, label, **options):
    return requires.enqueue(queue, JOBS, label, **options.pop("kwargs", {}), job_timeout=60, result_ttl=0, failure_ttl=600, **options)


def test_a_held_job_is_in_no_queue_and_runs_after_the_job_it_requires(lane):
    connection, queue = lane
    first = queue_job(queue, "first", job_id="req-first")
    second = queue_job(queue, "second", job_id="req-second", requires=[first])
    assert second.get_status() == rq.job.JobStatus.DEFERRED
    assert queue.job_ids == ["req-first"]                       # the held job is in no queue, so no round robin can pick it
    assert requires.held_jobs(connection, queue.name) == [("req-second", ["req-first"])]
    burst(queue)
    assert order(connection) == ["start first", "end first", "start second", "end second"]


def test_a_requirement_that_failed_fails_what_waits_for_it_and_the_chain_after_it_with_a_line(lane):
    connection, queue = lane
    bad = queue_job(queue, "bad", job_id="req-bad", kwargs={"fail": True})
    waits = queue_job(queue, "waits", job_id="req-waits", requires=[bad])
    queue_job(queue, "chain", job_id="req-chain", requires=[waits])
    queue_job(queue, "free", job_id="req-free")
    burst(queue)
    assert order(connection).count("start waits") == 0 and "start chain" not in order(connection)
    assert "end free" in order(connection)                         # a job with no requirement was not held back
    registry = rq.registry.FailedJobRegistry(queue.name, connection=connection, serializer=queue.serializer)
    assert set(registry.get_job_ids()) >= {"req-bad", "req-waits", "req-chain"}
    for job_id, text in (("req-waits", "required job req-bad failed"), ("req-chain", "required job req-waits did not run: required job req-bad failed")):
        assert text in redisqueue.fetch_job(job_id, connection, queue.serializer)._exc_info
    assert requires.held_jobs(connection, queue.name) == []         # nothing waits for ever


def test_an_optional_requirement_that_failed_still_lets_the_job_run(lane):
    connection, queue = lane
    bad = queue_job(queue, "bad", job_id="req-bad", kwargs={"fail": True})
    queue_job(queue, "after", job_id="req-after", requires=[bad], optional=True)
    burst(queue)
    assert "end after" in order(connection)


def test_a_job_that_requires_a_job_already_done_and_gone_is_not_held(lane):
    connection, queue = lane
    job = queue_job(queue, "late", job_id="req-late", requires=["req-never-existed"])
    assert job.get_status() == rq.job.JobStatus.QUEUED and queue.job_ids == ["req-late"]


def test_two_jobs_that_declare_no_requirement_run_together(lane):
    connection, queue = lane
    queue_job(queue, "one", job_id="req-one", kwargs={"seconds": 3.0})
    queue_job(queue, "two", job_id="req-two", kwargs={"seconds": 3.0})
    burst(queue, workers=2)
    seen = order(connection)
    assert seen.index("start two") < seen.index("end one") and seen.index("start one") < seen.index("end two")      # they overlapped


def test_dropping_the_held_jobs_of_a_queue_leaves_none_waiting(lane):
    connection, queue = lane
    first = queue_job(queue, "first", job_id="req-first")
    queue_job(queue, "second", job_id="req-second", requires=[first])
    queue_job(queue, "third", job_id="req-third", requires=[first])
    assert requires.drop_held(connection, queue.name) == 2
    assert requires.held_jobs(connection, queue.name) == []
    burst(queue)
    assert order(connection) == ["start first", "end first"]


def test_a_requirement_must_share_the_serializer_of_the_job_that_requires_it(lane):
    connection, queue = lane
    other = redisqueue.queue(redisqueue.GENERATION_QUEUE, connection)
    if redisqueue.serializer_for(redisqueue.GENERATION_QUEUE) is queue.serializer:
        pytest.skip("both queues share one serializer here")
    held = queue_job(other, "pickled", job_id="req-pickled")
    with pytest.raises(ValueError, match="serializer"):
        queue_job(queue, "json", job_id="req-json", requires=[held])
    other.empty()


def test_a_job_cannot_require_itself_or_name_an_empty_requirement(lane):
    connection, queue = lane
    with pytest.raises(ValueError, match="itself"):
        queue_job(queue, "self", job_id="req-self", requires=["req-self"])
    with pytest.raises(ValueError, match="job id"):
        queue_job(queue, "empty", job_id="req-empty", requires=[""])


def test_the_steps_of_a_run_come_after_the_steps_they_require():
    def step(key, kind):
        return commands.Step(key, kind, key, {})
    ok = [step("phenomena", "scatter-phenomena"), step("phenomena-index", "finish-phenomena")]
    assert commands.validate_order(ok) == ok
    assert commands.required_stage_keys(ok, 1) == ["phenomena"] and commands.required_stage_keys(ok, 0) == []
    with pytest.raises(ValueError, match="comes after it"):
        commands.validate_order(list(reversed(ok)))
    assert commands.validate_order([step("phenomena-index", "finish-phenomena")])        # a requirement that is not in the run is nothing to wait for


def test_a_failed_run_says_which_stage_did_not_run_because_the_stage_it_requires_did_not_finish(redis_server, key_prefix, monkeypatch):
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    connection = redisqueue.connect(redis_server)
    steps = [commands.Step("phenomena", "scatter-phenomena", "Scatter the phenomena", {}),
             commands.Step("phenomena-index", "finish-phenomena", "Index the phenomena", {})]
    run_id = runrecord.start_run(connection, "Plan", [Stage("phenomena", "Scatter the phenomena", units=2), Stage("phenomena-index", "Index the phenomena", units=1)],
                                 steps=[step.as_dict() for step in steps], now=100)
    runrecord.begin_stage(connection, run_id, 0, now=101)
    runrecord.finish_run(connection, run_id, "failed", error="unit 1 of stage Scatter the phenomena failed: boom", now=102)
    lines = [fields["message"] for _id, fields in runrecord.read_events(connection, run_id)]
    assert any("dropped: stage Index the phenomena did not run" in line and "phenomena" in line and "failed" in line for line in lines)
