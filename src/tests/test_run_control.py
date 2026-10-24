"""PERF.99: what an admin does to a run in Redis (`queue/runcontrol.py`): cancel through the stop command, clear a run
that stopped answering, delete a finished one, and keep its log. The autouse `_run_queue_inline` fixture gives each test
a queue and key prefix of its own."""

import logging
import os
import time

import pytest
import rq

from planetgen.queue import commands, redisqueue, runcontrol, runrecord, units, wrapper
from planetgen.queue.runrecord import Stage

KIND = "control-test"


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture(autouse=True)
def _kind(monkeypatch):
    monkeypatch.setitem(commands.REGISTRY, KIND, commands.CommandSpec(KIND, "A test", unit=(__name__, "_noop"),
                                                                     expand=(__name__, "_cut")))


def _noop(args, report):
    return {}


def _cut(params):
    return commands.Expansion([({"n": n}, 1.0) for n in range(params["count"])], {})


@pytest.fixture
def jobs_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))
    return tmp_path / "jobs"


def _queued_run(conn, count=75, done=0, stages=1):
    """A live run whose first stage has `count` units waiting on the bulk queue, `done` of them finished."""
    run_id = runrecord.start_run(conn, "Control test", [Stage(f"s{n}", f"Stage {n}", units=count) for n in range(stages)],
                                 who="boss", now=100.0)
    runrecord.begin_stage(conn, run_id, 0, now=101.0)
    units.enqueue_units(conn, run_id, 0, KIND, [({"n": n}, 1.0) for n in range(count)])
    for unit in range(done):
        job = rq.job.Job.fetch(units.job_id(run_id, 0, unit), connection=conn, serializer=redisqueue.JsonSerializer)
        job.delete()
        runrecord.unit_done(conn, run_id, 0, unit, now=102.0)
    return run_id


def _waiting(conn):
    return redisqueue.queue(redisqueue.BULK_QUEUE, conn).count


def test_cancelling_a_run_cancels_its_waiting_units_and_frees_the_machine_in_under_a_second(conn, jobs_dir):
    run_id = _queued_run(conn, count=75, done=5)
    assert _waiting(conn) == 70 and runrecord.active_run(conn) == run_id
    began = time.monotonic()
    stopped = runcontrol.cancel_run(conn, run_id, who="boss")
    took = time.monotonic() - began
    assert stopped == {"cancelled": 70, "stopped": 0}
    assert took < 1.0, f"cancelling a 75 unit run took {took:.2f} s"
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "cancelled" and run["error"] == "cancelled by boss"
    assert _waiting(conn) == 0
    assert runrecord.active_run(conn) is None, "the machine is free"
    kinds = [fields["kind"] for _id, fields in runrecord.read_events(conn, run_id)]
    assert "run.cancelling" in kinds and kinds[-1] == "run.cancelled"


def test_a_started_unit_is_sent_the_stop_command(conn, jobs_dir):
    run_id = _queued_run(conn, count=4)
    running = units.job_id(run_id, 0, 2)
    job = rq.job.Job.fetch(running, connection=conn, serializer=redisqueue.JsonSerializer)
    queue = redisqueue.queue(redisqueue.BULK_QUEUE, conn)
    queue.remove(job)
    job.worker_name = "worker-for-the-test"
    job.set_status(rq.job.JobStatus.STARTED)
    job.save()
    conn.zadd(queue.started_job_registry.key, {f"{running}:execution": time.time() + 60})
    listener = conn.pubsub()
    listener.subscribe("rq:pubsub:worker-for-the-test")
    listener.get_message(timeout=1)
    assert runcontrol.cancel_run(conn, run_id, who="boss") == {"cancelled": 3, "stopped": 1}
    message = None
    for _ in range(20):
        message = listener.get_message(timeout=0.2)
        if message and message["type"] == "message":
            break
    assert message and message["type"] == "message" and b"stop-job" in message["data"] and running.encode() in message["data"]


def test_a_cancelled_run_queues_no_next_stage(conn, jobs_dir):
    run_id = _queued_run(conn, count=2, stages=2)
    runrecord.request_cancel(conn, run_id)
    wrapper.advance(conn, run_id, 1)
    assert runrecord.read_run(conn, run_id)["stages"][1]["started"] is None
    assert _waiting(conn) == 2, "only the first stage's units are queued"


def test_a_run_that_is_over_or_unknown_is_not_cancelled_again(conn, jobs_dir):
    run_id = _queued_run(conn, count=2)
    assert runcontrol.cancel_run(conn, run_id) is not None
    assert runcontrol.cancel_run(conn, run_id) is None
    assert runcontrol.cancel_run(conn, "20260101-000000-deadbeef") is None
    assert runrecord.read_run(conn, run_id)["error"] == "cancelled"        # no one named


def test_a_run_that_shows_life_is_not_cleared_and_a_stale_one_is(conn, jobs_dir):
    run_id = _queued_run(conn, count=3)
    runrecord.touch(conn, run_id)
    assert runcontrol.clear_run(conn, run_id, who="boss", now=105.0) is None
    assert runrecord.read_run(conn, run_id)["state"] == "running"
    later = runrecord.read_run(conn, run_id)["updated"] + 10_000
    assert runcontrol.clear_run(conn, run_id, who="boss", now=later) == {"cancelled": 3, "stopped": 0}
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "failed" and run["error"].startswith("cleared by boss")
    assert runrecord.active_run(conn) is None
    assert _waiting(conn) == 0


def test_a_run_whose_follower_is_gone_can_be_cleared_at_once(conn, jobs_dir):
    run_id = _queued_run(conn, count=1)
    conn.hset(runrecord.run_key(run_id), "owner", "some-other-host:1")
    conn.delete(runrecord.alive_key(run_id))
    assert runcontrol.clear_run(conn, run_id, who="boss") is not None
    assert runrecord.read_run(conn, run_id)["state"] == "failed"


def test_only_a_finished_run_is_deleted(conn, jobs_dir):
    live = _queued_run(conn, count=1)
    assert runrecord.delete_run(conn, live) is False
    runcontrol.cancel_run(conn, live)
    assert runrecord.delete_run(conn, live) is True
    assert runrecord.read_run(conn, live) is None
    assert not [key for key in conn.scan_iter(match=f"{runrecord.PREFIX}:run:{live}*")]
    assert runrecord.delete_run(conn, live) is False


def test_clearing_the_finished_runs_keeps_the_live_one(conn, jobs_dir):
    first = _queued_run(conn, count=1)
    runcontrol.cancel_run(conn, first)
    live = _queued_run(conn, count=1)
    assert runcontrol.clear_finished(conn) == 1
    assert runrecord.read_run(conn, live)["state"] == "running" and runrecord.read_run(conn, first) is None


def test_the_log_is_written_where_the_download_link_finds_it_and_the_newest_twenty_are_kept(conn, jobs_dir):
    run_ids = []
    for number in range(22):
        run_id = runrecord.start_run(conn, "Log test", [Stage("s", "S", units=1)], run_id=f"20260101-0000{number:02d}-aaaaaaaa",
                                     now=float(number))
        runrecord.log(conn, run_id, f"a line of run {number}", now=float(number))
        runrecord.finish_run(conn, run_id, "done", now=float(number) + 1)
        assert runcontrol.write_log(conn, run_id) == runcontrol.log_path(run_id)
        run_ids.append(run_id)
    logs = sorted(os.listdir(jobs_dir / "runs"))
    assert logs == sorted(f"{run_id}.log" for run_id in run_ids[2:])
    text = (jobs_dir / "runs" / f"{run_ids[-1]}.log").read_text()
    assert "a line of run 21" in text and "run started: Log test" in text and text.rstrip().endswith("event=run.done")


def test_a_run_started_by_a_job_is_found_from_the_job(conn, jobs_dir):
    first = runrecord.start_run(conn, "One", [Stage("s", "S", units=1)], job="20261024-010203-abcd")
    runrecord.finish_run(conn, first, "done")
    second = runrecord.start_run(conn, "Two", [Stage("s", "S", units=1)], job="20261024-010203-abcd")
    other = runrecord.read_run(conn, second)
    assert other["job"] == "20261024-010203-abcd"
    assert runrecord.job_runs(conn, "20261024-010203-abcd") == [first, second]
    assert runrecord.job_runs(conn, "20261024-999999-ffff") == []
    assert runrecord.read_run(conn, first)["argv"] == []


def test_what_the_follower_logs_joins_the_event_stream(conn, caplog):
    caplog.set_level(logging.INFO, logger="planetgen")        # whatever level an earlier test left on the logger
    run_id = runrecord.start_run(conn, "Tee", [Stage("s", "S", units=1)])
    handler = runrecord.EventHandler(conn, run_id)
    logger = logging.getLogger("planetgen")
    logger.addHandler(handler)
    try:
        logger.info("sector 7 saved")
        logger.debug("not a line for the log")
        logger.warning("slow disk")
    finally:
        logger.removeHandler(handler)
    lines = [(fields["level"], fields["message"]) for _id, fields in runrecord.read_events(conn, run_id)
             if fields["kind"] == "log"]
    assert lines == [("info", "sector 7 saved"), ("warning", "slow disk")]


def test_a_run_through_the_wrapper_records_its_job_logs_its_lines_and_leaves_a_log_file(conn, jobs_dir, monkeypatch, caplog):
    caplog.set_level(logging.INFO, logger="planetgen")
    monkeypatch.setenv(wrapper.JOB_ENV, "20261024-040506-beef")
    monkeypatch.setattr(commands, "unit_function", lambda kind: _noop)
    monkeypatch.setattr(commands, "expand_function", lambda kind: _cut)

    def result(index, unit, outcome):
        logging.getLogger("planetgen").info(f"unit {unit} came back")

    run_id = wrapper.run("Through the wrapper", [commands.Step("s", KIND, "A stage", {"count": 3})],
                         lambda index, plan: None, result)
    assert runrecord.job_runs(conn, "20261024-040506-beef") == [run_id]
    assert runrecord.read_run(conn, run_id)["state"] == "done"
    text = (jobs_dir / "runs" / f"{run_id}.log").read_text()
    assert text.count("came back") == 3 and 'msg="Stage 1 of 1 finished after ' in text and "event=stage.finish" in text
