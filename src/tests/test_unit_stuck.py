"""OPS.140: a stuck unit is dumped, stopped and queued again once; stuck again, it fails the run with both stacks; a unit whose counter moves is never touched."""

import json
import os
import subprocess
import sys
import textwrap


import pytest
import rq

from planetgen.queue import beatjudge, commands, pause, redisqueue, runrecord, unitbeat, units, unitstuck


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture(autouse=True)
def _no_hold_left(conn):
    pause.resume(conn)
    yield
    pause.resume(conn)


@pytest.fixture
def run_id(conn, monkeypatch):
    monkeypatch.setitem(commands.REGISTRY, "stuckkind", commands.CommandSpec("stuckkind", "A test", unit=(__name__, "_unit"), expand=(__name__, "_unit")))
    run = runrecord.start_run(conn, "T", [runrecord.Stage("stuckkind", "Stuck", units=2)])
    runrecord.begin_stage(conn, run, 0)
    units.enqueue_units(conn, run, 0, "stuckkind", [({"n": 0}, 1.0), ({"n": 1}, 1.0)])
    yield run
    redisqueue.queue(redisqueue.BULK_QUEUE, conn).empty()


def _unit(args, report):
    return {}


def _beat(conn, run_id, unit, pid=999_999_999, host="elsewhere"):
    now = unitbeat.server_time(conn)
    conn.hset(runrecord.progress_key(run_id, 0), unit, json.dumps({"done": 3.0, "estimate": 9.0, "beat": now, "started": now - 3000, "moved": now - 2900, "phase": "p",
                                                                    "waiting_on": "", "restarts": 0, "pid": pid, "host": host}))


def _quiet():
    return beatjudge.Verdict(beatjudge.QUIET, reason="the count has stood at 3 for 2900 s")


def _stuck(conn, run_id, errors=None, stopped=None, forgotten=None, stack="STACK"):
    return unitstuck.Stuck(conn, run_id, 0, dump=lambda raw: stack, stop=lambda job_id: (stopped if stopped is not None else []).append(job_id),
                           record_error=lambda message, **fields: (errors if errors is not None else []).append(message),
                           forget=lambda unit: (forgotten if forgotten is not None else []).append(unit))


def _job(conn, run_id, unit):
    return redisqueue.fetch_job(units.job_id(run_id, 0, unit), conn)


def test_a_stuck_unit_is_dumped_and_stopped_and_the_run_notes_it(conn, run_id):
    stopped, errors = [], []
    action = _stuck(conn, run_id, errors, stopped)
    _beat(conn, run_id, 1)
    action.step([1], {1: _quiet()})
    assert stopped == [units.job_id(run_id, 0, 1)]
    record = action.records()[1]
    assert record["phase"] == "stopping" and record["attempts"] == 1 and record["dump"] == "STACK"
    run = conn.hgetall(runrecord.run_key(run_id))
    assert run[b"killed"] == b"1" and b"quiet" in run[b"killed_reason"]
    assert runrecord.read_run(conn, run_id)["state"] == "running" and errors and "STACK" in errors[0]
    action.step([1], {1: _quiet()})
    assert stopped == [units.job_id(run_id, 0, 1)], "a stop under way is not repeated"


def test_a_stopped_unit_is_queued_again_once_it_has_left_its_horse(conn, run_id):
    forgotten = []
    action = _stuck(conn, run_id, forgotten=forgotten)
    _beat(conn, run_id, 1)
    action.step([1], {1: _quiet()})
    job = _job(conn, run_id, 1)
    queue = redisqueue.queue(redisqueue.BULK_QUEUE, conn)
    queue.remove(job.id)
    job.set_status(rq.job.JobStatus.STARTED)
    action.step([], {})
    assert action.records()[1]["phase"] == "stopping", "the horse still holds it: wait"
    job.set_status(rq.job.JobStatus.STOPPED)
    action.step([], {})
    assert action.records()[1]["phase"] == "requeued" and forgotten == [1]
    assert units.job_id(run_id, 0, 1) in queue.job_ids and _job(conn, run_id, 1).get_status(refresh=True) == rq.job.JobStatus.QUEUED
    assert not conn.hexists(runrecord.progress_key(run_id, 0), 1), "the beat of the stopped attempt is cleared"
    assert runrecord.read_run(conn, run_id)["state"] == "running"


def test_a_unit_stuck_on_its_second_attempt_fails_the_run_with_both_stacks(conn, run_id):
    errors = []
    action = _stuck(conn, run_id, errors, stack="FIRST")
    _beat(conn, run_id, 1)
    action.step([1], {1: _quiet()})
    job = _job(conn, run_id, 1)
    redisqueue.queue(redisqueue.BULK_QUEUE, conn).remove(job.id)
    job.set_status(rq.job.JobStatus.STOPPED)
    action.step([], {})
    _beat(conn, run_id, 1)
    action.dump = None
    second = unitstuck.Stuck(conn, run_id, 0, dump=lambda raw: "SECOND", stop=lambda job_id: None, record_error=lambda message, **f: errors.append(message))
    second.step([1], {1: _quiet()})
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "failed" and "stuck twice" in run["error"] and "FIRST" in run["error"] and "SECOND" in run["error"]
    assert run["killed"] and any("SECOND" in line for line in errors)


def test_a_stop_that_does_not_take_effect_fails_the_run(conn, run_id, monkeypatch):
    action = _stuck(conn, run_id)
    _beat(conn, run_id, 0)
    action.step([0], {0: _quiet()})
    job = _job(conn, run_id, 0)
    redisqueue.queue(redisqueue.BULK_QUEUE, conn).remove(job.id)
    job.set_status(rq.job.JobStatus.STARTED)
    record = action.records()[0]
    conn.hset(action.key, 0, json.dumps({**record, "at": record["at"] - unitstuck.STOP_WAIT_SECONDS - 5}))
    action.step([], {})
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "failed" and "could not be stopped" in run["error"]


def test_a_unit_whose_counter_moves_is_never_touched_however_far_over_its_estimate(conn, run_id):
    stopped = []
    action = _stuck(conn, run_id, stopped=stopped)
    now = unitbeat.server_time(conn)
    conn.hset(runrecord.progress_key(run_id, 0), 0, json.dumps({"done": 500.0, "estimate": 10.0, "beat": now, "started": now - 20_000, "moved": now - 1.0, "phase": "p",
                                                                "waiting_on": "", "restarts": 0}))
    verdict = beatjudge.judge(unitbeat.read(conn, run_id, 0)[0], beatjudge.History(), now)
    assert not verdict.actionable
    action.step([], {0: verdict})
    assert stopped == [] and action.records() == {} and runrecord.read_run(conn, run_id)["state"] == "running"


def test_a_unit_that_finished_meanwhile_is_left_alone(conn, run_id):
    stopped = []
    action = _stuck(conn, run_id, stopped=stopped)
    runrecord.unit_done(conn, run_id, 0, 1, seconds=1.0)
    action.step([1], {1: _quiet()})
    assert stopped == []


def test_the_time_limit_is_only_a_backstop():
    assert units.unit_timeout(None) == units.unit_timeout(1) == units.unit_timeout(10 ** 9) == redisqueue.MAX_TIMEOUT_SECONDS


def test_a_frozen_horse_writes_its_stack_when_it_is_sent_the_signal(tmp_path):
    path = str(tmp_path / "stack.txt")
    script = textwrap.dedent(f"""
        import sys, time
        sys.path.insert(0, {os.path.dirname(os.path.dirname(os.path.dirname(unitbeat.__file__)))!r})
        from planetgen.queue import unitbeat
        assert unitbeat.register_dump({path!r})
        def frozen_work_loop():
            print("ready", flush=True)
            time.sleep(60)
        frozen_work_loop()
    """)
    child = subprocess.Popen([sys.executable, "-I", "-c", script], stdout=subprocess.PIPE, text=True)
    try:
        assert child.stdout.readline().strip() == "ready"
        import socket
        text = unitstuck.request_dump({"stack_file": path, "pid": child.pid, "host": socket.gethostname()})
        assert "frozen_work_loop" in text, text
        assert child.poll() is None, "the signal dumps and does not end the process"
        gone = unitstuck.request_dump({"stack_file": path, "pid": 999_999_999, "host": socket.gethostname()})
        assert "gone" in gone
    finally:
        child.kill()
        child.wait()


def test_the_watch_acts_on_a_unit_only_after_two_checks_in_a_row_say_quiet(conn, run_id):
    from planetgen.queue import unitwatch
    stopped = []
    watch = unitwatch.UnitWatch(conn, run_id, 0, "stuckkind", record_error=lambda *a, **k: None, stuck=_stuck(conn, run_id, stopped=stopped))
    _beat(conn, run_id, 1)
    now = unitbeat.server_time(conn)
    conn.hset(runrecord.progress_key(run_id, 0), 1, json.dumps({"done": 3.0, "estimate": 9.0, "beat": now, "started": now - 3000, "moved": now - 2900, "phase": "p",
                                                                "waiting_on": "", "restarts": 0}))
    watch.check()
    assert stopped == [], "one check acts on nothing"
    watch.check()
    assert stopped == [units.job_id(run_id, 0, 1)]


def _stopped_and_requeued(conn, run_id, action, unit=1):
    _beat(conn, run_id, unit)
    action.step([unit], {unit: _quiet()})
    job = _job(conn, run_id, unit)
    redisqueue.queue(redisqueue.BULK_QUEUE, conn).remove(job.id)
    job.set_status(rq.job.JobStatus.STOPPED)
    action.step([], {})
    return job


def test_a_restarted_unit_is_queued_at_the_front_and_the_queue_is_held_for_it_alone(conn, run_id):
    action = _stuck(conn, run_id)
    job = _stopped_and_requeued(conn, run_id, action)
    assert redisqueue.queue(redisqueue.BULK_QUEUE, conn).job_ids[0] == job.id
    state = pause.state(conn)
    assert state["paused"] and state["solo"] == {"job": job.id, "run": run_id, "stage": 0, "unit": 1}
    assert state["reason"] == "solo run of unit 1 (stage 0) after a restart" and state["by"] == f"run {run_id}"
    assert "unit 1 runs alone" in runrecord.log_text(conn, run_id)


def test_the_hold_lifts_when_the_restarted_unit_completes(conn, run_id):
    action = _stuck(conn, run_id)
    _stopped_and_requeued(conn, run_id, action)
    action.step([], {})
    assert pause.is_paused(conn), "the unit has not completed"
    runrecord.unit_done(conn, run_id, 0, 1, seconds=1.0)
    action.step([], {})
    assert not pause.is_paused(conn)


def test_the_hold_lifts_when_the_restarted_unit_fails_the_run(conn, run_id):
    action = _stuck(conn, run_id)
    _stopped_and_requeued(conn, run_id, action)
    _beat(conn, run_id, 1)
    _stuck(conn, run_id).step([1], {1: _quiet()})
    assert runrecord.read_run(conn, run_id)["state"] == "failed" and not pause.is_paused(conn)


def test_the_hold_lifts_when_the_run_is_canceled_meanwhile(conn, run_id):
    action = _stuck(conn, run_id)
    _stopped_and_requeued(conn, run_id, action)
    runrecord.finish_run(conn, run_id, "canceled")
    action.step([], {})
    assert not pause.is_paused(conn)


def test_a_pause_by_hand_is_kept_and_the_unit_waits_for_the_resume(conn, run_id):
    pause.pause("boss", conn)
    action = _stuck(conn, run_id)
    _stopped_and_requeued(conn, run_id, action)
    action.step([], {})
    assert pause.state(conn)["by"] == "boss" and pause.state(conn)["solo"] is None
    runrecord.unit_done(conn, run_id, 0, 1, seconds=1.0)
    action.step([], {})
    assert pause.is_paused(conn), "the follower lifts only the hold it made"


def test_nothing_is_held_when_the_follower_is_told_not_to(conn, run_id):
    action = unitstuck.Stuck(conn, run_id, 0, dump=lambda raw: "S", stop=lambda job_id: None, record_error=lambda message, **f: None, solo=False)
    _stopped_and_requeued(conn, run_id, action)
    assert not pause.is_paused(conn)


def test_a_hold_whose_follower_is_gone_lifts_when_the_unit_job_has_ended(conn, run_id):
    action = _stuck(conn, run_id)
    job = _stopped_and_requeued(conn, run_id, action)
    assert pause.is_paused(conn)
    job.set_status(rq.job.JobStatus.FINISHED)
    assert not pause.is_paused(conn)


def test_a_unit_queued_again_in_the_pass_that_finds_it_stopped_is_not_acted_on_with_the_verdict_of_its_old_attempt(conn, run_id):
    """OPS.143: the unit is still dead or quiet when its stop takes effect, so it is in `agreed` and gone from the verdicts the watch forgot; that is not a crash."""
    stopped = []
    action = _stuck(conn, run_id, stopped=stopped)
    _beat(conn, run_id, 1)
    action.step([1], {1: _quiet()})
    job = _job(conn, run_id, 1)
    redisqueue.queue(redisqueue.BULK_QUEUE, conn).remove(job.id)
    job.set_status(rq.job.JobStatus.STOPPED)
    verdicts = {1: _quiet()}
    action = unitstuck.Stuck(conn, run_id, 0, dump=lambda raw: "STACK", stop=stopped.append, record_error=lambda *a, **k: None, forget=lambda unit: verdicts.pop(unit, None))
    action.step([1], verdicts)         # requeues, forgets the unit, and must not look the old verdict up
    assert action.records()[1]["phase"] == "requeued" and len(stopped) == 1
