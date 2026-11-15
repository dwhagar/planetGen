"""PERF.153, PERF.154: a run that finds the machine taken waits in the one waiting list (`planetgen/queue/runline.py`),
the list is refused when it holds `queue.max_waiting` entries, and an API job shows that it waits
(`wrapper.report_wait`, `api_jobs.status`)."""

import os
import socket
import time
from dataclasses import asdict

import rq
import pytest

from planetgen.queue import api_jobs, commands, redisqueue, runline, runrecord, supervisor, wrapper
from planetgen.queue.runrecord import Stage
from planetgen.util.settings import get_settings


@pytest.fixture
def conn(redis_server, key_prefix, monkeypatch):
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    return redisqueue.connect(redis_server)


def _start(title="Waiting"):
    """What a run asks the supervisor to start (OPS.211): the arguments of `runrecord.start_run`, as `wrapper.run` leaves them."""
    return {"title": title, "stages": [asdict(Stage("mass", "Mass"))]}


def _ask(conn):
    """Stands for the supervisor: one turn of its duty that starts runs, at every poll of the run that asked."""
    return lambda: supervisor.RunStarts(conn).turn()


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
    run_id = runline.start_when_free(conn, _start(), 3600, _ask(conn), on_wait=lambda *a: shown.append(a))
    assert runrecord.active_run(conn) == run_id and shown == [] and runline.entries(conn) == []


def test_a_run_waits_for_the_one_that_holds_the_machine_and_then_takes_it(conn):
    holder = runrecord.start_run(conn, "Holder", [Stage("mass", "Mass")])
    shown = []
    clock = Clock(lambda n: runrecord.finish_run(conn, holder, "done") if n == 3 else None)
    run_id = runline.start_when_free(conn, _start(), 3600, _ask(conn), on_wait=lambda *a: shown.append(a),
                                     clock=clock.time, sleep=clock.sleep)
    assert runrecord.active_run(conn) == run_id != holder
    assert shown == [(holder, 1), (None, None)] and runline.entries(conn) == []


def test_the_wait_ends_after_the_bound_and_names_the_run(conn):
    holder = runrecord.start_run(conn, "Holder", [Stage("mass", "Mass")])
    clock = Clock()
    with pytest.raises(runrecord.Busy) as refused:
        runline.start_when_free(conn, _start(), 120, _ask(conn), clock=clock.time, sleep=clock.sleep)
    assert holder in str(refused.value) and "waited 2 minutes" in str(refused.value)
    assert runline.entries(conn) == [] and runrecord.active_run(conn) == holder


def test_a_wait_of_zero_does_not_wait(conn):
    runrecord.start_run(conn, "Holder", [Stage("mass", "Mass")])
    clock = Clock()
    with pytest.raises(runrecord.Busy):
        runline.start_when_free(conn, _start(), 0, _ask(conn), clock=clock.time, sleep=clock.sleep)
    assert clock.sleeps == 0


def test_the_line_is_first_come_first_served_and_shows_the_place(conn):
    holder = runrecord.start_run(conn, "Holder", [Stage("mass", "Mass")])
    ahead = runline.join(conn, runline.RUN, "Ahead")
    shown = []

    def move(n):
        if n == 2:                                   # the run ahead gives up
            runline.leave(conn, ahead)
        if n == 4:
            runrecord.finish_run(conn, holder, "done")

    clock = Clock(move)
    run_id = runline.start_when_free(conn, _start(), 3600, _ask(conn), on_wait=lambda *a: shown.append(a),
                                     clock=clock.time, sleep=clock.sleep)
    assert [position for _holder, position in shown[:-1]] == [2, 1] and runrecord.active_run(conn) == run_id


def test_a_waiter_that_was_killed_drops_out_of_the_line(conn):
    ticket = runline.join(conn, runline.RUN, "Dead")
    conn.delete(runline._alive_key(ticket))                       # no sign of life
    assert runline.entries(conn) == []
    run_id = runline.start_when_free(conn, _start(), 60, _ask(conn))        # not held up by it
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


def test_a_ninth_entry_is_refused_with_the_line_named(conn):
    holder = runrecord.start_run(conn, "Holder", [Stage("mass", "Mass")])
    for number in range(8):
        runline.join(conn, runline.RUN, f"Run {number}", who="boss")
    clock = Clock()
    with pytest.raises(runline.LineFull) as refused:
        runline.start_when_free(conn, _start(), 3600, _ask(conn), title="Ninth", clock=clock.time, sleep=clock.sleep)
    message = str(refused.value)
    assert "8 are waiting" in message and '"Run 0" (boss)' in message and '"Run 7" (boss)' in message and "holds 8" in message
    assert clock.sleeps == 0 and len(runline.entries(conn)) == 8 and runrecord.active_run(conn) == holder
    assert isinstance(refused.value, runrecord.Busy), "the commands that refuse to wait answer it as they answer any busy machine"


def test_the_cap_is_a_setting_and_counts_jobs_and_runs_together(conn, monkeypatch):
    monkeypatch.setattr(get_settings().queue, "max_waiting", 2)
    runline.join(conn, runline.JOB, "A job", who="boss", fields={"root": "/jobs"})
    runline.join(conn, runline.RUN, "A run")
    with pytest.raises(runline.LineFull):
        runline.join(conn, runline.RUN, "Another run")


def test_a_run_a_job_started_gets_a_place_in_a_full_line(conn, monkeypatch):
    monkeypatch.setattr(get_settings().queue, "max_waiting", 1)
    runline.join(conn, runline.JOB, "A job", fields={"root": "/jobs"})
    holder = runrecord.start_run(conn, "Holder", [Stage("mass", "Mass")])
    clock = Clock(lambda n: runrecord.finish_run(conn, holder, "done") if n == 2 else None)
    run_id = runline.start_when_free(conn, _start(), 3600, _ask(conn), admitted=True, clock=clock.time, sleep=clock.sleep)
    assert runrecord.active_run(conn) == run_id


def test_jobs_and_runs_each_keep_their_own_order_and_the_list_shows_both_as_they_arrived(conn):
    first_job = runline.join(conn, runline.JOB, "Job 1", fields={"root": "/jobs"})
    first_run = runline.join(conn, runline.RUN, "Run 1")
    second_job = runline.join(conn, runline.JOB, "Job 2", fields={"root": "/jobs"})
    second_run = runline.join(conn, runline.RUN, "Run 2")
    assert [(e["ticket"], e["type"], e["position"]) for e in runline.entries(conn)] == [
        (first_job, "job", 1), (first_run, "run", 1), (second_job, "job", 2), (second_run, "run", 2)]


def test_a_run_is_not_held_up_by_a_job_waiting_ahead_of_it(conn):
    """A waiting job waits for the running job; the stage of that job must not wait behind it (no deadlock)."""
    runline.join(conn, runline.JOB, "Waiting job", fields={"root": "/jobs"})
    run_id = runline.start_when_free(conn, _start(), 60, _ask(conn))
    assert runrecord.active_run(conn) == run_id


def test_a_run_started_inside_a_run_of_the_same_process_fails_instead_of_waiting_for_itself(conn):
    """ADM.136: the second run would wait at position 1 behind the first, which its own process is following."""
    runrecord.start_run(conn, "Generate job", [runrecord.Stage("a", "A")])
    with pytest.raises(runrecord.Busy, match="would wait for itself"):
        runline.start_when_free(conn, _start("Nested"), 600, _ask(conn), title="Nested", nested_guard=True)


def _left_behind(conn, pid=None):
    """A run of a process that is still running but no longer follows it: its owner is a live process id (this test's parent, a
    stand-in for a worker that lives on), and nothing says it is followed (its `alive` key has run out)."""
    holder = runrecord.start_run(conn, "A New galaxy job's earlier run", [runrecord.Stage("a", "A")])
    conn.hset(runrecord.run_key(holder), "owner", f"{socket.gethostname()}:{pid or os.getppid()}")
    conn.delete(runrecord.alive_key(holder))
    return holder


def test_a_run_left_by_a_process_that_lives_on_is_freed_and_the_next_run_starts(conn):
    """ADM.136: the holder's process id is running (a long-lived worker, or a reused id), so the pid check alone called it live and a New
    galaxy job waited on it for ever. Nothing says it is followed (its `alive` key ran out), so it is dead."""
    holder = _left_behind(conn)
    assert runrecord._pid_alive(os.getppid()) and runrecord.holder_is_dead(conn, holder)
    clock = Clock()
    run_id = runline.start_when_free(conn, _start(), 3600, _ask(conn), clock=clock.time, sleep=clock.sleep, nested_guard=True)
    assert runrecord.active_run(conn) == run_id != holder
    left = runrecord.read_run(conn, holder)
    assert left["state"] == "failed" and "process was gone" in left["error"]


def test_a_run_that_is_followed_is_not_freed_however_long_it_takes(conn):
    """The `alive` key is there (its follower says so), and the follower is a live process of this host: the run waits for it."""
    holder = _left_behind(conn)
    runrecord.touch(conn, holder)
    assert not runrecord.holder_is_dead(conn, holder)
    clock = Clock()
    with pytest.raises(runrecord.Busy, match="waited"):
        runline.start_when_free(conn, _start(), 600, _ask(conn), clock=clock.time, sleep=clock.sleep)
    assert runrecord.read_run(conn, holder)["state"] == "running"


def test_a_run_of_a_process_id_that_is_gone_is_still_freed_at_once(conn):
    holder = runrecord.start_run(conn, "Killed", [runrecord.Stage("a", "A")])
    conn.hset(runrecord.run_key(holder), "owner", f"{socket.gethostname()}:2147483646")
    runrecord.touch(conn, holder)       # its key has not run out yet: the missing process is enough
    assert runrecord.holder_is_dead(conn, holder)


def test_the_heartbeat_keeps_a_busy_follower_alive_and_stops_with_it(conn):
    """A follower busy in a callback does not poll; the heartbeat says the run is followed meanwhile, and stops with the follower."""
    run_id = runrecord.start_run(conn, "Busy", [runrecord.Stage("a", "A")])
    conn.delete(runrecord.alive_key(run_id))
    with runrecord.Heartbeat(conn, run_id, interval=0.02):
        deadline = time.monotonic() + 5
        while not conn.exists(runrecord.alive_key(run_id)) and time.monotonic() < deadline:
            time.sleep(0.02)
        assert conn.exists(runrecord.alive_key(run_id)) and not runrecord.holder_is_dead(conn, run_id)
    conn.delete(runrecord.alive_key(run_id))
    time.sleep(0.2)
    assert not conn.exists(runrecord.alive_key(run_id)), "the heartbeat is over with the block"


def test_a_run_whose_set_up_fails_is_ended_not_left_holding_the_machine(conn, monkeypatch):
    """ADM.136: the set-up after the machine is taken (workers, relay, log) was outside the `try`, so a failure there left the run live
    for a process that went on; the next New galaxy waited for it."""
    kind = next(iter(commands.REGISTRY))

    def refuse(*args, **kwargs):
        raise RuntimeError("no worker could be started")

    monkeypatch.setattr(wrapper, "stage_ensure", refuse)
    with pytest.raises(RuntimeError, match="no worker could be started"):
        wrapper.run("Set-up fails", [commands.Step("a", kind, "A", {})], lambda index, plan: None, lambda index, unit, result: None)
    assert runrecord.active_run(conn) is None
    assert runrecord.read_run(conn, runrecord.recent_runs(conn, limit=1)[0])["state"] == "failed"
