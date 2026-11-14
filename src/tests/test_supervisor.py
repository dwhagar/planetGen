"""OPS.207: the supervisor process (`planetgen/queue/supervisor.py`): one per host by a file lock, a sign of life in Redis, its loop driven by a fake clock, store
and process table, and started on demand by a submit, a status poll and the warm-up."""

import socket
import subprocess
import sys

import pytest

from planetgen.queue import redisqueue, runrecord, supervisor, warmstart, workers


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


class FakeClock:
    def __init__(self):
        self.t = 1000.0

    def now(self):
        return self.t

    def sleep(self, seconds):
        self.t += seconds


class FakeStore:
    def __init__(self, busy=lambda now: False, clock=None):
        self.published, self.retracted, self._busy, self.clock = [], 0, busy, clock

    def publish(self, pid, since, now, turns):
        self.published.append((pid, since, now, turns))

    def retract(self):
        self.retracted += 1

    def busy(self):
        return self._busy(self.clock.t)


def _supervisor(clock, store, **kwargs):
    lines = []
    return supervisor.Supervisor(clock, store, supervisor.ProcessTable(), pid=77, lines=lines.append, **kwargs), lines


def test_an_idle_supervisor_publishes_each_turn_and_ends_after_its_idle_spell():
    clock = FakeClock()
    store = FakeStore(clock=clock)
    loop, lines = _supervisor(clock, store, idle_seconds=10, turn_seconds=2)
    assert loop.run() == "idle"
    assert [turn[3] for turn in store.published] == [1, 2, 3, 4, 5, 6], "a turn every 2 s, the sixth at 10 s idle"
    assert {turn[1] for turn in store.published} == {1000.0}, "`since` stays the time of the first turn"
    assert store.retracted == 1 and "ended (nothing to watch for 10 s, after 6 turn(s))" in lines[-1]


def test_something_to_watch_keeps_it_alive_and_the_idle_spell_counts_from_when_it_stopped():
    clock = FakeClock()
    store = FakeStore(busy=lambda now: now < 1040, clock=clock)
    loop, _ = _supervisor(clock, store, idle_seconds=10, turn_seconds=2)
    assert loop.run() == "idle"
    assert clock.t >= 1050 and clock.t < 1054, "it ended about 10 s after the work ended at 1040"


def test_a_stop_request_ends_it_at_once_and_max_turns_bounds_a_test():
    clock = FakeClock()
    loop, _ = _supervisor(clock, FakeStore(busy=lambda now: True, clock=clock), stop=lambda: clock.t >= 1006)
    assert loop.run() == "stopped"
    clock2 = FakeClock()
    loop, _ = _supervisor(clock2, FakeStore(busy=lambda now: True, clock=clock2))
    assert loop.run(max_turns=3) == "turns" and loop.turns == 3


def test_only_one_process_holds_the_host_lock_and_the_kernel_frees_it_with_the_process(tmp_path):
    path = str(tmp_path / "supervisor.lock")
    first = supervisor.take_lock(path)
    assert first is not None
    assert supervisor.take_lock(path) is None, "a second holder is refused while the first lives"
    first.close()
    second = supervisor.take_lock(path)
    assert second is not None
    second.close()
    code = ("import sys; from planetgen.queue import supervisor; h = supervisor.take_lock(sys.argv[1]); print('held' if h else 'refused'); sys.stdout.flush(); "
            "import os; os.kill(os.getpid(), 9)")
    out = subprocess.run([sys.executable, "-c", code, path], capture_output=True, text=True, check=False).stdout.strip()
    assert out == "held"
    assert supervisor.take_lock(path) is not None, "a supervisor killed with -9 does not block the next"


def test_the_redis_store_publishes_a_sign_of_life_that_expires_and_is_taken_down(conn):
    store = supervisor.RedisStore(conn)
    assert not supervisor.is_alive(conn)
    store.publish(77, 1000.0, 1002.0, 2)
    assert supervisor.is_alive(conn)
    fields = {runrecord._s(k): runrecord._s(v) for k, v in conn.hgetall(supervisor.alive_key()).items()}
    assert fields == {"pid": "77", "since": "1000.0", "turn_at": "1002.0", "turns": "2"}
    assert 0 < conn.ttl(supervisor.alive_key()) <= supervisor.ALIVE_SECONDS
    store.retract()
    assert not supervisor.is_alive(conn)


def test_the_redis_store_sees_a_waiting_run_as_something_to_watch(conn):
    store = supervisor.RedisStore(conn)
    assert not store.busy()
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("a", "A", units=1)])
    assert store.busy()
    runrecord.finish_run(conn, run_id, "done")
    assert not store.busy()


def test_a_caller_starts_one_supervisor_and_the_callers_after_it_do_not(conn):
    started = []
    assert supervisor.ensure_started(conn, start=started.append) is True
    assert supervisor.ensure_started(conn, start=started.append) is False, "the first caller's claim holds off the second"
    assert len(started) == 1 and started[0][1:3] == ["-m", "planetgen.cli.supervisor"]
    conn.delete(supervisor.starting_key())
    supervisor.RedisStore(conn).publish(77, 1.0, 2.0, 1)
    assert supervisor.ensure_started(conn, start=started.append) is False, "one is alive: nothing to start"
    assert len(started) == 1


def test_a_supervisor_that_could_not_be_started_never_fails_the_caller(conn):
    def refuse(argv):
        raise OSError("no fork")
    assert supervisor.ensure_started(conn, start=refuse) is False


@pytest.mark.real_queue
def test_a_submit_a_status_poll_and_the_warm_up_each_call_for_a_supervisor(conn, monkeypatch):
    calls = []
    monkeypatch.setattr(supervisor, "ensure_started", lambda *args, **kwargs: calls.append(1) or False)
    monkeypatch.setattr(workers.updateguard, "is_updating", lambda: False)
    workers.ensure_workers(redisqueue.BULK_QUEUE, conn)
    assert calls, "ensure_workers (every submit)"
    calls.clear()
    monkeypatch.setattr(warmstart, "enabled", lambda: True)
    monkeypatch.setattr(warmstart, "_client", lambda: conn)
    monkeypatch.setattr(warmstart, "acquire", lambda database, client=None: False)
    assert warmstart.request("galaxy", "sector") == "running"
    assert calls, "the warm-up's start"


def test_the_supervisors_acts_are_kept_newest_first_with_their_reasons_and_shown_with_its_sign_of_life(conn):
    assert supervisor.overview(conn) == {"alive": False, "recent": []}
    store = supervisor.RedisStore(conn)
    clock = FakeClock()
    loop = supervisor.Supervisor(clock, store, supervisor.ProcessTable(), pid=77, lines=lambda text: None, idle_seconds=4, turn_seconds=2)
    loop.run()
    assert [(act["what"], act["at"]) for act in supervisor.overview(conn)["recent"]] == [("ended", 1004.0), ("started", 1000.0)]
    assert supervisor.overview(conn)["recent"][1]["reason"] == "something asked for one"
    assert supervisor.overview(conn)["alive"] is False, "it retracted its sign of life when it ended"
    store.publish(77, 1000.0, 1002.0, 2)
    assert supervisor.overview(conn)["alive"] and supervisor.overview(conn)["turns"] == 2 and supervisor.overview(conn)["since"] == 1000.0


def test_the_list_of_acts_keeps_only_the_newest_twenty(conn):
    store = supervisor.RedisStore(conn)
    for number in range(30):
        store.note(float(number), f"act {number}")
    recent = supervisor.overview(conn)["recent"]
    assert len(recent) == supervisor.RECENT_KEPT and recent[0]["what"] == "act 29"


class FakeStart:
    def __init__(self, fail_on=None):
        self.calls, self.fail_on = [], fail_on

    def __call__(self, queue_name, connection, url, mysql_host, count):
        if queue_name == self.fail_on:
            raise OSError("no fork")
        self.calls.append((queue_name, url, mysql_host, count))
        return 2


def test_an_ask_for_workers_is_left_for_the_supervisor_which_acts_on_it_once_and_the_later_ask_wins(conn):
    supervisor.nudge(conn, "bulk", url="redis://x", count=4)
    supervisor.nudge(conn, "bulk", url="redis://x", count=6)
    supervisor.nudge(conn, "maintenance", url="redis://x")
    assert conn.hlen(supervisor.nudges_key()) == 2
    assert supervisor.RedisStore(conn).busy(), "an ask not acted on yet is something to do: the supervisor does not end"
    start = FakeStart()
    starts = supervisor.WorkerStarts(conn, start=start)
    assert starts.turn() == 4
    assert sorted(start.calls) == [("bulk", "redis://x", None, 6), ("maintenance", "redis://x", None, None)]
    assert starts.turn() == 0 and conn.hlen(supervisor.nudges_key()) == 0, "taken and cleared: not acted on twice"


def test_one_queues_failed_start_leaves_the_others_and_never_stops_the_turn(conn):
    supervisor.nudge(conn, "bulk")
    supervisor.nudge(conn, "maintenance")
    start = FakeStart(fail_on="bulk")
    assert supervisor.WorkerStarts(conn, start=start).turn() == 2
    assert [call[0] for call in start.calls] == ["maintenance"]


@pytest.mark.real_queue
def test_ensure_workers_only_leaves_an_ask_and_start_workers_is_what_the_supervisor_calls(conn, monkeypatch):
    called = []
    monkeypatch.setattr(workers, "start_workers", lambda *args, **kwargs: called.append(args) or 3)
    assert workers.ensure_workers("planetgen-bulk", conn, count=2) == 0
    assert called == [], "a caller starts nothing itself"
    assert conn.hlen(supervisor.nudges_key()) == 1
    assert supervisor.WorkerStarts(conn).turn() == 3 and called and called[0][0] == "planetgen-bulk"
    conn.delete(supervisor.nudges_key())


def test_the_supervisor_runs_each_of_its_duties_every_turn_and_one_failing_does_not_stop_the_others():
    clock = FakeClock()
    seen = []

    class Duty:
        def __init__(self, name, fail=False):
            self.name, self.fail = name, fail

        def turn(self):
            seen.append(self.name)
            if self.fail:
                raise ConnectionError("down")

        def clear(self):
            seen.append("clear " + self.name)

    loop, _ = _supervisor(clock, FakeStore(busy=lambda now: True, clock=clock), duties=(Duty("a", fail=True), Duty("b")))
    loop.run(max_turns=2)
    assert seen == ["a", "b", "a", "b", "clear a", "clear b"]


def _finished_stage(conn, ended_at):
    """A run whose stage 0 was completed by its only unit at `ended_at` and whose end never queued the next opening (the unit's process was killed in it)."""
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("a", "A", units=1), runrecord.Stage("b", "B", units=1)])
    runrecord.begin_stage(conn, run_id, 0)
    runrecord.set_stage_plan(conn, run_id, 0, 1, 1.0, {"units": 1, "detail": {}})
    assert runrecord.unit_done(conn, run_id, 0, 0, now=ended_at).advance
    return run_id


def test_the_reconciler_ends_a_stage_that_finished_and_never_advanced_after_its_grace_and_says_so(conn):
    run_id = _finished_stage(conn, ended_at=1000.0)
    repaired, acts = [], []
    now = [1010.0]

    def repair(connection, rid):
        repaired.append(rid)
        return [0]

    reconciler = supervisor.Reconciler(conn, clock=lambda: now[0], repair=repair)
    loop = supervisor.Supervisor(FakeClock(), supervisor.RedisStore(conn), supervisor.ProcessTable(), pid=77, lines=acts.append, duties=(reconciler,))
    assert reconciler.turn() == [] and repaired == [], "10 s after the last unit: the unit's own process is still at it"
    now[0] = 1031.0
    assert reconciler.turn() == [(run_id, 0)] and repaired == [run_id]
    assert any("ended stage 1" in line and "had not queued the next stage" in line for line in acts)
    assert loop.duties == (reconciler,)


def test_the_reconciler_leaves_a_stage_a_process_is_ending_now_and_marks_every_live_run_watched(conn):
    run_id = _finished_stage(conn, ended_at=1000.0)
    conn.set(runrecord.ending_key(run_id, 0), 1)
    repaired = []
    reconciler = supervisor.Reconciler(conn, clock=lambda: 9999.0, repair=lambda c, r: repaired.append(r) or [0])
    assert reconciler.turn() == [] and repaired == []
    assert conn.exists(runrecord.supervised_key(run_id)), "the run has a follower that is not its starter"


def test_a_run_the_supervisor_watches_is_not_taken_over_when_its_starting_process_is_gone(conn):
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("a", "A", units=1)])
    conn.hset(runrecord.run_key(run_id), "owner", f"{socket.gethostname()}:999999999")          # a process id that is not running
    conn.delete(runrecord.alive_key(run_id))
    assert runrecord.holder_is_dead(conn, run_id), "without the supervisor the run is taken for a dead one, as before"
    conn.set(runrecord.supervised_key(run_id), "1", ex=15)
    assert not runrecord.holder_is_dead(conn, run_id), "OPS.210: with it, the run carries on"
    runrecord.finish_run(conn, run_id, "done")
    assert runrecord.holder_is_dead(conn, run_id), "a run that is over is never held"
