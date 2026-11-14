"""OPS.207: the supervisor process (`planetgen/queue/supervisor.py`): one per host by a file lock, a sign of life in Redis, its loop driven by a fake clock, store
and process table, and started on demand by a submit, a status poll and the warm-up."""

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
    assert store.retracted == 1 and "ends after 6 turn(s): idle" in lines[-1]


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
