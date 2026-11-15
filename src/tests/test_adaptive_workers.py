"""Adaptive worker control: the bulk worker limit follows the host's memory (queue/adaptive.py)."""

import os
import uuid

import pytest
import rq

from planetgen.queue import adaptive, memorygate, redisqueue, runrecord, workers
from planetgen.queue.memorygate import Claim, Holder, Host

pytestmark = pytest.mark.real_queue


@pytest.fixture
def lane(redis_server, key_prefix, monkeypatch):
    """A bulk queue of this test's own, the gate's keys under its prefix, workers registered by hand, starts and shutdowns recorded."""
    token = uuid.uuid4().hex[:12]
    monkeypatch.setattr(redisqueue, "BULK_QUEUE", f"planetgen-bulk-test-{token}")
    monkeypatch.setattr(memorygate, "PREFIX", f"{key_prefix}:gate")
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    connection = redisqueue.connect()
    started, told = [], []
    monkeypatch.setattr(workers.redisqueue, "start_detached", lambda argv, cwd: started.append(argv))
    monkeypatch.setattr(workers.command, "send_shutdown_command", lambda conn, name: told.append(name))
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 6)
    monkeypatch.setattr(workers, "host_worker_processes", lambda: 0)
    connection.delete(workers.WANTED_KEY)
    made = []

    def register(state="busy"):
        worker = rq.Worker([redisqueue.queue(redisqueue.BULK_QUEUE, connection)], connection=connection)
        worker.register_birth()
        connection.hset(worker.key, "pid", os.getpid())
        worker.set_state(state)
        made.append(worker)
        return worker

    yield connection, started, told, register
    for worker in made:
        worker.register_death()
        connection.delete(workers.LENT_PREFIX + worker.name)
    redisqueue.queue(redisqueue.BULK_QUEUE, connection).empty()
    connection.delete(workers.STARTING_PREFIX + redisqueue.BULK_QUEUE, workers.DECIDING_KEY)


def host(available=8000.0):
    return Host(16000.0, available, 8, 0.5)


def test_a_memory_wait_under_thirty_seconds_changes_nothing(lane):
    connection, started, told, register = lane
    [register() for _ in range(4)]
    assert adaptive.on_memory_wait(connection, "w", "unit:x", 29.0) is None
    assert adaptive.limit_of(connection) is None and not told


def test_a_memory_wait_over_thirty_seconds_lowers_the_limit_by_one_and_asks_a_worker_to_leave(lane):
    connection, started, told, register = lane
    made = [register() for _ in range(4)]
    assert adaptive.on_memory_wait(connection, made[0].name, "unit:x", 31.0, now=1000.0) == 3
    assert adaptive.limit_of(connection) == 3
    assert told == [made[0].name]                                      # the one that waits at the gate goes first
    recent = adaptive.recent(connection)
    assert recent[0]["from"] == 4 and recent[0]["to"] == 3 and "31 s" in recent[0]["why"]


def test_two_lowerings_are_a_cooldown_apart(lane):
    connection, started, told, register = lane
    made = [register() for _ in range(5)]
    assert adaptive.on_memory_wait(connection, made[0].name, "unit:x", 40.0, now=1000.0) == 4
    assert adaptive.on_memory_wait(connection, made[1].name, "unit:x", 45.0, now=1005.0) is None
    assert adaptive.on_memory_wait(connection, made[1].name, "unit:x", 70.0, now=1000.0 + adaptive.LOWER_COOLDOWN_SECONDS) == 3


def test_the_limit_never_goes_below_one(lane):
    connection, started, told, register = lane
    only = register()
    assert adaptive.on_memory_wait(connection, only.name, "unit:x", 500.0, now=1000.0) is None
    assert adaptive.limit_of(connection) is None and not told


def test_ensure_workers_starts_no_more_than_the_limit(lane):
    connection, started, told, register = lane
    [register() for _ in range(4)]
    for _ in range(10):
        redisqueue.queue(redisqueue.BULK_QUEUE, connection).enqueue(os.getpid, job_timeout=600)
    runrecord.start_run(connection, "Fill", [runrecord.Stage("mass", "Mass", units=200)])
    connection.hset(adaptive._key("state"), mapping={"limit": 4, "lowered_at": adaptive.time.time()})
    assert workers.ensure_workers(redisqueue.BULK_QUEUE, connection) == 0         # six wanted, but the limit holds it at the four it has
    connection.delete(adaptive._key("state"))
    assert workers.ensure_workers(redisqueue.BULK_QUEUE, connection) == 2


def test_a_host_under_pressure_lowers_the_limit_before_any_job_has_waited(lane):
    connection, started, told, register = lane
    [register() for _ in range(4)]
    low = host(available=memorygate.reserve_mb(host()) * 0.5)
    assert adaptive.adapt(connection, 6, 3, host=low, held=[], now=1000.0) == 3
    assert adaptive.adapt(connection, 6, 3, host=host(), held=[], now=1001.0) == 3        # not raised at once


def test_the_limit_rises_when_memory_has_been_free_and_jobs_are_queued(lane):
    connection, started, told, register = lane
    [register() for _ in range(3)]
    connection.hset(adaptive._key("state"), mapping={"limit": 3, "lowered_at": 1000.0})
    assert adaptive.adapt(connection, 6, 9, host=host(), held=[], now=1000.0 + adaptive.RAISE_QUIET_SECONDS - 1) == 3     # still quiet for too short a time
    assert adaptive.adapt(connection, 6, 9, host=host(), held=[], now=1000.0 + adaptive.RAISE_QUIET_SECONDS) == 4
    assert adaptive.adapt(connection, 6, 9, host=host(), held=[], now=1001.0 + adaptive.RAISE_QUIET_SECONDS) == 4         # one step per RAISE_STEP_SECONDS
    assert adaptive.recent(connection)[0]["to"] == 4


def test_the_limit_does_not_rise_without_room_or_without_queued_jobs(lane):
    connection, started, told, register = lane
    [register() for _ in range(3)]
    connection.hset(adaptive._key("state"), mapping={"limit": 3, "lowered_at": 0.0})
    claim = Claim("unit:x", 2100.0, 1.0, memorygate.LANE_BULK)
    full = host(available=memorygate.reserve_mb(host()) + 100.0)
    assert adaptive.adapt(connection, 6, 9, host=full, held=[], now=1000.0) == 3                      # no room for another worker
    assert adaptive.adapt(connection, 6, 0, host=host(), held=[Holder("a", claim, 1, "h", 0.0, 0.0)][:0], now=1000.0) == 3      # nothing queued that lacks a worker


def test_the_limit_is_dropped_once_it_reaches_the_cap(lane):
    connection, started, told, register = lane
    [register() for _ in range(5)]
    connection.hset(adaptive._key("state"), mapping={"limit": 5, "lowered_at": 0.0})
    assert adaptive.adapt(connection, 6, 9, host=host(), held=[], now=1000.0) == 6
    assert adaptive.limit_of(connection) is None


def test_the_control_is_off_with_the_gate(lane, monkeypatch):
    connection, started, told, register = lane
    monkeypatch.setenv(memorygate.DISABLE_ENV, "0")
    assert adaptive.adapt(connection, 6, 9, host=host(available=10.0), held=[], now=1000.0) == 6


def test_a_gate_wait_for_memory_reaches_the_control_and_a_wait_for_a_core_does_not(lane, monkeypatch):
    connection, started, told, register = lane
    me = register()
    [register() for _ in range(3)]
    seen = []
    monkeypatch.setattr(adaptive, "on_memory_wait", lambda conn, name, kind, waited, now=None: seen.append((name, kind, waited)))
    claim = Claim("unit:x", 2100.0, 1.0, memorygate.LANE_BULK)
    other = Holder("o", Claim("unit:y", 7000.0, 1.0, memorygate.LANE_BULK), os.getpid(), memorygate.socket.gethostname(), 7000.0, 7000.0)
    memorygate._write(connection, "o", other.claim, os.getpid(), "j", memorygate.time.time())
    monkeypatch.setattr(memorygate, "read_host", lambda: host(available=1500.0))
    monkeypatch.setattr(memorygate, "held_mb", lambda pid: 7000.0)
    now = [0.0]
    ticks = []

    def sleep(seconds):
        now[0] += 15.0
        ticks.append(now[0])
        if len(ticks) >= 3:
            raise StopIteration

    with pytest.raises(StopIteration):
        memorygate.acquire(connection, claim, me.name, clock=lambda: now[0], sleep=sleep)
    assert seen and seen[-1][1] == "unit:x" and seen[-1][2] >= 30.0
    seen.clear()
    monkeypatch.setattr(memorygate, "read_host", lambda: host(available=14000.0, ) if False else Host(16000.0, 14000.0, 8, 3.0))      # memory is free; the load is high
    now[0] = 0.0
    ticks.clear()
    with pytest.raises(StopIteration):
        memorygate.acquire(connection, claim, me.name, clock=lambda: now[0], sleep=sleep)
    assert not seen
