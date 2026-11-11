"""OPS.145: after a stuck unit is stopped and queued again, the other tasks are held so the unit runs alone: the pause flag carries the hold and names it, the workers take no
other job, the gate gives the unit the host, the scheduler does not charge its time to its kind, and the hold lifts by itself."""

import json
import os
import subprocess
import sys
import time
import types
import uuid

import pytest
import rq

from planetgen.queue import api_jobs, memorygate, pause, redisqueue, runrecord, scheduler
from planetgen.queue.memorygate import Claim, Holder, Host


@pytest.fixture
def alive_job(monkeypatch):
    """The unit's job reads as queued (the tests name jobs that were never made)."""
    monkeypatch.setattr(pause, "_solo_is_over", lambda connection, solo: False)


@pytest.fixture
def conn(redis_server, key_prefix, monkeypatch):
    monkeypatch.setattr(memorygate, "PREFIX", f"{key_prefix}:gate")
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    return redisqueue.connect(redis_server)


def host(available=8000.0, total=16000.0, cores=8, load=0.5):
    return Host(total, available, cores, load)


def holder(name, mb, held=0.0, lane="bulk", kind="unit:x", solo=False):
    return Holder(name, Claim(kind, mb, 1.0, lane, solo=solo), 1, "h", held, held)


def solo_claim(kind="unit:scatter"):
    return Claim(kind, 500.0, 1.0, "bulk", solo=True)


# --- the hold is the pause flag ------------------------------------------------------------------------------

def test_a_solo_hold_is_the_pause_flag_with_the_reason_and_the_unit(conn, alive_job):
    assert pause.hold_solo("run-1", 2, 7, "job-7", conn)
    state = pause.state(conn)
    assert state["paused"] and state["by"] == "run run-1" and state["at"] is not None
    assert state["reason"] == "solo run of unit 7 (stage 2) after a restart"
    assert state["solo"] == {"job": "job-7", "run": "run-1", "stage": 2, "unit": 7}


def test_a_pause_by_hand_is_not_replaced_by_a_solo_hold_and_a_hold_not_by_a_pause(conn, alive_job):
    pause.pause("boss", conn)
    assert not pause.hold_solo("run-1", 0, 1, "job-1", conn)
    assert pause.state(conn)["by"] == "boss" and pause.state(conn)["solo"] is None
    pause.resume(conn)
    assert pause.hold_solo("run-1", 0, 1, "job-1", conn)
    assert not pause.pause("boss", conn) and not pause.hold_solo("run-2", 0, 2, "job-2", conn), "the first hold is the one kept"
    assert pause.solo_of(conn)["unit"] == 1


def test_a_solo_hold_lets_only_its_unit_through_and_a_pause_by_hand_lets_none(conn, monkeypatch):
    monkeypatch.setattr(pause, "_solo_is_over", lambda connection, solo: False)
    assert not pause.blocks(conn) and not pause.blocks(conn, "any")
    pause.hold_solo("run-1", 0, 1, "job-1", conn)
    assert not pause.blocks(conn), "a worker with no job yet must be able to reach the unit in the queue"
    assert not pause.blocks(conn, "job-1") and pause.blocks(conn, "job-2")
    pause.resume(conn)
    pause.pause("boss", conn)
    assert pause.blocks(conn) and pause.blocks(conn, "job-1")


def test_the_follower_lifts_only_the_hold_it_made(conn, monkeypatch):
    monkeypatch.setattr(pause, "_solo_is_over", lambda connection, solo: False)
    pause.hold_solo("run-1", 0, 1, "job-1", conn)
    assert not pause.lift_solo("run-1", 0, 2, conn) and not pause.lift_solo("run-2", 0, 1, conn) and pause.is_paused(conn)
    assert pause.lift_solo("run-1", 0, 1, conn) and not pause.is_paused(conn)
    pause.pause("boss", conn)
    assert not pause.lift_solo("run-1", 0, 1, conn) and pause.is_paused(conn), "a pause by hand is Resume's"


@pytest.mark.parametrize("status,stays", [("queued", True), ("started", True), ("finished", False), ("failed", False), ("canceled", False), ("stopped", False), (None, False)])
def test_the_hold_lifts_by_itself_when_its_job_has_ended_or_is_gone(conn, monkeypatch, status, stays):
    job = None if status is None else types.SimpleNamespace(get_status=lambda refresh=False: status)
    monkeypatch.setattr(redisqueue, "fetch_job", lambda job_id, connection=None, serializer=None: job)
    pause.hold_solo("run-1", 0, 1, "job-1", conn)
    assert pause.is_paused(conn) is stays
    assert (conn.get(pause.key()) is not None) is stays


# --- the gate ---------------------------------------------------------------------------------------------------

def test_the_solo_claim_asks_for_the_whole_budget_less_the_hosts_and_the_sites_reserve():
    asked = memorygate.with_solo_mb(solo_claim(), host(total=16000.0), places=0)
    assert asked.mb == pytest.approx(16000.0 - memorygate.reserve_mb(host()))
    assert asked.solo and memorygate.with_solo_mb(Claim("unit:x", 500.0), host()).mb == 500.0
    kept = memorygate.solo_mb(host(), places=2)
    assert kept == pytest.approx(16000.0 - memorygate.reserve_mb(host()) - 2 * memorygate.KEEP_MB_PER_PLACE)


def test_a_solo_claim_waits_for_every_bulk_or_maintenance_claim_and_not_for_the_sites():
    claim = memorygate.with_solo_mb(solo_claim(), host(), places=0)
    refused = memorygate.decide(claim, host(), [holder("w1", 500.0), holder("w2", 500.0, lane="maintenance", kind="unit:warm")], places=0)
    assert not refused.admit and "solo run waits" in refused.reason and "unit:warm" in refused.reason
    assert memorygate.decide(claim, host(available=2000.0), [holder("v", 700.0, lane="interactive", kind="api:x")], places=0).admit
    assert memorygate.decide(claim, host(), [], places=0).admit


def test_while_a_solo_claim_is_held_no_bulk_claim_is_admitted_and_the_site_still_is():
    claim = memorygate.with_solo_mb(solo_claim(), host(), places=0)
    solo = holder("solo", claim.mb, held=300.0, kind="unit:scatter", solo=True)
    other = memorygate.decide(Claim("unit:y", 200.0, 1.0, "bulk"), host(), [solo], places=0)
    assert not other.admit and "solo run holds the host" in other.reason
    assert memorygate.decide(Claim("api:x", 300.0, 1.0, "interactive"), host(), [solo], places=0).admit, "the visitor's page is not held by the unit"


def test_the_claim_survives_the_claims_table(conn):
    memorygate.acquire(conn, solo_claim(), "w1", pid=os.getpid())
    (held,) = memorygate.holders(conn)
    assert held.claim.solo and held.claim.mb > 8000.0


# --- the scheduler ------------------------------------------------------------------------------------------------

def ticket(name, claim, since=0.0):
    return scheduler.Ticket(name, claim, since, 1, "h")


def test_the_solo_ticket_goes_before_every_other_and_holds_the_bulk_lanes_back():
    solo = ticket("s", memorygate.with_solo_mb(solo_claim(), host(), places=0), since=100.0)
    interactive = ticket("i", Claim("api:x", 100.0, 1.0, "interactive"), since=0.0)
    unit = ticket("u", Claim("unit:y", 100.0, 1.0, "bulk"), since=0.0)
    assert [t.name for t in scheduler.rank([unit, interactive, solo], now=100.0, used={})][0] == "s"
    held = scheduler.judge(unit, [solo], host(), [], 100.0, {}, places=0)
    assert not held.admit and "kept for the solo run" in held.reason
    assert scheduler.judge(interactive, [solo], host(), [], 100.0, {}, places=0).admit


def test_a_solo_run_is_not_a_turn_of_its_kind_so_nothing_is_behind_it_afterwards(conn):
    memorygate.acquire(conn, Claim("unit:scatter", 500.0, 2.0, "bulk"), "plain", pid=os.getpid())
    conn.hset(memorygate._key("claims"), "plain", json.dumps({**json.loads(conn.hget(memorygate._key("claims"), "plain")), "started": time.time() - 100}))
    memorygate.release(conn, "plain")
    plain = scheduler.usage(conn)["unit:scatter"]
    assert plain > 150.0
    memorygate.acquire(conn, solo_claim("unit:grid"), "solo", pid=os.getpid())
    conn.hset(memorygate._key("claims"), "solo", json.dumps({**json.loads(conn.hget(memorygate._key("claims"), "solo")), "started": time.time() - 5000}))
    memorygate.release(conn, "solo")
    used = scheduler.usage(conn)
    assert "unit:grid" not in used and used["unit:scatter"] == pytest.approx(plain, rel=0.01)


def test_the_jobs_held_during_a_solo_run_go_at_once_when_it_ends_and_keep_the_waiting_they_did(conn, monkeypatch):
    monkeypatch.setattr(memorygate, "read_host", lambda: host())
    monkeypatch.setattr(memorygate, "kept_mb", lambda held, places=None: 0.0)
    now = [1000.0]
    sched = scheduler.Scheduler(places=0, clock=lambda: now[0])
    claim = memorygate.with_solo_mb(solo_claim(), host(), places=0)
    solo = Holder("solo", claim, 1, "h", 300.0, 300.0)
    waiting = Claim("unit:y", 400.0, 1.0, "bulk")
    first = sched.decide(conn, waiting, "w9", host(), [solo], 0.0, waited=500.0, pid=9)
    assert not first.admit
    since = json.loads(conn.hget(memorygate._key("tickets"), "w9"))["since"]
    now[0] += 3000.0
    later = sched.decide(conn, waiting, "w9", host(), [solo], 0.0, waited=3500.0, pid=9)
    assert not later.admit and json.loads(conn.hget(memorygate._key("tickets"), "w9"))["since"] == since
    assert sched.decide(conn, waiting, "w9", host(), [], 0.0, waited=3500.0, pid=9).admit, "the solo claim ended: the first job that waited goes"


# --- the workers ----------------------------------------------------------------------------------------------------



def _sleeping_job(queue, seconds):
    return queue.enqueue(api_jobs.execute, args=(api_jobs.function_reference(api_jobs.run_command), [["sleep", str(seconds)], dict(os.environ), "/", 30, False]), job_timeout=60)


def _quick_job(queue, at_front=False):
    return queue.enqueue(api_jobs.execute, args=(api_jobs.function_reference(api_jobs.run_command), [["true"], dict(os.environ), "/", 30, False]), job_timeout=60,
                         at_front=at_front)


def _worker(name, redis_url):
    argv = redisqueue.worker_argv([name], redis_url, python=sys.executable, burst=False)
    argv[argv.index(name):argv.index(name)] = ["--idle-seconds", "8"]
    return subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


@pytest.mark.real_queue
def test_a_restarted_unit_runs_alone_and_the_job_a_worker_took_beside_it_waits_for_it(redis_server):
    connection = redisqueue.connect()
    pause.resume(connection)
    name = f"planetgen-bulk-test-{uuid.uuid4().hex[:12]}"
    queue = redisqueue.queue(name, connection)
    other = _quick_job(queue)
    unit = _sleeping_job(queue, 6)
    queue.remove(unit.id)
    redisqueue.requeue_job(unit, connection, at_front=True)
    assert queue.job_ids[0] == unit.id
    assert pause.hold_solo("run-t", 0, 1, unit.id, connection)
    workers = [_worker(name, redis_server), _worker(name, redis_server)]
    try:
        time.sleep(4)
        assert unit.get_status(refresh=True) == "started" and other.get_status(refresh=True) == "queued", "the other job is turned away, not run beside the unit"
        assert pause.is_paused(connection)
        deadline = time.monotonic() + 40
        while other.get_status(refresh=True) != "finished" and time.monotonic() < deadline:
            time.sleep(0.5)
        assert unit.get_status(refresh=True) == "finished" and other.get_status(refresh=True) == "finished"
        assert not pause.is_paused(connection), "the hold lifted by itself when the unit's job ended"
    finally:
        for process in workers:
            process.terminate()
            process.wait(30)
        pause.resume(connection)
        queue.empty()
