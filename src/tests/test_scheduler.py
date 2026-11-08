"""OPS.130: the fair-share scheduler on top of the memory gate (queue/scheduler.py): priority, aging and its bound, fair share, the reserved places,
and that a large claim is not starved by a stream of small ones."""

import json
import os
import time
import types

import pytest

from planetgen.queue import memorygate, redisqueue, runrecord, scheduler
from planetgen.queue.memorygate import Claim, Holder, Host
from planetgen.queue.scheduler import Ticket


@pytest.fixture
def conn(redis_server, key_prefix, monkeypatch):
    monkeypatch.setattr(memorygate, "PREFIX", f"{key_prefix}:gate")
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    return redisqueue.connect(redis_server)


@pytest.fixture
def fake_host(monkeypatch):
    state = types.SimpleNamespace(host=Host(16000.0, 8000.0, 8, 0.5), held={}, dead=set())
    monkeypatch.setattr(memorygate, "read_host", lambda: state.host)
    monkeypatch.setattr(memorygate, "held_mb", lambda pid: state.held.get(pid, 0.0))
    monkeypatch.setattr(memorygate, "alive", lambda pid: pid not in state.dead)
    monkeypatch.setattr(memorygate, "kept_mb", lambda held, places=None: 0.0)
    return state


def host(available=8000.0, total=16000.0, cores=8, load=0.5):
    return Host(total, available, cores, load)


def holder(name, mb, held=0.0, lane="bulk", cpu=1.0, kind="unit:x"):
    return Holder(name, Claim(kind, mb, cpu, lane), 1, "h", held, held)


def ticket(name, kind="unit:x", mb=500.0, lane="bulk", since=0.0, cpu=1.0, exclusive=False):
    return Ticket(name, Claim(kind, mb, cpu, lane, exclusive), since, 1, "h")


# --- the order --------------------------------------------------------------------------------------------------------

def test_the_site_goes_before_a_run_and_a_run_before_the_upkeep():
    interactive = ticket("a", "api:x", lane="interactive")
    long_job = ticket("b", "api:y", lane="long")
    prefetch = ticket("c", "fn:planetgen.web.lib.tileprefetch.build_nodes", lane="maintenance")
    unit = ticket("d", "unit:scatter-mass", lane="bulk")
    upkeep = ticket("e", "unit:warm-search", lane="maintenance")
    order = [t.name for t in scheduler.rank([upkeep, unit, prefetch, long_job, interactive], now=0.0, used={})]
    assert order == ["a", "c", "b", "d", "e"]


def test_a_job_that_has_waited_the_bound_outranks_a_fresh_job_of_any_class_whatever_the_fair_share():
    """The worst case the order adds: no job waits longer than `AGING_BOUND_SECONDS` behind jobs that arrived after it, even the lightest class with the
    heaviest use against the heaviest class with none."""
    now = 10_000.0
    used = {"unit:warm-search": 1000.0}                                    # all of the use is the waiter's kind: the full penalty
    waiter = ticket("w", "unit:warm-search", lane="maintenance", since=now - scheduler.AGING_BOUND_SECONDS)
    fresh = ticket("f", "api:x", lane="interactive", since=now)
    assert scheduler.priority(waiter, now, used) >= scheduler.priority(fresh, now, used) - 1e-9
    just_short = ticket("s", "unit:warm-search", lane="maintenance", since=now - scheduler.AGING_BOUND_SECONDS + 30.0)
    assert scheduler.priority(just_short, now, used) < scheduler.priority(fresh, now, used)


def test_within_a_class_the_kind_that_has_used_less_goes_first():
    used = {"unit:warm-search": 900.0, "unit:check-math": 100.0}
    heavy = ticket("heavy", "unit:warm-search", lane="maintenance", since=0.0)
    light = ticket("light", "unit:check-math", lane="maintenance", since=0.0)
    assert [t.name for t in scheduler.rank([heavy, light], 0.0, used)] == ["light", "heavy"]
    assert [t.name for t in scheduler.rank([heavy, light], 0.0, {})] == ["heavy", "light"]               # no use known: first come, by name


def test_use_is_charged_when_a_claim_is_given_up_and_halves_over_the_half_life(conn, fake_host):
    memorygate.acquire(conn, Claim("unit:x", 500.0, cpu=2.0), "w1", pid=os.getpid())
    raw = json.loads(conn.hget(memorygate._key("claims"), "w1"))
    raw["started"] -= 100.0
    conn.hset(memorygate._key("claims"), "w1", json.dumps(raw))
    memorygate.release(conn, "w1")
    now = time.time()
    assert scheduler.usage(conn, now)["unit:x"] == pytest.approx(200.0, rel=0.05)           # 100 s on 2 cores
    assert scheduler.usage(conn, now + scheduler.USAGE_HALF_LIFE)["unit:x"] == pytest.approx(100.0, rel=0.05)
    scheduler.charge(conn, "unit:x", 100.0, now=now)
    assert scheduler.usage(conn, now)["unit:x"] == pytest.approx(300.0, rel=0.05)


# --- the rule ---------------------------------------------------------------------------------------------------------

def test_a_lower_job_waits_for_a_higher_one_that_fits_only_when_the_two_do_not_fit_together():
    busy = [holder("w1", 100.0, held=100.0)]
    high, low = ticket("high", "api:x", mb=1200.0, lane="interactive"), ticket("low", "unit:x", mb=600.0, lane="bulk")
    assert scheduler.judge(low, [high], host(available=8000.0), busy, 0.0, {}, places=0).admit              # room for both: nobody is held up
    tight = host(available=2000.0)                                                                           # either fits, not both
    blocked = scheduler.judge(low, [high], tight, busy, 0.0, {}, places=0)
    assert not blocked.admit and "goes first" in blocked.reason
    assert scheduler.judge(high, [low], tight, busy, 0.0, {}, places=0).admit                                # and the higher one is not held up by the lower


def test_a_higher_job_that_does_not_fit_does_not_hold_up_a_lower_one_that_does_until_it_has_waited():
    busy = [holder("w1", 100.0, held=100.0)]
    small_room = host(available=2000.0)                      # 2,500 MB for the interactive job does not fit; 500 MB for the bulk one does
    big = ticket("big", "api:big", mb=2500.0, lane="interactive", since=0.0)
    small = ticket("small", "unit:x", mb=500.0, lane="bulk", since=0.0)
    assert scheduler.judge(small, [big], small_room, busy, 10.0, {}, places=0).admit                 # backfill: the host is not left idle for order
    kept = scheduler.judge(small, [big], small_room, busy, scheduler.RESERVE_AFTER_SECONDS + 1.0, {}, places=0)
    assert not kept.admit and "kept for" in kept.reason                                              # the room is kept for the one that has waited
    over = scheduler.judge(small, [big], small_room, busy, memorygate.MAX_WAIT_SECONDS + 1.0, {}, places=0)
    assert over.admit                                                                                # past the gate's own bound it is no longer kept


def test_a_large_claim_is_not_starved_by_a_stream_of_small_ones(monkeypatch):
    """A 4,000 MB job waits for room while a small one arrives every 5 s and holds its 800 MB for 30 s. Without the kept room the large job never fits;
    with it, the room is kept once it has waited `RESERVE_AFTER_SECONDS`, the small ones in memory end, and it goes before the gate's own bound."""

    def simulate(seconds):
        plenty = host(available=6000.0)
        admitted = {f"early{n}": -25 + 5 * n for n in range(3)}                  # the stream is under way: three small jobs are in memory already
        held = [Holder(name, Claim("unit:small", 800.0), 1, "h", 0.0, 0.0) for name in admitted]
        waiting = {}
        big = ticket("big", "unit:big", mb=4000.0, since=0.0)
        waiting["big"] = big
        for now in range(0, seconds, 5):
            held = [h for h in held if admitted[h.name] + 30 > now]
            waiting["small%d" % now] = ticket("small%d" % now, "unit:small", mb=800.0, since=now)
            for name, current in sorted(waiting.items(), key=lambda item: -scheduler.priority(item[1], now, {})):
                others = [t for t in waiting.values() if t.name != name]
                if scheduler.judge(current, others, plenty, held, float(now), {}, places=0).admit:
                    held.append(Holder(name, current.claim, 1, "h", 0.0, 0.0))
                    admitted[name] = now
                    del waiting[name]
                    break               # one decision at a time, as under the gate's lock
        return admitted

    admitted = simulate(300)
    assert "big" in admitted and admitted["big"] < memorygate.MAX_WAIT_SECONDS
    monkeypatch.setattr(scheduler, "RESERVE_AFTER_SECONDS", 1e9)
    assert "big" not in simulate(300)                    # the control: without the kept room the stream starves it


def test_the_weighted_share_between_kinds_follows_the_use(monkeypatch):
    """Two kinds of upkeep wait for one place at a time; the one that has used more of the host goes later, and goes at once when the other has used as much."""
    now = 0.0
    one_place = host(available=memorygate.reserve_mb(host()) + 1500.0)
    busy = [holder("w1", 100.0, held=100.0)]
    used = {"unit:warm-search": 600.0, "unit:check-math": 0.0}
    heavy = ticket("heavy", "unit:warm-search", mb=900.0, lane="maintenance")
    light = ticket("light", "unit:check-math", mb=900.0, lane="maintenance")
    assert scheduler.judge(light, [heavy], one_place, busy, now, used, places=0).admit
    assert not scheduler.judge(heavy, [light], one_place, busy, now, used, places=0).admit
    used["unit:check-math"] = 600.0
    assert scheduler.judge(ticket("heavy", "unit:warm-search", mb=900.0, lane="maintenance", since=-1.0), [light], one_place, busy, now, used, places=0).admit


# --- the places kept for the site -------------------------------------------------------------------------------------

def test_the_places_kept_for_the_short_lanes_keep_cores_against_the_other_lanes():
    eight = host(cores=8)
    assert memorygate.kept_cores([], eight, places=3) == pytest.approx(1.5)                           # 3 places at half a core
    assert memorygate.kept_cores([], eight, places=40) == pytest.approx(2.0)                          # never more than a quarter of the cores
    assert memorygate.kept_cores([holder("i", 100.0, lane="interactive", cpu=0.5)], eight, places=3) == pytest.approx(1.0)
    assert memorygate.kept_cores([holder("i", 100.0, lane="long", cpu=4.0)], eight, places=3) == 0.0       # the places are in use


def test_bulk_work_cannot_take_the_cores_a_visitor_job_needs_but_the_visitor_job_can():
    busy = [holder(f"w{n}", 100.0, held=100.0) for n in range(6)]                                       # six of eight cores in use by the bulk lane
    eight = host(cores=8)
    assert memorygate.decide(Claim("unit:x", 100.0), eight, busy, places=0).admit                          # no places kept: the seventh core is free
    refused = memorygate.decide(Claim("unit:x", 100.0), eight, busy, places=3)
    assert not refused.admit and "core" in refused.reason                                              # 1.5 of the 8 are kept for the site
    assert memorygate.decide(Claim("api:x", 100.0, cpu=0.5, lane="interactive"), eight, busy, places=3).admit


def test_a_warm_up_unit_waits_for_a_visitors_job_that_asked_after_it(conn, fake_host):
    """End to end in Redis: a warm-up unit and a visitor's job both wait for room; when it frees, only one fits, and the visitor's job goes first."""
    memorygate.acquire(conn, Claim("unit:big", 100.0, lane="bulk"), "w0", pid=1000)
    sched = scheduler.Scheduler(places=0)
    warm = Claim("unit:warm-search", 900.0, lane="maintenance")
    visitor = Claim("api:planetgen.queue.pagejobs:build_answer", 900.0, lane="interactive")
    fake_host.host = host(available=500.0)                                                                # nothing fits yet
    assert not sched.decide(conn, warm, "warm", fake_host.host, memorygate.holders(conn), 0.0, waited=60.0, pid=2001).admit
    assert not sched.decide(conn, visitor, "visitor", fake_host.host, memorygate.holders(conn), 0.0, waited=1.0, pid=2002).admit
    fake_host.host = host(available=2300.0)                                                               # room for either, not both
    turned_away = sched.decide(conn, warm, "warm", fake_host.host, memorygate.holders(conn), 0.0, waited=65.0, pid=2001)
    assert not turned_away.admit and "goes first" in turned_away.reason                                   # it asked first, and still waits
    assert sched.decide(conn, visitor, "visitor", fake_host.host, memorygate.holders(conn), 0.0, waited=6.0, pid=2002).admit
    memorygate._write(conn, "visitor", visitor, 2002, "", time.time())
    assert not sched.decide(conn, warm, "warm", fake_host.host, memorygate.holders(conn), 0.0, waited=70.0, pid=2001).admit
    memorygate.release(conn, "visitor")
    assert sched.decide(conn, warm, "warm", fake_host.host, memorygate.holders(conn), 0.0, waited=75.0, pid=2001).admit


def test_a_ticket_is_kept_while_the_job_waits_and_dropped_when_it_goes_or_its_worker_is_gone(conn, fake_host):
    fake_host.host = host(available=3000.0)
    memorygate.acquire(conn, Claim("unit:big", 2500.0), "w1", pid=1001)
    sched = scheduler.Scheduler(places=0, clock=lambda: 1000.0)
    claim = Claim("unit:big", 2500.0)
    decision = sched.decide(conn, claim, "w2", fake_host.host, memorygate.holders(conn), 0.0, waited=40.0, pid=1002)
    assert not decision.admit
    (waiting,) = scheduler.tickets(conn, now=1000.0)
    assert waiting.name == "w2" and waiting.since == pytest.approx(960.0)                                # waiting since 40 s before the first poll
    sched.decide(conn, claim, "w2", fake_host.host, memorygate.holders(conn), 0.0, waited=45.0, pid=1002)
    assert scheduler.tickets(conn, now=1000.0)[0].since == pytest.approx(960.0)                          # the wait is not restarted by a poll
    assert scheduler.tickets(conn, now=1000.0 + scheduler.TICKET_SECONDS + 1.0) == []                    # a worker that stopped polling is forgotten
    sched.decide(conn, claim, "w2", fake_host.host, memorygate.holders(conn), 0.0, waited=50.0, pid=1002)
    fake_host.dead.add(1002)
    assert scheduler.tickets(conn, now=1000.0) == []                                                     # and so is one whose process is gone
    sched.decide(conn, claim, "w3", fake_host.host, memorygate.holders(conn), 0.0, waited=0.0, pid=1003)
    sched.done(conn, "w3")
    assert scheduler.tickets(conn, now=1000.0) == []


def test_acquire_with_the_scheduler_admits_in_order_and_leaves_no_ticket(conn, fake_host):
    fake_host.host = host(available=3000.0)
    memorygate.acquire(conn, Claim("unit:big", 2500.0), "w1", pid=1001)
    polls = []

    def on_wait(decision, waited):
        polls.append(waited)
        if len(polls) >= 3:
            memorygate.release(conn, "w1")

    admission = memorygate.acquire(conn, Claim("unit:big", 2500.0), "w2", pid=1002, sleep=lambda s: None, on_wait=on_wait, scheduler=scheduler.Scheduler(places=0))
    assert admission.decision.admit and len(polls) >= 3
    assert scheduler.tickets(conn) == [] and [h.name for h in memorygate.holders(conn)] == ["w2"]


def test_hold_uses_the_scheduler_and_charges_the_kind(conn, fake_host):
    with memorygate.hold(conn, Claim("unit:x", 300.0)):
        assert scheduler.tickets(conn) == []
    assert "unit:x" in scheduler.usage(conn)


def test_the_snapshot_lists_the_waiting_in_the_order_they_would_go(conn, fake_host):
    fake_host.host = host(available=3000.0)
    memorygate.acquire(conn, Claim("unit:big", 2500.0), "w1", pid=1001)
    for name, claim in (("a", Claim("unit:warm", 2500.0, lane="maintenance")), ("b", Claim("api:x", 2500.0, lane="interactive"))):
        scheduler.Scheduler(places=0).decide(conn, claim, name, fake_host.host, memorygate.holders(conn), 0.0, waited=0.0, pid=os.getpid())
    shown = scheduler.snapshot(conn)
    assert [t.name for t, _ in shown["waiting"]] == ["b", "a"]
