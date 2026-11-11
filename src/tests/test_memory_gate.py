"""PERF.254: the memory and load gate every job of the work queue passes (queue/memorygate.py, queue/claims.py)."""

import json
import os
import signal
import types

import pytest

from planetgen.queue import allowlist, claims, commands, memorygate, redisqueue, runrecord, units, workers
from planetgen.queue.memorygate import Claim, Host, Holder


@pytest.fixture
def conn(redis_server, key_prefix, monkeypatch):
    """Redis under this test's own prefix, for the gate and the run records."""
    monkeypatch.setattr(memorygate, "PREFIX", f"{key_prefix}:gate")
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    return redisqueue.connect(redis_server)


def host(available=8000.0, total=16000.0, cores=8, load=0.5):
    return Host(total, available, cores, load)


def holder(name, mb, held=0.0, lane="bulk", cpu=1.0, kind="unit:x"):
    return Holder(name, Claim(kind, mb, cpu, lane), 1, "h", held, held)


# --- every kind of job has a claim ------------------------------------------------------------------------------------

def test_every_function_a_queue_may_run_has_a_claim_or_runs_inside_another_job():
    """Nothing is exempt from the gate (Boss, 3:34 PM EST on 2026-11-08): a new kind of job fails here until it has a claim."""
    for function in allowlist.FUNCTIONS:
        assert function in claims.FUNCTION_CLAIMS or function in claims.SPECIAL_FUNCTIONS or function in claims.CALLBACKS, function
    assert set(claims.FUNCTION_CLAIMS) <= allowlist.FUNCTIONS
    assert claims.CALLBACKS <= allowlist.FUNCTIONS and claims.SPECIAL_FUNCTIONS <= allowlist.FUNCTIONS
    assert not set(claims.FUNCTION_CLAIMS) & (claims.SPECIAL_FUNCTIONS | claims.CALLBACKS)


def test_every_function_the_api_dispatcher_may_call_has_a_claim():
    assert set(claims.API_CLAIMS) == set(allowlist.API_FUNCTIONS)
    for reference in allowlist.API_FUNCTIONS:
        mb, cpu = claims.api_table(reference)
        assert mb >= memorygate.MIN_CLAIM_MB and cpu > 0


def test_every_command_with_units_has_a_memory_figure():
    for kind, spec in commands.REGISTRY.items():
        if not spec.inline:
            assert claims.unit_table_mb(kind) >= memorygate.MIN_CLAIM_MB, kind
    assert not set(claims.UNIT_BY_NATURE_MB) & set(workers.STAGE_WORKER_MEMORY_MB)       # one figure for a command
    assert set(claims.UNIT_BY_NATURE_MB) <= set(commands.REGISTRY)


def test_a_kind_with_no_figure_is_refused_by_the_table_and_gets_the_default_at_the_worker():
    with pytest.raises(claims.NoClaim):
        claims.function_table("planetgen.nothing.here")
    with pytest.raises(claims.NoClaim):
        claims.unit_table_mb("link-neighbors-but-not")
    with pytest.raises(claims.NoClaim):
        claims.api_table("planetgen.api.nothing:job")
    job = types.SimpleNamespace(func_name="planetgen.nothing.here", origin=redisqueue.INTERACTIVE_QUEUE)
    assert claims.default_claim(job) == Claim("fn:planetgen.nothing.here", float(claims.DEFAULT_MB), 1.0, "interactive")


def test_the_lane_of_each_queue():
    assert claims.lane_of(redisqueue.BULK_QUEUE) == "bulk" and claims.lane_of(redisqueue.GENERATION_QUEUE) == "bulk"
    for name in redisqueue.maintenance_lane():
        assert claims.lane_of(name) == "maintenance"
    assert claims.lane_of(redisqueue.INTERACTIVE_QUEUE) == "interactive" and claims.lane_of(redisqueue.LONG_QUEUE) == "long"
    assert claims.lane_of("something-else") == "other"


def test_claim_for_a_unit_an_api_job_a_page_job_and_a_cut(conn):
    unit = types.SimpleNamespace(func_name=claims.UNIT_FUNCTION, args=("r", 3, 0, "scatter-mass", {}, 1.0), origin=redisqueue.BULK_QUEUE, id="r-s3-u0")
    assert claims.claim_for(unit, conn) == Claim("unit:scatter-mass", float(workers.SCATTER_WORKER_MEMORY_MB), 1.0, "bulk")
    api = types.SimpleNamespace(func_name=claims.EXECUTE_FUNCTION, args=("planetgen.queue.api_jobs:regenerate_sector", {}), origin=redisqueue.LONG_QUEUE, id="j")
    assert claims.claim_for(api, conn) == Claim("api:planetgen.queue.api_jobs:regenerate_sector", float(workers.SECTOR_WORKER_MEMORY_MB), 1.0, "long")
    page = types.SimpleNamespace(func_name="planetgen.web.lib.tileprefetch.build_nodes", args=(), origin=redisqueue.NODES_QUEUE, id="n")
    assert claims.claim_for(page, conn).lane == "maintenance"
    steps = [commands.Step("a", "link-neighbors", "A").as_dict(), commands.Step("b", "scatter-backfill", "B").as_dict()]
    conn.hset(runrecord.run_key("r"), mapping={"stages": json.dumps([{"key": "a", "label": "A"}, {"key": "b", "label": "B"}]), "steps": json.dumps(steps)})
    cut = types.SimpleNamespace(func_name=claims.OPEN_FUNCTION, args=("r", 1), origin=redisqueue.BULK_QUEUE, id="r-open1")
    assert claims.claim_for(cut, conn) == Claim("cut:scatter-backfill", float(claims.CUT_MB), 1.0, "bulk", True)
    memorygate.learn(conn, "cut:scatter-backfill", 2000.0)
    assert claims.claim_for(cut, conn).mb == pytest.approx(2000.0 * memorygate.LEARN_MARGIN)       # what the first cut held is what the next claims
    gone = types.SimpleNamespace(func_name=claims.OPEN_FUNCTION, args=("nobody", 4), origin=redisqueue.BULK_QUEUE, id="x")
    assert claims.claim_for(gone, conn).kind == "cut:unknown"


# --- the rule ---------------------------------------------------------------------------------------------------------

def test_a_job_alone_is_admitted_whatever_it_claims():
    decision = memorygate.decide(Claim("unit:x", 99999.0), host(available=100.0), [])
    assert decision.admit and "alone" in decision.reason


def test_a_claim_is_held_against_what_is_available_less_the_reserve_and_the_claims_not_yet_in_memory():
    reserve = memorygate.reserve_mb(host())
    assert reserve == max(memorygate.RESERVE_MIN_MB, memorygate.RESERVE_SHARE * 16000.0)
    other = holder("w1", 3000.0, held=1000.0)                          # 2,000 MB of its claim are still to come
    free = 8000.0 - reserve - memorygate.kept_mb([other], places=0) - 2000.0
    fits = memorygate.decide(Claim("unit:x", free - 1.0), host(), [other], places=0)
    assert fits.admit and fits.free_mb == pytest.approx(free)
    refused = memorygate.decide(Claim("unit:x", free + 1.0), host(), [other], places=0)
    assert not refused.admit and "needs" in refused.reason and "still to be taken" in refused.reason


def test_a_claim_that_is_in_memory_already_is_not_counted_twice():
    grown = holder("w1", 3000.0, held=3500.0)                         # past its claim: nothing more to come
    assert memorygate.decide(Claim("unit:x", 6000.0), host(available=8000.0), [grown], places=0).admit


def test_the_worker_asking_is_credited_with_what_it_holds_already():
    other = holder("w1", 2000.0, held=2000.0)
    claim = Claim("unit:scatter", 2100.0)
    assert not memorygate.decide(claim, host(available=3000.0), [other], places=0).admit
    assert memorygate.decide(claim, host(available=3000.0), [other], own_mb=1400.0, places=0).admit


def test_the_places_kept_for_the_short_lanes_hold_memory_against_the_bulk_lane_only():
    kept = 3 * memorygate.KEEP_MB_PER_PLACE
    other = holder("w1", 100.0, held=100.0)
    tight = host(available=memorygate.reserve_mb(host()) + kept + 500.0)
    assert not memorygate.decide(Claim("unit:x", 600.0, lane="bulk"), tight, [other], places=3).admit
    assert memorygate.decide(Claim("api:x", 600.0, lane="interactive"), tight, [other], places=3).admit
    in_use = holder("w2", kept, held=kept, lane="interactive")
    assert memorygate.kept_mb([in_use], places=3) == 0.0                 # the places are in use: nothing more to keep


def test_cores_are_claimed_too():
    others = [holder(f"w{n}", 100.0, held=100.0) for n in range(8)]
    decision = memorygate.decide(Claim("unit:x", 100.0), host(cores=8), others, places=0)
    assert not decision.admit and "core" in decision.reason
    assert memorygate.decide(Claim("api:x", 100.0, cpu=0.25, lane="interactive"), host(cores=8), others, places=0).admit       # a quarter over for the short lanes


def test_a_loaded_host_makes_an_ordinary_job_wait_and_an_interactive_one_less_so():
    other = [holder("w1", 100.0, held=100.0)]
    busy = host(load=1.5)
    assert not memorygate.decide(Claim("unit:x", 100.0), busy, other, places=0).admit
    assert memorygate.decide(Claim("api:x", 100.0, lane="interactive"), busy, other, places=0).admit
    assert not memorygate.decide(Claim("api:x", 100.0, lane="interactive"), host(load=2.5), other, places=0).admit


def test_the_cut_waits_for_the_bulk_lane_to_be_clear():
    cut = Claim("cut:scatter-backfill", 1500.0, lane="bulk", exclusive=True)
    assert not memorygate.decide(cut, host(), [holder("w1", 700.0, held=700.0)], places=0).admit
    assert memorygate.decide(cut, host(), [holder("w1", 700.0, held=700.0, lane="interactive")], places=0).admit
    assert memorygate.decide(cut, host(), [], places=0).admit


# --- the gate in Redis ------------------------------------------------------------------------------------------------

@pytest.fixture
def fake_host(monkeypatch):
    state = types.SimpleNamespace(host=host(), held={}, dead=set())
    monkeypatch.setattr(memorygate, "read_host", lambda: state.host)
    monkeypatch.setattr(memorygate, "held_mb", lambda pid: state.held.get(pid, 0.0))
    monkeypatch.setattr(memorygate, "alive", lambda pid: pid not in state.dead)
    monkeypatch.setattr(memorygate, "kept_mb", lambda held, places=None: 0.0)
    return state


def test_a_claim_is_taken_held_and_given_up(conn, fake_host):
    admission = memorygate.acquire(conn, Claim("unit:x", 500.0), "w1", job="j1", pid=os.getpid())
    assert admission.waited < 1 and admission.decision.admit and not admission.timed_out
    assert [h.name for h in memorygate.holders(conn)] == ["w1"]
    fake_host.held[os.getpid()] = 640.0
    assert memorygate.refresh(conn, "w1") == 640.0
    assert memorygate.release(conn, "w1") == 640.0
    assert memorygate.holders(conn) == [] and memorygate.learned(conn) == {"unit:x": 640.0}
    assert memorygate.release(conn, "w1") == 0.0                          # a second release is nothing


def test_a_job_waits_until_the_claim_that_stands_in_its_way_is_given_up(conn, fake_host):
    fake_host.host = host(available=3000.0)
    memorygate.acquire(conn, Claim("unit:big", 2500.0), "w1", pid=1001)
    waits = []

    def sleep(_seconds):
        waits.append(1)
        if len(waits) >= 3:
            memorygate.release(conn, "w1")

    admission = memorygate.acquire(conn, Claim("unit:big", 2500.0), "w2", pid=1002, sleep=sleep, on_wait=lambda decision, waited: waits.append(decision.reason))
    assert admission.decision.admit and len(waits) >= 4 and any(isinstance(item, str) and "needs" in item for item in waits)


def test_the_wait_is_bounded_and_the_table_cap_lets_the_job_in(conn, fake_host, monkeypatch):
    fake_host.host = host(available=3000.0)
    memorygate.acquire(conn, Claim("unit:big", 2500.0), "w1", pid=1001)
    ticks = iter(range(0, 10000, 100))
    clock = lambda: next(ticks)
    caps = {"unit:big": 2}
    admission = memorygate.acquire(conn, Claim("unit:big", 2500.0), "w2", pid=1002, clock=clock, sleep=lambda s: None, max_wait=300.0, cap_of=caps.get)
    assert admission.timed_out and admission.waited >= 300 and "table cap" in admission.decision.reason
    assert {h.name for h in memorygate.holders(conn)} == {"w1", "w2"}


def test_a_job_the_table_cap_would_still_hold_goes_on_waiting(conn, fake_host):
    fake_host.host = host(available=3000.0)
    memorygate.acquire(conn, Claim("unit:big", 2500.0), "w1", pid=1001)
    ticks = iter(range(0, 10000, 100))
    polls = []

    def on_wait(decision, waited):
        polls.append(waited)
        if waited > 900:
            memorygate.release(conn, "w1")          # the cap is 1 and w1 holds it: only its leaving lets w2 in

    admission = memorygate.acquire(conn, Claim("unit:big", 2500.0), "w2", pid=1002, clock=lambda: next(ticks), sleep=lambda s: None, max_wait=300.0,
                                   cap_of={"unit:big": 1}.get, on_wait=on_wait)
    assert not admission.timed_out and max(polls) > 900


def test_a_claim_of_a_dead_process_or_one_that_ran_out_is_dropped(conn, fake_host):
    memorygate.acquire(conn, Claim("unit:x", 100.0), "w1", pid=2001)
    memorygate.acquire(conn, Claim("unit:x", 100.0), "w2", pid=2002)
    fake_host.dead.add(2001)
    assert [h.name for h in memorygate.holders(conn)] == ["w2"]
    assert memorygate.holders(conn, now=10 ** 12) == []
    assert conn.hlen(memorygate._key("claims")) == 0


def test_the_gate_off_admits_at_once(conn, fake_host, monkeypatch):
    fake_host.host = host(available=10.0)
    memorygate.acquire(conn, Claim("unit:x", 5000.0), "w1", pid=3001)
    monkeypatch.setenv(memorygate.DISABLE_ENV, "0")
    admission = memorygate.acquire(conn, Claim("unit:x", 5000.0), "w2", pid=3002)
    assert admission.decision.admit and admission.decision.reason == "the gate is off"


def test_a_claim_is_the_table_figure_raised_by_what_was_seen_and_by_a_retry(conn):
    assert memorygate.claim_mb(conn, "unit:x", 700.0) == 700.0
    memorygate.learn(conn, "unit:x", 500.0)
    assert memorygate.claim_mb(conn, "unit:x", 700.0) == 700.0               # the table's figure stands while it is the larger
    memorygate.learn(conn, "unit:x", 1000.0)
    assert memorygate.claim_mb(conn, "unit:x", 700.0) == pytest.approx(1000.0 * memorygate.LEARN_MARGIN)
    assert memorygate.learn(conn, "unit:x", 300.0) == 1000.0                  # the largest peak is kept
    assert memorygate.claim_mb(conn, "unit:new", None) == memorygate.MIN_CLAIM_MB
    memorygate.learn(conn, "unit:new", 400.0)
    assert memorygate.claim_mb(conn, "unit:new", None) == pytest.approx(440.0)
    bigger = memorygate.grow(conn, "unit:x", "job-9", 1200.0, peak_mb=1300.0)
    assert bigger == pytest.approx(1800.0, rel=1e-3) and memorygate.claim_mb(conn, "unit:x", 700.0, job_id="job-9") == pytest.approx(1800.0, rel=1e-3)
    assert memorygate.claim_mb(conn, "unit:x", 700.0, job_id="job-10") == pytest.approx(1800.0, rel=1e-3)       # the kind is learned from it too
    assert conn.get(memorygate._key("override:job-10")) is None
    memorygate.forget(conn)
    assert memorygate.learned(conn) == {} and memorygate.claim_mb(conn, "unit:x", 700.0, job_id="job-9") == 700.0


def test_hold_is_a_claim_around_work_that_is_not_a_job(conn, fake_host):
    with memorygate.hold(conn, Claim("test-gate", 900.0), pid=os.getpid()) as name:
        assert [h.name for h in memorygate.holders(conn)] == [name]
    assert memorygate.holders(conn) == []


def test_the_snapshot_says_what_the_gate_holds(conn, fake_host):
    memorygate.acquire(conn, Claim("unit:x", 900.0), "w1", pid=4001)
    fake_host.held[4001] = 300.0
    snapshot = memorygate.snapshot(conn)
    assert snapshot["claimed_mb"] == 900.0 and snapshot["held_mb"] == 300.0 and snapshot["reserve_mb"] == memorygate.reserve_mb(host())


def test_the_real_host_reads_and_a_process_is_measured():
    reading = memorygate.read_host()
    assert reading.total_mb > 0 and 0 < reading.available_mb <= reading.total_mb and reading.cores >= 1 and reading.load_per_core >= 0
    assert memorygate.held_mb(os.getpid()) > 0 and memorygate.held_mb(2 ** 22 + 1) == 0.0


# --- a job killed for memory ------------------------------------------------------------------------------------------

def _killed_job(retries_left, func=units.UNIT_FUNCTION, conn=None):
    args = ("r1", 1, 101, "scatter-mass", {}, 1.0) if func == units.UNIT_FUNCTION else ("r1", 1)
    return types.SimpleNamespace(args=args, connection=conn, func_name=func, id="r1-s1-u101", worker_name="w1", retries_left=retries_left)


@pytest.fixture
def running(conn, fake_host):
    stages = [{"label": "Plan the galaxy"}, {"label": "Scatter the massive stars"}]
    conn.hset(runrecord.run_key("r1"), mapping={"state": "running", "kind": runrecord.GENERATION, "stages": json.dumps(stages)})
    memorygate.acquire(conn, Claim("unit:scatter-mass", 2100.0), "w1", job="r1-s1-u101", pid=5001)
    fake_host.held[5001] = 2600.0
    memorygate.refresh(conn, "w1")
    return conn


def test_a_unit_killed_with_sigkill_is_retried_once_with_a_larger_claim(running):
    units.horse_killed(_killed_job(1, conn=running), 4242, -signal.SIGKILL, None)
    assert runrecord._s(running.hget(runrecord.run_key("r1"), "state")) == "running"          # the run goes on
    assert memorygate.claim_mb(running, "unit:scatter-mass", 2100.0, job_id="r1-s1-u101") == pytest.approx(max(2100.0 * 1.5, 2600.0 * 1.1), rel=1e-3)
    assert memorygate.learned(running)["unit:scatter-mass"] > 2100.0


def test_a_unit_killed_again_stops_the_run_naming_the_unit_and_what_it_claimed(running):
    units.horse_killed(_killed_job(0, conn=running), 4242, -signal.SIGKILL, None)
    assert runrecord._s(running.hget(runrecord.run_key("r1"), "state")) == "failed"
    error = runrecord._s(running.hget(runrecord.run_key("r1"), "error"))
    assert error.startswith("unit 101 of stage Scatter the massive stars was killed: its worker process 4242 ended with")
    assert "it claimed 2,100 MB on the memory gate and the most seen was 2,600 MB" in error


def test_a_kill_that_is_not_sigkill_is_not_retried(running):
    units.horse_killed(_killed_job(1, conn=running), 4242, -signal.SIGTERM, None)
    assert runrecord._s(running.hget(runrecord.run_key("r1"), "state")) == "failed"


def test_an_opening_killed_with_sigkill_is_retried_once_then_the_run_stops_naming_the_stage(running):
    job = _killed_job(1, func=units.OPEN_FUNCTION, conn=running)
    units.horse_killed(job, 4242, -signal.SIGKILL, None)
    assert runrecord._s(running.hget(runrecord.run_key("r1"), "state")) == "running"
    job.retries_left = 0
    units.horse_killed(job, 4242, -signal.SIGKILL, None)
    error = runrecord._s(running.hget(runrecord.run_key("r1"), "error"))
    assert error.startswith("the opening of stage Scatter the massive stars was killed: its worker process 4242")


def test_a_job_that_is_not_part_of_a_run_is_only_learned_from(conn, fake_host):
    memorygate.acquire(conn, Claim("api:x", 700.0), "w1", job="j", pid=6001)
    fake_host.held[6001] = 900.0
    memorygate.refresh(conn, "w1")
    job = types.SimpleNamespace(args=(), connection=conn, func_name=claims.EXECUTE_FUNCTION, id="j", worker_name="w1", retries_left=0)
    units.horse_killed(job, 4242, -signal.SIGKILL, None)
    assert memorygate.claim_mb(conn, "api:x", 700.0) == pytest.approx(700.0 * 1.5, rel=1e-3)       # the next job of the kind claims more
