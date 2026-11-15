"""OPS.232: named resources (`planetgen/queue/tokens.py`). A claim holds tokens besides memory; the gate lets a claim go only when no claim held holds one of its names, and the restarted unit that runs
alone (OPS.145) is the claim that holds `alone`."""

import dataclasses

import pytest

from planetgen.queue import claims, commands, memorygate, scheduler, tokens
from planetgen.queue.memorygate import Claim, Holder, Host

HOST = Host(16000.0, 12000.0, 8, 0.1)


def claim(kind, *held, lane=memorygate.LANE_BULK, mb=500.0):
    return Claim(kind, mb, 1.0, lane, tokens=tuple(held))


def holder(name, *the_claim):
    return Holder(name, the_claim[0], 1, "h", the_claim[0].mb, the_claim[0].mb)


def test_a_token_is_held_by_one_claim_at_a_time_and_the_waiter_is_told_who_holds_it():
    writer = claim("unit:bulk-writer", tokens.MYSQL_WRITE)
    other = claim("unit:other-writer", tokens.MYSQL_WRITE)
    plain = claim("unit:plain")
    held = [holder("a", writer)]
    refused = memorygate.decide(other, HOST, held, places=2, name="b")
    assert not refused.admit and refused.reason == "waits for unit:bulk-writer to release mysql_write"
    assert memorygate.decide(plain, HOST, held, places=2, name="c").admit, "a claim without the token is not held back"
    assert memorygate.decide(other, HOST, [], places=2, name="b").admit
    assert memorygate.decide(writer, HOST, held, places=2, name="a").admit, "a claim never waits for its own old holding"


def test_the_token_of_a_site_job_holds_too_and_each_warm_set_is_its_own_name():
    held = [holder("a", claim("fn:build-tables", tokens.warm_set("bright-star-tables"), lane=memorygate.LANE_OTHER))]
    assert tokens.warm_set("bright-star-tables") == "warm:bright-star-tables"
    same = claim("fn:build-tables-again", tokens.warm_set("bright-star-tables"), lane=memorygate.LANE_OTHER)
    different = claim("fn:build-shape", tokens.warm_set("shape"), lane=memorygate.LANE_OTHER)
    assert not memorygate.decide(same, HOST, held, places=2, name="b").admit
    assert memorygate.decide(different, HOST, held, places=2, name="c").admit
    visitor = claim("fn:page", tokens.MYSQL_WRITE, lane=memorygate.LANE_INTERACTIVE, mb=300.0)
    assert not memorygate.decide(visitor, HOST, [holder("w", claim("unit:w", tokens.MYSQL_WRITE))], places=2, name="v").admit


def test_the_alone_token_is_the_solo_claim_and_the_gate_treats_it_so():
    solo = claim("unit:scatter", tokens.ALONE)
    assert solo.solo and not claim("unit:scatter").solo
    busy = memorygate.decide(solo, HOST, [holder("a", claim("unit:other"))], places=2, name="s")
    assert not busy.admit and busy.reason == "a solo run waits for unit:other to end"
    assert memorygate.decide(solo, HOST, [holder("v", claim("fn:page", lane=memorygate.LANE_INTERACTIVE, mb=300.0))], places=2, name="s").admit
    behind = memorygate.decide(claim("unit:x"), HOST, [holder("s", solo)], places=2, name="x")
    assert not behind.admit and behind.reason == "a solo run holds the host"
    assert scheduler.weight(solo) == scheduler.SOLO_WEIGHT


def test_a_claim_with_tokens_round_trips_through_the_gates_rows():
    original = claim("unit:x", tokens.MYSQL_WRITE, tokens.warm_set("shape"))
    row = {"kind": original.kind, "mb": original.mb, "cpu": original.cpu, "lane": original.lane, "exclusive": False, "tokens": list(original.tokens)}
    assert memorygate.claim_of(row) == original
    assert memorygate.claim_of({k: v for k, v in row.items() if k != "tokens"}).tokens == (), "a row without tokens holds none"


def test_a_command_names_its_tokens_and_the_claim_of_its_unit_holds_them(monkeypatch):
    spec = commands.spec("scatter-mass")
    assert spec.tokens == (), "none of the commands asks for a token yet"
    monkeypatch.setitem(commands.REGISTRY, "scatter-mass", dataclasses.replace(spec, tokens=(tokens.MYSQL_WRITE,)))
    assert claims.unit_tokens("scatter-mass") == (tokens.MYSQL_WRITE,)
    assert claims.unit_tokens("no-such-command") == ()


def test_the_solo_hold_gives_the_unit_the_alone_token(monkeypatch):
    class Job:
        id, func_name, origin, args = "job-1", claims.UNIT_FUNCTION, "bulk", ("run", 0, 0, "scatter-mass", {}, 1.0)

    monkeypatch.setattr(claims.pause, "solo_of", lambda conn: {"job": "job-1"})
    monkeypatch.setattr(claims.memorygate, "claim_mb", lambda conn, kind, mb, job_id=None, token=None: mb)
    monkeypatch.setattr(claims, "lane_of", lambda origin: memorygate.LANE_BULK)
    assert tokens.ALONE in claims.claim_for(Job(), None).tokens
    monkeypatch.setattr(claims.pause, "solo_of", lambda conn: {"job": "another"})
    assert tokens.ALONE not in claims.claim_for(Job(), None).tokens


# --- the simulator -----------------------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("seed", [0, 1])
def test_in_the_simulator_a_token_is_never_held_twice_and_the_job_that_runs_alone_has_the_bulk_and_maintenance_lanes_to_itself(seed):
    from tests import enginesim
    measures, broken = enginesim.run(seed, 120, named=True)
    assert not broken, broken
    assert measures["makespan"] < 100000.0, "every job finished: a token never held a job back for ever"
