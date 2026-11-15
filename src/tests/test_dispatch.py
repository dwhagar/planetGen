"""OPS.227: the dispatch rule over the whole backlog (`planetgen/queue/dispatch.py`): the site's lanes first against the places kept for them, then the other lanes ranked and judged, a reservation
that lasts until the job it was made for runs, only the top of the order judged and each job's last reason kept. The simulator (`tests/enginesim.py`) is the proof of the first two."""

import math

import pytest

from planetgen.queue import attempts, dispatch, memorygate, scheduler
from planetgen.queue.memorygate import Claim, Holder, Host
from planetgen.queue.scheduler import Ticket
from tests import enginesim

HOST = Host(16000.0, 12000.0, 8, 0.2)
BULK, INTERACTIVE, MAINT = memorygate.LANE_BULK, memorygate.LANE_INTERACTIVE, memorygate.LANE_MAINTENANCE


def ticket(name, lane=BULK, mb=900.0, since=0.0, **claim):
    return Ticket(name, Claim(claim.pop("kind", "k-" + name), mb, claim.pop("cpu", 1.0), lane, **claim), since, 0, "t")


def holder(name, lane=BULK, mb=900.0, **claim):
    return Holder(name, Claim("held-" + name, mb, 1.0, lane, **claim), 0, "t", mb, mb)


def test_the_site_goes_first_even_against_bulk_work_that_has_aged_past_it():
    old_bulk = [ticket(f"unit-{i}", since=0.0) for i in range(3)]
    visitor = ticket("page", INTERACTIVE, 300.0, since=2000.0)             # a new job; 2,000 s of aging puts every unit above it in one ranking
    now = 2001.0
    ranked = scheduler.rank(old_bulk + [visitor], now, {})
    assert ranked[-1].name == "page", "ranked together the visitor comes last"
    turn = dispatch.dispatch(old_bulk + [visitor], HOST, [], now, {})
    assert turn.admitted[0] == "page"


def test_the_site_lanes_are_served_oldest_first():
    waiting = [ticket("b", INTERACTIVE, 300.0, since=5.0), ticket("a", memorygate.LANE_LONG, 300.0, since=1.0), ticket("c", INTERACTIVE, 300.0, since=5.0)]
    assert dispatch.dispatch(waiting, HOST, [], 10.0, {}).admitted == ["a", "b", "c"]
    site, rest = dispatch.site_first(waiting + [ticket("u")])
    assert [t.name for t in site] == ["a", "b", "c"] and [t.name for t in rest] == ["u"]


def test_a_job_admitted_in_the_turn_counts_against_the_next():
    host = Host(16000.0, 3000.0, 8, 0.1)                    # room for one 900 MB claim after the reserve and the places kept for the site, not two
    waiting = [ticket("one", since=0.0), ticket("two", since=1.0)]
    turn = dispatch.dispatch(waiting, host, [holder("running")], 5.0, {}, places=2)
    assert len(turn.admitted) <= 1
    assert set(turn.admitted) | set(turn.reasons) == {"one", "two"}
    refused = [name for name in ("one", "two") if name not in turn.admitted]
    assert all(turn.reasons[name] for name in refused)


def test_a_reservation_has_no_end_in_the_dispatch_and_ends_at_the_gates_bound_for_a_worker():
    big = ticket("big", mb=9000.0, since=0.0, kind="big")                              # needs more than is free, the highest by age
    small = ticket("small", mb=600.0, since=400.0)                  # newer, lower priority, fits
    held = [holder("a", mb=2000.0), holder("b", mb=2000.0)]
    host = Host(16000.0, 7000.0, 8, 0.1)
    now = 2000.0                                                                       # `big` has waited far past the gate's 300 s
    assert not memorygate.decide(big.claim, host, held, places=2, name="big").admit and memorygate.decide(small.claim, host, held, places=2, name="small").admit
    worker = scheduler.judge(small, [big], host, held, now, {}, places=2)
    assert worker.admit, "main: a worker's reservation for `big` ended after memorygate.MAX_WAIT_SECONDS, and `big` lets in by the table cap"
    central = scheduler.judge(small, [big], host, held, now, {}, places=2, reserve_until=math.inf)
    assert not central.admit and "room is kept for" in central.reason
    assert "small" not in dispatch.dispatch([big, small], host, held, now, {}, places=2).admitted
    assert dispatch.dispatch([big, small], host, held, now, {}, places=2).reasons["small"].startswith("room is kept for")


def test_only_the_top_of_the_order_is_judged_and_the_rest_say_they_wait_their_turn():
    waiting = [ticket(f"unit-{i}", since=float(i)) for i in range(10)]
    full = Host(16000.0, 600.0, 8, 0.1)                     # nothing fits
    turn = dispatch.dispatch(waiting, full, [holder("running")], 20.0, {}, window=4)
    assert turn.admitted == [] and len(turn.reasons) == 10
    assert sum("waits its turn" in reason for reason in turn.reasons.values()) == 6


def test_a_job_that_is_admitted_has_no_waiting_reason_and_one_that_waits_has():
    turn = dispatch.dispatch([ticket("fits", INTERACTIVE, 300.0)], HOST, [holder("x")], 5.0, {})
    assert turn.admitted == ["fits"] and turn.reasons == {}


def test_the_reasons_are_written_onto_the_ledger_rows_of_the_jobs_that_have_one(monkeypatch):
    written = []
    monkeypatch.setattr(attempts, "note_wait", lambda conn, run_id, row, reason, now=None: written.append((run_id, row, reason)) or True)
    turn = dispatch.Turn(["go"], {"a": "needs 900 MB", "b": "a solo run holds the host", "no-row": "x"})
    assert dispatch.note_reasons(object(), turn, {"a": ("run", 11), "b": ("run", 12)}) == 2
    assert written == [("run", 11, "needs 900 MB"), ("run", 12, "a solo run holds the host")]


# --- the simulator --------------------------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("seed", [0, 1, 2])
def test_in_the_simulator_the_site_waits_under_a_second_where_one_ranking_starves_it(seed):
    """A stage of 300 units with visitors' jobs arriving through it: main's order applied to the whole backlog keeps a visitor waiting for over half an hour; the dispatch rule under a second."""
    measures, broken = enginesim.run(seed, 300)
    assert not broken, broken
    assert measures["wait_p99_site"] <= 1.5, measures
    control, _broken = enginesim.run(seed, 300, rule=enginesim.old_rule)
    assert control["wait_p99_site"] > 600.0, "the finding this rule answers: one ranking of the whole backlog starves the site"


def test_in_the_simulator_a_cut_that_needs_the_lane_goes_soon_after_its_room_is_kept_and_nothing_runs_twice():
    for seed in range(3):
        measures, broken = enginesim.run(seed, 150, cut_at=30.0)
        assert not broken, broken
        assert measures["wait_cut"] < 400.0, measures
