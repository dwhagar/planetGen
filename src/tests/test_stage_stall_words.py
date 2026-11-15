"""OPS.253: a stage that finishes no unit for hours says what it waits on, in the run line, in the log and on its bar.

Boss's New galaxy run sat 4 h 17 m in the bright-star backfill with "waiting on progress" and nothing more: the five sky cells draw one at a time, the others wait on the lock, and a wait is not stuck
(`beatjudge`), so no line said so, and the stage's bar only moved when a whole unit ended.
"""

import json

import pytest

from planetgen.generation import run_galaxy
from planetgen.queue import redisqueue, runrecord, runview, units
from planetgen.web import queue_page


def test_the_cause_names_what_the_running_units_do_the_largest_group_first():
    views = {n: {"state": "waiting", "waiting_on": "sky lock"} for n in range(4)}
    views[4] = {"state": "working", "waiting_on": ""}
    assert runview.stall_cause(views) == "4 waiting on sky lock, 1 working"
    assert runview.stall_cause({}) == ""


def test_a_stalled_estimate_shows_its_cause_and_a_moving_one_does_not():
    assert queue_page._eta_text({"seconds": 5.0, "stalled": True, "cause": "4 waiting on sky lock, 1 working"}, True) == "waiting on progress (4 waiting on sky lock, 1 working)"
    assert queue_page._eta_text({"seconds": 5.0, "stalled": True}, True) == "waiting on progress"


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


def test_a_stage_with_no_unit_finished_says_so_every_half_hour_and_names_what_the_others_do(conn, monkeypatch):
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=2)])
    runrecord.begin_stage(conn, run_id, 0)
    key = runrecord.unitstate_key(run_id, 0)
    conn.hset(key, 0, json.dumps({"state": "waiting", "waiting_on": "sky lock"}))
    conn.hset(key, 1, json.dumps({"state": "working", "waiting_on": ""}))
    now = [0.0]
    polls = [0]
    said = []
    monkeypatch.setattr(units.log, "normal", lambda message, *a, **k: said.append(message))

    def sleep(_seconds):
        now[0] += 700.0
        polls[0] += 1
        if polls[0] == 8:                           # 5,600 s in: both units end
            conn.sadd(runrecord.done_key(run_id, 0), 0, 1)
            conn.hset(runrecord.results_key(run_id, 0), mapping={0: "{}", 1: "{}"})

    units.follow_stage(conn, run_id, 0, 2, lambda unit, result: None, clock=lambda: now[0], sleep=sleep)
    assert len(said) == 2 and said[0].startswith("Stage 1: no unit has finished for 35 m; 0 of 2 units are done, the others: ")
    assert "1 waiting on sky lock" in said[0] and "1 working" in said[0] and "1 h 10 m" in said[1]


class _Bar:
    def __init__(self):
        self.total, self.done = 0.0, 0.0

    def update(self, advance=None, total=None):
        self.total = total if total is not None else self.total
        self.done += advance or 0.0


def _follower(sky_costs, chunks):
    follower = run_galaxy.BackfillFollower(None, "2026-11-14")
    follower.hold = lambda context: context
    run_galaxy.steps.Step = lambda *a, **k: _Bar()
    follower.on_plan({"units": chunks + len(sky_costs), "detail": {"sectors": 200, "sky_cells": len(sky_costs), "sky_costs": sky_costs}})
    return follower


def test_the_backfill_bar_counts_a_sky_cell_by_its_own_counter_and_never_goes_back(monkeypatch):
    saved = run_galaxy.steps.Step
    try:
        follower = _follower([1000.0, 3000.0], chunks=2)
        bar = follower.bar
        assert bar.total == 4200.0
        follower.on_progress(0, 50.0, 100.0)             # a sector chunk: 50 of its 100 sectors
        follower.on_progress(2, 25.0, 100.0)             # a sky cell a quarter through: a quarter of 1000
        assert bar.done == 300.0
        follower.on_progress(2, 10.0, 100.0)             # a counter that reads lower does not move the bar back
        assert bar.done == 300.0
        follower.on_result(2, {"sectors": 0, "stars": 5, "visited": 0, "shells": {}})     # the cell ends: its whole cost
        follower.on_result(0, {"sectors": 10, "stars": 1, "visited": 100, "shells": {}})
        assert bar.done == 1000.0 + 100.0
    finally:
        run_galaxy.steps.Step = saved


@pytest.mark.parametrize("route", ["a backfill run of its own (the API's sectors, planetgen backfill)", "a stage of a whole Generate job (planetgen run-job)"])
def test_every_way_a_run_gets_its_backfill_gets_the_bar_that_counts_a_sky_cell(route):
    """OPS.253: `planetgen galaxy` lists `BackfillFollower` itself (run_galaxy `tail`), the backfill run of its own is `planned_backfill`, and a Generate job's stage is `job._after_sectors("backfill")`: one class,
    so the bar is the same on every route, and a cell with a sky draw counts by its own counter on each."""
    import types
    from planetgen.generation import job
    args = types.SimpleNamespace(sky_limit=None)
    if route.startswith("a backfill run"):
        follower = run_galaxy.planned_backfill(args, "2026-11-14 00:00:00").followers[0]
    else:
        follower = job._after_sectors("backfill")(args, 0, "2026-11-14 00:00:00").followers[0]
    assert isinstance(follower, run_galaxy.BackfillFollower)
    saved, run_galaxy.steps.Step = run_galaxy.steps.Step, (lambda *a, **k: _Bar())
    try:
        follower.hold = lambda context: context
        follower.on_plan({"units": 3, "detail": {"sectors": 100, "sky_cells": 1, "sky_costs": [900.0]}})
        follower.on_progress(2, 50.0, 100.0)
        assert follower.bar.total == 1000.0 and follower.bar.done == 450.0
    finally:
        run_galaxy.steps.Step = saved


def test_the_run_line_and_the_generate_page_both_read_the_cause_of_a_stalled_stage(conn, monkeypatch):
    """The Queue page's list and detail (`api/workqueue.py`) read it by `summary` with `cause_of_stalled`; the Generate page's status reads it by `runview.read`: the same words on both."""
    from planetgen.queue import runeta
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("test", "T", units=2)])
    runrecord.begin_stage(conn, run_id, 0)
    conn.hset(runrecord.unitstate_key(run_id, 0), 0, json.dumps({"state": "waiting", "waiting_on": "sky lock"}))
    monkeypatch.setattr(runeta, "estimate", lambda run, now: runeta.Eta(seconds=5.0, low=1.0, high=9.0, stalled=True))
    run = runrecord.read_run(conn, run_id)
    by_list = runview.summary(run, 100.0, runview.cause_of_stalled(conn, run, 100.0))["eta"]
    by_read = runview.read(conn, run_id, 100.0)["eta"]
    assert by_list["cause"] == by_read["cause"] == "1 waiting on sky lock"
    assert queue_page._eta_text(by_list, True) == "waiting on progress (1 waiting on sky lock)"
    monkeypatch.setattr(runeta, "estimate", lambda run, now: runeta.Eta(seconds=5.0, low=1.0, high=9.0, stalled=False))
    assert "cause" not in runview.read(conn, run_id, 100.0)["eta"]
