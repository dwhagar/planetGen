"""PERF.213: a step parameter that is worked out when its stage opens (`wrapper.opened_by`) and not when the run is built. This
file drives the wrapper with commands of the test's own, like `test_finish_run.py`: the clock a stage records when it opens, the
stage that takes it as its parameter, what `finish_run` keeps of it, and the followers of the stages that work on a run's sectors.
The autouse `_run_queue_inline` fixture gives each test a queue of its own and runs the waiting units in this process."""

import datetime
import json

import pytest

from planetgen.generation import finish_run, run_galaxy
from planetgen.queue import commands, redisqueue, runrecord, units, wrapper

CALLS = []
"""list: What the test commands did, in order: `("cut", stage name, the `since` it was given)`."""

FAULTS = {}
"""dict: Stage name -> the unit that raises while it is there."""

SEED, BUILD, DATABASE = "AB" * 16, "77F01600909B06030C0300", "galaxy"


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture(autouse=True)
def _calls(monkeypatch):
    CALLS.clear()
    FAULTS.clear()
    clocks = iter(f"2026-10-30 06:00:{second:02d}" for second in range(1, 60))
    monkeypatch.setattr(wrapper, "_database_clock", lambda params: next(clocks))


def _unit(args, report):
    if FAULTS.get(args["stage"]) == args["n"]:
        raise RuntimeError(f"{args['stage']} stopped")
    return {"n": args["n"]}


def _cut(params):
    CALLS.append(("cut", params["name"], params.get("since"), params.get("centers")))
    detail = {"name": params["name"]}
    if params.get("planned") == "random-start":
        detail["centers"] = [[4, 5, 6], [7, 8, 9]]
    return commands.Expansion([({"stage": params["name"], "n": n}, 1.0) for n in range(params["count"])], detail)


@pytest.fixture(autouse=True)
def _commands(monkeypatch):
    monkeypatch.setitem(commands.REGISTRY, "plain", commands.CommandSpec("plain", "A test", unit=(__name__, "_unit"), expand=(__name__, "_cut")))
    monkeypatch.setattr(finish_run, "FINISHABLE", finish_run.FINISHABLE | {"plain"})


def _step(name, count=2, **params):
    return commands.Step(name, "plain", name.capitalize(), {"name": name, "count": count, **params})


def _run(steps, finishes=None, **kwargs):
    return wrapper.run("Stages", steps, lambda index, plan: None, lambda index, unit, result: None, database=DATABASE,
                       galaxy_seed=SEED, build=BUILD, finishes=finishes, **kwargs)


def _newest(conn):
    return max(runrecord.recent_runs(conn, limit=50), key=lambda run_id: runrecord.read_run(conn, run_id)["created"])


def _since_of(name):
    return next(call[2] for call in CALLS if call[1] == name)


def test_a_stage_takes_the_clock_of_the_stage_it_names_when_it_opens(conn):
    steps = [_step("sectors"), _step("middle"), _step("link", since=wrapper.opened_by("sectors"))]
    run_id = _run(steps)
    assert _since_of("link") == "2026-10-30 06:00:01" and _since_of("sectors") is None
    assert runrecord.opened_clock(conn, run_id, 0) == "2026-10-30 06:00:01"


def test_only_a_stage_some_step_refers_to_records_a_clock(conn):
    run_id = _run([_step("sectors"), _step("middle"), _step("link", since=wrapper.opened_by("sectors"))])
    assert [runrecord.opened_clock(conn, run_id, index) for index in range(3)] == ["2026-10-30 06:00:01", None, None]
    plain = _run([_step("a"), _step("b")])
    assert [runrecord.opened_clock(conn, plain, index) for index in range(2)] == [None, None]


def test_a_second_stage_that_refers_to_the_same_one_gets_the_same_clock(conn):
    _run([_step("sectors"), _step("link", since=wrapper.opened_by("sectors")), _step("settle", since=wrapper.opened_by("sectors"))])
    assert _since_of("link") == _since_of("settle") == "2026-10-30 06:00:01"


def test_the_recorded_step_keeps_the_deferred_parameter(conn):
    run_id = _run([_step("sectors"), _step("link", since=wrapper.opened_by("sectors"))])
    link = next(item for item in runrecord.read_steps(conn, run_id) if item["key"] == "link")
    assert link["params"]["since"] == wrapper.opened_by("sectors")


def test_a_deferred_parameter_that_names_no_stage_of_the_run_is_an_error(conn):
    with pytest.raises(wrapper.Unresolved):
        wrapper.resolve_params(conn, _run([_step("a")]), {"since": wrapper.opened_by("nowhere")})
    assert wrapper.resolve_params(conn, _run([_step("a")]), {"since": wrapper.opened_by("nowhere")}, strict=False) == {
        "since": wrapper.opened_by("nowhere")}


def test_a_stage_that_has_not_opened_has_no_clock_to_take(conn):
    run_id = _run([_step("a"), _step("b", since=wrapper.opened_by("a"))])
    with pytest.raises(wrapper.Unresolved):
        wrapper.resolve_params(conn, run_id, {"since": wrapper.opened_by("b")})


def test_scatter_settings_of_a_stage_that_has_not_opened_are_left_for_when_it_opens_when_not_strict(conn):
    """PERF.203 step 5: finishing a run plans its stages from the record with `strict=False`; settings that only the scatter's own
    stage works out (from the galaxy as the stages before it leave it) must not be worked out then."""
    marker = {"deferred": "scatter-settings", "scope": "bright-stars", "stage": "mass"}
    run_id = _run([_step("a")])
    assert wrapper.resolve_params(conn, run_id, {"settings": marker}, strict=False) == {"settings": marker}
    conn.hset(runrecord.run_key(run_id), f"deferred:{wrapper.scatter_settings_name(marker)}", json.dumps({"mass": {"limit": 8.0}}))
    assert wrapper.resolve_params(conn, run_id, {"settings": marker}, strict=False) == {"settings": {"limit": 8.0}}


def test_other_parameters_pass_through_unchanged(conn):
    run_id = _run([_step("a")])
    params = {"since": "2026-01-01 00:00:00", "count": 3, "items": [1, 2], "nested": {"deferred_not": 1}}
    assert wrapper.resolve_params(conn, run_id, params) == params


def _stopped(conn, steps, stage, unit):
    FAULTS[stage] = unit
    with pytest.raises(units.RunStopped):
        _run(steps)
    FAULTS.clear()
    CALLS.clear()
    return _newest(conn)


def test_a_run_carried_on_keeps_the_clock_the_stopped_run_took(conn):
    steps = [_step("sectors"), _step("link", since=wrapper.opened_by("sectors"))]
    run_id = _stopped(conn, steps, "sectors", 1)       # the sectors stage opened, saved some sectors and stopped
    found = finish_run.plan(conn, run_id, SEED, BUILD, DATABASE)
    assert found.steps[1].params["since"] == "2026-10-30 06:00:01"
    _run(found.steps, finishes=run_id)
    assert _since_of("link") == "2026-10-30 06:00:01", "the sectors the first run saved are still the run's"


def test_a_run_carried_on_from_after_the_stage_keeps_its_clock_too(conn):
    steps = [_step("sectors"), _step("link", since=wrapper.opened_by("sectors"))]
    run_id = _stopped(conn, steps, "link", 0)
    found = finish_run.plan(conn, run_id, SEED, BUILD, DATABASE)
    assert found.steps[0].params[wrapper.FINISHED_BY] == run_id
    assert found.steps[1].params["since"] == "2026-10-30 06:00:01"
    _run(found.steps, finishes=run_id)
    assert _since_of("link") == "2026-10-30 06:00:01"


def test_a_stage_that_never_opened_leaves_the_parameter_for_the_run_that_opens_it(conn):
    steps = [_step("first"), _step("sectors"), _step("link", since=wrapper.opened_by("sectors"))]
    run_id = _stopped(conn, steps, "first", 0)
    found = finish_run.plan(conn, run_id, SEED, BUILD, DATABASE)
    assert found.steps[2].params["since"] == wrapper.opened_by("sectors")
    _run(found.steps, finishes=run_id)
    assert _since_of("link") == "2026-10-30 06:00:01"


# --- the followers of the stages that work on a run's sectors --------------------------------------------------------------------

def _args():
    import argparse
    return argparse.Namespace()


@pytest.mark.parametrize("follower", (run_galaxy.LinkFollower, run_galaxy.BackfillFollower, run_galaxy.SettleFollower,
                                      run_galaxy.RunBlocksFollower))
def test_a_follower_hands_the_wrapper_the_deferred_clock_or_the_clock_it_was_given(follower):
    marker = wrapper.opened_by("sectors")
    assert follower(_args(), marker).params()["since"] == marker
    started = datetime.datetime(2026, 10, 30, 6, 0, 1)
    assert follower(_args(), started).params()["since"] == "2026-10-30 06:00:01"


def test_the_galaxy_run_gives_every_stage_that_works_on_its_sectors_the_clock_of_the_first_sectors_stage(monkeypatch):
    seen = {}
    monkeypatch.setattr(run_galaxy, "_run_stages", lambda args, title, tail: seen.update(tail=tail))
    monkeypatch.setattr(run_galaxy.run_plan, "bright_star_followers", lambda args: ([], None))
    monkeypatch.setattr(run_galaxy.run_plan, "phenomena_followers", lambda args: ([], None))
    run_galaxy.run_galaxy(_args())
    marker = wrapper.opened_by("sectors")
    followers = seen["tail"](marker)
    sincers = [f for f in followers if "since" in f.params()]
    assert {f.key for f in sincers} == {"link", "backfill", "settle", "run-blocks"}
    assert all(f.params()["since"] == marker for f in sincers)


# --- a parameter that reads the plan of an earlier stage (PERF.214) --------------------------------------------------------------

def _centers_of(name):
    return next(call[3] for call in CALLS if call[1] == name)


def test_a_stage_reads_what_the_stage_before_it_planned(conn):
    _run([_step("start", planned="random-start"), _step("neighborhood", centers=wrapper.plan_detail("start", "centers"))])
    assert _centers_of("neighborhood") == [[4, 5, 6], [7, 8, 9]]


def test_a_plan_field_the_stage_did_not_plan_is_an_error(conn):
    run_id = _run([_step("start"), _step("neighborhood")])
    with pytest.raises(wrapper.Unresolved):
        wrapper.resolve_params(conn, run_id, {"centers": wrapper.plan_detail("start", "centers")})
    with pytest.raises(wrapper.Unresolved):
        wrapper.resolve_params(conn, run_id, {"centers": wrapper.plan_detail("nowhere", "centers")})


def test_a_run_carried_on_keeps_the_random_starts_the_stopped_run_picked(conn):
    steps = [_step("start", planned="random-start", items=[]), _step("neighborhood", centers=wrapper.plan_detail("start", "centers"))]
    run_id = _stopped(conn, steps, "start", 1)
    found = finish_run.plan(conn, run_id, SEED, BUILD, DATABASE)
    assert found.steps[0].params["planned"] == "" and found.steps[0].params["items"] == [4, 5, 6, 7, 8, 9]
    assert found.steps[1].params["centers"] == [[4, 5, 6], [7, 8, 9]]


def test_a_stage_stopped_before_it_planned_its_starts_plans_them_again(conn):
    steps = [_step("first"), _step("start", planned="random-start", items=[])]
    run_id = _stopped(conn, steps, "first", 0)
    found = finish_run.plan(conn, run_id, SEED, BUILD, DATABASE)
    assert found.steps[1].params["planned"] == "random-start"
