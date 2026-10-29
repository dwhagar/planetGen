"""PERF.182: carrying on a generation run that stopped partway (`generation/finish_run.py`). This file drives the wrapper with
commands of the test's own: which stages the record says finished, what `plan` carries over and runs again, and why it
refuses. `test_finish_run_generation.py` stops real runs of the generator and checks that the run carried on ends with the
galaxy an uninterrupted run makes. The autouse `_run_queue_inline` fixture gives each test a queue of its own and runs the
waiting units in this process."""

import pytest

from planetgen.generation import finish_run
from planetgen.queue import commands, redisqueue, runrecord, units, wrapper

CALLS = []
"""list: What the test commands did, in order."""

FAULTS = {}
"""dict: Stage name -> the unit that raises while it is there (`_unit`), until the test clears it: the fault that stops a run."""

SEED, BUILD, DATABASE = "AB" * 16, "77F01600909B06030C0300", "galaxy"


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture(autouse=True)
def _calls():
    CALLS.clear()
    FAULTS.clear()


def _unit(args, report):
    CALLS.append(("unit", args["stage"], args["n"]))
    if FAULTS.get(args["stage"]) == args["n"]:
        raise RuntimeError(f"{args['stage']} stopped")
    return {"n": args["n"]}


def _cut(params):
    return commands.Expansion([({"stage": params["name"], "n": n}, 1.0) for n in range(params["count"])],
                              {"name": params["name"]})


def _prepare(params):
    CALLS.append(("prepare", params["name"]))


def _conclude(params, plan, results, seconds):
    CALLS.append(("conclude", params["name"]))
    return {"stage": params["name"], "units": len(results)}


@pytest.fixture(autouse=True)
def _commands(monkeypatch):
    for kind in ("hooked", "scatter-mass", "scatter-luminosity", "scatter-phenomena", "finish-phenomena"):
        monkeypatch.setitem(commands.REGISTRY, kind, commands.CommandSpec(
            kind, "A test", unit=(__name__, "_unit"), expand=(__name__, "_cut"), prepare=(__name__, "_prepare"),
            conclude=(__name__, "_conclude")))
    monkeypatch.setitem(commands.REGISTRY, "plain", commands.CommandSpec(
        "plain", "A test", unit=(__name__, "_unit"), expand=(__name__, "_cut")))
    monkeypatch.setattr(finish_run, "FINISHABLE", finish_run.FINISHABLE | {"hooked", "plain"})


def _step(name, kind="hooked", count=2, **params):
    return commands.Step(name, kind, name.capitalize(), {"name": name, "count": count, **params})


def _run(steps, galaxy_seed=SEED, build=BUILD, finishes=None, database=DATABASE, **kwargs):
    return wrapper.run("Stages", steps, kwargs.pop("on_plan", lambda index, plan: None),
                       kwargs.pop("on_result", lambda index, unit, result: None), database=database,
                       galaxy_seed=galaxy_seed, build=build, finishes=finishes, **kwargs)


def _newest(conn):
    """The run that started last: the ids sort by the second they started in, and two runs of a test can share one."""
    return max(runrecord.recent_runs(conn, limit=50), key=lambda run_id: runrecord.read_run(conn, run_id)["created"])


def _stopped(conn, steps, stage, unit, **kwargs):
    """The id of the run of `steps` that stopped in `stage` at unit `unit` (the fault is gone after it)."""
    FAULTS[stage] = unit
    with pytest.raises(units.RunStopped):
        _run(steps, **kwargs)
    FAULTS.clear()
    CALLS.clear()
    return _newest(conn)


def _three(conn, stage="second", unit=1):
    return _stopped(conn, [_step("first"), _step("second"), _step("third")], stage, unit)


def _plan(conn, run_id, galaxy_seed=SEED, build=BUILD, database=DATABASE):
    return finish_run.plan(conn, run_id, galaxy_seed, build, database)


def test_a_run_stopped_in_its_second_stage_is_carried_on_from_there(conn):
    run_id = _three(conn)
    found = _plan(conn, run_id)
    assert (found.first, found.stopped) == (1, 1)
    assert [step.params.get(wrapper.FINISHED_BY) for step in found.steps] == [run_id, None, None]
    assert [step.params.get(finish_run.AGAIN) for step in found.steps] == [None, True, None]
    new_id = _run(found.steps, finishes=run_id)
    assert CALLS == [("prepare", "second"), ("unit", "second", 0), ("unit", "second", 1), ("conclude", "second"),
                     ("prepare", "third"), ("unit", "third", 0), ("unit", "third", 1), ("conclude", "third")], CALLS
    new = runrecord.read_run(conn, new_id)
    assert new["state"] == "done" and new["finishes"] == run_id and new["galaxy_seed"] == SEED and new["build"] == BUILD
    assert [stage["units"] for stage in new["stages"]] == [0, 2, 2]


def test_a_stage_the_run_before_finished_is_carried_over_with_its_conclusion(conn):
    run_id = _three(conn)
    planned, concluded = {}, {}
    new_id = _run(_plan(conn, run_id).steps, finishes=run_id, on_plan=lambda index, plan: planned.update({index: plan}),
                  on_stage_done=lambda index, conclusion: concluded.update({index: conclusion}))
    assert planned[0]["units"] == 0
    assert planned[0]["detail"] == {"skip": f"finished by run {run_id}", wrapper.FINISHED_BY: run_id}
    assert concluded[0] == {"stage": "first", "units": 2} and concluded[1] == {"stage": "second", "units": 2}
    assert not [call for call in CALLS if call[1] == "first"], "a stage carried over is neither prepared, run nor concluded"
    assert runrecord.read_run(conn, new_id)["state"] == "done"


def test_a_finish_that_stops_is_finished_in_turn(conn):
    run_id = _three(conn)
    third = _stopped(conn, _plan(conn, run_id).steps, "third", 0, finishes=run_id)
    assert third != run_id
    found = _plan(conn, third)
    assert (found.first, found.stopped) == (2, 2)
    assert [step.params.get(wrapper.FINISHED_BY) for step in found.steps] == [third, third, None]
    assert [step.params.get(finish_run.AGAIN) for step in found.steps] == [None, None, True]
    done = _run(found.steps, finishes=third)
    assert CALLS == [("prepare", "third"), ("unit", "third", 0), ("unit", "third", 1), ("conclude", "third")]
    assert runrecord.read_run(conn, done)["state"] == "done"


def test_a_stage_cut_to_no_unit_is_finished_when_the_next_one_began(conn):
    run_id = _stopped(conn, [_step("first", count=0), _step("second")], "second", 0)
    queued = [commands.Step.from_dict(item) for item in runrecord.read_steps(conn, run_id)]
    assert finish_run.stage_finished(conn, run_id, 0, queued[0])
    assert not finish_run.stage_finished(conn, run_id, 1, queued[1])
    assert finish_run.first_unfinished(conn, run_id, queued) == 1


def test_a_stage_whose_units_are_done_but_not_concluded_is_not_finished(conn):
    run_id = _three(conn)
    queued = [commands.Step.from_dict(item) for item in runrecord.read_steps(conn, run_id)]
    assert finish_run.stage_finished(conn, run_id, 0, queued[0])
    conn.hdel(runrecord.stage_key(run_id, 0), "concluded")
    conn.hdel(runrecord.stage_key(run_id, 1), "started")         # the next stage never began
    assert not finish_run.stage_finished(conn, run_id, 0, queued[0])
    assert finish_run.first_unfinished(conn, run_id, queued) == 0


def test_the_special_phenomena_are_rewound_with_the_scatter_they_were_written_into(conn):
    steps = [_step("start", "plain"), _step("phenomena", "scatter-phenomena"), _step("phenomena-index", "finish-phenomena"),
             _step("blocks", "plain")]
    run_id = _stopped(conn, steps, "phenomena-index", 0)
    found = _plan(conn, run_id)
    assert (found.first, found.stopped) == (1, 2)
    assert [step.params.get(wrapper.FINISHED_BY) for step in found.steps] == [run_id, None, None, None]
    assert [step.params.get(finish_run.AGAIN) for step in found.steps] == [None, None, True, None]


def test_the_bright_stars_of_a_run_that_draws_both_passes_are_rewound_together(conn):
    run_id = _stopped(conn, [_step("mass", "scatter-mass", clear="all"), _step("luminosity", "scatter-luminosity")],
                      "luminosity", 0)
    found = _plan(conn, run_id)
    assert (found.first, found.stopped) == (0, 1)


def test_the_bright_stars_of_a_run_that_clears_its_own_pass_are_not(conn):
    steps = [_step("mass", "scatter-mass", clear="all"), _step("luminosity", "scatter-luminosity", clear="luminosity")]
    run_id = _stopped(conn, steps, "luminosity", 0)
    assert _plan(conn, run_id).first == 1


def test_a_run_already_carried_on_is_not_carried_on_twice(conn):
    run_id = _three(conn)
    carried = _stopped(conn, _plan(conn, run_id).steps, "third", 0, finishes=run_id)
    with pytest.raises(finish_run.CannotFinish, match=f"already carried on by run {carried}"):
        _plan(conn, run_id)
    assert _plan(conn, carried).first == 2


def test_the_page_leaves_the_build_unchecked(conn):
    run_id = _three(conn)
    assert _plan(conn, run_id, build=None).first == 1


@pytest.mark.parametrize("change, message", [
    ({"galaxy_seed": "CD" * 16}, "planned again"),
    ({"build": "FFFFFFFFFFFFFFFFFFFFFF"}, "PlanetGen changed"),
    ({"database": "elsewhere"}, "database"),
])
def test_a_run_is_not_carried_on_on_another_galaxy_or_build_or_database(conn, change, message):
    run_id = _three(conn)
    with pytest.raises(finish_run.CannotFinish, match=message):
        _plan(conn, run_id, **change)


def test_a_run_that_recorded_no_galaxy_is_not_carried_on(conn):
    run_id = _stopped(conn, [_step("first")], "first", 0, galaxy_seed=None, build=None)
    with pytest.raises(finish_run.CannotFinish, match="did not record"):
        _plan(conn, run_id)


def test_a_stage_with_no_way_to_carry_on_from_its_middle_is_refused(conn, monkeypatch):
    run_id = _stopped(conn, [_step("first"), _step("second", "hooked")], "second", 0)
    monkeypatch.setattr(finish_run, "FINISHABLE", frozenset({"plain"}))
    with pytest.raises(finish_run.CannotFinish, match="Second.*cannot be carried on from the middle"):
        _plan(conn, run_id)


def test_a_stage_the_run_finished_may_be_one_with_no_way_to_carry_on(conn, monkeypatch):
    run_id = _stopped(conn, [_step("first", "hooked"), _step("second", "plain")], "second", 0)
    monkeypatch.setattr(finish_run, "FINISHABLE", frozenset({"plain"}))
    assert _plan(conn, run_id).first == 1


def test_a_run_that_was_done_has_nothing_to_finish(conn):
    run_id = _run([_step("first")])
    with pytest.raises(finish_run.CannotFinish, match="nothing to finish"):
        _plan(conn, run_id)
    assert finish_run.last_incomplete(conn) is None


def test_the_last_incomplete_run_is_the_newest_generation_run_when_it_stopped(conn):
    first = _run([_step("first")])
    assert finish_run.last_incomplete(conn) is None
    stopped = _three(conn)
    assert stopped != first and finish_run.last_incomplete(conn) == stopped
    _run([_step("again")])
    assert finish_run.last_incomplete(conn) is None


def test_a_run_whose_process_is_gone_is_taken_for_stopped(conn):
    step = _step("first")
    run_id = runrecord.start_run(conn, "Gone", [runrecord.Stage("first", "First", units=1)], steps=[step.as_dict()],
                                 galaxy_seed=SEED, build=BUILD, database=DATABASE)
    runrecord.begin_stage(conn, run_id, 0)
    conn.delete(runrecord.alive_key(run_id))
    conn.hset(runrecord.run_key(run_id), "owner", "elsewhere:1")
    assert finish_run.last_incomplete(conn) == run_id
    assert runrecord.read_run(conn, run_id)["state"] == "failed"
    found = _plan(conn, run_id)
    assert found.first == 0 and found.steps[0].params[finish_run.AGAIN] is True
    assert runrecord.active_run(conn) is None, "the dead run no longer holds the machine"


def test_a_run_still_running_is_not_carried_on(conn):
    step = _step("first")
    run_id = runrecord.start_run(conn, "Live", [runrecord.Stage("first", "First", units=1)], steps=[step.as_dict()],
                                 galaxy_seed=SEED, build=BUILD, database=DATABASE)
    try:
        with pytest.raises(finish_run.CannotFinish, match="still running"):
            finish_run.last_incomplete(conn)
        with pytest.raises(finish_run.CannotFinish, match="still running"):
            _plan(conn, run_id)
    finally:
        runrecord.finish_run(conn, run_id, "cancelled")
