"""PERF.172: the work of a stage that is not a unit -- a command's `prepare` when the stage opens and `conclude` when its last
unit is done -- and what happens to a run when one of them, or the opening of a stage, fails. The autouse
`_run_queue_inline` fixture gives each test a queue of its own and runs the waiting units in this process."""

import pytest

from planetgen.queue import commands, redisqueue, runrecord, units, wrapper

CALLS = []
"""list: What the test commands did, in order."""


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture(autouse=True)
def _calls():
    CALLS.clear()


def _unit(args, report):
    CALLS.append(("unit", args["stage"], args["n"]))
    return {"n": args["n"]}


def _cut(params):
    if params.get("explode") == "expand":
        raise RuntimeError("the galaxy cannot be read")
    return commands.Expansion([({"stage": params["name"], "n": n}, 1.0) for n in range(params["count"])], {"name": params["name"]})


def _prepare(params):
    CALLS.append(("prepare", params["name"]))
    if params.get("explode") == "prepare":
        raise RuntimeError("cannot empty the table")


def _conclude(params, plan, results, seconds):
    CALLS.append(("conclude", params["name"], [result["n"] for result in results]))
    if params.get("explode") == "conclude":
        raise RuntimeError("cannot store the record")
    return {"stage": params["name"], "units": len(results)}


@pytest.fixture(autouse=True)
def _commands(monkeypatch):
    monkeypatch.setitem(commands.REGISTRY, "hooked", commands.CommandSpec(
        "hooked", "A test", unit=(__name__, "_unit"), expand=(__name__, "_cut"), prepare=(__name__, "_prepare"),
        conclude=(__name__, "_conclude")))
    monkeypatch.setitem(commands.REGISTRY, "plain", commands.CommandSpec(
        "plain", "A test", unit=(__name__, "_unit"), expand=(__name__, "_cut")))


def _step(name, kind="hooked", count=2, **params):
    return commands.Step(name, kind, name.capitalize(), {"name": name, "count": count, **params})


def _run(steps, done=None):
    planned, concluded = [], {}

    def on_stage_done(index, conclusion):
        concluded[index] = conclusion
        CALLS.append(("done", index))

    run_id = wrapper.run("Hooks", steps, lambda index, plan: planned.append(index), lambda index, unit, result: None,
                         on_stage_done=on_stage_done)
    return run_id, planned, concluded


def test_a_stage_is_prepared_before_it_is_cut_and_concluded_before_the_next_is_opened(conn):
    run_id, planned, concluded = _run([_step("first"), _step("second")])
    assert CALLS.index(("prepare", "first")) < CALLS.index(("unit", "first", 0))
    assert CALLS.index(("conclude", "first", [0, 1])) < CALLS.index(("prepare", "second"))
    assert planned == [0, 1] and concluded == {0: {"stage": "first", "units": 2}, 1: {"stage": "second", "units": 2}}
    assert runrecord.read_run(conn, run_id)["state"] == "done"


def test_a_command_with_no_hooks_has_no_conclusion(conn):
    _run_id, _planned, concluded = _run([_step("only", kind="plain")])
    assert concluded == {0: None}
    assert [call for call in CALLS if call[0] in ("prepare", "conclude")] == []


def test_a_stage_that_cannot_be_cut_fails_the_run_at_once_naming_the_stage(conn):
    with pytest.raises(units.RunStopped, match=r"stage Second could not start: RuntimeError: the galaxy cannot be read"):
        _run([_step("first"), _step("second", explode="expand")])
    run = runrecord.read_run(conn, runrecord.recent_runs(conn)[0])
    assert run["state"] == "failed" and run["error"].startswith("stage Second could not start")
    assert runrecord.active_run(conn) is None, "the machine is free"


def test_a_prepare_that_fails_fails_the_run_naming_the_stage(conn):
    with pytest.raises(units.RunStopped, match="stage Second could not start: RuntimeError: cannot empty the table"):
        _run([_step("first"), _step("second", explode="prepare")])
    assert ("unit", "second", 0) not in CALLS, "nothing of the stage ran"


def test_a_conclusion_that_fails_fails_the_run_before_the_next_stage_opens(conn):
    with pytest.raises(units.RunStopped, match="the end of stage First failed: RuntimeError: cannot store the record"):
        _run([_step("first", explode="conclude"), _step("second")])
    assert ("prepare", "second") not in CALLS


def test_a_prepare_that_fails_in_the_first_stage_raises_in_the_starter(conn):
    with pytest.raises(RuntimeError, match="cannot empty the table"):
        _run([_step("first", explode="prepare")])
    assert runrecord.read_run(conn, runrecord.recent_runs(conn)[0])["state"] == "failed"


def test_a_stage_cut_into_no_unit_is_still_prepared_and_concluded(conn):
    run_id, planned, concluded = _run([_step("empty", count=0), _step("after")])
    assert [call for call in CALLS if call[0] in ("prepare", "conclude")] == [
        ("prepare", "empty"), ("conclude", "empty", []), ("prepare", "after"), ("conclude", "after", [0, 1])]
    assert ("done", 0) in CALLS and ("done", 1) in CALLS
    assert concluded[0] == {"stage": "empty", "units": 0} and runrecord.read_run(conn, run_id)["state"] == "done"


def test_the_plan_of_a_stage_with_nothing_to_do_stores_its_row_as_skipped(conn, monkeypatch):
    stored = {}
    from planetgen.generation import stages
    monkeypatch.setattr(stages, "write_skipped_row", lambda config, key, row, reason: stored.update(key=key, reason=reason))
    monkeypatch.setattr(stages, "write_queued_row", lambda *args: stored.update(ran=True))
    monkeypatch.setattr(wrapper, "_mysql_config", lambda named: None)
    monkeypatch.setitem(commands.REGISTRY, "idle", commands.CommandSpec(
        "idle", "A test", unit=(__name__, "_unit"), expand=(__name__, "_cut_nothing")))
    row = {"command": "test", "n": 1, "total": 1, "label": "Idle", "settings": {}}
    step = commands.Step("idle", "idle", "Idle", {"name": "idle", "database": {"database": "d"}})
    wrapper.run("Skipped", [step], lambda index, plan: None, lambda index, unit, result: None, rows={"idle": row})
    assert stored == {"key": "idle", "reason": "nothing was saved"}


def _cut_nothing(params):
    return commands.Expansion([({"stage": "idle", "n": 0}, 1.0)], {"skip": "nothing was saved"})
