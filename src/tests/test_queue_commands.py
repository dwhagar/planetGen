"""PERF.118, PERF.121, PERF.122: every stage of a `galaxy`, `plan`, `scatter` or fill run is a registered command, and the wrapper (`planetgen/queue/wrapper.py`)
turns the arguments into the ordered steps that will run, writes the run record, queues the first stage and lets the done set
queue the next (`planetgen/queue/commands.py`)."""

import json

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.generation import stages
from planetgen.queue import commands, redisqueue, runrecord, units, wrapper

GALAXY_RUNS = (["galaxy"], ["galaxy", "--shell", "2"], ["sectors"], ["sectors", "--shell", "2"], ["link"], ["backfill"],
               ["settle"])
PLAN_RUNS = (
    ["plan"], ["scatter", "--mass"], ["scatter", "--luminosity"], ["scatter", "--phenomena"],
    ["scatter", "--nebula-field"], ["scatter", "--band", "--down-to", "10"], ["fill-map-luminosity"],
    ["rescore-habitability"], ["clean-nebula-field"], ["clean-nebula-field", "--apply"],
    ["clear-scatter", "--mass"], ["clear-scatter", "--luminosity"], ["clear-scatter", "--phenomena"],
    ["clear-scatter", "--nebula-field"])


def _args(argv):
    parser, _parsers = generate_cli.build_parser()
    return parser.parse_args(argv)


@pytest.mark.parametrize("argv", GALAXY_RUNS + PLAN_RUNS)
def test_every_step_of_a_galaxy_or_plan_run_is_a_registered_command(argv):
    args = _args(argv)
    steps = wrapper.steps_for(args.command, args)
    assert steps, "a run does something"
    for step in steps:
        assert step.kind in commands.REGISTRY, f"{step.key} runs the unregistered command {step.kind}"
    assert [step.key for step in steps] == [stage.key for stage in stages.for_command(args.command, args)
                                            if stage.skip is None], "exactly the stages the request selects"


@pytest.mark.parametrize("argv", GALAXY_RUNS + PLAN_RUNS)
def test_the_same_arguments_give_the_same_list(argv):
    args = _args(argv)
    first = wrapper.steps_for(args.command, args)
    assert first == wrapper.steps_for(args.command, _args(argv))
    assert [commands.Step.from_dict(json.loads(json.dumps(step.as_dict()))) for step in first] == first


def test_a_step_the_request_does_not_select_is_not_in_the_list():
    kinds = lambda argv: [step.kind for step in wrapper.steps_for(*(lambda a: (a.command, a))(_args(argv)))]  # noqa: E731
    assert "scatter-mass" in kinds(["galaxy"]) and "scatter-mass" not in kinds(["sectors"])
    assert kinds(["galaxy", "--shell", "2"])[0] == "generate-sectors" and kinds(["link"]) == ["link-neighbours"]
    assert kinds(["plan"]) == ["plan-skeleton"]
    assert kinds(["scatter", "--mass"]) == ["scatter-mass", "refresh-block-lights"]
    assert "scatter-backfill" in kinds(["galaxy"]) and kinds(["backfill"]) == ["scatter-backfill"]
    assert kinds(["scatter", "--phenomena"]) == [
        "scatter-phenomena", "finish-phenomena", "draw-nebula-field", "refresh-block-lights"]
    assert kinds(["fill-map-luminosity"]) == ["fill-map-luminosity", "refresh-block-lights"]
    assert kinds(["rescore-habitability"]) == ["rescore-habitability"]
    assert kinds(["clean-nebula-field"]) == ["clean-nebula-field"]
    assert kinds(["clear-scatter", "--luminosity"]) == ["clear-scatter"]


def test_the_habitability_re_score_makes_the_block_totals_in_its_conclude_hook():
    """PERF.202: the tail the starter used to do after the last range is the command's `conclude`, run by the unit that
    completes the stage."""
    assert commands.conclude_function("rescore-habitability").__name__ == "conclude_rescore_habitability"


def test_the_star_scatters_are_the_stages_the_wrapper_queues_first_in_the_order_of_their_passes():
    args = _args(["galaxy"])
    assert [step.kind for step in wrapper.queued(wrapper.steps_for(args.command, args))
            if step.kind in ("scatter-mass", "scatter-luminosity")] == ["scatter-mass", "scatter-luminosity"]
    for argv, kinds in ((["scatter", "--mass"], ["scatter-mass", "refresh-block-lights"]),
                        (["scatter", "--luminosity"], ["scatter-luminosity", "refresh-block-lights"]),
                        (["scatter", "--phenomena"], ["scatter-phenomena", "finish-phenomena", "draw-nebula-field",
                                                      "refresh-block-lights"])):
        args = _args(argv)
        assert [step.kind for step in wrapper.queued(wrapper.steps_for(args.command, args))] == kinds
    args = _args(["scatter", "--band", "--down-to", "10"])
    assert [step.kind for step in wrapper.queued(wrapper.steps_for(args.command, args))] == ["scatter-band"]


def test_only_a_registered_command_with_units_can_be_queued_and_every_command_has_a_summary(monkeypatch):
    for kind, command in commands.REGISTRY.items():
        assert command.kind == kind and command.summary
        assert (command.unit is None) == (command.expand is None), "a command has both its unit and its cut, or neither"
    with pytest.raises(ValueError):
        commands.build_steps("phenomenon", _args(["galaxy"]))
    monkeypatch.setitem(commands.REGISTRY, "an-inline-command", commands.CommandSpec("an-inline-command", "Runs here"))
    with pytest.raises(ValueError, match="cannot be cut"):
        commands.expand_function("an-inline-command")


def test_the_step_the_wrapper_adds_when_the_stored_scatter_is_incomplete_comes_before_the_luminosity_pass():
    args = _args(["scatter", "--luminosity"])
    step, listing = wrapper.step_for(args, "mass", "Massive stars", {"seed": 5})
    assert step.params["seed"] == 5 and step.kind == "scatter-mass"
    keys = [item.key for item in listing]
    assert keys.index("mass") < keys.index("luminosity") and keys.count("mass") == 1


def test_arguments_from_neither_command_give_a_list_of_the_one_step():
    class Args:
        workers = 2

    step, listing = wrapper.step_for(Args(), "mass", "Massive stars", {"seed": 1})
    assert listing == [step] and step.params["workers"] == 2 and step.params["seed"] == 1


# --- The wrapper runs a list of stages: the done set queues the next one ---------------------------------------

@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


def _register(monkeypatch, **cuts):
    """Test commands `kind=cut(params) -> Expansion`, whose units record where they ran in `RAN`."""
    for kind in cuts:
        monkeypatch.setitem(commands.REGISTRY, kind, commands.CommandSpec(kind, "A test", unit=(__name__, "_unit"),
                                                                          expand=(__name__, "_cut")))
    monkeypatch.setattr(commands, "unit_function", lambda kind: lambda args, report: RAN.append((kind, args["n"])) or {
        "kind": kind, "n": args["n"]})
    monkeypatch.setattr(commands, "expand_function", lambda kind: cuts[kind])
    RAN.clear()


RAN = []


def _units(count):
    return lambda params: commands.Expansion([({"n": n}, 1.0) for n in range(count)], {"count": count})


def test_the_unit_that_completes_a_stage_queues_the_next_in_the_order_of_the_list(monkeypatch, conn):
    _register(monkeypatch, first=_units(3), second=_units(2), third=_units(1))
    monkeypatch.setitem(commands.REGISTRY, "an-inline-command", commands.CommandSpec("an-inline-command", "Runs here"))
    plans, results = [], {}
    listing = [commands.Step("inline", "an-inline-command", "An inline step"),
               commands.Step("a", "first", "First"), commands.Step("b", "second", "Second"),
               commands.Step("c", "third", "Third")]
    run_id = wrapper.run("Three stages", listing[1:], lambda index, plan: plans.append((index, plan["kind"], plan["units"])),
                         lambda index, unit, result: results.__setitem__((index, unit), result), listing=listing)
    assert plans == [(0, "first", 3), (1, "second", 2), (2, "third", 1)]
    assert [kind for kind, _n in RAN[:3]] == ["first"] * 3, "a stage is queued only when the one before it is done"
    assert sorted(RAN) == sorted([("first", n) for n in range(3)] + [("second", n) for n in range(2)] + [("third", 0)])
    assert [kind for kind, _n in RAN] == ["first"] * 3 + ["second"] * 2 + ["third"]
    assert results[(1, 1)] == {"kind": "second", "n": 1} and len(results) == 6
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "done" and [stage["key"] for stage in run["stages"]] == ["a", "b", "c"]
    assert [stage["units"] for stage in run["stages"]] == [3, 2, 1] == [stage["done"] for stage in run["stages"]]
    assert [step["key"] for step in runrecord.read_steps(conn, run_id)] == ["inline", "a", "b", "c"], \
        "the run record keeps the whole list, the inline step too"
    assert runrecord.active_run(conn) is None


def test_a_stage_with_no_units_is_passed_over_and_the_run_still_ends(monkeypatch, conn):
    _register(monkeypatch, first=_units(2), empty=_units(0), last=_units(1), none=_units(0))
    seen = []
    steps = [commands.Step("a", "first", "A"), commands.Step("b", "empty", "B"), commands.Step("c", "last", "C"),
             commands.Step("d", "none", "D")]
    run_id = wrapper.run("Gaps", steps, lambda index, plan: seen.append((index, plan["units"])),
                         lambda index, unit, result: None)
    assert seen == [(0, 2), (1, 0), (2, 1), (3, 0)]
    assert [kind for kind, _n in RAN] == ["first", "first", "last"]
    assert runrecord.read_run(conn, run_id)["state"] == "done"


def test_the_plan_of_a_stage_is_in_redis_for_any_process_following_the_run(monkeypatch, conn):
    _register(monkeypatch, first=_units(2))
    run_id = wrapper.run("One", [commands.Step("a", "first", "A", {"x": 1})], lambda index, plan: None,
                         lambda index, unit, result: None)
    plan = runrecord.read_plan(conn, run_id, 0)
    assert plan == {"kind": "first", "units": 2, "cost": 2.0, "detail": {"count": 2}}
    assert runrecord.read_plan(conn, run_id, 1) is None
    assert runrecord.read_steps(conn, run_id)[0]["params"] == {"x": 1}


def test_a_unit_that_fails_in_a_later_stage_fails_the_run(monkeypatch, conn):
    _register(monkeypatch, first=_units(1), second=_units(2))

    def explode(args, report):
        if args["n"] == 1:
            raise ValueError("no stars here")
        return {}

    monkeypatch.setattr(commands, "unit_function", lambda kind: (lambda args, report: {}) if kind == "first" else explode)
    with pytest.raises(units.RunStopped, match="unit 1 of stage B failed: ValueError: no stars here"):
        wrapper.run("Fails", [commands.Step("a", "first", "A"), commands.Step("b", "second", "B")],
                    lambda index, plan: None, lambda index, unit, result: None)
    assert runrecord.read_run(conn, runrecord.recent_runs(conn)[0])["state"] == "failed"


def test_a_full_redis_stops_the_run_with_one_plain_error(monkeypatch, conn):
    import redis
    _register(monkeypatch, first=_units(1))

    def full(*args, **kwargs):
        raise redis.exceptions.OutOfMemoryError("command not allowed when used memory > 'maxmemory'.")

    monkeypatch.setattr(units, "follow_stage", full)
    with pytest.raises(RuntimeError, match="Redis is full.*Raise maxmemory"):
        wrapper.run("Full", [commands.Step("a", "first", "A")], lambda index, plan: None, lambda index, unit, result: None)
    assert runrecord.read_run(conn, runrecord.recent_runs(conn)[0])["state"] == "failed"


# --- The commands run alone ------------------------------------------------------------------------------------

def _alone(args, report):
    report(1, 1)
    return {"squared": args["n"] ** 2}


def test_a_unit_runs_alone_and_the_command_line_lists_and_runs_commands(monkeypatch, tmp_path):
    from planetgen.util import log
    monkeypatch.setitem(commands.REGISTRY, "alone", commands.CommandSpec("alone", "Runs alone", unit=(__name__, "_alone")))
    monkeypatch.setitem(commands.REGISTRY, "an-inline-command", commands.CommandSpec("an-inline-command", "Runs here"))
    lines = []
    monkeypatch.setattr(log, "normal", lines.append)
    assert commands.run_alone("alone", {"n": 4}) == {"squared": 16}
    generate_cli.run_commands(_args(["commands"]))
    assert any(line.startswith("scatter-mass") and "units" in line for line in lines)
    assert any(line.startswith("an-inline-command") and "inline" in line for line in lines)
    assert len(lines) == len(commands.REGISTRY), "one line per command"
    path = tmp_path / "unit.json"
    path.write_text(json.dumps({"n": 5}))
    lines.clear()
    generate_cli.run_command_unit(_args(["command-run", "alone", f"@{path}"]))
    generate_cli.run_command_unit(_args(["command-run", "alone", '{"n": 3}']))
    assert [json.loads(line) for line in lines] == [{"squared": 25}, {"squared": 9}]
