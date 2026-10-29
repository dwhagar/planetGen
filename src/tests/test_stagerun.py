"""PERF.202: the stages of a generating verb as parts of one wrapper run (`planetgen/generation/stagerun.py`)."""

import datetime

import pytest

from planetgen.generation import run_galaxy, stagerun
from planetgen.queue import commands, redisqueue, runrecord, wrapper


def _count(args, report):
    report(1, 1)
    return {"counted": args["n"]}


def _cut_two(params):
    return commands.Expansion([({"n": 1}, 1.0), ({"n": 2}, 1.0)], {"sectors": 2})


def _cut_refused(params):
    raise RuntimeError(f"The {params['what']} is already stored -- run 'planetgen clear-scatter' first.")


def _unit_refused(args, report):
    raise RuntimeError("The bright stars are split from the massive stars at the mass limit.")


def _cut_one(params):
    return commands.Expansion([({}, 1.0)], {"sectors": 1})


@pytest.fixture
def args(mysql_config):
    return run_galaxy._default_generation_args(mysql_config)


def _register(monkeypatch, **named):
    """Registers a test command for each `kind=(expand, unit)`."""
    for kind, (expand, unit) in named.items():
        monkeypatch.setitem(commands.REGISTRY, kind, commands.CommandSpec(
            kind, "A test", unit=(__name__, unit.__name__), expand=(__name__, expand.__name__), whole=True))
    monkeypatch.setattr(commands, "unit_function", lambda kind: named[kind][1])
    monkeypatch.setattr(commands, "expand_function", lambda kind: named[kind][0])


def _part(mysql_config, kind, **kwargs):
    step = commands.Step(kind, kind, f"Stage {kind}", {"database": wrapper.database_args(mysql_config), **kwargs.pop("params", {})})
    return stagerun.Part(step, **kwargs)


def test_two_stages_are_followed_in_one_run_each_with_its_own_totals(monkeypatch, args, mysql_config):
    _register(monkeypatch, first=(_cut_two, _count), second=(_cut_two, _count))
    first = _part(mysql_config, "first", sums=("counted",), bar="First", advance=lambda result: result["counted"])
    second = _part(mysql_config, "second", sums=("counted",))
    run_id = stagerun.run_parts(args, "Two stages", [first, second])
    assert first.totals == {"counted": 3} and second.totals == {"counted": 3}
    conn = redisqueue.connect()
    record = runrecord.read_run(conn, run_id)
    assert record["state"] == "done" and [stage["key"] for stage in record["stages"]] == ["first", "second"]
    assert [step["kind"] for step in runrecord.read_steps(conn, run_id)] == ["first", "second"]


def test_a_stage_that_refuses_stops_the_run_with_its_own_words_before_the_next_one(monkeypatch, args, mysql_config):
    _register(monkeypatch, refused=(_cut_refused, _count), after=(_cut_two, _count))
    refusing = _part(mysql_config, "refused", params={"what": "massive-star scatter"})
    after = _part(mysql_config, "after", sums=("counted",))
    with pytest.raises(RuntimeError, match="The massive-star scatter is already stored -- run 'planetgen clear-scatter'"):
        stagerun.run_parts(args, "Refused first", [refusing, after])
    assert after.totals == {"counted": 0}, "the stage after a refusal never ran"


def test_two_refusals_fire_in_the_order_the_stages_are_listed(monkeypatch, args, mysql_config):
    _register(monkeypatch, expansion=(_cut_refused, _count), unit=(_cut_one, _unit_refused), fine=(_cut_two, _count))
    first_in_order = [_part(mysql_config, "unit"), _part(mysql_config, "expansion", params={"what": "bright-star scatter"})]
    with pytest.raises(RuntimeError, match="split from the massive stars"):
        stagerun.run_parts(args, "Unit refuses first", first_in_order)
    other_order = [_part(mysql_config, "expansion", params={"what": "bright-star scatter"}), _part(mysql_config, "unit")]
    with pytest.raises(RuntimeError, match="The bright-star scatter is already stored"):
        stagerun.run_parts(args, "Expansion refuses first", other_order)
    ok_then_refusal = [_part(mysql_config, "fine", sums=("counted",)), _part(mysql_config, "unit")]
    with pytest.raises(RuntimeError, match="split from the massive stars"):
        stagerun.run_parts(args, "Both stages", ok_then_refusal)
    assert ok_then_refusal[0].totals == {"counted": 3}, "the stage before the refusal was done and counted"


def test_the_link_settle_and_backfill_parts_are_one_run_of_three_stages(args, mysql_config):
    since = datetime.datetime(2000, 1, 1)
    parts = [run_galaxy.link_part(args, since), run_galaxy.backfill_part(args, since)[0],
             run_galaxy.settle_part(args, since)]
    run_id = stagerun.run_parts(args, "Three stages", parts)
    conn = redisqueue.connect()
    assert [step["kind"] for step in runrecord.read_steps(conn, run_id)][-3:] == [
        "link-neighbours", "scatter-backfill", "settle-sector-paths"]
    assert parts[0].totals == {"linked": 0} and parts[2].totals == {"saved": 0}
