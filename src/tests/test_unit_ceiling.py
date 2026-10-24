"""PERF.155: no unit of a stage is expected to run longer than `queue.unit_ceiling_seconds` unless the stage is one unit by
nature. A test holds every registered command to it, and a unit that runs longer is written to the run's log."""

import time

import pytest

from planetgen.queue import commands, redisqueue, runrecord, wrapper
from planetgen.util import settings


def test_the_ceiling_is_a_setting_that_defaults_to_five_minutes():
    option = settings.Settings.model_fields["queue"].annotation.model_fields["unit_ceiling_seconds"]
    assert option.default == 300 and settings.get_settings().queue.unit_ceiling_seconds > 0


def test_every_queued_command_is_one_unit_by_nature_or_declares_a_longest_unit_under_the_ceiling():
    ceiling = settings.Settings.model_fields["queue"].annotation.model_fields["unit_ceiling_seconds"].default
    for kind, command in commands.REGISTRY.items():
        if command.inline:
            assert not command.whole and command.unit_seconds is None, f"{kind} is inline: it has no units to bound"
        elif command.whole:
            assert command.unit_seconds is None, f"{kind} is one unit by nature: a longest unit has no meaning"
        else:
            assert command.unit_seconds is not None, f"{kind} is cut into units and must declare how long its longest runs"
            assert 0 < command.unit_seconds <= ceiling, \
                f"{kind}: its longest unit is expected at {command.unit_seconds} s, over the {ceiling} s ceiling; cut it finer"


def _unit(args, report):
    time.sleep(args["pause"])
    return {}


def _cut(params):
    return commands.Expansion([({"pause": params["pause"]}, 1.0)] * 2, {})


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


def _long_lines(conn, run_id):
    return [fields for _id, fields in runrecord.read_events(conn, run_id) if fields.get("kind") == "unit.long"]


def _run(monkeypatch, conn, whole, pause, ceiling):
    monkeypatch.setitem(commands.REGISTRY, "ceiling-test", commands.CommandSpec(
        "ceiling-test", "A test", unit=(__name__, "_unit"), expand=(__name__, "_cut"), whole=whole,
        unit_seconds=None if whole else 1))
    monkeypatch.setattr(commands, "unit_function", lambda kind: _unit)
    monkeypatch.setattr(commands, "expand_function", lambda kind: _cut)
    monkeypatch.setattr(settings.get_settings().queue, "unit_ceiling_seconds", ceiling)
    steps = [commands.Step("only", "ceiling-test", "Only stage", {"pause": pause})]
    run_id = wrapper.run("Ceiling", steps, lambda index, plan: None, lambda index, unit, result: None, listing=steps)
    return _long_lines(conn, run_id)


def test_a_unit_that_runs_over_the_ceiling_is_written_to_the_runs_log_naming_the_stage(monkeypatch, conn):
    lines = _run(monkeypatch, conn, whole=False, pause=0.3, ceiling=0.1)
    assert len(lines) == 1 and lines[0]["level"] == "warning"
    assert "Only stage" in lines[0]["message"] and "unit ceiling" in lines[0]["message"] and "2 units" in lines[0]["message"]


def test_a_unit_under_the_ceiling_writes_nothing(monkeypatch, conn):
    assert _run(monkeypatch, conn, whole=False, pause=0.01, ceiling=30) == []


def test_a_stage_that_is_one_unit_by_nature_is_not_held_to_the_ceiling(monkeypatch, conn):
    assert _run(monkeypatch, conn, whole=True, pause=0.3, ceiling=0.1) == []
