# tests/test_cut_neutral.py

"""
PERF.260: a stage may be flagged "dynamic split is byte identical" (`commands.SPLIT_BYTE_IDENTICAL`, so the queue engine may pick any cut for it) only if it stores the same bytes at two cuts.

* `test_every_command_flagged_byte_identical_stores_the_same_bytes_at_two_cuts` is the merge gate's test (Boss, 2026-11-12, 3:09 PM EST: "if a flag is ever set when something is not byte
  identical then it should fail"). It finds the flagged commands in `commands.REGISTRY` itself, with no list kept by hand, runs each at two cut sizes from the same stored state and fails on
  any difference. It fails, never warns or skips: a flagged command with no recipe in `tests/cut_neutral_recipes.py` cannot be run at two cuts and fails until its flag is removed.
* The probe tests show the harness catches what it must: a stage that is neutral passes, one that stores a different thing at a different cut fails, as do one that does not reproduce
  itself, one whose two cuts make the same number of units and one that stores nothing.
* `test_each_recipe_cuts_its_stage_two_ways` runs every recipe and prints its verdict (the report of PERF.260, `-s`); it fails only if the recipe itself proves nothing (the two cuts make the
  same number of units, or the stage stores nothing), not when the stage turns out to be cut dependent: a stage that is not flagged may be.
"""

import dataclasses

import pytest

from planetgen.db import store
from planetgen.queue import commands
from tests import cut_neutral_recipes as recipes
from tests.cut_neutral_support import (Context, Recipe, check, dump, flagged_commands,
                                       gate_failures, restore)

pytestmark = pytest.mark.db

BYTE_IDENTICAL = getattr(commands, "SPLIT_BYTE_IDENTICAL", "byte-identical")


@pytest.fixture
def context(mysql_config, monkeypatch, redis_server, large_spiral_fallback):
    return Context(mysql_config, monkeypatch)


# --- the gate ------------------------------------------------------------------------------------

def test_every_command_flagged_byte_identical_stores_the_same_bytes_at_two_cuts(context):
    pristine = dump(context.mysql_config)

    def fresh(_recipe):
        restore(context.mysql_config, pristine)      # each recipe starts from the empty database, not from what the last one left
        return context

    failures = gate_failures(recipes.RECIPES, fresh)
    assert not failures, "a flag \"dynamic split is byte identical\" is set on a command that is not:\n" + "\n".join(failures)


# --- every recipe, as a report ---------------------------------------------------------------------

@pytest.mark.parametrize("kind", sorted(recipes.RECIPES))
def test_each_recipe_cuts_its_stage_two_ways(context, kind):
    verdict = check(recipes.RECIPES[kind], context)
    print(verdict.line(), "| stored:", verdict.changed_tables, "|", [line[:90] for line in verdict.differences])
    assert verdict.distinct_cuts and verdict.changed_tables, "the recipe proves nothing: " + verdict.line()


# --- the probe: a small stage of its own, to show the harness catches what it must ---------------------------

PROBE_CUT = 1
PROBE_MODE = "neutral"
PROBE_ITEMS = 6
_calls = []


def probe_expand(params):
    pieces = [list(range(start, min(start + PROBE_CUT, PROBE_ITEMS))) for start in range(0, PROBE_ITEMS, PROBE_CUT)]
    return commands.Expansion([({"database": params["database"], "items": piece, "piece": len(piece)}, float(len(piece))) for piece in pieces], {})


def probe_unit(params, report):
    conn = store.get_connection(_config_of(params))
    try:
        for item in params["items"]:
            value = {"neutral": float(item) * 1.5,
                     "cut dependent": float(item) * 1.5 + params["piece"],       # the size of the piece leaks into the bytes
                     "unstable": float(item) * 1.5 + len(_calls),                # a different number every run
                     "nothing": None}[PROBE_MODE]
            _calls.append(item)
            if value is not None:
                conn.execute("INSERT INTO cut_probe (item, value) VALUES (?, ?)", (item, value))
        conn.commit()
        return {"items": len(params["items"])}
    finally:
        conn.close()


_config = {}


def _config_of(params):
    return _config["config"]


@dataclasses.dataclass(frozen=True)
class FlaggedSpec(commands.CommandSpec):
    """A command spec that carries the flag as `OPS.161` will (the field is the registry's own once that lands)."""

    split: str = BYTE_IDENTICAL


def _probe_recipe(cuts=((1,), (6,))):
    def prepare(context):
        conn = store.get_connection(context.mysql_config)
        try:
            conn.execute("CREATE TABLE IF NOT EXISTS cut_probe (item INT, value DOUBLE NULL)")
            conn.execute("TRUNCATE TABLE cut_probe")
            conn.commit()
        finally:
            conn.close()

    return Recipe("cut-probe", (("tests.test_cut_neutral", "PROBE_CUT"),), cuts, prepare=prepare,
                  run_params=lambda context: {"database": "x"})


def _register(context, monkeypatch, mode, flagged=True, recipe=True):
    _config["config"] = context.mysql_config
    monkeypatch.setattr("tests.test_cut_neutral.PROBE_MODE", mode)
    spec_class = FlaggedSpec if flagged else commands.CommandSpec
    monkeypatch.setitem(commands.REGISTRY, "cut-probe", spec_class("cut-probe", "A probe stage", unit=("tests.test_cut_neutral", "probe_unit"),
                                                                    expand=("tests.test_cut_neutral", "probe_expand")))
    return {"cut-probe": _probe_recipe()} if recipe else {}


def _probe_only():
    """The registry of the probe alone: the commands the program flags are held to their own recipes by the gate's test, and these tests are about the probe."""
    return {"cut-probe": commands.REGISTRY["cut-probe"]}


def test_a_flagged_stage_that_stores_the_same_bytes_at_any_cut_passes(context, monkeypatch):
    found = _register(context, monkeypatch, "neutral")
    assert gate_failures(found, lambda _recipe: context, _probe_only()) == []
    verdict = check(found["cut-probe"], context)
    assert verdict.passed and verdict.units == (6, 1) and verdict.line().startswith("IDENTICAL cut-probe")


def test_a_flagged_stage_that_stores_a_different_thing_at_a_different_cut_fails_the_gate(context, monkeypatch):
    found = _register(context, monkeypatch, "cut dependent")
    failures = gate_failures(found, lambda _recipe: context, _probe_only())
    assert len(failures) == 1 and failures[0].startswith("FAILS cut-probe") and "the cuts differ in the stored bytes" in failures[0]
    assert "cut_probe" in failures[0], "the message names the table that differs"


def test_a_flagged_stage_that_does_not_reproduce_itself_fails_the_gate(context, monkeypatch):
    found = _register(context, monkeypatch, "unstable")
    failures = gate_failures(found, lambda _recipe: context, _probe_only())
    assert len(failures) == 1 and "the same cut does not reproduce itself" in failures[0]


def test_a_flagged_stage_that_stores_nothing_fails_the_gate(context, monkeypatch):
    found = _register(context, monkeypatch, "nothing")
    failures = gate_failures(found, lambda _recipe: context, _probe_only())
    assert len(failures) == 1 and "stored nothing" in failures[0]


def test_two_cuts_that_make_the_same_number_of_units_fail_the_gate(context, monkeypatch):
    found = _register(context, monkeypatch, "neutral")
    found["cut-probe"] = _probe_recipe(cuts=((6,), (10,)))      # both are one unit of six
    failures = gate_failures(found, lambda _recipe: context, _probe_only())
    assert len(failures) == 1 and "the two cuts make the same number of units" in failures[0]


def test_a_flagged_command_with_no_recipe_fails_the_gate_until_the_flag_is_removed(context, monkeypatch):
    found = _register(context, monkeypatch, "neutral", recipe=False)
    failures = gate_failures(found, lambda _recipe: context, _probe_only())
    assert len(failures) == 1 and "cut-probe is flagged byte identical but has no recipe" in failures[0]
    monkeypatch.setitem(commands.REGISTRY, "cut-probe", commands.CommandSpec("cut-probe", "A probe stage"))      # the flag is removed
    assert gate_failures(found, lambda _recipe: context, _probe_only()) == []


def test_the_flagged_commands_are_found_in_the_registry_with_no_list_and_an_unflagged_one_is_not_run(context, monkeypatch):
    found = _register(context, monkeypatch, "cut dependent", flagged=False)
    assert "cut-probe" not in flagged_commands()
    assert gate_failures(found, lambda _recipe: context, _probe_only()) == [], "a stage that is not flagged may be cut dependent"
    _register(context, monkeypatch, "cut dependent", flagged=True)
    assert "cut-probe" in flagged_commands()
    assert flagged_commands() == [kind for kind, spec in commands.REGISTRY.items() if getattr(spec, "split", None) == BYTE_IDENTICAL]


def test_every_command_the_registry_flags_has_a_recipe_that_is_not_a_stand_in():
    """The recipes are kept by command kind; one for a kind that is not registered is a stale recipe."""
    assert set(recipes.RECIPES) <= set(commands.REGISTRY)
    assert all(recipe.kind == kind for kind, recipe in recipes.RECIPES.items())
