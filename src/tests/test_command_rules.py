"""PERF.120: the rules every command follows (docs/design/architecture/command-inventory.md, "Rules for a new
command"), held for every kind in the registry so a new command cannot skip one."""

import inspect
import pathlib

import pytest

from planetgen.queue import commands

DOCS = pathlib.Path(__file__).resolve().parents[2] / "docs"


@pytest.fixture(scope="module")
def documentation():
    return "\n".join(path.read_text(encoding="utf-8") for path in (DOCS / "design").rglob("*.md")) \
        + "\n".join(path.read_text(encoding="utf-8") for path in (DOCS / "reference").rglob("*.md"))


@pytest.mark.parametrize("kind", sorted(commands.REGISTRY))
def test_a_command_is_a_registered_kind_with_json_arguments_and_a_result(kind):
    """One path: a unit that takes `(params, report)` and an expansion that takes `(params)`, both found by name."""
    spec = commands.spec(kind)
    assert spec.kind == kind and spec.summary.strip(), "a command has a name and says what it does"
    assert not spec.inline, f"{kind} is inline: give it a unit and an expand (PERF.120 leaves no second way)"
    unit, expand = commands.unit_function(kind), commands.expand_function(kind)
    assert len(inspect.signature(unit).parameters) == 2, f"{kind}'s unit takes (params, report)"
    assert len(inspect.signature(expand).parameters) == 1, f"{kind}'s expand takes (params)"


@pytest.mark.parametrize("kind", sorted(commands.REGISTRY))
def test_a_command_declares_how_long_its_longest_unit_runs_or_that_it_is_one_unit(kind):
    spec = commands.spec(kind)
    assert spec.whole or spec.unit_seconds, f"{kind}: set whole=True, or unit_seconds (PERF.155)"


@pytest.mark.parametrize("kind", sorted(commands.REGISTRY))
def test_a_command_has_a_docs_line(kind, documentation):
    assert kind in documentation, f"{kind} is in no document under docs/design or docs/reference"


def test_the_command_line_names_only_registered_commands():
    for table in (commands.KIND_OF_STAGE, commands.KIND_OF_WARM_FAMILY, commands.KIND_OF_CHECK):
        assert set(table.values()) <= set(commands.REGISTRY)
