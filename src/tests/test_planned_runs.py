# tests/test_planned_runs.py

"""PERF.203 step 1: a verb's run is a `followed.Planned`, built apart from starting it, and the followers it holds are the stages the verb lists."""

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.generation import followed, run_plan, run_population, stages

PARSERS = generate_cli.build_parser()[1]


def _args(command, *extra):
    namespace = PARSERS[command].parse_args(list(extra))
    namespace.command = command
    return namespace


PLANNED = (
    ("rescore-habitability", run_plan.planned_habitability, ()),
    ("fill-map-luminosity", run_plan.planned_map_luminosity, ()),
    ("clean-nebula-field", run_plan.planned_clean_nebula_field, ()),
    ("clear-scatter", run_plan.planned_clear_scatter, ("--mass",)),
    ("population", run_population.planned_population, ()),
    ("territories", run_population.planned_territories, ()),
)


@pytest.mark.parametrize("command, planned, extra", PLANNED, ids=[row[0] for row in PLANNED])
def test_the_followers_of_a_planned_run_are_the_stages_of_its_verb(command, planned, extra):
    args = _args(command, *extra)
    found = planned(args)
    assert isinstance(found, followed.Planned) and found.title
    assert [follower.key for follower in found.followers] == [stage.key for stage in stages.for_command(command, args)]


def test_run_planned_hands_the_followers_the_title_and_the_options_to_run_followed_and_returns_the_result(monkeypatch):
    calls = []
    monkeypatch.setattr(followed, "run_followed", lambda args, followers, title, **options: calls.append((args, followers, title, options)))
    planned = followed.Planned(["a", "b"], "A run", lambda: 7, seed=5, galaxy=False)
    assert followed.run_planned("args", planned) == 7
    assert calls == [("args", ["a", "b"], "A run", {"seed": 5, "galaxy": False})]


def test_a_planned_run_without_a_result_returns_none(monkeypatch):
    monkeypatch.setattr(followed, "run_followed", lambda *args, **options: None)
    assert followed.run_planned("args", followed.Planned(["a"], "A run")) is None
