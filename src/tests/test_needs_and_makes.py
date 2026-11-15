"""OPS.238, first part: what each command needs and makes (`queue/needs.py`) is checked against the commands, against each other and against the order a galaxy run has today."""

import pytest

from planetgen.generation import stages
from planetgen.queue import commands, needs


def test_every_command_named_is_registered_and_every_artifact_is_in_the_vocabulary():
    unknown = sorted({kind for kind in (*needs.MAKES, *needs.NEEDS) if kind not in commands.REGISTRY})
    assert not unknown, f"declared for commands that are not registered: {unknown}"
    outside = sorted({name for table in (needs.MAKES, needs.NEEDS) for names in table.values() for name in names if name not in needs.ARTIFACTS})
    assert not outside, f"artifacts that are not in ARTIFACTS: {outside}"


def test_every_artifact_is_made_by_a_command_and_needed_by_one_or_made_as_an_end_product():
    unmade = sorted(name for name in needs.ARTIFACTS if not needs.makers_of(name))
    assert not unmade, f"artifacts no command makes: {unmade}"
    needed_unmade = sorted({name for names in needs.NEEDS.values() for name in names if not needs.makers_of(name)})
    assert not needed_unmade, f"needed and never made: {needed_unmade}"


def test_a_command_never_needs_what_it_makes():
    own = sorted(kind for kind in needs.NEEDS if set(needs.NEEDS[kind]) & set(needs.makes_of(kind)))
    assert not own, own


def test_the_graph_has_no_cycle():
    graph = {}
    for needer, maker in needs.edges():
        if needer != maker:
            graph.setdefault(needer, set()).add(maker)
    done, active = set(), []

    def visit(kind):
        assert kind not in active, f"a cycle: {active[active.index(kind):] + [kind]}"
        if kind in done:
            return
        active.append(kind)
        for maker in graph.get(kind, ()):
            visit(maker)
        active.pop()
        done.add(kind)

    for kind in graph:
        visit(kind)


def _galaxy_kinds(**flags):
    """The command of each stage of a galaxy run, in its order (`stages.galaxy_stages` for a random start, the unplaced and the system kinds left out)."""
    args = type("Args", (), flags)()
    return [commands.KIND_OF_STAGE[stage.key] for stage in stages.galaxy_stages(args)]


def test_the_order_of_a_galaxy_run_meets_every_need_a_command_of_the_run_can_meet():
    kinds = _galaxy_kinds()
    assert kinds[0] == "generate-start-sector" and kinds[-1] == "refresh-block-lights" and "run-territories" in kinds
    for place, kind in enumerate(kinds):
        for artifact in needs.needs_of(kind):
            makers = [maker for maker in needs.makers_of(artifact) if maker in kinds]
            assert not makers or any(kinds.index(maker) < place for maker in makers), f"{kind} (stage {place + 1}) needs {artifact}, made only after it by {makers}"


@pytest.mark.parametrize("scatter", ["mass", "luminosity", "phenomena", "nebula_field", "band"])
def test_the_order_of_a_scatter_run_meets_the_needs_too(scatter):
    args = type("Args", (), {scatter: True, "prepare_limits": scatter == "mass"})()
    kinds = [commands.KIND_OF_STAGE[stage.key] for stage in stages.scatter_stages(args)]
    assert kinds
    for place, kind in enumerate(kinds):
        for artifact in needs.needs_of(kind):
            makers = [maker for maker in needs.makers_of(artifact) if maker in kinds]
            assert not makers or any(kinds.index(maker) < place for maker in makers), f"{kind} needs {artifact}, made only after it by {makers}"


def test_what_a_command_requires_today_is_a_need_the_declarations_have():
    """`commands.REQUIRES` goes when the planner reads the declarations: until then every requirement it holds is one they state."""
    for kind, required in commands.REQUIRES.items():
        for other in required:
            assert set(needs.needs_of(kind)) & set(needs.makes_of(other)), f"{kind} requires {other}, which no declared need of {kind} is made by"
