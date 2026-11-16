"""OPS.267: the plan's closing line and warning, which read `args.edge` only where the Edge dial was given."""

import types

from planetgen.generation import run_plan

SUMMARY = {"shell_count": 4, "outer_shell_index": 3, "reach_shell_index": 3, "total_candidate_sectors": 1234, "elapsed_s": 0.5, "edge_confirmed": False}


def test_a_plan_cut_by_its_reach_with_no_edge_says_so_and_does_not_fail(caplog):
    """`planetgen plan --reach-shells 3` has no `--edge` (`args.edge` is `None`): the warning names the reach."""
    run_plan.log_skeleton_summary(types.SimpleNamespace(edge=None, reach_shells=3), SUMMARY)
    assert "--reach-shells" in caplog.text and "shell 3" in caplog.text
    assert "--edge" not in caplog.text


def test_a_plan_cut_by_its_edge_still_names_the_edge(caplog):
    run_plan.log_skeleton_summary(types.SimpleNamespace(edge=20000.0, reach_shells=None), SUMMARY)
    assert "at the --edge (" in caplog.text and "Re-run with a larger --edge" in caplog.text


def test_a_plan_that_reached_its_edge_warns_of_nothing(caplog):
    run_plan.log_skeleton_summary(types.SimpleNamespace(edge=None, reach_shells=3), {**SUMMARY, "edge_confirmed": True})
    assert "true edge was not reached" not in caplog.text
