# tests/test_plan_outline.py

"""
GEN.617: the line that says how dense the planned galaxy is and how much of it the outline holds. The plan prints it,
the Generate page shows it under the dials, and both come from `generation.plan_outline`, so the numbers are the same.
An outline with no sector in it is refused with the reason.
"""

import json

import pytest

from planetgen.db import store
from planetgen.generation import plan_outline, run_plan
from planetgen.queue import commands
from tests.test_galaxy_gen import _PLAN_SHAPE_ARGV
from tests.test_work_units import _skeleton_args, control_database  # noqa: F401 -- the fixture

SEED_HEX = "00112233445566778899AABBCCDDEEFF"
SEED = bytes.fromhex(SEED_HEX)
EMPTY = ["--mass", "3e5", "--disk-scale-length", "500pc"]


def test_the_line_gives_the_crest_the_share_and_the_size_of_the_outline():
    made, refusal = plan_outline.preview(["--mass", "1e8", "--disk-scale-length", "2000pc", "--flatness", "20",
                                          "--bulge-share", "0.001"], SEED)
    assert refusal is None and not made.empty
    assert made.line().startswith("Outline: the disc's crest is 3.32 times the one-star density, about ")
    assert f"{made.shells:,} shells out to " in made.line() and f"{made.sectors:,} candidate sectors" in made.line()
    assert made.shells == len(plan_outline.build(plan_outline._parser().parse_args(
        ["--mass", "1e8", "--disk-scale-length", "2000pc", "--flatness", "20", "--bulge-share", "0.001"]),
        SEED, 4.0).extents)


def test_a_crest_under_twenty_is_warned_of_and_a_dense_disc_is_not():
    thin, _ = plan_outline.preview(["--mass", "1e8", "--disk-scale-length", "2000pc", "--flatness", "20",
                                    "--bulge-share", "0.001"], SEED)
    assert any("times the one-star density (under 20)" in note for note in thin.warnings)
    dense, _ = plan_outline.preview(["--mass", "6e9"], SEED)
    assert dense.crest >= 20 and dense.warnings == []


def test_an_outline_with_no_sector_is_refused_with_the_reason():
    made, refusal = plan_outline.preview(EMPTY, SEED)
    assert made.empty and made.shells == 0 and made.share == 0.0
    assert made.line().startswith("Outline: empty") and "one-star-per-sector density" in refusal


def test_a_flag_that_is_not_valid_is_a_value_error():
    with pytest.raises(ValueError):
        plan_outline.preview(["--mass", "nonsense"], SEED)


@pytest.mark.fresh_database
def test_the_plan_prints_the_line_the_preview_gives(control_database):  # noqa: F811
    args = _skeleton_args(control_database, "--seed", SEED_HEX)
    params = json.loads(json.dumps(run_plan.skeleton_params(args, control_database)))
    result = commands.run_alone("plan-skeleton", params)
    printed = [line for line in result["lines"] if line.startswith("Outline: ")]
    made, _ = plan_outline.preview(_PLAN_SHAPE_ARGV, SEED)
    assert printed == [made.line() + "."]
    assert result["summary"]["shell_count"] == made.shells
    assert result["summary"]["total_candidate_sectors"] == made.sectors


@pytest.mark.fresh_database
def test_the_plan_refuses_an_empty_outline_and_stores_nothing(control_database):  # noqa: F811
    args = _skeleton_args(control_database, "--seed", SEED_HEX, *EMPTY)
    params = json.loads(json.dumps(run_plan.skeleton_params(args, control_database)))
    with pytest.raises(ValueError, match="nothing to plan"):
        commands.run_alone("plan-skeleton", params)
    conn = store.get_connection(control_database)
    try:
        assert store.get_galaxy_shape(conn) is None
    finally:
        conn.close()


@pytest.mark.fresh_database
def test_the_plan_says_where_its_time_went_step_by_step(control_database):  # noqa: F811
    """PERF.163: a table of the unit's steps, the costliest first, whose seconds add up to the total it names."""
    args = _skeleton_args(control_database, "--seed", SEED_HEX)
    params = json.loads(json.dumps(run_plan.skeleton_params(args, control_database)))
    result = commands.run_alone("plan-skeleton", params)
    start = next(i for i, line in enumerate(result["lines"]) if line.startswith("Where the plan's time went ("))
    rows = result["lines"][start + 1:]
    seconds = [float(row.split("s ", 1)[0]) for row in rows]
    assert seconds == sorted(seconds, reverse=True) and len(rows) >= 8
    named = float(result["lines"][start].split("(")[1].split("s)")[0])
    assert sum(seconds) == pytest.approx(named, abs=0.01 * len(rows))
    text = "\n".join(rows)
    for step in ("the globular clusters", "the satellite galaxies", "the component field", "storing the shells",
                 "storing the clusters and satellites", "the outline"):
        assert step in text
