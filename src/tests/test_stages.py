# tests/test_stages.py

"""
UX.89: every staged action numbers and announces its stages, the optional ones
as skipped with a reason, and the Generate page counts the stages of a whole
job instead of its processes.
"""

import argparse
import json

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.generation import stages
from planetgen.queue import progress_file
from planetgen.util import log
from planetgen.web import generate_page, jobs


def _galaxy(*flags):
    _parser, parsers = generate_cli.build_parser()
    return parsers["galaxy"].parse_args(list(flags))


def _plan(*flags):
    _parser, parsers = generate_cli.build_parser()
    return parsers["plan"].parse_args(list(flags))


def _labels(found):
    return [(stage.key, stage.skip is not None) for stage in found]


def test_a_random_start_galaxy_run_has_nine_stages_with_the_unasked_ones_skipped():
    found = stages.galaxy_stages(_galaxy("--then-scatter"))
    assert _labels(found) == [
        ("start", False), ("neighborhood", False), ("link", False), ("mass", False), ("luminosity", False),
        ("phenomena", False), ("backfill", False), ("population", True), ("settle", False)]
    assert "population" in found[7].skip


def test_without_then_scatter_the_three_scatter_stages_are_skipped_and_still_numbered():
    found = stages.galaxy_stages(_galaxy("--backfill-from", "none", "--no-settle"))
    assert len(found) == 9
    assert {stage.key for stage in found if stage.skip} == {"mass", "luminosity", "phenomena", "backfill",
                                                            "population", "settle"}


def test_a_galaxy_run_with_an_address_has_one_sector_stage():
    found = stages.galaxy_stages(_galaxy("--ring", "3", "--layer", "0"))
    assert [stage.key for stage in found][:2] == ["sectors", "link"] and len(found) == 8


def test_plan_stages_follow_the_flags():
    assert [s.key for s in stages.plan_stages(_plan("--no-bright-stars"))] == ["skeleton"]
    assert [s.key for s in stages.plan_stages(_plan())] == ["skeleton", "phenomena", "mass", "luminosity"]
    assert [s.key for s in stages.plan_stages(_plan("--bright-stars-only"))] == ["mass", "luminosity"]
    assert [s.key for s in stages.plan_stages(_plan("--phenomena-only"))] == ["phenomena"]


def test_entering_a_stage_announces_the_skipped_ones_before_it_and_the_rest_at_the_end(monkeypatch, tmp_path):
    lines = []
    monkeypatch.setattr(log, "normal", lambda message, *a, **k: lines.append(message))
    monkeypatch.setenv(progress_file.ENV_VAR, str(tmp_path / "progress.json"))
    stages.begin(stages.galaxy_stages(_galaxy("--then-scatter")))
    stages.enter("start")
    stages.enter("link")           # "neighborhood" never ran
    stages.skip("mass", "no layers")
    stages.enter("settle")         # the ones between are skipped
    stages.finish()
    text = "\n".join(lines)
    assert text.startswith("9 stages: 1. Generate the starting sector;")
    assert "Stage 2 of 9: Generate the neighborhood -- skipped: not needed." in text
    assert "Stage 4 of 9: Scatter the massive stars -- skipped: no layers." in text
    assert "Stage 8 of 9: Run the population pass -- skipped: the population pass was not asked for." in text
    assert "Stage 9 of 9: Save the sector paths." in text
    assert json.loads((tmp_path / "progress.json").read_text())["stage"]["index"] == 9


def test_a_stage_key_the_run_does_not_have_is_ignored():
    stages.begin(stages.plan_stages(_plan("--phenomena-only")))
    stages.enter("mass")
    stages.finish()


FORM = {"confirm": "db", "max_ring": "", "skip_bright_stars": ""}


def _steps(action, form):
    return generate_page.build_job(action, form, "db")[2]


def test_the_new_galaxy_job_counts_every_stage_of_every_step():
    steps = _steps("new_galaxy", {"confirm": "db"})
    found = jobs.job_stages(steps)
    assert [s["label"] for s in found][:3] == [generate_page.MATH_CHECK_LABEL, "Reset the galaxy", "Plan the galaxy"]
    assert len(found) == 12 and found[-1]["n"] == 12
    assert [s["skipped"] is not None for s in found].count(True) == 1     # the population pass


def test_ticking_skip_the_scatter_lists_those_stages_as_skipped_with_that_reason():
    found = jobs.job_stages(_steps("new_galaxy", {"confirm": "db", "skip_bright_stars": "on"}))
    scatter = [s for s in found if s["label"].startswith("Scatter")]
    assert len(scatter) == 3 and all("box is ticked" in s["skipped"] for s in scatter)


def test_a_plan_job_lists_its_scatter_steps_stages():
    found = jobs.job_stages(_steps("plan", {"confirm": "db"}))
    assert [s["label"] for s in found][-3:] == ["Scatter the massive stars", "Scatter the bright stars",
                                                  "Scatter the phenomena"]


def test_the_stage_view_names_the_running_stage_by_its_place_in_the_whole_job():
    steps = _steps("new_galaxy", {"confirm": "db"})
    job = {"stages": jobs.job_stages(steps), "step": 4, "finished": False,
           "progress": {"stage": {"index": 4, "total": 9, "label": "Scatter the massive stars", "skipped": None}}}
    view = generate_page.stage_view(job)
    assert view["stage_text"] == "Stage 7 of 12: Scatter the massive stars"
    states = [s["state"] for s in view["stage_list"]]
    assert states[:6] == ["done"] * 6 and states[6] == "running" and "waiting" in states
    assert [s["state"] for s in view["stage_list"] if s["state"] == "skipped"] == ["skipped"]
