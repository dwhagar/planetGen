# tests/test_job_steps.py

"""
PERF.203 step 4: the generation steps of a Generate job run as the one step `planetgen run-job`. The page builds and reviews the job as
its list of commands; when it starts, they are one run with one stage list, which the page numbers and estimates as the commands'
stages were.
"""

import json

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.generation import job, stages
from planetgen.web import generate_page, job_runner, job_steps, jobs

DB = "planetgen_job_steps"

JOBS = [
    ("new_galaxy", {"confirm": DB}),
    ("new_galaxy", {"confirm": DB, "skip_scatters": "1"}),
    ("plan", {}),
    ("redo_scatters", {"redo_mass": "1", "redo_phenomena": "1"}),
    ("galaxy", {"mode": "random"}),
]


def _built(action, form):
    kind, title, steps = generate_page.build_job(action, form, DB)
    return title, steps


def _entries(step):
    _parser, parsers = generate_cli.build_parser()
    start = len(jobs.GENERATE_COMMAND) + 1
    shared = parsers["run-job"].parse_args(step["argv"][start + 1:])
    return job.parse_entries(shared.spec, shared, parsers, generate_cli.validate_command)


@pytest.mark.parametrize("action, form", JOBS)
def test_the_generating_steps_of_a_job_are_one_run_job_step(action, form):
    title, steps = _built(action, form)
    ran = job_steps.collapse(steps, title)
    work = [step for step in ran if step["argv"][1:3] == jobs.GENERATE_COMMAND]
    assert len(work) == 1 and work[0]["argv"][3] == "run-job" and work[0]["label"] == title
    assert [step["label"] for step in ran if step not in work] == [
        step["label"] for step in steps if step["argv"][1:3] != jobs.GENERATE_COMMAND]
    spec = json.loads(work[0]["argv"][4])
    assert [(entry["command"], entry["argv"]) for entry in spec] == [
        (step["argv"][3], [word for word in step["argv"][4:] if word not in job_steps.SHARED_FLAGS])
        for step in steps if step["argv"][1:3] == jobs.GENERATE_COMMAND]


@pytest.mark.parametrize("action, form", JOBS)
def test_the_stage_list_of_the_step_is_the_one_the_run_makes(action, form):
    """The page numbers the job's stages by the list the step carries and the run announces its stages by its own: they are the same,
    stage for stage, with the block light totals once at the last scatter."""
    title, steps = _built(action, form)
    ran = job_steps.collapse(steps, title)
    (step,) = [step for step in ran if step["argv"][1:3] == jobs.GENERATE_COMMAND]
    entries = _entries(step)
    args = generate_cli.build_parser()[1]["run-job"].parse_args([step["argv"][4]])
    args.job = entries
    listed = job.stage_list(job.build_parts(entries, args))
    assert [stage.label for stage in listed] == [stage["label"] for stage in step["stages"]]
    assert [stage.skip for stage in listed] == [stage["skipped"] for stage in step["stages"]]
    assert [stage["label"] for stage in step["stages"]].count("Make the block light totals") <= 1


def test_a_new_galaxy_keeps_its_data_check_and_reset_steps_apart():
    title, steps = _built("new_galaxy", {"confirm": DB})
    ran = job_steps.collapse(steps, title)
    assert [step["label"] for step in ran] == [generate_page.DATA_CHECK_LABEL, "Reset the galaxy", "New galaxy"]
    assert ran[0] is steps[0] and ran[1] is steps[1]


def test_the_flags_the_page_adds_to_every_step_are_given_once():
    title, steps = _built("plan", {})
    for step in steps[1:]:
        step["argv"] = step["argv"] + ["--disregard-missing-data", "--debug"]
    (step,) = [step for step in job_steps.collapse(steps, title) if step["argv"][3:4] == ["run-job"]]
    assert step["argv"][5:] == ["--disregard-missing-data", "--debug"]
    assert all(not set(job_steps.SHARED_FLAGS) & set(entry["argv"]) for entry in json.loads(step["argv"][4]))


@pytest.mark.parametrize("action, form", [("bright_band", {"down_to": "100"}), ("reset", {"confirm": DB}),
                                          ("check_db", {})])
def test_a_job_with_nothing_a_run_can_hold_is_left_as_it_is(action, form):
    title, steps = _built(action, form)
    assert job_steps.collapse(steps, title) == steps


def test_a_group_of_one_command_is_already_one_run_and_stays_a_step_of_its_own():
    title, steps = _built("plan", {"skip_scatters": "1"})
    assert job_steps.collapse(steps, title) == steps


def test_a_step_that_is_not_a_stage_of_a_run_ends_the_group_it_stands_in():
    python = "python"
    generate = [python, *jobs.GENERATE_COMMAND]
    steps = [{"label": "a", "argv": generate + ["plan"]}, {"label": "b", "argv": generate + ["scatter", "--mass"]},
             {"label": "band", "argv": generate + ["scatter", "--band", "--down-to", "100"]},
             {"label": "c", "argv": generate + ["link"]}, {"label": "d", "argv": generate + ["settle"]}]
    ran = job_steps.collapse(steps, "T")
    assert [step["label"] for step in ran] == ["Rest of the job", "band", "Rest of the job"]
    assert [entry["command"] for entry in json.loads(ran[0]["argv"][4])] == ["plan", "scatter"]
    assert [entry["command"] for entry in json.loads(ran[2]["argv"][4])] == ["link", "settle"]


def test_a_started_job_holds_the_collapsed_steps(tmp_path):
    title, steps = _built("plan", {})
    job_id = jobs.start_job("plan", title, steps, root=str(tmp_path), spawn=False)
    saved = json.loads((tmp_path / job_id / "job.json").read_text())
    assert [step["label"] for step in saved["steps"]] == [generate_page.DATA_CHECK_LABEL, title]
    assert saved["steps"][1]["argv"][3] == "run-job"
    jobs.release_lock(job_id)


def test_the_estimate_of_a_run_job_step_is_one_figure_for_each_stage(monkeypatch):
    from planetgen.generation import stats
    monkeypatch.setattr(stats, "stage_seconds", lambda conn, key, settings, limit=50: 2.0)
    title, steps = _built("plan", {})
    ran = job_steps.collapse(steps, title)
    per_step, per_stage = job_runner.step_estimates(None, ran, {generate_page.DATA_CHECK_LABEL: 1.0})
    (step,) = [step for step in ran if step["argv"][3:4] == ["run-job"]]
    assert per_step == [1.0, 2.0 * len(step["stages"])]
    assert per_stage == [1.0] + [2.0] * len(step["stages"])


def test_a_stage_with_no_stored_time_makes_the_step_unestimated(monkeypatch):
    from planetgen.generation import stats
    monkeypatch.setattr(stats, "stage_seconds", lambda conn, key, settings, limit=50: None)
    title, steps = _built("plan", {})
    ran = job_steps.collapse(steps, title)
    per_step, per_stage = job_runner.step_estimates(None, ran, {})
    assert per_step[1] is None and None in per_stage
