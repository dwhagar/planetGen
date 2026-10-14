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
    scatter = [s for s in found if s["label"].startswith("Scatter") and "neighborhood" not in s["label"]]
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


# --- PERF.56: each stage's time, stored with the settings it ran with ---------------------

@pytest.fixture
def control_config(mysql_config, monkeypatch):
    """The test database doubling as the control database (as `test_job_tree.py` does)."""
    from planetgen.db import store
    store.get_control_connection(mysql_config, ensure_schema=True).close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    return mysql_config


def _control(mysql_config):
    from planetgen.db import store
    return store.get_control_connection(mysql_config)


def test_each_stage_is_stored_with_its_time_settings_and_what_it_did(control_config, monkeypatch):
    mysql_config = control_config
    from planetgen.db import store
    from planetgen.generation import stats as generation_stats
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "1")
    args = _galaxy("--then-scatter", "--phenomenon-min-mass", "16", "--bright-star-min-luminosity", "9000",
                   "--workers", "2", "--no-settle", "--mysql-host", mysql_config.host,
                   "--mysql-port", str(mysql_config.port), "--mysql-user", mysql_config.user,
                   "--mysql-password", mysql_config.password, "--mysql-database", mysql_config.database)
    stages.begin(stages.galaxy_stages(args), args, "galaxy")
    stages.enter("start")
    stages.enter("mass")
    stages.note(layers=40, layers_modified=7, objects=120)
    stages.enter("luminosity")
    stages.finish()
    conn = _control(mysql_config)
    try:
        history = {run["stage_key"]: run for run in generation_stats.stage_history(conn)}
    finally:
        conn.close()
    assert history["mass"]["settings"] == {"workers": 2, "mass_limit_sol": 16.0}
    assert history["mass"]["metrics"] == {"layers": 40, "layers_modified": 7, "objects": 120}
    assert history["luminosity"]["settings"] == {"workers": 2, "mass_limit_sol": 16.0, "luminosity_floor_sol": 9000.0}
    assert history["mass"]["seconds"] >= 0 and not history["mass"]["skipped"]
    assert history["population"]["skipped"] and "population pass" in history["population"]["skip_reason"]
    assert history["settle"]["skipped"] and history["link"]["skipped"] is True
    assert history["mass"]["stage_n"] == 4 and history["mass"]["stage_total"] == 9


def test_an_estimate_prefers_the_runs_with_the_same_settings(control_config, monkeypatch):
    mysql_config = control_config
    from planetgen.generation import stats as generation_stats
    conn = _control(mysql_config)
    try:
        for seconds, mass in ((10.0, 8.0), (20.0, 8.0), (100.0, 14.0)):
            conn.execute(
                "INSERT INTO generation_stage_runs (database_name, command, stage_key, stage_n, stage_total, label,"
                " skipped, started_at, finished_at, seconds, workers, settings, metrics, version_key)"
                " VALUES ('d', 'galaxy', 'mass', 4, 9, 'x', 0, NOW(6), NOW(6), ?, 1, ?, '{}', ?)",
                (seconds, json.dumps({"workers": 1, "mass_limit_sol": mass}),
                 generation_stats.current_version_key()))
        conn.commit()
        assert generation_stats.stage_seconds(conn, "mass", {"workers": 1, "mass_limit_sol": 8.0}) == 15.0
        assert generation_stats.stage_seconds(conn, "mass", {"workers": 1, "mass_limit_sol": 14.0}) == 100.0
        assert generation_stats.stage_seconds(conn, "mass", {"workers": 4, "mass_limit_sol": 20.0}) is not None
        assert generation_stats.stage_seconds(conn, "settle", {"workers": 1}) is None
    finally:
        conn.close()

def test_the_backfill_stage_is_the_mass_scatter_from_the_neighborhood():
    # GEN.187: the backfill is now the mass rings around the generated sectors.
    (found,) = [stage for stage in stages.galaxy_stages(_galaxy("--then-scatter")) if stage.key == "backfill"]
    assert found.label == "Scatter the massive stars from the neighborhood"
    (off,) = [stage for stage in stages.galaxy_stages(_galaxy("--backfill-from", "none")) if stage.key == "backfill"]
    assert off.skip and "--backfill-from none" in off.skip


def test_a_run_is_estimated_from_its_stages_stored_times_by_their_settings(control_config):
    """Boss: an ETA looks up the stored stats by the settings the stage runs with (PERF.56, PERF.33)."""
    from planetgen.generation import stats as generation_stats
    args = _plan("--phenomenon-min-mass", "8", "--workers", "1")
    conn = _control(control_config)
    try:
        assert stages.estimate_seconds(conn, "plan", args) is None   # nothing recorded yet
        for key, seconds, settings in (("skeleton", 5.0, {"max_ring": 60}), ("phenomena", 7.0, {"mass_limit_sol": 8.0}),
                                       ("mass", 11.0, {"mass_limit_sol": 8.0}),
                                       ("luminosity", 13.0, {"mass_limit_sol": 8.0})):
            conn.execute(
                "INSERT INTO generation_stage_runs (database_name, command, stage_key, stage_n, stage_total, label,"
                " skipped, started_at, finished_at, seconds, workers, settings, metrics, version_key)"
                " VALUES ('d', 'plan', ?, 1, 4, 'x', 0, NOW(6), NOW(6), ?, 1, ?, '{}', ?)",
                (key, seconds, json.dumps(settings), generation_stats.current_version_key()))
        conn.commit()
        assert stages.estimate_seconds(conn, "plan", args) == 36.0
        # A stage that is skipped counts nothing.
        assert stages.estimate_seconds(conn, "plan", _plan("--no-bright-stars")) == 5.0
    finally:
        conn.close()


# --- PERF.55: the whole job's layers per second, and the command line's overall bar -------

def test_the_whole_job_stores_its_layers_per_second_over_all_its_layers(control_config, monkeypatch):
    mysql_config = control_config
    from planetgen.generation import stats as generation_stats
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "1")
    args = _galaxy("--then-scatter", "--phenomenon-min-mass", "16", "--workers", "2", "--mysql-host", mysql_config.host,
                   "--mysql-port", str(mysql_config.port), "--mysql-user", mysql_config.user,
                   "--mysql-password", mysql_config.password, "--mysql-database", mysql_config.database)
    stages.begin(stages.galaxy_stages(args), args, "galaxy")
    stages.enter("mass")
    stages.note(layers=40, layers_modified=0, objects=0)    # a pass that placed nothing still visited its layers
    stages.enter("luminosity")
    stages.note(layers=60, layers_modified=5, objects=9)
    stages.finish()
    conn = _control(mysql_config)
    try:
        (job,) = generation_stats.stage_history(conn, stages.JOB_KEY)
        rate = generation_stats.job_layers_per_second(conn, {"workers": 2, "mass_limit_sol": 16.0})
        fallback = generation_stats.job_layers_per_second(conn, {"workers": 2})
    finally:
        conn.close()
    assert job["metrics"]["layers"] == 100 and job["seconds"] > 0
    assert job["metrics"]["layers_per_second"] == pytest.approx(100 / job["seconds"])
    assert job["settings"]["mass_limit_sol"] == 16.0 and rate == pytest.approx(job["metrics"]["layers_per_second"])
    assert fallback == pytest.approx(rate)


def test_a_run_that_visited_no_layers_stores_no_whole_job_row(control_config, monkeypatch):
    mysql_config = control_config
    from planetgen.generation import stats as generation_stats
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "1")
    args = _galaxy("--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port), "--mysql-user",
                   mysql_config.user, "--mysql-password", mysql_config.password, "--mysql-database",
                   mysql_config.database)
    stages.begin(stages.galaxy_stages(args), args, "galaxy")
    stages.enter("start")
    stages.finish()
    conn = _control(mysql_config)
    try:
        assert generation_stats.stage_history(conn, stages.JOB_KEY) == []
        assert generation_stats.job_layers_per_second(conn) is None
    finally:
        conn.close()


def test_the_overall_bar_is_the_time_so_far_against_the_stored_time_of_every_stage(monkeypatch):
    monkeypatch.delenv(progress_file.ENV_VAR, raising=False)
    found = stages.plan_stages(_plan())
    stages.begin(found)
    try:
        assert stages._state["overall"] is not None
        clock = stages._state["overall"]["t0"]
        description, done, total = stages.overall_progress(clock + 30.0)
        assert description == "Whole job (stage 1 of 4)" and done == pytest.approx(30.0) and total is None
        stages._state["overall"]["estimate"] = 100.0
        assert stages.overall_progress(clock + 30.0)[1:] == (pytest.approx(30.0), 100.0)
        stages.enter("mass")
        assert stages.overall_progress(clock + 30.0)[0] == "Whole job (stage 3 of 4)"
        # Past its estimate the bar waits just short of full rather than claiming to be done.
        _d, done, total = stages.overall_progress(clock + 200.0)
        assert done == pytest.approx(200.0) and total == pytest.approx(200.0 * stages.OVERALL_MARGIN)
        assert done / total < 1.0
    finally:
        stages.finish()
    assert stages.overall_progress() is None


def test_a_single_stage_run_and_a_web_run_have_no_command_line_overall_bar(monkeypatch):
    monkeypatch.delenv(progress_file.ENV_VAR, raising=False)
    stages.begin(stages.plan_stages(_plan("--phenomena-only")))
    assert stages.overall_progress() is None
    stages.finish()
    monkeypatch.setenv(progress_file.ENV_VAR, "/tmp/progress-from-the-web.json")
    stages.begin(stages.plan_stages(_plan()))
    assert stages.overall_progress() is None
    stages.finish()


def test_the_generation_display_carries_the_overall_bar_and_stops_updating_it(monkeypatch):
    from planetgen.generation import run_common
    monkeypatch.delenv(progress_file.ENV_VAR, raising=False)
    stages.begin(stages.plan_stages(_plan()))
    try:
        progress = run_common._generation_progress()
        with progress:
            tasks = [task for task in progress.tasks if task.description.startswith("Whole job")]
            assert len(tasks) == 1 and tasks[0].total is None
            stop = progress.overall_stop
            assert not stop.is_set()
        assert stop.is_set()
    finally:
        stages.finish()
    with run_common._generation_progress() as plain:
        assert plain.tasks == [] and plain.overall_stop is None


# --- PERF.65 / PERF.66: one total for steps and stages, and the whole job's time left --------

def test_a_web_jobs_stages_are_numbered_across_the_whole_job(monkeypatch):
    lines = []
    monkeypatch.setattr(log, "normal", lambda message, *a, **k: lines.append(message))
    monkeypatch.setenv(stages.STAGE_OFFSET_ENV, "3")
    monkeypatch.setenv(stages.STAGE_TOTAL_ENV, "12")
    stages.begin(stages.galaxy_stages(_galaxy("--then-scatter")))
    stages.enter("start")
    stages.enter("settle")
    stages.finish()
    text = "\n".join(lines)
    assert text.startswith("9 stages (stages 4 to 12 of 12): 4. Generate the starting sector;")
    assert "Stage 4 of 12: Generate the starting sector." in text
    assert "Stage 5 of 12: Generate the neighborhood -- skipped: not needed." in text
    assert "Stage 12 of 12: Save the sector paths." in text


def test_a_run_on_its_own_counts_its_own_stages(monkeypatch):
    lines = []
    monkeypatch.setattr(log, "normal", lambda message, *a, **k: lines.append(message))
    monkeypatch.delenv(stages.STAGE_OFFSET_ENV, raising=False)
    monkeypatch.delenv(stages.STAGE_TOTAL_ENV, raising=False)
    stages.begin(stages.plan_stages(_plan("--phenomena-only")))
    stages.enter("phenomena")
    stages.finish()
    assert lines[0] == "1 stage: 1. Scatter the phenomena." and "Stage 1 of 1:" in lines[1]


def test_a_total_too_small_for_the_runs_stages_is_ignored(monkeypatch):
    monkeypatch.setenv(stages.STAGE_OFFSET_ENV, "8")
    monkeypatch.setenv(stages.STAGE_TOTAL_ENV, "9")
    stages.begin(stages.galaxy_stages(_galaxy("--then-scatter")))
    try:
        assert stages._numbering() == (0, 9)
    finally:
        stages.finish()


def test_the_runner_and_the_stages_agree_on_the_variables_that_number_a_jobs_stages():
    from planetgen.web import job_runner
    assert (job_runner.STAGE_OFFSET_ENV, job_runner.STAGE_TOTAL_ENV) == (stages.STAGE_OFFSET_ENV, stages.STAGE_TOTAL_ENV)


def test_the_new_galaxy_log_counts_twelve_tasks_not_four(tmp_path):
    """Boss's sample: a New galaxy job printed "Step 1 of 4" for a job of 12 stages."""
    import sys
    from planetgen.web import job_runner
    steps = _steps("new_galaxy", {"confirm": "db"})
    shown = [{"label": step["label"], "stages": step["stages"],
              "argv": [sys.executable, "-c", "import os; print(os.environ['PLANETGEN_STAGE_OFFSET'], os.environ['PLANETGEN_STAGE_TOTAL'])"]}
             for step in steps]
    counts = job_runner._step_stage_counts(shown)
    assert counts == [1, 1, 1, 9] and sum(counts) == 12
    headings = [job_runner._step_heading(sum(counts[:i]) + 1, count, sum(counts)) for i, count in enumerate(counts)]
    assert headings == ["Step 1 of 12", "Step 2 of 12", "Step 3 of 12", "Steps 4 to 12 of 12"]
    job_dir = tmp_path / "jobs" / "20261010-120000-abcd"
    job_dir.mkdir(parents=True)
    (job_dir / "job.json").write_text(json.dumps({
        "id": job_dir.name, "kind": "new_galaxy", "title": "New galaxy", "database": "db", "created_at": 0,
        "cwd": str(tmp_path), "env": {}, "steps": shown}))
    assert job_runner.run(str(job_dir)) == 0
    log_text = (job_dir / "output.log").read_text()
    assert "=== Step 1 of 12: " in log_text and "=== Step 3 of 12: " in log_text
    assert "=== Steps 4 to 12 of 12: " in log_text and "of 4" not in log_text
    assert "=== Steps 4 to 12 exited with status 0" in log_text
    assert "\n3 12\n" in log_text      # the last step's process was told its stages start at 4 of 12


def test_the_command_line_bar_never_shows_the_running_stages_time_as_the_whole_jobs_until_the_last(monkeypatch):
    monkeypatch.delenv(progress_file.ENV_VAR, raising=False)
    stages.begin(stages.plan_stages(_plan()))
    try:
        clock = stages._state["overall"]["t0"]
        assert stages.overall_progress(clock + 30.0)[2] is None          # nothing finished, nothing to average
        stages.enter("skeleton")
        stages.enter("phenomena")
        stages._state["job_seconds"], stages._state["ran"] = 20.0, 1       # one stage took 20 s
        # Stage 2 of 4 is 25 s old: the time left counts 2 more stages at the 20 s average.
        total = stages.overall_progress(clock + 45.0)[2]
        assert total == pytest.approx(20.0 + 25.0 + 2 * 20.0)
        stages.enter("luminosity")
        stages._state["job_seconds"], stages._state["ran"] = 60.0, 3
        # The last stage, 10 s in: nothing after it, so the time left is its own.
        total = stages.overall_progress(clock + 70.0)[2]
        assert total == pytest.approx(max(60.0 + max(20.0, 10.0), 70.0 * stages.OVERALL_MARGIN))
    finally:
        stages.finish()
