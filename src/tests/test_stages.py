# tests/test_stages.py

"""
UX.89: every staged action numbers and announces its stages, the optional ones
as skipped with a reason, and the Generate page counts the stages of a whole
job instead of its processes.
"""

import json
import sys
import types

import time

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.generation import stages
from planetgen.queue import progress_file
from planetgen.util import log
from planetgen.web import generate_page, jobs


def _galaxy(*flags):
    _parser, parsers = generate_cli.build_parser()
    return parsers["galaxy"].parse_args(list(flags))


def _sectors(*flags):
    _parser, parsers = generate_cli.build_parser()
    return parsers["sectors"].parse_args(list(flags))


def _plan(*flags):
    _parser, parsers = generate_cli.build_parser()
    return parsers["plan"].parse_args(list(flags))


def _scatter(*flags):
    _parser, parsers = generate_cli.build_parser()
    return parsers["scatter"].parse_args(list(flags))


def _four_stages():
    """The skeleton and the three scatters, the list a `plan` run once had (the stage code is the same for any list)."""
    return [stages.plan_stages(_plan())[0], stages.SCATTER_STAGE["phenomena"], stages.SCATTER_STAGE["mass"],
            stages.SCATTER_STAGE["luminosity"]]


def _labels(found):
    return [(stage.key, stage.skip is not None) for stage in found]


def test_a_random_start_galaxy_run_has_fourteen_stages_and_skips_none_of_them():
    found = stages.galaxy_stages(_galaxy())
    assert _labels(found) == [
        ("start", False), ("neighborhood", False), ("link", False), ("mass", False), ("luminosity", False),
        ("phenomena", False), ("phenomena-index", False), ("nebula-draw", False), ("block-lights", False),
        ("backfill", False), ("population", False), ("territories", False), ("settle", False),
        ("run-blocks", False)]


def test_a_galaxy_run_of_one_sector_and_its_neighborhood_has_both_stages():
    found = stages.galaxy_stages(_galaxy("--shell", "3", "--band", "0", "--wedge", "1", "--radius-pc", "20"))
    assert [stage.key for stage in found][:3] == ["sectors", "neighborhood", "link"] and len(found) == 14
    assert [stage.key for stage in stages.galaxy_stages(_galaxy("--center-sector", "1", "--radius-pc", "20"))][:2] \
        == ["sectors", "link"]


def test_a_sectors_run_has_the_sector_stages_alone():
    assert [stage.key for stage in stages.for_command("sectors", _sectors())] == ["start", "neighborhood"]


def test_a_galaxy_run_with_an_address_has_one_sector_stage():
    found = stages.galaxy_stages(_galaxy("--shell", "3", "--band", "0"))
    assert [stage.key for stage in found][:2] == ["sectors", "link"] and len(found) == 13


def test_plan_is_the_skeleton_and_each_scatter_and_fill_is_one_stage_of_its_own():
    assert [s.key for s in stages.plan_stages(_plan())] == ["skeleton"]
    for flag, keys in (("--mass", ["mass", "block-lights"]), ("--luminosity", ["luminosity", "block-lights"]),
                       ("--phenomena", ["phenomena", "phenomena-index", "nebula-draw", "block-lights"]),
                       ("--nebula-field", ["nebula-draw"]), ("--band", ["band"])):
        assert [s.key for s in stages.scatter_stages(_scatter(flag))] == keys, flag
        assert [s.key for s in stages.for_command("scatter", _scatter(flag))] == keys
    assert stages.scatter_stages(_scatter("--mass")) != stages.scatter_stages(_scatter("--luminosity"))
    for command, keys in (("fill-map-luminosity", ["map-luminosity", "block-lights"]), ("rescore-habitability", ["habitability"]),
                          ("clean-nebula-field", ["nebula-field"])):
        assert [s.key for s in stages.for_command(command, None)] == keys
    for kind in stages.CLEAR_KINDS:
        found = stages.clear_stages(types.SimpleNamespace(**{k.replace("-", "_"): k == kind for k in stages.CLEAR_KINDS}))
        assert [s.key for s in found] == [f"clear-{kind}"]
        assert stages.stage_settings(f"clear-{kind}", None) == {"part": kind}
    assert stages.for_command("benchmark", None) == []


def test_entering_a_stage_announces_the_skipped_ones_before_it_and_the_rest_at_the_end(monkeypatch, tmp_path):
    lines = []
    monkeypatch.setattr(log, "normal", lambda message, *a, **k: lines.append(message))
    monkeypatch.setenv(progress_file.ENV_VAR, str(tmp_path / "progress.json"))
    stages.begin(stages.galaxy_stages(_galaxy()))
    stages.enter("start")
    stages.enter("link")           # "neighborhood" never ran
    stages.skip("mass", "no layers")
    stages.enter("settle")         # the ones between are skipped
    stages.finish()
    text = "\n".join(lines)
    assert text.startswith("Stage 1 of 14: Generate the starting sector")
    assert "stages:" not in text and "stage:" not in text      # PERF.145: no list of the stages up front
    assert "Stage 2 of 14: Generate the neighborhood -- skipped: not needed" in text
    assert "Stage 4 of 14: Scatter the massive stars -- skipped: no layers" in text
    assert "Stage 11 of 14: Run the population pass -- skipped: not needed" in text
    assert "Stage 12 of 14: Recompute the territories -- skipped: not needed" in text
    assert "Stage 13 of 14: Save the sector paths" in text
    assert "Stage 14 of 14: Update the map's block totals -- skipped: not needed" in text
    assert json.loads((tmp_path / "progress.json").read_text())["stage"]["index"] == 14


def test_a_stage_key_the_run_does_not_have_is_ignored():
    stages.begin(stages.scatter_stages(_scatter("--phenomena")))
    stages.enter("mass")
    stages.finish()


FORM = {"confirm": "db", "max_ring": "", "skip_scatters": ""}


def _steps(action, form):
    return generate_page.build_job(action, form, "db")[2]


def test_the_new_galaxy_job_counts_every_stage_of_every_step():
    steps = _steps("new_galaxy", {"confirm": "db"})
    found = jobs.job_stages(steps)
    assert [s["label"] for s in found][:3] == [generate_page.DATA_CHECK_LABEL, "Reset the galaxy", "Plan the galaxy"]
    assert len(found) == 16 and found[-1]["n"] == 16
    assert not any(s["skipped"] for s in found)


def test_ticking_skip_the_scatter_leaves_the_scatter_steps_out():
    found = jobs.job_stages(_steps("new_galaxy", {"confirm": "db", "skip_scatters": "on"}))
    scatter = [s for s in found if s["label"].startswith("Scatter") and "neighborhood" not in s["label"]]
    assert scatter == [] and len(found) == 8


def test_a_plan_job_lists_its_scatter_steps_stages():
    found = jobs.job_stages(_steps("plan", {"confirm": "db"}))
    assert [s["label"] for s in found][-8:] == [
        "Scatter the massive stars", "Make the block light totals", "Scatter the bright stars",
        "Make the block light totals", "Scatter the phenomena", "Index the phenomena", "Draw the nebula field",
        "Make the block light totals"]


def test_the_stage_view_names_the_running_stage_by_its_place_in_the_whole_job():
    steps = _steps("new_galaxy", {"confirm": "db"})
    job = {"stages": jobs.job_stages(steps), "step": 6, "finished": False,
           "progress": {"stage": {"index": 1, "total": 2, "label": "Scatter the massive stars", "skipped": None}}}
    view = generate_page.stage_view(job)
    assert view["stage_text"] == "Stage 7 of 16: Scatter the massive stars"
    states = [s["state"] for s in view["stage_list"]]
    assert states[:6] == ["done"] * 6 and states[6] == "running" and "waiting" in states
    assert "skipped" not in states


# --- PERF.56: each stage's time, stored with the settings it ran with ---------------------

@pytest.fixture
def control_config(mysql_config, monkeypatch):
    """The test database doubling as the control database (as `test_job_runner_step_times.py` does)."""
    from planetgen.db import store
    store.get_control_connection(mysql_config, ensure_schema=True).close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    return mysql_config


def _control(mysql_config):
    from planetgen.db import store
    return store.get_control_connection(mysql_config)


def test_each_stage_is_stored_with_its_time_settings_and_what_it_did(control_config, monkeypatch):
    mysql_config = control_config
    from planetgen.generation import stats as generation_stats
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "1")
    args = _galaxy("--phenomenon-min-mass", "16", "--bright-star-min-luminosity", "9000",
                   "--workers", "2", "--mysql-host", mysql_config.host,
                   "--mysql-port", str(mysql_config.port), "--mysql-user", mysql_config.user,
                   "--mysql-password", mysql_config.password, "--mysql-database", mysql_config.database)
    stages.begin(stages.galaxy_stages(args), args, "galaxy")
    stages.enter("start")
    stages.enter("mass")
    stages.note(shells=40, shells_modified=7, objects=120)
    stages.enter("luminosity")
    stages.finish()
    conn = _control(mysql_config)
    try:
        history = {run["stage_key"]: run for run in generation_stats.stage_history(conn)}
    finally:
        conn.close()
    assert history["mass"]["settings"] == {"workers": 2, "mass_limit_sol": 16.0}
    metrics = dict(history["mass"]["metrics"])
    assert metrics.pop("seconds_per_shell") == pytest.approx(history["mass"]["seconds"] / 40)
    assert metrics == {"shells": 40, "shells_modified": 7, "objects": 120}
    assert history["luminosity"]["settings"] == {"workers": 2, "mass_limit_sol": 16.0, "luminosity_floor_sol": 9000.0}
    assert history["mass"]["seconds"] >= 0 and not history["mass"]["skipped"]
    assert history["population"]["skipped"] and history["territories"]["skipped"] and history["settle"]["skipped"] and history["link"]["skipped"] is True
    assert history["mass"]["stage_n"] == 4 and history["mass"]["stage_total"] == 14


def test_a_stage_handed_to_the_queue_is_not_stored_by_this_process(control_config, monkeypatch):
    mysql_config = control_config
    from planetgen.generation import stats as generation_stats
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "1")
    args = _scatter("--mass", "--workers", "2", "--mysql-host", mysql_config.host,
                 "--mysql-port", str(mysql_config.port), "--mysql-user", mysql_config.user,
                 "--mysql-password", mysql_config.password, "--mysql-database", mysql_config.database)
    stages.begin([stages.SCATTER_STAGE["mass"], stages.SCATTER_STAGE["luminosity"]], args, "scatter")
    stages.enter("mass")
    row = stages.hand_over_open_stage()
    assert row["command"] == "scatter" and row["settings"]["workers"] == 2 and (row["n"], row["total"]) == (1, 2)
    stages.note(shells=10, shells_modified=1, objects=3)
    stages.enter("luminosity")
    stages.finish()
    conn = _control(mysql_config)
    try:
        assert [run["stage_key"] for run in generation_stats.stage_history(conn)
                if run["stage_key"] in ("mass", "luminosity")] == ["luminosity"]
        (job,) = generation_stats.stage_history(conn, stages.JOB_KEY)
    finally:
        conn.close()
    assert job["metrics"]["shells"] == 10, "the whole job still counts the handed-over stage"


def test_a_stage_with_nothing_recorded_hands_over_no_row(monkeypatch):
    from planetgen.generation import stats as generation_stats
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "0")
    args = _scatter("--mass")
    stages.begin(stages.for_command("scatter", args), args, "scatter")
    stages.enter("mass")
    assert stages.hand_over_open_stage() is None
    stages.finish()


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
    (found,) = [stage for stage in stages.galaxy_stages(_galaxy()) if stage.key == "backfill"]
    assert found.label == "Scatter the massive stars from the neighborhood" and found.skip is None


def test_a_run_is_estimated_from_its_stages_stored_times_by_their_settings(control_config):
    """Boss: an ETA looks up the stored stats by the settings the stage runs with (PERF.56, PERF.33)."""
    from planetgen.generation import stats as generation_stats
    runs = (("plan", _plan("--workers", "1")),
            ("scatter", _scatter("--phenomena", "--workers", "1")),
            ("scatter", _scatter("--mass", "--phenomenon-min-mass", "8", "--workers", "1")),
            ("scatter", _scatter("--luminosity", "--workers", "1")))
    conn = _control(control_config)
    try:
        for command, args in runs:
            assert stages.estimate_seconds(conn, command, args) is None   # nothing recorded yet
        for key, seconds, settings in (("skeleton", 5.0, {"max_ring": 60}), ("phenomena", 7.0, {"mass_limit_sol": 8.0}),
                                       ("phenomena-index", 1.0, {}), ("nebula-draw", 2.0, {}), ("block-lights", 3.0, {}),
                                       ("mass", 11.0, {"mass_limit_sol": 8.0}),
                                       ("luminosity", 13.0, {"mass_limit_sol": 8.0})):
            conn.execute(
                "INSERT INTO generation_stage_runs (database_name, command, stage_key, stage_n, stage_total, label,"
                " skipped, started_at, finished_at, seconds, workers, settings, metrics, version_key)"
                " VALUES ('d', 'plan', ?, 1, 4, 'x', 0, NOW(6), NOW(6), ?, 1, ?, '{}', ?)",
                (key, seconds, json.dumps(settings), generation_stats.current_version_key()))
        conn.commit()
        assert [stages.estimate_seconds(conn, command, args) for command, args in runs] == [5.0, 13.0, 14.0, 16.0]
        # A stage with no time stored makes the estimate unknown.
        assert stages.estimate_seconds(conn, "scatter", _scatter("--band")) is None
    finally:
        conn.close()


# --- PERF.55: the whole job's layers per second, and the command line's overall bar -------

def test_the_whole_job_stores_its_shells_per_second_over_all_its_layers(control_config, monkeypatch):
    mysql_config = control_config
    from planetgen.generation import stats as generation_stats
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "1")
    args = _galaxy("--phenomenon-min-mass", "16", "--workers", "2", "--mysql-host", mysql_config.host,
                   "--mysql-port", str(mysql_config.port), "--mysql-user", mysql_config.user,
                   "--mysql-password", mysql_config.password, "--mysql-database", mysql_config.database)
    stages.begin(stages.galaxy_stages(args), args, "galaxy")
    stages.enter("mass")
    stages.note(shells=40, shells_modified=0, objects=0)    # a pass that placed nothing still visited its layers
    stages.enter("luminosity")
    stages.note(shells=60, shells_modified=5, objects=9)
    stages.finish()
    conn = _control(mysql_config)
    try:
        (job,) = generation_stats.stage_history(conn, stages.JOB_KEY)
        rate = generation_stats.job_shells_per_second(conn, {"workers": 2, "mass_limit_sol": 16.0})
        fallback = generation_stats.job_shells_per_second(conn, {"workers": 2})
    finally:
        conn.close()
    assert job["metrics"]["shells"] == 100 and job["seconds"] > 0
    assert job["metrics"]["shells_per_second"] == pytest.approx(100 / job["seconds"])
    assert job["settings"]["mass_limit_sol"] == 16.0 and rate == pytest.approx(job["metrics"]["shells_per_second"])
    assert fallback == pytest.approx(rate)


def test_a_run_that_visited_no_layers_still_stores_its_whole_process_row(control_config, monkeypatch):
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
        (job,) = generation_stats.stage_history(conn, stages.JOB_KEY)
        assert job["metrics"] == {"stages_run": 1, "shells": 0} and job["seconds"] >= 0
        assert generation_stats.stage_history(conn, stages.SCATTER_KEY) == []
        assert generation_stats.job_shells_per_second(conn) is None
    finally:
        conn.close()


def test_the_scatters_are_stored_on_their_own_and_added_up_for_a_scatter_with_no_record(control_config, monkeypatch):
    mysql_config = control_config
    from planetgen.generation import stats as generation_stats
    monkeypatch.setenv(generation_stats.STATS_ENV_VAR, "1")
    args = _galaxy("--workers", "2", "--mysql-host", mysql_config.host,
                   "--mysql-port", str(mysql_config.port), "--mysql-user", mysql_config.user,
                   "--mysql-password", mysql_config.password, "--mysql-database", mysql_config.database)
    stages.begin(stages.galaxy_stages(args), args, "galaxy")
    stages.enter("mass")
    stages.note(shells=40, objects=10)
    time.sleep(0.02)
    stages.enter("luminosity")
    stages.note(shells=60, objects=5)
    time.sleep(0.02)
    stages.finish()
    conn = _control(mysql_config)
    try:
        (scatter,) = generation_stats.stage_history(conn, stages.SCATTER_KEY)
        mass = generation_stats.stage_seconds(conn, "mass", {"workers": 2})
        luminosity = generation_stats.stage_seconds(conn, "luminosity", {"workers": 2})
        phenomena = generation_stats.stage_seconds(conn, "phenomena", {"workers": 2})
        link = generation_stats.stage_seconds(conn, "link", {"workers": 2})
    finally:
        conn.close()
    assert scatter["metrics"]["operations"] == 2 and scatter["metrics"]["shells"] == 100
    assert scatter["metrics"]["objects"] == 15
    assert scatter["metrics"]["seconds_per_shell"] == pytest.approx(scatter["seconds"] / 100)
    assert mass != luminosity                      # each kind keeps its own figure
    assert phenomena == pytest.approx(scatter["seconds"] / 2)    # a scatter never run takes the common one
    assert link is None                            # any other stage has no common figure


def test_the_overall_bar_is_the_time_so_far_against_the_stored_time_of_every_stage(monkeypatch):
    monkeypatch.delenv(progress_file.ENV_VAR, raising=False)
    found = _four_stages()
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
    stages.begin(stages.scatter_stages(_scatter("--nebula-field")))
    assert stages.overall_progress() is None
    stages.finish()
    monkeypatch.setenv(progress_file.ENV_VAR, "/tmp/progress-from-the-web.json")
    stages.begin(_four_stages())
    assert stages.overall_progress() is None
    stages.finish()


def test_the_generation_display_carries_the_overall_bar_and_stops_updating_it(monkeypatch):
    from planetgen.generation import run_common
    monkeypatch.delenv(progress_file.ENV_VAR, raising=False)
    stages.begin(_four_stages())
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
    monkeypatch.setenv(stages.STAGE_TOTAL_ENV, "17")
    stages.begin(stages.galaxy_stages(_galaxy()))
    stages.enter("start")
    stages.enter("settle")
    stages.finish()
    text = "\n".join(lines)
    assert text.startswith("Stage 4 of 17: Generate the starting sector")      # a step of 14 stages: the banner names the step
    assert "stages (stages" not in text
    assert "Stage 5 of 17: Generate the neighborhood -- skipped: not needed" in text
    assert "Stage 16 of 17: Save the sector paths" in text


def test_a_run_on_its_own_counts_its_own_stages(monkeypatch):
    lines = []
    monkeypatch.setattr(log, "normal", lambda message, *a, **k: lines.append(message))
    monkeypatch.delenv(stages.STAGE_OFFSET_ENV, raising=False)
    monkeypatch.delenv(stages.STAGE_TOTAL_ENV, raising=False)
    stages.begin(stages.scatter_stages(_scatter("--nebula-field")))
    stages.enter("nebula-draw")
    stages.finish()
    assert lines == ["Stage 1 of 1: Draw the nebula field",         # no banner here, so these two are the mentions
                     "Stage 1 of 1 finished after 0 s"]               # UX.122: and the stage says it finished


def test_a_step_of_one_stage_under_the_runners_banner_logs_nothing_for_it_but_a_skipped_stage_says_why(monkeypatch):
    """PERF.145: "=== Stage 4 of 13: Plan the galaxy ===" is the one mention of the stage."""
    lines = []
    monkeypatch.setattr(log, "normal", lambda message, *a, **k: lines.append(message))
    monkeypatch.setenv(stages.STAGE_OFFSET_ENV, "3")
    monkeypatch.setenv(stages.STAGE_TOTAL_ENV, "13")
    stages.begin(stages.scatter_stages(_scatter("--nebula-field")))
    stages.enter("nebula-draw")
    stages.finish()
    assert lines == []
    stages.begin(stages.scatter_stages(_scatter("--nebula-field")))
    stages.skip("nebula-draw", "nothing to scatter")
    stages.finish()
    assert lines == ["Stage 4 of 13: Draw the nebula field -- skipped: nothing to scatter"]


def test_every_stage_of_every_command_is_logged_at_most_once(monkeypatch):
    """PERF.145: no stage appears twice in a run's log, whichever command, with or without the runner's numbering."""
    for command, make in (("galaxy", lambda: stages.galaxy_stages(_galaxy())),
                          ("plan", lambda: stages.plan_stages(_plan()))):
        for env in ({}, {stages.STAGE_OFFSET_ENV: "2", stages.STAGE_TOTAL_ENV: "30"}):
            lines = []
            monkeypatch.setattr(log, "normal", lambda message, *a, **k: lines.append(message))
            for name in (stages.STAGE_OFFSET_ENV, stages.STAGE_TOTAL_ENV):
                monkeypatch.delenv(name, raising=False)
            for name, value in env.items():
                monkeypatch.setenv(name, value)
            listed = make()
            stages.begin(listed, command=command)
            for stage in listed:
                stages.enter(stage.key)
            stages.finish()
            ends = [line for line in lines if " finished after " in line]
            lines = [line for line in lines if line not in ends]
            numbers = [line.split(":")[0] for line in lines]
            assert len(numbers) == len(set(numbers)) and all(line.startswith("Stage ") for line in lines)
            # a one-stage command under the runner's numbering is named by the runner's banner alone
            assert len(lines) == (0 if len(listed) == 1 and env else len(listed))
            assert len(ends) == len(lines)          # UX.122: and every stage it names also says it finished, once


def test_a_total_too_small_for_the_runs_stages_is_ignored(monkeypatch):
    monkeypatch.setenv(stages.STAGE_OFFSET_ENV, "8")
    monkeypatch.setenv(stages.STAGE_TOTAL_ENV, "9")
    stages.begin(stages.galaxy_stages(_galaxy()))
    try:
        assert stages._numbering() == (0, 14)
    finally:
        stages.finish()


def test_the_runner_and_the_stages_agree_on_the_variables_that_number_a_jobs_stages():
    from planetgen.web import job_runner
    assert (job_runner.STAGE_OFFSET_ENV, job_runner.STAGE_TOTAL_ENV) == (stages.STAGE_OFFSET_ENV, stages.STAGE_TOTAL_ENV)


def test_the_new_galaxy_log_counts_every_stage_not_every_step(tmp_path):
    """Boss's sample: a New galaxy job printed "Step 1 of 4" for a job of 16 stages (16 stages and 10 steps with the data check, GEN.590, PERF.172)."""
    import sys
    from planetgen.web import job_runner
    steps = _steps("new_galaxy", {"confirm": "db"})
    shown = [{"label": step["label"], "stages": step["stages"],
              "argv": [sys.executable, "-c", "import os; print(os.environ['PLANETGEN_STAGE_OFFSET'], os.environ['PLANETGEN_STAGE_TOTAL'])"]}
             for step in steps]
    counts = job_runner._step_stage_counts(shown)
    assert counts == [1, 1, 1, 2, 1, 2, 2, 4, 1, 1] and sum(counts) == 16
    job_dir = tmp_path / "jobs" / "20261010-120000-abcd"
    job_dir.mkdir(parents=True)
    (job_dir / "job.json").write_text(json.dumps({
        "id": job_dir.name, "kind": "new_galaxy", "title": "New galaxy", "database": "db", "created_at": 0,
        "cwd": str(tmp_path), "env": {}, "steps": shown}))
    assert job_runner.run(str(job_dir)) == 0
    log_text = (job_dir / "output.log").read_text()
    assert "=== Stage 1 of 16: " in log_text and "=== Stage 3 of 16: " in log_text
    assert "=== Stages 4 to 5 of 16: " in log_text and "of 10" not in log_text
    assert "=== Stages 4 to 5 of 16 finished after " in log_text
    from planetgen.util import logrecord
    records = [logrecord.parse_file_line(line) for line in (job_dir / "records.log").read_text().splitlines()]
    # UX.123: the page and the download hold the display line (with its time) of each stage mark, the records file the record
    assert [logrecord.display_line(r) for r in records if r.event] == [line for line in log_text.splitlines() if " job: === " in line]
    assert sum(1 for r in records if r.event == logrecord.STAGE_START) == 10
    assert "\n3 16\n" in log_text      # the sectors step's process was told its stages start at 4 of 16


def test_the_command_line_bar_never_shows_the_running_stages_time_as_the_whole_jobs_until_the_last(monkeypatch):
    monkeypatch.delenv(progress_file.ENV_VAR, raising=False)
    stages.begin(_four_stages())
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


def test_every_stage_marker_comes_from_one_formatter():
    """UX.98: one wording for a single stage and for a run, in the section headers and the closing line."""
    from planetgen.util import stagemark
    assert stagemark.marker(4, 12) == "Stage 4 of 12"
    assert stagemark.marker(4, 12, 9) == "Stages 4 to 12 of 12"
    assert stagemark.section(1, 12, 1, "Reset") == "=== Stage 1 of 12: Reset ==="
    assert stagemark.section(4, 12, 9, "Generate") == "=== Stages 4 to 12 of 12: Generate ==="
    from planetgen.util import logrecord
    done = logrecord.stage_record(logrecord.STAGE_FINISH, "job", 4, 12, seconds=30.6, count=9)
    failed = logrecord.stage_record(logrecord.STAGE_FAILED, "job", 2, 12, seconds=4, reason="exited with status 1")
    assert logrecord.plain_line(done) == "=== Stages 4 to 12 of 12 finished after 31 s ==="
    assert logrecord.plain_line(failed) == "=== Stage 2 of 12 failed after 4 s: exited with status 1 ==="


def test_a_stage_announcement_is_one_divider_on_the_console_and_one_record_in_the_record_file(tmp_path, monkeypatch, capsys):
    """UX.122: the same stage mark is the plain divider the console shows and the record the record file holds."""
    from planetgen.util import logrecord
    monkeypatch.delenv(stages.STAGE_OFFSET_ENV, raising=False)
    monkeypatch.delenv(stages.STAGE_TOTAL_ENV, raising=False)
    path = tmp_path / "records.log"
    log.set_record_file(str(path))
    log.configure(log.NORMAL)
    try:
        stages.begin(stages.scatter_stages(_scatter("--nebula-field")))
        stages.enter("nebula-draw")
        stages.finish()
        for handler in log._logger.handlers:
            handler.flush()
    finally:
        log.set_record_file(None)
        log.configure(log.NORMAL)
    records = [logrecord.parse_file_line(line) for line in path.read_text().splitlines()]
    assert [(r.event, r.message) for r in records] == [("stage.start", "Stage 1 of 1: Draw the nebula field"),
                                                       ("stage.finish", "Stage 1 of 1 finished after 0 s")]
    assert records[0].fields == {"event": "stage.start", "stage": 1, "stages": 1}
    assert records[1].fields == {"event": "stage.finish", "stage": 1, "stages": 1, "seconds": 0.0}
    assert capsys.readouterr().out.splitlines() == ["=== Stage 1 of 1: Draw the nebula field ===",
                                                    "=== Stage 1 of 1 finished after 0 s ==="]


def test_a_stage_that_fails_says_so_in_the_console_and_the_record_file(tmp_path, monkeypatch, capsys):
    """UX.122: an error ends the open stage with its failure line, as the job runner's banner does for a web job's step."""
    from planetgen.util import logrecord
    monkeypatch.delenv(stages.STAGE_OFFSET_ENV, raising=False)
    monkeypatch.delenv(stages.STAGE_TOTAL_ENV, raising=False)
    path = tmp_path / "records.log"
    log.set_record_file(str(path))
    log.configure(log.NORMAL)
    try:
        stages.begin(stages.scatter_stages(_scatter("--nebula-field")))
        stages.enter("nebula-draw")
        stages.fail("MySQL went away")
        stages.fail("again")                                  # nothing is open any more: nothing more is said
        for handler in log._logger.handlers:
            handler.flush()
    finally:
        log.set_record_file(None)
        log.configure(log.NORMAL)
    records = [logrecord.parse_file_line(line) for line in path.read_text().splitlines()]
    assert [(r.event, r.severity, r.message) for r in records] == [
        ("stage.start", "info", "Stage 1 of 1: Draw the nebula field"),
        ("stage.failed", "error", "Stage 1 of 1 failed after 0 s: MySQL went away")]
    assert capsys.readouterr().out.splitlines() == ["=== Stage 1 of 1: Draw the nebula field ===",
                                                    "=== Stage 1 of 1 failed after 0 s: MySQL went away ==="]
    assert stages._state["stages"] == [] and stages._state["open"] is None


def test_a_command_that_raises_ends_its_open_stage_as_failed(monkeypatch, capsys):
    """UX.122: `planetgen` closes the run's open stage with the error's name and message."""
    from planetgen.cli import generate

    def handler(args):
        stages.begin(stages.scatter_stages(_scatter("--nebula-field")))
        stages.enter("nebula-draw")
        raise RuntimeError("the sky fell")

    monkeypatch.setitem(generate._COMMAND_HANDLERS, "check-math", handler)
    monkeypatch.setattr(sys, "argv", ["planetgen", "check-math"])
    monkeypatch.setattr(generate, "_wait_for_enter", lambda: None)
    monkeypatch.delenv(stages.STAGE_OFFSET_ENV, raising=False)
    monkeypatch.delenv(stages.STAGE_TOTAL_ENV, raising=False)
    try:
        with pytest.raises(RuntimeError):
            generate.main()
    finally:
        log.configure(log.NORMAL)
    assert "=== Stage 1 of 1 failed after 0 s: RuntimeError: the sky fell ===" in capsys.readouterr().out
