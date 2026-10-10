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
