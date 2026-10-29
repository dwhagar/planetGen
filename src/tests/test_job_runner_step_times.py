# tests/test_job_runner_step_times.py

"""
What a Generate page job keeps about its steps now that the job tree is gone (PERF.98, PERF.100): the time each step
took, by its label, in `step-times.json` of the jobs folder, which the next job plans its estimates by; and the command
line a run records (`generate._run_argv`, without the database and debug options).

Tests that take the `mysql_config` fixture (see `conftest.py`) are skipped, not failed, when no MySQL test server is
configured/reachable.
"""

import json
import sys
import time

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.web import job_runner as jobRunner


@pytest.fixture
def control_config(mysql_config, monkeypatch):
    store.get_control_connection(mysql_config, ensure_schema=True).close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    return mysql_config


def test_the_run_argv_leaves_out_the_database_and_debug_options():
    argv = ["galaxy", "--shell", "3", "--mysql-password", "secret", "--debug", "f.log", "--band=1",
            "--mysql-user=u", "--mysql-database", "db", "--debug", "--yes"]
    assert generate_cli._run_argv(argv) == ["galaxy", "--shell", "3", "--band=1", "--yes"]


def test_a_web_job_keeps_the_time_of_each_step_in_the_jobs_directory(control_config, tmp_path):
    jobs_dir = tmp_path / "jobs"
    job_dir = jobs_dir / "20261001-120000-abcd"
    job_dir.mkdir(parents=True)
    (job_dir / "job.json").write_text(json.dumps({
        "id": job_dir.name, "kind": "galaxy", "title": "Generate sectors", "database": control_config.database,
        "created_at": time.time(), "cwd": str(tmp_path), "env": {},
        "steps": [{"label": "Step one", "argv": [sys.executable, "-c", "pass"]},
                  {"label": "Step two", "argv": [sys.executable, "-c", "pass"]}],
    }))
    assert jobRunner.run(str(job_dir)) == 0
    state = json.loads((job_dir / "state.json").read_text())
    assert state["step_estimates"] == [None, None] and state["step_started_at"] > 0
    recorded = jobRunner.recorded_step_seconds(str(jobs_dir))
    assert set(recorded) == {"Step one", "Step two"} and all(seconds > 0 for seconds in recorded.values())
    # A second run of the same steps now has their recorded times to plan by.
    job_dir.joinpath("state.json").unlink()
    assert jobRunner.run(str(job_dir)) == 0
    estimates = json.loads((job_dir / "state.json").read_text())["step_estimates"]
    assert all(seconds is not None and seconds > 0 for seconds in estimates)


def test_no_record_of_step_times_is_no_estimate(tmp_path):
    assert jobRunner.recorded_step_seconds(str(tmp_path)) == {}
    (tmp_path / jobRunner.STEP_TIMES_NAME).write_text("not json")
    assert jobRunner.recorded_step_seconds(str(tmp_path)) == {}
