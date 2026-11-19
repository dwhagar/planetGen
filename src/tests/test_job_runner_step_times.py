# tests/test_job_runner_step_times.py

"""
What a Generate page job keeps about its steps now that the job tree is gone (PERF.98, PERF.100): the time each step
took, by its label, in `step-times.json` of the jobs folder, which the next job plans its estimates by; and the command
line a run records (`generate._run_argv`, without the database and debug options).

Tests that take the `mysql_config` fixture (see `conftest.py`) are skipped, not failed, when no MySQL test server is
configured/reachable.
"""

import sys

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.queue import pagerun
from planetgen.web import job_runner as jobRunner, jobs
from tests import page_jobs


def test_the_run_argv_leaves_out_the_database_and_debug_options():
    argv = ["galaxy", "--shell", "3", "--mysql-password", "secret", "--debug", "f.log", "--band=1",
            "--mysql-user=u", "--mysql-database", "db", "--debug", "--yes"]
    assert generate_cli._run_argv(argv) == ["galaxy", "--shell", "3", "--band=1", "--yes"]


def test_a_web_job_keeps_the_time_of_each_step_in_the_jobs_directory(control_config, tmp_path, monkeypatch):
    jobs_dir = tmp_path / "jobs"
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(jobs_dir))
    page_jobs.forget_all_jobs()
    steps = [{"label": "Step one", "argv": [sys.executable, "-c", "pass"]}, {"label": "Step two", "argv": [sys.executable, "-c", "pass"]}]
    try:
        job_id = jobs.start_job("galaxy", "Generate sectors", steps, database=control_config.database, root=str(jobs_dir), spawn=False)
        assert page_jobs.run_to_end(job_id) == 0
        with page_jobs.ledger_rows() as control:
            job = pagerun.read(control, job_id)
        assert job["step_estimates"] == [None, None] and job["step_started_at"] > 0
        recorded = jobRunner.recorded_step_seconds(str(jobs_dir))
        assert set(recorded) == {"Step one", "Step two"} and all(seconds > 0 for seconds in recorded.values())
        # A second job of the same steps now has their recorded times to plan by.
        second = jobs.start_job("galaxy", "Generate sectors", steps, database=control_config.database, root=str(jobs_dir), spawn=False)
        assert page_jobs.run_to_end(second) == 0
        with page_jobs.ledger_rows() as control:
            estimates = pagerun.read(control, second)["step_estimates"]
        assert all(seconds is not None and seconds > 0 for seconds in estimates)
    finally:
        page_jobs.forget_all_jobs()


def test_no_record_of_step_times_is_no_estimate(tmp_path):
    assert jobRunner.recorded_step_seconds(str(tmp_path)) == {}
    (tmp_path / jobRunner.STEP_TIMES_NAME).write_text("not json")
    assert jobRunner.recorded_step_seconds(str(tmp_path)) == {}
