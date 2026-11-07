# tests/test_error_log.py

"""
OPS.116: one `error.log` in the log folder that records every error and crash, from every source, while each error still goes to its own log:
every logger at ERROR and above (the program's, another library's, Flask's), an uncaught exception in the main thread or another, a fatal
signal, a worker or a work horse that was killed or crashed (written by the parent), a web job's failed step, a one-off command that failed,
the 5xx answers of the web app, and `update.sh`, `install.sh` and the warm-up stages.
"""

import json
import logging
import os
import pathlib
import re
import signal
import subprocess
import sys
import threading
import time

import pytest

from planetgen.queue import api_jobs, redisqueue, units
from planetgen.util import joblog, log, logpaths, logrecord
from tests.test_api import client  # noqa: F401 -- the web app over a throwaway database

SRC = pathlib.Path(__file__).resolve().parents[1]
REPO = SRC.parent
LOGRECORD_SH = REPO / "scripts" / "logrecord.sh"
LINE = re.compile(r"^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d [A-Za-z0-9+\-:]+ (ERROR|CRITICAL)\s+\w[\w-]*: .* \(.*log=[\w.]+ pid=\d+\)$")


@pytest.fixture
def folder(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_LOG_DIR", str(tmp_path))
    log.configure(log.NORMAL, console=False)
    yield tmp_path
    log.configure(log.NORMAL, console=False)


def _lines(folder):
    path = folder / "error.log"
    return path.read_text(encoding="utf-8").splitlines() if path.exists() else []


def _python(code, folder, **env):
    environment = dict(os.environ, PLANETGEN_LOG_DIR=str(folder), PYTHONPATH=str(SRC), **env)
    return subprocess.run([sys.executable, "-c", code], env=environment, capture_output=True, text=True, timeout=60, cwd=str(folder))


# --- the logging layer: every logger at ERROR and above -----------------------------------------------------------------

def test_the_error_log_is_in_the_log_folder_and_the_folders_rotation_rule_covers_it(folder):
    assert logpaths.ERROR_LOG_NAME == "error.log" and "error.log" in logpaths.JOB_LOG_NAMES
    assert logpaths.error_log_path() == str(folder / "error.log")
    assert logpaths.job_log_path("error.log").endswith("error.log")


def test_an_error_from_the_programs_logger_is_in_the_error_log_with_the_log_it_also_went_to(folder):
    log.error("the sector could not be saved")
    log.get("warm").error("a family failed", family="stars")
    lines = _lines(folder)
    assert len(lines) == 2 and all(LINE.match(line) for line in lines)
    assert "generate: the sector could not be saved (log=console pid=" in lines[0]
    assert "warm: a family failed (family=stars log=warmup.log pid=" in lines[1]


@pytest.mark.parametrize("logger_name, component", [
    ("planetgen.queue.workers", "queue"),
    ("planetgen_app", "web"),                      # Flask's logger for the web app
    ("rq.worker", "lib"),
    ("urllib3.connectionpool", "lib"),
    ("some.other.library", "lib"),
])
def test_an_error_from_any_logger_is_in_the_error_log(folder, logger_name, component):
    logging.getLogger(logger_name).error("it failed")
    assert len(_lines(folder)) == 1 and f" {component}: it failed" in _lines(folder)[0]


def test_nothing_at_error_or_above_is_missing_and_nothing_below_is_there(folder):
    sources = [log.error, log.get("tiles").error, logging.getLogger("planetgen.web.x").error, logging.getLogger("flask.app").error,
               logging.getLogger().error, logging.getLogger("planetgen_app").error, log.get("orbits").critical,
               logging.getLogger("rq.worker").critical]
    for number, emit in enumerate(sources):
        emit(f"error number {number}")
    log.normal("information")
    log.warning("a warning")
    logging.getLogger("flask.app").warning("another warning")
    logging.getLogger("some.library").info("more information")
    lines = _lines(folder)
    assert [line.split(": ", 1)[1].split(" (")[0] for line in lines] == [f"error number {n}" for n in range(len(sources))]


def test_an_error_with_an_exception_carries_its_traceback(folder):
    try:
        raise ValueError("bad value")
    except ValueError:
        log.exception("could not read the file")
    text = "\n".join(_lines(folder))
    assert "could not read the file" in text and "Traceback (most recent call last)" in text and "ValueError: bad value" in text


def test_the_file_is_made_when_the_first_error_comes_with_the_mode_the_web_server_can_write(folder):
    assert not (folder / "error.log").exists()
    log.normal("nothing wrong")
    assert not (folder / "error.log").exists()
    log.error("something wrong")
    assert (folder / "error.log").stat().st_mode & 0o777 == 0o660


def test_a_folder_that_cannot_be_written_is_skipped_not_an_error(tmp_path, monkeypatch):
    blocker = tmp_path / "file"
    blocker.write_text("not a folder")
    monkeypatch.setenv("PLANETGEN_LOG_DIR", str(blocker / "logs"))
    log.configure(log.NORMAL, console=False)
    log.error("still logged where it goes")                      # no exception from the unwritable error log
    log.error("and again")


def test_the_program_rotates_the_file_itself_where_the_system_does_not(folder, monkeypatch):
    monkeypatch.setattr(logpaths, "log_rotation_mode", lambda config=None: "app")
    log.configure(log.NORMAL, console=False)
    log.error("rotated by the program")
    assert type(log._error_handler.inner).__name__ == "_RotatingHandler"


def test_an_error_is_still_written_to_the_log_it_would_have_gone_to(folder):
    with joblog.JobLog("tiles.log", "tiles", "Whole-map build"):
        log.get("tiles").error("a level could not be built")
    assert any("a level could not be built" in line for line in (folder / "tiles.log").read_text().splitlines())
    assert any("a level could not be built (log=tiles.log" in line for line in _lines(folder))


# --- a job's own log -----------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("raises", [SystemExit(3), RuntimeError("boom")])
def test_a_job_that_ends_badly_leaves_its_closing_line_in_the_error_log_too(folder, raises):
    with pytest.raises(type(raises)):
        with joblog.JobLog("orbits.log", "orbits", "Orbital update"):
            raise raises
    lines = _lines(folder)
    assert len(lines) == 1 and "orbits: Orbital update " in lines[0] and "log=orbits.log" in lines[0]


def test_a_job_that_finishes_leaves_nothing(folder):
    with joblog.JobLog("orbits.log", "orbits", "Orbital update"):
        log.normal("done")
    assert _lines(folder) == []


# --- an uncaught exception, a fatal signal -------------------------------------------------------------------------------

def test_an_uncaught_exception_in_the_main_thread_is_a_critical_line_with_its_traceback(tmp_path):
    result = _python("from planetgen.util import log\nlog.configure(console=False)\nraise RuntimeError('main crashed')\n", tmp_path)
    assert result.returncode == 1 and "RuntimeError: main crashed" in result.stderr        # the usual traceback still reaches stderr
    text = "\n".join(_lines(tmp_path))
    assert "CRITICAL" in text and "Uncaught exception" in text and "RuntimeError: main crashed" in text


def test_an_uncaught_exception_in_another_thread_is_recorded_with_the_thread(tmp_path):
    code = ("import threading\nfrom planetgen.util import log\nlog.configure(console=False)\n"
            "def work():\n    raise KeyError('thread crashed')\n"
            "thread = threading.Thread(target=work, name='the-helper')\nthread.start()\nthread.join()\n")
    result = _python(code, tmp_path)
    assert result.returncode == 0
    text = "\n".join(_lines(tmp_path))
    assert "CRITICAL" in text and "Uncaught exception in thread the-helper" in text and "KeyError: 'thread crashed'" in text


def test_a_ctrl_c_is_not_an_error(tmp_path):
    result = _python("from planetgen.util import log\nlog.configure(console=False)\nraise KeyboardInterrupt\n", tmp_path)
    assert result.returncode != 0 and _lines(tmp_path) == []


def test_a_fatal_signal_leaves_a_traceback_in_the_error_log(tmp_path):
    (tmp_path / "error.log").write_text("")                      # the setup scripts make the file; the dump goes only into one that is there
    code = ("import os, signal\nfrom planetgen.util import log\nlog.configure(console=False)\n"
            "os.kill(os.getpid(), signal.SIGSEGV)\ntime.sleep(5)\n")
    result = _python(code.replace("import os, signal", "import os, signal, time"), tmp_path)
    assert result.returncode == -signal.SIGSEGV
    assert "Fatal Python error: Segmentation fault" in (tmp_path / "error.log").read_text()


# --- a worker or a work horse that was killed or crashed, written by the parent ---------------------------------------

def _start(code, tmp_path):
    proc = redisqueue.start_detached([sys.executable, "-c", code, "--burst", "planetgen-test-queue"], str(tmp_path))
    deadline = time.monotonic() + 20
    thread = next((t for t in threading.enumerate() if t.name == f"planetgen-worker-{proc.pid}"), None)
    if thread is not None:                       # a reaper that has already ended has already written its line
        thread.join(max(0.1, deadline - time.monotonic()))
    return proc


def test_a_worker_killed_by_a_signal_is_logged_by_its_parent_with_the_signal(folder):
    proc = _start("import os, signal; os.kill(os.getpid(), signal.SIGKILL)", folder)
    lines = _lines(folder)
    assert len(lines) == 1 and f"Worker process {proc.pid} " in lines[0] and "was killed by SIGKILL" in lines[0] and "status=-9" in lines[0]


def test_a_worker_that_crashed_is_logged_with_its_exit_status(folder):
    _start("import sys; sys.exit(3)", folder)
    assert "exited with status 3" in "\n".join(_lines(folder))


@pytest.mark.parametrize("code", ["pass", "import sys; sys.exit(75)", "import os, signal; os.kill(os.getpid(), signal.SIGTERM)"])
def test_a_worker_that_ended_well_or_was_asked_to_stop_is_not_an_error(folder, code):
    _start(code, folder)
    assert _lines(folder) == []


class _Job:
    id = "job-1"
    func_name = "planetgen.queue.units.something_else"
    args = ()
    connection = None


def test_a_work_horse_that_was_killed_is_logged_with_the_signal(folder):
    units.horse_killed(_Job(), 4242, signal.SIGKILL, None)          # the raw wait status of a process killed by SIGKILL is the signal number
    lines = _lines(folder)
    assert len(lines) == 1 and "process running job job-1" in lines[0] and "pid 4242" in lines[0] and "was killed by SIGKILL" in lines[0]


def test_a_unit_whose_horse_was_killed_is_logged_and_fails_its_run(folder, monkeypatch):
    from planetgen.queue import runrecord, wrapper
    failed = []
    monkeypatch.setattr(runrecord, "finish_run", lambda conn, run, state, error=None: failed.append((run, state, error)))
    monkeypatch.setattr(wrapper, "stage_name", lambda conn, run, stage: "Scatter")
    job = _Job()
    job.func_name, job.args = units.UNIT_FUNCTION, ("run-9", 2, 7)
    units.horse_killed(job, 4243, signal.SIGKILL, None)
    assert failed and failed[0][1] == "failed"
    assert "unit 7 of run run-9, pid 4243, was killed by SIGKILL" in "\n".join(_lines(folder))


def test_a_work_horse_ended_by_sigterm_is_logged_because_rq_only_calls_the_handler_for_an_end_nobody_asked_for(folder):
    """OPS.118: a unit's horse ended with status 15 on Boss's test server; the next one names its signal, its unit and its run in the error log."""
    units.horse_killed(_Job(), 4244, signal.SIGTERM, None)
    assert "was killed by SIGTERM" in "\n".join(_lines(folder))


def test_a_decoded_wait_status_is_the_exit_status_or_minus_the_signal():
    assert log.decode_wait_status(signal.SIGKILL) == -9
    assert log.decode_wait_status(3 << 8) == 3
    assert log.decode_wait_status(0) == 0


# --- a one-off command, the lines a failing program printed ----------------------------------------------------------

def test_a_failing_command_is_logged_with_its_exit_status_and_the_error_lines_it_printed(folder):
    code = "import sys; print('working'); print('ValueError: no such sector', file=sys.stderr); sys.exit(2)"
    out = api_jobs.run_command([sys.executable, "-c", code], dict(os.environ), str(folder), 30)
    assert out["returncode"] == 2
    lines = _lines(folder)
    assert len(lines) == 1 and "exited with status 2; its last error lines: ValueError: no such sector" in lines[0]


def test_a_command_that_timed_out_or_could_not_start_is_logged(folder):
    api_jobs.run_command([sys.executable, "-c", "import time; time.sleep(30)"], dict(os.environ), str(folder), 1)
    api_jobs.run_command(["/nonexistent/planetgen-command"], dict(os.environ), str(folder), 5)
    text = "\n".join(_lines(folder))
    assert "timed out after 1 s" in text and "could not be started" in text


def test_a_command_that_succeeded_is_not_an_error(folder):
    api_jobs.run_command([sys.executable, "-c", "print('fine')"], dict(os.environ), str(folder), 30)
    assert _lines(folder) == []


def test_the_error_lines_of_a_program_are_picked_out_of_what_it_printed():
    text = "starting\nTraceback (most recent call last):\n  File \"x.py\", line 1\nValueError: bad\nwarning: fine\nERROR: lost it\n"
    assert log.error_lines(text) == ["Traceback (most recent call last):", "ValueError: bad", "ERROR: lost it"]
    assert log.error_lines(text, 1) == ["ERROR: lost it"] and log.error_lines("") == []


# --- a web job's failed step -----------------------------------------------------------------------------------------------

def test_a_web_job_step_that_failed_is_logged_with_its_error_lines(folder, mysql_config, monkeypatch):
    from planetgen.db import store
    from planetgen.web import job_runner
    store.get_control_connection(mysql_config, ensure_schema=True).close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    job_dir = folder / "jobs" / "20261001-120000-abcd"
    job_dir.mkdir(parents=True)
    failing = "import sys; print('ValueError: the sector is gone'); sys.exit(3)"
    (job_dir / "job.json").write_text(json.dumps({
        "id": job_dir.name, "kind": "galaxy", "title": "Generate sectors", "database": mysql_config.database,
        "created_at": time.time(), "cwd": str(folder), "env": {},
        "steps": [{"label": "Fill", "argv": [sys.executable, "-c", failing]}]}))
    assert job_runner.run(str(job_dir)) == 1
    lines = _lines(folder)
    assert len(lines) == 1 and "job: Step 'Fill' of web job 20261001-120000-abcd exited with status 3; its last error lines: " \
        "ValueError: the sector is gone" in lines[0]


# --- the web app: a 5xx answer ------------------------------------------------------------------------------------------

def _route(client, path):
    rule = next(rule for rule in client.application.url_map.iter_rules() if rule.rule == path)
    return rule.endpoint


def test_an_unhandled_error_in_a_route_is_logged_with_its_route_and_traceback(folder, client):
    app = client.application
    endpoint = _route(client, "/api/health")

    def broken():
        raise RuntimeError("route exploded")

    app.view_functions[endpoint] = broken
    app.config["PROPAGATE_EXCEPTIONS"] = False
    try:
        assert client.get("/api/health").status_code == 500
    finally:
        app.config["PROPAGATE_EXCEPTIONS"] = None
    text = "\n".join(_lines(folder))
    assert "GET /api/health" in text and f"route {endpoint}" in text and "RuntimeError: route exploded" in text
    assert len([line for line in _lines(folder) if "Unhandled exception" in line]) == 1       # once, not again for the 500 answer


def test_a_5xx_answer_without_an_exception_is_logged_with_its_route(folder, client):
    app = client.application
    endpoint = _route(client, "/api/health")
    app.view_functions[endpoint] = lambda: ("bad gateway", 502)
    assert client.get("/api/health").status_code == 502
    lines = _lines(folder)
    assert len(lines) == 1 and "web: 502 on GET /api/health" in lines[0] and f"route={endpoint}" in lines[0]


def test_the_deliberate_503_of_a_closed_site_or_a_rebuilding_cache_is_not_an_error(folder, client):
    from flask import g
    app = client.application
    endpoint = _route(client, "/api/health")

    def closed():
        g.expected_unavailable = True
        return ("closed for an update", 503)

    app.view_functions[endpoint] = closed
    assert client.get("/api/health").status_code == 503
    assert _lines(folder) == []


def test_a_4xx_answer_is_not_an_error_of_the_program(folder, client):
    assert client.get("/api/no-such-route").status_code == 404
    assert _lines(folder) == []


# --- the shell scripts ------------------------------------------------------------------------------------------------------

def _bash(script, tmp_path, **env):
    environment = dict(os.environ, **env)
    return subprocess.run(["bash", "-c", script], env=environment, capture_output=True, text=True, timeout=120, cwd=str(tmp_path))


def test_a_script_error_is_in_the_error_log_with_the_log_it_also_went_to_and_a_warning_is_not(tmp_path):
    target = tmp_path / "error.log"
    result = _bash(f'source "{LOGRECORD_SH}"; log_error_file_begin update.log; log_error "disk full"; log_warning "slow disk"; '
                   f'log_error_file_to "{target}"; log_error "again"; log_error_only "warm-up did not finish"', tmp_path)
    assert result.returncode == 0
    lines = target.read_text().splitlines()
    assert len(lines) == 3 and all(LINE.match(line) for line in lines)
    assert [line.split(": ", 1)[1].split(" (")[0] for line in lines] == ["disk full", "again", "warm-up did not finish"]
    assert all("(log=update.log pid=" in line for line in lines)
    assert "WARNING" in result.stderr and "slow disk" in result.stderr           # the screen shows the warning, the error log does not


def test_a_script_that_stops_on_a_failing_command_says_which_command_and_which_step(tmp_path):
    target = tmp_path / "error.log"
    script = (f'set -euo pipefail; source "{LOGRECORD_SH}"; log_error_file_begin console; log_error_file_to "{target}"; log_track_failures; '
              'trap \'log_exit_failure "The update" $?\' EXIT; log_stage 4 9 "Checking the math"; ls /nonexistent-planetgen-folder; echo not reached')
    result = _bash(script, tmp_path)
    assert result.returncode != 0 and "not reached" not in result.stdout
    text = target.read_text()
    assert "The update failed with exit status" in text and "during stage 4 of 9, checking the math" in text
    assert "The command was: ls /nonexistent-planetgen-folder" in text


def test_a_script_that_ends_well_leaves_nothing_in_the_error_log(tmp_path):
    target = tmp_path / "error.log"
    script = (f'set -euo pipefail; source "{LOGRECORD_SH}"; log_error_file_begin console; log_error_file_to "{target}"; log_track_failures; '
              'trap \'log_exit_failure "The update" $?\' EXIT; log_stage 1 9 "Pulling"; true')
    assert _bash(script, tmp_path).returncode == 0
    assert not target.exists() or target.read_text() == ""


def test_a_warm_up_stage_that_does_not_finish_is_in_the_error_log_while_the_screen_only_warns(tmp_path):
    fake = tmp_path / "fakepython"
    fake.write_text("#!/bin/bash\nexit 3\n")
    fake.chmod(0o755)
    target = tmp_path / "error.log"
    script = (f'SCRIPT_DIR="{tmp_path}"; source "{REPO}/scripts/deploy-common.sh"; PYTHON="{fake}"; log_error_file_begin update.log; '
              f'log_error_file_to "{target}"; run_warm_stage 2 2 "the full warm-up"; echo "returned $?"')
    env = {key: value for key, value in os.environ.items() if key != "WARM_STAGE_SECONDS"}
    result = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True, timeout=120)
    assert "returned 0" in result.stdout                                            # a failed warm-up never fails the update
    assert "WARNING" in result.stderr
    lines = target.read_text().splitlines()
    assert len(lines) == 1 and "the cache warm-up (part 2 of 2, the full warm-up) did not finish: it failed with status 3" in lines[0]


def test_update_sh_and_install_sh_start_the_error_log_and_record_what_stopped_them():
    update, install = (REPO / "update.sh").read_text(), (REPO / "install.sh").read_text()
    assert "log_error_file_begin update.log" in update and "log_track_failures" in update and "\nstart_error_log\n" in update
    assert "log_error_file_begin console" in install and "log_track_failures" in install and "\nstart_error_log\n" in install
    assert 'log_exit_failure "The install"' in install
    assert 'log_exit_failure "The update"' in (REPO / "scripts" / "deploy-common.sh").read_text()


# --- the component and the setup ------------------------------------------------------------------------------------------

def test_the_log_of_each_component_is_the_one_its_lines_also_go_to():
    assert logrecord.origin_log("warm") == "warmup.log" and logrecord.origin_log("tiles") == "tiles.log"
    assert logrecord.origin_log("orbits") == "orbits.log" and logrecord.origin_log("update") == "update.log"
    assert logrecord.origin_log("web") == "apache" and logrecord.origin_log("generate") == "console"
    assert set(logrecord.COMPONENT_LOGS) <= set(logrecord.COMPONENTS)


def test_the_setup_makes_the_file_for_the_web_servers_user_and_macos_rotates_it():
    setup = (REPO / "examples" / "apache" / "setup-debug-log.sh").read_text()
    assert "$ACTIVITY_DIR/error.log  $APACHE_USER:$APACHE_GROUP  660" in setup
    assert "error_log(" in (REPO / "examples" / "apache" / "log-locations.py").read_text()
