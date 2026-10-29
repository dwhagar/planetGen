# tests/test_test_runner.py

"""
Tests for the test suite's own organisation (TEST.139): the area and heavy
markers (`tests/areas.py`, added by `conftest.py`) and the slice picker
`scripts/test.py`.
"""

import importlib.util
import json
import os
import sys
import threading

import pytest

from tests import areas

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
TESTS_DIR = os.path.dirname(__file__)

_spec = importlib.util.spec_from_file_location("test_runner_script", os.path.join(REPO_ROOT, "scripts", "test.py"))
runner = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(runner)


@pytest.fixture(autouse=True)
def _progress_logs_in_a_throwaway_directory(tmp_path, monkeypatch):
    """`runner.main` starts a progress log (TEST.150); these tests must not leave logs in the checkout's own."""
    monkeypatch.setattr(runner, "PROGRESS_DIR", str(tmp_path / "progress"))
    # ... and must not clear the test Redis, which the tests running beside this one are using (TEST.152).
    monkeypatch.setattr(runner, "clear_test_redis", lambda: "")
    # ... and must not look at, or restart, the services the tests running beside this one stand on (TEST.376); the
    # tests of the check itself are in test_lane_services.py.
    monkeypatch.setattr(runner, "_services", lambda: _QuietServices)
    # ... nor end the queue workers of the tests beside this one (TEST.393).
    monkeypatch.setattr(runner, "reap_strays", lambda: "")


class _QuietWatchdog:
    """A watchdog that never trips."""

    def __init__(self, environ, on_trip, **options):
        self.events, self.what, self.failed, self.tripped = [], "", False, threading.Event()

    def start(self):
        pass

    def stop(self):
        pass


class _QuietServices:
    Watchdog = _QuietWatchdog

    @staticmethod
    def preflight(environ, browser=True, say=print):
        return []


@pytest.mark.parametrize("name, area", [
    ("test_web_pages.py", "web"),
    ("test_fuzz_web_routes.py", "web"),
    ("test_api_auth_sweep.py", "api"),
    ("test_bughunt_api_gaps.py", "api"),
    ("test_db_old_schemas.py", "db"),
    ("test_name_generation.py", "names"),
    ("test_orbital_motion.py", "physics"),
    ("test_galaxy_gen.py", "generation"),
    ("test_fuzz_cli.py", "ops"),
])
def test_a_test_file_gets_the_area_its_name_says(name, area):
    assert areas.area_of(name) == area


def test_a_name_no_rule_matches_is_other():
    assert areas.area_of("test_zzz_unheard_of.py") == "other"


def test_every_area_the_rules_give_is_one_the_runner_and_pytest_ini_know():
    known = set(runner.AREAS)
    assert {area for area, _pattern in areas.AREA_RULES} | {"other"} <= known
    with open(os.path.join(REPO_ROOT, "pytest.ini"), encoding="utf-8") as f:
        ini = f.read()
    for area in known:
        assert f"area_{area}:" in ini, f"area_{area} is not registered in pytest.ini"


def test_the_heavy_files_exist():
    """A renamed or deleted heavy file leaves a dead entry in `HEAVY_FILES`."""
    missing = sorted(name for name in areas.HEAVY_FILES if not os.path.exists(os.path.join(TESTS_DIR, name)))
    assert not missing, f"update HEAVY_FILES in tests/areas.py: {missing}"


def test_quick_leaves_out_database_browser_and_heavy_tests():
    [(label, marker, workers)] = runner.build_phases("quick", [], 4, 2)
    assert (label, workers) == ("quick", 4)
    assert marker == "not db and not browser and not heavy"


def test_full_runs_the_real_queue_tests_alone_and_the_browser_tests_last_on_fewer_workers():
    phases = runner.build_phases("full", [], 8, 2)
    assert phases == [("non-browser", "not browser and not real_queue", 8), ("real-queue", "real_queue and not browser", 1),
                      ("browser", "browser", 2)]


def test_budget_runs_only_the_budget_tests_on_one_worker_and_switches_them_on(monkeypatch):
    assert runner.build_phases("budget", [], 8, 2) == [("budget", "budget", 1)]
    seen = []
    monkeypatch.setattr(runner, "run_phase", lambda command, environment, **options: seen.append(environment) or (0, {"passed": 1}))
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    assert runner.main(["budget"]) == 0
    assert seen[0]["PLANETGEN_BUDGET"] == "1"


def test_a_run_ends_the_queue_workers_an_earlier_run_left_before_its_first_phase(monkeypatch, capsys):
    """TEST.398: a worker left by an earlier run picks up the units of this one, so it is ended before the first phase."""
    order = []
    monkeypatch.setattr(runner, "reap_strays", lambda: order.append("reap") or "2 queue worker process(es) a test started and did not stop")
    monkeypatch.setattr(runner, "run_phase", lambda command, environment, **options: order.append("phase") or (0, {"passed": 1}))
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    assert runner.main(["budget"]) == 0
    assert order[0] == "reap" and "phase" in order and order.index("reap") < order.index("phase")
    assert "before the run, ended 2 queue worker process(es)" in capsys.readouterr().out


def test_only_the_end_of_version_run_switches_on_the_version_gate(monkeypatch):
    seen = []
    monkeypatch.setattr(runner, "run_phase", lambda command, environment, **options: seen.append(environment) or (0, {"passed": 1}))
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    monkeypatch.delenv("PLANETGEN_VERSION_GATE", raising=False)
    assert runner.main(["full"]) == 0
    assert seen and all(environment["PLANETGEN_VERSION_GATE"] == "1" for environment in seen)
    seen.clear()
    assert runner.main(["quick"]) == 0
    assert seen and all("PLANETGEN_VERSION_GATE" not in environment for environment in seen)


def test_an_area_run_adds_the_browser_phase_only_for_the_web_area():
    assert [p[0] for p in runner.build_phases("area", ["db", "api"], 4, 2)] == ["non-browser", "real-queue"]
    web = runner.build_phases("area", ["web", "api"], 4, 2)
    assert web[0][1] == "(area_web or area_api) and not browser and not real_queue"
    assert web[1] == ("real-queue", "(area_web or area_api) and real_queue and not browser", 1)
    assert web[2] == ("browser", "(area_web or area_api) and browser", 2)


def test_pytest_command_uses_xdist_only_for_more_than_one_worker():
    assert "-n" not in runner.pytest_command("db", 1, [])
    command = runner.pytest_command("db", 4, ["-x"])
    assert command[command.index("-n") + 1] == "4" and command[-1] == "-x"
    assert command[command.index("-m", 3) + 1] == "db"


def test_a_run_is_remembered_for_the_next_one(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "RUNS_FILE", str(tmp_path / "runs.json"))
    assert runner.load_runs() == {}
    runner.save_run("full", 3312.4, True)
    assert runner.load_runs()["full"]["seconds"] == 3312
    assert runner.minutes(3312) == "55 min 12 s"


def test_everything_after_a_double_dash_goes_to_pytest(monkeypatch):
    calls = []
    monkeypatch.setattr(runner, "run_phase", lambda command, environment, **options: calls.append(command) or (0, {"passed": 1}))
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    assert runner.main(["quick", "--", "-x", "-k", "login", "--durations=5"]) == 0
    assert calls[0][-4:] == ["-x", "-k", "login", "--durations=5"]


def test_a_failing_phase_is_kept_with_its_full_text(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from tests.failure_capture import FailureCapture

    log = tmp_path / "failures.jsonl"
    capture = FailureCapture(str(log))
    failed = SimpleNamespace(outcome="failed", nodeid="src/tests/test_x.py::test_a", when="call", duration=1.234,
                             longrepr="Traceback ... AssertionError: boom")
    capture.pytest_runtest_logreport(failed)
    capture.pytest_runtest_logreport(SimpleNamespace(outcome="passed", nodeid="src/tests/test_x.py::test_b",
                                                     when="call", duration=0.1, longrepr=None))
    monkeypatch.setattr(runner, "FAILURES_FILE", str(log))
    found = runner.failed_tests(str(log))
    assert list(found) == ["src/tests/test_x.py::test_a"]
    record = found["src/tests/test_x.py::test_a"]
    assert record["when"] == "call" and record["seconds"] == 1.23 and "boom" in record["text"]


def test_a_missing_failure_log_means_no_failures(tmp_path):
    assert runner.failed_tests(str(tmp_path / "none.jsonl")) == {}


def test_changed_selects_and_full_refreshes_the_map():
    assert runner.depmap_options("changed", "", False) == ["--changed", ""]
    assert runner.depmap_options("changed", "abc123", False) == ["--changed", "abc123"]
    assert runner.depmap_options("full", "", False) == ["--record-deps"]
    assert runner.depmap_options("full", "abc123", True) == ["--record-deps", "--deps-check", "abc123"]
    assert runner.depmap_options("quick", "", False) == []


def test_a_check_with_no_change_to_compare_is_refused_and_says_what_base_to_pass(monkeypatch, capsys):
    """TEST.294: `full --check` on a checkout that is main compares it with itself, and called a failing test missed."""
    from tests import depmap

    ran = []
    monkeypatch.setattr(runner, "run_phase", lambda *args, **kwargs: ran.append(args) or (0, {}))
    monkeypatch.setattr(depmap, "changed_files", lambda base="": set())
    assert runner.main(["full", "--check"]) == 2
    said = capsys.readouterr().out
    assert "--base" in said and "Nothing was run" in said and not ran
    assert runner.check_without_changes("abc123").startswith("no file differs from abc123")

    monkeypatch.setattr(depmap, "changed_files", lambda base="": {"src/planetgen/generation/system.py"})
    assert runner.check_without_changes("") == ""
    monkeypatch.setattr(depmap, "changed_files", lambda base="": None)   # git cannot say: pytest reports that itself
    assert runner.check_without_changes("") == ""


def test_changed_runs_the_same_phases_as_full():
    assert runner.build_phases("changed", [], 4, 2) == runner.build_phases("full", [], 4, 2)


def test_changed_options_reach_pytest(monkeypatch):
    calls = []
    monkeypatch.setattr(runner, "run_phase", lambda command, environment, **options: calls.append(command) or (5, {}))
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    assert runner.main(["changed", "--base", "HEAD~1"]) == 0
    assert all(call[-2:] == ["--changed", "HEAD~1"] for call in calls) and len(calls) == 3


def test_the_summary_line_of_pytest_gives_the_counts():
    assert runner.summary_counts("=== 3 failed, 15419 passed, 22 skipped, 2 errors in 1081.44s (0:18:01) ===") == {
        "failed": 3, "passed": 15419, "skipped": 22, "errors": 2}
    assert runner.summary_counts("12 passed in 0.5s") == {}  # only the `=` banner line
    assert runner.summary_counts("FAILED src/tests/test_x.py::test_a - 3 passed") == {}


def test_the_last_line_says_pass_fail_or_interrupted(monkeypatch, capsys):
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    monkeypatch.setattr(runner, "rerun_failures", lambda environment: None)
    monkeypatch.setattr(runner, "run_phase", lambda command, environment, **options: (0, {"passed": 7, "skipped": 1}))
    assert runner.main(["quick"]) == 0
    assert "quick: PASS in" in capsys.readouterr().out.splitlines()[-1]
    monkeypatch.setattr(runner, "run_phase", lambda command, environment, **options: (1, {"passed": 5, "failed": 2}))
    assert runner.main(["quick"]) == 1
    last = capsys.readouterr().out.splitlines()[-1]
    assert "quick: FAIL in" in last and "5 passed, 2 failed" in last

    def stopped(command, environment, **options):
        raise KeyboardInterrupt("SIGTERM")

    monkeypatch.setattr(runner, "run_phase", stopped)
    assert runner.main(["quick"]) == 130
    last = capsys.readouterr().out.splitlines()[-1]
    assert "quick: INTERRUPTED (SIGTERM) in the quick phase" in last and "not finished" in last and "with --resume continues it" in last


def test_run_phase_passes_the_output_through_and_returns_the_counts(capfd):
    code, counts = runner.run_phase([sys.executable, "-c", "print('hello'); print('=== 2 passed, 1 skipped in 0.10s ===')"],
                                    dict(os.environ))
    assert code == 0 and counts == {"passed": 2, "skipped": 1}
    assert "hello" in capfd.readouterr().out


def test_timed_out_tests_are_named_from_the_failure_log(tmp_path):
    log = tmp_path / "failures.jsonl"
    log.write_text(json.dumps({"nodeid": "t.py::test_slow", "when": "call", "text": "Failed: Timeout (>900.0s) from pytest-timeout."}) + "\n"
                   + json.dumps({"nodeid": "t.py::test_wrong", "when": "call", "text": "AssertionError"}) + "\n")
    assert runner.timed_out_tests(str(log)) == ["t.py::test_slow"]
    assert runner.timed_out_tests(str(tmp_path / "missing.jsonl")) == []


def test_the_last_line_counts_the_tests_that_timed_out(monkeypatch, capsys):
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    monkeypatch.setattr(runner, "rerun_failures", lambda environment: None)
    monkeypatch.setattr(runner, "timed_out_tests", lambda: ["t.py::test_slow"])
    monkeypatch.setattr(runner, "run_phase", lambda command, environment, **options: (1, {"passed": 5, "failed": 1}))
    assert runner.main(["quick"]) == 1
    out = capsys.readouterr().out
    assert "TIMED OUT: 1 test(s)" in out and "  t.py::test_slow" in out
    assert out.splitlines()[-1].endswith("5 passed, 1 failed (1 TIMED OUT)")


def test_pytest_ini_gives_every_test_a_time_limit_that_stops_a_hang():
    import configparser
    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "pytest.ini"))
    assert 600 <= int(config["pytest"]["timeout"]) <= 1800  # the slowest test recorded took 479 s under load
    assert config["pytest"]["timeout_method"].strip() == "signal"  # stops the test in a worker without losing the worker


def test_a_test_that_hangs_is_stopped_and_named(tmp_path):
    import subprocess
    (tmp_path / "test_hang.py").write_text("import time\ndef test_hangs():\n    time.sleep(30)\n")
    result = subprocess.run([sys.executable, "-m", "pytest", "-p", "no:cacheprovider", "-c", os.devnull, "--rootdir", str(tmp_path),
                             "--timeout", "2", "--timeout-method", "signal", "-q", str(tmp_path / "test_hang.py")],
                            cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode == 1 and "Timeout (>2.0s) from pytest-timeout" in result.stdout


def test_the_worker_count_is_cut_to_the_memory_and_the_idle_cores():
    assert runner.cap_workers(8, 20000, 0.5, 8, 1000) == (8, "")
    assert runner.cap_workers(8, 3500, 0.5, 8, 1000) == (3, "3500 MB free at 1000 MB a worker")
    assert runner.cap_workers(8, 20000, 5.4, 8, 1000) == (3, "load 5.4 on 8 cores")
    assert runner.cap_workers(8, 20000, 12.0, 8, 1000) == (1, "load 12.0 on 8 cores")
    assert runner.cap_workers(2, 1200, 0.0, 4, 2500) == (1, "1200 MB free at 2500 MB a worker")
    assert runner.cap_workers(4, None, None, 8, 1000) == (4, "")  # no reading limits nothing
    assert runner.cap_workers(1, 100, 30.0, 8, 1000) == (1, "")  # never below one, and nothing to say about one worker


def test_a_phase_says_when_the_machine_cut_its_workers(monkeypatch, capsys):
    monkeypatch.setattr(runner, "available_memory_mb", lambda: 2200)
    monkeypatch.setattr(runner.os, "getloadavg", lambda: (0.0, 0.0, 0.0))
    seen = []
    monkeypatch.setattr(runner, "run_phase", lambda command, environment, **options: seen.append(command) or (0, {"passed": 1}))
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    assert runner.main(["full", "-n", "6"]) == 0
    out = capsys.readouterr().out
    assert "2 worker(s) (capped: 2200 MB free at 1000 MB a worker)" in out
    assert seen[0][seen[0].index("-n") + 1] == "2" and "-n" not in seen[2]  # browser: 2200 MB / 2500 -> one worker
    seen.clear()
    assert runner.main(["full", "-n", "6", "--no-worker-cap"]) == 0
    assert seen[0][seen[0].index("-n") + 1] == "6" and seen[2][seen[2].index("-n") + 1] == "2"
