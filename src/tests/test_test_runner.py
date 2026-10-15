# tests/test_test_runner.py

"""
Tests for the test suite's own organisation (TEST.139): the area and heavy
markers (`tests/areas.py`, added by `conftest.py`) and the slice picker
`scripts/test.py`.
"""

import importlib.util
import os
import sys

import pytest

from tests import areas

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
TESTS_DIR = os.path.dirname(__file__)

_spec = importlib.util.spec_from_file_location("test_runner_script", os.path.join(REPO_ROOT, "scripts", "test.py"))
runner = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(runner)


@pytest.mark.parametrize("name, area", [
    ("test_web_pages.py", "web"),
    ("test_fuzz_web_routes.py", "web"),
    ("test_api_auth_sweep.py", "api"),
    ("test_bughunt_api_gaps.py", "api"),
    ("test_db_old_schemas.py", "db"),
    ("test_wikiclient_wikijs.py", "wiki"),
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


def test_full_runs_the_browser_tests_after_the_rest_on_fewer_workers():
    phases = runner.build_phases("full", [], 8, 2)
    assert phases == [("non-browser", "not browser", 8), ("browser", "browser", 2)]


def test_budget_runs_only_the_budget_tests_on_one_worker_and_switches_them_on(monkeypatch):
    assert runner.build_phases("budget", [], 8, 2) == [("budget", "budget", 1)]
    seen = []
    monkeypatch.setattr(runner, "run_phase", lambda command, environment: seen.append(environment) or (0, {"passed": 1}))
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    assert runner.main(["budget"]) == 0
    assert seen[0]["PLANETGEN_BUDGET"] == "1"


def test_an_area_run_adds_the_browser_phase_only_for_the_web_area():
    assert [p[0] for p in runner.build_phases("area", ["db", "api"], 4, 2)] == ["non-browser"]
    web = runner.build_phases("area", ["web", "api"], 4, 2)
    assert web[0][1] == "(area_web or area_api) and not browser"
    assert web[1] == ("browser", "(area_web or area_api) and browser", 2)


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
    monkeypatch.setattr(runner, "run_phase", lambda command, environment: calls.append(command) or (0, {"passed": 1}))
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


def test_changed_runs_the_same_two_phases_as_full():
    assert runner.build_phases("changed", [], 4, 2) == runner.build_phases("full", [], 4, 2)


def test_changed_options_reach_pytest(monkeypatch):
    calls = []
    monkeypatch.setattr(runner, "run_phase", lambda command, environment: calls.append(command) or (5, {}))
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    assert runner.main(["changed", "--base", "HEAD~1"]) == 0
    assert all(call[-2:] == ["--changed", "HEAD~1"] for call in calls) and len(calls) == 2


def test_the_summary_line_of_pytest_gives_the_counts():
    assert runner.summary_counts("=== 3 failed, 15419 passed, 22 skipped, 2 errors in 1081.44s (0:18:01) ===") == {
        "failed": 3, "passed": 15419, "skipped": 22, "errors": 2}
    assert runner.summary_counts("12 passed in 0.5s") == {}  # only the `=` banner line
    assert runner.summary_counts("FAILED src/tests/test_x.py::test_a - 3 passed") == {}


def test_the_last_line_says_pass_fail_or_interrupted(monkeypatch, capsys):
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    monkeypatch.setattr(runner, "rerun_failures", lambda environment: None)
    monkeypatch.setattr(runner, "run_phase", lambda command, environment: (0, {"passed": 7, "skipped": 1}))
    assert runner.main(["quick"]) == 0
    assert "quick: PASS in" in capsys.readouterr().out.splitlines()[-1]
    monkeypatch.setattr(runner, "run_phase", lambda command, environment: (1, {"passed": 5, "failed": 2}))
    assert runner.main(["quick"]) == 1
    last = capsys.readouterr().out.splitlines()[-1]
    assert "quick: FAIL in" in last and "5 passed, 2 failed" in last

    def stopped(command, environment):
        raise KeyboardInterrupt("SIGTERM")

    monkeypatch.setattr(runner, "run_phase", stopped)
    assert runner.main(["quick"]) == 130
    last = capsys.readouterr().out.splitlines()[-1]
    assert "quick: INTERRUPTED (SIGTERM) in the quick phase" in last and "not finished" in last


def test_run_phase_passes_the_output_through_and_returns_the_counts(capfd):
    code, counts = runner.run_phase([sys.executable, "-c", "print('hello'); print('=== 2 passed, 1 skipped in 0.10s ===')"],
                                    dict(os.environ))
    assert code == 0 and counts == {"passed": 2, "skipped": 1}
    assert "hello" in capfd.readouterr().out
