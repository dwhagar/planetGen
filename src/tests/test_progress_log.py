# tests/test_progress_log.py

"""TEST.150: the progress log of a test run and the id of its build (the first step of resuming a run, TEST.148)."""

import datetime
import json
import os
import subprocess
import sys
from types import SimpleNamespace

import pytest

from tests import buildid, progress_log

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SERVERS = lambda: {"mysql": "10.11.0-MariaDB", "redis": "7.0.0"}  # noqa: E731
NOW = datetime.datetime(2026, 10, 15, 12, 0, tzinfo=datetime.timezone.utc)


def _report(nodeid, when, outcome, duration=0.1, **extra):
    return SimpleNamespace(nodeid=nodeid, when=when, outcome=outcome, duration=duration, **extra)


def _feed(log, nodeid, setup="passed", call="passed", teardown="passed", **extra):
    for when, outcome in (("setup", setup), ("call", call), ("teardown", teardown)):
        if outcome is not None:
            log.pytest_runtest_logreport(_report(nodeid, when, outcome, **(extra if when == "call" else {})))


def test_a_line_is_written_when_a_test_finishes_not_before(tmp_path):
    path = tmp_path / "log.jsonl"
    log = progress_log.ProgressLog(str(path), "gw1")
    log.pytest_runtest_logreport(_report("t.py::test_a", "setup", "passed"))
    log.pytest_runtest_logreport(_report("t.py::test_a", "call", "passed", 0.25))
    assert not path.exists()  # killed here, the test has no line and runs again on a resume
    log.pytest_runtest_logreport(_report("t.py::test_a", "teardown", "passed", 0.05))
    [line] = [json.loads(text) for text in path.read_text().splitlines()]
    assert (line["kind"], line["nodeid"], line["outcome"], line["worker"]) == ("test", "t.py::test_a", "passed", "gw1")
    assert line["seconds"] == pytest.approx(0.4)


@pytest.mark.parametrize("phases, extra, outcome", [
    (dict(), {}, "passed"),
    (dict(call="failed"), {}, "failed"),
    (dict(teardown="failed"), {}, "failed"),
    (dict(setup="failed", call=None), {}, "failed"),
    (dict(setup="skipped", call=None), {}, "skipped"),
    (dict(call="skipped"), {"wasxfail": "known"}, "xfailed"),
    (dict(call="passed"), {"wasxfail": "known"}, "xpassed"),
    (dict(call="failed"), {"wasxfail": "known"}, "failed"),
])
def test_the_outcome_is_the_worst_of_the_phases(tmp_path, phases, extra, outcome):
    path = tmp_path / "log.jsonl"
    log = progress_log.ProgressLog(str(path))
    _feed(log, "t.py::test_a", **phases, **extra)
    assert json.loads(path.read_text())["outcome"] == outcome


def test_reading_a_log_gives_the_header_and_each_tests_last_line_and_skips_a_torn_line(tmp_path):
    path = tmp_path / "log.jsonl"
    path.write_text(json.dumps({"kind": "header", "build": "abc"}) + "\n"
                    + json.dumps({"kind": "test", "nodeid": "a", "outcome": "failed"}) + "\n"
                    + json.dumps({"kind": "test", "nodeid": "a", "outcome": "passed"}) + "\n"
                    + '{"kind": "test", "nodeid": "b", "outco')
    header, tests = progress_log.read(str(path))
    assert header["build"] == "abc"
    assert {nodeid: line["outcome"] for nodeid, line in tests.items()} == {"a": "passed"}
    assert progress_log.read(str(tmp_path / "missing.jsonl")) == ({}, {})


def test_a_real_run_with_workers_writes_one_line_per_test(tmp_path):
    (tmp_path / "test_sample.py").write_text(
        "import pytest\n"
        "def test_passes(): pass\n"
        "def test_fails(): assert False\n"
        "def test_skips(): pytest.skip('no')\n"
        "@pytest.mark.xfail\n"
        "def test_known(): assert False\n")
    log = tmp_path / "progress.jsonl"
    environment = dict(os.environ, PLANETGEN_PROGRESS_LOG=str(log), PYTHONPATH=SRC)
    subprocess.run([sys.executable, "-m", "pytest", "-p", "tests.progress_log", "-p", "no:cacheprovider", "-n", "2",
                    "-q", "--rootdir", str(tmp_path), "-c", os.devnull, str(tmp_path / "test_sample.py")],
                   cwd=tmp_path, env=environment, capture_output=True, text=True)
    _, tests = progress_log.read(str(log))
    assert {nodeid.split("::")[-1]: line["outcome"] for nodeid, line in tests.items()} == {
        "test_passes": "passed", "test_fails": "failed", "test_skips": "skipped", "test_known": "xfailed"}
    assert {line["worker"] for line in tests.values()} <= {"gw0", "gw1"}


# --- the build id ------------------------------------------------------------------------------------------------

class Repo(os.PathLike):
    """A throwaway git repository: a path that also runs git in itself."""

    def __init__(self, path):
        self.path = path

    def __fspath__(self):
        return str(self.path)

    __str__ = __fspath__

    def __truediv__(self, name):
        return self.path / name

    def git(self, *args):
        subprocess.run(["git", "-c", "user.email=t@example.com", "-c", "user.name=t", *args], cwd=self.path, check=True,
                       capture_output=True)


@pytest.fixture
def repo(tmp_path):
    def git(*args):
        subprocess.run(["git", "-c", "user.email=t@example.com", "-c", "user.name=t", *args], cwd=tmp_path, check=True,
                       capture_output=True)

    git("init", "-q", "-b", "main")
    (tmp_path / "code.py").write_text("x = 1\n")
    (tmp_path / ".gitignore").write_text("ignored.txt\n")
    git("add", "-A")
    git("commit", "-q", "-m", "first")
    return Repo(tmp_path)


def _id(repo, command="full", args=(), base="", environ=None, now=NOW):
    parts = buildid.build_parts(str(repo), command, list(args), base, environ=environ or {}, now=now, servers=SERVERS)
    return buildid.build_id(parts), parts


def test_the_same_inputs_give_the_same_id(repo):
    assert _id(repo)[0] == _id(repo)[0]
    assert len(_id(repo)[0]) == 16 and _id(repo)[1]["dirty"] == "clean"


def test_a_change_to_a_tracked_file_changes_the_id_and_putting_it_back_restores_it(repo):
    clean, _ = _id(repo)
    (repo / "code.py").write_text("x = 2\n")
    edited, parts = _id(repo)
    assert edited != clean and parts["dirty"] != "clean"
    (repo / "code.py").write_text("x = 1\n")
    assert _id(repo)[0] == clean


def test_a_new_file_changes_the_id_unless_git_ignores_it(repo):
    clean, _ = _id(repo)
    (repo / "ignored.txt").write_text("noise")
    assert _id(repo)[0] == clean
    (repo / "new.py").write_text("y = 1\n")
    added, _ = _id(repo)
    assert added != clean
    (repo / "new.py").write_text("y = 2\n")
    assert _id(repo)[0] != added


def test_a_commit_of_the_same_content_keeps_the_tree_and_a_different_content_changes_it(repo):
    before = _id(repo)
    repo.git("commit", "-q", "--allow-empty", "-m", "empty")
    assert _id(repo)[0] == before[0]  # the same files, so the same build
    (repo / "code.py").write_text("x = 3\n")
    repo.git("commit", "-q", "-am", "edit")
    assert _id(repo)[1]["tree"] != before[1]["tree"]


@pytest.mark.parametrize("change", [
    dict(command="changed"),
    dict(args=["-x"]),
    dict(environ={"PLANETGEN_TEST_SOMETHING": "1"}),
    dict(now=NOW + datetime.timedelta(days=1)),
])
def test_the_command_its_arguments_the_environment_and_the_day_change_the_id(repo, change):
    assert _id(repo, **change)[0] != _id(repo)[0]


def test_the_log_paths_and_the_worker_name_do_not_change_the_id(repo):
    base = _id(repo)[0]
    assert _id(repo, environ={"PLANETGEN_FAILURE_LOG": "/x", "PLANETGEN_PROGRESS_LOG": "/y", "PYTEST_XDIST_WORKER": "gw3",
                              "HOME": "/elsewhere"})[0] == base


def test_the_servers_versions_change_the_id(repo):
    base = _id(repo)[0]
    parts = buildid.build_parts(str(repo), "full", [], environ={}, now=NOW,
                                servers=lambda: {"mysql": "8.4.0", "redis": "7.0.0"})
    assert buildid.build_id(parts) != base
    assert buildid.differences(_id(repo)[1], parts) == ["mysql"]


def test_no_part_shows_a_value_from_the_environment(repo):
    secret = "hunter2-secret"
    parts = _id(repo, environ={"PLANETGEN_TEST_MYSQL_PASSWORD": secret})[1]
    assert secret not in json.dumps(parts)


# --- the runner -------------------------------------------------------------------------------------------------

def test_the_runner_starts_a_log_for_the_build_and_hands_its_path_to_pytest(tmp_path, monkeypatch):
    import importlib.util
    spec = importlib.util.spec_from_file_location("planetgen_test_runner", os.path.join(os.path.dirname(SRC), "scripts", "test.py"))
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    monkeypatch.setattr(runner, "PROGRESS_DIR", str(tmp_path / "progress"))
    seen = []
    monkeypatch.setattr(runner, "run_phase", lambda command, environment: seen.append(environment) or (0, {"passed": 1}))
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    old = tmp_path / "progress"
    old.mkdir()
    (old / "stale.jsonl").write_text("{}\n")
    os.utime(old / "stale.jsonl", (0, 0))
    (old / "recent.jsonl").write_text("{}\n")
    assert runner.main(["quick"]) == 0
    path = seen[0]["PLANETGEN_PROGRESS_LOG"]
    assert os.path.dirname(path) == str(old) and path.endswith(".jsonl")
    header, tests = progress_log.read(path)
    assert header["build"] == os.path.basename(path)[:-len(".jsonl")] and header["command"] == "quick"
    assert set(header["parts"]) == {"tree", "dirty", "packages", "python", "mysql", "redis", "command", "selection", "env", "date"}
    assert not (old / "stale.jsonl").exists() and (old / "recent.jsonl").exists()
    # a new run of the same build starts the log again
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps({"kind": "test", "nodeid": "a", "outcome": "passed"}) + "\n")
    assert runner.main(["quick"]) == 0
    assert progress_log.read(path)[1] == {}
