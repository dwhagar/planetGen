# tests/test_resume_across_builds.py

"""TEST.208: `scripts/test.py --resume` reuses a pass from an earlier run of another build when nothing the test
depends on has changed since, and runs the test when a file it used, a package, a server or the day has."""

import datetime
import importlib.util
import json
import os
from types import SimpleNamespace

import pytest

from tests import buildid, progress_log, reuse

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(SRC)
NOW = datetime.datetime(2026, 10, 15, 12, 0, tzinfo=datetime.timezone.utc)
PARTS = {"tree": "t1", "dirty": "clean", "packages": "p1", "python": "3.13.0", "mysql": "8.0", "redis": "7.0",
         "command": "c1", "selection": "s1", "env": "e1", "date": "2026-10-15"}
FILES = {"src/planetgen/a.py": "1", "src/planetgen/b.py": "1", "src/planetgen/templates/page.html": "1",
         "src/tests/test_a.py": "1", "src/tests/test_b.py": "1", "src/tests/conftest.py": "1", "docs/notes.md": "1"}


def _log(directory, name, parts=None, manifest=None, tests=None, started="2026-10-15T10:00:00Z", command="full"):
    """An earlier run's progress log: `tests` is `{node id: (files used, opaque)}`, all of them passed."""
    path = os.path.join(directory, name + ".jsonl")
    lines = [{"kind": "header", "build": name, "command": command, "parts": parts or PARTS, "started": started,
              "manifest": manifest if manifest is not None else FILES}]
    for nodeid, (files, opaque) in (tests or {}).items():
        lines.append({"kind": "test", "nodeid": nodeid, "outcome": "passed"})
        lines.append({"kind": "deps", "nodeid": nodeid, "files": files, "opaque": opaque})
    with open(path, "w", encoding="utf-8") as f:
        f.write("".join(json.dumps(line) + "\n" for line in lines))
    return path


TESTS = {
    "src/tests/test_a.py::test_uses_a": (["src/tests/test_a.py", "src/planetgen/a.py"], False),
    "src/tests/test_a.py::test_uses_page": (["src/tests/test_a.py", "src/planetgen/templates/page.html"], False),
    "src/tests/test_b.py::test_uses_b": (["src/tests/test_b.py", "src/planetgen/b.py"], False),
    "src/tests/test_b.py::test_runs_a_subprocess": (["src/tests/test_b.py"], True),
}


@pytest.fixture
def directory(tmp_path):
    path = tmp_path / "progress"
    path.mkdir()
    return str(path)


def _reusable(directory, now_manifest, parts=PARTS):
    found, looked, blocked = reuse.reusable(directory, parts, now_manifest)
    return set(found), blocked


def test_an_unchanged_tree_reuses_every_pass(directory):
    _log(directory, "old", tests=TESTS)
    assert _reusable(directory, dict(FILES))[0] == set(TESTS)


def test_the_tree_may_differ_from_the_earlier_runs_where_no_test_looked(directory):
    _log(directory, "old", tests=TESTS)
    reused, blocked = _reusable(directory, dict(FILES, **{"docs/notes.md": "2", "docs/new.md": "1"}))
    assert reused == set(TESTS) and not blocked


def test_a_changed_file_a_test_used_reruns_that_test_only(directory):
    _log(directory, "old", tests=TESTS)
    opaque = "src/tests/test_b.py::test_runs_a_subprocess"      # reruns for any change under src/planetgen/
    reused, _ = _reusable(directory, dict(FILES, **{"src/planetgen/a.py": "2"}))
    assert reused == set(TESTS) - {"src/tests/test_a.py::test_uses_a", opaque}
    reused, _ = _reusable(directory, dict(FILES, **{"src/planetgen/templates/page.html": "2"}))
    assert reused == set(TESTS) - {"src/tests/test_a.py::test_uses_page", opaque}


def test_a_changed_test_file_reruns_every_test_in_it(directory):
    _log(directory, "old", tests=TESTS)
    reused, _ = _reusable(directory, dict(FILES, **{"src/tests/test_b.py": "2"}))
    assert {nodeid.split("::")[0] for nodeid in reused} == {"src/tests/test_a.py"}


def test_a_test_that_starts_a_subprocess_reruns_for_a_change_to_the_code_but_not_to_the_docs(directory):
    _log(directory, "old", tests=TESTS)
    opaque = "src/tests/test_b.py::test_runs_a_subprocess"
    assert opaque in _reusable(directory, dict(FILES, **{"docs/notes.md": "2"}))[0]
    assert opaque not in _reusable(directory, dict(FILES, **{"src/planetgen/b.py": "2"}))[0]


def test_a_change_to_the_test_setup_reruns_everything(directory):
    _log(directory, "old", tests=TESTS)
    assert _reusable(directory, dict(FILES, **{"src/tests/conftest.py": "2"}))[0] == set()


def test_a_new_package_server_or_day_reuses_nothing_and_names_the_part(directory):
    _log(directory, "old", tests=TESTS)
    for part, value in (("packages", "p2"), ("python", "3.14.0"), ("mysql", "8.4"), ("redis", "8.0"), ("env", "e2"),
                        ("date", "2026-10-16")):
        reused, blocked = _reusable(directory, dict(FILES), dict(PARTS, **{part: value}))
        assert reused == set() and blocked == [part]


def test_the_tree_the_command_and_the_selection_do_not_keep_a_pass_out(directory):
    _log(directory, "old", tests=TESTS)
    parts = dict(PARTS, tree="t2", dirty="d2", command="c2", selection="s2")
    assert _reusable(directory, dict(FILES), parts)[0] == set(TESTS)


def test_the_newest_run_is_the_one_a_pass_is_taken_from(directory):
    _log(directory, "older", tests={"src/tests/test_a.py::test_uses_a": (["src/tests/test_a.py"], False)},
         started="2026-10-15T09:00:00Z")
    _log(directory, "newer", tests={"src/tests/test_a.py::test_uses_a": (["src/tests/test_a.py"], False)},
         started="2026-10-15T11:00:00Z")
    found, _looked, _blocked = reuse.reusable(directory, PARTS, dict(FILES))
    assert found["src/tests/test_a.py::test_uses_a"]["build"] == "newer"


def test_a_pass_without_a_deps_line_or_a_log_without_a_manifest_is_never_reused(directory):
    path = _log(directory, "old", tests=TESTS)
    lines = [json.loads(text) for text in open(path)]
    with open(path, "w") as f:   # the deps line of one test never got written (the run was killed between the two)
        f.write("".join(json.dumps(line) + "\n" for line in lines
                        if not (line["kind"] == "deps" and line["nodeid"].endswith("test_uses_b"))))
    _log(directory, "ancient", manifest=None, tests=TESTS)
    os.remove(os.path.join(directory, "ancient.jsonl"))
    with open(os.path.join(directory, "legacy.jsonl"), "w") as f:
        f.write(json.dumps({"kind": "header", "build": "legacy", "parts": PARTS}) + "\n")
        f.write(json.dumps({"kind": "test", "nodeid": "src/tests/test_a.py::test_legacy", "outcome": "passed"}) + "\n")
    reused = _reusable(directory, dict(FILES))[0]
    assert "src/tests/test_b.py::test_uses_b" not in reused and "src/tests/test_a.py::test_legacy" not in reused
    assert len(reused) == len(TESTS) - 1


def test_the_manifest_lists_tracked_and_untracked_files_and_notices_an_edit(tmp_path):
    import subprocess
    run = lambda *args: subprocess.run(["git", *args], cwd=tmp_path, check=True, capture_output=True)  # noqa: E731
    run("init", "-q")
    (tmp_path / "a.py").write_text("one")
    (tmp_path / "b.py").write_text("two")
    run("add", "a.py", "b.py")
    run("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "x")
    clean = reuse.manifest(str(tmp_path))
    assert set(clean) == {"a.py", "b.py"}
    (tmp_path / "a.py").write_text("three")
    (tmp_path / "c.py").write_text("new")
    edited = reuse.manifest(str(tmp_path))
    assert reuse.changed_between(clean, edited) == {"a.py", "c.py"}
    (tmp_path / "a.py").write_text("one")
    (tmp_path / "c.py").unlink()
    assert reuse.changed_between(clean, reuse.manifest(str(tmp_path))) == set()   # back to its committed text


def test_a_passed_test_writes_a_deps_line_and_a_failed_one_does_not(tmp_path):
    path = tmp_path / "log.jsonl"
    log = progress_log.ProgressLog(str(path))
    config = SimpleNamespace(pluginmanager=SimpleNamespace(get_plugin=lambda name: log))
    progress_log.record_deps(config, "t.py::test_a", {"src/b.py", "src/a.py"}, False)
    [line] = [json.loads(text) for text in path.read_text().splitlines()]
    assert line == {"kind": "deps", "nodeid": "t.py::test_a", "files": ["src/a.py", "src/b.py"], "opaque": False}
    assert progress_log.read(str(path)) == ({}, {})   # a deps line is neither the header nor a finished test
    progress_log.record_deps(SimpleNamespace(pluginmanager=SimpleNamespace(get_plugin=lambda name: None)), "x", [], False)


# --- the runner ---------------------------------------------------------------------------------------------------

def _runner(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("planetgen_test_runner_reuse", os.path.join(ROOT, "scripts", "test.py"))
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    monkeypatch.setattr(runner, "PROGRESS_DIR", str(tmp_path / "progress"))
    monkeypatch.setattr(runner, "clear_test_redis", lambda: "")
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    monkeypatch.setattr(runner, "rerun_failures", lambda environment: None)
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    return runner


def _earlier_run_of_another_build(runner, directory, changed_file=None):
    """A log of `quick` as this checkout's environment would write it, from a tree one file apart."""
    parts = buildid.build_parts(runner.ROOT, "quick", [], "", now=NOW)
    manifest = reuse.manifest(runner.ROOT)
    some_file = "src/planetgen/_version.py"
    manifest[some_file] = "other"                      # the earlier tree differed in this file
    os.makedirs(directory, exist_ok=True)
    tests = {"src/tests/test_x.py::test_one": ([f"src/tests/test_x.py", changed_file or "src/planetgen/_other.py"], False),
             "src/tests/test_x.py::test_two": (["src/tests/test_x.py", some_file], False)}
    return _log(directory, "earlierbuild", parts=dict(parts, tree="earlier", command="old"), manifest=manifest, tests=tests,
                started="2026-10-15T09:30:00Z", command="quick")


def test_a_resume_after_a_merge_skips_the_tests_whose_files_did_not_change_and_carries_them_on(tmp_path, monkeypatch):
    runner = _runner(tmp_path, monkeypatch)
    directory = str(tmp_path / "progress")
    _earlier_run_of_another_build(runner, directory)
    resumed = runner.start_progress_log("quick", [], "", directory=directory, now=NOW, resume=True)
    assert resumed.problem is None and resumed.passed == 1
    assert json.load(open(resumed.skip_file)) == ["src/tests/test_x.py::test_one"]    # test_two used the changed file
    assert resumed.reused == [("earlierbuild", "2026-10-15T09:30:00Z", 1)]
    header, tests = progress_log.read(resumed.path)
    assert header["build"] == resumed.identifier and set(tests) == {"src/tests/test_x.py::test_one"}
    assert tests["src/tests/test_x.py::test_one"]["reused_from"] == "earlierbuild"
    assert "src/tests/test_x.py::test_one" in reuse.read_deps(resumed.path)
    # The carried pass counts for the next resume of this build too.
    again = runner.start_progress_log("quick", [], "", directory=directory, now=NOW, resume=True)
    assert again.identifier == resumed.identifier and json.load(open(again.skip_file)) == ["src/tests/test_x.py::test_one"]


def test_a_resume_in_another_environment_runs_every_test_and_says_why(tmp_path, monkeypatch):
    runner = _runner(tmp_path, monkeypatch)
    directory = str(tmp_path / "progress")
    _earlier_run_of_another_build(runner, directory)
    tomorrow = runner.start_progress_log("quick", [], "", directory=directory, now=NOW + datetime.timedelta(days=1), resume=True)
    assert tomorrow.skip_file is None and tomorrow.passed == 0 and "date" in tomorrow.note


def test_main_says_how_many_tests_were_reused_and_from_which_run(tmp_path, monkeypatch, capsys):
    runner = _runner(tmp_path, monkeypatch)
    commands = []
    monkeypatch.setattr(runner, "run_phase", lambda command, environment: commands.append(command) or (0, {"passed": 3}))
    monkeypatch.setattr(runner.datetime, "datetime", type("D", (datetime.datetime,), {
        "now": classmethod(lambda cls, tz=None: NOW)}))
    _earlier_run_of_another_build(runner, str(tmp_path / "progress"))
    assert runner.main(["quick", "--resume"]) == 0
    out = capsys.readouterr().out.splitlines()
    assert any("RESUMING build" in line and "the 1 tests that passed earlier are skipped" in line
               and "1 in the run started 2026-10-15T09:30:00Z of build earlierbuild" in line for line in out)
    assert "--skip-passed" in commands[0]
    assert "(RESUMED build" in out[-1] and "1 more passed earlier" in out[-1] and "of build earlierbuild" in out[-1]
