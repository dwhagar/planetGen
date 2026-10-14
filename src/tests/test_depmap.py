# tests/test_depmap.py

"""The change-based test selection (TEST.139): the map, the selection rules, and a real recording."""

import json
import subprocess
import sys
import textwrap

import pytest

from tests import depmap


def _map(tests, files):
    return {"format": depmap.FORMAT, "files": files,
            "tests": {nodeid: {"f": [files.index(f) for f in used], "o": opaque, "s": 1.0}
                      for nodeid, (used, opaque) in tests.items()}}


FILES = ["src/planetgen/a.py", "src/planetgen/b.py", "src/tests/test_a.py", "src/tests/test_b.py", "src/html/x.html"]
DATA = _map({
    "src/tests/test_a.py::test_one": (["src/tests/test_a.py", "src/planetgen/a.py"], False),
    "src/tests/test_b.py::test_two": (["src/tests/test_b.py", "src/planetgen/b.py", "src/html/x.html"], False),
    "src/tests/test_b.py::test_cli": (["src/tests/test_b.py"], True),
}, FILES)
ALL = list(DATA["tests"])


def test_a_change_selects_the_tests_that_used_the_file():
    chosen, _ = depmap.select(DATA, {"src/planetgen/a.py"}, ALL)
    assert chosen == {"src/tests/test_a.py::test_one", "src/tests/test_b.py::test_cli"}  # the second is opaque


def test_a_template_a_test_opened_selects_it():
    chosen, _ = depmap.select(DATA, {"src/html/x.html"}, ALL)
    assert "src/tests/test_b.py::test_two" in chosen and "src/tests/test_a.py::test_one" not in chosen


def test_a_changed_test_file_selects_all_its_tests():
    chosen, _ = depmap.select(DATA, {"src/tests/test_b.py"}, ALL)
    assert chosen == {"src/tests/test_b.py::test_two", "src/tests/test_b.py::test_cli"}


def test_a_test_that_starts_a_subprocess_depends_on_the_program():
    chosen, _ = depmap.select(DATA, {"src/planetgen/unrelated.py"}, ALL)
    assert chosen == {"src/tests/test_b.py::test_cli"}
    assert depmap.select(DATA, {"docs/readme.md"}, ALL)[0] == set()


def test_a_test_the_map_has_not_seen_is_always_selected():
    chosen, _ = depmap.select(DATA, {"docs/readme.md"}, ALL + ["src/tests/test_new.py::test_x"])
    assert chosen == {"src/tests/test_new.py::test_x"}


def test_a_change_to_the_test_setup_or_no_map_selects_everything():
    for name in ("src/tests/conftest.py", "pytest.ini", "src/tests/fixtures/old_schemas/v60.sql", "src/tests/db_support.py"):
        assert depmap.select(DATA, {name}, ALL)[0] == set(ALL), name
    assert depmap.select(None, {"docs/readme.md"}, ALL)[0] == set(ALL)


def test_merging_keeps_unrun_tests_and_drops_deleted_files():
    old = _map({"src/tests/test_a.py::test_one": (["src/planetgen/a.py"], False),
                "src/tests/gone.py::test_old": (["src/planetgen/b.py"], False)}, FILES)
    merged = depmap.merge_map(old, {"src/tests/test_b.py::test_two": ({"src/planetgen/c.py"}, True, 2.5)},
                              lambda name: name != "src/tests/gone.py")
    assert set(merged["tests"]) == {"src/tests/test_a.py::test_one", "src/tests/test_b.py::test_two"}
    entry = merged["tests"]["src/tests/test_b.py::test_two"]
    assert entry["o"] is True and [merged["files"][i] for i in entry["f"]] == ["src/planetgen/c.py"]


def test_the_importers_of_a_changed_module_are_found_through_import_statements():
    found = depmap.importers_of({"src/planetgen/util/draw.py"})
    assert any(name.startswith("src/planetgen/") for name in found["src/planetgen/util/draw.py"])
    assert all(name.endswith(".py") for name in found["src/planetgen/util/draw.py"])


def test_a_changed_package_init_reaches_whatever_imports_inside_the_package():
    found = depmap.importers_of({"src/planetgen/names/__init__.py"})
    assert "src/tests/test_names.py" in found["src/planetgen/names/__init__.py"]


def test_a_recorded_run_lists_the_files_the_test_used(tmp_path):
    """Runs a three-test project (one fails) with --record-deps in a child pytest and reads the map it wrote: the
    failing test is not in it (a file outside the repository is not listed, so this proves that the plugin
    loads, records and writes, not what it lists)."""
    if not hasattr(sys, "monitoring"):
        pytest.skip("recording needs Python 3.12 or later")
    (tmp_path / "helper.py").write_text("def twice(x):\n    return 2 * x\n")
    (tmp_path / "data.txt").write_text("hi")
    (tmp_path / "test_sample.py").write_text(textwrap.dedent("""
        import helper
        def test_calls_the_helper():
            assert helper.twice(2) == 4
        def test_reads_a_file():
            assert open(__file__.replace("test_sample.py", "data.txt")).read() == "hi"
        def test_fails():
            assert helper.twice(2) == 5
        """))
    (tmp_path / "conftest.py").write_text(
        "import sys\nsys.path.insert(0, %r)\npytest_plugins = ('tests.depmap',)\n" % str(depmap.ROOT + "/src"))
    out = tmp_path / "map.json"
    subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", str(tmp_path), "--record-deps",
                    "--deps-file", str(out), "--rootdir", str(tmp_path)], cwd=tmp_path, check=False,
                   capture_output=True)
    data = json.loads(out.read_text())
    assert data["format"] == depmap.FORMAT
    assert sorted(name.split("::")[1] for name in data["tests"]) == ["test_calls_the_helper", "test_reads_a_file"]
