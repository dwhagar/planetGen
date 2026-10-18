"""
The global documents (DOC.69, `scripts/build_global_docs.py`): one script rebuilds or checks the architecture map, the
plan pages and the docs index, in that order; the tests that compare each with its source are the `globaldocs`
tests, left out of a run unless PLANETGEN_GLOBAL_DOCS=1 (`scripts/test.py --global-docs`), so a stale global document
never fails a lane's gate.
"""

import importlib.util
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_spec = importlib.util.spec_from_file_location("build_global_docs", os.path.join(ROOT, "scripts", "build_global_docs.py"))
build = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(build)


class _Script:
    def __init__(self, calls, name, code):
        self.calls, self.name, self.code = calls, name, code

    def main(self, argv):
        self.calls.append((self.name, list(argv)))
        return self.code


def _load(calls, stale=()):
    return lambda name: _Script(calls, name, 1 if name in stale else 0)


def test_the_index_is_built_last_as_it_lists_the_other_pages():
    assert [name for _written, name in build.GENERATORS] == [
        "build_architecture_map", "build_todo_docs", "build_docs_index"]


def test_rebuilding_runs_every_generator_to_write():
    calls = []
    assert build.run(load=_load(calls)) == []
    assert calls == [(name, []) for _written, name in build.GENERATORS]


def test_checking_names_every_stale_document_and_runs_them_all():
    calls = []
    failed = build.run(check=True, load=_load(calls, stale={"build_architecture_map", "build_docs_index"}))
    assert failed == ["docs/design/architecture-map.html", "docs/index.html"]
    assert [name for name, _argv in calls] == [name for _written, name in build.GENERATORS]
    assert all(argv == ["--check"] for _name, argv in calls)


def test_the_command_exits_one_when_a_document_is_stale(monkeypatch, capsys):
    monkeypatch.setattr(build, "run", lambda check=False: ["docs/plan/*.html"] if check else [])
    assert build.main(["--check"]) == 1
    assert "docs/plan/*.html" in capsys.readouterr().out
    assert build.main([]) == 0


def test_the_freshness_tests_are_the_globaldocs_tests():
    marked = []
    for name in ("test_architecture_map.py", "test_docs_index.py", "test_todo_docs.py"):
        with open(os.path.join(ROOT, "src", "tests", name), encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        marked += [lines[i + 1].split("(")[0].removeprefix("def ") for i, line in enumerate(lines)
                   if line == "@pytest.mark.globaldocs"]
    assert {"test_pages_are_up_to_date", "test_the_index_is_up_to_date",
            "test_the_page_lists_every_module_of_the_package"} <= set(marked)
