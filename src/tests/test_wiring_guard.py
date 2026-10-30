"""
The wiring guard (GEN.754, scripts/wiring_guard.py): nothing built or downloaded stays unread without an open item that will
wire it, and the guard finds each kind of unread thing.
"""

import pathlib
import sys
import textwrap

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import wiring_guard as guard  # noqa: E402


@pytest.fixture(scope="module")
def code():
    return guard.Code()


def test_the_register_and_the_code_agree(code):
    todo = (ROOT / "docs" / "TODO.md").read_text(encoding="utf-8")
    assert guard.problems(guard.unread(code), guard.load_register(), todo) == []


def test_the_register_parses_every_line_it_has():
    entries = guard.load_register()
    assert entries and all(kind in guard.KINDS for kind, _ in entries)
    assert all(item == guard.BY_DESIGN or item.split(".")[0].isalpha() for item, _ in entries.values())


def test_the_guard_reaches_the_modules_the_app_reaches(code):
    assert "planetgen.physics.main_sequence" in code.reachable
    assert "planetgen.physics.line_opacity" not in code.reachable


def test_a_module_named_by_a_string_is_reached(code):
    assert "planetgen.generation.phenomena.builders" in code.reachable


TODO = "- [ ] **GEN.1 an open item**\n- [x] **GEN.2 a closed item**\n  text\n- [ ] **TEST.3b another**\n"


def test_the_open_items_are_the_unchecked_lines():
    assert guard.open_items(TODO) == {"GEN.1", "TEST.3b"}


def test_something_unread_and_not_in_the_register_is_named():
    found = {("module", "planetgen.physics.lonely")}
    assert any("planetgen.physics.lonely" in message and "nothing reads it" in message
               for message in guard.problems(found, {}, TODO))


def test_a_register_line_for_something_read_now_is_stale():
    entries = {("module", "planetgen.physics.wired"): ("GEN.1", "was waiting")}
    assert any("read now" in message for message in guard.problems(set(), entries, TODO))


@pytest.mark.parametrize("item,ok", [("GEN.1", True), ("GEN.1 and GEN.2", True), ("GEN.2", False), ("GEN.9", False),
                                     (guard.BY_DESIGN, True)])
def test_the_item_must_be_open_or_the_thing_by_design(item, ok):
    entries = {("data", "x.json"): (item, "a reason")}
    assert (guard.problems({("data", "x.json")}, entries, TODO) == []) is ok


def test_a_register_line_without_a_reason_is_refused():
    entries = {("data", "x.json"): ("GEN.1", "")}
    assert any("no reason" in message for message in guard.problems({("data", "x.json")}, entries, TODO))


def _write(root, relative, text=""):
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text), encoding="utf-8")


@pytest.fixture
def tiny(tmp_path, monkeypatch):
    """A small package: an app module that uses one physics module, a dataset with a reader nothing calls, a data file
    nothing names and a physics module nothing imports."""
    package = tmp_path / "src" / "planetgen"
    for folder in ("", "cli", "physics", "datasets"):
        _write(package, f"{folder}/__init__.py" if folder else "__init__.py")
    _write(package, "cli/run.py", """
        from planetgen.physics import used
        from planetgen.datasets import catalogue, fed
        value = used.VALUE + fed.size() + len(catalogue.ALL)
        """)
    _write(package, "physics/used.py", "VALUE = 1\nSTAND_IN = 2\nOTHER = 3\n")
    _write(package, "physics/lonely.py", "STAND_IN_OF_LONELY = 4\n")
    _write(package, "datasets/fetch.py", "def locate(dataset):\n    return dataset\n")
    _write(package, "datasets/fed.py", """
        from planetgen.datasets import fetch
        DATASET = 'fed'
        def size():
            return 1
        def read_table():
            return fetch.locate(DATASET)
        """)
    _write(package, "datasets/idle.py", "from planetgen.datasets import fetch\nDATASET = 'idle'\n")
    _write(package, "datasets/catalogue.py", "from planetgen.datasets import fed, fetch, idle\nALL = (fed.DATASET, idle.DATASET)\n")
    _write(package, "physics/data/named.json", "{}")
    _write(package, "physics/data/orphan.json", "{}")
    _write(package, "cli/reads.py", "NAME = 'named.json'\n")
    monkeypatch.setattr(guard, "SRC", str(tmp_path / "src"))
    monkeypatch.setattr(guard, "PACKAGE", str(package))
    monkeypatch.setattr(guard, "DATA_FOLDER", str(package / "physics" / "data"))
    return guard.Code()


def test_an_unimported_library_module_is_unread(tiny):
    assert tiny.unread_modules() == ["planetgen.physics.lonely"]


def test_a_dataset_nothing_imports_and_a_reader_nothing_calls_are_unread(tiny):
    datasets, readers = tiny.unread_datasets()
    assert datasets == ["idle"]
    assert readers == ["fed.read_table"]


def test_a_data_file_no_module_names_is_unread(tiny):
    assert tiny.unread_data_files() == ["orphan.json"]


def test_a_row_whose_constants_nothing_reads_is_unread(tiny):
    tables = {"6.1": (["fill_id", "constant"], [["A1", "used.VALUE and used.OTHER"], ["A2", "lonely.STAND_IN_OF_LONELY"],
                                                ["A3", "a law with no constant"], ["A4", "used.STAND_IN"]])}
    assert tiny.unread_rows(tables) == ["6.1:A2", "6.1:A4"]
