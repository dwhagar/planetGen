# tests/test_allocate_numbers.py

"""
OPS.71: `scripts/allocate_numbers.py` gives a branch's epochs and migrations the next free numbers after main's, in one pass, with no hand
edits: two branches that took the same epoch and migration both merge in turn. Each test builds a small git repository with main and a
branch that took the numbers main has since taken.
"""

import gzip
import importlib.util
import json
import os
import subprocess

import pytest

from planetgen.util import generator_epoch

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "allocate_numbers.py")
VERSIONS = "src/planetgen/db/migrations/versions"
EPOCHS = "src/tests/golden/epochs.json"
FIXTURES = "src/tests/fixtures/old_schemas"
SCHEMA = "src/planetgen/db/schema.sql"


@pytest.fixture(scope="module")
def tool():
    spec = importlib.util.spec_from_file_location("allocate_numbers", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _git(root, *args, check=True):
    done = subprocess.run(["git", *args], cwd=root, capture_output=True, text=True)
    if check:
        assert done.returncode == 0, done.stderr
    return done.stdout


def _write(root, path, text):
    full = os.path.join(root, path)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w", encoding="utf-8") as handle:
        handle.write(text)


def _migration(number, name, previous, doc):
    return (f'"""{doc}\n\nSchema v{number}. Epoch lines in the docs say {number}."""\n\nrevision = "{number:04d}"\n'
            f'down_revision = "{previous:04d}"\n\n\ndef upgrade():\n    pass\n')


def _epochs(root, entries):
    data = {"stages": ["plan", "scatter", "fill", "population"], "stage_note": "note", "epochs": [
        {"epoch": index + 1, "first_release": None, "reason": reason, "affects": ["fill"]} for index, reason in enumerate(entries)]}
    _write(root, EPOCHS, generator_epoch.format_epochs(data))


def _schema(root, columns):
    _write(root, SCHEMA, "-- a comment\nCREATE TABLE t (\n    id INT" + "".join(f",\n    {c} INT" for c in columns) + "\n);\n\n-- more\n")


def _fixture(root, version, text):
    full = os.path.join(root, FIXTURES, f"schema_v{version}.sql.gz")
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with gzip.open(full, "wt", encoding="utf-8") as handle:
        handle.write(text)


@pytest.fixture
def repo(tmp_path):
    """main with migration 0001, epochs 1-2; a branch that adds 0002 and epoch 3; main meanwhile adds its own 0002 and epoch 3."""
    root = str(tmp_path)
    _git(root, "init", "-q", "-b", "main")
    _git(root, "config", "user.email", "t@example.org")
    _git(root, "config", "user.name", "t")
    _write(root, f"{VERSIONS}/0001_first.py", _migration(1, "first", 0, "First."))
    _epochs(root, ["first", "second"])
    _schema(root, [])
    _write(root, f"{FIXTURES}/README.md", "# Old schemas\n\nCommit each file came from:\n\n- v0: aaa (start)\n")
    _write(root, "docs/notes.md", "A note about main.\n")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "common")
    _git(root, "checkout", "-q", "-b", "branch")
    _write(root, f"{VERSIONS}/0002_own_change.py", _migration(2, "own_change", 1, "Own change."))
    _epochs(root, ["first", "second", "own reason"])
    _schema(root, ["own"])
    _fixture(root, 1, "-- old own\nCREATE TABLE t (id INT);\n")
    _write(root, f"{FIXTURES}/README.md", "# Old schemas\n\nCommit each file came from:\n\n- v0: aaa (start)\n- v1: bbb (before own change)\n")
    _write(root, "docs/notes.md", "A note about main.\nThe own change is migration 0002 (schema v2) and epoch 3; v1 stays.\n")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "own")
    _git(root, "checkout", "-q", "main")
    _write(root, f"{VERSIONS}/0002_main_change.py", _migration(2, "main_change", 1, "Main change."))
    _epochs(root, ["first", "second", "main reason"])
    _schema(root, ["main"])
    _fixture(root, 1, "-- old main\nCREATE TABLE t (id INT);\n")
    _write(root, f"{FIXTURES}/README.md", "# Old schemas\n\nCommit each file came from:\n\n- v0: aaa (start)\n- v1: ccc (before main change)\n")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "main")
    _git(root, "checkout", "-q", "branch")
    return root


def _epoch_reasons(root):
    with open(os.path.join(root, EPOCHS), encoding="utf-8") as handle:
        return [(entry["epoch"], entry["reason"]) for entry in json.load(handle)["epochs"]]


def test_a_branch_that_took_main_numbers_is_told_so(tool, repo):
    ok, message = tool.check(repo, "main")
    assert not ok
    assert "migration 0002 -> 0003" in message and "epoch 3 -> 4" in message


def test_the_branch_numbers_move_after_main_before_the_merge(tool, repo, capsys):
    assert tool.main(["--root", repo, "--base", "main"]) == 0
    names = sorted(os.listdir(os.path.join(repo, VERSIONS)))
    assert names == ["0001_first.py", "0003_own_change.py"]
    with open(os.path.join(repo, VERSIONS, "0003_own_change.py"), encoding="utf-8") as handle:
        text = handle.read()
    assert 'revision = "0003"' in text and 'down_revision = "0002"' in text and "Schema v3." in text
    assert _epoch_reasons(repo) == [(1, "first"), (2, "second"), (3, "main reason"), (4, "own reason")]
    with open(os.path.join(repo, "docs/notes.md"), encoding="utf-8") as handle:
        note = handle.read()
    assert "migration 0003 (schema v3) and epoch 4; v1 stays" in note and note.startswith("A note about main.\n")
    assert tool.check(repo, "main")[0]
    assert tool.main(["--root", repo, "--base", "main", "--check"]) == 0
    assert "tests.golden_update" in capsys.readouterr().out


def test_the_first_migration_is_tested_from_the_schema_main_ships(tool, repo):
    tool.main(["--root", repo, "--base", "main"])
    names = sorted(name for name in os.listdir(os.path.join(repo, FIXTURES)) if name.endswith(".gz"))
    assert names == ["schema_v1.sql.gz", "schema_v2.sql.gz"]                  # main's v1, and the new one from main's schema.sql
    with gzip.open(os.path.join(repo, FIXTURES, "schema_v2.sql.gz"), "rt", encoding="utf-8") as handle:
        assert handle.read() == "CREATE TABLE t (\n    id INT,\n    main INT\n);\n"
    with open(os.path.join(repo, FIXTURES, "README.md"), encoding="utf-8") as handle:
        readme = handle.read()
    assert "- v2: " in readme and "(before own change)" in readme and "- v1: bbb" not in readme


def test_a_merge_with_the_epoch_file_in_conflict_is_numbered_the_same(tool, repo):
    done = subprocess.run(["git", "merge", "main"], cwd=repo, capture_output=True, text=True)
    assert done.returncode != 0
    assert "<<<<<<<" in open(os.path.join(repo, EPOCHS), encoding="utf-8").read()
    _git(repo, "checkout", "--theirs", "--", f"{FIXTURES}/README.md", f"{FIXTURES}/schema_v1.sql.gz", SCHEMA)
    _git(repo, "add", "-A", "--", f"{FIXTURES}/README.md", f"{FIXTURES}/schema_v1.sql.gz", SCHEMA)
    assert tool.main(["--root", repo, "--base", "main"]) == 0
    assert _epoch_reasons(repo) == [(1, "first"), (2, "second"), (3, "main reason"), (4, "own reason")]
    assert sorted(n for n in os.listdir(os.path.join(repo, VERSIONS)) if n.endswith(".py")) == [
        "0001_first.py", "0002_main_change.py", "0003_own_change.py"]


def test_a_dry_run_changes_nothing(tool, repo):
    before = _git(repo, "status", "--porcelain")
    assert tool.main(["--root", repo, "--base", "main", "--dry-run"]) == 0
    assert _git(repo, "status", "--porcelain") == before
    assert sorted(os.listdir(os.path.join(repo, VERSIONS)))[1] == "0002_own_change.py"


def test_numbers_already_following_main_are_left_alone(tool, repo):
    tool.main(["--root", repo, "--base", "main"])
    snapshot = _git(repo, "status", "--porcelain")
    assert tool.main(["--root", repo, "--base", "main"]) == 0
    assert _git(repo, "status", "--porcelain") == snapshot
    assert _epoch_reasons(repo)[-1] == (4, "own reason")


def test_a_branch_with_two_migrations_moves_both_in_a_chain(tool, repo):
    _write(repo, f"{VERSIONS}/0003_own_more.py", _migration(3, "own_more", 2, "More of the own change."))
    _fixture(repo, 2, "-- between\nCREATE TABLE t (id INT, own INT);\n")
    assert tool.main(["--root", repo, "--base", "main"]) == 0
    names = sorted(name for name in os.listdir(os.path.join(repo, VERSIONS)) if name.endswith(".py"))
    assert names == ["0001_first.py", "0003_own_change.py", "0004_own_more.py"]
    with open(os.path.join(repo, VERSIONS, "0004_own_more.py"), encoding="utf-8") as handle:
        text = handle.read()
    assert 'revision = "0004"' in text and 'down_revision = "0003"' in text and "Schema v4." in text
    fixtures = sorted(name for name in os.listdir(os.path.join(repo, FIXTURES)) if name.endswith(".gz"))
    assert fixtures == ["schema_v1.sql.gz", "schema_v2.sql.gz", "schema_v3.sql.gz"]
    with gzip.open(os.path.join(repo, FIXTURES, "schema_v3.sql.gz"), "rt", encoding="utf-8") as handle:
        assert "own INT" in handle.read()                                      # the one made between the two moved with them


def test_two_migrations_numbered_alike_are_refused(tool, repo, capsys):
    _write(repo, f"{VERSIONS}/0002_other.py", _migration(2, "other", 1, "Other."))
    assert tool.main(["--root", repo, "--base", "main"]) == 2
    assert "numbered alike" in capsys.readouterr().err
