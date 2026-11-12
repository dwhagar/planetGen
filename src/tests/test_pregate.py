"""TEST.720: the pre-gate quick lane (`scripts/pregate.py`) finds what it is there to find."""

import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import pregate  # noqa: E402


def _git(repo, *args):
    return subprocess.run(["git", *args], cwd=repo, check=True, capture_output=True, text=True).stdout


def _repo(tmp_path):
    _git(tmp_path, "init", "-q", "-b", "main")
    _git(tmp_path, "config", "user.email", "t@example.com")
    _git(tmp_path, "config", "user.name", "t")
    (tmp_path / "a.txt").write_text("one\n")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "base")


def test_a_function_defined_twice_in_a_changed_file_is_named(tmp_path, monkeypatch):
    (tmp_path / "m.py").write_text("def _duration(x):\n    return x\n\n\ndef other():\n    pass\n\n\ndef _duration(y):\n    return y\n")
    monkeypatch.setattr(pregate, "ROOT", tmp_path)
    found = pregate.duplicate_definitions(["m.py"])
    assert found == ["m.py: _duration at lines 1 and 9"]


def test_overloads_property_setters_and_registered_functions_are_not_duplicates(tmp_path, monkeypatch):
    (tmp_path / "m.py").write_text(
        "import typing\n\n\n@typing.overload\ndef f(a): ...\n\n\n@typing.overload\ndef f(a, b): ...\n\n\n"
        "class C:\n    @property\n    def v(self):\n        return 1\n\n    @v.setter\n    def v(self, x):\n        pass\n"
    )
    monkeypatch.setattr(pregate, "ROOT", tmp_path)
    assert pregate.duplicate_definitions(["m.py"]) == []


def test_a_conflict_with_the_base_is_found_and_a_clean_merge_is_not(tmp_path, monkeypatch):
    _repo(tmp_path)
    monkeypatch.setattr(pregate, "ROOT", tmp_path)
    _git(tmp_path, "checkout", "-q", "-b", "feature")
    (tmp_path / "a.txt").write_text("feature\n")
    _git(tmp_path, "commit", "-q", "-am", "feature")
    _git(tmp_path, "checkout", "-q", "main")
    assert pregate.merges_cleanly("main") == (True, [])
    (tmp_path / "a.txt").write_text("main\n")
    _git(tmp_path, "commit", "-q", "-am", "main moves")
    _git(tmp_path, "checkout", "-q", "feature")
    clean, files = pregate.merges_cleanly("main")
    assert clean is False and files == ["a.txt"]


def test_the_set_files_name_real_test_files_and_the_prs_own_tests_are_added():
    fixed = pregate.test_files([])
    assert fixed and all((ROOT / name).exists() for name in fixed)
    assert "src/tests/test_american_spelling.py" in fixed
    added = pregate.test_files(["src/tests/test_pregate.py", "src/planetgen/util/log.py", "src/tests/test_does_not_exist.py"])
    assert "src/tests/test_pregate.py" in added and "src/tests/test_does_not_exist.py" not in added
    assert len(added) == len(set(added))
