# tests/test_regenerate_generated.py

"""TEST.721: `scripts/regenerate_generated.py` knows every generated file's builder, reports stale ones, and repairs a merge conflict in a generated file by rebuilding it."""

import importlib.util
import pathlib
import subprocess
import sys

import pytest

SCRIPT = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "regenerate_generated.py"


@pytest.fixture(scope="module")
def tool():
    spec = importlib.util.spec_from_file_location("regenerate_generated", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["regenerate_generated"] = module
    spec.loader.exec_module(module)
    return module


BUILDER = '''import sys
root = __import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
source = open(root + "/srca.txt").read() + open(root + "/srcb.txt").read()
built = "".join(f"built: {line}\\n" for line in source.splitlines())
if "--check" in sys.argv:
    sys.exit(0 if open(root + "/out.txt").read() == built else 1)
open(root + "/out.txt", "w").write(built)
'''


def _git(repo, *args):
    return subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, check=False)


@pytest.fixture
def repo(tmp_path, tool, monkeypatch):
    """A throwaway repository with one source, one generated file and its builder; two branches that each add a source line and so conflict on the generated file."""
    (tmp_path / "scripts").mkdir()
    (tmp_path / "scripts" / "build_fake.py").write_text(BUILDER, encoding="utf-8")
    (tmp_path / "srca.txt").write_text("a\n", encoding="utf-8")
    (tmp_path / "srcb.txt").write_text("b\n", encoding="utf-8")
    (tmp_path / "other.txt").write_text("one\n", encoding="utf-8")
    for args in (("init", "-q", "-b", "main"), ("config", "user.email", "t@example.com"), ("config", "user.name", "t"), ("config", "commit.gpgsign", "false")):
        _git(tmp_path, *args)
    _run_build(tmp_path)
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-qm", "base")
    for branch, name, text in (("a", "srca.txt", "a changed\n"), ("b", "srcb.txt", "b changed\n")):
        _git(tmp_path, "checkout", "-qb", branch, "main")
        (tmp_path / name).write_text(text, encoding="utf-8")
        _run_build(tmp_path)
        _git(tmp_path, "commit", "-qam", branch)
    _git(tmp_path, "checkout", "-q", "a")
    monkeypatch.setattr(tool, "ROOT", str(tmp_path))
    monkeypatch.setattr(tool, "GENERATED", (tool.Generated("fake", "build_fake.py", ("out.txt",)),))
    return tmp_path


def _run_build(repo):
    subprocess.run([sys.executable, str(repo / "scripts" / "build_fake.py")], check=True)


def test_every_generated_file_has_a_builder_that_exists_and_takes_check(tool):
    root = SCRIPT.parent.parent
    for entry in tool.GENERATED:
        script = root / "scripts" / entry.script
        assert script.exists(), entry
        assert "--check" in script.read_text(encoding="utf-8"), entry
        assert tool.owner_of(entry.paths[0] + ("x" if entry.paths[0].endswith("/") else "")) is entry


def test_the_global_documents_are_not_in_the_set(tool):
    assert tool.owner_of("docs/generated/architecture-map.html") is None and tool.owner_of("docs/site/index.html") is None


def test_a_conflict_in_a_generated_file_is_rebuilt_from_the_merged_sources(tool, repo):
    assert _git(repo, "merge", "b").returncode != 0
    assert tool.unmerged_paths() == ["out.txt"]
    assert tool.command_resolve(list(tool.GENERATED)) == 0
    assert tool.unmerged_paths() == []
    out = (repo / "out.txt").read_text(encoding="utf-8")
    assert "<<<<" not in out and "built: a changed" in out and "built: b changed" in out
    assert tool.command_check(list(tool.GENERATED)) == 0
    assert _git(repo, "commit", "-qm", "merge b").returncode == 0


def test_resolve_refuses_while_a_source_conflict_is_left(tool, repo):
    (repo / "other.txt").write_text("a\n", encoding="utf-8")
    _git(repo, "commit", "-qam", "a other")
    _git(repo, "checkout", "-q", "b")
    (repo / "other.txt").write_text("b\n", encoding="utf-8")
    _git(repo, "commit", "-qam", "b other")
    _git(repo, "checkout", "-q", "a")
    assert _git(repo, "merge", "b").returncode != 0
    assert set(tool.unmerged_paths()) == {"out.txt", "other.txt"}
    assert tool.command_resolve(list(tool.GENERATED)) == 2
    assert "out.txt" in tool.unmerged_paths()               # nothing was touched


def test_check_names_a_stale_file(tool, repo, capsys):
    (repo / "srca.txt").write_text("changed\n", encoding="utf-8")
    assert tool.command_check(list(tool.GENERATED)) == 1
    assert "out of date: fake" in capsys.readouterr().out
    assert tool.command_regenerate(list(tool.GENERATED)) == 0
    assert tool.command_check(list(tool.GENERATED)) == 0
