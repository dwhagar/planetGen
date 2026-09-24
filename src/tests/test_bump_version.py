# tests/test_bump_version.py

"""
Tests for `scripts/bump_version.py`, which stamps `changes/` release notes
into a real version after merge (see `changes/README.md`). Each test builds
a tiny stand-in repo in `tmp_path` holding just the three places a version
lives, so nothing here touches the real checkout's files.
"""

import importlib.util
import os
import shutil
import subprocess

import pytest

# Same two-levels-up convention test_version_sync.py uses (src layout).
REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")

_spec = importlib.util.spec_from_file_location(
    "bump_version", os.path.join(REPO_ROOT, "scripts", "bump_version.py"))
bump_version = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bump_version)

FIXED_NOTE = "### Fixed\n- Fixed a thing.\n"
ADDED_NOTE = "### Added\n- Added a thing.\n"


def _make_repo(root, version="1.2.3"):
    os.makedirs(os.path.join(root, "src", "stellarObjects"))
    os.makedirs(os.path.join(root, "changes"))
    with open(os.path.join(root, bump_version.VERSION_FILE), "w") as f:
        f.write(f'"""Docstring."""\n\n__version__ = "{version}"\n\nREPO_URL = "x"\n')
    with open(os.path.join(root, "README.md"), "w") as f:
        f.write(f"# Title\n\n**Version:** {version} &middot; [Changelog](CHANGELOG.md)\n")
    with open(os.path.join(root, "CHANGELOG.md"), "w") as f:
        f.write(f"# Changelog\n\n## [{version}] - 2026-01-01\n\n### Fixed\n- Old fix.\n")
    with open(os.path.join(root, "changes", "README.md"), "w") as f:
        f.write("# Pending release notes\n")
    return root


def _note(root, name, body):
    with open(os.path.join(root, "changes", name), "w") as f:
        f.write(body)


def _git(root, *args):
    subprocess.run(["git", *args], cwd=root, check=True, capture_output=True)


def _init_git(root):
    _git(root, "init", "-q", "-b", "main")
    _git(root, "config", "user.email", "t@example.com")
    _git(root, "config", "user.name", "t")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "base")


def _read(root, path):
    with open(os.path.join(root, path)) as f:
        return f.read()


needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="git not installed")


@pytest.mark.parametrize("level, expected", [
    ("patch", "1.2.4"), ("minor", "1.3.0"), ("major", "2.0.0"),
])
def test_next_version(level, expected):
    assert bump_version.next_version("1.2.3", level) == expected


def test_stamp_updates_all_three_places_and_removes_the_note(tmp_path):
    root = _make_repo(str(tmp_path))
    _note(root, "thing.patch.md", FIXED_NOTE)

    assert bump_version.stamp(root, date="2026-02-02") == ["1.2.4"]

    assert '__version__ = "1.2.4"' in _read(root, bump_version.VERSION_FILE)
    assert "**Version:** 1.2.4 &middot;" in _read(root, "README.md")
    changelog = _read(root, "CHANGELOG.md")
    assert changelog.startswith(
        "# Changelog\n\n## [1.2.4] - 2026-02-02\n\n### Fixed\n- Fixed a thing.\n\n## [1.2.3]")
    assert not os.path.exists(os.path.join(root, "changes", "thing.patch.md"))
    assert os.path.exists(os.path.join(root, "changes", "README.md"))
    assert bump_version.current_versions(root) == {
        bump_version.VERSION_FILE: "1.2.4", "README.md": "1.2.4", "CHANGELOG.md": "1.2.4"}


def test_each_note_gets_its_own_release(tmp_path):
    root = _make_repo(str(tmp_path))
    _note(root, "a-fix.patch.md", FIXED_NOTE)
    _note(root, "b-feature.minor.md", ADDED_NOTE)

    # No git here, so filename order decides.
    assert bump_version.stamp(root, date="2026-02-02") == ["1.2.4", "1.3.0"]
    changelog = _read(root, "CHANGELOG.md")
    assert changelog.index("## [1.3.0]") < changelog.index("## [1.2.4]") < changelog.index("## [1.2.3]")


def test_no_notes_is_a_no_op(tmp_path):
    root = _make_repo(str(tmp_path))
    before = _read(root, "CHANGELOG.md")
    assert bump_version.stamp(root) == []
    assert _read(root, "CHANGELOG.md") == before


def test_dry_run_changes_nothing(tmp_path):
    root = _make_repo(str(tmp_path))
    _note(root, "thing.patch.md", FIXED_NOTE)
    assert bump_version.stamp(root, dry_run=True) == ["1.2.4"]
    assert '__version__ = "1.2.3"' in _read(root, bump_version.VERSION_FILE)
    assert os.path.exists(os.path.join(root, "changes", "thing.patch.md"))


@pytest.mark.parametrize("name, body, complaint", [
    ("thing.md", FIXED_NOTE, "isn't named"),
    ("thing.bugfix.md", FIXED_NOTE, "isn't named"),
    ("thing.patch.md", "", "is empty"),
    ("thing.patch.md", "- no heading\n", "must start with a section heading"),
    ("thing.patch.md", "### Fixes\n- typo'd section\n", "use one of"),
    ("thing.patch.md", "## [9.9.9] - 2026-01-01\n### Fixed\n- x\n", "'#'/'##' heading"),
])
def test_bad_notes_are_rejected_before_anything_is_stamped(tmp_path, name, body, complaint):
    root = _make_repo(str(tmp_path))
    _note(root, "good.patch.md", FIXED_NOTE)
    _note(root, name, body)
    with pytest.raises(bump_version.BumpError, match=complaint):
        bump_version.stamp(root)
    assert '__version__ = "1.2.3"' in _read(root, bump_version.VERSION_FILE)


def test_refuses_to_stamp_when_the_three_places_already_disagree(tmp_path):
    root = _make_repo(str(tmp_path))
    with open(os.path.join(root, "README.md"), "w") as f:
        f.write("**Version:** 1.2.2\n")
    _note(root, "thing.patch.md", FIXED_NOTE)
    with pytest.raises(bump_version.BumpError, match="disagree"):
        bump_version.stamp(root)


@needs_git
def test_notes_are_released_in_merge_order_and_committed(tmp_path):
    root = _make_repo(str(tmp_path))
    _init_git(root)
    # "z-..." merges first, so it must get the lower version despite sorting last by name.
    _note(root, "z-first.patch.md", FIXED_NOTE)
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "first")
    _note(root, "a-second.minor.md", ADDED_NOTE)
    _git(root, "add", "-A")
    env_date = {**os.environ, "GIT_COMMITTER_DATE": "2030-01-01T00:00:00Z"}
    subprocess.run(["git", "commit", "-q", "-m", "second"], cwd=root, check=True, env=env_date)

    assert bump_version.stamp(root, date="2026-02-02", commit=True) == ["1.2.4", "1.3.0"]

    log = subprocess.run(["git", "log", "--format=%s", "-2"], cwd=root,
                         check=True, capture_output=True, text=True).stdout.split("\n")
    assert log[:2] == ["Release 1.3.0 (a-second)", "Release 1.2.4 (z-first)"]
    status = subprocess.run(["git", "status", "--porcelain"], cwd=root,
                            check=True, capture_output=True, text=True).stdout
    assert status == ""


@needs_git
def test_check_pr_passes_a_pr_that_only_adds_a_note(tmp_path):
    root = _make_repo(str(tmp_path))
    _init_git(root)
    _git(root, "checkout", "-q", "-b", "feature")
    _note(root, "feature.patch.md", FIXED_NOTE)
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "feature")
    assert bump_version.check_pr(root, "main") == []


@needs_git
def test_check_pr_rejects_a_hand_bumped_version(tmp_path):
    root = _make_repo(str(tmp_path))
    _init_git(root)
    _git(root, "checkout", "-q", "-b", "feature")
    _note(root, "feature.patch.md", FIXED_NOTE)
    bump_version.apply_release(root, "1.2.4", "2026-02-02", FIXED_NOTE)
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "old-style bump")
    problems = bump_version.check_pr(root, "main")
    assert len(problems) == 3
    assert all("no longer bump the version" in p for p in problems)


@needs_git
def test_check_pr_requires_a_note_unless_allowed(tmp_path):
    root = _make_repo(str(tmp_path))
    _init_git(root)
    _git(root, "checkout", "-q", "-b", "feature")
    with open(os.path.join(root, "other.txt"), "w") as f:
        f.write("x")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "no note")
    assert any("adds no release note" in p for p in bump_version.check_pr(root, "main"))
    assert bump_version.check_pr(root, "main", allow_no_fragment=True) == []


def test_real_pending_notes_are_valid():
    """Catches a malformed note in a PR during the normal test run, before merge."""
    bump_version.find_fragments(REPO_ROOT)
