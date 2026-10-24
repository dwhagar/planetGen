# tests/test_bump_version.py

"""
Tests for `scripts/bump_version.py`, which stamps `changes/` release notes
into a real version after merge (see `changes/README.md`). Each test builds
a tiny stand-in repo in `tmp_path` holding just the three places a version
lives, so nothing here touches the real checkout's files.
"""

import importlib.util
import io
import json
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

@pytest.fixture(autouse=True)
def _no_revision_hold(monkeypatch):
    """The tests below describe the normal scheme; the hold has its own."""
    monkeypatch.setattr(bump_version, "REVISION_HOLD", False)


HEADLINE = "**Headline:** A thing works better.\n\n"
FIXED_NOTE = HEADLINE + "### Fixed\n- Fixed a thing.\n"
ADDED_NOTE = HEADLINE + "### Added\n- Added a thing.\n"


def _make_repo(root, version="1.2.3"):
    os.makedirs(os.path.join(root, os.path.dirname(bump_version.VERSION_FILE)))
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
    ("patch", "1.3.4"), ("minor", "1.3.4"), ("major", "2.0.4"),
])
def test_next_version(level, expected):
    assert bump_version.next_version("1.2.3", level) == expected


def test_stamp_updates_all_three_places_and_removes_the_note(tmp_path):
    root = _make_repo(str(tmp_path))
    _note(root, "thing.patch.md", FIXED_NOTE)

    assert bump_version.stamp(root, date="2026-02-02") == ["1.3.4"]

    assert '__version__ = "1.3.4"' in _read(root, bump_version.VERSION_FILE)
    assert "**Version:** 1.3.4 &middot;" in _read(root, "README.md")
    changelog = _read(root, "CHANGELOG.md")
    assert changelog.startswith(
        "# Changelog\n\n## [1.3.4] - 2026-02-02\n\n**Headline:** A thing works better.\n\n### Fixed\n- Fixed a thing.\n\n## [1.2.3]")
    assert not os.path.exists(os.path.join(root, "changes", "thing.patch.md"))
    assert os.path.exists(os.path.join(root, "changes", "README.md"))
    assert bump_version.current_versions(root) == {
        bump_version.VERSION_FILE: "1.3.4", "README.md": "1.3.4", "CHANGELOG.md": "1.3.4"}


def test_each_note_gets_its_own_release(tmp_path):
    root = _make_repo(str(tmp_path))
    _note(root, "a-fix.patch.md", FIXED_NOTE)
    _note(root, "b-feature.minor.md", ADDED_NOTE)

    # No git here, so filename order decides.
    assert bump_version.stamp(root, date="2026-02-02") == ["1.3.4", "1.4.5"]
    changelog = _read(root, "CHANGELOG.md")
    assert changelog.index("## [1.4.5]") < changelog.index("## [1.3.4]") < changelog.index("## [1.2.3]")


def test_no_notes_is_a_no_op(tmp_path):
    root = _make_repo(str(tmp_path))
    before = _read(root, "CHANGELOG.md")
    assert bump_version.stamp(root) == []
    assert _read(root, "CHANGELOG.md") == before


def test_dry_run_changes_nothing(tmp_path):
    root = _make_repo(str(tmp_path))
    _note(root, "thing.patch.md", FIXED_NOTE)
    assert bump_version.stamp(root, dry_run=True) == ["1.3.4"]
    assert '__version__ = "1.2.3"' in _read(root, bump_version.VERSION_FILE)
    assert os.path.exists(os.path.join(root, "changes", "thing.patch.md"))


@pytest.mark.parametrize("name, body, complaint", [
    ("thing.md", FIXED_NOTE, "isn't named"),
    ("thing.bugfix.md", FIXED_NOTE, "isn't named"),
    ("thing.patch.md", "", "is empty"),
    ("thing.patch.md", "- no heading\n", "must start with a Headline"),
    ("thing.patch.md", "### Fixed\n- no headline\n", "must start with a Headline"),
    ("thing.patch.md", "**Headline:**\n### Fixed\n- empty headline\n", "must start with a Headline"),
    ("thing.patch.md", HEADLINE + "- no heading\n", "must start with a section heading"),
    ("thing.patch.md", HEADLINE + "### Fixes\n- typo'd section\n", "use one of"),
    ("thing.patch.md", HEADLINE + "## [9.9.9] - 2026-01-01\n### Fixed\n- x\n", "'#'/'##' heading"),
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

    assert bump_version.stamp(root, date="2026-02-02", commit=True) == ["1.3.4", "1.4.5"]

    log = subprocess.run(["git", "log", "--format=%s", "-2"], cwd=root,
                         check=True, capture_output=True, text=True).stdout.split("\n")
    assert log[:2] == ["Release 1.4.5 (a-second)", "Release 1.3.4 (z-first)"]
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
def test_check_pr_rejects_a_major_note_while_the_revision_is_on_hold(tmp_path, monkeypatch):
    monkeypatch.setattr(bump_version, "REVISION_HOLD", True)
    root = _make_repo(str(tmp_path))
    _init_git(root)
    _git(root, "checkout", "-q", "-b", "feature")
    _note(root, "feature.major.md", ADDED_NOTE)
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "feature")
    problems = bump_version.check_pr(root, "main")
    assert any("'major' note" in p and "on hold" in p for p in problems)


@needs_git
def test_check_pr_rejects_a_hand_bumped_version(tmp_path):
    root = _make_repo(str(tmp_path))
    _init_git(root)
    _git(root, "checkout", "-q", "-b", "feature")
    _note(root, "feature.patch.md", FIXED_NOTE)
    bump_version.apply_release(root, "1.3.4", "2026-02-02", FIXED_NOTE)
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


# -- BUILD counts releases; the TODO counters only allocate IDs ------------------------------

def _todo_map(root, next_free, todo_ids=()):
    """A stand-in docs/plan/todo-number-map.md (and docs/TODO.md)."""
    os.makedirs(os.path.join(root, "docs", "plan"), exist_ok=True)
    rows = "".join(f"| {cat} | {cat}.{n} |\n" for cat, n in next_free.items())
    with open(os.path.join(root, bump_version.TODO_MAP_FILE), "w") as f:
        f.write("# Map\n\n## Next free IDs\n\n| Category | Next free ID |\n|---|---|\n" + rows)
    with open(os.path.join(root, bump_version.TODO_FILE), "w") as f:
        f.write("# TODO\n\n" + "".join(f"- [ ] **{i} Something**\n" for i in todo_ids))


def _all_next_free(**overrides):
    table = {cat: 1 for cat in bump_version.TODO_CATEGORIES}
    table.update(overrides)
    return table


def test_build_goes_up_by_one_on_every_release_and_never_resets(tmp_path):
    root = _make_repo(str(tmp_path), version="1.2.3")
    _todo_map(root, _all_next_free(UX=5, API=9, DB=13))  # the counters no longer matter
    _note(root, "a-fix.patch.md", FIXED_NOTE)
    _note(root, "b-feature.minor.md", ADDED_NOTE)
    _note(root, "c-break.major.md", ADDED_NOTE)
    assert bump_version.stamp(root, date="2026-02-02") == ["1.3.4", "1.4.5", "2.0.6"]
    assert '__version__ = "2.0.6"' in _read(root, bump_version.VERSION_FILE)


def test_a_stale_next_free_table_is_rejected(tmp_path):
    root = _make_repo(str(tmp_path))
    _todo_map(root, _all_next_free(MAP=5), ["MAP.5"])
    _note(root, "thing.patch.md", FIXED_NOTE)
    with pytest.raises(bump_version.BumpError, match="next free MAP ID is MAP.5"):
        bump_version.todo_counters(root)


def test_a_missing_category_row_is_rejected(tmp_path):
    root = _make_repo(str(tmp_path))
    table = _all_next_free()
    del table["POP"]
    _todo_map(root, table)
    with pytest.raises(bump_version.BumpError, match="no row for POP"):
        bump_version.todo_counters(root)


def test_real_todo_counters_are_consistent():
    """The real 'Next free IDs' table covers every category and is ahead of docs/TODO.md."""
    assert sum(bump_version.todo_counters(REPO_ROOT).values()) > 0


def test_revision_hold_keeps_the_revision_and_treats_major_as_minor(tmp_path, monkeypatch):
    monkeypatch.setattr(bump_version, "REVISION_HOLD", True)
    assert bump_version.next_version("8.0.711", "patch") == "8.0.712"
    assert bump_version.next_version("8.0.711", "minor") == "8.0.712"
    assert bump_version.next_version("8.0.711", "major") == "8.0.712"

    root = _make_repo(str(tmp_path), version="8.0.711")
    _note(root, "a-feature.minor.md", ADDED_NOTE)
    _note(root, "b-fix.patch.md", FIXED_NOTE)
    assert bump_version.stamp(root, date="2026-02-02") == ["8.0.712", "8.0.713"]
    changelog = _read(root, "CHANGELOG.md")
    assert changelog.index("## [8.0.713]") < changelog.index("## [8.0.712]") < changelog.index("## [8.0.711]")
    assert '__version__ = "8.0.713"' in _read(root, bump_version.VERSION_FILE)


def test_merge_bodies_puts_new_bullets_first_in_each_section():
    merged = bump_version.merge_bodies(
        "### Added\n- New.\n### Fixed\n- New fix.\n",
        "\n\n### Fixed\n- Old fix.\n### Changed\n- Old change.\n")
    assert merged == ("### Added\n- New.\n\n### Fixed\n- New fix.\n- Old fix.\n\n"
                      "### Changed\n- Old change.")


def test_a_headline_line_leads_the_note_and_lands_under_the_release_heading(tmp_path):
    root = _make_repo(str(tmp_path))
    _note(root, "thing.minor.md", "**Headline:** The map opens on the biggest group.\n\n### Added\n- Added a thing.\n")
    bump_version.stamp(root, date="2026-02-02")
    changelog = _read(root, bump_version.CHANGELOG_FILE)
    assert ("## [1.3.4] - 2026-02-02\n\n**Headline:** The map opens on the biggest group.\n\n### Added"
            in changelog)


@pytest.mark.parametrize("body, complaint", [
    ("**Headline:** Only a headline.\n", "no section"),
    ("**Headline:** Then bullets.\n- no heading\n", "must start with a section heading"),
    ("**Headline:**\n### Fixed\n- x\n", "must start with a Headline"),
])
def test_a_headline_still_needs_a_section_after_it(tmp_path, body, complaint):
    root = _make_repo(str(tmp_path))
    _note(root, "thing.patch.md", body)
    with pytest.raises(bump_version.BumpError, match=complaint):
        bump_version.stamp(root)


def test_merge_bodies_keeps_the_newest_headline():
    merged = bump_version.merge_bodies(
        "**Headline:** New one.\n\n### Added\n- New.\n",
        "\n\n**Headline:** Old one.\n\n### Added\n- Old.\n")
    assert merged == "**Headline:** New one.\n\n### Added\n- New.\n- Old."
    assert bump_version.merge_bodies("### Added\n- New.\n", "\n**Headline:** Old one.\n### Added\n- Old.\n") \
        .startswith("**Headline:** Old one.\n\n### Added")


# -- the generator epoch (OPS.56) ----------------------------------------------------------------

GENERATION_NOTE = FIXED_NOTE.rstrip("\n") + "\ngeneration-output: changed\n"


def _epochs(root, count, pending=False):
    os.makedirs(os.path.join(root, bump_version.GOLDEN_DIR), exist_ok=True)
    entries = [{"epoch": n, "first_release": None, "reason": "r"} for n in range(1, count + 1)]
    if pending:
        entries[-1]["pending"] = True
    with open(os.path.join(root, bump_version.EPOCHS_FILE), "w") as f:
        json.dump({"epochs": entries}, f, indent=2)


def _golden(root, text):
    with open(os.path.join(root, bump_version.GOLDEN_DIR, "sector.txt"), "w") as f:
        f.write(text)


def _epoch_repo(tmp_path):
    root = _make_repo(str(tmp_path))
    _epochs(root, 1)
    _golden(root, "one")
    _init_git(root)
    _git(root, "checkout", "-q", "-b", "feature")
    return root


def _commit(root, note=None):
    if note:
        _note(root, "feature.patch.md", note)
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "change")


@needs_git
def test_check_pr_passes_golden_epoch_and_note_line_together(tmp_path):
    root = _epoch_repo(tmp_path)
    _golden(root, "two")
    _epochs(root, 2, pending=True)
    _commit(root, GENERATION_NOTE)
    assert bump_version.check_pr(root, "main") == []


@needs_git
def test_check_pr_passes_a_pr_that_touches_none_of_them(tmp_path):
    root = _epoch_repo(tmp_path)
    _commit(root, FIXED_NOTE)
    assert bump_version.check_pr(root, "main") == []


@needs_git
@pytest.mark.parametrize("golden, epoch, line, missing", [
    (True, False, False, "epoch bumped"),
    (True, True, False, "generation-output: changed"),
    (True, False, True, "epoch bumped"),
    (False, True, True, "golden files changed"),
    (False, True, False, "golden files changed"),
    (False, False, True, "golden files changed"),
])
def test_check_pr_rejects_a_partial_epoch_change(tmp_path, golden, epoch, line, missing):
    root = _epoch_repo(tmp_path)
    if golden:
        _golden(root, "two")
    if epoch:
        _epochs(root, 2, pending=True)
    _commit(root, GENERATION_NOTE if line else FIXED_NOTE)
    problems = bump_version.check_pr(root, "main")
    assert len(problems) == 1 and "go together" in problems[0] and missing in problems[0], problems


@needs_git
def test_check_pr_rejects_a_jump_of_two_epochs(tmp_path):
    root = _epoch_repo(tmp_path)
    _golden(root, "two")
    _epochs(root, 3, pending=True)
    _commit(root, GENERATION_NOTE)
    assert any("by one at most" in p for p in bump_version.check_pr(root, "main"))


def _no_golden_repo(tmp_path):
    root = _make_repo(str(tmp_path))
    _epochs(root, 1)
    _init_git(root)
    _git(root, "checkout", "-q", "-b", "feature")
    return root


@needs_git
def test_with_no_golden_file_checked_in_the_epoch_and_the_note_go_together(tmp_path):
    root = _no_golden_repo(tmp_path)
    _epochs(root, 2, pending=True)
    _commit(root, GENERATION_NOTE)
    assert bump_version.check_pr(root, "main") == []


@needs_git
@pytest.mark.parametrize("epoch, line, missing", [(True, False, "generation-output"), (False, True, "epoch bumped")])
def test_with_no_golden_file_a_lone_epoch_or_note_line_is_still_refused(tmp_path, epoch, line, missing):
    root = _no_golden_repo(tmp_path)
    if epoch:
        _epochs(root, 2, pending=True)
    _commit(root, GENERATION_NOTE if line else FIXED_NOTE)
    problems = bump_version.check_pr(root, "main")
    assert len(problems) == 1 and missing in problems[0] and "golden files changed" not in problems[0].split("but not")[1], problems


@needs_git
def test_the_pr_that_adds_the_counter_is_not_held_to_the_rule(tmp_path):
    root = _make_repo(str(tmp_path))
    _init_git(root)
    _git(root, "checkout", "-q", "-b", "feature")
    _epochs(root, 1)
    _golden(root, "one")
    _commit(root, FIXED_NOTE)
    assert bump_version.check_pr(root, "main") == []


def test_stamping_a_generation_note_names_the_release_of_the_pending_epoch(tmp_path):
    root = _make_repo(str(tmp_path), version="1.2.3")
    _epochs(root, 2, pending=True)
    _note(root, "a-fix.patch.md", FIXED_NOTE)
    _note(root, "b-change.patch.md", GENERATION_NOTE)
    bump_version.stamp(root, date="2026-02-02", out=io.StringIO())
    epochs = json.load(open(os.path.join(root, bump_version.EPOCHS_FILE)))["epochs"]
    assert epochs[1]["first_release"] == "1.4.5" and "pending" not in epochs[1]   # the second note's release
    assert epochs[0]["first_release"] is None


def test_the_real_epoch_file_is_well_formed_and_the_accessor_reads_it():
    from planetgen.util import generator_epoch

    assert generator_epoch.generator_epoch() == len(generator_epoch.epoch_history()) >= 1
    assert generator_epoch.first_release(1) is None
    with pytest.raises(ValueError):
        generator_epoch.first_release(0)
