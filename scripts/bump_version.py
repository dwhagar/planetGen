#!/usr/bin/env python3
"""
Stamps pending release notes into a real version number.

Why this exists
===============
Every PR used to pick "main's version + 1" on its own and write that number
into three places at once (`src/stellarObjects/_version.py`'s `__version__`,
`README.md`'s `**Version:**` badge, and a new top `CHANGELOG.md` entry).
Two PRs open at the same time always claimed the same number, and whichever
merged second needed a hand-renumbering merge (e.g. 5.41.0/5.41.1 ->
5.43.0/5.43.1, or 5.46.15 -> 5.46.16).

Now a PR never touches those three places. It adds one note file under
`changes/` instead (see `changes/README.md`), named uniquely per PR so two
PRs can't collide, and `.github/workflows/stamp-version.yml` runs this
script on `main` after every merge: each pending note becomes its own
release (oldest merge first), bumping all three places together.

Usage
=====
    python scripts/bump_version.py              # stamp every pending note
    python scripts/bump_version.py --dry-run    # show what stamping would do
    python scripts/bump_version.py --commit     # stamp, one git commit per release
    python scripts/bump_version.py --check      # validate pending note files only
    python scripts/bump_version.py --check-pr origin/main
                                                # CI's PR guard (see check_pr)

Deliberately standard-library only, so it runs anywhere without installing
the package first.
"""

from __future__ import annotations

import argparse
import datetime
import os
import re
import subprocess
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

VERSION_FILE = os.path.join("src", "stellarObjects", "_version.py")
README_FILE = "README.md"
CHANGELOG_FILE = "CHANGELOG.md"
CHANGES_DIR = "changes"

BUMP_LEVELS = ("patch", "minor", "major")

# The section headings CHANGELOG.md already uses (Keep a Changelog's set).
SECTIONS = ("Added", "Changed", "Deprecated", "Removed", "Fixed", "Security")

FRAGMENT_RE = re.compile(r"^(?P<slug>[A-Za-z0-9][A-Za-z0-9_-]*)\.(?P<level>patch|minor|major)\.md$")
VERSION_RE = re.compile(r"^__version__\s*=\s*(['\"])(?P<version>[^'\"]+)\1", re.M)
BADGE_RE = re.compile(r"(\*\*Version:\*\*\s*)(?P<version>\S+)")
CHANGELOG_TOP_RE = re.compile(r"^## \[(?P<version>[^\]]+)\]", re.M)


class BumpError(Exception):
    """A problem a person needs to fix before anything can be stamped."""


class Fragment:
    """One pending release note: `changes/<slug>.<level>.md`."""

    def __init__(self, filename, level, body):
        self.filename = filename
        self.level = level
        self.body = body

    def __repr__(self):
        return f"Fragment({self.filename!r})"


# -- small file helpers ------------------------------------------------------

def _read(root, relative_path):
    with open(os.path.join(root, relative_path), encoding="utf-8") as f:
        return f.read()


def _write(root, relative_path, contents):
    with open(os.path.join(root, relative_path), "w", encoding="utf-8", newline="\n") as f:
        f.write(contents)


def _git(root, *args):
    return subprocess.run(
        ["git", *args], cwd=root, check=True, capture_output=True, text=True
    ).stdout


# -- version parsing ---------------------------------------------------------

def parse_version(version):
    parts = version.split(".")
    if len(parts) != 3 or not all(p.isdigit() for p in parts):
        raise BumpError(f"version {version!r} isn't in MAJOR.MINOR.PATCH form")
    return tuple(int(p) for p in parts)


def next_version(version, level):
    major, minor, patch = parse_version(version)
    if level == "major":
        return f"{major + 1}.0.0"
    if level == "minor":
        return f"{major}.{minor + 1}.0"
    return f"{major}.{minor}.{patch + 1}"


def version_in(text, pattern, what):
    match = pattern.search(text)
    if not match:
        raise BumpError(f"couldn't find the version in {what}")
    return match.group("version")


def current_versions(root):
    """The version as each of the three places currently states it."""
    return {
        VERSION_FILE: version_in(_read(root, VERSION_FILE), VERSION_RE, VERSION_FILE),
        README_FILE: version_in(_read(root, README_FILE), BADGE_RE, f"{README_FILE}'s **Version:** badge"),
        CHANGELOG_FILE: version_in(_read(root, CHANGELOG_FILE), CHANGELOG_TOP_RE, f"{CHANGELOG_FILE}'s top entry"),
    }


# -- fragments ---------------------------------------------------------------

def validate_fragment_body(filename, body):
    """Returns a list of problems with one note's contents (empty if fine)."""
    problems = []
    lines = [line for line in body.splitlines() if line.strip()]
    if not lines:
        return [f"{CHANGES_DIR}/{filename} is empty"]
    first = lines[0].strip()
    allowed = ", ".join(f"'### {s}'" for s in SECTIONS)
    if not first.startswith("### "):
        problems.append(f"{CHANGES_DIR}/{filename} must start with a section heading ({allowed})")
    for line in lines:
        if line.startswith("### ") and line[4:].strip() not in SECTIONS:
            problems.append(f"{CHANGES_DIR}/{filename} has heading {line.strip()!r}; use one of {allowed}")
        elif re.match(r"^#{1,2} ", line):
            problems.append(
                f"{CHANGES_DIR}/{filename} has a '#'/'##' heading ({line.strip()!r}); "
                "the release heading is added for you, so only use '###' sections"
            )
    return problems


def find_fragments(root):
    """
    Reads every pending note under `changes/`, validating each one.

    Raises:
        BumpError: listing every problem found, if any note is misnamed or
                   malformed, so nothing is stamped half-way.
    """
    directory = os.path.join(root, CHANGES_DIR)
    if not os.path.isdir(directory):
        return []
    fragments, problems = [], []
    for filename in sorted(os.listdir(directory)):
        if filename in ("README.md", ".gitkeep") or filename.startswith("."):
            continue
        match = FRAGMENT_RE.match(filename)
        if not match:
            problems.append(
                f"{CHANGES_DIR}/{filename} isn't named '<short-name>.<patch|minor|major>.md'"
            )
            continue
        body = _read(root, os.path.join(CHANGES_DIR, filename))
        problems.extend(validate_fragment_body(filename, body))
        fragments.append(Fragment(filename, match.group("level"), body))
    if problems:
        raise BumpError("\n".join(problems))
    return fragments


def order_fragments(root, fragments):
    """
    Oldest merge first, so versions follow the order PRs landed on main.

    Uses when each note first appeared along main's first-parent history
    (i.e. the merge that brought it in). Falls back to filename order for
    anything git can't date, such as an uncommitted note.
    """
    def merged_at(fragment):
        path = f"{CHANGES_DIR}/{fragment.filename}"
        try:
            out = _git(root, "log", "--first-parent", "--diff-filter=A", "--format=%ct", "--", path)
        except (OSError, subprocess.CalledProcessError):
            return float("inf")
        stamps = [int(s) for s in out.split()]
        return min(stamps) if stamps else float("inf")

    return sorted(fragments, key=lambda f: (merged_at(f), f.filename))


# -- stamping ----------------------------------------------------------------

def apply_release(root, version, date, body):
    """Writes one release into all three places at once."""
    text = _read(root, VERSION_FILE)
    _write(root, VERSION_FILE, VERSION_RE.sub(lambda m: f"__version__ = {m.group(1)}{version}{m.group(1)}", text, count=1))

    text = _read(root, README_FILE)
    _write(root, README_FILE, BADGE_RE.sub(lambda m: f"{m.group(1)}{version}", text, count=1))

    text = _read(root, CHANGELOG_FILE)
    top = CHANGELOG_TOP_RE.search(text)
    entry = f"## [{version}] - {date}\n\n{body.strip()}\n\n"
    if top:
        text = text[:top.start()] + entry + text[top.start():]
    else:
        text = text.rstrip("\n") + "\n\n" + entry
    _write(root, CHANGELOG_FILE, text)


def stamp(root, date=None, dry_run=False, commit=False, out=sys.stdout):
    """
    Turns each pending note into its own release, oldest first.

    Returns:
        list[str]: the versions released, in order (empty if nothing was
                   pending).
    """
    date = date or datetime.datetime.now(datetime.timezone.utc).date().isoformat()
    fragments = order_fragments(root, find_fragments(root))
    if not fragments:
        print("No pending release notes in changes/; nothing to stamp.", file=out)
        return []

    versions = current_versions(root)
    if len(set(versions.values())) != 1:
        details = ", ".join(f"{path} says {v}" for path, v in versions.items())
        raise BumpError(f"the three version places disagree before stamping ({details}); fix that first")

    version = versions[VERSION_FILE]
    released = []
    for fragment in fragments:
        version = next_version(version, fragment.level)
        released.append(version)
        print(f"{CHANGES_DIR}/{fragment.filename} -> {version} ({fragment.level})", file=out)
        if dry_run:
            continue
        apply_release(root, version, date, fragment.body)
        os.remove(os.path.join(root, CHANGES_DIR, fragment.filename))
        if commit:
            _git(root, "add", "-A", "--", VERSION_FILE, README_FILE, CHANGELOG_FILE,
                 f"{CHANGES_DIR}/{fragment.filename}")
            slug = FRAGMENT_RE.match(fragment.filename).group("slug")
            _git(root, "commit", "-q", "-m", f"Release {version} ({slug})")
    return released


# -- PR guard ----------------------------------------------------------------

def _versions_at(root, ref):
    def show(path):
        return _git(root, "show", f"{ref}:{path}")
    return {
        VERSION_FILE: version_in(show(VERSION_FILE), VERSION_RE, f"{ref}:{VERSION_FILE}"),
        README_FILE: version_in(show(README_FILE), BADGE_RE, f"{ref}:{README_FILE}"),
        CHANGELOG_FILE: version_in(show(CHANGELOG_FILE), CHANGELOG_TOP_RE, f"{ref}:{CHANGELOG_FILE}"),
    }


def check_pr(root, base, allow_no_fragment=False):
    """
    CI's guard for a PR branch checked out at `root`, compared against
    `base` (e.g. `origin/main`).

    Returns:
        list[str]: problems found (empty means the PR is fine).
    """
    problems = []
    before, after = _versions_at(root, base), _versions_at(root, "HEAD")
    for path in (VERSION_FILE, README_FILE, CHANGELOG_FILE):
        if before[path] != after[path]:
            problems.append(
                f"{path} changes the version ({before[path]} -> {after[path]}). PRs no longer "
                f"bump the version; add a note under {CHANGES_DIR}/ instead and it's stamped "
                "after merge (see changes/README.md)."
            )

    try:
        find_fragments(root)
    except BumpError as e:
        problems.extend(str(e).splitlines())

    merge_base = _git(root, "merge-base", base, "HEAD").strip()
    added = _git(root, "diff", "--name-only", "--diff-filter=A", merge_base, "HEAD", "--", f"{CHANGES_DIR}/")
    added_notes = [p for p in added.split() if FRAGMENT_RE.match(os.path.basename(p))]
    if not added_notes and not allow_no_fragment:
        problems.append(
            f"This PR adds no release note under {CHANGES_DIR}/. Add one "
            "(see changes/README.md), or label the PR 'no-release' if it "
            "really shouldn't produce a version."
        )
    return problems


# -- CLI ---------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="only validate pending note files")
    mode.add_argument("--check-pr", metavar="BASE", help="validate a PR branch against BASE (e.g. origin/main)")
    mode.add_argument("--dry-run", action="store_true", help="print what stamping would do, change nothing")
    mode.add_argument("--commit", action="store_true", help="stamp and make one git commit per release")
    parser.add_argument("--allow-no-fragment", action="store_true",
                        help="with --check-pr: don't require the PR to add a note")
    parser.add_argument("--date", help="release date to write (default: today, UTC)")
    parser.add_argument("--root", default=REPO_ROOT, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    try:
        if args.check:
            fragments = find_fragments(args.root)
            print(f"{len(fragments)} pending release note(s), all valid.")
        elif args.check_pr:
            problems = check_pr(args.root, args.check_pr, args.allow_no_fragment)
            if problems:
                print("Release-note check failed:\n" + "\n".join(f"  - {p}" for p in problems), file=sys.stderr)
                return 1
            print("Release-note check passed.")
        else:
            stamp(args.root, date=args.date, dry_run=args.dry_run, commit=args.commit)
    except BumpError as e:
        print(f"bump_version: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
