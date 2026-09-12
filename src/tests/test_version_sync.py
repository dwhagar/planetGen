# tests/test_version_sync.py

"""
Guards against `README.md`'s version badge and `CHANGELOG.md`'s top entry
drifting out of sync with `stellarObjects._version.__version__` -- the
single source of truth `setup.py`/the CLI `--version` banners already read
(see `_version.py`'s own docstring). Nothing enforces either file getting
bumped alongside a real release; this caught README.md sitting 3 releases
stale (5.24.0 while `__version__` was already 5.27.0) with no CI signal at
all until a deploy-readiness check noticed it by hand.
"""

import os
import re

from stellarObjects._version import __version__

# This file lives at src/tests/, two levels under the repo root (src
# layout), same relative-path convention test_examples.py uses.
REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")


def _read(relative_path):
    with open(os.path.join(REPO_ROOT, relative_path), encoding="utf-8") as f:
        return f.read()


def test_readme_version_badge_matches_version_py():
    readme = _read("README.md")
    match = re.search(r"\*\*Version:\*\*\s*(\S+)", readme)
    assert match, "README.md's '**Version:** ...' badge line is missing or reformatted"
    assert match.group(1) == __version__, (
        f"README.md's version badge says {match.group(1)!r}, but "
        f"stellarObjects/_version.py's __version__ is {__version__!r} -- "
        "bump the README badge alongside any version release."
    )


def test_changelog_top_entry_matches_version_py():
    changelog = _read("CHANGELOG.md")
    match = re.search(r"^## \[(\S+)\]", changelog, re.M)
    assert match, "CHANGELOG.md's top '## [x.y.z] - ...' entry heading is missing or reformatted"
    assert match.group(1) == __version__, (
        f"CHANGELOG.md's top entry is [{match.group(1)}], but "
        f"stellarObjects/_version.py's __version__ is {__version__!r} -- "
        "add a new top entry (or bump this one) alongside any version release."
    )
