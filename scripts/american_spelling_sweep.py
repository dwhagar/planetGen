#!/usr/bin/env python3
"""
The American spelling of the whole tree (DOC.103): code, identifiers, comments, tests and documents.

    python scripts/american_spelling_sweep.py apply     # rewrite the tree in place (and rename the files whose names hold a British word)
    python scripts/american_spelling_sweep.py check     # list what is left; exit 1 if anything is

`apply` can be run on any tree at any time: it only touches a word that is still British, so a second run changes nothing. It is how DOC.103 is
merged on a train: merge the branch taking the train's side of every conflict, then run `apply` and `check`, rebuild the generated files
(`scripts/build_reference_pages.py`, `scripts/build_api_compatibility.py`) and commit.

The words: centre (centres, centred, centring), cancelled, cancelling, vapour, flavour, ionising, ionised, ionisation, oxidised, oxidising, disc
(discs, as a whole word: `nuclear_disc` counts, `discrete` does not) and link-neighbours. A proper noun that holds "Centre" (Centre National, Centre for, ...)
stays. Left alone (`SKIP_PREFIX`, `SKIP_FILES`, `SKIP_LINES`): the design documents, the plan and the TODO (and `docs/site/`, the web pages built from them), the changelog and the release notes, the
generated pages, the migrations and the old-schema fixtures (they name the stored words as they were), the reference data (paper titles), the
golden data, the licence text, the vendor scripts, the files that hold the British word list on purpose, and the lines that tell what a word was.
"""

import os
import re
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

SKIP_PREFIX = (
    "docs/TODO.md", "docs/plan/", "docs/design/", "docs/generated/", "docs/site/", "docs/references/", "docs/REFERENCES.md", "docs/REFERENCES.html", "docs/reference/references-data.json",
    "docs/reference/api-versions.json", "docs/reference/api-compatibility.md", "CHANGELOG.md", "changes/", "DATA-LICENSE.md", ".github/",
    "src/planetgen/db/migrations/versions/", "src/tests/fixtures/", "src/tests/golden/", "src/html/static/vendor/",
)
"""tuple: Paths (files or folders, from the repository root) that keep their words."""

SKIP_FILES = (
    "src/planetgen/util/plainwords.py", "src/tests/test_plain_words.py", "src/tests/test_american_spelling.py",
    "src/tests/test_migration_american_spelling.py", "scripts/american_spelling_sweep.py",
)
"""tuple: The files that hold the British words on purpose (the word list, its tests, this script)."""

SKIP_LINES = re.compile(r"old spelling|^- \*\*13\*\* \((?:GEN\.1212, )?DOC\.103\)|unified-system-barycentre")
"""re.Pattern: A line that tells what a word was, or names a file that keeps its name, is left as it is."""

EXTENSIONS = (".py", ".js", ".json", ".html", ".md", ".sql", ".txt", ".sh", ".css", ".yml", ".yaml", ".cfg", ".toml", ".ini", ".mjs")
PROPER = r"(?!\s+(?:National|Nationale|for|de|d'|of)\b)"

RULES = [(re.compile(pattern), replacement) for pattern, replacement in (
    (r"\bCentre" + PROPER, "Center"),
    (r"centred", "centered"), (r"Centred", "Centered"), (r"CENTRED", "CENTERED"),
    (r"centring", "centering"), (r"Centring", "Centering"),
    (r"centre(s?)(?![a-z])", r"center\1"), (r"CENTRE(S?)(?![A-Z])", r"CENTER\1"),
    (r"(?<=[a-z])Centre(s?)(?![a-z])", r"Center\1"),
    (r"cancelled", "canceled"), (r"Cancelled", "Canceled"), (r"CANCELLED", "CANCELED"),
    (r"cancelling", "canceling"), (r"Cancelling", "Canceling"), (r"CANCELLING", "CANCELING"),
    (r"vapour", "vapor"), (r"Vapour", "Vapor"), (r"VAPOUR", "VAPOR"),
    (r"flavour", "flavor"), (r"Flavour", "Flavor"), (r"FLAVOUR", "FLAVOR"),
    (r"ionising", "ionizing"), (r"Ionising", "Ionizing"), (r"IONISING", "IONIZING"),
    (r"ionised", "ionized"), (r"Ionised", "Ionized"), (r"IONISED", "IONIZED"),
    (r"ionisation", "ionization"), (r"Ionisation", "Ionization"),
    (r"oxidised", "oxidized"), (r"Oxidised", "Oxidized"), (r"OXIDISED", "OXIDIZED"),
    (r"oxidising", "oxidizing"), (r"oxidisation", "oxidization"),
    (r"(?<![A-Za-z])disc(s?)(?![A-Za-z])", r"disk\1"), (r"(?<![A-Za-z])Disc(s?)(?![A-Za-z])", r"Disk\1"),
    (r"(?<![A-Za-z])DISC(S?)(?![A-Za-z])", r"DISK\1"),
    (r"link-neighbours", "link-neighbors"), (r"link_neighbours", "link_neighbors"),
)]
"""list: `(pattern, replacement)`, applied in order to every line that is not skipped."""


def _tracked():
    out = subprocess.check_output(["git", "ls-files"], cwd=ROOT, text=True)
    return [name for name in out.split("\n") if name]


def _in_scope(name):
    return not (name in SKIP_FILES or name.startswith(SKIP_PREFIX))


def _fix(text):
    """`text` with every British word of the rules changed, line by line (a line of `SKIP_LINES` stays)."""
    lines = text.split("\n")
    for number, line in enumerate(lines):
        if SKIP_LINES.search(line):
            continue
        for rule, replacement in RULES:
            line = rule.sub(replacement, line)
        lines[number] = line
    return "\n".join(lines)


def _read(name):
    try:
        with open(os.path.join(ROOT, name), encoding="utf-8") as handle:
            return handle.read()
    except (OSError, UnicodeDecodeError):
        return None


def _renamed(name):
    """The name of a file with its British words changed (`test_disc_x.py` becomes `test_disk_x.py`)."""
    folder, base = os.path.split(name)
    return os.path.join(folder, _fix(base))


def apply():
    """Rewrite the tree; returns the names of the files changed or renamed."""
    changed = []
    for name in _tracked():
        if not _in_scope(name):
            continue
        if name.endswith(EXTENSIONS):
            text = _read(name)
            if text is not None and _fix(text) != text:
                with open(os.path.join(ROOT, name), "w", encoding="utf-8") as handle:
                    handle.write(_fix(text))
                changed.append(name)
        new = _renamed(name)
        if new != name:
            subprocess.check_call(["git", "mv", name, new], cwd=ROOT)
            changed.append(f"{name} -> {new}")
    return changed


def left():
    """`[str]` the lines and file names that still hold a British word of the rules."""
    found = []
    for name in _tracked():
        if not _in_scope(name):
            continue
        if _renamed(name) != name:
            found.append(f"{name}: the file name")
        if name.endswith(EXTENSIONS):
            text = _read(name)
            for number, line in enumerate((text or "").split("\n"), 1):
                if not SKIP_LINES.search(line) and _fix(line) != line:
                    found.append(f"{name}:{number}: {line.strip()[:100]}")
    return found


def main(argv):
    if len(argv) != 2 or argv[1] not in ("apply", "check"):
        print(__doc__)
        return 2
    if argv[1] == "apply":
        changed = apply()
        print(f"{len(changed)} file(s) changed")
    remaining = left()
    for line in remaining:
        print(line)
    print("ok: no British word is left" if not remaining else f"{len(remaining)} place(s) still hold a British word")
    return 1 if remaining else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
