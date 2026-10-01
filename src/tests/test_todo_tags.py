"""TODO(...) tags in the code name an item in docs/TODO.md.

Items in docs/TODO.md have permanent category IDs (for example UX.2 or
MAP.2.1), each written as a "**ID " heading. A tag in the code must use
one of those IDs, as TODO(CAT.N) or TODO(CAT.N.M...), so it keeps pointing
at the right item. The old running numbers ("TODO(installers #50)") were
renumbered with the file and no longer mean anything on their own.

docs/TODO.md is read when the test runs, so the check follows the file as
it changes. With no tags in the code, there is nothing to check.
"""

import os
import re

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

CATEGORIES = (
    "UX", "MAP", "NAV", "GEN", "PERF", "DB", "API",
    "ADM", "SEC", "USR", "OPS", "DOC", "VIEW", "POP",
)

# Where tags are checked: these directories (recursively) and top-level files.
SCAN_DIRS = ("src", "scripts", "examples")
SCAN_FILES = ("generate.py", "install.sh", "update.sh", "install.ps1", "update.ps1")

# This file and the installer test name tags on purpose (the latter checks
# that the old "TODO(installers #50)" tag stays gone).
SKIP_FILES = {"test_todo_tags.py", "test_install_python_deps.py"}
SKIP_DIRS = {"__pycache__", "vendor", "node_modules", ".git"}

TAG_RE = re.compile(r"\bTODO\(([^)]*)\)")
ID_RE = re.compile(r"(?:%s)\.\d+(?:\.\d+)*" % "|".join(CATEGORIES))
HEADING_RE = re.compile(r"\*\*((?:%s)\.\d+(?:\.\d+)*) " % "|".join(CATEGORIES))


def _source_files():
    for name in SCAN_FILES:
        path = os.path.join(ROOT, name)
        if os.path.isfile(path):
            yield path
    for top in SCAN_DIRS:
        for dirpath, dirnames, filenames in os.walk(os.path.join(ROOT, top)):
            dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
            for name in filenames:
                if name in SKIP_FILES or name.endswith((".pyc", ".min.js")):
                    continue
                yield os.path.join(dirpath, name)


def _tags():
    """Every TODO(...) tag in the scanned files, as (path:line, inside)."""
    found = []
    for path in _source_files():
        try:
            with open(path, encoding="utf-8") as f:
                lines = f.readlines()
        except (UnicodeDecodeError, OSError):
            continue
        for number, line in enumerate(lines, 1):
            for match in TAG_RE.finditer(line):
                where = "%s:%d" % (os.path.relpath(path, ROOT), number)
                found.append((where, match.group(1)))
    return found


def _todo_ids():
    path = os.path.join(ROOT, "docs", "TODO.md")
    with open(path, encoding="utf-8") as f:
        return set(HEADING_RE.findall(f.read()))


def test_tags_use_category_ids():
    bad = [(where, inside) for where, inside in _tags() if not ID_RE.fullmatch(inside)]
    assert not bad, "TODO tags must be TODO(CAT.N[.M...]): %r" % bad


def test_tags_name_items_in_todo_md():
    tags = [(where, inside) for where, inside in _tags() if ID_RE.fullmatch(inside)]
    if not tags:
        return
    ids = _todo_ids()
    missing = [(where, inside) for where, inside in tags if inside not in ids]
    assert not missing, "TODO tags with no '**ID ' heading in docs/TODO.md: %r" % missing


def test_tag_patterns():
    assert ID_RE.fullmatch("MAP.2.1")
    assert ID_RE.fullmatch("PERF.2")
    assert not ID_RE.fullmatch("installers #50")
    assert not ID_RE.fullmatch("FOO.1")
    assert not ID_RE.fullmatch("MAP")
    assert HEADING_RE.findall("- **MAP.2.1 Drill-down stages** text") == ["MAP.2.1"]
