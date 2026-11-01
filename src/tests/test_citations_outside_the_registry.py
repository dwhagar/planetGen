# tests/test_citations_outside_the_registry.py

"""
DOC.99: a citation is written only in the registry (`docs/reference/references-data.json`, `docs/REFERENCES.md` and the pages built from them).

Everywhere else a source is its citation ID in square brackets, `[R0123]` (`docs/design/interface/references.md`, "The in-line citation
style"; rule 17 of `docs/design/interface/text-style.md`). This test reads every tracked text file, runs the rewrite of
`scripts/rewrite_citations.py` over it and fails when the rewrite would change anything: an author and year, a DOI or an arXiv number
written out, and fails on a bibliography entry in APA form (`Surname, I. (2011). Title.`, `scripts/rewrite_bibliographies.py`). The failure names the file, the line and the ID to write instead. To cite a new source, `python scripts/build_references.py add
<doi or arXiv number>` lists it and prints its ID; to fix a place, `python scripts/rewrite_citations.py apply --all-candidates --strings
--scope <file>`.

The files that keep names are `rewrite_citations.KEEP` (the registry and its builders, the data license, the changelog history, the download
catalogue that needs its addresses). `docs/TODO.md`, `docs/plan/`, `docs/generated/` and `docs/archive/` follow the sweep through the TODO
thread and are checked once that has run (`TODO_OWNED`).
"""

import difflib
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import rewrite_bibliographies as rb  # noqa: E402
import rewrite_citations as rc  # noqa: E402

TODO_OWNED = ()
"""tuple: Text kept out of the lock until the TODO thread has swept it. Empty now: the TODO, plan, generated and archive text is swept (2026-11-01), so the lock covers all of it."""

SHOWN = 25
"""int: How many places a failure lists."""


def _places(name, reg):
    path = ROOT / name
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return []
    entries = ["%s:%d: %s\n      write: its ID, [R0123] (python scripts/rewrite_bibliographies.py apply ties an entry to the list)" % (name, n + 1, line.strip()[:110])
               for n, line in enumerate(text.split("\n")) if name.endswith(".md") and rb.APA.match(line)]
    new = rc.rewrite_file(name, text, reg, rc.Counter(), strings=True)
    if new == text:
        return entries
    old_lines, new_lines = text.split("\n"), new.split("\n")
    places = []
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, old_lines, new_lines, autojunk=False).get_opcodes():
        if tag != "equal":
            places.append("%s:%d: %s\n      write: %s" % (name, i1 + 1, old_lines[i1].strip()[:110], new_lines[j1].strip()[:110] if j1 < len(new_lines) else ""))
    return entries + places


def test_no_citation_is_written_outside_the_registry():
    reg = rc.Registry(all_candidates=True)
    found = []
    for name in rc.tracked([]):
        if not name.startswith(TODO_OWNED):
            found += _places(name, reg)
    assert not found, (
        "%d places write a source by name; write its ID, [R0123] (python scripts/build_references.py add <doi> lists a new source):\n  " % len(found)
        + "\n  ".join(found[:SHOWN])
    )
