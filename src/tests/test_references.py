"""
docs/REFERENCES.md and docs/REFERENCES.html are written from one cache by scripts/build_references.py and must match it;
every listed paper carries a URL; and every author-year citation in the repository text has an outcome (DOC.95): tied to a paper in the
list, labelled as not a paper, or listed for review with what was tried.
"""

import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_references as builder  # noqa: E402

CACHE = json.loads((ROOT / "docs" / "reference" / "references-data.json").read_text(encoding="utf-8"))


def test_pages_match_the_cache():
    result = subprocess.run([sys.executable, str(ROOT / "scripts" / "build_references.py"), "render", "--check"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_every_entry_has_a_url():
    missing = [key for key, entry in CACHE["entries"].items() if builder.entry_urls(entry)[0] is None]
    assert missing == []


def test_the_opening_paragraph_thanks_the_researchers():
    text = (ROOT / "docs" / "REFERENCES.md").read_text(encoding="utf-8")
    assert "would not be possible without the dedicated work of the researchers" in text


def _outcome(key):
    status = CACHE.get("authoryear_status", {}).get(key, {})
    if status.get("matched"):
        return "tied"
    if status.get("not_a_paper"):
        return "not a paper"
    if status.get("review"):
        return "review"
    return None


def test_every_author_year_citation_has_an_outcome():
    """A new author-and-year citation written anywhere in the repository must be tied to a paper, labelled as not a paper, or put on the review list
    (run `build_references.py harvest`, then `settle_references.py triage local search`); it may not simply be unresolved."""
    *_, ays, _ctx, _more = builder.harvest([])
    open_keys = sorted("|".join(k) for k in ays if _outcome("|".join(k)) is None)
    assert open_keys == [], "%d citations with no outcome, first: %s" % (len(open_keys), open_keys[:10])


def test_every_outcome_is_complete():
    entries = CACHE["entries"]
    for key, status in CACHE.get("authoryear_status", {}).items():
        if status.get("matched"):
            for ref in status["matched"]:
                assert ref in entries, "%s is tied to %s, which is not in the list" % (key, ref)
                assert builder.entry_urls(entries[ref])[0] is not None, "%s: %s has no link" % (key, ref)
        elif status.get("not_a_paper"):
            assert status["not_a_paper"].get("what"), "%s is not a paper but says nothing about what it is" % key
        elif status.get("review"):
            assert status["review"].get("why"), "%s is for review but says nothing about what was tried" % key


def test_the_overview_shows_the_three_counts():
    text = (ROOT / "docs" / "REFERENCES.md").read_text(encoding="utf-8")
    groups = builder.citation_groups(CACHE)
    assert groups["open"] == []
    line = next(ln for ln in text.splitlines() if "author-year citations in the notes" in ln)
    for group in ("tied", "not", "review"):
        assert format(len(groups[group]), ",") in line
