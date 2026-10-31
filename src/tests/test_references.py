"""
docs/REFERENCES.md and docs/REFERENCES.html are written from one cache by scripts/build_references.py and must match it;
every listed paper carries a URL; and every author-year citation in the repository text has an outcome (DOC.95): tied to a paper in the
list, labelled as not a paper, or listed for review with what was tried.
"""

import json
import pathlib
import subprocess
import sys

import pytest

from planetgen.web.lib import reference_urls

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


# --- citation IDs (DOC.100) -------------------------------------------------------------------------------------------------------

def _tiny(*keys):
    return {"entries": {k: {"authors": [["Smith", "A."]], "year": 2000, "title": "Paper " + k, "doi": k if k.startswith("10.") else ""} for k in keys}}


def test_every_source_has_one_id_and_no_two_share_it():
    assert builder.assign_ids(json.loads(json.dumps(CACHE))) == 0, "a source has no ID: run build_references.py render"
    ids = CACHE["ids"]
    for key, entry in CACHE["entries"].items():
        if entry.get("authors") or entry.get("title"):
            assert builder.id_of(CACHE, key).startswith("R")
    by_source = {}
    for name, sid in ids.items():
        by_source.setdefault(sid, set())
    assert len(by_source) == len(set(ids.values()))
    assert all(reference_urls.is_paper_id(sid) for sid in ids.values())
    shown = [builder.id_of(CACHE, k) for k, e in CACHE["entries"].items() if e.get("authors") or e.get("title")]
    assert len(shown) == len(set(shown)), "two sources share an ID"


def test_ids_are_numbered_once_and_never_change():
    cache = _tiny("10.1/a", "10.1/b")
    assert builder.assign_ids(cache) == 2
    first = dict(cache["ids"])
    assert builder.assign_ids(cache) == 0
    assert cache["ids"] == first
    cache["entries"].update(_tiny("10.1/c")["entries"])
    assert builder.assign_ids(cache) == 1
    assert all(cache["ids"][k] == v for k, v in first.items())


def test_an_id_is_never_reused_after_its_source_is_removed():
    cache = _tiny("10.1/a", "10.1/b")
    builder.assign_ids(cache)
    gone = builder.id_of(cache, "10.1/b")
    del cache["entries"]["10.1/b"]
    cache["entries"].update(_tiny("10.1/z")["entries"])
    builder.assign_ids(cache)
    assert builder.id_of(cache, "10.1/z") != gone
    assert builder.resolve_id(cache, gone) == gone


def test_a_preprint_that_gets_a_doi_keeps_its_id():
    cache = {"entries": {"arXiv:1234.5678": {"authors": [["Smith", "A."]], "year": 2000, "title": "P", "arxiv": "1234.5678", "doi": ""}}}
    builder.assign_ids(cache)
    before = builder.id_of(cache, "arXiv:1234.5678")
    del cache["entries"]["arXiv:1234.5678"]
    cache["entries"]["10.1/p"] = {"authors": [["Smith", "A."]], "year": 2000, "title": "P", "arxiv": "1234.5678", "doi": "10.1/p"}
    assert builder.assign_ids(cache) == 0
    assert builder.id_of(cache, "10.1/p") == before


def test_two_ids_that_turn_out_to_be_one_paper_keep_the_older_as_the_id_and_the_newer_as_an_alias():
    cache = {"entries": {"arXiv:1.1": {"authors": [["S", "A."]], "year": 1, "title": "P", "arxiv": "1.1", "doi": ""},
                         "10.1/p": {"authors": [["S", "A."]], "year": 1, "title": "P", "arxiv": "", "doi": "10.1/p"}}}
    builder.assign_ids(cache)
    older, newer = sorted((builder.id_of(cache, "arXiv:1.1"), builder.id_of(cache, "10.1/p")))
    cache["entries"]["10.1/p"]["arxiv"] = "1.1"
    del cache["entries"]["arXiv:1.1"]
    builder.assign_ids(cache)
    assert builder.id_of(cache, "10.1/p") == older
    assert builder.resolve_id(cache, newer) == older


def test_the_reference_list_prints_each_id_and_its_papers_page_anchor():
    text = (ROOT / "docs" / "REFERENCES.md").read_text(encoding="utf-8")
    for sid in sorted(set(CACHE["ids"].values()))[:50]:
        assert 'id="paper-%s"' % sid.lower() in text and "**%s**" % sid in text


def test_the_id_helpers_of_the_site_agree_with_the_builder():
    assert reference_urls.paper_anchor("R0123") == "paper-r0123"
    assert builder.format_id(7) == "R0007" and builder.format_id(12345) == "R12345"
    with pytest.raises(ValueError):
        reference_urls.paper_anchor("10.1086/151682")
