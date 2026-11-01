# tests/test_citations_links.py

"""
DOC.97: a citation on a References page is a link, and only there.

`planetgen.web.lib.citations.link` turns `[R0123]` into the source's name as a link to its entry on the reference list;
`docs/references/citations.json` (built by `scripts/build_reference_pages.py`) holds the names, so a link cannot differ from the list.
A bracket that names no source, and an author and year written out in the text, fail instead of leaving a dead link.
"""

import json
import pathlib
import re
import sys

import pytest

from planetgen.web.lib import citations, reference_urls, references_site

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_references as refs  # noqa: E402

DATA = {
    "sources": {
        "R0001": {"kind": "paper", "label": "Moe et al. 2019", "title": "Mind your Ps & Qs"},
        "R0002": {"kind": "book", "label": "Salaris and Cassisi 2005", "title": "Evolution of Stars"},
        "R0010": {"kind": "data", "label": "Gaia DR3", "title": "Gaia DR3"},
    },
    "aliases": {"R0003": "R0001"},
}


def test_a_bracket_becomes_a_link_to_the_entry():
    out = citations.link("Binaries are common [R0001].", DATA)
    assert out == 'Binaries are common <a class="citation" href="/references/papers#paper-r0001" title="Mind your Ps &amp; Qs">Moe et al. 2019</a>.'


def test_several_sources_and_a_place_in_the_source():
    out = citations.link("See [R0001, R0002: Table 2].", DATA)
    assert out.count("<a ") == 2 and out.endswith("</a>, Table 2.")
    assert "#paper-r0001" in out and "#paper-r0002" in out


def test_a_dataset_links_to_the_data_page():
    assert 'href="/references/data#paper-r0010"' in citations.link("[R0010]", DATA)


def test_an_alias_follows_to_the_surviving_source():
    assert 'href="/references/papers#paper-r0001"' in citations.link("[R0003]", DATA)


def test_text_outside_the_brackets_is_escaped():
    assert citations.link("a < b [R0001]", DATA).startswith("a &lt; b ")


def test_an_unknown_source_fails():
    with pytest.raises(citations.CitationError, match="R0999"):
        citations.link("[R0999]", DATA)


def test_an_author_and_year_written_out_fails():
    for text in ("Moe et al. 2019 found it.", "Moe and Di Stefano (2017) found it.", "Moe, 2019"):
        with pytest.raises(citations.CitationError):
            citations.link(text, DATA)
    assert citations.link("Moe et al. 2019 [R0001]", DATA, strict=False)


def test_ids_in_lists_each_source_once_in_order():
    assert citations.ids_in("[R0002] then [R0001, R0002: Fig. 1] and [R0001]") == ["R0002", "R0001"]


def test_the_page_helper_reads_the_built_names():
    cache = json.loads((ROOT / "docs" / "reference" / "references-data.json").read_text(encoding="utf-8"))
    sid = next(iter(cache["ids"].values()))
    assert f"#paper-{sid.lower()}" in references_site.cite(f"[{sid}]")


def test_every_source_of_the_list_has_a_name_and_a_working_address():
    built = json.loads((ROOT / "docs" / "references" / "citations.json").read_text(encoding="utf-8"))
    cache = json.loads((ROOT / "docs" / "reference" / "references-data.json").read_text(encoding="utf-8"))
    listed = {sid for sid in cache["ids"].values() if sid not in cache["id_aliases"]}
    assert set(built["sources"]) == listed
    for sid, entry in built["sources"].items():
        assert entry["label"].strip() and entry["kind"] in ("paper", "book", "data")
        assert re.fullmatch(r"/references/(papers|data)#paper-r\d{4,}", reference_urls.reference_path("data" if entry["kind"] == "data" else "papers", sid))
    for old, new in built["aliases"].items():
        assert new in built["sources"]
