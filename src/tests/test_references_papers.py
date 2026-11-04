# tests/test_references_papers.py

"""
DOC.89 (References R3): the Papers and Data pages of the References section.

* The build (`scripts/build_reference_pages.py`): `papers.json` has a row for every source of the reference list, each with a link a reader can
  follow (a DOI page, else an arXiv page, else a search), a field read from the journal, and totals that are counted here a second way.
  `data.json` has a row for every data set, archive, catalog, and web source.
* Where we use each source (`web/lib/references_usage.py`): the area comes from the name of the file that cites it, a source cited in two areas
  is in both, and the registry rows are read as rows.
* The pages: `/references/papers` holds the anchor `paper-r0123` of every source exactly once, so a citation link always lands; the views by
  area, field, year, and type list a group; the search finds a source by title, author, DOI, or ID; the totals are the build's own.
  `/references/data` holds the anchor of every data source. The sentences keep the plain-words guide.
"""

import collections
import json
import pathlib
import re
import sys

import pytest

from planetgen.util import plainwords
from planetgen.web.lib import references_site, references_usage, reference_urls

from tests.test_class_pages import _html, app, client  # noqa: F401

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_reference_pages as builder  # noqa: E402

CACHE = json.loads((ROOT / "docs" / "reference" / "references-data.json").read_text(encoding="utf-8"))
PAPERS = references_site.load("papers.json")
DATA = references_site.load("data.json")
CITATIONS = references_site.load("citations.json")
ANCHOR = re.compile(r'id="(paper-r\d{4,})"')


def _usage(**uses):
    """A usage dict for a page test: `_usage(R0001=("stars",))` is a source cited in one place of the area stars."""
    return {"uses": {sid: {"values": [f"src/planetgen/physics/{area}.py" for area in areas], "notes": [], "stand_ins": [], "areas": list(areas)}
                     for sid, areas in uses.items()}, "files": 0}


# --- The build -------------------------------------------------------------------------------------------------------------------

def test_every_source_of_the_list_has_one_row_and_the_ids_are_the_lists_own():
    rows = PAPERS["rows"]
    assert len(rows) == len(CACHE["entries"]) > 4000
    ids = [row["id"] for row in rows]
    assert len(set(ids)) == len(ids)
    assert all(reference_urls.is_paper_id(sid) for sid in ids)
    assert set(ids) == {CACHE["ids"][key.lower()] for key in CACHE["entries"]}


def test_every_row_has_a_link_a_reader_can_follow_and_the_first_is_the_best_one():
    for row in PAPERS["rows"]:
        links = references_site.paper_links(row)
        assert links and links[0]["url"].startswith(("https://", "http://")), row["id"]
        if row.get("doi"):
            assert links[0]["kind"] == "doi" and links[0]["url"].startswith("https://doi.org/"), row["id"]
        elif row.get("arxiv"):
            assert links[0]["kind"] == "arxiv" and links[0]["url"].startswith("https://arxiv.org/abs/"), row["id"]
        assert len({link["url"] for link in links}) == len(links), row["id"]


def test_the_links_come_in_the_order_doi_arxiv_ads_and_a_search_is_the_last_resort():
    both = references_site.paper_links({"doi": "10.1/ab(c)", "arxiv": "1234.5678", "title": "T", "who": "Smith, J."})
    assert [link["kind"] for link in both] == ["doi", "arxiv", "ads"]
    assert both[0]["url"] == "https://doi.org/10.1/ab(c)" and both[2]["url"].endswith("10.1/ab(c)/abstract")
    only_arxiv = references_site.paper_links({"arxiv": "astro-ph/0508219"})
    assert [link["kind"] for link in only_arxiv] == ["arxiv", "ads"] and "arXiv:astro-ph/0508219" in only_arxiv[1]["url"]
    own = references_site.paper_links({"link": "https://example.org/a"})
    assert [link["kind"] for link in own] == ["link"]
    search = references_site.paper_links({"title": "Some Title", "who": "Smith, J. & Doe, A."})
    assert [link["kind"] for link in search] == ["search"] and "Some%20Title%20Smith" in search[0]["url"]
    assert references_site.paper_links({})[0]["kind"] == "search"


@pytest.mark.parametrize("journal, field", [
    ("Icarus", "Planetary science"), ("MNRAS", "Astronomy and astrophysics"), ("GeCoA", "Earth and geoscience"),
    ("JAtS", "Atmosphere and climate"), ("AsBio", "Life sciences"), ("JPCA", "Chemistry"), ("PhRvL", "Physics"),
    ("Math. Comp.", "Mathematics and computing"), ("Nature", "General science"), ("Arch Intern Med (Chic)", "Medicine and health"),
    ("The American Historical Review", "Other"),
])
def test_the_field_is_read_from_the_journal(journal, field):
    assert builder.field_of({"journal": journal}) == field


def test_a_source_with_no_journal_has_no_field_unless_its_arxiv_id_says_astronomy():
    assert builder.field_of({}) == "Not stated"
    assert builder.field_of({"arxiv": "astro-ph/0508219"}) == "Astronomy and astrophysics"
    assert builder.field_of({"book": "Handbook of Exoplanets"}) == "Planetary science"


def test_the_fields_are_the_listed_ones():
    names = {name for name, _pattern in builder.FIELDS} | {builder.NO_JOURNAL, builder.OTHER_FIELD}
    assert {row["field"] for row in PAPERS["rows"]} <= names
    assert collections.Counter(row["field"] for row in PAPERS["rows"])["Astronomy and astrophysics"] > 1500


def test_the_totals_are_the_lists_own_counted_a_second_way():
    totals, entries = PAPERS["totals"], list(CACHE["entries"].values())
    assert totals["papers"] == len(entries)
    assert totals["with_doi"] == sum(1 for e in entries if e.get("doi"))
    assert totals["with_arxiv"] == sum(1 for e in entries if e.get("arxiv"))
    assert totals["datasets"] == sum(1 for e in entries if e.get("type") == "dataset")
    assert totals["books"] == sum(1 for e in entries if e.get("type") in ("book", "monograph", "edited-book"))
    years = [int(e["year"]) for e in entries if e.get("year")]
    assert totals["oldest"]["year"] == min(years) and totals["newest"]["year"] == max(years)
    for key in ("by_field", "by_decade", "by_type"):
        assert sum(number for _name, number in totals[key]) == totals["papers"], key
    assert [name for name, _n in totals["by_decade"]] == sorted(name for name, _n in totals["by_decade"])


def test_the_decade_of_a_year():
    assert builder.decade_of(1994) == "1990s" and builder.decade_of(2000) == "2000s" and builder.decade_of(0) == "Not stated"


def test_the_data_file_has_a_row_for_each_source_of_the_data_page():
    assert {row["id"] for row in DATA["setup"]} == {CACHE["ids"]["setup:" + d["url"].lower()] for d in CACHE["setup_datasets"]}
    assert {row["id"] for row in DATA["archives"]} == {CACHE["ids"]["archive:" + a["url"].lower()] for a in CACHE["archives"]}
    assert {row["id"] for row in DATA["web"]} == {CACHE["ids"]["web:" + w["url"].lower()] for w in CACHE["web_sources"]}
    assert all(row["link"].startswith(("https://", "http://")) for key in ("setup", "archives", "web") for row in DATA[key])


def test_a_title_has_no_markup_and_no_character_references():
    assert builder.clean("A &amp; B <i>in</i>   space") == "A & B in space"
    assert builder.clean("a &amp;lt; 200 au and 0.3 <z< 4") == "a < 200 au and 0.3 <z< 4"
    assert not [row["id"] for row in PAPERS["rows"] if re.search(r"</?(i|b|sub|sup|em|mml:\w+|span)\b|&[a-z]+;", row.get("title", "") + row.get("venue", ""))]


def test_the_rows_are_written_one_to_a_line_and_the_build_is_stable():
    text = builder.render({"a": 1, "rows": [{"x": 1}, {"x": 2}]}, "rows")
    assert text == '{\n  "a": 1,\n  "rows": [\n    {"x": 1},\n    {"x": 2}\n  ]\n}\n'
    assert json.loads(text) == {"a": 1, "rows": [{"x": 1}, {"x": 2}]}
    built = builder.build()
    assert built["papers.json"].count("\n") == len(PAPERS["rows"]) + 12 or built["papers.json"].count("\n") > len(PAPERS["rows"])
    assert builder.build() == built


# --- Where we use each source ------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("path, area", [
    ("src/planetgen/physics/ozone.py", "atmospheres"), ("src/planetgen/galaxy/potential.py", "galaxy"),
    ("src/planetgen/physics/tidal_heating.py", "orbits"), ("docs/design/habitability/life-classes.md", "life"),
    ("src/planetgen/physics/rings.py", "rings-magnetism"), ("src/planetgen/physics/solar_particles.py", "radiation"),
    ("src/planetgen/physics/belt_mass.py", "moons-belts"), ("src/planetgen/physics/saturation.py", "climate"),
    ("src/planetgen/physics/body_colour.py", "appearance"), ("src/planetgen/physics/main_sequence.py", "stars"),
    ("src/planetgen/physics/bulk_structure.py", "bodies"), ("docs/design/README.md", "other"),
])
def test_the_area_is_read_from_the_name_of_the_file(path, area):
    assert references_usage.area_of(path) == area


def test_every_area_of_the_rules_is_an_area_of_the_method_pages_and_has_a_name():
    assert {area for area, _pattern in references_usage.RULES} <= set(reference_urls.AREAS)
    assert set(references_usage.AREA_TITLES) == set(reference_urls.AREAS) | {references_usage.OTHER}


def test_the_registry_rows_are_read_as_rows_with_their_section():
    text = ("## 5.3 Atmosphere\n\n| Id | Kind | Fill |\n|---|---|---|\n| AT.a1, a2 | D | 1 the law [R0001, R0002: Table 2] |\n"
            "| BC1 | D | nothing cited |\n| GR6 and GR3 thick disc | x | [R0003] |\n\n## 5.1 Galaxy laws\n\n| Id | Fill |\n|---|---|\n| GR1 | [R0003] |\n")
    assert references_usage.registry_rows(text) == [("AT.a1, a2", "5.3 Atmosphere", ["R0001", "R0002"]), ("GR1", "5.1 Galaxy laws", ["R0003"])]


def test_a_scan_finds_the_places_of_a_source_by_kind_and_area(tmp_path):
    def write(relative, text):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    write("src/planetgen/physics/ozone.py", "# Cross sections [R0001: Table 1] and [R0002].\n")
    write("src/planetgen/physics/data/gases.json", '{"source": "[R0001]"}\n')
    write("src/planetgen/web/lib/citations.py", "# an example [R0001]\n")
    write("src/planetgen/web/static/x.js", "// [R0001]\n")
    write("docs/design/habitability/life-classes.md", "Life [R0001, R0003].\n")
    write("docs/design/research/gap-filling-plan.md", "## 5.1 Galaxy laws\n\n| Id | Fill |\n|---|---|\n| GR1 | [R0001] |\n")
    found = references_usage.scan(str(tmp_path))
    assert found["files"] == 3
    one = found["uses"]["R0001"]
    assert one["values"] == ["src/planetgen/physics/data/gases.json", "src/planetgen/physics/ozone.py"]
    assert one["notes"] == ["docs/design/habitability/life-classes.md"] and one["stand_ins"] == ["GR1"]
    assert one["areas"] == ["galaxy", "atmospheres", "life"] or sorted(one["areas"]) == ["atmospheres", "galaxy", "life"]
    assert found["uses"]["R0003"]["areas"] == ["life"] and references_usage.cited_by(one) == 4


def test_the_real_scan_reads_the_project_and_cites_most_sources_somewhere():
    found = references_usage.uses()
    assert found["files"] > 500 and len(found["uses"]) > 1500
    ids = {row["id"] for row in PAPERS["rows"]}
    assert not [sid for sid in found["uses"] if sid not in ids and sid not in CITATIONS["sources"]][:5]
    assert all(entry["areas"] for entry in found["uses"].values())


# --- The Papers page ----------------------------------------------------------------------------------------------------------------

def test_the_papers_page_holds_every_anchor_exactly_once_so_a_citation_always_lands(client):
    html = _html(client, "/references/papers")
    anchors = ANCHOR.findall(html)
    assert len(anchors) == len(set(anchors)) == len(PAPERS["rows"])
    assert 'id="out-of-universe-banner"' in html and "<title>Papers" in html
    papers = {sid for sid, source in CITATIONS["sources"].items() if source["kind"] in ("paper", "book")}
    assert {reference_urls.paper_anchor(sid) for sid in papers} <= set(anchors)


def test_the_data_page_holds_every_anchor_of_a_data_or_web_source(client):
    html = _html(client, "/references/data")
    anchors = ANCHOR.findall(html)
    assert len(anchors) == len(set(anchors))
    data = {sid for sid, source in CITATIONS["sources"].items() if source["kind"] in ("data", "web")}
    assert data and {reference_urls.paper_anchor(sid) for sid in data} == set(anchors)


def test_every_citation_link_names_a_page_that_holds_its_anchor(client):
    pages = {"papers": set(ANCHOR.findall(_html(client, "/references/papers"))), "data": set(ANCHOR.findall(_html(client, "/references/data")))}
    for sid, source in CITATIONS["sources"].items():
        page = "data" if source["kind"] in ("data", "web") else "papers"
        assert reference_urls.paper_anchor(sid) in pages[page], sid


def test_the_a_to_z_list_has_a_letter_index_and_the_links_of_each_source(client):
    html = _html(client, "/references/papers")
    assert 'href="#letter-a"' in html and 'id="letter-a"' in html
    assert html.count('class="reference-entry reference-compact"') == len(PAPERS["rows"])
    assert "https://doi.org/10.1093/mnras/276.1.268" in html
    assert 'aria-current="page">A to Z</a>' in html


def test_a_view_by_area_lists_the_groups_and_a_group_lists_its_sources(client):
    html = _html(client, "/references/papers?view=area")
    for title in references_usage.AREA_TITLES.values():
        assert f">{title}</a>" in html, title
    group = _html(client, "/references/papers?view=area&group=stars")
    assert "<h2 id=\"group-stars-heading\">Stars</h2>" in group
    ids = ANCHOR.findall(group)
    assert ids and len(ids) == len(set(ids))
    stars = {sid for sid, entry in references_usage.uses()["uses"].items() if "stars" in entry["areas"]}
    assert {reference_urls.paper_anchor(sid) for sid in stars if sid in {row["id"] for row in PAPERS["rows"]}} == set(ids)


def test_a_source_cited_in_two_areas_is_in_both_and_has_one_anchor_in_a_view():
    row = PAPERS["rows"][0]
    usage = _usage(**{row["id"]: ("stars", "climate")})
    groups = references_site.paper_groups([row], usage, "area")
    assert sorted(groups) == ["climate", "stars"]
    assert references_site.paper_groups([row], _usage(), "area")["none"][0] == references_site.UNUSED_TITLE
    page = references_site.papers_page({"rows": [row, row], "totals": PAPERS["totals"]}, usage, view="area", group="stars")
    assert [bool(entry["anchor"]) for entry in page["blocks"][0]["entries"]] == [True, False]


@pytest.mark.parametrize("view, group", [("field", "planetary-science"), ("year", "1990s"), ("type", "book")])
def test_the_other_views_list_a_group(client, view, group):
    html = _html(client, f"/references/papers?view={view}&group={group}")
    ids = ANCHOR.findall(html)
    assert ids and len(ids) == len(set(ids))
    expected = {
        "field": lambda row: row["field"] == "Planetary science",
        "year": lambda row: row.get("year", 0) // 10 * 10 == 1990,
        "type": lambda row: row["type"] == "Book",
    }[view]
    assert len(ids) == sum(1 for row in PAPERS["rows"] if expected(row))


def test_a_group_that_does_not_exist_shows_the_groups_and_a_view_that_does_not_exist_shows_the_a_to_z_list(client):
    assert "Choose a group to list its sources." in _html(client, "/references/papers?view=field&group=nothing")
    assert len(ANCHOR.findall(_html(client, "/references/papers?view=nothing"))) == len(PAPERS["rows"])


def test_the_totals_view_shows_the_builds_numbers_and_the_sources_we_cite_most(client):
    html = _html(client, "/references/papers?view=totals")
    totals = PAPERS["totals"]
    for number in (totals["papers"], totals["with_doi"], totals["with_arxiv"]):
        assert f"{number:,}" in html
    for heading in ("The reference list in numbers", "Most cited by us", "Where we use them", "The field of the journal", "The decade of the paper"):
        assert f">{heading}</h2>" in html, heading
    page = references_site.papers_page(PAPERS, references_usage.uses(), view="totals")
    most = next(section for section in page["sections"] if section["id"] == "most-cited")["tables"][0]["rows"]
    counts = [int(row["cells"][2].replace(",", "")) for row in most]
    assert len(counts) == references_site.MOST_CITED and counts == sorted(counts, reverse=True)
    by_area = next(section for section in page["sections"] if section["id"] == "by-area")["tables"][0]["rows"]
    assert {row["cells"][0]["text"] for row in by_area} >= {"Stars", "Atmospheres", "Life"}


@pytest.mark.parametrize("query, ident", [("abundances at high redshift", "R0001"), ("R0006", "R0006"), ("10.1093/mnras/276.1.268", "R0001"),
                                          ("Abbot Switzer steppenwolf", "R0010"), ("2508.18083", "R0006")])
def test_the_search_finds_a_source_by_title_author_doi_or_id(client, query, ident):
    html = _html(client, "/references/papers?q=" + query.replace(" ", "+"))
    assert f'id="{reference_urls.paper_anchor(ident)}"' in html and "sources match" in html or "1 sources match" in html


def test_a_search_with_no_match_says_so_and_a_broad_one_is_cut(client):
    assert "No source matches." in _html(client, "/references/papers?q=zzzzzzqqqq")
    page = references_site.papers_page(PAPERS, _usage(), q="a")
    assert len(page["blocks"][0]["entries"]) == references_site.SEARCH_LIMIT and f"The first {references_site.SEARCH_LIMIT} are shown." in page["note"]


def test_a_full_entry_lists_the_places_that_cite_it(client):
    sid = next(sid for sid, entry in references_usage.uses()["uses"].items() if entry["values"] and entry["stand_ins"])
    html = _html(client, f"/references/papers?q={sid}")
    assert "Code and data:" in html and "Stand-ins:" in html and "Used in:" in html
    entry = references_site.paper_entry(next(row for row in PAPERS["rows"] if row["id"] == sid), references_usage.uses())
    assert entry["cited_by"] == references_usage.cited_by(references_usage.uses()["uses"][sid]) >= 2


def test_a_source_that_nothing_cites_says_so():
    row = PAPERS["rows"][0]
    entry = references_site.paper_entry(row, _usage())
    assert entry["areas"] == [] and entry["cited_by"] == 0


def test_more_places_than_fit_are_counted_not_listed():
    row = PAPERS["rows"][0]
    usage = {"uses": {row["id"]: {"values": [f"src/planetgen/physics/m{n}.py" for n in range(12)], "notes": [], "stand_ins": [], "areas": ["bodies"]}}}
    entry = references_site.paper_entry(row, usage)
    assert len(entry["code"]["names"]) == references_site.SHOWN and entry["code"]["more"] == 4
    assert entry["code"]["names"][0] == "physics/m0.py"


def test_the_entry_text_is_escaped(client, monkeypatch):
    row = dict(PAPERS["rows"][0], title="<script>alert(1)</script> & more", who="A <b>", venue="V\"q")
    monkeypatch.setattr(references_site, "load", lambda name, folder=None: {"rows": [row], "totals": PAPERS["totals"]} if name == "papers.json" else None)
    html = _html(client, "/references/papers")
    assert "<script>alert" not in html and "&lt;script&gt;alert(1)&lt;/script&gt; &amp; more" in html


def test_a_site_that_has_not_built_the_papers_or_the_data_says_so_and_keeps_the_banner(client, monkeypatch):
    monkeypatch.setattr(references_site, "load", lambda name, folder=None: None)
    for path in ("/references/papers", "/references/data"):
        html = _html(client, path)
        assert 'id="out-of-universe-banner"' in html and "The pages are not built" in html and "build_reference_pages.py" in html


def test_the_references_pages_get_no_help_links_in_their_text(client):
    assert 'class="help-mark"' not in _html(client, "/references/papers?q=planet") and 'class="help-mark"' not in _html(client, "/references/data")


# --- The Data page ----------------------------------------------------------------------------------------------------------------

def test_the_licence_tables_are_read_from_the_file_and_each_row_has_a_cell_for_each_column():
    tables = references_site.licence_tables(references_site.licence_text())
    assert len(tables) > 10 and sum(len(table["rows"]) for table in tables) > 100
    assert all(table["head"] and all(len(row) == len(table["head"]) for row in table["rows"]) for table in tables)
    assert not [cell for table in tables for row in table["rows"] for cell in row if "`" in cell]


def test_the_licence_reader_handles_a_small_text_and_a_missing_file(tmp_path):
    text = "# Title\n\n## Part `x`\n\n| A | B |\n|---|---|\n| `one` | [two](https://e.org) |\n| short |\n\nText.\n"
    assert references_site.licence_tables(text) == [{"heading": "Part x", "head": ["A", "B"], "rows": [["one", "two"], ["short", ""]]}]
    assert references_site.licence_text(str(tmp_path / "nothing.md")) == "" and references_site.licence_tables("") == []


def test_the_terms_of_the_setup_data_sets_are_matched_by_key(client):
    terms = references_site._setup_terms(references_site.licence_tables(references_site.licence_text()))
    assert {"isochrones", "ozone_uv"} <= set(terms)
    page = references_site.data_page(DATA, references_site.licence_text())
    setup = next(section for section in page["sections"] if section["id"] == "setup")["tables"][0]["rows"]
    by_name = {row["cells"][1]["text"]: row["cells"][3] for row in setup}
    assert any(text for text in by_name.values()) and any("CC BY" in text or "Open" in text for text in by_name.values())


def test_the_data_page_has_its_sections_and_the_terms_of_use(client):
    html = _html(client, "/references/data")
    for heading in ("Data sets downloaded at setup", "Archives, catalogs, and tools", "Standards and web pages", "Data sets that have a DOI", "Terms of use",
                    "With thanks"):
        assert f">{heading}</h2>" in html, heading
    assert 'id="out-of-universe-banner"' in html and "NIST Chemistry WebBook" in html
    assert 'href="/references/papers?view=type&amp;group=data-set"' in html


def test_the_overview_links_to_the_papers_and_the_data_page(client):
    html = _html(client, "/references")
    assert 'href="/references/papers"' in html and 'href="/references/data"' in html


# --- The words --------------------------------------------------------------------------------------------------------------------

def _prose(page):
    texts = list(page.get("intro", [])) + [p for section in page["sections"] for p in section["paragraphs"]]
    heads = [section["heading"] for section in page["sections"]] + [head for section in page["sections"] for table in section["tables"] for head in table["head"]]
    return texts, heads


def test_the_sentences_of_the_papers_and_the_data_pages_keep_the_guide():
    pages = [references_site.papers_page(PAPERS, references_usage.uses(), view=view) for view in ("az", "area", "totals")]
    pages.append(references_site.data_page(DATA, ""))
    pages.append(references_site.data_page(DATA, references_site.licence_text()))
    for page in pages:
        texts, heads = _prose(page)
        broken = [(sentence, rule) for text in texts for sentence in plainwords.sentences(text)
                  for rule, _detail in plainwords.check_sentence(sentence, references=True)]
        assert not broken, broken
        bad = [(label, plainwords.check_label(label)) for label in heads if plainwords.check_label(label)]
        assert not bad, bad
    assert all(page["title"] and page["views"] is not None for page in pages)
