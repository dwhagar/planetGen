# tests/test_references_stand_ins.py

"""
DOC.90 (References R4): the stand-in registry pages.

* The registry is read whole (`scripts/build_stand_in_values.py`): two tables of the plan carry a section number another table already has
  (6.127 and 6.134); each is read once, under its own key, so no row is read twice and none is left out.
* The build (`scripts/build_reference_pages.py`): `docs/references/stand-ins.json` holds every row once, with a page address no other row has,
  and `fitted.json` holds the tuning targets.
* The pages (`web/references_pages.py`, `web/lib/references_site.py`): `/references/stand-ins` filters the rows, each row has its own page with
  its citations linked, and `/references/fitted` lists the rows at tweak levels L2 to L4 and the targets.
"""

import json
import pathlib
import re
import sys

import pytest

from planetgen.util import plainwords
from planetgen.web.lib import reference_urls, references_site

from tests.test_class_pages import _html, app, client  # noqa: F401

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_reference_pages as builder  # noqa: E402
import build_stand_in_values as registry  # noqa: E402

PLAN = ROOT / "docs" / "design" / "research" / "gap-filling-plan.md"
BANNER_ID = 'id="out-of-universe-banner"'


@pytest.fixture(scope="module")
def data():
    return json.loads((ROOT / "docs" / "references" / "stand-ins.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def overview():
    return json.loads((ROOT / "docs" / "references" / "overview.json").read_text(encoding="utf-8"))


# --- The registry is read whole --------------------------------------------------------------------------------------------------

def test_two_tables_with_one_section_number_are_each_read_once_under_their_own_title():
    tables = registry.registry_tables(PLAN.read_text(encoding="utf-8"))
    assert len(tables) == len(registry.SECTIONS) == len({key for key, *_rest in registry.sections()})
    giant, ocean = (tables[key][1] for key in ("6.134 Stand-ins of the structure of a giant",
                                               "6.134 Stand-ins of the ocean cover law and the mantle's water store"))
    assert {"GS2", "GS8"} <= {row[0] for row in giant} and {"HY1", "HY2"} <= {row[0] for row in ocean}
    assert not {row[0] for row in giant} & {row[0] for row in ocean}
    high_gravity = tables["6.127 Stand-ins of the high-gravity class line"][1]
    assert "GV1" in {row[0] for row in high_gravity}
    assert "EJ1" in {row[0] for row in tables["6.127 Stand-ins of the architecture's losses"][1]}


def test_every_row_of_the_plan_is_in_the_list_once():
    tables = registry.registry_tables(PLAN.read_text(encoding="utf-8"))
    keys = [f"{number}:{row[0]}" for key, number, _title, _columns in registry.sections() for row in tables[key][1]]
    assert len(keys) == len(set(keys)) == registry._counts(tables)[0]


# --- The build -------------------------------------------------------------------------------------------------------------------

def test_the_data_holds_every_row_once_and_no_two_rows_share_an_address(data):
    rows = references_site.rows_of(data)
    total = builder.stand_in_counts(PLAN.read_text(encoding="utf-8"))["rows"]
    assert len(rows) == total > 1000
    assert len({row["key"] for _part, row in rows}) == len({row["slug"] for _part, row in rows}) == total
    assert all(row["slug"] == reference_urls.anchor(row["key"]) for _part, row in rows)


def test_every_row_has_a_tier_a_tweak_level_and_something_to_say(data, overview):
    tiers, tweaks = {row[0] for row in overview["tiers"]}, {row[0] for row in overview["tweaks"]}
    for part, row in references_site.rows_of(data):
        assert row["tier_code"] in tiers and row["tweak"] in tweaks, row["key"]
        assert row["what"] and row["id"], row["key"]
        assert row["fields"] or row["kind"], row["key"]
        assert row["kind"] in data["kinds"], (row["key"], row["kind"])


def test_the_counts_by_tier_and_tweak_match_the_overview(data, overview):
    rows = [row for _part, row in references_site.rows_of(data)]
    for code, number in overview["stand_ins"]["by_tier"].items():
        assert number == sum(1 for row in rows if row["tier_code"] == code), code
    for code, number in overview["stand_ins"]["by_tweak"].items():
        assert number == sum(1 for row in rows if row["tweak"] == code), code


def test_the_unread_and_burn_down_lines_come_from_their_registers(data):
    by_key = {row["key"]: row for _part, row in references_site.rows_of(data)}
    assert by_key["6.7:AM64"]["burn_down"][:2] == ["replace", "GEN.804"]
    assert by_key["5.3:AT.b5"]["burn_down"][0] == "no data"
    assert any(row["unread"] for row in by_key.values())


def test_a_bracket_that_names_no_source_stops_the_build_and_names_the_row(monkeypatch):
    real = builder.stand_ins.registry_tables

    def broken(text):
        tables = real(text)
        header, rows = tables["6.7"]
        rows[0][header.index("evidence")] = "read in [R9999]"
        return tables

    monkeypatch.setattr(builder.stand_ins, "registry_tables", broken)
    with pytest.raises(ValueError, match=r"stand-in 6\.7:AM1: R9999 is not a source of the reference list"):
        builder.build_stand_ins()


def test_the_tuning_targets_are_all_in_the_fitted_data():
    fitted = json.loads((ROOT / "docs" / "references" / "fitted.json").read_text(encoding="utf-8"))
    table = json.loads((ROOT / "src" / "planetgen" / "physics" / "data" / "tuning_targets.json").read_text(encoding="utf-8"))
    assert [t["id"] for t in fitted["targets"]] == [t["id"] for t in table["targets"]]
    assert all(t["statistic"] and t["label"] for t in fitted["targets"])


def test_the_committed_files_are_what_the_registries_give_now():
    for name, build in (("stand-ins.json", builder.build_stand_ins), ("fitted.json", builder.build_fitted)):
        assert (ROOT / "docs" / "references" / name).read_text(encoding="utf-8") == builder.render(build()), name


# --- The address of a stand-in ---------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("key, path", [
    ("6.7:AM1", "/references/stand-ins/6-7-am1"), ("5.3:AT.a1, a2", "/references/stand-ins/5-3-at-a1-a2"),
    ("6.3a:--shape", "/references/stand-ins/6-3a-shape"),
])
def test_the_key_of_a_row_is_its_section_and_its_id(key, path):
    assert reference_urls.reference_path("stand-in", key) == path


@pytest.mark.parametrize("key", ["AM1", "6.7", "6.7:", "seven:AM1", ""])
def test_an_id_without_its_section_is_an_error_and_not_a_guess(key):
    with pytest.raises(ValueError):
        reference_urls.reference_path("stand-in", key)


# --- The pages -------------------------------------------------------------------------------------------------------------------

def test_the_index_lists_every_row_with_a_link_and_the_banner(client, data):
    html = _html(client, "/references/stand-ins")
    total = len(references_site.rows_of(data))
    assert BANNER_ID in html and html.count(BANNER_ID) == 1 and '<body class="app references-mode">' in html
    assert f"{total:,} of {total:,} stand-ins" in html
    assert html.count('href="/references/stand-ins/') == total
    assert 'href="/references/stand-ins/5-1-gr1"' in html


def test_the_filters_narrow_the_rows_and_name_how_many_are_left(client, data):
    rows = [row for _part, row in references_site.rows_of(data)]
    total = len(rows)
    for query, wanted in (("tier=T3", [r for r in rows if r["tier_code"] == "T3"]), ("tweak=L4", [r for r in rows if r["tweak"] == "L4"]),
                          ("tier=T9&tweak=L3", [r for r in rows if r["tier_code"] == "T9" and r["tweak"] == "L3"])):
        html = _html(client, f"/references/stand-ins?{query}")
        assert f"{len(wanted):,} of {total:,} stand-ins" in html, query
        assert html.count('href="/references/stand-ins/') == len(wanted), query
    html = _html(client, "/references/stand-ins?section=6-7-stand-ins-of-the-atmosphere-mixture-draw")
    assert 'href="/references/stand-ins/6-7-am1"' in html and 'href="/references/stand-ins/5-1-gr1"' not in html
    html = _html(client, "/references/stand-ins?text=zzzznothingzzzz")
    assert "0 of" in html and "No stand-in matches these filters." in html


def test_a_page_for_a_row_says_what_it_is_how_sure_it_is_and_what_replaces_it(client):
    html = _html(client, "/references/stand-ins/6-7-am64")
    assert BANNER_ID in html and "<title>Stand-in AM64" in html
    assert "Its source tier is" in html and "Its tweak level is" in html
    assert "The item GEN.804 replaces it with data." in html
    assert "How it was made and decided" in html and "<dt>What it is</dt>" in html


def test_a_citation_on_a_row_page_is_a_link_to_the_papers_page(client, data):
    row = next(row for _part, row in references_site.rows_of(data) if any("[R" in text for _label, text in row["fields"]))
    html = _html(client, f"/references/stand-ins/{row['slug']}")
    assert re.search(r'<a class="citation" href="/references/(?:papers|data)#paper-r\d{4,}"', html), row["key"]
    assert not re.search(r"\[R\d{4,}", html)


def test_every_stand_in_has_a_page_and_a_sample_of_them_render(client, data, overview):
    rows = references_site.rows_of(data)
    for _part, row in rows:
        page = references_site.standin(data, overview, row["slug"])
        assert page["row"]["key"] == row["key"] and page["paragraphs"], row["key"]
    for _part, row in rows[::60]:
        response = client.get(f"/references/stand-ins/{row['slug']}")
        assert response.status_code == 200, row["key"]
        assert BANNER_ID in response.get_data(as_text=True), row["key"]


def test_a_row_the_registry_does_not_hold_is_not_found(client):
    assert client.get("/references/stand-ins/9-9-nothing").status_code == 404
    assert client.get("/references/stand-ins/am1").status_code == 404


def test_a_page_links_to_the_last_and_next_row_and_back_to_its_section(client):
    html = _html(client, "/references/stand-ins/5-1-gr2")
    assert 'href="/references/stand-ins/5-1-gr1">Last: GR1</a>' in html and 'href="/references/stand-ins/5-1-gr3">Next: GR3</a>' in html
    assert 'href="/references/stand-ins#section-' in html


def test_the_fitted_page_lists_every_row_at_l2_to_l4_and_every_target(client, data):
    html = _html(client, "/references/fitted")
    rows = [row for _part, row in references_site.rows_of(data) if row["tweak"] in ("L2", "L3", "L4")]
    assert BANNER_ID in html and len(rows) > 100
    assert sum(html.count(f'href="/references/stand-ins/{row["slug"]}"') for row in rows) == len(rows)
    table = json.loads((ROOT / "src" / "planetgen" / "physics" / "data" / "tuning_targets.json").read_text(encoding="utf-8"))
    assert f"The tuning targets ({len(table['targets'])})" in html
    assert all(f"<th scope=\"row\">{target['id']}</th>" in html for target in table["targets"])
    assert re.search(r'<a class="citation" href="/references/papers#paper-r\d{4,}"', html)


def test_a_site_that_has_not_built_the_data_says_so_on_every_page(client, monkeypatch):
    monkeypatch.setattr(references_site, "load", lambda name, folder=None: None)
    for path in ("/references/stand-ins", "/references/stand-ins/5-1-gr1", "/references/fitted"):
        html = _html(client, path)
        assert BANNER_ID in html and "The pages are not built" in html, path


def test_the_sentences_of_the_pages_keep_the_guide(data, overview):
    targets = json.loads((ROOT / "docs" / "references" / "fitted.json").read_text(encoding="utf-8"))
    index = references_site.standins(data, overview)
    fitted = references_site.fitted(data, targets, overview)
    texts = list(index["intro"]) + list(fitted["intro"])
    for _part, row in list(references_site.rows_of(data))[::40]:
        texts += references_site.standin(data, overview, row["slug"])["paragraphs"][:3]
    broken = [(sentence, rule) for text in texts for sentence in plainwords.sentences(text)
              for rule, _detail in plainwords.check_sentence(sentence, references=True)]
    assert not broken, broken[:5]


# --- The Papers page links the stand-ins that cite a source ------------------------------------------------------------------------

def test_a_paper_entry_links_the_stand_in_rows_that_cite_it_when_their_pages_exist():
    row = {"id": "R0001", "label": "A 2000", "type": "Journal article", "field": "Other", "title": "t", "who": "A."}
    usage = {"uses": {"R0001": {"values": [], "notes": [], "stand_ins": ["5.1:GR1", "9.9:ZZ9", "AB1"], "areas": ["galaxy"]}}}
    names = references_site.paper_entry(row, usage)["stand_ins"]["names"]
    assert names == [{"text": "GR1", "url": "/references/stand-ins/5-1-gr1"}, {"text": "ZZ9", "url": ""}, {"text": "AB1", "url": ""}]


def test_the_scan_keeps_the_section_with_the_row_id(tmp_path):
    from planetgen.web.lib import references_usage
    plan = tmp_path / "docs" / "design" / "research"
    plan.mkdir(parents=True)
    (plan / "gap-filling-plan.md").write_text("## 5.1 Galaxy laws\n\n| Id | Fill |\n|---|---|\n| GR1 | [R0001] |\n", encoding="utf-8")
    assert references_usage.scan(str(tmp_path))["uses"]["R0001"]["stand_ins"] == ["5.1:GR1"]
