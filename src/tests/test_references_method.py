# tests/test_references_method.py

"""
DOC.91 (References R5): the method pages (`/references/method/<area>`), what we have and lack (`/references/gaps`) and the computed
scorecard (`/references/scorecard`). The properties page left `/help/properties` and its rows are on the method pages. Every number on
these pages is counted from the stand-in registry (`docs/references/stand-ins.json`), the sources, and the property sheet.
"""

import re

import pytest

from planetgen.web.lib import property_help, reference_urls, references_site, references_usage
from planetgen.web.lib.reference_urls import AREAS

from tests.test_class_pages import _html, app, client  # noqa: F401
from tests.test_references_crawl import BANNER


@pytest.fixture(scope="module")
def data():
    return references_site.load("stand-ins.json")


def test_every_row_of_the_property_sheet_is_on_exactly_one_method_page():
    by_area = property_help.by_area()
    assert set(by_area) == set(AREAS)
    listed = [(row["group"], row["label"], row["source"]) for rows in by_area.values() for row in rows]
    every = [(row["group"], row["label"], row["source"]) for row in property_help.rows()]
    assert sorted(listed) == sorted(every) and len(listed) == len(set(listed))
    assert all(row["applies"] in ("Stars", "Neutron stars") for row in by_area["stars"]) and by_area["stars"]


def test_every_group_of_the_property_sheet_names_an_area():
    from planetgen.describe import sheet
    assert set(sheet.GROUPS.values()) <= set(property_help.GROUP_AREA)
    assert set(property_help.GROUP_AREA.values()) <= set(AREAS)


@pytest.mark.parametrize("area", AREAS)
def test_a_method_page_opens_with_the_banner_and_lists_the_stand_ins_of_its_area(client, data, area):  # noqa: F811
    html = _html(client, f"/references/method/{area}")
    assert BANNER in html and f"<title>Method: {references_usage.AREA_TITLES[area]}" in html
    keys = [row["key"] for section, row in references_site.rows_of(data) if references_usage.area_of(section["title"]) == area]
    for key in keys[:20]:
        assert f'href="{reference_urls.reference_path("stand-in", key)}"' in html
    for row in property_help.by_area()[area]:
        assert f'<tr id="{reference_urls.anchor(row["label"])}"' in html, row["label"]


def test_a_value_link_lands_on_the_row_of_its_name(client):  # noqa: F811
    html = _html(client, "/references/method/orbits")
    assert 'id="eccentricity"' in html
    path, fragment = reference_urls.reference_path("method", "orbits/eccentricity").split("#")
    assert path == "/references/method/orbits" and fragment == "eccentricity"


def test_an_unknown_area_is_not_found(client):  # noqa: F811
    assert client.get("/references/method/nothing").status_code == 404


def test_the_old_properties_page_is_gone_with_no_redirect(client):  # noqa: F811
    response = client.get("/help/properties")
    assert response.status_code == 404 and "Location" not in response.headers


def test_the_body_pages_link_to_the_method_pages_not_the_old_page():
    from pathlib import Path
    html = (Path(__file__).resolve().parents[1] / "planetgen" / "web" / "templates" / "item.html").read_text(encoding="utf-8")
    assert "property_help_page" not in html and 'reference_link("method"' in html


def test_the_scorecard_counts_come_from_the_registry(client, data):  # noqa: F811
    page = references_site.scorecard_page(data, references_usage.uses())
    rows = page["sections"][0]["table"]["rows"]
    assert len(rows) == len(AREAS) + 2 and rows[-1][0] == "All areas"
    total = sum(len(section["rows"]) for section in data["sections"])
    assert rows[-1][1] == f"{total:,}" and sum(int(row[1].replace(",", "")) for row in rows[:-1]) == total
    html = _html(client, "/references/scorecard")
    assert BANNER in html and 'id="life"' in html and 'id="all"' in html


def test_the_gaps_page_lists_the_weak_values_of_every_area(client, data):  # noqa: F811
    page = references_site.gaps_page(data, references_usage.uses())
    ids = [section["id"] for section in page["sections"]]
    assert ids[0] == "summary" and set(AREAS) <= set(ids)
    summary = page["sections"][0]["table"]["rows"]
    weak = sum(int(row[-1].replace(",", "")) for row in summary)
    expected = sum(1 for _section, row in references_site.rows_of(data) if row["tier_code"] in references_site.LACKING_TIERS or row["tweak"] == "L4")
    assert weak == expected
    html = _html(client, "/references/gaps")
    assert BANNER in html and re.search(r'<section class="panel" id="life"', html)


def test_the_start_page_leads_to_every_new_page(client):  # noqa: F811
    html = _html(client, "/references")
    for area in AREAS:
        assert f'href="/references/method/{area}"' in html
    assert 'href="/references/gaps"' in html and 'href="/references/scorecard"' in html


def test_a_body_page_links_a_value_to_its_method_row_and_keeps_no_source_line():
    """DOC.93: where a body's page printed the one-line source it shows the quiet link; the line itself is on the method page."""
    assert property_help.method_ref("planet", "Eccentricity") == "orbits/Eccentricity"
    assert property_help.method_ref("star", "Rotation period").startswith("stars/")
    assert property_help.method_ref("planet", "No such value") is None
    from planetgen.web.helpers import value_source
    assert str(value_source("planet", "Eccentricity")) == (
        '<a class="reference-link" href="/references/method/orbits#eccentricity">More information</a>')
    assert str(value_source("planet", "No such value")) == ""
    for kind, label in (("planet", "Eccentricity"), ("moon", "Eccentricity"), ("star", "Rotation period")):
        path, _fragment = reference_urls.reference_path("method", property_help.method_ref(kind, label)).split("#")
        assert any(row["label"] == label for row in property_help.by_area()[path.rsplit("/", 1)[1]])
