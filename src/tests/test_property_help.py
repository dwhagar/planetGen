# tests/test_property_help.py

"""
UX.187: the properties help page (`/help/properties`) says, for each physics value of a star, planet or moon page, its unit,
which bodies have it and the one-line source its page shows under it. It is built from the property sheet and its source
table, so the page and the pages it explains cannot disagree; it is in the Help menu and linked from a body's page.
"""

import re

from markupsafe import escape

from planetgen.describe import sheet
from planetgen.describe.sheet_sources import SOURCES
from planetgen.web.lib import property_help
from tests.test_class_pages import _html, app, client  # noqa: F401


def test_every_source_line_is_on_the_page_once_with_its_label_unit_and_who_has_it():
    rows = property_help.rows()
    lines = {line for kind in SOURCES.values() for line in kind.values()}
    assert {row["source"] for row in rows} == lines
    for row in rows:
        assert row["label"] and row["unit"] and row["applies"]
        assert "—" not in row["unit"] + row["applies"]
    seen = [(row["label"], row["source"]) for row in rows]
    assert len(seen) == len(set(seen))


def test_a_value_that_planets_and_moons_share_is_one_row_and_a_planet_only_one_says_so():
    by_label = {(row["group"], row["label"], row["applies"]): row for row in property_help.rows()}
    assert ("Orbit", "Eccentricity", "Planets, Moons") in by_label
    assert ("Orbit", "Inclination", "Planets") in by_label
    assert ("Composition", "Ring outer edge", "Planets") in by_label
    assert ("Rotation", "Rotation period", "Stars") in by_label
    assert ("Composition", "Alpha enhancement [alpha/Fe]", "Stars") in by_label    # UX.254: a row worked out from the stored values


def test_the_unit_column_reads_as_a_unit_a_noun_or_words():
    units = {row["label"]: row["unit"] for row in property_help.rows()}
    assert units["Tectonic regime"] == property_help.WORDS and units["Life stage"] == property_help.WORDS
    assert units["Rotation period"] == "h" and units["Magnetic moment"] == "A·m²"
    assert units["Eccentricity"].startswith("a number")


def test_the_page_lists_every_row_under_its_sheet_heading(client):  # noqa: F811
    html = _html(client, "/help/properties")
    assert "<title>Properties" in html
    headings = re.findall(r'<h2 id="properties-\d+-heading">([^<]*)</h2>', html)
    assert headings == [group["label"] for group in property_help.page()["groups"]]
    for row in property_help.rows():
        assert f"<td>{escape(row['source'])}</td>" in html, row["label"]
        assert f'<th scope="row">{escape(row["label"])}</th><td>{escape(row["unit"])}</td>' in html


def test_the_help_menu_opens_the_page_and_marks_it_current(client):  # noqa: F811
    html = _html(client, "/help/properties")
    assert 'data-href="/help/properties" aria-current="page">Properties</sl-menu-item>' in html
    assert 'class="site-help-trigger" aria-haspopup="menu" aria-current="page"' in html


def test_the_sheet_rows_the_page_explains_are_rows_of_the_sheet():
    for kind, name in property_help.KINDS:
        keys = {spec.key for spec in sheet.SPECS[kind]} | (set(sheet.CHEMISTRY_KEYS) if kind == "star" else set())    # UX.254
        assert set(SOURCES[kind]) <= keys, name
