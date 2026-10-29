# tests/test_object_page_map_symbols.py

"""
UX.233: the point-of-interest symbols beside the bodies of an object page's tree, and the collapsed "Map symbols" key that follows them
(`web/lib/poi_icons.mark_of_body`, `web/lib/objecttree.py`, `partials/poi_key.html`). The key is shown only on a page that draws a
symbol, its lines are the Galaxy Map's (one macro draws both), and the symbols are the ones the map's system level gives a body
(`static/bodymarkers.js`).

The page tests are skipped without the MySQL test server.
"""

import pathlib
import re

import pytest
from flask import render_template_string

from planetgen.db import store
from planetgen.physics import habitability
from planetgen.physics.habitability_world import DOMAINS
from planetgen.web import item_pages, system_pages
from planetgen.web.lib import objecttree, poi_icons
from tests.publicids import pid
from tests.test_web_a11y import site_app, site_db  # noqa: F401 -- fixtures

TEMPLATES = pathlib.Path(item_pages.__file__).parent / "templates"

CONTEXT = objecttree.TreeContext(item_url=lambda ref: "/object/" + ref.replace(":", "-"),
                                 tree_url=lambda ref: "/object/" + ref.replace(":", "-") + "/tree",
                                 static_url=lambda path: "/static/" + path)


SCORES = {"Blue": 1.0, "Green": 0.8, "Yellow": 0.5, "Red": 0.2}
"""dict: A stored PHI-6 score in each colour (`habitability.tier`: Blue at 1, Green from 0.75, Yellow from 0.4)."""


def _score(tier_name):
    assert habitability.tier(SCORES[tier_name]) == tier_name
    return SCORES[tier_name]


@pytest.mark.parametrize("kind,summary,expected", [
    ("belt", {}, ("asteroid-field", None)),
    ("planet", {"body_type": "g"}, ("gas-giant", None)),
    ("planet", {"body_type": "t", "phi4": _score("Blue")}, ("habitable-world", 0)),
    ("moon", {"body_type": "t", "phi4": _score("Green")}, ("habitable-world", 1)),
    ("planet", {"body_type": "t", "phi4": _score("Yellow")}, ("habitable-world", 2)),
    ("planet", {"body_type": "t", "phi4": _score("Red")}, ("uninhabitable-world", 3)),
    ("planet", {"body_type": "t", "phi4": None}, None),
    ("star", {"phi4": 0.9}, None),
    ("comet", {}, None),
])
def test_a_body_gets_the_symbol_the_maps_system_level_gives_it(kind, summary, expected):
    assert poi_icons.mark_of_body(kind, summary) == expected


def test_a_row_draws_its_symbol_first_and_a_row_without_one_draws_none():
    habitable = {"ref": "planet:A", "kind": "planet", "id": "A", "name": "Garden", "child_count": 0, "counts": {},
                 "summary": {"body_type": "t", "phi4": _score("Green")}}
    star = {"ref": "star:B", "kind": "star", "id": "B", "name": "Sol", "child_count": 0, "counts": {}, "summary": {}}
    html = objecttree.nodes_html([habitable, star], CONTEXT)
    assert html.count('data-poi="habitable-world"') == 1 and "poi-tier-1" in html
    assert html.index("tree-poi") < html.index("Class") if "Class" in html else True
    assert 'data-poi' not in objecttree.nodes_html([star], CONTEXT)
    assert objecttree.has_poi_mark([habitable, star]) and not objecttree.has_poi_mark([star])


def test_a_tree_with_no_static_url_draws_no_symbol():
    bare = objecttree.TreeContext(item_url=CONTEXT.item_url, tree_url=CONTEXT.tree_url)
    gas = {"ref": "planet:G", "kind": "planet", "id": "G", "name": "Big", "child_count": 0, "counts": {}, "summary": {"body_type": "g"}}
    assert "data-poi" not in objecttree.nodes_html([gas], bare)


def test_the_galaxy_page_and_the_object_pages_draw_the_key_with_the_same_macro():
    galaxy = (TEMPLATES / "galaxy.html").read_text(encoding="utf-8")
    key = (TEMPLATES / "partials" / "poi_key.html").read_text(encoding="utf-8")
    assert "{{ poi_key() }}" in galaxy
    assert re.search(r"macro poi_key_panel\(.*?\{\{ poi_key\(\) \}\}", key, re.S)
    for name in ("system.html", "item.html"):
        assert "poi_key_panel(" in (TEMPLATES / name).read_text(encoding="utf-8"), name


@pytest.fixture(scope="module")
def client(site_app):
    return site_app.test_client()


@pytest.fixture(scope="module")
def system_ref(site_db):
    conn = store.get_connection(site_db)
    try:
        row = conn.execute("SELECT s.id FROM star_systems s JOIN planets p ON p.star_system_id = s.id ORDER BY s.id LIMIT 1").fetchone()
        assert row, "the sample galaxy holds no system with a planet"
        return pid("system", row["id"], conn)
    finally:
        conn.close()


def _key(page):
    match = re.search(r'<dl class="poi-key".*?</dl>', page, re.S)
    return match.group(0) if match else None


def _patched(monkeypatch, change):
    original = item_pages.apiclient.get_object_children

    def children(db, ref, limit=None, after=None):
        page = original(db, ref, limit=limit, after=after)
        for item in page["items"]:
            change(item)
        return page

    monkeypatch.setattr(item_pages.apiclient, "get_object_children", children)


def test_a_system_page_with_a_habitable_world_shows_the_collapsed_key_with_the_maps_lines(client, system_ref, monkeypatch, site_app):
    def make_habitable(item):
        if item["kind"] == "planet":
            item["summary"].update(body_type="t", phi4=_score("Green"), **{f"tier_{domain}": 1 for domain in DOMAINS})

    _patched(monkeypatch, make_habitable)
    page = client.get(f"/system/{system_ref}").get_data(as_text=True)
    assert 'data-poi="habitable-world"' in page
    panel = re.search(r'<details class="panel" id="system-map-symbols">\s*<summary>Map symbols</summary>', page)
    assert panel, "the key is a collapsed panel under the tree"
    with site_app.test_request_context():
        galaxy_key = render_template_string('{% from "partials/poi_key.html" import poi_key with context %}{{ poi_key() }}')
    assert _key(page) == _key(galaxy_key) and _key(page).count("<dt ") == len(poi_icons.ICONS)


def test_a_system_page_with_no_symbol_in_its_tree_has_no_key(client, system_ref, monkeypatch):
    def make_barren_unscored(item):
        if item["kind"] == "planet":
            item["summary"].update(body_type="t", phi4=None)
        if item["kind"] == "belt":
            item["kind"] = "comet"

    _patched(monkeypatch, make_barren_unscored)
    page = client.get(f"/system/{system_ref}").get_data(as_text=True)
    assert "data-poi=" not in page and "system-map-symbols" not in page and 'class="poi-key"' not in page
