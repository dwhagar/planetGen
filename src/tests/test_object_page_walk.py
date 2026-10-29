# tests/test_object_page_walk.py

"""
UX.177: every kind of object page, walked on the small generated database `test_web_a11y.py` builds: each star, planet,
moon, belt, comet and facility the database holds opens its page (property sheet, places, tree), each tree row's link opens the
page it names, a planet with no moons has no tree, and an object that does not exist is a 404.

Skipped without the MySQL test server. The browser behaviour of the tree is in `test_web_browser_object_tree.py`, the
accessibility of the pages in `test_web_a11y.py`.
"""

import re

import pytest

from planetgen.db import store
from planetgen.web import item_pages
from tests.publicids import pid
from tests.test_number_guard import bare_cells
from tests.test_web_a11y import site_app, site_db  # noqa: F401 -- fixtures

PUBLIC_ID = r"[0-9A-F]{10}-\d{7}-\d{3}"


WALKED = 6
"""int: How many objects of one kind the walk opens at most (the tests below slice at 6 or fewer)."""

WALKED_TABLES = {"star": ("stars", "1"), "planet": ("planets", "1"), "moon": ("moons", "planet_id IS NOT NULL"),
                 "belt": ("asteroid_belts", "1"), "comet": ("comets", "1"), "facility": ("facilities", "1")}
"""dict: The kinds the walk opens and where their rows are (a moon of a free planet has no system, so it is left out)."""


def _bodies(site_db):
    """`{kind: [ref, ...]}` for the first `WALKED` objects of each kind the database holds, read from the tables. It
    makes no API call: a loop over the systems' API answers spent the API's per-IP limit and then read a 429's body
    as a system with no "stars" (TEST.349, TEST.359)."""
    found = {kind: [] for kind in item_pages.ITEM_PAGE_KINDS}
    conn = store.get_connection(site_db)
    try:
        for kind, (table, where) in WALKED_TABLES.items():
            rows = conn.execute(f"SELECT id FROM {table} WHERE {where} ORDER BY id LIMIT {WALKED}").fetchall()
            found[kind] = [f"{kind}:{pid(kind, row['id'], conn)}" for row in rows]
    finally:
        conn.close()
    assert found["star"] and found["planet"], "the sample galaxy holds no star or no planet to walk"
    return found


@pytest.fixture(scope="module")
def client(site_app):
    return site_app.test_client()


@pytest.fixture(scope="module")
def bodies(site_db):
    return _bodies(site_db)


def test_every_kind_of_body_in_the_database_has_a_page_with_its_sheet_and_places(client, bodies):
    walked = set()
    for kind, refs in bodies.items():
        for ref in refs[:6]:
            response = client.get(f"/object/{ref}")
            assert response.status_code == 200, ref
            page = response.get_data(as_text=True)
            assert 'id="item-properties"' in page and 'id="location"' in page, ref
            assert 'href="/system/' in page, ref
            walked.add(kind)
    assert {"star", "planet"} <= walked, f"the database held no {set(item_pages.ITEM_PAGE_KINDS) - walked}"


def test_no_page_of_a_body_shows_a_number_with_no_unit_scale_or_noun(client, bodies):
    """GEN.772: the property sheet and the places table of every walked body hold no cell that is a number alone."""
    for refs in bodies.values():
        for ref in refs[:6]:
            page = client.get(f"/object/{ref}").get_data(as_text=True)
            assert not bare_cells(page), f"{ref}: {bare_cells(page)}"


def test_every_row_of_a_stars_tree_opens_the_page_it_names(client, bodies):
    checked = 0
    for ref in bodies["star"][:5]:
        page = client.get(f"/object/{ref}").get_data(as_text=True)
        for href in set(re.findall(r'class="tree-name" href="([^"]+)"', page)):
            assert client.get(href).status_code == 200, href
            checked += 1
            for below in client.get(href + "/tree").get_json().get("html", "").split('class="tree-name" href="')[1:]:
                assert client.get(below.split('"', 1)[0]).status_code == 200
    assert checked, "no star in the database holds a body"


def test_a_level_of_the_tree_is_served_for_every_parent_kind(client, bodies):
    for kind in ("star", "planet"):
        for ref in bodies[kind][:4]:
            response = client.get(f"/object/{ref}/tree")
            assert response.status_code == 200, ref
            assert set(response.get_json()) == {"html", "next"}
            assert response.headers["Cache-Control"] == "no-store"


def test_a_body_with_nothing_below_it_has_no_tree_panel(client, bodies):
    for ref in bodies["moon"][:3] + bodies["comet"][:3]:
        page = client.get(f"/object/{ref}").get_data(as_text=True)
        assert 'id="item-members"' not in page, ref
    assert client.get(f"/object/{bodies['comet'][0]}/tree").get_json() == {"html": "", "next": None} if bodies["comet"] else True


@pytest.mark.parametrize("ref", ["planet:FE00000000-0000000-001", "star:FE00000000-0000000-000", "moon:nonsense", "planet:",
                                 "gizmo:FE81000A2B-0000005-001"])
def test_an_object_that_does_not_exist_is_a_404(client, ref):
    assert client.get(f"/object/{ref}").status_code == 404
    assert client.get(f"/object/{ref}/tree").status_code in (404, 400)
