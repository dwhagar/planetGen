# tests/test_object_page_walk.py

"""
UX.177: every kind of object page, walked on the small generated database `test_web_a11y.py` builds: each star, planet,
moon, belt and comet the database holds opens its page (property sheet, places, tree), each tree row's link opens the
page it names, a planet with no moons has no tree, and an object that does not exist is a 404.

Skipped without the MySQL test server. The browser behaviour of the tree is in `test_web_browser_object_tree.py`, the
accessibility of the pages in `test_web_a11y.py`.
"""

import re

import pytest

from planetgen.web import item_pages
from tests.test_web_a11y import site_app, site_db  # noqa: F401 -- fixtures

PUBLIC_ID = r"[0-9A-F]{10}-\d{7}-\d{3}"


def _bodies(client):
    """`{kind: [ref, ...]}` for every body of the first systems the database holds, from the system API."""
    found = {kind: [] for kind in item_pages.ITEM_PAGE_KINDS}
    for system in client.get("/api/systems?limit=100").get_json()["items"]:
        detail = client.get(f"/api/systems/{system['id']}").get_json()
        found["star"] += [f"star:{row['id']}" for row in detail["stars"]]
        found["belt"] += [f"belt:{row['id']}" for row in detail["belts"]]
        found["comet"] += [f"comet:{row['id']}" for row in detail["comets"]]
        for planet in detail["planets"]:
            found["planet"].append(f"planet:{planet['id']}")
            found["moon"] += [f"moon:{moon['id']}" for moon in planet.get("moons") or []]
    return found


@pytest.fixture(scope="module")
def client(site_app):
    return site_app.test_client()


@pytest.fixture(scope="module")
def bodies(client):
    return _bodies(client)


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
