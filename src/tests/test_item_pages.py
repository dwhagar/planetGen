# tests/test_item_pages.py

"""
UX.143: every kind of item a reference can name has a page, from one registry, and a star, planet, moon, belt or
comet gets the item page at `/object/<ref>`.
"""

import pytest

from planetgen.describe import model
from planetgen.galaxy import objectref
from planetgen.web import item_pages
from planetgen.web.app import create_app
from planetgen.web.lib import apiclient
from tests.test_web_pages import _FakeConfig


PLANET = "planet:FE81000A2B-0000005-001"
MOON = "moon:FE81000A2B-0000005-002"


def _document():
    moon = model.section("moon", MOON, "Tethys", [model.facts("Moon", [model.fact("mass", "Mass", model.qty("mass", 7.3e22))])])
    planet = model.section("planet", PLANET, "Kepler b", [
        model.facts("Planet", [model.fact("mass", "Mass", model.qty("mass", 5.97e24))]),
        model.prose("description", "A rocky world."),
    ], [moon])
    system = model.section("system", "system:9", "Kepler", [], [planet])
    return model.document("system", "system:9", "Kepler", [system])


@pytest.fixture
def client(monkeypatch):
    objects = {
        PLANET: {"ref": PLANET, "name": "Kepler b", "parents": [
            {"ref": "galaxy", "kind": "galaxy"},
            {"ref": "sector:4", "kind": "sector", "name": "Sector 004"},
            {"ref": "system:9", "kind": "system", "name": "Kepler"}]},
        MOON: {"ref": MOON, "name": "Tethys", "parents": [
            {"ref": "galaxy", "kind": "galaxy"},
            {"ref": "sector:4", "kind": "sector", "name": "Sector 004"},
            {"ref": "system:9", "kind": "system", "name": "Kepler"},
            {"ref": PLANET, "kind": "planet", "name": "Kepler b"}]},
        "star:FE81000A2B-0000005-000": {"ref": "star:FE81000A2B-0000005-000", "name": "Kepler", "parents": []},
    }

    def get_object(db, ref):
        if ref not in objects:
            raise apiclient.NotFoundError("no such object")
        return objects[ref]

    monkeypatch.setattr(apiclient, "get_object", get_object)
    monkeypatch.setattr(apiclient, "get_system_description", lambda db, system_id: _document())
    monkeypatch.setattr(apiclient, "get_system", lambda db, system_id: {"inside": None})
    app = create_app(_FakeConfig)
    app.testing = True
    return app.test_client()


def test_every_kind_a_reference_can_name_has_a_page():
    assert sorted(item_pages.page_kinds()) == sorted(objectref.KINDS)
    assert not set(item_pages.OWN_PAGE_KINDS) & set(item_pages.ITEM_PAGE_KINDS)
    assert set(item_pages.KIND_LABELS) == set(item_pages.ITEM_PAGE_KINDS)


def test_a_planet_has_a_page_with_its_facts_places_and_members(client):
    response = client.get(f"/object/{PLANET}")
    assert response.status_code == 200
    page = response.get_data(as_text=True)
    assert "<title>" in page and "Kepler b" in page
    assert "A rocky world." in page
    assert "Mass" in page
    assert 'href="/sector/4"' in page          # the places it is in
    assert 'href="/system/9"' in page
    assert f"/object/{MOON}" in page            # the bodies that belong to it, on pages of their own
    assert "#planet-FE81000A2B-0000005-001" in page   # back to its row in the system's list


def test_a_moon_names_the_planet_it_circles(client):
    page = client.get(f"/object/{MOON}").get_data(as_text=True)
    assert f"/object/{PLANET}" in page


def test_a_body_with_no_system_is_not_found(client):
    assert client.get("/object/star:FE81000A2B-0000005-000").status_code == 404


def test_an_unknown_body_is_not_found(client):
    assert client.get("/object/planet:FE81000A2B-0000005-099").status_code == 404


def test_the_page_has_a_location_block_from_the_parent_chain(client):
    page = client.get(f"/object/{MOON}").get_data(as_text=True)
    assert 'id="location"' in page
    for url in ('href="/galaxy"', 'href="/sector/4"', 'href="/system/9"', f'href="/object/{PLANET}"'):
        assert url in page
    assert "Not inside any nebula or supernova remnant." in page
