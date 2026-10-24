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
FREE_PLANET = "rogue_planet:FE81000A2B-0000009"
FREE_MOON = "moon:FE81000A2B-0000009-001"
FACILITY = "facility:FE81000A2B-0000005-009"
LONE_FACILITY = "facility:FE81000A2B-000000A"


def _document():
    moon = model.section("moon", MOON, "Tethys", [model.facts("Moon", [model.fact("mass", "Mass", model.qty("mass", 7.3e22))])])
    planet = model.section("planet", PLANET, "Kepler b", [
        model.facts("Planet", [model.fact("mass", "Mass", model.qty("mass", 5.97e24, "kg"))]),
        model.prose("description", "A rocky world."),
    ], [moon])
    system = model.section("system", "system:9", "Kepler", [], [planet])
    return model.document("system", "system:9", "Kepler", [system])


STAR = "star:FE81000A2B-0000005-000"


def _detail():
    star = {"id": "FE81000A2B-0000005-000", "name": "Kepler", "role": "single", "star_type": "G2V"}
    moon = {"id": "FE81000A2B-0000005-002", "name": "Tethys", "star_id": None}
    other = {"id": "FE81000A2B-0000005-003", "name": "Kepler c", "star_id": star["id"], "moons": []}
    planet = {"id": "FE81000A2B-0000005-001", "name": "Kepler b", "star_id": star["id"], "moons": [moon]}
    return {"inside": None, "stars": [star], "planets": [planet, other], "belts": [], "comets": []}


def _children(ref):
    moon = {"ref": MOON, "kind": "moon", "id": "FE81000A2B-0000005-002", "name": "Tethys", "child_count": 0, "counts": {},
            "summary": {"planet_class": "C", "distance_km": 294660.0, "period_years": 0.0052}}
    planet = {"ref": PLANET, "kind": "planet", "id": "FE81000A2B-0000005-001", "name": "Kepler b", "child_count": 1,
              "counts": {"moon": 1}, "summary": {"planet_class": "J", "distance_km": 1.5e8, "period_years": 1.0}}
    items = {PLANET: [moon], STAR: [planet]}.get(ref, [])
    return {"ref": ref, "kind": ref.partition(":")[0], "total": len(items), "items": items, "next": None}


def _sheet(ref):
    mass = model.qty("mass", 5.97e24, "kg")
    return {"kind": ref.partition(":")[0], "groups": [
        {"key": "physical", "label": "Physical", "rows": [
            {"key": "mass_kg", "label": "Mass", "value": [mass], "text": "1 Earth mass", "stored": True},
            {"key": "ocean_depth_km", "label": "Ocean depth", "value": ["not stored"], "text": "not stored",
             "stored": False}]}]}


@pytest.fixture
def client(monkeypatch):
    objects = {
        PLANET: {"ref": PLANET, "name": "Kepler b", "phi4": 0.8, "tier_pressure": 0, "tier_temperature": 1, "tier_air": 2,
                 "tier_water": 0, "tier_radiation": 3, "tier_suitability": 1, "parents": [
            {"ref": "galaxy", "kind": "galaxy"},
            {"ref": "sector:4", "kind": "sector", "name": "Sector 004"},
            {"ref": "system:9", "kind": "system", "name": "Kepler"}]},
        MOON: {"ref": MOON, "name": "Tethys", "parents": [
            {"ref": "galaxy", "kind": "galaxy"},
            {"ref": "sector:4", "kind": "sector", "name": "Sector 004"},
            {"ref": "system:9", "kind": "system", "name": "Kepler"},
            {"ref": PLANET, "kind": "planet", "name": "Kepler b"}]},
        STAR: {"ref": STAR, "name": "Kepler", "parents": []},
        FACILITY: {"ref": FACILITY, "name": "New Hope", "parents": [
            {"ref": "galaxy", "kind": "galaxy"},
            {"ref": "sector:4", "kind": "sector", "name": "Sector 004"},
            {"ref": "system:9", "kind": "system", "name": "Kepler"},
            {"ref": PLANET, "kind": "planet", "name": "Kepler b"}]},
        LONE_FACILITY: {"ref": LONE_FACILITY, "name": "Lone Light", "parents": [
            {"ref": "galaxy", "kind": "galaxy"},
            {"ref": "sector:4", "kind": "sector", "name": "Sector 004"}]},
        FREE_MOON: {"ref": FREE_MOON, "name": "Lone Moon", "phi4": 0.5, "tier_pressure": 1, "tier_temperature": 1,
                    "tier_air": 1, "tier_water": 1, "tier_radiation": 1, "tier_suitability": 1, "parents": [
            {"ref": "galaxy", "kind": "galaxy"},
            {"ref": FREE_PLANET, "kind": "rogue_planet", "name": "Drifter"}]},
    }

    def get_object(db, ref):
        if ref not in objects:
            raise apiclient.NotFoundError("no such object")
        return objects[ref]

    def get_facility(db, facility_id):
        return {"id": facility_id, "name": "New Hope", "kind": "colony", "type_id": None, "affiliation_name": "Velarix Concord",
                "description": "A first colony.", "ratings": {"crime": "yellow", "housing": None, "resources": None,
                                                              "maintenance": None, "health": None},
                "rating_words": {"crime": "Elevated"}}

    monkeypatch.setattr(apiclient, "get_object", get_object)
    monkeypatch.setattr(apiclient, "get_facility", get_facility)
    monkeypatch.setattr(apiclient, "get_system_description", lambda db, system_id: _document())
    monkeypatch.setattr(apiclient, "get_system", lambda db, system_id: _detail())
    monkeypatch.setattr(apiclient, "get_system_facilities", lambda db, system_id: [])
    monkeypatch.setattr(apiclient, "get_object_sheet", lambda db, ref: _sheet(ref))
    monkeypatch.setattr(apiclient, "get_object_children", lambda db, ref, limit=None, after=None: _children(ref))
    app = create_app(_FakeConfig)
    app.testing = True
    return app.test_client()


def test_every_kind_a_reference_can_name_has_a_page():
    assert sorted(item_pages.page_kinds()) == sorted(objectref.KINDS)
    assert not set(item_pages.OWN_PAGE_KINDS) & set(item_pages.ITEM_PAGE_KINDS)
    assert set(item_pages.KIND_LABELS) == set(item_pages.ITEM_PAGE_KINDS)


def test_a_planet_has_a_page_with_its_sheet_text_places_and_tree(client):
    response = client.get(f"/object/{PLANET}")
    assert response.status_code == 200
    page = response.get_data(as_text=True)
    assert "<title>" in page and "Kepler b" in page
    assert "A rocky world." in page                      # the description's paragraph
    assert 'id="item-properties"' in page and "Ocean depth" in page and "sheet-missing" in page   # the sheet, a gap shown
    assert 'href="/sector/4"' in page and 'href="/system/9"' in page          # the places it is in
    assert f'href="/object/{STAR}"' in page               # the star it orbits is on the chain
    assert 'data-tree-root' in page and f"/object/{MOON}" in page            # the bodies that belong to it, each a page
    assert page.index('id="location"') < page.index('id="item-properties"') < page.index('id="item-facts"') \
        < page.index('id="item-members"')
    assert "#planet-" not in page, "no more anchors on the system's page"


def test_the_sheet_shows_a_quantity_in_the_readers_units_and_the_tree_row_says_what_it_is(client):
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert 'data-q="mass"' in page                       # a quantity the unit script can rewrite
    assert 'class="tree-name"' in page and "Class" in page


def test_a_scored_planet_shows_its_symbols_under_the_title_and_other_kinds_do_not(client):
    """UX.100: the item page of a scored planet or moon carries the seven PHI-6 symbols; a star or unscored body none."""
    body = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert body.count("hab-symbol ") == 7 and 'aria-label="Radiation: Red"' in body
    assert body.index("hab-symbols") < body.index("About Kepler b")
    # UX.188: the symbols are the first row of the At a glance table, above the property sheet
    assert body.index('id="item-glance"') < body.index("hab-symbols") < body.index('id="item-properties"')
    assert 'class="glance-habitability"' in body
    assert "hab-symbol" not in client.get(f"/object/{MOON}").get_data(as_text=True)  # the fake moon was not scored
    assert "hab-symbol" not in client.get(f"/object/{STAR}").get_data(as_text=True)


def test_a_moon_names_the_planet_it_circles_and_its_star(client):
    page = client.get(f"/object/{MOON}").get_data(as_text=True)
    assert f"/object/{PLANET}" in page and f"/object/{STAR}" in page
    assert page.index(f"/object/{STAR}") < page.index(f"/object/{PLANET}"), "the star comes before the planet"
    assert 'data-tree-root' not in page, "a moon holds nothing, so there is no tree"


def test_a_planet_links_to_the_next_planet_out(client):
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert 'rel="next"' in page and "Kepler c" in page and 'rel="prev"' not in page


def test_a_body_with_no_system_is_not_found(client):
    assert client.get(f"/object/{STAR}").status_code == 404


def test_a_facility_has_a_page_with_its_host_its_record_and_no_tree(client):
    """UX.183: where it stands, what its record says and its sheet; it holds nothing, so there is no tree."""
    page = client.get(f"/object/{FACILITY}").get_data(as_text=True)
    assert "New Hope" in page and "Object ID FE81000A2B-0000005-009" in page
    assert f'href="/object/{PLANET}"' in page and "Built on" in page          # the planet it stands on
    assert "A first colony." in page and "affiliated with Velarix Concord" in page and "Crime" in page
    assert 'id="item-properties"' in page and "Ocean depth" in page            # the sheet, from the fake
    assert "item-members" not in page and "Navigate" not in page
    assert 'href="/system/9"' in page


def test_a_standalone_facility_has_a_page_with_no_system(client):
    page = client.get(f"/object/{LONE_FACILITY}").get_data(as_text=True)
    assert "Lone Light" in page and 'href="/sector/4"' in page
    assert "Part of the" not in page and "free planet" not in page


def test_an_unknown_body_is_not_found(client):
    assert client.get("/object/planet:FE81000A2B-0000005-099").status_code == 404


def test_the_page_has_a_location_block_from_the_parent_chain(client):
    page = client.get(f"/object/{MOON}").get_data(as_text=True)
    assert 'id="location"' in page
    for url in ('href="/galaxy"', 'href="/sector/4"', 'href="/system/9"', f'href="/object/{PLANET}"'):
        assert url in page
    assert "Not inside any nebula or supernova remnant." in page


def test_the_tree_route_answers_the_rows_of_a_level_as_json(client):
    answer = client.get(f"/object/{PLANET}/tree")
    assert answer.status_code == 200 and answer.headers["Cache-Control"] == "no-store"
    body = answer.get_json()
    assert body["next"] is None and f'href="/object/{MOON}"' in body["html"] and "tree-leaf" in body["html"]
    assert client.get("/object/ship:1/tree").status_code == 404


def test_a_moon_of_a_free_planet_has_a_page_that_links_to_its_planet(client):
    response = client.get(f"/object/{FREE_MOON}")
    assert response.status_code == 200
    page = response.get_data(as_text=True)
    assert "Lone Moon" in page and "Orbits the free planet Drifter" in page
    assert "Open the free planet Drifter" in page and "/phenomenon/rogue_planet/" in page
    assert 'id="item-properties"' in page
    assert "system</span>" not in page


def test_a_planets_page_lists_the_facilities_it_hosts_and_its_chips(client, monkeypatch):
    facility = {"id": "FE81000A2B-0000005-009", "host_type": "planet", "host_id": "FE81000A2B-0000005-001",
                "name": "Outpost <9>", "kind": "mining-colony",
                "placement": "terrestrial", "orbit_distance_km": None}
    other = dict(facility, host_id="FE81000A2B-0000005-003", name="Elsewhere")
    monkeypatch.setattr(apiclient, "get_system_facilities", lambda db, system_id: [facility, other])
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert "Outpost &lt;9&gt;" in page and "Mining colony" in page and "Elsewhere" not in page
    assert f'<a href="/object/{FACILITY}">Outpost &lt;9&gt;</a>' in page          # UX.183: each is a page of its own
    star = client.get(f"/object/{STAR}").get_data(as_text=True)
    assert "Outpost" not in star and 'class="item-chips"' not in star


def test_a_life_worlds_page_names_its_dominant_species(client, monkeypatch):
    """UX.176: the dominant species of a life world is on the world's own page; a lifeless body never looks one up."""
    detail = _detail()
    detail["planets"][0]["life_chemical"] = "chlorophyll"
    asked = []

    def get_planet_species(db, planet_id):
        asked.append(planet_id)
        return {"id": 3, "name": "Vel<ar>"}

    monkeypatch.setattr(apiclient, "get_system", lambda db, system_id: detail)
    monkeypatch.setattr(apiclient, "get_planet_species", get_planet_species)
    monkeypatch.setattr(item_pages, "population_status", lambda: {"species": True})
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert 'Dominant species: <a href="/species/3">Vel&lt;ar&gt;</a>' in page
    assert asked == ["FE81000A2B-0000005-001"]
    assert "Dominant species" not in client.get(f"/object/{STAR}").get_data(as_text=True)
    assert asked == ["FE81000A2B-0000005-001"]


def test_no_link_to_a_body_goes_to_an_anchor_on_the_systems_page():
    """UX.176: a body has a page of its own; the system page's row anchors, their script and their helper are gone."""
    import pathlib
    import re

    root = pathlib.Path(__file__).resolve().parents[1]
    sources = [path for folder in ("planetgen/web", "html/static") for path in (root / folder).rglob("*")
               if path.suffix in (".py", ".js", ".html", ".css") and "__pycache__" not in path.parts]
    anchor = re.compile(r"body_anchor|bodyanchor|body-targeted|[\"'/]system/[^\"'\s]*#(?:star|planet|moon|belt|comet)-")
    found = [f"{path.relative_to(root)}: {match.group(0)}" for path in sources
             for match in [anchor.search(path.read_text(encoding="utf-8"))] if match]
    assert not found, found


def test_a_page_lists_a_row_only_when_it_holds_something_or_its_absence_matters():
    """UX.188: none and 0 are left out of a page unless the row is one a captain wants to know is absent."""
    from planetgen.web.item_pages import glance_view, sheet_context

    def row(key, value, absent, show):
        return {"key": key, "label": key, "value": [str(value)], "text": str(value), "stored": value is not None,
                "absent": absent, "show_absent": show}

    sheet = {"kind": "planet", "groups": [
        {"key": "atmosphere", "label": "Atmosphere", "rows": [
            row("atmospheric_pressure_pa", 0, True, True), row("p_h2s_kpa", 0, True, False),
            row("p_co2_kpa", 3, False, False)]},
        {"key": "climate", "label": "Climate", "rows": [row("surface_temperature_k", 288, False, False)]}]}
    shown = [group_row["label"] for group in sheet_context(sheet)["sheet"] for group_row in group["rows"]]
    assert shown == ["atmospheric_pressure_pa", "p_co2_kpa", "surface_temperature_k"]
    assert sheet_context(sheet)["hidden_rows"] == 1
    assert [r["label"] for r in glance_view(sheet)] == ["surface_temperature_k", "atmospheric_pressure_pa"]
