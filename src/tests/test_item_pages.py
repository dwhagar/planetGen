# tests/test_item_pages.py

"""
UX.143: every kind of item a reference can name has a page, from one registry, and a star, planet, moon, belt or
comet gets the item page at `/object/<ref>`.
"""

import re

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
    monkeypatch.setattr(apiclient, "get_life_air", lambda db, kind, body_id: {"id": body_id, "life_class_id": None, "life_air": None})
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
    assert 'class="summary inline-links"' in page and "Kepler b is a world" in page      # UX.220: the summary is written from the data
    assert "A rocky world." not in page                  # and the stored class description no longer opens the page
    assert 'id="item-properties"' in page and "Mass" in page and "Ocean depth" not in page   # the sheet; a value never stored is not listed (UX.199)
    assert 'href="/sector/4"' in page and 'href="/system/9"' in page          # the places it is in
    assert f'href="/object/{STAR}"' in page               # the star it orbits is on the chain
    assert 'data-tree-root' in page and f"/object/{MOON}" in page            # the bodies that belong to it, each a page
    assert page.index('id="item-facts"') < page.index('id="location"') < page.index('id="item-properties"') \
        < page.index('id="item-members"')      # UX.215: the Location row is the first row of the At a glance table
    assert "#planet-" not in page, "no more anchors on the system's page"


def test_the_sheet_shows_a_quantity_in_the_readers_units_and_the_tree_row_says_what_it_is(client):
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert 'data-q="mass"' in page                       # a quantity the unit script can rewrite
    assert 'class="tree-name"' in page and "Class" in page


def test_a_scored_planet_shows_its_symbols_under_the_title_and_other_kinds_do_not(client):
    """UX.100: the item page of a scored planet or moon carries the seven PHI-6 symbols; a star or unscored body none."""
    body = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert body.count("hab-symbol ") == 7 and 'aria-label="Radiation: Red"' in body
    assert body.index("About Kepler b") < body.index("hab-symbols")        # UX.215: the symbols are a row of the table inside About
    # UX.188: the symbols are the first row of the At a glance table, above the property sheet
    assert body.index('id="item-glance-heading"') < body.index("hab-symbols") < body.index('id="item-properties"')
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
    assert 'id="item-properties"' in page and "Mass" in page            # the sheet, from the fake
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
    for url in ('href="/sector/4"', 'href="/system/9"', f'href="/object/{PLANET}"'):
        assert url in page
    assert 'href="/galaxy"' not in page.split('id="location"')[1].split("</section>")[0]    # UX.200: no Galaxy link
    assert "Not inside any nebula" not in page


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
    detail["planets"][0]["life_class"] = 1
    detail["planets"][0]["life_class_id"] = "1a"
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


def test_the_glance_panel_flags_a_very_high_swing_of_starlight_only(monkeypatch):
    """GEN.1135: Class Q is retired; the panel says "very high" at 0.75 peak to peak or more and shows nothing below it."""
    from planetgen.describe import sheet as sheet_registry
    from planetgen.web.item_pages import glance_view

    def sheet_with(swing):
        return {"kind": "planet", "groups": [{"key": "climate", "label": "Climate", "rows": [
            {"key": "flux_variation", "label": "Seasonal variation of that energy, peak to peak", "value": [f"{swing * 100:.0f} %"],
             "text": f"{swing * 100:.0f} %", "raw": swing, "stored": True, "absent": swing == 0, "show_absent": False}]}]}

    assert "flux_variation" in sheet_registry.GLANCE["planet"] and "flux_variation" in sheet_registry.GLANCE["moon"]
    assert glance_view(sheet_with(0.74)) == []
    assert glance_view(sheet_with(0.0)) == []
    flagged = glance_view(sheet_with(0.9))
    assert [r["label"] for r in flagged] == [sheet_registry.SWING_FLAG_LABEL]
    assert str(flagged[0]["html"]).startswith("Very high")
    assert glance_view(sheet_with(sheet_registry.SWING_FLAG_MIN))


def test_the_symbols_on_a_planets_page_go_to_its_own_habitability_group_and_that_group_links_the_help_page(client):
    """UX.197: the icons link to this page's Habitability section, whose help link goes on to the levels page."""
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert f'href="/object/{PLANET}#habitability"' in page


def test_a_barren_airless_dry_world_says_so_in_sentences_and_leaves_the_sections_out():
    """UX.203: no life, no atmosphere and no water become sentences and a "None", not tables of zeros."""
    from planetgen.describe import sheet as registry
    from planetgen.web.item_pages import sheet_context

    barren = {"life_microbial": 0, "life_flora": 0, "life_fauna": 0, "life_sentient": 0, "atmospheric_pressure_pa": 0,
              "hydrosphere": "dry", "water_mass_fraction": 0, "ocean_fraction": 0}
    assert registry.missing("planet", barren) == ["life", "atmosphere", "water"]
    assert registry.missing("star", barren) == []
    # a trace of air, of water or of life, and a value that is not stored, never count as none
    assert registry.missing("moon", {**barren, "atmospheric_pressure_pa": 5.0}) == ["life", "water"]
    assert registry.missing("moon", {**barren, "water_mass_fraction": 1e-3}) == ["life", "atmosphere"]
    assert registry.missing("moon", {**barren, "hydrosphere": "vapour"}) == ["life", "atmosphere"]
    assert registry.missing("moon", {**barren, "life_fauna": 1}) == ["atmosphere", "water"]
    assert registry.missing("planet", {**barren, "life_fauna": None, "atmospheric_pressure_pa": None,
                                        "hydrosphere": None}) == []

    def row(key):
        return {"key": key, "label": key, "value": ["0"], "text": "0", "stored": True, "absent": True, "show_absent": True}
    groups = [{"key": key, "label": key.title(), "rows": [row(key + "_row")]} for key in ("atmosphere", "climate", "water", "life")]
    context = sheet_context({"kind": "planet", "groups": groups, "missing": ["life", "atmosphere", "water"]})
    assert [group["key"] for group in context["sheet"]] == ["atmosphere", "climate"]
    assert [r["html"] for r in context["sheet"][0]["rows"]] == ["None"]
    assert context["absences"] == [registry.MISSING_SENTENCES[key] for key in ("life", "atmosphere", "water")]
    assert [group["key"] for group in sheet_context({"kind": "planet", "groups": groups})["sheet"]] == [
        "atmosphere", "climate", "water", "life"]


def _day_sheet(**shown):
    """A sheet whose rotation group holds the stored day and night rows named in `shown` (a key to its text and number)."""
    rows = [{"key": key, "label": key, "value": [text], "text": text, "raw": number, "stored": True, "absent": number == 0}
            for key, (text, number) in shown.items()]
    return {"kind": "planet", "groups": [{"key": "rotation", "label": "Rotation", "rows": rows}]}


def test_the_summary_and_the_rotation_group_give_the_day_and_night_shares(client, monkeypatch):
    """UX.205: an Earth-like world gets the share of the time its star is up, a tidally locked one says day side,
    terminator and night side, and a moon in its planet's shadow gets the extra night with the longest eclipse."""
    earth = _day_sheet(day_share_equator=("50.2 %", 0.502), day_share_mean=("50.1 %", 0.501), solar_day_hours=("24 h", 24.0),
                       daylight_hours=("12 h", 12.0))
    locked = _day_sheet(day_share_mean=("50 %", 0.5))
    shadowed = _day_sheet(**{**{"day_share_equator": ("47 %", 0.47), "day_share_mean": ("47 %", 0.47),
                                "solar_day_hours": ("42 h", 42.0), "daylight_hours": ("20 h", 20.0)},
                             "parent_shadow_share": ("4.6 %", 0.046), "longest_eclipse_hours": ("2.3 h", 2.3)})
    for ref, sheet in ((PLANET, earth), (PLANET, locked), (MOON, shadowed)):
        monkeypatch.setattr(apiclient, "get_object_sheet", lambda db, r, sheet=sheet: {**sheet, "kind": ref.partition(":")[0]})
        page = client.get(f"/object/{ref}").get_data(as_text=True)
        summary = page[page.index('id="item-summary-heading"'):page.index('id="item-properties"')]
        if sheet is earth:
            assert "Its star is up 50.2 % of the time at the equator, 12 h of daylight in each solar day of 24 h" in summary
            assert "keeps one face" not in summary and "shadow" not in summary
        elif sheet is locked:
            assert "keeps one face to its star" in summary and "terminator" in summary and "solar day" not in summary
        else:
            assert "in the shadow 4.6 % of the time, the longest eclipse lasting 2.3 h" in summary
        assert "longest_eclipse_hours" in page if sheet is shadowed else True     # and the row is in the rotation group


def test_an_object_with_no_day_and_night_shares_says_nothing_of_them(client):
    page = client.get(f"/object/{STAR}").get_data(as_text=True)
    assert "star is up" not in page and "one face" not in page


def test_the_breadcrumb_calls_the_system_link_system_and_the_star_link_by_the_stars_name(client):
    """The system link of the trail reads "System" and the star's its own name, so a shared name is not shown twice."""
    page = client.get(f"/object/{MOON}").get_data(as_text=True)
    trail = re.search(r'<nav[^>]*aria-label="Breadcrumb".*?</nav>', page, re.S)
    assert trail, "no breadcrumb trail"
    assert re.search(r'<a [^>]*href="/system/9"[^>]*>System</a>', trail.group(0))


def test_the_properties_follow_the_about_block_right_after_at_a_glance(client):
    """Boss, 10-24: header and render, At a glance, About, then the property list."""
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    order = [page.index(f'id="{section}"') for section in ("item-facts", "item-glance-heading", "item-properties")]
    assert order == sorted(order)


def test_a_stored_key_reads_as_words_on_the_sheet():
    """A text value such as the life stage technological_civilization shows no underscore."""
    from planetgen.describe import sheet as registry

    spec = next(spec for spec in registry.SPECS["planet"] if spec.key == "life_stage")
    assert registry._value(spec, "planets", "technological_civilization") == ["Technological civilization"]


def test_every_row_the_atmosphere_section_lists_has_a_label_and_a_value():
    """UX.198: gases the air does not hold (partial pressure 0) are not listed, and the rest are never blank."""
    from planetgen.web.item_pages import sheet_context

    def row(key, label, value, absent):
        return {"key": key, "label": label, "value": [value], "text": value, "stored": True, "absent": absent,
                "show_absent": False}
    gases = [row("p_n2", "Nitrogen partial pressure", "78 kPa", False)] + [
        row(f"p_gas{i}", f"Gas {i} partial pressure", "0 Pa", True) for i in range(14)]
    context = sheet_context({"kind": "planet", "groups": [{"key": "atmosphere", "label": "Atmosphere", "rows": gases}]})
    listed = context["sheet"][0]["rows"]
    assert [r["label"] for r in listed] == ["Nitrogen partial pressure"]
    assert all(r["label"] and str(r["html"]) for r in listed) and context["hidden_rows"] == 14


def _rich_sheet(ref):
    q = model.qty
    row = lambda key, label, value, raw=None: {"key": key, "label": label, "value": value, "text": label, "stored": True, "raw": raw}
    return {"kind": "planet", "groups": [
        {"key": "identity", "label": "Identity", "rows": [
            row("planet_class", "Class", ["M"]),
            row("zone", "Zone", ["Ecosphere / Habitable Zone, ", q("length", 1.4e8, "km"), " to ", q("length", 2e8, "km"), " from its star"])]},
        {"key": "atmosphere", "label": "Atmosphere", "rows": [
            row("atmospheric_pressure_pa", "Surface pressure", [q("pressure", 101325, "Pa")], 101325.0),
            row("p_N2", "Nitrogen partial pressure", [q("pressure", 79000, "Pa")], 79000.0),
            row("p_O2", "Oxygen partial pressure", [q("pressure", 22325, "Pa")], 22325.0)]},
        {"key": "composition", "label": "Composition", "rows": [
            row("metal_mass_fraction", "Metal", ["0.3"], 0.3), row("rock_mass_fraction", "Rock", ["0.7"], 0.7),
            row("mantle_redox", "Mantle oxidation", ["Reduced"])]},
        {"key": "habitability", "label": "Habitability", "rows": [row("phi4", "Habitability score (PHI-6)", ["0.8"], 0.8)]},
        {"key": "physical", "label": "Physical", "rows": [row("mass_kg", "Mass", [q("mass", 5.97e24, "kg")])]}]}


def _rich_document():
    document = _document()
    planet = document["sections"][0]["children"][0]
    planet["blocks"] += [
        model.ratings({"name": "PHI-6", "score": 0.8, "tier": "Green"}, [], []),
        model.facts("Radiation", [model.fact("dose", "Total dose", "1 mSv")]),
        model.prose("life", "Life based on carbon lives on this world."), model.prose("flavour", "Sensors show a green haze.")]
    return document


def test_a_planet_page_is_laid_out_about_habitability_atmosphere_composition_then_properties(client, monkeypatch):
    """UX.215 to UX.222: About (a table with the Location row first, a summary, Life, Civilization), then PHI-6, the Atmosphere
    and Composition tables, each a heading, then the rest of the properties."""
    monkeypatch.setattr(apiclient, "get_object_sheet", lambda db, ref: _rich_sheet(ref))
    monkeypatch.setattr(apiclient, "get_system_description", lambda db, system_id: _rich_document())
    detail = _detail()
    detail["planets"][0].update(planet_class="M", zone="e", body_type="t", phi4=0.8, tier_suitability=1, life_microbial=2, life_flora=0,
                                life_fauna=0, life_sentient=0, atmospheric_pressure_pa=101325.0, inhabited=True,
                                equipment={"required": ["a gas-mix mask"], "items": []})
    monkeypatch.setattr(apiclient, "get_system", lambda db, system_id: detail)
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    headings = re.findall(r'<h([23])[^>]*>(.*?)</h\1>', page)
    names = [re.sub(r"<[^>]+>", "", text) for _level, text in headings]
    wanted = ["Sensor view", "About Kepler b", "At a glance", "Summary", "Life", "Civilization and politics", "Habitability (PHI-6)",
              "Atmosphere", "Composition", "Properties"]
    positions = [names.index(name) for name in wanted]
    assert positions == sorted(positions), names
    assert "Rendered example" not in page
    assert page.index('id="location"') < page.index("hab-symbols") < page.index('class="summary inline-links"')
    assert page.index("Life based on carbon") < page.index('id="habitability"') < page.index('id="sheet-atmosphere"') \
        < page.index('id="sheet-composition"') < page.index('id="item-properties"')
    assert 'id="location"' in page.split('id="item-facts"')[1].split("</section>")[0]     # the Location row is inside About
    assert 'href="/classes/planet/M"' in page and "Ecosphere / Habitable Zone, " in page and " from its star" in page
    assert "Mole fraction" in page and "Partial pressure" in page and "Share by mass" in page
    assert "Metal (iron and nickel core)" in page and "70 percent of its mass" in page
    assert "Dominant species" not in page and "A technological civilization or a colony lives here." in page


def test_gen1108_a_planets_page_has_the_climate_by_latitude_table_between_the_atmosphere_and_the_composition(client, monkeypatch):
    from planetgen.describe import sheet as sheet_registry
    from planetgen.physics import regional_temperature
    earth = regional_temperature.compute(dict(surface_temperature_k=288.0, atmospheric_pressure_pa=101325.0, gravity_g=1.0, rotation_period_hours=24.0,
                                              axial_tilt_deg=23.4, period_years=1.0))
    monkeypatch.setattr(apiclient, "get_object_sheet", lambda db, ref: {**_rich_sheet(ref), "bands": sheet_registry.bands_table(earth)})
    monkeypatch.setattr(apiclient, "get_system_description", lambda db, system_id: _rich_document())
    detail = _detail()
    detail["planets"][0].update(planet_class="M", zone="e", body_type="t", phi4=0.8, tier_suitability=1, life_microbial=2, life_flora=0,
                                life_fauna=0, life_sentient=0, atmospheric_pressure_pa=101325.0, inhabited=True,
                                equipment={"required": ["a gas-mix mask"], "items": []})
    monkeypatch.setattr(apiclient, "get_system", lambda db, system_id: detail)
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert page.index('id="sheet-atmosphere"') < page.index('id="sheet-bands"') < page.index('id="sheet-composition"')
    table = page.split('id="sheet-bands"')[1].split("</section>")[0]
    assert "Climate by latitude" in table and "North polar" in table and "South polar" in table and "Equator" in table
    assert "Annual mean temperature" in table and "Longest day of the year" in table and "Polar day and night" in table
    assert "no day and night range stored" not in page


def test_a_planet_with_no_life_or_people_has_no_life_or_civilization_section(client):
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert 'id="life-paragraph-heading"' not in page and 'id="civilization-heading"' not in page


def test_a_planet_with_life_lists_the_classes_it_supports_and_what_the_class_lives_on(client, monkeypatch):
    detail = _detail()
    detail["planets"][0].update(planet_class="M", zone="e", body_type="t", phi4=0.8, tier_suitability=1, life_microbial=2, life_flora=3,
                                life_fauna=1, life_sentient=0, life_class=1, life_class_variant="a", life_class_id="1a",
                                life_class_best="1a", life_class_best_score=7, life_class_scores={"1a": 7, "2": 4},
                                life_stage="multicellularity", life_animal_grade=2)
    monkeypatch.setattr(apiclient, "get_system", lambda db, system_id: detail)
    found = apiclient.get_object
    monkeypatch.setattr(apiclient, "get_object", lambda db, ref: {**found(db, ref), "life_microbial": 2, "life_flora": 3,
                                                                  "life_fauna": 1, "life_sentient": 0})
    readings = [
        {"kind": "microbial", "label": "Microbes", "score": 2, "abundance": "Rare", "trend": "Stable", "security": "Secure",
         "years_left_gyr": None, "text": "Microbes: rare, stable, secure."},
        {"kind": "flora", "label": "Flora", "score": 3, "abundance": "Occasional", "trend": "Declining", "security": "Precarious",
         "years_left_gyr": 0.4, "text": "Flora: occasional, declining, precarious."},
        {"kind": "fauna", "label": "Fauna", "score": 1, "abundance": "Rare", "trend": "Stable", "security": "Secure",
         "years_left_gyr": None, "text": "Fauna: rare, stable, secure."},
        {"kind": "sentient", "label": "Sentient life", "score": 0, "abundance": "Absent", "trend": None, "security": None,
         "years_left_gyr": None, "text": "Sentient life: absent."}]
    monkeypatch.setattr(apiclient, "get_life_readings", lambda db, kind, body_id: {"id": body_id, "life_class_id": "1a", "readings": readings})
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert re.search(r"Conditions here suit <a [^>]*>Class 1a</a> life \(oxygenic light-eaters\) best, suitability 7 of 9", page)
    assert re.search(r"they also suit <a [^>]*>Class 2</a> \(anoxygenic light-eaters[^<]*, suitability 4 of 9\)\.", page)
    assert re.search(r"Life here is <a [^>]*>Class 1a</a>, at the stage of small mobile animals \(animal grade A2\)\. Conditions here suit", page)
    glance = page.split('class="life-glance"')[1].split("</p>")[0]
    assert "Microbes: rare, stable, secure. Flora: occasional, declining, precarious. Fauna: rare, stable, secure. Conditions suit <a" in glance
    assert "Sentient" not in glance and "of 9" not in glance
    assert page.index('id="item-glance-heading"') < page.index('class="life-glance"') < page.index('id="item-summary-heading"')
    table = page.split('class="life-readings"')[1].split("</table>")[0]
    assert "<th scope=\"row\">Flora</th><td>Occasional (3 of 9)</td><td>Declining</td><td>Precarious</td><td>about 0.4 billion years</td>" in table
    assert "Its life scores are" not in page and "latest milestone" not in page and "evolving at" not in page
    assert '<a class="help-link" href="/help/life/classes/1#class-1a">Class 1a</a> (oxygenic light-eaters) best, then <a class="help-link" href="/help/life/classes/2#class-2">Class 2</a>' in glance
    assert 'href="/help/life">What these words mean</a>' in glance
    for topic in ("abundance", "trend", "security"):
        assert f'<a class="help-link help-mark" href="/help/life/{topic}"' in table
    assert 'href="/help/life/classes/1#class-1a">Class 1a</a>, at the stage of small mobile animals' in page
    assert '/help/life/stages">stages and animal grades</a>' in page and '/help/life/suitability">suitability</a>' in page
    assert page.index('id="life-paragraph-heading"') < page.index("Life here is <a") < page.index('id="life-class-words"')
    block = page.split('id="life-class-words"')[1].split("</details>")[0]
    assert "eats, drinks, breathes and gives off" in page and "<dt>Drinks</dt><dd>water</dd>" in block
    assert re.search(r'<details[^>]*id="life-class-words"(?![^>]*\bopen\b)', page), "the block starts collapsed"
    assert 'class="hab-life-class"' in page


def _planet_with_life(monkeypatch, history):
    detail = _detail()
    detail["planets"][0].update(planet_class="M", zone="e", body_type="t", phi4=0.8, tier_suitability=1, life_microbial=2, life_flora=3,
                                life_fauna=1, life_sentient=0, life_class=1, life_class_variant="a", life_class_id="1a",
                                life_class_best="1a", life_class_best_score=7, life_class_scores={"1a": 7},
                                life_stage="multicellularity", life_animal_grade=2)
    monkeypatch.setattr(apiclient, "get_system", lambda db, system_id: detail)
    asked = []

    def get_life_air(db, kind, body_id):
        asked.append((kind, body_id))
        return {"id": body_id, "life_class_id": "1a", "life_air": history}

    monkeypatch.setattr(apiclient, "get_life_air", get_life_air)
    return asked


def test_a_planet_whose_life_changed_its_air_has_a_collapsed_life_and_air_chart_in_the_life_section(client, monkeypatch):
    from tests.test_life_air_chart import HISTORY
    asked = _planet_with_life(monkeypatch, {**HISTORY, "speculative": {"tier": "T7", "tweak": "L4"}})
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert asked == [("planet", "FE81000A2B-0000005-001")]
    assert re.search(r'<details[^>]*id="life-and-air"(?![^>]*\bopen\b)', page), "the chart starts collapsed"
    block = page.split('id="life-and-air"')[1].split("</details>")[0]
    assert "Life and air through the ages" in page and '<svg class="life-air-chart"' in block and block.count("<polyline") == 3
    assert '<li value="1">' in block and "Speculative chemistry" in block and 'class="life-air-table"' in block and "29 kPa" in block
    assert page.index('id="life-paragraph-heading"') < page.index('id="life-and-air"') < page.index('id="item-properties"')


def test_a_living_planet_whose_air_the_life_did_not_change_has_no_chart(client, monkeypatch):
    asked = _planet_with_life(monkeypatch, None)
    assert 'id="life-and-air"' not in client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert asked == [("planet", "FE81000A2B-0000005-001")]


def test_a_planet_with_no_living_class_does_not_ask_for_a_chart(client, monkeypatch):
    monkeypatch.setattr(apiclient, "get_life_air", lambda *args: pytest.fail("the chart is asked for only for a body with a living class"))
    page = client.get(f"/object/{PLANET}").get_data(as_text=True)
    assert 'id="life-and-air"' not in page


def test_gen1108_the_page_draws_the_seven_bands_as_a_table_of_columns():
    """GEN.1108: the sheet's bands are a table with a column for each band and a row for each figure, a band with no such figure a dash; a body with no bands has no table."""
    from planetgen.describe import sheet as sheet_registry
    from planetgen.physics import latitude_bands, regional_temperature
    earth = regional_temperature.compute(dict(surface_temperature_k=288.0, atmospheric_pressure_pa=101325.0, gravity_g=1.0, rotation_period_hours=24.0,
                                              axial_tilt_deg=23.4, period_years=1.0))
    view = item_pages.bands_view({"bands": sheet_registry.bands_table(earth)})
    assert view["columns"] == [sheet_registry.BAND_LABELS[name] for name in latitude_bands.NAMES]
    assert [row["label"] for row in view["rows"]][:2] == ["Annual mean temperature", "Day and night swing (hottest to coldest)"]
    assert all(len(row["cells"]) == 7 for row in view["rows"])
    moon = regional_temperature.compute(dict(surface_temperature_k=250.0, atmospheric_pressure_pa=0.0, gravity_g=0.165, rotation_period_hours=655.7,
                                             axial_tilt_deg=1.5, period_years=0.0748))
    rings = item_pages.bands_view({"bands": sheet_registry.bands_table(moon)})
    assert "point under its star" in rings["caption"] and [row["label"] for row in rings["rows"]] == ["Annual mean temperature",
                                                                                                      "Day and night swing (hottest to coldest)", "State"]
    assert item_pages.bands_view({"bands": None}) is None and item_pages.bands_view({}) is None
