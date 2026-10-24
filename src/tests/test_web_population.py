# tests/test_web_population.py

"""
The population pages (`web/population_pages.py`, POP.1 to POP.4): Species,
a species, Polities and a polity, "Dominant species" on a life world's
planet row and "Territory of ..." on an owned system -- and all of them
hidden until a population pass has made species.

Most tests fake the population calls on top of `test_web_system_phen`'s
fake data layer; the last runs a real population pass.
"""

import re

import markupsafe
import pytest

from planetgen.web.app import create_app
from planetgen.api.config import Config

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.web.population_pages import format_years  # noqa: E402
from planetgen.db import store  # noqa: E402
from planetgen.population import model
from tests.test_population import _set_age, galaxy  # noqa: E402,F401
from tests.test_web_facilities import _planet as _facility_planet  # noqa: E402
from tests.test_web_system_phen import app, client, fake  # noqa: E402,F401
from tests.publicids import pid

SPECIES = {
    "id": 3, "name": "Vel<ar>", "homeworld_planet_id": 11, "homeworld_name": "Home", "star_system_id": 5,
    "system_name": "Kepler", "life_chemical": "chlorophyll", "life_stage": "technological_civilization",
    "build": "gracile", "climate": "temperate", "size": "medium", "civilization_age_years": 3.4e6,
    "era": "ancient", "spacefaring": True, "polity_id": 2, "polity_name": "Velar Concord",
    "energy_index": 5, "materials_index": 4, "information_index": 6, "medical_index": 3, "propulsion_index": 5,
    "defense_index": 2, "tech_level": 4.35, "tech_source": "generated", "tech_band": "Fusion",
    "limiting_domain": "materials", "lopsided": True,
}
POLITY = {
    "id": 2, "name": "Velar Concord", "government": "Concord", "color": "#3366cc", "reach_ly": 80.0,
    "species_id": 3, "species_name": "Vel<ar>", "era": "ancient", "civilization_age_years": 3.4e6,
    "capital_system_id": 5, "capital_name": "Kepler", "system_count": 2,
}


class FakePopulation:
    def __init__(self, status):
        self.status = status
        self.calls = []
        self.asked = []

    def get_population_status(self, db):
        self.calls.append("status")
        return dict(self.status)

    def get_species_list(self, db, spacefaring=None, limit=None, offset=None, sort=None, descending=False,
                         eras=(), facets=False):
        self.calls.append(("species", spacefaring))
        self.asked.append(("species", sort, descending, list(eras)))
        items = [SPECIES] if spacefaring in (None, True) and (not eras or SPECIES["era"] in eras) else []
        body = {"items": items, "total": len(items), "limit": limit, "offset": offset}
        if facets:
            body["facets"] = {"spacefaring": [{"value": "yes", "count": 1}, {"value": "no", "count": 4}],
                              "era": [{"value": "ancient", "count": 1}, {"value": "young_age", "count": 4}]}
        return body

    def get_species(self, db, species_id):
        if species_id != SPECIES["id"]:
            raise apiclient.NotFoundError("no such species")
        return SPECIES

    def get_polities(self, db, limit=None, offset=None, sort=None, descending=False, governments=(), eras=(),
                     facets=False):
        self.asked.append(("polities", sort, descending, list(governments), list(eras)))
        body = {"items": [POLITY], "total": 1, "limit": limit, "offset": offset}
        if facets:
            body["facets"] = {"government": [{"value": "Concord", "count": 1}],
                              "era": [{"value": "ancient", "count": 1}]}
        return body

    def get_polity(self, db, polity_id, limit=None, offset=None, sort=None, descending=False):
        if polity_id != POLITY["id"]:
            raise apiclient.NotFoundError("no such polity")
        self.asked.append(("polity", sort, descending))
        systems = [{"id": 5, "name": "Kepler", "distance_ly": 0.0}, {"id": 6, "name": "Far<b>", "distance_ly": 7.25}]
        if descending:
            systems.reverse()
        return {**POLITY, "systems": systems[offset:offset + limit], "limit": limit, "offset": offset}

    def get_planet_species(self, db, planet_id):
        self.calls.append(("planet_species", planet_id))
        return SPECIES if planet_id == SPECIES["homeworld_planet_id"] else None

    habitable = []

    def get_system_habitable_species(self, db, system_id):
        self.calls.append(("habitable", system_id))
        return self.habitable

    def get_system_owner(self, db, system_id):
        self.calls.append(("owner", system_id))
        return {"polity_id": 2, "polity_name": "Velar Concord", "color": "#3366cc", "distance_ly": 0.0}


ALL = {"generated": True, "species": True, "polities": True, "territories": True}
NONE = dict(apiclient.POPULATION_NONE)


def _install(monkeypatch, status):
    data = FakePopulation(status)
    for name in ("get_population_status", "get_species_list", "get_species", "get_polities", "get_polity",
                 "get_planet_species", "get_system_owner", "get_system_habitable_species"):
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    return data


def _nav(html):
    return re.search(r'<nav class="site-sections" aria-label="Main">.*?</nav>', html, re.S).group(0)


def _planet(planet_id, name, life=None):
    planet = _facility_planet(planet_id, name, "t", float(planet_id))
    planet["life_chemical"] = life
    return planet


def test_format_years():
    assert format_years(None) == ""
    assert format_years(342) == "about 340 years"
    assert format_years(12_345) == "about 12,300 years"
    assert format_years(3.4e6) == "about 3.4 million years"
    assert format_years(2.2e9) == "about 2.2 billion years"


def test_without_population_nothing_shows(monkeypatch, client, fake):
    data = _install(monkeypatch, NONE)
    fake.system["planets"] = [_planet(11, "Home", life="chlorophyll")]
    page = client.get("/system/5").get_data(as_text=True)
    assert "Species" not in _nav(page)
    assert "Dominant species" not in page and "Territory of" not in page
    assert not [call for call in data.calls if call != "status"]  # no per-planet or owner lookups
    for url in ("/species/3", "/polities/2"):
        assert client.get(url).status_code == 404
    assert "peoples" not in client.get("/classes").get_data(as_text=True)


@pytest.mark.parametrize("url, kind", [("/species", "species"), ("/polities", "polities")])
def test_the_lists_before_the_population_pass_are_empty_pages_not_404s(monkeypatch, client, fake, url, kind):
    """UX.133: the lists answer 200 with a line that says the pass has not run; no row is looked up."""
    data = _install(monkeypatch, NONE)
    response = client.get(url)
    assert response.status_code == 200
    page = response.get_data(as_text=True)
    assert f"The population pass has not run on this galaxy yet, so it has no {kind}." in page
    assert "datatable" not in page.split("<main", 1)[-1]
    assert not [call for call in data.calls if call != "status"]
    assert "Species" not in _nav(page)


def test_status_failure_hides_the_pages(monkeypatch, client, fake):
    def boom(db):
        raise apiclient.TransportUnreachable("down")
    monkeypatch.setattr(apiclient, "get_population_status", boom)
    page = client.get("/classes").get_data(as_text=True)
    assert "Species" not in _nav(page)


def test_species_without_polities(monkeypatch, client, fake):
    _install(monkeypatch, {"generated": True, "species": True, "polities": False, "territories": False})
    page = client.get("/species").get_data(as_text=True)
    assert "Species" in _nav(page) and "Polities</a>" not in page
    polities = client.get("/polities")
    assert polities.status_code == 200 and "No species in this galaxy has founded a polity." in polities.get_data(as_text=True)
    assert "Territory of" not in client.get("/system/5").get_data(as_text=True)


def test_species_list_and_filter(monkeypatch, client, fake):
    data = _install(monkeypatch, ALL)
    page = client.get("/species").get_data(as_text=True)
    assert '<a href="/species/3">Vel&lt;ar&gt;</a>' in page
    assert '<a href="/system/5">Kepler</a>' in page and '<a href="/polities/2">Velar Concord</a>' in page
    assert 'href="/polities">Polities</a>' in page
    assert '<span class="datatable-option-label">Not spacefaring</span> <span class="datatable-count">(4)</span>' in page
    filtered = client.get("/species?spacefaring=no").get_data(as_text=True)
    assert ("species", False) in data.calls
    assert 'name="spacefaring" value="no" checked' in filtered
    assert "<em>None</em>" in filtered
    # Both choices, or a stray one, filter nothing.
    client.get("/species?spacefaring=yes&spacefaring=no")
    assert data.calls[-1] == ("species", None)
    client.get("/species?spacefaring=bogus")
    assert data.calls[-1] == ("species", None)


def test_species_table_sorts_and_filters_by_era(monkeypatch, client, fake):
    data = _install(monkeypatch, ALL)
    page = client.get("/species?species_sort=polity&species_order=desc&era=young_age").get_data(as_text=True)
    assert ("species", "polity", True, ["young_age"]) in data.asked
    assert '<th scope="col" data-col="polity" aria-sort="descending">' in page
    assert '<span class="datatable-option-label">Young age</span>' in page


def test_population_tables_serve_json(monkeypatch, client, fake):
    data = _install(monkeypatch, ALL)
    body = client.get("/table/species?sort=era&facets=1").get_json()
    assert body["rows"][0][0] == {"text": "Vel<ar>", "href": "/species/3"}
    assert body["rows"][0][4] == {"text": "4.35 Fusion"}
    assert body["rows"][0][6] == {"text": "Velar Concord", "href": "/polities/2"}
    assert {o["value"]: o["label"] for o in body["facets"]["spacefaring"]} == {
        "yes": "Spacefaring", "no": "Not spacefaring"}
    body = client.get("/table/polities?government=Concord&era=ancient&facets=1").get_json()
    assert body["rows"][0][0] == {"text": "Velar Concord", "href": "/polities/2", "swatch": "#3366cc"}
    assert body["rows"][0][4] == {"text": "2"} and body["rows"][0][5]["text"] == "80.0 ly"
    assert ("polities", "name", False, ["Concord"], ["ancient"]) in data.asked
    body = client.get("/table/polity-systems?polity=2&sort=distance&order=desc").get_json()
    assert [row[0]["text"] for row in body["rows"]] == ["Far<b>", "Kepler"]
    assert body["total"] == 2
    assert client.get("/table/polity-systems?polity=9").status_code == 404


def test_species_page(monkeypatch, client, fake):
    _install(monkeypatch, ALL)
    page = client.get("/species/3").get_data(as_text=True)
    assert "<h1>Vel&lt;ar&gt;</h1>" in page or "Vel&lt;ar&gt;</h1>" in page
    assert "Technological civilization" in page and "about 3.4 million years" in page
    assert '<a href="/system/5">Kepler</a>' in page and '<a href="/polities/2">Velar Concord</a>' in page
    # POP.9: the level with its band, the six indices, what limits it and the lopsided flag.
    assert "4.35 Fusion" in page and "Energy index" in page and "Defense index" in page
    assert "5 of 7" in page and "2 of 7" in page
    assert re.search(r"Limited by.*?Materials", page, re.S) and "Lopsided" in page
    assert client.get("/species/4").status_code == 404


def test_polity_pages(monkeypatch, client, fake):
    _install(monkeypatch, ALL)
    listing = client.get("/polities").get_data(as_text=True)
    assert '<a href="/polities/2">Velar Concord</a>' in listing and "80.0 ly" in listing
    assert "Territories on the Galaxy Map" in listing
    page = client.get("/polities/2").get_data(as_text=True)
    assert '<a href="/system/6">Far&lt;b&gt;</a>' in page and "7.3 ly" in page
    assert page.index("/system/5") < page.index("/system/6")
    assert 'data-col="distance" aria-sort="ascending"' in page
    reversed_page = client.get("/polities/2?systems_sort=distance&systems_order=desc").get_data(as_text=True)
    assert reversed_page.index("/system/6") < reversed_page.rindex("/system/5")
    assert client.get("/polities/9").status_code == 404


def test_system_page_names_species_and_territory(monkeypatch, client, fake):
    data = _install(monkeypatch, ALL)
    fake.system["planets"] = [_planet(11, "Home", life="chlorophyll"), _planet(12, "Rock")]
    page = client.get("/system/5").get_data(as_text=True)
    assert 'Dominant species: <a href="/species/3">Vel&lt;ar&gt;</a>' in page
    assert "Species: Vel&lt;ar&gt;" in page
    assert 'Territory of <a href="/polities/2">Velar Concord</a>' in page
    assert ("planet_species", 12) not in data.calls  # lifeless planets are never looked up


def test_classes_link_the_population_pages(monkeypatch, client, fake):
    _install(monkeypatch, ALL)
    page = client.get("/classes").get_data(as_text=True)
    assert '<a href="/species">Species</a>' in page and '<a href="/polities">Polities</a>' in page
    assert 'href="/species"' in _nav(page)


# --- A real database ---------------------------------------------------------------------

@pytest.fixture
def real_client(mysql_config):
    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        WEB_DATABASE = ""
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    monkeypatch_env = pytest.MonkeyPatch()
    monkeypatch_env.setenv("PLANETGEN_PAGE_CACHE", "off")
    application = create_app(RealConfig)
    application.testing = True
    yield application.test_client()
    monkeypatch_env.undo()


def test_real_population(mysql_config, galaxy, real_client):  # noqa: F811
    before = real_client.get("/classes").get_data(as_text=True)
    assert 'href="/species"' not in _nav(before)
    for url in ("/species", "/polities"):   # UX.133: an empty page that says so, not a 404
        empty = real_client.get(url)
        assert empty.status_code == 200 and "The population pass has not run" in empty.get_data(as_text=True)
    assert real_client.get("/species/1").status_code == 404

    conn = store.get_connection(mysql_config)
    try:
        model.run_pass(conn)
        with conn:
            _set_age(conn, galaxy["a"], 1e6)
        model.run_pass(conn)
        homeworld = conn.execute("SELECT homeworld_planet_id FROM species WHERE star_system_id = ?",
                                 (galaxy["a"],)).fetchone()["homeworld_planet_id"]
        name = conn.execute("SELECT name FROM species WHERE homeworld_planet_id = ?", (homeworld,)).fetchone()["name"]
    finally:
        conn.close()

    listing = real_client.get("/species").get_data(as_text=True)
    # The page escapes the name; some draws have an apostrophe ("Epteyn'Ska").
    assert 'href="/species"' in _nav(listing) and str(markupsafe.escape(name)) in listing
    assert real_client.get("/polities").status_code == 200
    system = real_client.get(f"/system/{pid('system', galaxy['a'])}").get_data(as_text=True)
    assert "Dominant species:" in system and "Territory of" in system


LIMITS = [
    {"name": "comfort_high_c", "unit": "C", "meaning": "greatest temperature with no gear (Blue)", "human": 31.0,
     "value": 55.0, "changed": True},
    {"name": "pressure_no_gear_kpa", "unit": "kPa", "meaning": "least pressure with no gear at all (Blue)",
     "human": 50.0, "value": 50.0, "changed": False},
]


def test_species_page_lists_the_environment_limits(monkeypatch, client, fake):
    data = _install(monkeypatch, ALL)
    monkeypatch.setattr(data, "get_species", lambda db, species_id: {**SPECIES, "limits": LIMITS})
    monkeypatch.setattr(apiclient, "get_species", data.get_species)
    page = client.get("/species/3").get_data(as_text=True)
    assert "Environment Limits" in page
    assert re.search(r"Comfort high</th><td>55 \u00b0C</td><td>31 \u00b0C</td><td>Yes</td>", page)
    assert re.search(r"Pressure no gear</th><td>50 kPa</td><td>50 kPa</td><td>No</td>", page)


def _habitable(name, ids):
    domains = {d: {"score": 1.0, "tier": "Blue"} for d in ("pressure", "temperature", "air", "water", "radiation")}
    domains["suitability"] = {"score": 0.6, "tier": "Yellow"}
    return {"id": ids, "name": name, "homeworld_planet_id": 11, "tech_level": 4.0, "equipment": "ideal",
            "unaided": True, "equipped": False, "failing": [], "phi4": 0.8, "domains": domains}


def test_system_page_lists_the_species_that_could_live_on_each_world(monkeypatch, client, fake):
    data = _install(monkeypatch, ALL)
    data.habitable = [{"kind": "planet", "id": 11, "name": "Home<b>", "unaided": [_habitable("Vel<ar>", 3)],
                       "equipped": [_habitable("Other", 4)]}]
    page = client.get("/system/5").get_data(as_text=True)
    assert "Species that could live here" in page and "Home&lt;b&gt;" in page
    assert '<a href="/species/3">Vel&lt;ar&gt;</a>' in page and '<a href="/species/4">Other</a>' in page
    assert "Suitability: Yellow" in page and "Pressure: Blue" in page
    assert "<caption>Unaided</caption>" in page and "<caption>With equipment</caption>" in page


def test_system_page_has_no_species_panel_when_none_could_live_there(monkeypatch, client, fake):
    _install(monkeypatch, ALL)
    assert "Species that could live here" not in client.get("/system/5").get_data(as_text=True)
