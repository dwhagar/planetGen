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

SPECIES = {
    "id": 3, "name": "Vel<ar>", "homeworld_planet_id": 11, "homeworld_name": "Home", "star_system_id": 5,
    "system_name": "Kepler", "life_chemical": "chlorophyll", "life_stage": "technological_civilization",
    "build": "gracile", "climate": "temperate", "size": "medium", "civilization_age_years": 3.4e6,
    "era": "ancient", "spacefaring": True, "polity_id": 2, "polity_name": "Velar Concord",
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

    def get_population_status(self, db):
        self.calls.append("status")
        return dict(self.status)

    def get_species_list(self, db, spacefaring=None, limit=None, offset=None):
        self.calls.append(("species", spacefaring))
        items = [SPECIES] if spacefaring in (None, True) else []
        return {"items": items, "total": len(items), "limit": limit, "offset": offset}

    def get_species(self, db, species_id):
        if species_id != SPECIES["id"]:
            raise apiclient.NotFoundError("no such species")
        return SPECIES

    def get_polities(self, db, limit=None, offset=None):
        return {"items": [POLITY], "total": 1, "limit": limit, "offset": offset}

    def get_polity(self, db, polity_id, limit=None, offset=None):
        if polity_id != POLITY["id"]:
            raise apiclient.NotFoundError("no such polity")
        systems = [{"id": 5, "name": "Kepler", "distance_ly": 0.0}, {"id": 6, "name": "Far<b>", "distance_ly": 7.25}]
        return {**POLITY, "systems": systems[offset:offset + limit], "limit": limit, "offset": offset}

    def get_planet_species(self, db, planet_id):
        self.calls.append(("planet_species", planet_id))
        return SPECIES if planet_id == SPECIES["homeworld_planet_id"] else None

    def get_system_owner(self, db, system_id):
        self.calls.append(("owner", system_id))
        return {"polity_id": 2, "polity_name": "Velar Concord", "color": "#3366cc", "distance_ly": 0.0}


ALL = {"generated": True, "species": True, "polities": True, "territories": True}
NONE = dict(apiclient.POPULATION_NONE)


def _install(monkeypatch, status):
    data = FakePopulation(status)
    for name in ("get_population_status", "get_species_list", "get_species", "get_polities", "get_polity",
                 "get_planet_species", "get_system_owner"):
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
    assert format_years(12_345) == "about 1.23 × 10⁴ years"
    assert format_years(3.4e6) == "about 3.4 million years"
    assert format_years(2.2e9) == "about 2.2 billion years"


def test_without_population_nothing_shows(monkeypatch, client, fake):
    data = _install(monkeypatch, NONE)
    fake.system["planets"] = [_planet(11, "Home", life="chlorophyll")]
    page = client.get("/system/5").get_data(as_text=True)
    assert "Species" not in _nav(page)
    assert "Dominant species" not in page and "Territory of" not in page
    assert not [call for call in data.calls if call != "status"]  # no per-planet or owner lookups
    for url in ("/species", "/species/3", "/polities", "/polities/2"):
        assert client.get(url).status_code == 404
    assert "peoples" not in client.get("/classes").get_data(as_text=True)


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
    assert client.get("/polities").status_code == 404
    assert "Territory of" not in client.get("/system/5").get_data(as_text=True)


def test_species_list_and_filter(monkeypatch, client, fake):
    data = _install(monkeypatch, ALL)
    page = client.get("/species").get_data(as_text=True)
    assert '<a href="/species/3">Vel&lt;ar&gt;</a>' in page
    assert '<a href="/system/5">Kepler</a>' in page and '<a href="/polities/2">Velar Concord</a>' in page
    assert 'href="/species" aria-current="page">All</a>' in page
    assert 'href="/polities">Polities</a>' in page
    filtered = client.get("/species?spacefaring=0").get_data(as_text=True)
    assert ("species", False) in data.calls
    assert 'href="/species?spacefaring=0" aria-current="page">Not spacefaring</a>' in filtered
    assert "<em>None</em>" in filtered
    client.get("/species?spacefaring=bogus")
    assert data.calls[-1] == ("species", None)


def test_species_page(monkeypatch, client, fake):
    _install(monkeypatch, ALL)
    page = client.get("/species/3").get_data(as_text=True)
    assert "<h1>Vel&lt;ar&gt;</h1>" in page or "Vel&lt;ar&gt;</h1>" in page
    assert "Technological civilization" in page and "about 3.4 million years" in page
    assert '<a href="/system/5">Kepler</a>' in page and '<a href="/polities/2">Velar Concord</a>' in page
    assert client.get("/species/4").status_code == 404


def test_polity_pages(monkeypatch, client, fake):
    _install(monkeypatch, ALL)
    listing = client.get("/polities").get_data(as_text=True)
    assert '<a href="/polities/2">Velar Concord</a>' in listing and "80.0 ly" in listing
    assert "Territories on the Galaxy Map" in listing
    page = client.get("/polities/2").get_data(as_text=True)
    assert '<a href="/system/6">Far&lt;b&gt;</a>' in page and "7.2 ly" in page
    assert page.index("/system/5") < page.index("/system/6")
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
    assert real_client.get("/species").status_code == 404

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
    system = real_client.get(f"/system/{galaxy['a']}").get_data(as_text=True)
    assert "Dominant species:" in system and "Territory of" in system
