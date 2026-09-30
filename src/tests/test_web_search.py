# tests/test_web_search.py

"""
The Flask-served search page (`/search`, `web/searchpage.py`) and the
`search.py` CGI shim that now redirects to it. Same fake data layer as
`test_web_pages.py`; the tests at the bottom use a real database.
"""

import re
from urllib.parse import urlencode

import pytest

import web  # noqa: F401 -- puts src/html/lib on sys.path
import apiclient  # noqa: E402
from stellarObjects import _db  # noqa: E402
from stellarObjects.config import SystemConfig  # noqa: E402
from stellarObjects.spaceSector import SpaceSector  # noqa: E402
from stellarObjects.systemData import StarSystem  # noqa: E402

from tests.test_web_pages import DB, FakeData, app, client, db_client  # noqa: F401 -- fixtures
from tests.webpage_support import run_page  # noqa: E402

_FACETS = {
    "type": [
        {"value": "star", "label": "Stars", "count": 3, "tooltip": None},
        {"value": "planet", "label": "Planets", "count": 2, "tooltip": None},
    ],
    "spectral": [
        {"value": "G", "label": "G", "count": 2, "tooltip": "Yellow <dwarf>"},
        {"value": "K", "label": "K", "count": 1, "tooltip": None},
    ],
    "luminosity": [], "class": [], "body": [], "life": [],
    "moon_class": [], "moon_body": [], "moon_life": [], "density": [],
}


def _result(rows, total=None, offset=0, limit=50):
    return {"rows": rows, "total": len(rows) if total is None else total, "limit": limit,
            "offset": offset, "truncated": False}


class SearchFake(FakeData):
    def __init__(self):
        super().__init__()
        self.search_calls = []
        self.results = None

    def get_search(self, db, texts, tags, sizes=None, limit=None, offsets=None):
        self.search_calls.append({"db": db, "texts": dict(texts), "tags": {k: set(v) for k, v in tags.items()},
                                  "sizes": dict(sizes or {}), "limit": limit, "offsets": dict(offsets or {})})
        results = self.results or {
            "sectors": _result([{"id": 5, "name": "Kepler <Reach>", "edge_mpc": 3.07}]) if texts.get("sector_q") else None,
            "systems": _result([{"id": 7, "name": "Kepler-42", "sector_id": None, "is_binary": 0,
                                 "star_summary": "M5V"}]) if texts.get("system_q") else None,
            "stars": _result([{"name": "Kepler-42 A", "role": "primary", "star_type": "M5V", "radius_km": 118000.4,
                               "star_system_id": 7, "system_name": "Kepler-42", "sector_id": 5}])
            if texts.get("star_q") or tags.get("spectral") else None,
            "planets": None, "moons": None, "belts": None,
        }
        return {
            "facets": _FACETS,
            "autocomplete": {"sectors": ["Kepler <Reach>"], "systems": ["Kepler-42"], "stars": [], "planets": [],
                             "moons": []},
            "facet_labels": {"type:star": "Stars", "spectral:G": "G", "spectral:K": "K"},
            "results": results,
        }


@pytest.fixture
def fake(monkeypatch):
    data = SearchFake()
    for name in ("get_sectors", "get_systems", "auth_me", "get_search"):
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    return data


def _panel(html, panel):
    part = html[html.index(f'id="search-{panel}"'):]
    return part[:part.index("</section>")]


# --- Rendering ---------------------------------------------------------------------

def test_empty_search_shows_form_tags_and_hint(client, fake):
    resp = client.get("/search")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert "<title>Search - " in html
    assert '<span aria-current="page">Search</span>' in html
    assert '<form method="get" action="/search" class="search-form"' in html
    assert "Enter a name above, or pick a tag below" in html
    assert 'id="search-sectors"' not in html
    assert fake.search_calls[0]["db"] == DB
    assert fake.search_calls[0]["texts"] == {k: "" for k in ("sector_q", "system_q", "star_q", "planet_q", "moon_q")}
    # Tags are plain GET links; the tooltip and names are escaped.
    assert '<a class="tag" href="/search?spectral=G" title="Yellow &lt;dwarf&gt;">' in html
    assert '<option value="Kepler &lt;Reach&gt;">' in html
    # No POST forms, no db in any search URL.
    assert 'method="post"' not in html
    assert "/search?db=" not in html and "&amp;db=" not in html


def test_q_searches_every_name(client, fake):
    html = client.get("/search?q=Kepler").get_data(as_text=True)
    call = fake.search_calls[-1]
    assert call["texts"] == {k: "Kepler" for k in ("sector_q", "system_q", "star_q", "planet_q", "moon_q")}
    assert call["offsets"] == {p: 0 for p in ("sectors", "systems", "stars", "planets", "moons", "belts")}
    assert call["limit"] == 50
    assert 'name="q" value="Kepler"' in html
    assert "Kepler &lt;Reach&gt;" in _panel(html, "sectors")
    # Not-yet-moved pages are plain GET links.
    assert f'href="/sector.py?db={DB}&amp;id=5"' in _panel(html, "sectors")
    stars = _panel(html, "stars")
    assert f'href="/system.py?db={DB}&amp;id=7"' in stars
    assert "118,000 km" in stars
    assert "Standalone" in _panel(html, "systems")
    # Chip to remove the name search.
    assert '<span class="filter-chip">Name: “Kepler”<a href="/search"' in html


def test_object_field_overrides_q_and_type_tag_narrows_it(client, fake):
    client.get("/search?q=Kep&star_q=Kepler-42+A&type=planet")
    texts = fake.search_calls[-1]["texts"]
    # With a type tag, q skips sectors/systems; star_q keeps its own term.
    assert texts == {"sector_q": "", "system_q": "", "star_q": "Kepler-42 A", "planet_q": "Kep", "moon_q": "Kep"}
    assert fake.search_calls[-1]["tags"]["type"] == {"planet"}


def test_repeated_tag_params_and_toggle_links(client, fake):
    html = client.get("/search?spectral=G&spectral=K&type=bogus").get_data(as_text=True)
    call = fake.search_calls[-1]
    assert call["tags"]["spectral"] == {"G", "K"}
    assert call["tags"]["type"] == set()  # unknown values dropped
    # Active tags link to the search without them; others add themselves.
    assert re.search(r'<a class="tag active" href="/search\?spectral=K"[^>]*>'
                     r'<span class="sr-only">Selected: </span>G ', html)
    assert '<a class="tag" href="/search?type=star&amp;spectral=G&amp;spectral=K">' in html
    # The form carries active tags as hidden fields.
    assert '<input type="hidden" name="spectral" value="G">' in html
    assert '<input type="hidden" name="spectral" value="K">' in html
    assert '<a class="clear-filters" href="/search">Clear all</a>' in html


def test_sizes_are_parsed_and_bad_numbers_ignored(client, fake):
    html = client.get("/search?planet_min_radius_km=1000&planet_max_radius_km=7000&star_max_radius_km=abc"
                      ).get_data(as_text=True)
    sizes = fake.search_calls[-1]["sizes"]
    assert sizes == {"star": None, "planet": (1000.0, 7000.0), "moon": None}
    assert "Planet size: 1,000–7,000 km" in html
    # The per-object fields are unfolded when one is in use.
    assert '<details class="search-more" open>' in html


def test_empty_form_fields_redirect_to_short_url(client, fake):
    resp = client.get("/search?q=Kepler&sector_q=&system_q=&star_min_radius_km=&spectral=G&db=other&stars_page=2")
    assert resp.status_code == 302
    assert resp.headers["Location"] == "/search?q=Kepler&spectral=G&stars_page=2"
    assert not fake.search_calls


def test_each_panel_pages_on_its_own(client, fake):
    fake.results = {
        "sectors": _result([{"id": i, "name": f"S{i}", "edge_mpc": 1.0} for i in range(50)], total=120, offset=50),
        "systems": None,
        "stars": _result([], total=70, offset=50),
        "planets": None, "moons": None, "belts": None,
    }
    html = client.get("/search?q=S&sectors_page=2&stars_page=9").get_data(as_text=True)
    assert fake.search_calls[-1]["offsets"]["sectors"] == 50
    assert fake.search_calls[-1]["offsets"]["stars"] == 400
    sectors = _panel(html, "sectors")
    assert "Sectors <span class=\"count\">(120)</span>" in sectors
    # The pager keeps the search and the other panel's page (as returned).
    assert 'href="/search?q=S&amp;stars_page=2&amp;sectors_page=3#search-sectors"' in sectors
    assert 'aria-label="Sectors result pages"' in sectors
    assert 'href="/search?q=S&amp;sectors_page=2&amp;stars_page=1#search-stars"' in _panel(html, "stars")


def test_no_results_message(client, fake):
    fake.results = {p: None for p in ("sectors", "systems", "stars", "planets", "moons", "belts")}
    html = client.get("/search?type=planet").get_data(as_text=True)
    assert "No matching objects." in html


def test_header_search_box_submits_to_search(client, fake):
    html = client.get("/").get_data(as_text=True)
    assert '<form class="site-search" role="search" method="get" action="/search">' in html
    assert 'name="q" placeholder="Search names"' in html


# --- CGI shim ---------------------------------------------------------------------------

def test_cgi_shim_redirects_get_with_every_filter():
    result = run_page("http://127.0.0.1:9/api", "search.py", query=urlencode([
        ("db", "x"), ("system_q", "Kepler 42"), ("spectral", "G"), ("spectral", "K"), ("sector_q", ""),
        ("planet_min_radius_km", "100"), ("stars_page", "3"), ("junk", "1"),
    ]))
    assert result.status_code == 301
    assert result.headers["Location"] == \
        "/search?system_q=Kepler+42&planet_min_radius_km=100&spectral=G&spectral=K&stars_page=3"


def test_cgi_shim_redirects_post():
    result = run_page("http://127.0.0.1:9/api", "search.py", method="POST",
                      body={"db": "x", "type": "belt", "belts_page": "2"})
    assert result.status_code == 301
    assert result.headers["Location"] == "/search?type=belt&belts_page=2"


def test_cgi_shim_without_params():
    result = run_page("http://127.0.0.1:9/api", "search.py", query={"db": "x"})
    assert result.status_code == 301
    assert result.headers["Location"] == "/search"


# --- Real database, in-process ---------------------------------------------------------

def test_real_search_pages_a_result_panel(db_client, mysql_config, monkeypatch):
    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")
    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    for i in range(60):
        _db.save_sector(SpaceSector(f"Pager Sector {i:03d}", edge_ly=10.0), config=mysql_config)
    _db.save_sector(SpaceSector("Unrelated", edge_ly=10.0), config=mysql_config)

    page1 = db_client.get("/search?sector_q=Pager")
    assert page1.status_code == 200
    panel1 = _panel(page1.get_data(as_text=True), "sectors")
    assert "(60)" in panel1
    assert "Pager Sector 049" in panel1 and "Pager Sector 050" not in panel1
    assert "Showing 1&ndash;50 of 60" in panel1
    assert 'href="/search?sector_q=Pager&amp;sectors_page=2#search-sectors"' in panel1

    panel2 = _panel(db_client.get("/search?sector_q=Pager&sectors_page=2").get_data(as_text=True), "sectors")
    assert "Pager Sector 050" in panel2 and "Pager Sector 059" in panel2
    assert "Pager Sector 049" not in panel2 and "Unrelated" not in panel2
    assert "Showing 51&ndash;60 of 60" in panel2


def test_real_q_and_tags(db_client, mysql_config):
    sector = SpaceSector("Tagged Sector", edge_ly=10.0)
    for star in ("G2V", "M5V"):
        cfg = SystemConfig()
        cfg.STAR_TYPE = star
        cfg.PLANETS = False
        cfg.BINARY_SYSTEM = False
        sector.add_system(StarSystem(system_config=cfg), position=(1.0, 1.0, 1.0), system_config=cfg)
    _db.save_sector(sector, config=mysql_config)

    html = db_client.get("/search?q=Tagged").get_data(as_text=True)
    assert "Tagged Sector" in _panel(html, "sectors")

    html = db_client.get("/search?spectral=G").get_data(as_text=True)
    stars = _panel(html, "stars")
    assert "(1)" in stars and "G2V" in stars and "M5V" not in stars
    assert '<a class="tag active" href="/search"' in html
