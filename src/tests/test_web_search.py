# tests/test_web_search.py

"""
The Flask-served search page (`/search`, `web/searchpage.py`) and the
old `search.py` URL that now redirects to it. Same fake data layer as
`test_web_pages.py`; the tests at the bottom use a real database.
"""

import re

import pytest

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.db import store  # noqa: E402
from planetgen.generation.config import SystemConfig  # noqa: E402
from planetgen.galaxy.sector import SpaceSector  # noqa: E402
from planetgen.generation.system import StarSystem  # noqa: E402

from tests.test_web_pages import DB, FakeData, app, client, db_client  # noqa: F401 -- fixtures

_FACETS = {
    "type": [
        {"value": "star", "label": "Stars", "count": 3, "tooltip": None},
        {"value": "planet", "label": "Planets", "count": 2, "tooltip": None},
    ],
    "spectral": [
        {"value": "G", "label": "G", "count": 2, "tooltip": "Yellow <dwarf>"},
        {"value": "K", "label": "K", "count": 1, "tooltip": None},
    ],
    "luminosity": [], "class": [], "body": [], "life": [], "equipment": [],
    "moon_class": [], "moon_body": [], "moon_life": [], "moon_equipment": [], "density": [],
}


def _result(rows, total=None, offset=0, limit=50):
    return {"rows": rows, "total": len(rows) if total is None else total, "limit": limit,
            "offset": offset, "truncated": False}


class SearchFake(FakeData):
    def __init__(self):
        super().__init__()
        self.search_calls = []
        self.results = None

    def get_search(self, db, texts, tags, sizes=None, limit=None, offsets=None, panels=None, sorts=None, stars=None, ranges=None):
        self.search_calls.append({"stars": stars, "ranges": dict(ranges or {}), "db": db, "texts": dict(texts), "tags": {k: set(v) for k, v in tags.items()},
                                  "sizes": dict(sizes or {}), "limit": limit, "offsets": dict(offsets or {}),
                                  "panels": panels, "sorts": dict(sorts or {})})
        results = self.results or {
            "sectors": _result([{"id": 5, "name": "Kepler <Reach>", "edge_mpc": 3.07}]) if texts.get("sector_q") else None,
            "systems": _result([{"id": 7, "name": "Kepler-42", "sector_id": None, "is_binary": 0,
                                 "star_summary": "M5V", "star_count": 1, "object_count": 6, "has_habitable_world": True}]) if texts.get("system_q") or stars else None,
            "stars": _result([{"id": 3, "ref": "star:FE81000A2B-0000005-001", "name": "Kepler-42 A", "role": "primary", "star_type": "M5V", "radius_km": 118000.4,
                               "star_system_id": 7, "system_name": "Kepler-42", "sector_id": 5}])
            if texts.get("star_q") or tags.get("spectral") else None,
            "planets": None, "moons": None, "belts": None,
        }
        if panels is not None:
            results = {panel: result if panel in panels else None for panel, result in results.items()}
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
    # UX.211: a tag's count says what it counts.
    counts = re.findall(r'<span class="tag-count">([^<]+)</span>', html)
    assert counts and all(re.fullmatch(r"\d+ results?", text) for text in counts), counts
    # No POST forms, no db in any search URL.
    assert 'method="post"' not in html
    assert "/search?db=" not in html and "&amp;db=" not in html


def test_q_searches_every_name(client, fake):
    html = client.get("/search?q=Kepler").get_data(as_text=True)
    call = fake.search_calls[0]
    assert call["texts"] == {k: "Kepler" for k in ("sector_q", "system_q", "star_q", "planet_q", "moon_q")}
    assert call["offsets"] == {p: 0 for p in ("sectors", "systems", "stars", "planets", "moons", "belts", "phenomena")}
    assert call["limit"] == 50
    # UX.97: the first request runs no panel; then every panel is asked for on its own.
    assert call["panels"] == []
    assert [c["panels"] for c in fake.search_calls[1:]] == [[p] for p in
                                                           ("sectors", "systems", "stars", "planets", "moons", "belts", "phenomena")]
    assert 'name="q" value="Kepler"' in html
    assert "Kepler &lt;Reach&gt;" in _panel(html, "sectors")
    # Not-yet-moved pages are plain GET links.
    assert 'href="/sector/5"' in _panel(html, "sectors")
    stars = _panel(html, "stars")
    assert 'href="/system/7"' in stars
    assert "118,000 km" in stars
    assert "Standalone" in _panel(html, "systems")
    # "Show on Galaxy Map" (MAP.25): a sector always, a standalone system never.
    assert 'href="/sector/5/galaxy" aria-label="Show Kepler &lt;Reach&gt; on the Galaxy Map">Show</a>' in _panel(html, "sectors")
    assert "/galaxy" not in _panel(html, "systems")
    # The name is not echoed as a chip, and the page has no second search box (UX.54).
    assert "Name: “Kepler”" not in html and 'class="site-search"' not in html


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
    # Tag groups fold away; only a group with a selected tag starts open.
    groups = re.findall(r'<details class="tag-group"( open)?>\s*<summary><h3>([^<]+)</h3>'
                        r'<span class="tag-group-meta">([^<]+)</span>', html)
    assert groups and all(bool(is_open) == ("selected" in meta) for is_open, _t, meta in groups)
    assert any(is_open and "2 selected" in meta for is_open, _t, meta in groups)


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
    assert fake.search_calls[0]["offsets"]["sectors"] == 50
    assert fake.search_calls[0]["offsets"]["stars"] == 400
    # The past-the-end page is asked for again on its own; the API answers with the last real page.
    stars_calls = [c for c in fake.search_calls if c["panels"] == ["stars"]]
    assert stars_calls[-1]["offsets"] == {"stars": 400}
    sectors = _panel(html, "sectors")
    assert "Sectors <span class=\"count\">(120)</span>" in sectors
    # The pager keeps the search and the other panel's page (as asked; it shows the last real one).
    assert 'href="/search?q=S&amp;stars_page=9&amp;sectors_page=3#search-sectors"' in sectors
    assert 'aria-label="Sectors pages"' in sectors
    assert 'href="/search?q=S&amp;sectors_page=2&amp;stars_page=1#search-stars"' in _panel(html, "stars")


def test_a_panel_that_times_out_says_so_and_the_others_show(client, fake, monkeypatch):
    """UX.97: each panel runs under its own limit; one that hits it is a note, not the busy page."""
    answer = fake.get_search

    def one_slow(db, texts, tags, sizes=None, limit=None, offsets=None, panels=None, sorts=None, stars=None, ranges=None):
        if panels == ["systems"]:
            raise apiclient.ApiError("QUERY_TIMEOUT", status_code=504)
        return answer(db, texts, tags, sizes=sizes, limit=limit, offsets=offsets, panels=panels, sorts=sorts, stars=stars)

    monkeypatch.setattr(apiclient, "get_search", one_slow)
    resp = client.get("/search?q=Kepler")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert "This took too long; try again or narrow the search." in _panel(html, "systems")
    assert 'href="/search?q=Kepler"' in _panel(html, "systems")
    assert "Kepler &lt;Reach&gt;" in _panel(html, "sectors")

    def broken(db, texts, tags, sizes=None, limit=None, offsets=None, panels=None, sorts=None, stars=None, ranges=None):
        if panels == ["systems"]:
            raise apiclient.ApiError("down", status_code=500)
        return answer(db, texts, tags, sizes=sizes, limit=limit, offsets=offsets, panels=panels, sorts=sorts, stars=stars)

    monkeypatch.setattr(apiclient, "get_search", broken)
    assert client.get("/search?q=Kepler").status_code == 502


def test_a_panels_table_route_pages_one_panel_alone(client, fake):
    fake.results = {
        "sectors": _result([{"id": i, "name": f"S{i}", "edge_mpc": 1.0} for i in range(50)], total=120, offset=50),
        "systems": None, "stars": None, "planets": None, "moons": None, "belts": None,
    }
    data = client.get("/table/search-sectors?q=S&offset=50&limit=50").get_json()
    assert data["total"] == 120 and len(data["rows"]) == 50
    assert data["rows"][0][0] == {"text": "S0", "href": "/sector/0"}
    call = fake.search_calls[-1]
    assert call["panels"] == ["sectors"] and call["offsets"] == {"sectors": 50} and call["texts"]["sector_q"] == "S"


def test_no_results_message(client, fake):
    fake.results = {p: None for p in ("sectors", "systems", "stars", "planets", "moons", "belts", "phenomena")}
    html = client.get("/search?type=planet").get_data(as_text=True)
    assert "No sectors, systems, stars, planets or moons match." in html


def test_the_header_search_icon_links_to_search(client, fake):
    """UX.160: the header's search icon is a link to this page, and this page leaves it out."""
    assert 'class="site-search-link" name="search" label="Search" title="Search" href="/search"' in client.get("/").get_data(as_text=True)
    assert "site-search-link" not in client.get("/search").get_data(as_text=True)


# --- Real database, in-process ---------------------------------------------------------

def test_real_search_pages_a_result_panel(db_client, mysql_config, monkeypatch):
    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")
    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    for i in range(60):
        store.save_sector(SpaceSector(f"Pager Sector {i:03d}", edge_ly=10.0), config=mysql_config)
    store.save_sector(SpaceSector("Unrelated", edge_ly=10.0), config=mysql_config)

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
    store.save_sector(sector, config=mysql_config)

    html = db_client.get("/search?q=Tagged").get_data(as_text=True)
    assert "Tagged Sector" in _panel(html, "sectors")

    html = db_client.get("/search?spectral=G").get_data(as_text=True)
    stars = _panel(html, "stars")
    assert "(1)" in stars and "G2V" in stars and "M5V" not in stars
    assert '<a class="tag active" href="/search"' in html


def test_real_phenomenon_and_class_tags(db_client, mysql_config):
    from planetgen.generation.phenomena.nebula import Nebula
    from planetgen.generation.phenomena.compact_remnant import BlackHole

    for name, nebula_class in (("Crab Mist", "D"), ("Faint Veil", "A")):
        nebula = Nebula(SystemConfig(), name=name, nebula_class=nebula_class)
        store.save_phenomenon(nebula, SystemConfig(), "nebula", config=mysql_config)
    store.save_phenomenon(BlackHole(SystemConfig()), SystemConfig(), "black-hole", config=mysql_config)

    html = db_client.get("/search").get_data(as_text=True)
    assert "Phenomenon Class" in html and "D: Classical H II region" in html

    html = db_client.get("/search?phenomenon=nebula").get_data(as_text=True)
    panel = _panel(html, "phenomena")
    assert "(2)" in panel and "Crab Mist" in panel and "Faint Veil" in panel
    assert "/phenomenon/nebula/" in panel

    panel = _panel(db_client.get("/search?phenomenon_class=nebula:D").get_data(as_text=True), "phenomena")
    assert "(1)" in panel and "Crab Mist" in panel and "Faint Veil" not in panel

    panel = _panel(db_client.get("/search?phenomenon=black_hole").get_data(as_text=True), "phenomena")
    assert "(1)" in panel and "Black Hole" in panel

    # UX.44: a class tag only narrows its own type, so the black hole
    # stays beside the class D nebula.
    html = db_client.get("/search?phenomenon=black_hole&phenomenon=nebula&phenomenon_class=nebula:D")
    panel = _panel(html.get_data(as_text=True), "phenomena")
    assert "(2)" in panel and "Crab Mist" in panel and "Black Hole" in panel and "Faint Veil" not in panel

    panel = _panel(db_client.get("/search?phenomenon_class=nebula:D&phenomenon_class=nebula:A")
                   .get_data(as_text=True), "phenomena")
    assert "(2)" in panel and "Crab Mist" in panel and "Faint Veil" in panel


def test_tags_in_one_group_combine_with_or(db_client, mysql_config):
    # UX.44: a star has one spectral type, so two spectral tags mean
    # either; a second group (luminosity) still narrows with AND.
    sector = SpaceSector("Either Sector", edge_ly=10.0)
    for star in ("G2V", "M5V", "K1III"):
        cfg = SystemConfig()
        cfg.STAR_TYPE = star
        cfg.PLANETS = False
        cfg.BINARY_SYSTEM = False
        sector.add_system(StarSystem(system_config=cfg), position=(1.0, 1.0, 1.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config)

    stars = _panel(db_client.get("/search?spectral=G&spectral=M").get_data(as_text=True), "stars")
    assert "(2)" in stars and "G2V" in stars and "M5V" in stars and "K1III" not in stars

    stars = _panel(db_client.get("/search?spectral=G&spectral=K&luminosity=V").get_data(as_text=True), "stars")
    assert "(1)" in stars and "G2V" in stars and "K1III" not in stars


def test_groups_with_no_matches_are_left_out(client, fake):
    """UX.54: only groups that have results get a card; with none, one line."""
    empty = {"rows": [], "total": 0, "total_capped": False}
    fake.results = {p: dict(empty) for p in ("sectors", "systems", "stars", "planets", "moons", "belts", "phenomena")}
    html = client.get("/search?q=zzz").get_data(as_text=True)
    assert 'id="search-sectors"' not in html and 'id="search-stars"' not in html
    assert "No sectors, systems, stars, planets or moons match." in html


def test_a_star_result_links_to_its_own_page(client, fake):
    """NAV.8: a body result opens the body's own page (UX.173 and UX.176 retired the row on the system page)."""
    html = client.get("/search?star_q=Kepler").get_data(as_text=True)
    assert 'href="/object/star:FE81000A2B-0000005-001"' in html and "#star-3" not in html


def _body_row(id_, name, **extra):
    return {"id": id_, "ref": f"planet:FE81000A2B-0000005-{id_:03d}", "name": name, "planet_class": "M", "body_type": "t", "life_class_id": None, "radius_km": 6371.0,
            "star_system_id": 7, "system_name": "Kepler-42", "sector_id": 5, "equipment": None, "phi4": None,
            **{f"tier_{d}": None for d in ("pressure", "temperature", "air", "water", "radiation", "suitability")}, **extra}


def test_planet_results_show_habitability_symbols_not_notes(client, fake):
    """UX.100: a scored planet's row carries the overall symbol and the six domains; an unscored one a dash."""
    scored = _body_row(11, "Scored", equipment="Required: gas-mix mask.", phi4=0.8, tier_pressure=0, tier_temperature=1, tier_air=2,
                       tier_water=0, tier_radiation=3, tier_suitability=1)
    fake.results = {"planets": _result([scored, _body_row(12, "Unscored", body_type="g")]),
                    **{p: None for p in ("sectors", "systems", "stars", "moons", "belts", "phenomena")}}
    html = client.get("/search?planet_q=Sc").get_data(as_text=True)
    panel = _panel(html, "planets")
    assert "data-col=\"phi4\"" in panel and ">Habitability" in panel
    for label in ("Overall habitability: Green", "Pressure: Blue", "Air chemistry: Yellow", "Radiation: Red"):
        assert f'aria-label="{label}"' in panel, label
    assert panel.count("hab-symbol ") == 7
    assert "<em>—</em>" in panel[panel.index("Unscored"):]


_SCORED = {"tier_pressure": 0, "tier_temperature": 1, "tier_air": 2, "tier_water": 0, "tier_radiation": 1, "tier_suitability": 1,
           "life_microbial": 2, "life_flora": 0, "life_fauna": 0, "life_sentient": 0}


def test_score_filters_show_as_chips_and_pass_to_the_api(client, fake):
    """UX.214: a score filter in the URL is a removable chip named in words and goes to the API as its facet."""
    fake.results = {"planets": _result([_body_row(11, "Scored", phi4=0.8, **_SCORED)]),
                    **{p: None for p in ("sectors", "systems", "stars", "moons", "belts", "phenomena")}}
    html = client.get("/search?habitability=green&domain=suitability:blue&life_level=flora:5&living=civilization"
                      "&moon_living=alive").get_data(as_text=True)
    sent = fake.search_calls[-1]["tags"]
    assert sent["habitability"] == {"green"} and sent["domain"] == {"suitability:blue"}
    assert sent["life_level"] == {"flora:5"} and sent["living"] == {"civilization"} and sent["moon_living"] == {"alive"}
    assert html.count('<span class="filter-chip">') == 5      # the labels come from the API's facet_labels
    assert "2 of 9" in _panel(html, "planets") and "Green" in _panel(html, "planets")


def test_a_score_column_sorts_and_the_sort_goes_to_the_api(client, fake):
    """UX.214: the Habitability, Suitability and life columns of the planet table sort; best first is the descending order."""
    fake.results = {"planets": _result([_body_row(11, "Scored", phi4=0.8, **_SCORED)]),
                    **{p: None for p in ("sectors", "systems", "stars", "moons", "belts", "phenomena")}}
    html = client.get("/search?planet_q=Sc").get_data(as_text=True)
    panel = _panel(html, "planets")
    for key in ("phi4", "suitability", "microbial", "flora", "fauna", "sentient"):
        assert f'data-col="{key}"' in panel and "aria-sort" in panel
    client.get("/search?planet_q=Sc&planets_sort=flora")
    assert fake.search_calls[-1]["sorts"] == {"planet": "flora"}
    client.get("/search?planet_q=Sc&planets_sort=phi4&planets_order=desc")
    assert fake.search_calls[-1]["sorts"] == {"planet": "-phi4"}
    client.get("/search?planet_q=Sc&planets_sort=nonsense")
    assert fake.search_calls[-1]["sorts"] == {}


# --- Ranges on physical values (UX.227) ------------------------------------------------------

def _planet_only(fake, rows):
    fake.results = {"planets": _result(rows), **{p: None for p in ("sectors", "systems", "stars", "moons", "belts", "phenomena")}}


def test_value_ranges_show_as_chips_and_go_to_the_api_in_their_unit(client, fake):
    """UX.227: a range on a physical value is a removable chip with its unit, and the API gets the bounds and the unit named."""
    _planet_only(fake, [_body_row(11, "Scored", phi4=0.8, **_SCORED)])
    html = client.get("/search?planet_min_gravity=0.8&planet_max_gravity=1.2&planet_max_surface_dose=2"
                      "&moon_min_surface_pressure=50&moon_surface_pressure_unit=kPa&star_min_luminosity=0.5").get_data(as_text=True)
    assert fake.search_calls[-1]["ranges"] == {
        "planet": {"gravity": (0.8, 1.2, None), "surface_dose": (None, 2.0, None)},
        "moon": {"surface_pressure": (50.0, None, "kPa")}, "star": {"luminosity": (0.5, None, None)}}
    assert 'Planet surface gravity: 0.8–1.2 g₀' in html and 'Planet surface radiation dose: ≤ 2 mSv/yr' in html
    assert 'Moon surface pressure: ≥ 50 kPa' in html and 'Star luminosity: ≥ 0.5 L☉' in html
    # a chip's cross takes that one range away and keeps the rest
    cross = re.search(r'Planet surface gravity: 0\.8–1\.2 g₀<a href="([^"]+)"', html).group(1)
    assert "planet_min_gravity" not in cross and "planet_max_gravity" not in cross and "planet_max_surface_dose=2" in cross
    assert '<details class="search-more search-values" open>' in html


def test_the_form_has_a_row_with_two_boxes_and_a_unit_menu_for_each_value(client, fake):
    html = client.get("/search?moon_min_surface_pressure=50&moon_surface_pressure_unit=kPa").get_data(as_text=True)
    assert "<summary>Search by physical values</summary>" in html
    assert 'name="planet_min_gravity"' in html and 'name="planet_max_gravity"' in html and 'name="star_max_age"' in html
    # the unit menu offers the units of the kind and shows the chosen one
    menu = re.search(r'<select name="moon_surface_pressure_unit"[^>]*>(.*?)</select>', html, re.S).group(1)
    assert '<option value="kPa" selected>kPa</option>' in menu and '<option value="bar">bar</option>' in menu
    assert '<select name="planet_surface_pressure_unit"' in html and 'value="bar" selected' in html
    # gases: the short list has rows, any other gas has the "another gas" row, and a gas in use gets its own row
    assert 'name="planet_min_gas_O2"' in html and 'name="planet_min_gas_C6H6"' not in html and 'name="planet_other_gas"' in html
    assert '<option value="C6H6">C6H6</option>' in html
    html = client.get("/search?planet_min_gas_C6H6=0.001").get_data(as_text=True)
    assert 'name="planet_min_gas_C6H6"' in html and '<option value="C6H6">C6H6</option>' not in html.split('name="planet_other_gas"')[1].split("</select>")[0]
    # a number of moons has no unit to choose and no menu
    assert 'name="planet_moon_count_unit"' not in html and 'name="planet_min_moon_count"' in html


def test_a_bad_value_range_is_ignored_like_a_bad_size(client, fake):
    _planet_only(fake, [_body_row(11, "Scored", phi4=0.8, **_SCORED)])
    client.get("/search?planet_min_gravity=abc&planet_max_surface_dose=nan&planet_min_surface_temperature=300"
               "&planet_max_surface_temperature=200&planet_min_surface_pressure=1&planet_surface_pressure_unit=km&moon_max_gravity=inf")
    assert fake.search_calls[-1]["ranges"] == {"planet": {"surface_pressure": (1.0, None, None)}}


def test_a_unit_that_says_nothing_is_dropped_from_the_url(client, fake):
    """The form sends every row's unit: one with no bound beside it, or the key's own unit, redirects to the short form."""
    resp = client.get("/search?planet_min_gravity=1&planet_gravity_unit=g0&planet_surface_pressure_unit=kPa&planet_other_gas=&planet_other_gas_unit=bar")
    assert resp.status_code == 302 and resp.headers["Location"] == "/search?planet_min_gravity=1"
    resp = client.get("/search?planet_min_gravity=1&planet_gravity_unit=m_s2")
    assert resp.status_code == 200


def test_the_another_gas_row_turns_into_the_key_of_that_gas(client, fake):
    resp = client.get("/search?moon_other_gas=C6H6&moon_other_gas_min=0.001&moon_other_gas_unit=mbar&moon_other_gas_max=")
    assert resp.status_code == 302 and resp.headers["Location"] == "/search?moon_min_gas_C6H6=0.001&moon_gas_C6H6_unit=mbar"
    client.get(resp.headers["Location"])
    assert fake.search_calls[-1]["ranges"] == {"moon": {"gas_C6H6": (0.001, None, "mbar")}}
    resp = client.get("/search?planet_other_gas=Unobtainium&planet_other_gas_min=1")
    assert resp.status_code == 302 and resp.headers["Location"] == "/search"


def test_an_order_by_a_value_goes_to_the_api_unless_a_table_header_orders_the_panel(client, fake):
    _planet_only(fake, [_body_row(11, "Scored", phi4=0.8, **_SCORED)])
    html = client.get("/search?planet_q=Sc&planet_value_sort=-gravity&star_value_sort=age&moon_value_sort=luminosity").get_data(as_text=True)
    assert fake.search_calls[-1]["sorts"] == {"planet": "-gravity", "star": "age"}       # a moon has no luminosity: ignored
    assert "Planets ordered by surface gravity, high to low" in html and "Stars ordered by age, low to high" in html
    assert '<option value="-gravity" selected>' in html
    client.get("/search?planet_q=Sc&planet_value_sort=-gravity&planets_sort=flora")
    assert fake.search_calls[-1]["sorts"] == {"planet": "flora"}
    # an order alone is not a search
    html = client.get("/search?planet_value_sort=gravity").get_data(as_text=True)
    assert "Enter a name above, or pick a tag below" in html


def test_the_values_a_search_used_are_shown_with_their_units(client, fake):
    row = _body_row(11, "Scored", phi4=0.8, **_SCORED)
    row["values"] = {"gravity": {"value": 0.93417, "unit": "g0", "symbol": "g₀"},
                     "surface_pressure": {"value": 101.3, "unit": "kPa", "symbol": "kPa"}}
    _planet_only(fake, [row])
    html = client.get("/search?planet_min_gravity=0.5").get_data(as_text=True)
    assert "Surface gravity 0.9342 g₀; Surface pressure 101.3 kPa" in _panel(html, "planets")
    _planet_only(fake, [_body_row(12, "Plain", phi4=0.8, **_SCORED)])
    assert "Surface gravity" not in _panel(client.get("/search?planet_q=Pl").get_data(as_text=True), "planets")


# --- The stars filter (UX.239) ---------------------------------------------------------------

def test_the_stars_filter_is_a_menu_in_the_form_and_goes_to_the_api(client, fake):
    html = client.get("/search?stars=3").get_data(as_text=True)
    assert fake.search_calls[0]["stars"] == "3"
    assert '<select name="stars" id="search-star-count">' in html and '<option value="3" selected>Trinary</option>' in html
    for label in ("Single", "Binary", "Trinary", "Quaternary", "Quintenary and beyond"):
        assert f">{label}</option>" in html
    assert "Stars: Trinary" in html                                    # the removable chip
    assert 'id="search-systems"' in html                               # the Systems panel runs with no name
    assert 'href="/search"' in html


def test_stars_typed_in_the_search_box_is_the_same_filter_and_the_rest_stays_a_name(client, fake):
    client.get("/search?q=Kepler+stars:3%2B")
    call = fake.search_calls[0]
    assert call["stars"] == "3+" and call["texts"]["system_q"] == "Kepler"
    client.get("/search?q=stars:Binary")
    call = fake.search_calls[-1]
    assert call["stars"] == "binary" and call["texts"]["system_q"] == ""


def test_a_stars_filter_that_is_not_one_is_ignored_on_the_page(client, fake):
    client.get("/search?stars=lots&q=Kepler")
    assert fake.search_calls[0]["stars"] is None


def test_no_filter_sends_none_and_the_menu_says_any(client, fake):
    html = client.get("/search?q=Kepler").get_data(as_text=True)
    assert fake.search_calls[0]["stars"] is None and '<option value="">Any number</option>' in html


def test_real_search_finds_a_known_binary_and_a_known_triple(db_client, mysql_config):
    from planetgen.util import draw
    for name, count in (("Lone Omicron", 1), ("Pair Omicron", 2), ("Triple Omicron", 3)):
        config = SystemConfig()
        config.PLANETS = False
        config.STARS = count
        with draw.bound(f"stars-filter-page:{name}"):
            system = StarSystem(system_config=config)
        assert len(system.stars) == count
        system.name = name
        store.save_system(system, config, config=mysql_config)

    binary = _panel(db_client.get("/search?stars=binary").get_data(as_text=True), "systems")
    assert "Pair Omicron" in binary and "Triple Omicron" not in binary and "Lone Omicron" not in binary
    triple = _panel(db_client.get("/search?q=stars:3").get_data(as_text=True), "systems")
    assert "Triple Omicron" in triple and "Pair Omicron" not in triple
    more = _panel(db_client.get("/search?stars=2%2B&q=Omicron").get_data(as_text=True), "systems")
    assert "Pair Omicron" in more and "Triple Omicron" in more and "Lone Omicron" not in more
