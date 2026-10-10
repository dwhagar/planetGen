# tests/test_web_galaxy.py

"""
The Flask-served Galaxy Map (`/galaxy`) and its tile JSON
(`/galaxy/tiles`), `planetgen/web/galaxy_views.py`, plus the old `galaxy.py`/
`galaxy_tiles.py` URLs that redirect there.

Most tests fake the data layer (the `apiclient` functions the view calls,
and the two `tilecache` uses for tiles), the same way
`test_web_pages.py` does. The tests at the bottom run against a real
throwaway database in-process and are skipped without a MySQL server.
"""

import json
import os
import re

import pytest
from flask import url_for

from planetgen.web.app import create_app
from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.api.config import Config

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.web.lib import tilecache  # noqa: E402
from planetgen.db import store as _db  # noqa: E402
from planetgen.galaxy.sector import SpaceSector  # noqa: E402
from planetgen.web import csrf  # noqa: E402
from planetgen.web.helpers import page_url  # noqa: E402


DB = "planetgen_galaxy_test"
STAMP = "00000000000000aa"


class _FakeConfig(Config):
    WEB_DATABASE = DB
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


def _placed(i, x=100.0, y=50.0, name=None, radius=None):
    return {
        "id": i, "name": name or f"Placed {i:03d}", "x": x, "y": y, "z": 0.0,
        "galactic_radius_pc": radius if radius is not None else 1000.0 + i,
        "ring_index": 10 + i, "layer_index": 0, "ring_slot_index": 0, "system_count": i % 4,
    }


class FakeData:
    """Stands in for the API: the page's two reads and tilecache's two."""

    def __init__(self):
        self.sectors = [_placed(1), _placed(2, x=-100.0), _placed(3, x=-100.0, y=-50.0)]
        self.shape = {"edge_pc": 3.526, "outer_ring_index": 400}
        self.dbs = set()
        self.tile_calls = []
        self.fail_tiles = None
        self.stage_calls = []
        self.polity_total = 2

    def get_polities(self, db, limit=None, offset=None):
        self.dbs.add(db)
        return {"items": [], "total": self.polity_total, "limit": limit, "offset": offset}

    def get_galaxy_sectors(self, db):
        self.dbs.add(db)
        return self.sectors

    def get_galaxy_shape(self, db):
        self.dbs.add(db)
        return self.shape

    def get_galaxy_changes(self, db, since=None):
        self.dbs.add(db)
        return {"stamp": STAMP, "state": "s", "full": since is None, "tiles": []}

    def get_galaxy_tiles(self, db, tile_keys):
        self.dbs.add(db)
        if self.fail_tiles:
            raise self.fail_tiles
        self.tile_calls.append(list(tile_keys))
        return {
            "tiles": {key: {"placed": [{"id": 1, "name": "T", "x": 1.0, "y": 0.0, "z": 0.0}], "planned": []}
                      for key in tile_keys},
            "edge_pc": 3.526, "has_shape": True,
        }

    def get_galaxy_stage(self, db, at=None):
        self.dbs.add(db)
        self.stage_calls.append(at)
        return {"at": at, "child_m": 243 if at is None else 27, "children": [], "sectors": None}

    def auth_me(self, cookie_header):
        return None


@pytest.fixture
def fake(monkeypatch, tmp_path):
    data = FakeData()
    for name in ("get_galaxy_sectors", "get_galaxy_shape", "get_polities", "auth_me"):
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    monkeypatch.setattr(tilecache, "get_galaxy_changes", data.get_galaxy_changes)
    monkeypatch.setattr(tilecache, "get_galaxy_tiles", data.get_galaxy_tiles)
    monkeypatch.setattr(tilecache, "get_galaxy_stage", data.get_galaxy_stage)
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    monkeypatch.delenv("PLANETGEN_TILE_CACHE_MAX_MB", raising=False)
    data.cache_root = tmp_path / "tiles"
    return data


@pytest.fixture
def app():
    application = create_app(_FakeConfig)
    application.testing = True
    return application


@pytest.fixture
def client(app):
    return app.test_client()


def _scene(html):
    match = re.search(r'<script type="application/json" id="galaxymap3d-data">(.*?)</script>', html, re.S)
    return json.loads(match.group(1))


# --- /galaxy ------------------------------------------------------------------------

def test_galaxy_page_renders_map_and_quadrant_summary(client, fake, app):
    resp = client.get("/galaxy")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert "<title>Galaxy Map - " in html
    assert re.search(r'<script type="module" src="/static/galaxymap3d.js\?v=[^"]+"></script>', html)
    # The Bookmarks menu (MAP.23), for this database's list.
    assert re.search(r'<script type="module" src="/static/bookmarks.js\?v=[^"]+"></script>', html)
    assert f'data-bookmarks-menu data-bookmarks-keys="map" data-bookmark-db="{DB}"' in html
    assert 'id="galaxymap3d-canvas"' in html
    assert "3 placed sectors" in html
    # The Galaxy section is current, with a breadcrumb.
    assert re.search(r'<a href="/galaxy" aria-current="page">Galaxy</a>', html)
    assert '<li><span aria-current="page">Galaxy</span></li>' in html
    # Quadrant rows are plain GET links.
    for label in ("I", "II", "III", "IV"):
        assert f'<a href="/galaxy?quadrant={label}#galaxy-table">Quadrant {label}</a>' in html
    # The database never leaves the server.
    assert fake.dbs == {DB}
    assert not re.search(r'/galaxy[^"\s]*db=', html)


def test_galaxy_scene_data_points_at_new_urls(client, fake, app):
    scene = _scene(client.get("/galaxy").get_data(as_text=True))
    assert scene["fetchPath"] == "/galaxy/tiles"
    assert "db" not in scene
    assert scene["storageKey"] == DB  # the browser's tile cache keeps its old key
    assert "{id}" in scene["sectorUrl"]
    with app.test_request_context("/"):
        assert scene["sectorUrl"].replace("{id}", "5") == page_url("sector", sector_id=5)
        assert scene["phenomenonUrl"].replace("{type}", "supernova_remnant").replace("{id}", "7") == page_url(
            "phenomenon", phenomenon_type="supernova_remnant", phenomenon_id=7)
        assert scene["systemUrl"].replace("{id}", "9") == page_url("system", system_id=9)
        assert scene["stagePath"] == url_for("web.galaxy_stage")
    assert scene["initial"]["stamp"] == STAMP
    assert scene["initial"]["tiles"]
    # Visitors get no Generate buttons.
    assert scene["generate"] is None


def test_galaxy_scene_data_gives_an_admin_the_generate_target(client, fake, monkeypatch):
    monkeypatch.setattr(apiclient, "auth_me", lambda cookie_header: {"username": "admin"})
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")
    target = _scene(client.get("/galaxy").get_data(as_text=True))["generate"]
    assert target["url"] == "/admin/generate"
    assert target["csrfField"] == csrf.FIELD_NAME
    assert target["csrfToken"]


def test_galaxy_first_view_is_written_to_the_disk_cache(client, fake):
    client.get("/galaxy")
    assert len(fake.tile_calls) == 1
    files = [name for _root, _dirs, names in os.walk(fake.cache_root) for name in names]
    assert "stamp.json" in files and any(name.startswith("t") for name in files)
    client.get("/galaxy")
    assert len(fake.tile_calls) == 1  # second view served from disk


def test_galaxy_without_shape_shows_hint(client, fake):
    fake.shape = None
    html = client.get("/galaxy").get_data(as_text=True)
    assert "density skeleton hasn" in html
    assert _scene(html)["hasShape"] is True  # from the tiles' own has_shape


def test_galaxy_quadrant_lists_its_sectors_nearest_first(client, fake):
    fake.sectors = [_placed(1, radius=900.0, name="Far"), _placed(2, radius=10.0, name="Near <b>x</b>"),
                    _placed(3, x=-100.0, name="Elsewhere")]
    resp = client.get("/galaxy?quadrant=i")
    html = resp.get_data(as_text=True)
    assert "<title>Galaxy Map: Quadrant I - " in html
    assert "Sectors in Quadrant I" in html
    assert "Elsewhere" not in html
    assert html.index("Near &lt;b&gt;x&lt;/b&gt;") < html.index(">Far<")
    assert '<a href="/galaxy" aria-current="page">Galaxy</a>' in html
    crumbs = re.search(r'<nav class="breadcrumbs".*?</nav>', html, re.S).group(0)
    assert '<a href="/galaxy">Galaxy</a>' in crumbs and "Quadrant I</span>" in crumbs
    assert '<a href="/galaxy#galaxy-table">All Quadrants</a>' in html


def test_galaxy_quadrant_sector_links_follow_the_sector_page(client, fake, app):
    html = client.get("/galaxy?quadrant=I").get_data(as_text=True)
    with app.test_request_context("/"):
        expected = page_url("sector", sector_id=1).replace("&", "&amp;")
    assert f'<a href="{expected}">Placed 001</a>' in html


def test_galaxy_quadrant_pages(client, fake):
    fake.sectors = [_placed(i) for i in range(1, 56)]
    page1 = client.get("/galaxy?quadrant=I").get_data(as_text=True)
    assert "Placed 050" in page1 and "Placed 051" not in page1
    assert 'href="/galaxy?quadrant=I&amp;page=2#galaxy-table"' in page1
    page2 = client.get("/galaxy?quadrant=I&page=2").get_data(as_text=True)
    assert "Placed 055" in page2 and "Placed 050" not in page2


def test_galaxy_quadrant_table_route_sorts_and_filters(client, fake):
    fake.sectors = [_placed(1, radius=900.0, name="Far"), _placed(2, radius=10.0, name="Near"),
                    _placed(3, x=-100.0, name="Elsewhere")]
    data = client.get("/table/galaxy-quadrant?quadrant=I&sort=name&order=desc&facets=1").get_json()
    assert [row[0]["text"] for row in data["rows"]] == ["Near", "Far"] and data["total"] == 2
    assert data["rows"][0][0]["href"].startswith("/sector/")
    assert data["facets"]["zone"] and all(o["label"].startswith("Zone ") for o in data["facets"]["zone"])
    nearest = client.get("/table/galaxy-quadrant?quadrant=I").get_json()
    assert [row[0]["text"] for row in nearest["rows"]] == ["Near", "Far"]
    assert client.get("/table/galaxy-quadrant").get_json()["total"] == 0


def test_galaxy_empty_quadrant(client, fake):
    html = client.get("/galaxy?quadrant=IV").get_data(as_text=True)
    assert "0 sectors" in html and "<em>None</em>" in html


def test_galaxy_unknown_quadrant_falls_back_to_summary(client, fake):
    html = client.get("/galaxy?quadrant=IX&db=other").get_data(as_text=True)
    assert "<h2 id=\"galaxy-table-heading\">Quadrants</h2>" in html
    assert fake.dbs == {DB}


def test_galaxy_api_failure_is_a_502_page(client, fake, monkeypatch):
    def boom(db):
        raise apiclient.ApiError("down")
    monkeypatch.setattr(apiclient, "get_galaxy_sectors", boom)
    resp = client.get("/galaxy")
    assert resp.status_code == 502
    assert resp.mimetype == "text/html"


# --- /galaxy/tiles -----------------------------------------------------------------------

def test_tiles_endpoint_serves_and_caches(client, fake):
    query = "/galaxy/tiles?tiles=12/2048/2048/2048,1/1/1/1"
    first = client.get(query)
    assert first.status_code == 200
    assert first.mimetype == "application/json"
    assert first.headers["Cache-Control"] == "no-store"
    assert first.headers["Content-Security-Policy"] == "default-src 'none'"
    body = first.get_json()
    assert set(body["tiles"]) == {"12/2048/2048/2048", "1/1/1/1"}
    assert body["cached"] == 0 and body["stamp"] == STAMP
    assert "history" in body  # no stamp sent: the browser gets the history

    second = client.get(query + f"&stamp={STAMP}").get_json()
    assert second["cached"] == 2
    assert second["tiles"] == body["tiles"]
    assert "history" not in second
    assert len(fake.tile_calls) == 1
    assert fake.dbs == {DB}


def test_tiles_endpoint_ignores_db_param(client, fake):
    client.get("/galaxy/tiles?tiles=1/1/1/1&db=someone_else")
    assert fake.dbs == {DB}


def test_tiles_endpoint_rejects_bad_keys(client, fake):
    resp = client.get("/galaxy/tiles?tiles=99/0/0/0")
    assert resp.status_code == 400
    assert "error" in resp.get_json()
    too_many = ",".join(f"12/{i}/0/0" for i in range(tilecache.MAX_TILES_PER_REQUEST + 1))
    assert client.get(f"/galaxy/tiles?tiles={too_many}").status_code == 400


def test_tiles_endpoint_api_failure_is_json(client, fake):
    fake.fail_tiles = apiclient.ApiError("secret detail")
    resp = client.get("/galaxy/tiles?tiles=1/1/1/1")
    assert resp.status_code == 502
    assert resp.mimetype == "application/json"
    assert "secret detail" not in resp.get_data(as_text=True)


# --- /galaxy/stage -----------------------------------------------------------------------

def test_stage_endpoint_serves_and_caches(client, fake):
    first = client.get("/galaxy/stage?at=243.7.14.0")
    assert first.status_code == 200
    assert first.headers["Cache-Control"] == "no-store"
    assert first.get_json()["at"] == "243.7.14.0"
    assert first.get_json()["stamp"] == STAMP
    galaxy = client.get("/galaxy/stage").get_json()
    assert galaxy["child_m"] == 243
    client.get("/galaxy/stage?at=243.7.14.0")
    client.get("/galaxy/stage")
    assert fake.stage_calls == ["243.7.14.0", None]
    assert fake.dbs == {DB}


def test_stage_endpoint_rejects_bad_keys(client, fake):
    for bad in ("81.0.0.0", "243.0.99.0", "nonsense", "1.0.0.0"):
        resp = client.get(f"/galaxy/stage?at={bad}")
        assert resp.status_code == 400
        assert "error" in resp.get_json()
    assert fake.stage_calls == []


def test_galaxy_page_offers_territories_once_polities_exist(client, fake):
    html = client.get("/galaxy").get_data(as_text=True)
    assert 'data-action="territories"' in html
    assert _scene(html)["territoryPath"] == "/galaxy/territories"


@pytest.mark.parametrize("failure", [None, apiclient.ApiError("down")])
def test_galaxy_page_hides_territories_without_population(client, fake, monkeypatch, failure):
    """No polities generated (or the count can't be read): no button."""
    if failure is None:
        fake.polity_total = 0
    else:
        def fail(db, limit=None, offset=None):
            raise failure
        monkeypatch.setattr(apiclient, "get_polities", fail)
    html = client.get("/galaxy").get_data(as_text=True)
    assert 'data-action="territories"' not in html
    assert 'id="galaxymap3d-territories"' not in html
    assert _scene(html)["territoryPath"] is None


def test_galaxy_pick_mode_banner_and_sector_links(client, fake):
    """NAV's "Pick on Galaxy Map" (design doc section 9): a banner whose
    Cancel keeps the other endpoint, and sector links that carry the
    pick on to the Sector Map."""
    html = client.get("/galaxy?pick=to&from=system:12").get_data(as_text=True)
    banner = re.search(r'<p class="pick-banner".*?</p>', html, re.S).group(0)
    assert "Choosing a destination" in banner
    assert 'href="/nav?from=system:12"' in banner
    scene = _scene(html)
    assert scene["pick"] == "to"
    # The page's script (navpick.js) adds the pick to sector links and the
    # map's own URLs; the page hands it the other end and where Cancel goes.
    assert scene["sectorUrl"] == "/sector/{id}"
    assert scene["pickOther"] == "system:12" and scene["pickCancel"] == "/nav?from=system:12"
    # The Bookmarks menu keeps the pick too (NAV.40).
    menu = re.search(r'<details class="bookmarks-menu" data-bookmarks-menu[^>]*>', html).group(0)
    for attribute in ('data-pick="to"', 'data-keep-name="from"', 'data-keep-value="system:12"', 'data-nav-url="/nav"'):
        assert attribute in menu


def test_galaxy_pick_mode_for_a_start_without_a_destination(client, fake):
    html = client.get("/galaxy?pick=from").get_data(as_text=True)
    assert "Choosing a start" in html
    scene = _scene(html)
    assert scene["pick"] == "from" and scene["pickOther"] is None
    assert scene["sectorUrl"] == "/sector/{id}"


@pytest.mark.parametrize("query", ["", "?pick=sideways", "?pick=to&from=bogus"])
def test_galaxy_without_a_valid_pick_has_no_banner(client, fake, query):
    html = client.get("/galaxy" + query).get_data(as_text=True)
    assert "pick-banner" not in html
    scene = _scene(html)
    assert scene["pick"] is None and scene["pickOther"] is None and scene["pickCancel"] is None
    assert scene["sectorUrl"] == "/sector/{id}"
    assert "data-pick=" not in html


# --- /galaxy/territories -----------------------------------------------------------------

def test_territories_endpoint_names_each_polity(client, fake, monkeypatch):
    """The overlay's endpoint folds each polity's name, color and system
    count (from `/api/polities`) into `/api/territories`' ids."""
    monkeypatch.setattr(apiclient, "get_territories", lambda db: {
        "points": [{"id": 1, "polity_id": 7, "color": "#d94f4f", "x": 1.0, "y": 2.0, "z": 0.0}],
        "polities": [{"id": 7, "capital_pc": [1.0, 2.0, 0.0], "reach_ly": 40.0},
                     {"id": 8, "capital_pc": None, "reach_ly": 20.0}],
    })
    monkeypatch.setattr(apiclient, "get_polities", lambda db, limit=None, offset=None: {
        "items": [{"id": 7, "name": "The Union", "color": "#d94f4f", "government": "federation",
                   "system_count": 12}],
        "total": 1, "limit": limit, "offset": 0,
    })
    body = client.get("/galaxy/territories").get_json()
    assert body["points"][0]["id"] == 1
    union, unnamed = body["polities"]
    assert (union["name"], union["system_count"], union["reach_ly"]) == ("The Union", 12, 40.0)
    # A polity the names page didn't reach is still drawn, just unnamed.
    assert unnamed["id"] == 8 and unnamed["name"] is None


def test_territories_endpoint_reports_an_api_failure(client, fake, monkeypatch):
    def fail(db, **kwargs):
        raise apiclient.ApiError("down")

    monkeypatch.setattr(apiclient, "get_territories", fail)
    resp = client.get("/galaxy/territories")
    assert resp.status_code == 502
    assert "error" in resp.get_json()


# --- The NAV course overlay --------------------------------------------------------------

def test_galaxy_page_draws_a_course(client, fake, monkeypatch):
    """`?course=<from>,<to>` embeds the course `nav_page.galaxy_course`
    works out, and says what it is showing."""
    from planetgen.web import nav_page

    asked = []
    course = {
        "scope": "galaxy", "navUrl": "/nav?from=system:1&to=system:2", "sector": None,
        "points": [{"name": "Alpha", "role": "origin", "url": "/system/1", "x": 1.0, "y": 2.0, "z": 0.0},
                   {"name": "Omega", "role": "destination", "url": "/system/2", "x": 90.0, "y": 2.0, "z": 0.0}],
    }
    monkeypatch.setattr(nav_page, "galaxy_course", lambda f, t: asked.append((f, t)) or course)
    html = client.get("/galaxy?course=system:1,system:2").get_data(as_text=True)
    assert asked == [("system:1", "system:2")]
    assert _scene(html)["course"]["points"][1]["name"] == "Omega"
    assert "Showing the course from <strong>Alpha</strong> to <strong>Omega</strong>" in html
    assert 'href="/nav?from=system:1&amp;to=system:2"' in html


def test_galaxy_page_shows_the_course_legend_and_readout(client, fake, monkeypatch):
    """NAV.20: the page names both line styles, the course's numbers and the stops as links."""
    from planetgen.web import nav_page

    points = [{"name": "Alpha", "role": "origin", "url": "/system/1", "x": 1.0, "y": 2.0, "z": 0.0},
              {"name": "Waypoint", "role": "hop", "url": "/system/9", "x": 40.0, "y": 9.0, "z": 0.0},
              {"name": "Omega", "role": "destination", "url": "/system/2", "x": 90.0, "y": 2.0, "z": 0.0}]
    course = {
        "scope": "galaxy", "navUrl": "/nav", "sector": None, "points": points, "direct": [points[0], points[2]],
        "readout": {"distance": "290 ly", "course": "045 mark 000", "frame": "Galactic", "route_distance": "300 ly",
                    "times": [{"label": "Warp 1", "text": "290 years"}],
                    "stops": [{"name": "Waypoint", "url": "/system/9"}]},
    }
    monkeypatch.setattr(nav_page, "galaxy_course", lambda f, t: course)
    html = client.get("/galaxy?course=system:1,system:2").get_data(as_text=True)
    assert "Route via adjacent systems" in html and "Direct line" in html
    assert '<a href="/system/9">Waypoint</a>' in html
    assert "<dt>Course</dt><dd>045 mark 000</dd>" in html and "<dt>Warp 1</dt><dd>290 years</dd>" in html
    assert _scene(html)["course"]["direct"][1]["name"] == "Omega"
    course["direct"] = None
    assert "Direct line" not in client.get("/galaxy?course=system:1,system:2").get_data(as_text=True)


def test_galaxy_page_without_a_course(client, fake):
    assert _scene(client.get("/galaxy").get_data(as_text=True))["course"] is None
    # A malformed value asks nothing and draws nothing.
    assert _scene(client.get("/galaxy?course=system:1").get_data(as_text=True))["course"] is None


def test_galaxy_page_rings_the_sectors_a_job_made(client, fake, monkeypatch):
    """ADM.31: `?made=<since>,<until>` embeds the sectors a run made and says how many."""
    asked = []
    made = {"total": 3, "items": [
        {"id": n, "name": f"S{n}", "x": 10.0 * n, "y": 0.0, "z": 0.0, "ring_index": 5 + n, "layer_index": 0,
         "ring_slot_index": n} for n in (1, 2)]}
    monkeypatch.setattr(apiclient, "get_galaxy_made", lambda db, since, until=None: asked.append((since, until)) or made)
    html = client.get("/galaxy?made=100,200").get_data(as_text=True)
    assert asked == [(100.0, 200.0)]
    scene = _scene(html)["made"]
    assert scene["total"] == 3 and scene["sectors"] == [[6, 0, 1], [7, 0, 2]] and scene["points"][1]["x"] == 20.0
    assert "Showing the 3 sectors that run made" in html and "the first 2 are ringed" in html
    assert _scene(client.get("/galaxy").get_data(as_text=True))["made"] is None
    assert _scene(client.get("/galaxy?made=soon").get_data(as_text=True))["made"] is None


# --- /galaxy/locate ----------------------------------------------------------------------

def test_locate_endpoint_passes_the_name_through(client, fake, monkeypatch):
    calls = []

    def fake_locate(db, q):
        calls.append((db, q))
        return [{"kind": "sector", "id": 3, "name": "Belcana", "sector_id": 3, "sector_name": "Belcana",
                 "ring": 7, "layer": 0, "slot": 2}]

    monkeypatch.setattr(apiclient, "get_galaxy_locate", fake_locate)
    body = client.get("/galaxy/locate?q=  Belcana  ").get_json()
    assert [m["name"] for m in body["matches"]] == ["Belcana"]
    assert calls == [(DB, "Belcana")]
    # A blank query asks the API nothing.
    assert client.get("/galaxy/locate?q=   ").get_json() == {"matches": []}
    assert len(calls) == 1


def test_locate_endpoint_reports_an_api_failure(client, fake, monkeypatch):
    def fail(db, q):
        raise apiclient.ApiError("down")

    monkeypatch.setattr(apiclient, "get_galaxy_locate", fail)
    resp = client.get("/galaxy/locate?q=Belcana")
    assert resp.status_code == 502
    assert "error" in resp.get_json()


# --- /galaxy/nebula/<id>/shape -----------------------------------------------------------

def test_nebula_shape_endpoint_passes_the_mesh_through(client, fake, monkeypatch):
    calls = []

    def fake_shape(db, nebula_id, lod="low"):
        calls.append((db, nebula_id, lod))
        return {"id": nebula_id, "lod": lod, "vertices": [[0, 0, 0]], "faces": []}

    monkeypatch.setattr(apiclient, "get_nebula_shape", fake_shape)
    resp = client.get("/galaxy/nebula/5/shape?lod=full")
    assert resp.get_json()["lod"] == "full"
    assert client.get("/galaxy/nebula/5/shape").get_json()["lod"] == "low"
    assert calls == [(DB, 5, "full"), (DB, 5, "low")]


def test_nebula_shape_endpoint_reports_a_missing_nebula_and_an_api_failure(client, fake, monkeypatch):
    def missing(db, nebula_id, lod="low"):
        raise apiclient.NotFoundError("no such nebula: 9")

    monkeypatch.setattr(apiclient, "get_nebula_shape", missing)
    assert client.get("/galaxy/nebula/9/shape").status_code == 404

    def fail(db, nebula_id, lod="low"):
        raise apiclient.ApiError("down")

    monkeypatch.setattr(apiclient, "get_nebula_shape", fail)
    resp = client.get("/galaxy/nebula/9/shape")
    assert resp.status_code == 502
    assert "error" in resp.get_json()


def test_nebula_surroundings_endpoint_passes_the_stars_through(client, fake, monkeypatch):
    calls = []

    def fake_surroundings(db, nebula_id):
        calls.append((db, nebula_id))
        return {"radius_pc": 5.0, "half_width_pc": 30.0, "stars": []}

    monkeypatch.setattr(apiclient, "get_nebula_surroundings", fake_surroundings)
    assert client.get("/galaxy/nebula/6/surroundings").get_json()["half_width_pc"] == 30.0
    assert calls == [(DB, 6)]

    def missing(db, nebula_id):
        raise apiclient.NotFoundError("no such nebula: 9")

    monkeypatch.setattr(apiclient, "get_nebula_surroundings", missing)
    assert client.get("/galaxy/nebula/9/surroundings").status_code == 404


# --- Old CGI URLs ---------------------------------------------------------------------------

def test_old_galaxy_url_redirects(client):
    result = client.get("/galaxy.py?db=x&quadrant=II&page=3")
    assert result.status_code == 301
    assert result.headers["Location"] == "/galaxy?quadrant=II&page=3"


def test_old_galaxy_tiles_url_redirects(client):
    result = client.get("/galaxy_tiles.py", query_string={"db": "x", "tiles": "1/1/1/1,2/2/2/2", "stamp": STAMP})
    assert result.status_code == 301
    assert result.headers["Location"] == f"/galaxy/tiles?tiles=1%2F1%2F1%2F1%2C2%2F2%2F2%2F2&stamp={STAMP}"


# --- Real database, in-process ----------------------------------------------------------------

@pytest.fixture
def db_client(mysql_config, tmp_path, monkeypatch):
    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")
    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    application = create_app(RealConfig)
    application.testing = True
    return application.test_client()


def _place_sector(mysql_config, name, address=(5, 1, 20)):
    from planetgen.galaxy.geometry import galactic_radius_pc, sector_position_pc

    position = sector_position_pc(*address, 3.526)
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            return _db.insert_sector(conn, SpaceSector(name=name), galaxy_position={
                "center_x_pc": position[0], "center_y_pc": position[1], "center_z_pc": position[2],
                "galactic_radius_pc": galactic_radius_pc(position),
                "ring_index": address[0], "layer_index": address[1], "ring_slot_index": address[2],
            })
    finally:
        conn.close()


def test_real_galaxy_page_and_tiles(db_client, mysql_config, tmp_path):
    sector_id = _place_sector(mysql_config, "Real <Placed>")
    resp = db_client.get("/galaxy")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert "1 placed sector" in html
    scene = _scene(html)
    assert scene["fetchPath"] == "/galaxy/tiles"
    placed = [entry for tile in scene["initial"]["tiles"].values() for entry in tile["placed"]]
    assert any(entry["id"] == sector_id for entry in placed)

    quadrant = re.search(r'href="/galaxy\?quadrant=(I|II|III|IV)#galaxy-table">Quadrant \1</a></td>\s*<td data-label="Sectors">1<',
                         html).group(1)
    listing = db_client.get(f"/galaxy?quadrant={quadrant}").get_data(as_text=True)
    assert "Real &lt;Placed&gt;" in listing

    query = "/galaxy/tiles?tiles=12/2048/2048/2048,1/1/1/1"
    first = db_client.get(query).get_json()
    assert set(first["tiles"]) == {"12/2048/2048/2048", "1/1/1/1"}
    assert len(first["stamp"]) == 16
    second = db_client.get(query).get_json()
    assert second["cached"] == 2
    assert any(names for _root, _dirs, names in os.walk(tmp_path / "tiles"))
