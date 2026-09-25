# tests/test_web_galaxy.py

"""
The Flask-served Galaxy Map (`/galaxy`) and its tile JSON
(`/galaxy/tiles`), `html/web/galaxy_views.py`, plus the `galaxy.py`/
`galaxy_tiles.py` CGI shims that redirect there.

Most tests fake the data layer (the `apiclient` functions the view calls,
and the two `tilecache` uses for tiles), the same way
`test_web_pages.py` does. The tests at the bottom run against a real
throwaway database in-process and are skipped without a MySQL server.
"""

import json
import os
import re

import pytest

from api.app import create_app
from api.config import Config

import web  # noqa: F401 -- puts src/html/lib on sys.path
import apiclient  # noqa: E402
import tilecache  # noqa: E402
from stellarObjects import _db  # noqa: E402
from stellarObjects.spaceSector import SpaceSector  # noqa: E402
from web.helpers import page_url  # noqa: E402

from tests.webpage_support import run_page  # noqa: E402

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
        "shell_index": 10 + i, "shell_slot_index": 0, "system_count": i % 4,
    }


class FakeData:
    """Stands in for the API: the page's two reads and tilecache's two."""

    def __init__(self):
        self.sectors = [_placed(1), _placed(2, x=-100.0), _placed(3, x=-100.0, y=-50.0)]
        self.shape = {"edge_pc": 3.526, "outer_shell_index": 400}
        self.dbs = set()
        self.tile_calls = []
        self.fail_tiles = None

    def get_galaxy_sectors(self, db):
        self.dbs.add(db)
        return self.sectors

    def get_galaxy_shape(self, db):
        self.dbs.add(db)
        return self.shape

    def get_galaxy_changes(self, db, since=None):
        self.dbs.add(db)
        return {"stamp": STAMP, "state": "s", "full": since is None, "tiles": []}

    def get_galaxy_tiles(self, db, tile_keys, density_key=None):
        self.dbs.add(db)
        if self.fail_tiles:
            raise self.fail_tiles
        self.tile_calls.append((list(tile_keys), density_key))
        return {
            "tiles": {key: {"placed": [{"id": 1, "name": "T", "x": 1.0, "y": 0.0, "z": 0.0}], "planned": []}
                      for key in tile_keys},
            "density": {"key": density_key, "points": [{"x": 1.0, "y": 2.0, "z": 3.0}]} if density_key else None,
            "edge_pc": 3.526, "has_shape": True,
        }

    def auth_me(self, cookie_header):
        return None


@pytest.fixture
def fake(monkeypatch, tmp_path):
    data = FakeData()
    for name in ("get_galaxy_sectors", "get_galaxy_shape", "auth_me"):
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    monkeypatch.setattr(tilecache, "get_galaxy_changes", data.get_galaxy_changes)
    monkeypatch.setattr(tilecache, "get_galaxy_tiles", data.get_galaxy_tiles)
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
    assert scene["initial"]["stamp"] == STAMP
    assert scene["initial"]["tiles"]


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


def test_galaxy_empty_quadrant(client, fake):
    html = client.get("/galaxy?quadrant=IV").get_data(as_text=True)
    assert "No sectors placed in this Quadrant yet." in html


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
    query = "/galaxy/tiles?tiles=12/2048/2048/2048,1/1/1/1&density=1/1/1/1"
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
    assert second["cached"] == 3
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


# --- CGI shims ---------------------------------------------------------------------------

def test_galaxy_shim_redirects_get():
    result = run_page("http://127.0.0.1:9/api", "galaxy.py", query={"db": "x", "quadrant": "ii", "page": "3"})
    assert result.status_code == 301
    assert result.headers["Location"] == "/galaxy?quadrant=II&page=3"


def test_galaxy_shim_redirects_post_and_drops_junk():
    result = run_page("http://127.0.0.1:9/api", "galaxy.py", method="POST",
                      body={"db": "x", "quadrant": "<x>", "page": "2"})
    assert result.status_code == 301
    assert result.headers["Location"] == "/galaxy"


def test_galaxy_tiles_shim_redirects():
    result = run_page("http://127.0.0.1:9/api", "galaxy_tiles.py",
                      query={"db": "x", "tiles": "1/1/1/1,2/2/2/2", "stamp": STAMP})
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


def _place_sector(mysql_config, name, shell_index=5, shell_slot_index=100):
    from stellarObjects.galaxyGeometry import galactic_radius_pc, sector_position_pc

    position = sector_position_pc(shell_index, shell_slot_index, 3.526)
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            return _db.insert_sector(conn, SpaceSector(name=name), galaxy_position={
                "center_x_pc": position[0], "center_y_pc": position[1], "center_z_pc": position[2],
                "galactic_radius_pc": galactic_radius_pc(position),
                "shell_index": shell_index, "shell_slot_index": shell_slot_index,
                "vertices_pc": {"inner": [], "outer": []},
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

    quadrant = re.search(r'href="/galaxy\?quadrant=(I|II|III|IV)#galaxy-table">Quadrant \1</a></td>\s*<td>1<',
                         html).group(1)
    listing = db_client.get(f"/galaxy?quadrant={quadrant}").get_data(as_text=True)
    assert "Real &lt;Placed&gt;" in listing

    query = "/galaxy/tiles?tiles=12/2048/2048/2048,1/1/1/1&density=1/1/1/1"
    first = db_client.get(query).get_json()
    assert set(first["tiles"]) == {"12/2048/2048/2048", "1/1/1/1"}
    assert len(first["stamp"]) == 16
    second = db_client.get(query).get_json()
    assert second["cached"] == 3
    assert any(names for _root, _dirs, names in os.walk(tmp_path / "tiles"))
