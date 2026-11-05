# tests/test_web_galaxy.py

"""
The Flask-served Galaxy Map (`/galaxy`) and its node data
(`/galaxy/nodes`), `planetgen/web/galaxy_views.py`, plus the old `galaxy.py`
URL that redirects there.

Most tests fake the data layer (the `apiclient` functions the view calls,
and the `tilecache` uses for nodes), the same way
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
from tests.node_support import node_payload  # noqa: E402
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
        "shell_index": 10 + i, "band_index": 0, "wedge_index": 0, "system_count": i % 4,
    }


class FakeData:
    """Stands in for the API: the page's two reads and tilecache's two."""

    def __init__(self):
        self.sectors = [_placed(1), _placed(2, x=-100.0), _placed(3, x=-100.0, y=-50.0)]
        self.shape = {"edge_pc": 3.526, "outer_shell_index": 400}
        self.dbs = set()
        self.node_calls = []
        self.fail_nodes = None
        self.stage_calls = []
        self.window_calls = []
        self.cover_calls = []
        self.cover_cells = [[2, 1, 0, 0.25]]
        self.fail_cover = None
        self.stage_children = []
        self.stage_sectors = []
        self.polity_total = 2
        self.regions = [
            {"key": "9.5.1.7", "name": "Vang Zarh Region", "parent": "Vang Reach", "sector_count": 3, "capacity": 12, "share": 0.25},
            {"key": "9.5.1.8", "name": "Kell Oon Region", "parent": "Vang Reach", "sector_count": 1, "capacity": 12, "share": 1 / 12},
        ]

    def get_polities(self, db, limit=None, offset=None):
        self.dbs.add(db)
        return {"items": [], "total": self.polity_total, "limit": limit, "offset": offset}

    def get_superbubbles(self, db):
        self.dbs.add(db)
        return []

    def get_containers(self, db, kind=None):
        self.dbs.add(db)
        return []

    def get_galaxy_sectors(self, db):
        self.dbs.add(db)
        return self.sectors

    def get_galaxy_charted_regions(self, db):
        self.dbs.add(db)
        return {"regions": self.regions, "edge_pc": 3.526}

    def get_galaxy_shape(self, db):
        self.dbs.add(db)
        return self.shape

    def get_galaxy_changes(self, db, since=None):
        self.dbs.add(db)
        return {"stamp": STAMP, "state": "s", "full": since is None, "nodes": []}

    def get_galaxy_nodes(self, db, node_keys):
        self.dbs.add(db)
        if self.fail_nodes:
            raise self.fail_nodes
        self.node_calls.append(list(node_keys))
        return {"nodes": {key: {"aggregate": {"generated": 2}, "stars": [], "points": [], "clouds": []}
                          for key in node_keys}, "edge_pc": 3.526, "has_shape": True}

    def get_galaxy_stage(self, db, at=None):
        self.dbs.add(db)
        self.stage_calls.append(at)
        return {"at": at, "child_m": 243 if at is None else 27, "children": [dict(c, stats=dict(c["stats"]))
                                                                             for c in self.stage_children],
                "sectors": [dict(r, stats=dict(r.get("stats") or {})) for r in self.stage_sectors] if at and at.startswith("3.") else None,
                "edge_pc": 3.526}

    def get_galaxy_window(self, db, at, g):
        self.dbs.add(db)
        self.window_calls.append((at, g))
        return {"at": at, "g": g, "child_m": int(at.split(".")[0]), "cells": [[3, 5, 0]], "children": [],
                "sectors": None}

    def get_galaxy_nebula_cover(self, db, at=None, g=None):
        self.dbs.add(db)
        self.cover_calls.append((at, g))
        if self.fail_cover:
            raise self.fail_cover
        child_m = 243 if not at else int(at.split(".")[0]) // (1 if g else 3)
        return {"at": at, "child_m": child_m, "cells": [list(cell) for cell in self.cover_cells]}

    def auth_me(self, cookie_header):
        return None


@pytest.fixture
def fake(monkeypatch, tmp_path):
    data = FakeData()
    for name in ("get_galaxy_sectors", "get_galaxy_shape", "get_polities", "get_superbubbles", "get_containers",
                 "auth_me", "get_galaxy_nebula_cover", "get_galaxy_charted_regions"):
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    monkeypatch.setattr(tilecache, "get_galaxy_changes", data.get_galaxy_changes)
    monkeypatch.setattr(tilecache, "get_galaxy_nodes", data.get_galaxy_nodes)
    monkeypatch.setattr(tilecache, "get_galaxy_stage", data.get_galaxy_stage)
    monkeypatch.setattr(tilecache, "get_galaxy_window", data.get_galaxy_window)
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    monkeypatch.delenv("PLANETGEN_CACHE_MAX_MB", raising=False)
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

def test_galaxy_page_names_the_maps_modules_up_front(client, fake):
    """MAP.161: every module the map script imports is a modulepreload hint at the very URL the
    script's own import() uses, so the browser fetches the tree in parallel and each file once."""
    html = client.get("/galaxy").get_data(as_text=True)
    hinted = re.findall(r'<link rel="modulepreload" href="([^"]+)">', html)
    assert len(hinted) == len(set(hinted)) >= 20
    script = re.search(r'<script type="module" src="([^"]+galaxymap3d\.js[^"]*)"', html).group(1)
    query = script.split("?", 1)[1]
    assert all(url.endswith("?" + query) for url in hinted), "same ?v= as the page's script"
    assert any(url.startswith("/static/galaxystageview.js?") for url in hinted)
    assert html.index('rel="modulepreload"') < html.index(script)


def test_galaxy_page_renders_map_and_charted_regions(client, fake, app):
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
    # The regions where sectors are charted, with their names (quadrants are gone).
    assert '<h2 id="galaxy-table-heading">Charted regions</h2>' in html
    assert "Vang Zarh Region" in html and "Kell Oon Region" in html and "25.0%" in html and "Vang Reach" in html
    assert html.index("Vang Zarh Region") < html.index("Kell Oon Region")
    assert "Quadrant" not in html
    # The database never leaves the server.
    assert fake.dbs == {DB}
    assert not re.search(r'/galaxy[^"\s]*db=', html)


def test_galaxy_scene_data_points_at_new_urls(client, fake, app):
    scene = _scene(client.get("/galaxy").get_data(as_text=True))
    assert scene["fetchPath"] == "/galaxy/nodes"
    assert "db" not in scene
    assert scene["storageKey"] == DB  # the browser's node cache keeps its old key
    assert "{id}" in scene["sectorUrl"]
    with app.test_request_context("/"):
        assert scene["sectorUrl"].replace("{id}", "5") == page_url("sector", sector_id=5)
        assert scene["phenomenonUrl"].replace("{type}", "supernova_remnant").replace("{id}", "7") == page_url(
            "phenomenon", phenomenon_type="supernova_remnant", phenomenon_id=7)
        assert scene["systemUrl"].replace("{id}", "9") == page_url("system", system_id=9)
        assert scene["stagePath"] == url_for("web.galaxy_stage")
        assert scene["nearbyUrl"] == url_for("web.nearby")  # NAV.45: "What's within 10 pc" on the map's panels
    assert scene["initial"]["stamp"] == STAMP
    assert scene["initial"]["nodes"]
    # Visitors get no Generate buttons.
    assert scene["generate"] is None


def test_galaxy_scene_data_gives_an_admin_the_generate_target(client, fake, monkeypatch):
    monkeypatch.setattr(apiclient, "auth_me", lambda cookie_header: {"username": "admin"})
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")
    target = _scene(client.get("/galaxy").get_data(as_text=True))["generate"]
    assert target["url"] == "/admin/generate"
    assert target["csrfField"] == csrf.FIELD_NAME
    assert target["csrfToken"]


@pytest.fixture
def running(redis_server, key_prefix, monkeypatch):
    """A run changing the galaxy (PERF.149), started and ended by the test; the site reads Redis through it."""
    from planetgen.queue import redisqueue, runrecord
    from planetgen.queue.runrecord import Stage
    from planetgen.web.lib import galaxystatus
    conn = redisqueue.connect(redis_server)
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    monkeypatch.setattr(galaxystatus, "_client", lambda: conn)
    galaxystatus.clear_cache()

    class Run:
        run_id = None

        def start(self):
            self.run_id = runrecord.start_run(conn, "Galaxy", [Stage("mass", "Scatter the massive stars", units=4)])
            runrecord.begin_stage(conn, self.run_id, 0)
            galaxystatus.clear_cache()

        def end(self):
            runrecord.finish_run(conn, self.run_id, "done")
            galaxystatus.clear_cache()

    yield Run()
    galaxystatus.clear_cache()


def test_the_galaxy_map_is_closed_while_a_galaxy_generation_runs_and_opens_when_it_ends(client, fake, running):
    """PERF.149 (Boss, 4:55 PM EDT 2026-10-22: no galaxy view during generation): the page says a generation is running."""
    running.start()
    resp = client.get("/galaxy")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 503
    assert "The Galaxy Map is closed while the galaxy is charted" in html
    assert "The galaxy is being changed by a background job" in html and "Scatter the massive stars" in html
    assert 'id="galaxymap3d-canvas"' not in html and "galaxymap3d.js" not in html
    assert '<meta http-equiv="refresh" content="30">' in html
    assert resp.headers["Cache-Control"] == "no-store"
    assert "/admin/queue" not in html                              # the Queue page is for admins
    assert fake.node_calls == []                                    # nothing was read from the half-built galaxy
    running.end()
    reopened = client.get("/galaxy")
    assert reopened.status_code == 200 and 'id="galaxymap3d-canvas"' in reopened.get_data(as_text=True)


def test_an_admin_is_linked_to_the_queue_page_from_the_closed_galaxy_map(client, fake, running, monkeypatch):
    monkeypatch.setattr(apiclient, "auth_me", lambda cookie_header: {"username": "admin"})
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")
    running.start()
    html = client.get("/galaxy").get_data(as_text=True)
    assert '<a href="/admin/queue">See the job on the Queue page</a>' in html


def test_galaxy_first_view_is_written_to_the_disk_cache(client, fake):
    client.get("/galaxy")
    assert len(fake.node_calls) == 1
    files = [name for _root, _dirs, names in os.walk(fake.cache_root) for name in names]
    assert "stamp.json" in files and any(name.startswith("n") for name in files)
    client.get("/galaxy")
    assert len(fake.node_calls) == 1  # second view served from disk


def test_galaxy_without_shape_shows_hint(client, fake):
    fake.shape = None
    html = client.get("/galaxy").get_data(as_text=True)
    assert "density skeleton hasn" in html
    assert _scene(html)["hasShape"] is True  # from the nodes' own has_shape


def test_galaxy_quadrant_address_is_ignored(client, fake):
    html = client.get("/galaxy?quadrant=I").get_data(as_text=True)
    assert "<title>Galaxy Map - " in html and "Quadrant" not in html


def test_galaxy_with_no_charted_regions_says_so(client, fake):
    fake.regions = []
    html = client.get("/galaxy").get_data(as_text=True)
    assert "No sectors are charted yet." in html


def test_galaxy_api_failure_is_a_502_page(client, fake, monkeypatch):
    def boom(db):
        raise apiclient.ApiError("down")
    monkeypatch.setattr(apiclient, "get_galaxy_sectors", boom)
    resp = client.get("/galaxy")
    assert resp.status_code == 502
    assert resp.mimetype == "text/html"


# --- /galaxy/nodes (MAP.260) ----------------------------------------------------------

def test_nodes_endpoint_serves_and_caches(client, fake):
    query = "/galaxy/nodes?keys=243.4.1.2,1.5.0.3"
    first = client.get(query)
    assert first.status_code == 200 and first.mimetype == "application/vnd.planetgen.tiles"
    assert first.headers["Cache-Control"] == "no-store"
    assert first.headers["Content-Security-Policy"] == "default-src 'none'"
    body = node_payload(first)
    assert "history" in body  # no stamp sent: the browser gets the history
    assert set(body["nodes"]) == {"243.4.1.2", "1.5.0.3"} and body["cached"] == 0 and body["stamp"] == STAMP
    assert body["nodes"]["243.4.1.2"]["aggregate"] == {"generated": 2}
    second = node_payload(client.get(query + f"&stamp={STAMP}"))
    assert second["cached"] == 2 and second["nodes"] == body["nodes"] and "history" not in second
    assert len(fake.node_calls) == 1 and fake.dbs == {DB}


def test_nodes_endpoint_ignores_db_param_and_rejects_bad_keys(client, fake):
    client.get("/galaxy/nodes?keys=243.1.1.1&db=someone_else")
    assert fake.dbs == {DB}
    for keys in ("99.0.0.0", "243.1.1", "x"):
        resp = client.get(f"/galaxy/nodes?keys={keys}")
        assert resp.status_code == 400 and "error" in resp.get_json()
    many = ",".join(f"243.1.1.{i}" for i in range(tilecache.MAX_NODES_PER_REQUEST + 1))
    assert client.get(f"/galaxy/nodes?keys={many}").status_code == 400


def test_nodes_endpoint_api_failure_is_json(client, fake):
    fake.fail_nodes = apiclient.ApiError("secret detail")
    resp = client.get("/galaxy/nodes?keys=243.1.1.1")
    assert resp.status_code == 502 and resp.mimetype == "application/json"
    assert "secret detail" not in resp.get_data(as_text=True)


# --- /galaxy/stage -----------------------------------------------------------------------

def test_stage_endpoint_serves_and_caches(client, fake):
    first = client.get("/galaxy/stage?at=243.2.1.0")
    assert first.status_code == 200
    assert first.headers["Cache-Control"] == "no-store"
    assert first.get_json()["at"] == "243.2.1.0"
    assert first.get_json()["stamp"] == STAMP
    galaxy = client.get("/galaxy/stage").get_json()
    assert galaxy["child_m"] == 243
    client.get("/galaxy/stage?at=243.2.1.0")
    client.get("/galaxy/stage")
    assert fake.stage_calls == ["243.2.1.0", None]
    assert fake.dbs == {DB}
    # MAP.174: the stage is cached, the nebula cover is not (it changes without a sector changing).
    assert fake.cover_calls == [("243.2.1.0", None), (None, None), ("243.2.1.0", None), (None, None)]
    assert first.get_json()["nebula_cover"] == [[2, 1, 0, 0.25]]


def _child(shell, band, wedge):
    return {"shell": shell, "band": band, "wedge": wedge, "generated": 3,
            "stats": {"systems": 4, "expected_systems": 5.0, "stars": 6, "mean_age_gy": 2.0, "luminosity_sol": 7.0}}


def test_stage_carries_the_nebula_cover_of_filled_and_unfilled_cells(client, fake):
    """MAP.174: a generated child has its share in `stats`; a cell that holds no sector is in `nebula_cover` alone."""
    fake.stage_children = [_child(2, 1, 0), _child(5, 0, 1)]
    fake.cover_cells = [[2, 1, 0, 0.25], [7, 0, 0, 0.5]]
    payload = client.get("/galaxy/stage?at=243.2.1.0").get_json()
    stats = {(c["shell"], c["band"], c["wedge"]): c["stats"] for c in payload["children"]}
    assert stats[(2, 1, 0)]["nebula_cover"] == 0.25 and stats[(5, 0, 1)]["nebula_cover"] == 0.0
    assert stats[(2, 1, 0)]["stars"] == 6
    assert payload["nebula_cover"] == [[2, 1, 0, 0.25], [7, 0, 0, 0.5]]
    again = client.get("/galaxy/stage?at=243.2.1.0").get_json()
    assert again["nebula_cover"] == payload["nebula_cover"]
    fake.cover_cells = []
    assert client.get("/galaxy/stage?at=243.2.1.0").get_json()["nebula_cover"] == []


def test_stage_draws_without_the_cover_when_the_api_cannot_give_it(client, fake):
    fake.stage_children = [_child(2, 1, 0)]
    fake.fail_cover = apiclient.ApiError("secret detail")
    resp = client.get("/galaxy/stage?at=243.2.1.0")
    assert resp.status_code == 200
    payload = resp.get_json()
    assert payload["nebula_cover"] is None and "nebula_cover" not in payload["children"][0]["stats"]
    assert "secret detail" not in resp.get_data(as_text=True)


def test_window_carries_the_nebula_cover_and_a_window_of_sectors_none(client, fake):
    big = client.get("/galaxy/window?at=27.3.5.0&g=5").get_json()
    assert big["nebula_cover"] == [[2, 1, 0, 0.25]]
    assert fake.cover_calls[-1] == ("27.3.5.0", "5")
    fake.cover_cells = []
    resp = client.get("/galaxy/window?at=1.4.1.0&g=5")
    assert resp.status_code == 200, resp.get_data(as_text=True)
    assert resp.get_json()["nebula_cover"] == []


def test_stage_endpoint_rejects_bad_keys(client, fake):
    for bad in ("729.0.0.0", "243.0.99.0", "nonsense", "1.0.0.0"):
        resp = client.get(f"/galaxy/stage?at={bad}")
        assert resp.status_code == 400
        assert "error" in resp.get_json()
    assert fake.stage_calls == []


def test_window_endpoint_reads_big_cells_from_the_block_totals(client, fake):
    first = client.get("/galaxy/window?at=27.3.5.0&g=5")
    assert first.status_code == 200
    assert first.headers["Cache-Control"] == "no-store"
    window = first.get_json()
    assert window["at"] == "27.3.5.0" and window["g"] == 5 and window["child_m"] == 27
    assert window["stamp"] == STAMP and [3, 5, 0] in window["cells"] and window["children"] == []
    assert fake.window_calls == [("27.3.5.0", 5)] and fake.stage_calls == []
    top = client.get("/galaxy/window?at=243.2.1.0").get_json()
    assert top["g"] == 9 and fake.window_calls[1:] == [("243.2.1.0", 9)] and fake.stage_calls == []


def test_window_endpoint_builds_a_window_of_sectors_from_the_parent_stages(client, fake):
    first = client.get("/galaxy/window?at=1.4.1.0&g=5")
    assert first.status_code == 200
    window = first.get_json()
    assert window["at"] == "1.4.1.0" and window["g"] == 5 and window["child_m"] == 1
    assert window["stamp"] == STAMP and [4, 1, 0] in window["cells"] and window["children"] == []
    assert fake.stage_calls and all(call is None or call.startswith("3.") for call in fake.stage_calls)
    assert fake.window_calls == []
    calls = list(fake.stage_calls)
    client.get("/galaxy/window?at=1.4.1.0&g=5")
    assert fake.stage_calls == calls  # the stages are cached


def test_window_endpoint_rejects_bad_keys_and_widths(client, fake):
    for query in ("at=729.0.0.0", "at=243.0.99.0", "at=nonsense", "", "at=27.3.5.0&g=4", "at=27.3.5.0&g=x", "at=27.3.5.0&g=13"):
        resp = client.get(f"/galaxy/window?{query}")
        assert resp.status_code == 400, query
        assert "error" in resp.get_json()
    assert fake.stage_calls == [] and fake.window_calls == []


def _sector(shell, band, wedge, sector_id):
    return {"shell": shell, "band": band, "wedge": wedge, "id": sector_id, "name": f"Sector {sector_id}", "system_count": 2, "stats": {"stars": 5}}


def test_near_endpoint_lists_the_generated_sectors_round_a_point(client, fake):
    """MAP.335: the sectors a scene can be fetched for round the eye, cut from the cached stages like a window of sectors."""
    from planetgen.galaxy import spherical

    x, y, z = spherical.sector_position_pc(4, 1, 0, 3.526)
    fake.stage_sectors = [_sector(4, 1, 0, 11), _sector(5, 1, 0, 12), _sector(60, 0, 0, 13)]
    resp = client.get(f"/galaxy/near?x={x}&y={y}&z={z}&r=5")
    assert resp.status_code == 200, resp.get_data(as_text=True)
    assert resp.headers["Cache-Control"] == "no-store"
    near = resp.get_json()
    assert near["at"] == "1.4.1.0" and near["edge_pc"] == 3.526 and near["r"] == 5 and near["stamp"] == STAMP
    ids = {s["id"] for s in near["sectors"]}
    assert 11 in ids and 13 not in ids, "the sector holding the point is in, one far away is not"
    assert {"shell", "band", "wedge", "id", "name", "system_count", "stats"} <= set(near["sectors"][0])
    calls = list(fake.stage_calls)
    client.get(f"/galaxy/near?x={x}&y={y}&z={z}&r=5")
    assert fake.stage_calls == calls, "the stages are cached"
    assert fake.window_calls == []


def test_near_endpoint_holds_the_radius_and_says_nothing_outside_the_grid(client, fake):
    far = client.get("/galaxy/near?x=100&y=50&z=0&r=100000").get_json()
    assert far["r"] == pytest.approx(4.5 * 3.526)
    outside = client.get("/galaxy/near?x=1e15&y=0&z=0&r=5").get_json()
    assert outside["at"] is None and outside["sectors"] == []
    plain = client.get("/galaxy/near?x=100&y=50&z=0")
    assert plain.status_code == 200 and plain.get_json()["r"] == 0


def test_near_endpoint_rejects_bad_numbers(client, fake):
    for query in ("", "x=1&y=2", "x=a&y=2&z=3&r=1", "x=1&y=2&z=3&r=-1", "x=nan&y=2&z=3&r=1", "x=1&y=inf&z=3&r=1"):
        resp = client.get(f"/galaxy/near?{query}")
        assert resp.status_code == 400, query
        assert "error" in resp.get_json()
    assert fake.stage_calls == [] or all(call is None for call in fake.stage_calls)


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


def test_the_galaxy_map_draws_the_course_through_the_waypoints_in_its_address(client, fake, monkeypatch):  # NAV.49
    from planetgen.web import nav_page

    asked = []

    def through(raws):
        asked.append(raws)
        return {"scope": "galaxy", "points": [{"name": "A", "role": "origin", "url": "/system/1", "x": 1.0, "y": 0.0, "z": 0.0, "hopUnknown": False},
                                              {"name": "B", "role": "waypoint", "url": "/system/2", "x": 2.0, "y": 1.0, "z": 0.0, "hopUnknown": False},
                                              {"name": "C", "role": "destination", "url": "/system/3", "x": 3.0, "y": 0.0, "z": 0.0, "hopUnknown": False}],
                "sector": None, "via": True, "waypoints": [0, 1, 2], "navUrl": "/nav?from=system:1&to=system:3", "direct": None,
                "directs": [[{"x": 1.0, "y": 0.0, "z": 0.0}, {"x": 2.0, "y": 1.0, "z": 0.0}]], "adjusted": [],
                "readout": {"distance": "2 ly", "course": None, "frame": "Galactic", "route_distance": "2 ly", "times": [],
                            "stops": [{"name": "B", "url": "/system/2"}], "diagram": None}}

    monkeypatch.setattr(nav_page, "galaxy_course_through", through)
    page = client.get("/galaxy?via=system:1,system:2,system:3").get_data(as_text=True)
    assert asked == [["system:1", "system:2", "system:3"]]
    assert "Showing the course through your <strong>3 waypoints</strong>, from <strong>A</strong> to <strong>C</strong>" in page
    assert "course-readout" in page and "Direct line" in page
    asked.clear()
    assert "Showing the course" not in client.get("/galaxy?via=system:1").get_data(as_text=True) and asked == []
    # A course of its own wins over the waypoints.
    monkeypatch.setattr(nav_page, "galaxy_course", lambda a, b: None)
    client.get("/galaxy?course=system:1,system:2&via=system:1,system:2,system:3")
    assert asked == []


def test_galaxy_page_without_a_course(client, fake):
    assert _scene(client.get("/galaxy").get_data(as_text=True))["course"] is None
    # A malformed value asks nothing and draws nothing.
    assert _scene(client.get("/galaxy?course=system:1").get_data(as_text=True))["course"] is None


def test_galaxy_page_rings_the_sectors_a_job_made(client, fake, monkeypatch):
    """ADM.31: `?made=<since>,<until>` embeds the sectors a run made and says how many."""
    asked = []
    made = {"total": 3, "items": [
        {"id": n, "name": f"S{n}", "x": 10.0 * n, "y": 0.0, "z": 0.0, "shell_index": 5 + n, "band_index": 0,
         "wedge_index": n} for n in (1, 2)]}
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


# --- /galaxy/nebula/<id>/shape (MAP.142 part 4: gone; the maps draw a nebula from its points) -----------------------------

def test_the_site_no_longer_serves_a_nebula_shape_mesh(client, fake):
    assert client.get("/galaxy/nebula/5/shape?lod=full").status_code == 404


def test_nebula_volume_endpoint_serves_the_packed_bytes(client, fake, monkeypatch):
    """MAP.142: `/galaxy/nebula/<id>/volume` hands the browser the decoded bytes of the API's base64 volume, cacheable; a missing nebula is a
    404 and an API failure a 502, as JSON."""
    import base64

    calls = []

    def fake_volume(db, nebula_id):
        calls.append((db, nebula_id))
        return {"id": nebula_id, "format": "NVOL1", "level": 2, "count": 0, "data": base64.b64encode(b"NVOL\x01\x02").decode("ascii")}

    monkeypatch.setattr(apiclient, "get_nebula_volume", fake_volume)
    resp = client.get("/galaxy/nebula/5/volume")
    assert resp.status_code == 200 and resp.data == b"NVOL\x01\x02"
    assert "max-age" in resp.headers["Cache-Control"]
    assert calls == [(DB, "5")]

    def missing(db, nebula_id):
        raise apiclient.NotFoundError("no such nebula: 9")

    monkeypatch.setattr(apiclient, "get_nebula_volume", missing)
    assert client.get("/galaxy/nebula/9/volume").status_code == 404

    def fail(db, nebula_id):
        raise apiclient.ApiError("down")

    monkeypatch.setattr(apiclient, "get_nebula_volume", fail)
    resp = client.get("/galaxy/nebula/9/volume")
    assert resp.status_code == 502 and "error" in resp.get_json()


def test_nebula_surroundings_endpoint_passes_the_stars_through(client, fake, monkeypatch):
    calls = []

    def fake_surroundings(db, nebula_id):
        calls.append((db, nebula_id))
        return {"radius_pc": 5.0, "half_width_pc": 30.0, "stars": []}

    monkeypatch.setattr(apiclient, "get_nebula_surroundings", fake_surroundings)
    assert client.get("/galaxy/nebula/6/surroundings").get_json()["half_width_pc"] == 30.0
    assert calls == [(DB, "6")]

    def missing(db, nebula_id):
        raise apiclient.NotFoundError("no such nebula: 9")

    monkeypatch.setattr(apiclient, "get_nebula_surroundings", missing)
    assert client.get("/galaxy/nebula/9/surroundings").status_code == 404


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
    from planetgen.galaxy.frame import galactic_radius_pc
    from planetgen.galaxy.spherical import sector_position_pc

    position = sector_position_pc(*address, 3.526)
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            return _db.insert_sector(conn, SpaceSector(name=name), galaxy_position={
                "center_x_pc": position[0], "center_y_pc": position[1], "center_z_pc": position[2],
                "galactic_radius_pc": galactic_radius_pc(position),
                "shell_index": address[0], "band_index": address[1], "wedge_index": address[2],
            })
    finally:
        conn.close()


def test_real_galaxy_page_and_nodes(db_client, mysql_config, tmp_path):
    from planetgen.galaxy.drill import drill_chain_of, format_drill_key

    sector_id = _place_sector(mysql_config, "Real <Placed>")
    resp = db_client.get("/galaxy")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert "1 placed sector" in html
    scene = _scene(html)
    assert scene["fetchPath"] == "/galaxy/nodes"
    assert scene["initial"]["nodes"], "the opening nodes are embedded"
    assert not any(section in node for node in scene["initial"]["nodes"].values()
                   for section in ("planned", "filled")), "MAP.157: sections the page never reads"
    assert sector_id

    assert '<h2 id="galaxy-table-heading">Charted regions</h2>' in html
    assert "Quadrant" not in html

    root = format_drill_key(drill_chain_of(5, 1, 20)[0])
    query = f"/galaxy/nodes?keys={root},1.5.1.20"
    first = node_payload(db_client.get(query))
    assert set(first["nodes"]) == {root, "1.5.1.20"}
    assert len(first["stamp"]) == 16
    second = node_payload(db_client.get(query))
    assert second["cached"] == 2
    assert any(names for _root, _dirs, names in os.walk(tmp_path / "tiles"))


# --- /galaxy/uncharted/<shell>/<band>/<wedge>/scene (MAP.162) -----------------------------

def _uncharted_cell(**overrides):
    cell = {"shell_index": 3, "band_index": 0, "wedge_index": 1, "designation": "ABC", "center_pc": [100.0, 50.0, 0.0],
            "edge_pc": 4.0, "sector_id": None,
            "stars": [{"id": 1, "x": 101.0, "y": 50.0, "z": 0.0, "luminosity_sol": 900.0, "temperature_k": 20000.0,
                       "radius_sol": 5.0, "star_type": "B"}],
            "scattered": [{"id": 5, "kind": "black-hole", "subtype": "stellar", "type": "black_hole",
                           "x": 99.0, "y": 50.0, "z": 0.0}]}
    cell.update(overrides)
    return cell


def test_the_galaxy_page_gives_the_map_the_uncharted_scene_url(client, fake, app):
    scene = _scene(client.get("/galaxy").get_data(as_text=True))
    with app.test_request_context("/"):
        assert scene["unchartedSceneUrl"].replace("{shell}", "3").replace("{band}", "-2").replace("{wedge}", "7") == url_for(
            "web.galaxy_uncharted_scene", shell=3, band=-2, wedge=7)


def test_the_uncharted_scene_endpoint_serves_the_cells_objects_marked_uncharted(client, fake, monkeypatch):
    calls = []

    def cell(db, shell, band, wedge):
        calls.append((db, shell, band, wedge))
        return _uncharted_cell()

    monkeypatch.setattr(apiclient, "get_uncharted_sector", cell)
    resp = client.get("/galaxy/uncharted/3/0/1/scene")
    scene = resp.get_json()
    assert calls == [(DB, 3, 0, 1)]
    assert scene["uncharted"] is True and scene["designation"] == "ABC" and scene["halfEdgePc"] == 2.0
    assert len(scene["stars"]) == 1 and len(scene["clouds"]) == 1
    assert resp.headers["Cache-Control"] == "no-store"


def test_the_uncharted_scene_endpoint_sends_a_generated_cell_to_its_sectors_scene(client, fake, monkeypatch):
    monkeypatch.setattr(apiclient, "get_uncharted_sector", lambda *args: _uncharted_cell(sector_id=42))
    resp = client.get("/galaxy/uncharted/3/0/1/scene")
    assert resp.status_code == 302 and resp.headers["Location"].endswith("/sector/42/scene")


def test_the_uncharted_scene_endpoint_reports_a_bad_cell_and_an_api_failure(client, fake, monkeypatch):
    def missing(*args):
        raise apiclient.NotFoundError("no such cell")

    monkeypatch.setattr(apiclient, "get_uncharted_sector", missing)
    assert client.get("/galaxy/uncharted/3/0/1/scene").status_code == 404

    def off_the_grid(*args):
        raise apiclient.ApiError("wedge_index 3 out of range", status_code=400)

    monkeypatch.setattr(apiclient, "get_uncharted_sector", off_the_grid)
    assert client.get("/galaxy/uncharted/0/0/3/scene").status_code == 404

    def fail(*args):
        raise apiclient.ApiError("down")

    monkeypatch.setattr(apiclient, "get_uncharted_sector", fail)
    resp = client.get("/galaxy/uncharted/3/0/1/scene")
    assert resp.status_code == 502 and "error" in resp.get_json()


def test_a_page_whose_data_timed_out_asks_again_by_itself(client, fake, monkeypatch):
    """PERF.64: the 504 page of a busy database reloads itself and says when a client may retry."""
    def timed_out(*args, **kwargs):
        raise apiclient.ApiError("QUERY_TIMEOUT", status_code=504)

    monkeypatch.setattr(apiclient, "get_galaxy_sectors", timed_out)
    response = client.get("/galaxy")
    assert response.status_code == 504 and response.headers["Retry-After"] == "10"
    html = response.get_data(as_text=True)
    assert '<meta http-equiv="refresh" content="10">' in html and "database is busy" in html


def test_a_page_that_waited_on_a_locked_row_answers_503_and_asks_again(client, fake, monkeypatch):
    """ADM.53: a lock wait timeout (1205) is a 503 with `Retry-After`, not the generic 500."""
    import pymysql

    def locked(*args, **kwargs):
        raise pymysql.err.OperationalError(1205, "Lock wait timeout exceeded; try restarting transaction")

    monkeypatch.setattr(apiclient, "get_galaxy_sectors", locked)
    response = client.get("/galaxy")
    assert response.status_code == 503 and response.headers["Retry-After"] == "10"
    assert '<meta http-equiv="refresh" content="10">' in response.get_data(as_text=True)
