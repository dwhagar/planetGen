# tests/test_api.py

"""
End-to-end tests for the read-only Flask API (`src/html/api/`) against a real,
throwaway MySQL database (see `conftest.py`'s `mysql_config` fixture) --
every test here is skipped, not failed, when no MySQL test server is
configured/reachable.

Seeds a small, real sector (via `stellarObjects._db.save_sector`, the same
path `sectorGen.py` uses) rather than hand-building rows, so these tests
exercise the exact same read path (`queryDb.py`/`stellarObjects._db.load_sector`/
`load_star_system`) production traffic does.
"""

import math
import re
import time

import pytest

from api.app import create_app
from api.config import Config
from stellarObjects import _db, adminAuth
from stellarObjects._db import MySQLConfig
from stellarObjects.config import SystemConfig
from stellarObjects.galaxyViewport import tiles_intersecting_sphere
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem
from wikiClient import WikiClientPageExistsError, WikiPage

TEST_ADMIN_PASSWORD = "a-strong-test-password-123"


@pytest.fixture
def seeded_sector(mysql_config):
    """Saves one small sector (two systems, one G2V, one M5V) and returns
    (mysql_config, sector_id, system_ids) for tests to build requests
    against."""
    sector = SpaceSector("Test Sector", edge_ly=10.0)

    cfg_a = SystemConfig()
    cfg_a.STAR_TYPE = "G2V"
    cfg_a.PLANETS = False
    # Pinned False: left at its own default, `BINARY_SYSTEM` now rolls
    # real chance (StarSystem._should_generate_binary), and this fixture's
    # own tests rely on both systems staying single (is_binary == 0).
    cfg_a.BINARY_SYSTEM = False
    system_a = StarSystem(system_config=cfg_a)
    sector.add_system(system_a, position=(1.0, 1.0, 1.0), system_config=cfg_a)

    cfg_b = SystemConfig()
    cfg_b.STAR_TYPE = "M5V"
    cfg_b.PLANETS = False
    cfg_b.BINARY_SYSTEM = False
    system_b = StarSystem(system_config=cfg_b)
    sector.add_system(system_b, position=(-2.0, 0.5, 3.0), system_config=cfg_b)

    sector_id = _db.save_sector(sector, config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        rows = conn.execute(
            "SELECT id, name FROM star_systems WHERE sector_id = ? ORDER BY id", (sector_id,)
        ).fetchall()
    finally:
        conn.close()

    return mysql_config, sector_id, [row["id"] for row in rows]


@pytest.fixture
def client(mysql_config):
    # Subclasses the real Config (not a bare new class carrying only
    # MYSQL_CONFIG) so tests exercise the actual rate-limit configuration
    # production uses -- a from-scratch class here would silently leave
    # RATELIMIT_DEFAULT/RATELIMIT_STORAGE_URI unset, making the global
    # default limit untestable (confirmed by testing: Flask-Limiter falls
    # back to an unconfigured in-memory store with no default limit at
    # all rather than erroring, so this gap wouldn't show up as a
    # failure -- just as tests silently not covering what they claim to).
    #
    # WRITE_MYSQL_CONFIG/CONTROL_MYSQL_CONFIG both point at this same
    # throwaway database as MYSQL_CONFIG -- there's no privilege
    # separation to test here (that's a deployment/grants concern, not
    # application logic), and the control schema's tables never collide
    # with the content schema's own (see test_admin_auth.py's module
    # docstring), so one database serves for all three roles in tests.
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


@pytest.fixture
def default_admin_client(mysql_config, client):
    """
    A `client` already logged in as the seeded default `admin`/`password`
    admin -- `must_change_credentials` is still set, so this is the
    fixture for tests confirming that gate actually blocks write/admin
    routes (see `authz.require_admin(fresh=True)`), not for exercising
    the writes themselves (use `admin_client` for that).

    Flask's test client keeps its own cookie jar across requests made on
    the same `client` instance, so the session cookie `POST /api/auth/login`
    sets is carried automatically into every later request this fixture's
    caller makes -- no manual cookie plumbing needed in tests, unlike the
    real CGI admin pages (`html/login.py` etc.), which do that relay by
    hand precisely because a CGI script has no persistent client object.
    """
    adminAuth.bootstrap_control_schema(mysql_config)
    # The write endpoints intentionally run with ensure_schema=False (see
    # routes._write_conn) -- a real deployment's content schema already
    # exists by the time an admin account is in use. This throwaway test
    # database starts completely empty, so lay the content schema down
    # here the same way any first `sectorGen.py`/`migrateDb.py` run would.
    _db.get_connection(mysql_config).close()
    response = client.post("/api/auth/login", json={
        "username": adminAuth.DEFAULT_ADMIN_USERNAME, "password": adminAuth.DEFAULT_ADMIN_PASSWORD,
    })
    assert response.status_code == 200
    assert response.get_json()["must_change_credentials"] is True
    return client


@pytest.fixture
def admin_client(default_admin_client):
    """A `client` logged in and past the forced credential change -- ready
    to exercise real write/admin endpoints against."""
    response = default_admin_client.post("/api/auth/change-credentials", json={
        "current_password": adminAuth.DEFAULT_ADMIN_PASSWORD,
        "new_username": "test-admin",
        "new_password": TEST_ADMIN_PASSWORD,
    })
    assert response.status_code == 200
    assert response.get_json()["must_change_credentials"] is False
    return default_admin_client


FAKE_WIKI_CONFIG = {
    "wikijs": {"base_url": "http://wiki.invalid/wikijs", "api_token": "test-token", "configured": True},
    "mediawiki": {
        "base_url": "http://wiki.invalid/mediawiki", "username": "Bot@pg", "password": "pw", "configured": True,
    },
}
"""dict: A `Config.WIKI_CONFIG`-shaped value with both backends
"configured" -- an unreachable `http://wiki.invalid/...` `base_url` on
purpose, so a test that forgets to install `fake_wiki_client` fails loudly
(a real connection attempt/timeout) rather than silently passing."""


@pytest.fixture
def client_with_wiki(mysql_config):
    """Same as `client` above, but with both wiki backends reporting as
    configured (see `FAKE_WIKI_CONFIG`) -- for tests that exercise the
    `POST .../wiki` routes' happy path, which `client`'s own unconfigured
    `WIKI_CONFIG` (nothing set in this test run's environment/config.json)
    would otherwise always reject with 501 before ever reaching
    `wikiClient.WikiClient` at all."""
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        WIKI_CONFIG = FAKE_WIKI_CONFIG

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


@pytest.fixture
def admin_client_with_wiki(mysql_config, client_with_wiki):
    """`admin_client`'s own login/credential-change flow, replayed against
    `client_with_wiki` instead of `client` -- duplicated rather than
    parameterizing `admin_client` itself, since fixtures can't take
    runtime arguments; kept to these two extra lines beyond what
    `default_admin_client`/`admin_client` already do."""
    adminAuth.bootstrap_control_schema(mysql_config)
    _db.get_connection(mysql_config).close()
    response = client_with_wiki.post("/api/auth/login", json={
        "username": adminAuth.DEFAULT_ADMIN_USERNAME, "password": adminAuth.DEFAULT_ADMIN_PASSWORD,
    })
    assert response.status_code == 200
    response = client_with_wiki.post("/api/auth/change-credentials", json={
        "current_password": adminAuth.DEFAULT_ADMIN_PASSWORD,
        "new_username": "test-admin-wiki",
        "new_password": TEST_ADMIN_PASSWORD,
    })
    assert response.status_code == 200
    return client_with_wiki


class _FakeWikiClient:
    """Stand-in for `wikiClient.WikiClient` (monkeypatched over
    `api.routes.WikiClient`) so these tests never make a real network call
    against a wiki instance -- records every `create_page` call (backend/
    path/title/content) on the class itself and returns a canned
    `WikiPage`, matching the real `WikiClient.create_page` contract."""

    calls = []

    def __init__(self, backend, **kwargs):
        self.backend = backend

    def create_page(self, path, title, content, **kwargs):
        _FakeWikiClient.calls.append({"backend": self.backend, "path": path, "title": title, "content": content})
        return WikiPage(id=1, path=path, title=title, url=f"http://wiki.invalid/{self.backend}/{path}")


@pytest.fixture
def fake_wiki_client(monkeypatch):
    """Installs `_FakeWikiClient` in place of `api.routes.WikiClient` for
    one test, resetting its recorded calls first."""
    _FakeWikiClient.calls = []
    monkeypatch.setattr("api.routes.WikiClient", _FakeWikiClient)
    return _FakeWikiClient


def test_health_ok(mysql_config, client):
    # `client`'s own database starts completely empty (see `mysql_config`'s
    # docstring) -- lay the schema down first, the same way a real
    # deployment always has one by the time its API is ever queried (see
    # `default_admin_client`'s identical step), so this exercises the
    # "reachable and current" case the test's own name promises rather
    # than the "never migrated at all" case `test_health_ok_reports_an_
    # uninitialized_schema_without_erroring` below covers instead.
    _db.get_connection(mysql_config).close()

    response = client.get("/api/health")
    assert response.status_code == 200
    body = response.get_json()
    assert body["status"] == "ok"
    assert body["schema_current"] is True
    assert body["schema_version"] == _db.SCHEMA_VERSION


def test_health_reports_an_uninitialized_schema_without_erroring(client):
    """
    `client`'s own database starts completely empty -- no `schema_migrations`
    table at all, one step further back than "some migrations pending"
    (which `schema_row` coming back empty already handles). `/api/health`
    should still report 200 with `schema_current: False` and a helpful
    `detail` rather than the earlier bug of an uncaught table-doesn't-exist
    error surfacing as a bare 503 "unreachable" (which is meant for an
    actually-unreachable server, not a reachable one that's simply never
    been migrated -- see `test_unreachable_database_returns_503_not_a_crash`
    for that case).
    """
    response = client.get("/api/health")
    assert response.status_code == 200
    body = response.get_json()
    assert body["status"] == "ok"
    assert body["schema_current"] is False
    assert body["schema_version"] is None
    assert "not been initialized" in body["detail"]


def test_unreachable_database_returns_503_not_a_crash():
    """
    Regression guard for a fixed bug: `queryDb.open_readonly` raises a
    bare `SystemExit` (correct for its own CLI callers) when the
    configured MySQL server is unreachable. `SystemExit` is a
    `BaseException`, not an `Exception` -- left uncaught, it would
    propagate straight through Flask's request dispatch instead of
    becoming any HTTP response at all, from every route that opens a
    connection via `routes.get_db()`, not just `/health`'s own explicit
    check. `get_db()` now converts it into an `ApiError` (503) so it
    becomes an ordinary exception from that point on.

    Deliberately builds its own app against a config that can never
    connect (a closed port on localhost, so this fails fast and needs no
    real MySQL server at all -- unlike every other test in this module,
    this one always runs, never skipped).
    """
    unreachable = MySQLConfig(host="127.0.0.1", port=1, user="x", password="x", database="x")

    class UnreachableConfig(Config):
        MYSQL_CONFIG = unreachable
        WRITE_MYSQL_CONFIG = unreachable
        CONTROL_MYSQL_CONFIG = unreachable
        RATELIMIT_STORAGE_URI = "memory://"

    test_client = create_app(UnreachableConfig).test_client()

    health_response = test_client.get("/api/health")
    assert health_response.status_code == 503
    assert health_response.get_json()["status"] == "error"

    sectors_response = test_client.get("/api/sectors")
    assert sectors_response.status_code == 503
    assert "error" in sectors_response.get_json()


def test_sectors_lists_seeded_sector(client, seeded_sector):
    _config, sector_id, _system_ids = seeded_sector

    response = client.get("/api/sectors")
    assert response.status_code == 200
    body = response.get_json()
    assert body["total"] >= 1
    assert body["limit"] == 100
    assert body["offset"] == 0
    ids = {entry["id"] for entry in body["items"]}
    assert sector_id in ids


def test_sectors_pagination(client, seeded_sector):
    response = client.get("/api/sectors?limit=1&offset=0")
    assert response.status_code == 200
    body = response.get_json()
    assert len(body["items"]) == 1
    assert body["limit"] == 1
    assert body["offset"] == 0


def test_sectors_rejects_invalid_limit(client):
    response = client.get("/api/sectors?limit=not-a-number")
    assert response.status_code == 400
    assert "error" in response.get_json()


def test_sector_detail_found_and_not_found(client, seeded_sector):
    _config, sector_id, system_ids = seeded_sector

    response = client.get(f"/api/sectors/{sector_id}")
    assert response.status_code == 200
    body = response.get_json()
    assert body["name"] == "Test Sector"
    assert body["system_count"] == len(system_ids)
    assert {s["id"] for s in body["systems"]} == set(system_ids)
    assert body["systems"][0]["stars"]

    response = client.get("/api/sectors/999999999")
    assert response.status_code == 404
    assert "error" in response.get_json()


def test_systems_filters_by_star_type_and_sector(client, seeded_sector):
    _config, sector_id, system_ids = seeded_sector

    response = client.get(f"/api/systems?sector_id={sector_id}")
    assert response.status_code == 200
    body = response.get_json()
    assert body["total"] == len(system_ids)

    response = client.get(f"/api/systems?sector_id={sector_id}&star_type=G")
    assert response.status_code == 200
    body = response.get_json()
    assert body["total"] == 1
    assert body["items"][0]["is_binary"] == 0
    assert body["items"][0]["star_summary"].startswith("G2V")


def test_systems_rejects_invalid_sector_id(client):
    response = client.get("/api/systems?sector_id=not-an-int")
    assert response.status_code == 400
    assert "error" in response.get_json()


def test_systems_sector_id_none_matches_standalone_systems(client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector

    response = client.get("/api/systems?sector_id=none")
    assert response.status_code == 200
    body = response.get_json()
    assert not (set(system_ids) & {item["id"] for item in body["items"]})


def test_system_detail_found_and_not_found(client, seeded_sector):
    _config, sector_id, system_ids = seeded_sector

    response = client.get(f"/api/systems/{system_ids[0]}")
    assert response.status_code == 200
    body = response.get_json()
    assert body["id"] == system_ids[0]
    assert body["sector_id"] == sector_id
    assert "stars" in body and body["stars"]
    # No stored page text since schema v28 -- see /text and /sections.
    assert "markdown_content" not in body and "wikitext_content" not in body
    assert {s["id"] for s in body["sector_siblings"]} == set(system_ids)
    # Nearest same-sector systems come from live rows (ids, current
    # names), never the system itself, nearest first.
    neighbors = body["nearest_neighbors"]
    assert 0 < len(neighbors) <= min(3, len(system_ids) - 1)
    assert system_ids[0] not in {n["id"] for n in neighbors}
    assert {n["id"] for n in neighbors} <= set(system_ids)
    distances = [n["distance_ly"] for n in neighbors]
    assert distances == sorted(distances)
    # The fixture's two systems sit at (1, 1, 1) and (-2, 0.5, 3) ly.
    assert neighbors[0]["distance_ly"] == pytest.approx(math.sqrt(9 + 0.25 + 4), rel=1e-3)

    response = client.get("/api/systems/999999999")
    assert response.status_code == 404


def test_system_text_renders_both_formats_from_the_database(client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector
    name = client.get(f"/api/systems/{system_ids[0]}").get_json()["name"]

    wikitext = client.get(f"/api/systems/{system_ids[0]}/text").get_json()
    assert wikitext["format"] == "wikitext"
    assert wikitext["content"].startswith(f"= {name} =")
    assert "[[Category:Star Systems]]" in wikitext["content"]

    markdown = client.get(f"/api/systems/{system_ids[0]}/text?format=markdown").get_json()
    assert markdown["content"].startswith(f"# {name}")

    assert client.get(f"/api/systems/{system_ids[0]}/text?format=html").status_code == 400
    assert client.get("/api/systems/999999999/text").status_code == 404


def test_system_text_follows_a_rename(admin_client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector
    response = admin_client.patch(f"/api/systems/{system_ids[0]}", json={"name": "Renamed After Generation"})
    assert response.status_code == 200

    content = admin_client.get(f"/api/systems/{system_ids[0]}/text?format=markdown").get_json()["content"]
    assert content.startswith("# Renamed After Generation")


def test_system_sections_cover_every_body(client, mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.MOONS = True
    cfg.COMETS = True
    cfg.ASTEROID_BELT = True
    cfg.BINARY_SYSTEM = False
    system_id = _db.save_system(StarSystem(system_config=cfg), cfg, config=mysql_config)

    detail = client.get(f"/api/systems/{system_id}").get_json()
    sections = client.get(f"/api/systems/{system_id}/sections").get_json()
    assert detail["planets"] and detail["planets"][0]["moons"] and detail["belts"] and detail["comets"]

    assert "This system contains" in sections["overview"] or "no stellar objects" in sections["overview"]
    assert set(sections["stars"]) == {str(s["id"]) for s in detail["stars"]}
    assert set(sections["planets"]) == {str(p["id"]) for p in detail["planets"]}
    assert set(sections["moons"]) == {str(m["id"]) for p in detail["planets"] for m in p["moons"]}
    assert set(sections["belts"]) == {str(b["id"]) for b in detail["belts"]}
    assert set(sections["comets"]) == {str(c["id"]) for c in detail["comets"]}
    for planet in detail["planets"]:
        # A planet's own section stops before its moons' sections begin.
        assert "### " not in sections["planets"][str(planet["id"])]
        assert {"habitable", "inhabited", "life_stage"} <= set(planet)

    assert client.get("/api/systems/999999999/sections").status_code == 404


def test_systems_near_requires_radius(client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector

    response = client.get(f"/api/systems/{system_ids[0]}/near")
    assert response.status_code == 400

    response = client.get(f"/api/systems/{system_ids[0]}/near?radius=-5")
    assert response.status_code == 400

    response = client.get(f"/api/systems/{system_ids[0]}/near?radius=1000")
    assert response.status_code == 200
    body = response.get_json()
    assert isinstance(body, list)
    other_ids = {entry["id"] for entry in body}
    assert system_ids[1] in other_ids


def test_nav_returns_direct_course_and_route_for_same_sector(client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector

    response = client.get(f"/api/nav?from={system_ids[0]}&to={system_ids[1]}")
    assert response.status_code == 200
    body = response.get_json()
    assert body["scope"] == "sector"
    assert body["direct"]["distance_ly"] > 0
    assert set(body["direct"]) == {"distance_ly", "azimuth_deg", "altitude_deg"}
    assert [leg["warp_factor"] for leg in body["warp_times"]] == [1, 3, 6, 9]
    assert body["route"]["path"] == [system_ids[0], system_ids[1]]
    assert len(body["origin_position"]) == 3
    assert len(body["destination_position"]) == 3
    assert set(body["route"]["positions"]) == {str(system_ids[0]), str(system_ids[1])}


def test_nav_requires_from_and_to(client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector

    response = client.get(f"/api/nav?to={system_ids[0]}")
    assert response.status_code == 400

    response = client.get(f"/api/nav?from={system_ids[0]}")
    assert response.status_code == 400

    response = client.get(f"/api/nav?from=not-an-int&to={system_ids[0]}")
    assert response.status_code == 400


def test_nav_returns_404_for_unknown_system(client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector

    response = client.get(f"/api/nav?from={system_ids[0]}&to=999999999")
    assert response.status_code == 404
    assert "error" in response.get_json()


def test_nav_returns_direct_course_for_system_to_phenomenon(client, mysql_config):
    """
    Regression test for a real bug this endpoint's own manual testing
    caught: `route["positions"]` can have both a plain-int key (a system
    hop) and a `queryDb._phenomenon_nav_key`-shaped string key (a
    phenomenon endpoint) once phenomenon endpoints exist -- Flask's
    `jsonify` sorts dict keys by default, and comparing an int key
    against a string key mid-sort raised a 500 `TypeError` here before
    `api/routes._route_for_json` fixed it, confirmed directly against a
    running instance of this exact endpoint.
    """
    sector = SpaceSector("Phenomenon Nav Sector", edge_ly=10.0)
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)

    empty_vertices = {"inner": [], "outer": []}
    sector_id = _db.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": 0.0, "center_y_pc": 0.0, "center_z_pc": 0.0,
        "galactic_radius_pc": 0.0, "vertices_pc": empty_vertices,
    })

    from stellarObjects.nebulaData import Nebula
    nebula = Nebula(SystemConfig())
    conn = _db.get_connection(mysql_config)
    try:
        system_id = conn.execute(
            "SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)
        ).fetchone()["id"]
        nebula_id = _db.insert_nebula(conn, nebula, sector_id=sector_id, placement={
            "center_x_pc": 5.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 5.0,
        })
        conn.commit()
    finally:
        conn.close()

    response = client.get(
        f"/api/nav?from={system_id}&to={nebula_id}&to_kind=phenomenon&to_type=nebula"
    )
    assert response.status_code == 200
    body = response.get_json()
    assert body["scope"] == "galaxy"
    assert body["route"] is not None
    phenomenon_key = f"phenomenon:nebula:{nebula_id}"
    assert body["route"]["path"][-1] == phenomenon_key
    assert phenomenon_key in body["route"]["positions"]
    assert str(system_id) in body["route"]["positions"]


def test_nav_returns_400_when_system_has_no_sector(client, mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        config_cur = conn.execute("INSERT INTO system_configs (markdown) VALUES (0)")
        config_id = config_cur.lastrowid
        origin_cur = conn.execute(
            "INSERT INTO star_systems (sector_id, system_config_id, name, quadrant) VALUES (NULL, ?, ?, ?)",
            (config_id, "Adrift", "I"),
        )
        origin_id = origin_cur.lastrowid
        destination_cur = conn.execute(
            "INSERT INTO star_systems (sector_id, system_config_id, name, quadrant) VALUES (NULL, ?, ?, ?)",
            (config_id, "Alsoadrift", "I"),
        )
        destination_id = destination_cur.lastrowid
        conn.commit()
    finally:
        conn.close()

    response = client.get(f"/api/nav?from={origin_id}&to={destination_id}")
    assert response.status_code == 400
    assert "error" in response.get_json()


def test_databases_lists_the_seeded_schema(client, seeded_sector):
    mysql_config, sector_id, system_ids = seeded_sector

    response = client.get("/api/databases")
    assert response.status_code == 200
    body = response.get_json()
    entry = next((item for item in body["items"] if item["name"] == mysql_config.database), None)
    assert entry is not None
    assert entry["sector_count"] == 1
    assert entry["system_count"] == len(system_ids)


def test_db_query_param_selects_a_different_database(client, mysql_config):
    # `client`'s own Config is pinned to `mysql_config`'s database --
    # this creates a *second*, separate schema (sharing the "planetgen"
    # prefix `list_databases`/`resolve_database` filter by, same as
    # `conftest.py`'s own throwaway names) to prove `?db=` actually
    # switches connections rather than coincidentally matching the default.
    import uuid

    other_name = f"planetgen_test_{uuid.uuid4().hex[:16]}"
    other_config = _db.MySQLConfig(
        host=mysql_config.host, port=mysql_config.port,
        user=mysql_config.user, password=mysql_config.password, database=other_name,
    )
    admin_conn = _db.get_connection(
        _db.MySQLConfig(
            host=mysql_config.host, port=mysql_config.port,
            user=mysql_config.user, password=mysql_config.password, database="",
        ),
        ensure_schema=False,
    )
    try:
        admin_conn.execute(f"CREATE DATABASE `{other_name}`")

        sector = SpaceSector("Other DB Sector", edge_ly=5.0)
        cfg = SystemConfig()
        cfg.STAR_TYPE = "K5V"
        cfg.PLANETS = False
        system = StarSystem(system_config=cfg)
        sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
        other_sector_id = _db.save_sector(sector, config=other_config)

        response = client.get(f"/api/sectors/{other_sector_id}?db={other_name}")
        assert response.status_code == 200
        assert response.get_json()["name"] == "Other DB Sector"

        # The same id, without ?db=, resolves against client's own
        # default database (mysql_config's) -- schema-initialize it first
        # (conftest.py's mysql_config fixture only CREATEs the database,
        # it never applies the schema -- that normally happens the first
        # time something writes to it, e.g. seeded_sector's save_sector
        # call, which this test doesn't use) so the "not found there"
        # case below is a clean 404, not a missing-table error.
        _db.get_connection(mysql_config).close()
        response = client.get(f"/api/sectors/{other_sector_id}")
        assert response.status_code == 404
    finally:
        admin_conn.execute(f"DROP DATABASE IF EXISTS `{other_name}`")
        admin_conn.close()


def test_db_query_param_rejects_unknown_database(client):
    response = client.get("/api/sectors?db=not_a_real_schema")
    assert response.status_code == 404
    assert "error" in response.get_json()


def test_galaxy_sectors_excludes_unplaced_sectors(client, seeded_sector):
    # seeded_sector's own sector is never given a galaxy placement.
    response = client.get("/api/galaxy/sectors")
    assert response.status_code == 200
    assert response.get_json() == {"items": []}


def test_galaxy_view_placed_sector_includes_edge_ly(client, mysql_config):
    """`GET /api/galaxy/view`'s own "placed" tier (queryDb.
    galaxy_sectors_in_view) exposes each placed sector's own real
    edge_ly -- added so static/galaxymap3d.js can color a sector marker
    by its own true stellar density (system_count / edge_ly ** 3)
    relative to the real local average, rather than raw system count
    alone."""
    sector = SpaceSector("Density Test Sector", edge_ly=10.0)
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)

    empty_vertices = {"inner": [], "outer": []}
    _db.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": 0.0, "center_y_pc": 0.0, "center_z_pc": 0.0,
        "galactic_radius_pc": 0.0, "vertices_pc": empty_vertices,
    })

    response = client.get("/api/galaxy/view?cx=0&cy=0&cz=0&radius_pc=1000")
    assert response.status_code == 200
    placed = response.get_json()["placed"]
    assert len(placed) == 1
    assert placed[0]["edge_ly"] == pytest.approx(10.0)
    assert placed[0]["system_count"] == 1



def _place_sector(mysql_config, name, center_pc, edge_ly=10.0, shell_index=None, shell_slot_index=None):
    sector = SpaceSector(name, edge_ly=edge_ly)
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=cfg), position=(0.0, 0.0, 0.0), system_config=cfg)
    position = {
        "center_x_pc": center_pc[0], "center_y_pc": center_pc[1], "center_z_pc": center_pc[2],
        "galactic_radius_pc": math.dist(center_pc, (0.0, 0.0, 0.0)),
        "vertices_pc": {"inner": [], "outer": []},
    }
    if shell_index is not None:
        position["shell_index"] = shell_index
        position["shell_slot_index"] = shell_slot_index
    return _db.save_sector(sector, config=mysql_config, galaxy_position=position)


def test_galaxy_tiles_returns_placed_sectors_by_cube(client, mysql_config):
    """`GET /api/galaxy/tiles` puts each placed sector in exactly the tile
    whose half-open box holds its center, and lists planned slots only
    for small tiles."""
    _place_sector(mysql_config, "Near Core", (5.0, 5.0, 5.0))
    _place_sector(mysql_config, "Far Out", (9000.0, -3000.0, 20.0))

    near = tiles_intersecting_sphere(12, (5.0, 5.0, 5.0), 0.0)[0]
    far = tiles_intersecting_sphere(12, (9000.0, -3000.0, 20.0), 0.0)[0]
    whole = "0/0/0/0"
    response = client.get(f"/api/galaxy/tiles?tiles={near},{far},{whole}")
    assert response.status_code == 200
    body = response.get_json()
    tiles = body["tiles"]
    assert [s["name"] for s in tiles[near]["placed"]] == ["Near Core"]
    assert [s["name"] for s in tiles[far]["placed"]] == ["Far Out"]
    assert sorted(s["name"] for s in tiles[whole]["placed"]) == ["Far Out", "Near Core"]
    assert tiles[near]["placed"][0]["system_count"] == 1
    assert tiles[near]["placed"][0]["edge_ly"] == pytest.approx(10.0)
    # No galaxy shape here, so small tiles list every slot unfiltered;
    # the level-0 tile is far too big to list slots at all.
    assert tiles[near]["planned"]
    assert tiles[whole]["planned"] == []
    assert body["density"] is None
    assert body["has_shape"] is False


def test_galaxy_tiles_rejects_bad_requests(client, mysql_config):
    _db.get_connection(mysql_config).close()
    assert client.get("/api/galaxy/tiles?tiles=13/0/0/0").status_code == 400
    assert client.get("/api/galaxy/tiles?tiles=1/2/0/0").status_code == 400
    assert client.get("/api/galaxy/tiles?tiles=nonsense").status_code == 400
    too_many = ",".join(f"12/{i}/0/0" for i in range(129))
    assert client.get(f"/api/galaxy/tiles?tiles={too_many}").status_code == 400
    empty = client.get("/api/galaxy/tiles?tiles=&density=1/0/0/0")
    assert empty.status_code == 200
    assert empty.get_json()["density"] == {"key": "1/0/0/0", "points": []}


def test_galaxy_stamp_changes_when_a_sector_is_placed(client, mysql_config):
    _db.get_connection(mysql_config).close()
    first = client.get("/api/galaxy/stamp").get_json()["stamp"]
    assert re.fullmatch(r"[0-9a-f]{16}", first)
    assert client.get("/api/galaxy/stamp").get_json()["stamp"] == first
    _place_sector(mysql_config, "New Sector", (1.0, 2.0, 3.0))
    assert client.get("/api/galaxy/stamp").get_json()["stamp"] != first

def test_search_returns_facets_and_matches_a_class_tag(client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector

    response = client.get("/api/search")
    assert response.status_code == 200
    body = response.get_json()
    assert set(body["facets"]) == {
        "type", "spectral", "luminosity", "class", "body", "life",
        "moon_class", "moon_body", "moon_life", "density",
    }
    assert body["results"]["stars"] is None  # no filter active yet -- no reason to run

    response = client.get("/api/search?spectral=G")
    assert response.status_code == 200
    body = response.get_json()
    assert body["results"]["stars"] is not None
    assert all(row["star_system_id"] in system_ids for row in body["results"]["stars"]["rows"])


def test_search_text_queries(client, seeded_sector):
    _config, _sector_id, _system_ids = seeded_sector

    for query_param in ("sector_q=Test", "system_q=Test", "star_q=Test", "planet_q=Test", "moon_q=Test"):
        response = client.get(f"/api/search?{query_param}")
        assert response.status_code == 200
        body = response.get_json()
        assert "results" in body


def test_unmatched_route_returns_json_404(client):
    response = client.get("/api/no-such-route")
    assert response.status_code == 404
    assert response.get_json() == {"error": "not found"}


# ---------------------------------------------------------------------
# Authentication (see auth.py/authz.py) -- login, the forced default-
# credential change, and API keys.
# ---------------------------------------------------------------------

def test_login_wrong_password_is_generic_401(mysql_config, client):
    adminAuth.bootstrap_control_schema(mysql_config)
    response = client.post("/api/auth/login", json={"username": "admin", "password": "wrong"})
    assert response.status_code == 401
    unknown_user_response = client.post("/api/auth/login", json={"username": "no-such-admin", "password": "wrong"})
    assert unknown_user_response.status_code == 401
    assert response.get_json()["error"] == unknown_user_response.get_json()["error"]


def test_login_success_sets_cookie_and_reports_must_change_credentials(default_admin_client):
    response = default_admin_client.get("/api/auth/me")
    assert response.status_code == 200
    body = response.get_json()
    assert body["username"] == "admin"
    assert body["must_change_credentials"] is True


def test_login_sets_cookie_scoped_to_root_path_not_api(mysql_config, client):
    """
    Regression test: the session cookie's `Path` attribute must be `/`,
    not `/api`. The documented Apache deployment (docs/apache-
    deployment.md, examples/apache/planetgen.conf.example) serves the CGI
    admin pages (admin.py, login.py, etc.) at the site root under the
    same `DocumentRoot` the API is mounted under at `/api` -- per RFC 6265
    path matching, a cookie scoped to `/api` is never attached by the
    browser to a request for `/admin.py`, so a `Path=/api` cookie would
    lock a user out of every admin page immediately after a successful
    login. Flask's own test client can't reproduce that failure directly
    (it only ever talks to `/api/...` routes here, never `/admin.py`), so
    this asserts the cookie's scope directly instead.
    """
    adminAuth.bootstrap_control_schema(mysql_config)
    response = client.post("/api/auth/login", json={
        "username": adminAuth.DEFAULT_ADMIN_USERNAME, "password": adminAuth.DEFAULT_ADMIN_PASSWORD,
    })
    assert response.status_code == 200
    set_cookie_headers = response.headers.get_all("Set-Cookie")
    assert any("path=/;" in h.lower() or h.lower().endswith("path=/") for h in set_cookie_headers), set_cookie_headers


def test_me_requires_auth(client):
    response = client.get("/api/auth/me")
    assert response.status_code == 401


def test_write_route_requires_auth(client):
    response = client.post("/api/sectors", json={"name": "Test", "edge_ly": 10.0})
    assert response.status_code == 401


def test_write_route_blocked_until_default_credentials_are_changed(default_admin_client):
    response = default_admin_client.post("/api/sectors", json={"name": "Test", "edge_ly": 10.0})
    assert response.status_code == 403


def test_change_credentials_requires_current_password(default_admin_client):
    response = default_admin_client.post("/api/auth/change-credentials", json={
        "current_password": "wrong", "new_username": "someone", "new_password": TEST_ADMIN_PASSWORD,
    })
    assert response.status_code == 400


def test_change_credentials_rejects_weak_password(default_admin_client):
    response = default_admin_client.post("/api/auth/change-credentials", json={
        "current_password": adminAuth.DEFAULT_ADMIN_PASSWORD, "new_username": "someone", "new_password": "short",
    })
    assert response.status_code == 400


def test_change_credentials_success_clears_the_flag_and_reissues_session(admin_client):
    response = admin_client.get("/api/auth/me")
    assert response.status_code == 200
    body = response.get_json()
    assert body["username"] == "test-admin"
    assert body["must_change_credentials"] is False


def test_logout_ends_the_session(admin_client):
    response = admin_client.post("/api/auth/logout")
    assert response.status_code == 200
    response = admin_client.get("/api/auth/me")
    assert response.status_code == 401


def test_api_key_create_list_revoke(admin_client):
    response = admin_client.get("/api/auth/api-keys")
    assert response.status_code == 200
    assert response.get_json()["items"] == []

    response = admin_client.post("/api/auth/api-keys", json={"label": "ci key"})
    assert response.status_code == 201
    created = response.get_json()
    assert created["key"].startswith("pg_")

    response = admin_client.get("/api/auth/api-keys")
    items = response.get_json()["items"]
    assert len(items) == 1
    assert items[0]["id"] == created["id"]
    assert "key" not in items[0]
    assert items[0]["revoked_at"] is None

    # The raw key works as a Bearer credential against a write route.
    response = admin_client.post(
        "/api/sectors", json={"name": "Via API Key", "edge_ly": 8.0},
        headers={"Authorization": f"Bearer {created['key']}"},
    )
    assert response.status_code == 201

    response = admin_client.delete(f"/api/auth/api-keys/{created['id']}")
    assert response.status_code == 200

    response = admin_client.post(
        "/api/sectors", json={"name": "Should fail", "edge_ly": 8.0},
        headers={"Authorization": f"Bearer {created['key']}"},
    )
    assert response.status_code == 401

    response = admin_client.delete(f"/api/auth/api-keys/{created['id']}")
    assert response.status_code == 404


def test_login_is_rate_limited(mysql_config, client):
    adminAuth.bootstrap_control_schema(mysql_config)
    for _ in range(10):
        response = client.post("/api/auth/login", json={"username": "admin", "password": "wrong"})
        assert response.status_code == 401

    response = client.post("/api/auth/login", json={"username": "admin", "password": "wrong"})
    assert response.status_code == 429


# ---------------------------------------------------------------------
# Write endpoints -- real inserts/updates/deletes (see routes.py),
# authenticated via `admin_client`.
# ---------------------------------------------------------------------

def test_create_sector_then_get_update_delete(admin_client):
    response = admin_client.post("/api/sectors", json={"name": "Test", "edge_ly": 10.0})
    assert response.status_code == 201
    sector_id = response.get_json()["id"]

    response = admin_client.get(f"/api/sectors/{sector_id}")
    assert response.status_code == 200
    assert response.get_json()["name"] == "Test"
    assert response.get_json()["edge_ly"] == pytest.approx(10.0)

    response = admin_client.patch(f"/api/sectors/{sector_id}", json={"name": "Renamed"})
    assert response.status_code == 200
    assert admin_client.get(f"/api/sectors/{sector_id}").get_json()["name"] == "Renamed"

    response = admin_client.patch("/api/sectors/999999999", json={"name": "Nope"})
    assert response.status_code == 404

    response = admin_client.delete(f"/api/sectors/{sector_id}")
    assert response.status_code == 200
    assert admin_client.get(f"/api/sectors/{sector_id}").status_code == 404

    response = admin_client.delete(f"/api/sectors/{sector_id}")
    assert response.status_code == 404


def test_create_sector_validates_request_body(admin_client):
    response = admin_client.post("/api/sectors", json={"name": "Test"})
    assert response.status_code == 400

    response = admin_client.post("/api/sectors", json={"name": "", "edge_ly": 10.0})
    assert response.status_code == 400

    response = admin_client.post("/api/sectors", json={"name": "Test", "edge_ly": -1})
    assert response.status_code == 400

    response = admin_client.post("/api/sectors", json={"name": "Test", "edge_ly": 10.0, "bogus": 1})
    assert response.status_code == 400

    response = admin_client.post("/api/sectors", data="not json", content_type="text/plain")
    assert response.status_code == 400


def test_delete_sector_detaches_rather_than_deletes_its_systems(admin_client, seeded_sector):
    _config, sector_id, system_ids = seeded_sector

    response = admin_client.delete(f"/api/sectors/{sector_id}")
    assert response.status_code == 200

    response = admin_client.get(f"/api/systems/{system_ids[0]}")
    assert response.status_code == 200
    assert response.get_json()["sector_id"] is None


def test_delete_system_bumps_its_sectors_modified_at(admin_client, seeded_sector):
    # v27: a deleted row leaves no timestamp behind, so its sector's
    # modified_at is what tells a cache the sector changed.
    config, sector_id, system_ids = seeded_sector

    def sector_modified_at():
        conn = _db.get_connection(config)
        try:
            return conn.execute("SELECT modified_at FROM sectors WHERE id = ?", (sector_id,)).fetchone()["modified_at"]
        finally:
            conn.close()

    before = sector_modified_at()
    time.sleep(0.05)
    response = admin_client.delete(f"/api/systems/{system_ids[0]}")
    assert response.status_code == 200
    assert sector_modified_at() > before


def test_create_system_generates_and_persists_a_standalone_system(admin_client):
    response = admin_client.post("/api/systems", json={"star_type": "G2V", "planets": False})
    assert response.status_code == 201
    system_id = response.get_json()["id"]

    response = admin_client.get(f"/api/systems/{system_id}")
    assert response.status_code == 200
    body = response.get_json()
    assert body["sector_id"] is None
    assert body["stars"][0]["star_type"].startswith("G2V")


def test_create_system_rejects_unrecognized_and_invalid_fields(admin_client):
    response = admin_client.post("/api/systems", json={"sector_id": 1})
    assert response.status_code == 400

    response = admin_client.post("/api/systems", json={"age": "ancient"})
    assert response.status_code == 400

    response = admin_client.post("/api/systems", json={"num_orbits": -1})
    assert response.status_code == 400


def test_update_and_delete_system(admin_client):
    response = admin_client.post("/api/systems", json={"star_type": "M5V", "planets": False})
    system_id = response.get_json()["id"]

    response = admin_client.patch(f"/api/systems/{system_id}", json={"name": "Renamed World"})
    assert response.status_code == 200
    assert admin_client.get(f"/api/systems/{system_id}").get_json()["name"] == "Renamed World"

    response = admin_client.patch(f"/api/systems/{system_id}", json={"star_type": "G2V"})
    assert response.status_code == 400  # unrecognized field for PATCH

    response = admin_client.delete(f"/api/systems/{system_id}")
    assert response.status_code == 200
    assert admin_client.get(f"/api/systems/{system_id}").status_code == 404


def test_write_endpoints_are_rate_limited_more_tightly_than_the_default(admin_client):
    # WRITE_RATE_LIMIT is 10/minute -- the 11th write in one minute must
    # be rejected with 429, well before the 50/hour global default would
    # ever kick in.
    for _ in range(10):
        response = admin_client.delete("/api/sectors/999999999")
        assert response.status_code == 404

    response = admin_client.delete("/api/sectors/999999999")
    assert response.status_code == 429
    body = response.get_json()
    assert body["error"] == "rate limit exceeded"


# ---------------------------------------------------------------------
# Wiki publishing (schema.sql's "v22" header note) -- POST .../wiki, the
# GET /api/wiki-config it's gated on, and PATCH /api/sectors/<id>'s own
# manual `wiki_url` field (html/admin.py's manual-link admin section).
# ---------------------------------------------------------------------

def test_wiki_config_reports_unconfigured_by_default(client):
    # Neither backend has any PLANETGEN_WIKIJS_*/PLANETGEN_MEDIAWIKI_*
    # env var or config.json section set in this test run.
    response = client.get("/api/wiki-config")
    assert response.status_code == 200
    assert response.get_json() == {"wikijs": False, "mediawiki": False}


def test_wiki_config_reports_configured_backends(client_with_wiki):
    response = client_with_wiki.get("/api/wiki-config")
    assert response.status_code == 200
    assert response.get_json() == {"wikijs": True, "mediawiki": True}


def test_upload_system_wiki_requires_auth(seeded_sector, client):
    _config, _sector_id, system_ids = seeded_sector
    response = client.post(f"/api/systems/{system_ids[0]}/wiki", json={"backend": "wikijs", "path": "x"})
    assert response.status_code == 401


def test_upload_system_wiki_returns_501_when_backend_not_configured(admin_client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector
    response = admin_client.post(f"/api/systems/{system_ids[0]}/wiki", json={"backend": "wikijs", "path": "x"})
    assert response.status_code == 501


def test_upload_system_wiki_rejects_invalid_body(admin_client_with_wiki, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector
    system_id = system_ids[0]

    response = admin_client_with_wiki.post(f"/api/systems/{system_id}/wiki", json={"backend": "confluence"})
    assert response.status_code == 400

    # wikijs requires a path (it addresses a page separately from its title)
    response = admin_client_with_wiki.post(f"/api/systems/{system_id}/wiki", json={"backend": "wikijs"})
    assert response.status_code == 400

    response = admin_client_with_wiki.post(
        f"/api/systems/{system_id}/wiki", json={"backend": "wikijs", "path": "x", "bogus": 1}
    )
    assert response.status_code == 400


def test_upload_system_wiki_wikijs_persists_url_on_the_matching_column(
    admin_client_with_wiki, seeded_sector, fake_wiki_client
):
    _config, _sector_id, system_ids = seeded_sector
    system_id = system_ids[0]

    response = admin_client_with_wiki.post(
        f"/api/systems/{system_id}/wiki", json={"backend": "wikijs", "path": "systems/test-system"},
    )
    assert response.status_code == 201
    page = response.get_json()
    assert page["url"] == "http://wiki.invalid/wikijs/systems/test-system"

    call = fake_wiki_client.calls[0]
    assert call["backend"] == "wikijs"
    assert call["path"] == "systems/test-system"

    detail = admin_client_with_wiki.get(f"/api/systems/{system_id}").get_json()
    assert detail["wikijs_url"] == page["url"]
    assert detail["mediawiki_url"] is None
    # Markdown, not wikitext, is what wikijs got -- rendered fresh
    markdown = admin_client_with_wiki.get(f"/api/systems/{system_id}/text?format=markdown").get_json()
    assert call["content"] == markdown["content"]


def test_upload_system_wiki_mediawiki_uses_the_system_name_as_the_page_path(
    admin_client_with_wiki, seeded_sector, fake_wiki_client
):
    _config, _sector_id, system_ids = seeded_sector
    system_id = system_ids[0]
    name = admin_client_with_wiki.get(f"/api/systems/{system_id}").get_json()["name"]

    response = admin_client_with_wiki.post(f"/api/systems/{system_id}/wiki", json={"backend": "mediawiki"})
    assert response.status_code == 201

    call = fake_wiki_client.calls[0]
    assert call["backend"] == "mediawiki"
    assert call["path"] == name  # no path given -- mediawiki addresses by title/name, not a separate path

    detail = admin_client_with_wiki.get(f"/api/systems/{system_id}").get_json()
    assert detail["mediawiki_url"] is not None
    assert detail["wikijs_url"] is None
    wikitext = admin_client_with_wiki.get(f"/api/systems/{system_id}/text?format=wikitext").get_json()
    assert call["content"] == wikitext["content"]


def test_upload_system_wiki_page_exists_maps_to_409(admin_client_with_wiki, seeded_sector, monkeypatch):
    class RaisingClient:
        def __init__(self, backend, **kwargs):
            pass

        def create_page(self, **kwargs):
            raise WikiClientPageExistsError("a page already exists there")

    monkeypatch.setattr("api.routes.WikiClient", RaisingClient)
    _config, _sector_id, system_ids = seeded_sector

    response = admin_client_with_wiki.post(
        f"/api/systems/{system_ids[0]}/wiki", json={"backend": "wikijs", "path": "x"},
    )
    assert response.status_code == 409


def test_upload_system_wiki_returns_404_for_unknown_system(admin_client_with_wiki, fake_wiki_client):
    response = admin_client_with_wiki.post(
        "/api/systems/999999999/wiki", json={"backend": "wikijs", "path": "x"},
    )
    assert response.status_code == 404
    assert fake_wiki_client.calls == []


def test_upload_sector_wiki_persists_url_and_uses_generated_content(
    admin_client_with_wiki, seeded_sector, fake_wiki_client
):
    _config, sector_id, _system_ids = seeded_sector

    response = admin_client_with_wiki.post(f"/api/sectors/{sector_id}/wiki", json={"backend": "mediawiki"})
    assert response.status_code == 201
    page = response.get_json()

    call = fake_wiki_client.calls[0]
    assert call["backend"] == "mediawiki"
    assert "Test Sector" in call["content"]  # seeded_sector's own sector name

    detail = admin_client_with_wiki.get(f"/api/sectors/{sector_id}").get_json()
    assert detail["wiki_url"] == page["url"]


def test_update_sector_wiki_url_manually_sets_and_clears(admin_client, seeded_sector):
    """The admin "manually set the wiki link" affordance (`html/admin.py`)
    -- no upload, no wiki reachability required at all."""
    _config, sector_id, _system_ids = seeded_sector

    response = admin_client.patch(f"/api/sectors/{sector_id}", json={"wiki_url": "https://wiki.example.com/Sector"})
    assert response.status_code == 200
    assert admin_client.get(f"/api/sectors/{sector_id}").get_json()["wiki_url"] == "https://wiki.example.com/Sector"

    response = admin_client.patch(f"/api/sectors/{sector_id}", json={"wiki_url": None})
    assert response.status_code == 200
    assert admin_client.get(f"/api/sectors/{sector_id}").get_json()["wiki_url"] is None


def test_create_sector_rejects_wiki_url(admin_client):
    # A brand-new sector has never been uploaded anywhere -- wiki_url is
    # an update-only field (SECTOR_UPDATE_FIELDS), not accepted at create.
    response = admin_client.post(
        "/api/sectors", json={"name": "Test", "edge_ly": 10.0, "wiki_url": "https://wiki.example.com/x"}
    )
    assert response.status_code == 400
