# tests/test_api.py

"""
End-to-end tests for the read-only Flask API (`src/planetgen/api/`) against a real,
throwaway MySQL database (see `conftest.py`'s `mysql_config` fixture) --
every test here is skipped, not failed, when no MySQL test server is
configured/reachable.

Seeds a small, real sector (via `planetgen.db.store.save_sector`, the same
path `sectorGen.py` uses) rather than hand-building rows, so these tests
exercise the exact same read path (`planetgen.db.query`/`planetgen.db.store.load_sector`/
`load_star_system`) production traffic does.
"""

import math
import re
import time

import pytest

from planetgen.web.app import create_app
from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.api.config import Config
from planetgen.db import store as _db
from planetgen.admin import auth as adminAuth
from planetgen.galaxy import sector as spaceSector
from planetgen.db.store import MySQLConfig
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.viewport import tile_keys_containing, tiles_intersecting_sphere
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem
from wikiClient import WikiClientPageExistsError, WikiPage
from planetgen.db import query

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
def first_admin_password(mysql_config):
    """Seeds the control schema (`bootstrap_control_schema`) and returns
    the random first password it made -- there is no published default."""
    _username, password = adminAuth.bootstrap_control_schema(mysql_config)
    return password


@pytest.fixture
def default_admin_client(mysql_config, client, first_admin_password):
    """
    A `client` already logged in as the seeded default `admin`/`password`
    admin -- `must_change_credentials` is still set, so this is the
    fixture for tests confirming that gate actually blocks write/admin
    routes (see `authz.require_admin(fresh=True)`), not for exercising
    the writes themselves (use `admin_client` for that).

    Flask's test client keeps its own cookie jar across requests made on
    the same `client` instance, so the session cookie `POST /api/auth/login`
    sets is carried automatically into every later request this fixture's
    caller makes -- no manual cookie plumbing needed in tests.
    """
    # The write endpoints intentionally run with ensure_schema=False (see
    # routes._write_conn) -- a real deployment's content schema already
    # exists by the time an admin account is in use. This throwaway test
    # database starts completely empty, so lay the content schema down
    # here the same way any first `sectorGen.py`/`planetgen.cli.migrate` run would.
    _db.get_connection(mysql_config).close()
    response = client.post("/api/auth/login", json={
        "username": adminAuth.DEFAULT_ADMIN_USERNAME, "password": first_admin_password,
    })
    assert response.status_code == 200
    assert response.get_json()["must_change_credentials"] is True
    return client


@pytest.fixture
def admin_client(default_admin_client, first_admin_password):
    """A `client` logged in and past the forced credential change -- ready
    to exercise real write/admin endpoints against."""
    response = default_admin_client.post("/api/auth/change-credentials", json={
        "current_password": first_admin_password,
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
    _username, first_password = adminAuth.bootstrap_control_schema(mysql_config)
    _db.get_connection(mysql_config).close()
    response = client_with_wiki.post("/api/auth/login", json={
        "username": adminAuth.DEFAULT_ADMIN_USERNAME, "password": first_password,
    })
    assert response.status_code == 200
    response = client_with_wiki.post("/api/auth/change-credentials", json={
        "current_password": first_password,
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
    monkeypatch.setattr("planetgen.api.routes.WikiClient", _FakeWikiClient)
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
    Regression guard for a fixed bug: `query.open_readonly` raises a
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


@pytest.mark.parametrize("star_type", ["%", "_", "_2V", "%V", "G%"])
def test_systems_star_type_wildcards_match_only_themselves(client, seeded_sector, star_type):
    # `%`/`_` in ?star_type= are literal characters, not LIKE wildcards:
    # no star type contains them, so nothing matches.
    _config, sector_id, _system_ids = seeded_sector
    response = client.get("/api/systems", query_string={"sector_id": sector_id, "star_type": star_type})
    assert response.status_code == 200
    assert response.get_json()["total"] == 0


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
    # Generation is random: a system can still come out with no moons, belt
    # or comet despite the flags above, so draw until one has every kind.
    for _ in range(50):
        system = StarSystem(system_config=cfg)
        bodies = system.planets or []
        if (any(getattr(b, "moons", None) for b in bodies)
                and any(b.body_type == "a" for b in bodies) and system.comets):
            break
    system_id = _db.save_system(system, cfg, config=mysql_config)

    detail = client.get(f"/api/systems/{system_id}").get_json()
    sections = client.get(f"/api/systems/{system_id}/sections").get_json()
    assert any(p["moons"] for p in detail["planets"]) and detail["belts"] and detail["comets"]

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
    assert set(body["direct"]) == {"distance_ly", "bearing_deg", "mark_deg", "elevation_deg", "frame"}
    assert body["direct"]["frame"] == "sector"
    assert [leg["warp_factor"] for leg in body["warp_times"]] == [1, 2, 4, 8, 9, 9.5, 9.9, 9.995]
    assert [leg["fold_factor"] for leg in body["fold_times"]] == [4, 5, 6, 6.5, 7, 7.5, 8, 8.5]
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
    hop) and a `query._phenomenon_nav_key`-shaped string key (a
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

    sector_id = _db.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": 0.0, "center_y_pc": 0.0, "center_z_pc": 0.0,
        "galactic_radius_pc": 0.0,
    })

    from planetgen.generation.phenomena.nebula import Nebula
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


def test_the_control_database_is_never_listed_or_selectable(client, seeded_sector, monkeypatch):
    """The control schema (admin logins, sessions, API keys) shares the
    `planetgen` prefix by default (`planetgen_control`); it must not show
    up in `/api/databases` or be reachable with `?db=`."""
    mysql_config, _sector_id, _system_ids = seeded_sector
    monkeypatch.setenv(_db.CONTROL_DB_ENV_VAR, mysql_config.database)

    names = {item["name"] for item in client.get("/api/databases").get_json()["items"]}
    assert mysql_config.database not in names
    response = client.get(f"/api/sectors?db={mysql_config.database}")
    assert response.status_code == 404
    with pytest.raises(ValueError):
        _db.resolve_database(mysql_config, mysql_config.database)

    monkeypatch.delenv(_db.CONTROL_DB_ENV_VAR)
    real_load_config = _db.load_config
    monkeypatch.setattr(_db, "load_config", lambda: {**real_load_config(), "control_database": mysql_config.database})
    assert mysql_config.database not in {entry["name"] for entry in _db.list_databases(mysql_config)}


def test_list_databases_prefix_is_not_a_like_pattern(mysql_config):
    """`_` and `%` in the prefix match only themselves."""
    _db.get_connection(mysql_config).close()
    name = mysql_config.database
    assert name in {entry["name"] for entry in _db.list_databases(mysql_config, prefix=name)}
    # With `_` as a wildcard, name[:-1] + "_" would match `name` itself.
    assert _db.list_databases(mysql_config, prefix=name[:-1] + "_") == []
    assert _db.list_databases(mysql_config, prefix="%") == []
    assert _db.list_databases(mysql_config, prefix="planetgen%test") == []
    assert _db.escape_like("a_b%c\\d") == "a\\_b\\%c\\\\d"


def test_galaxy_sectors_excludes_unplaced_sectors(client, seeded_sector):
    # seeded_sector's own sector is never given a galaxy placement.
    response = client.get("/api/galaxy/sectors")
    assert response.status_code == 200
    assert response.get_json() == {"items": []}


def _place_sector(mysql_config, name, center_pc, edge_ly=10.0, address=None):
    sector = SpaceSector(name, edge_ly=edge_ly)
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=cfg), position=(0.0, 0.0, 0.0), system_config=cfg)
    position = {
        "center_x_pc": center_pc[0], "center_y_pc": center_pc[1], "center_z_pc": center_pc[2],
        "galactic_radius_pc": math.dist(center_pc, (0.0, 0.0, 0.0)),
    }
    if address is not None:
        position["ring_index"], position["layer_index"], position["ring_slot_index"] = address
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
    # Neither sector has a grid address, so neither is in a filled summary.
    assert tiles[near]["filled"] == {"g": 1, "cells": []}
    assert tiles[whole]["filled"]["cells"] == []
    assert tiles[near]["clouds"] == [] and tiles[whole]["clouds"] == []
    assert "density" not in body
    assert body["has_shape"] is False


def test_galaxy_filled_in_box_counts_every_placed_sector(mysql_config):
    """A tile's filled summary lists each addressed sector at g = 1 and
    past `max_cells` groups them into cells three sectors a side, with the
    cell's wedge by center angle -- nothing is sampled away."""
    from planetgen.galaxy.geometry import ring_sector_count, sector_position_pc

    addresses = [(0, 0, 0), (1, 0, 4), (2, 1, 0), (7, -1, 30), (7, 0, 30)]
    ids = []
    for n, address in enumerate(addresses):
        center = sector_position_pc(*address, 4.0)
        ids.append(_place_sector(mysql_config, f"Cell {n}", center, address=address))
    _place_sector(mysql_config, "No Address", (1.0, 1.0, 1.0))
    lo, hi = (-100.0, -100.0, -100.0), (100.0, 100.0, 100.0)
    conn = _db.get_connection(mysql_config)
    try:
        exact = query.galaxy_filled_in_box(conn, lo, hi, 200.0, 4.0)
        grouped = query.galaxy_filled_in_box(conn, lo, hi, 200.0, 4.0, max_cells=4)
        coarse = query.galaxy_filled_in_box(conn, lo, hi, 3 * 4.0 * query.GALAXY_TILE_FILLED_SCALE, 4.0)
    finally:
        conn.close()
    assert exact["g"] == 1
    assert [cell[:5] for cell in exact["cells"]] == [
        [ring, layer, slot, sector_id, 1] for (ring, layer, slot), sector_id in zip(addresses, ids)
    ]
    assert exact["cells"][0][5] == "Cell 0"
    assert ring_sector_count(7) > 30

    def expected_cells(g):
        cells = {}
        for ring, layer, slot in addresses:
            cell_ring = ring // g
            wedges = max(3, round(2 * math.pi * (cell_ring + 0.5)))
            angle = (slot + 0.5) / ring_sector_count(ring)
            key = (cell_ring, (layer + (g - 1) // 2) // g, int(angle * wedges))
            cells[key] = cells.get(key, 0) + 1
        return cells

    for summary in (grouped, coarse):
        assert summary["g"] == 3
        assert {tuple(cell[:3]): cell[3] for cell in summary["cells"]} == expected_cells(3)


def test_galaxy_clouds_in_box_lists_every_cloud_reaching_the_box(mysql_config):
    """A tile lists each nebula or supernova remnant whose sphere reaches
    into it, wherever its center is, largest first; past `max_clouds` the
    smallest are dropped. Point-like phenomena aren't clouds."""
    from planetgen.generation.phenomena.nebula import Nebula
    from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
    from stellarObjects.utils import pc_to_ly

    cfg = SystemConfig()
    sector_id = _place_sector(mysql_config, "Cloud Home", (0.0, 0.0, 0.0))
    placement = {"center_x_pc": 0.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 0.0}
    nebula = Nebula(cfg, name="Big Cloud")
    nebula.radius_ly = pc_to_ly(30.0)
    remnant = SupernovaRemnant(cfg, name="Small Shell")
    remnant.radius_ly = pc_to_ly(2.0)
    conn = _db.get_connection(mysql_config)
    try:
        nebula_id = _db.insert_nebula(conn, nebula, sector_id=sector_id, placement=placement)
        remnant_id = _db.insert_supernova_remnant(conn, remnant, sector_id=sector_id, placement=placement)
        conn.commit()
        home = query.galaxy_clouds_in_box(conn, (-10.0, -10.0, -10.0), (10.0, 10.0, 10.0))
        reached = query.galaxy_clouds_in_box(conn, (25.0, -10.0, -10.0), (45.0, 10.0, 10.0))
        corner = query.galaxy_clouds_in_box(conn, (25.0, 25.0, 25.0), (45.0, 45.0, 45.0))
        capped = query.galaxy_clouds_in_box(conn, (-10.0, -10.0, -10.0), (10.0, 10.0, 10.0), max_clouds=1)
    finally:
        conn.close()

    assert [(c["type"], c["id"]) for c in home] == [("nebula", nebula_id), ("supernova_remnant", remnant_id)]
    assert home[0]["radius_pc"] == pytest.approx(30.0)
    assert home[0]["name"] == "Big Cloud"
    assert home[0]["class"] and home[0]["descriptor"]
    # 25 pc from the center along x: the 30 pc nebula reaches, the shell doesn't.
    assert [c["id"] for c in reached] == [nebula_id]
    # The near corner is 25*sqrt(3) ~ 43 pc away: out of reach.
    assert corner == []
    assert [c["id"] for c in capped] == [nebula_id]


def _bright_row(position_pc, luminosity_sol, edge_pc):
    from planetgen.galaxy.geometry import sector_address_at

    ring, layer, slot = sector_address_at(position_pc, edge_pc)
    return (ring, layer, slot, *(int(round(v * 1000)) for v in position_pc), "young", "B", "III", 1e31, 1e7,
            15000.0, luminosity_sol * 3.82e26, 0.05, 0.1, 8.0, None, 1)


@pytest.mark.parametrize("outline, ranges_per_ring", [(False, False), (True, False), (True, True)])
def test_galaxy_bright_stars_in_box_matches_a_brute_force_search(mysql_config, monkeypatch, outline, ranges_per_ring):
    """Every box -- small (exact address ranges, or one range per ring),
    wrapped past +X, around the axis, or most of the galaxy (the
    luminosity index) -- lists exactly the most luminous stars inside it."""
    import random

    if ranges_per_ring:
        monkeypatch.setattr(query, "BRIGHT_STAR_MAX_EXACT_RANGES", 1)

    edge = 3.5
    rng = random.Random(7)
    stars = []
    for _ in range(1500):
        r, theta = 180.0 * math.sqrt(rng.random()), rng.uniform(0, 2 * math.pi)
        stars.append(((r * math.cos(theta), r * math.sin(theta), rng.uniform(-20.0, 20.0)), rng.uniform(500, 5e5)))
    if outline:
        rings = int(200 / edge) + 1
        _db.replace_galaxy_layers([(j, rings) for j in range(-7, 8)], config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        _db.insert_bright_stars(conn, (_bright_row(p, lum, edge) for p, lum in stars))
        conn.commit()
        ids = [row["id"] for row in conn.execute("SELECT id FROM bright_stars ORDER BY id").fetchall()]
        boxes = [((-300, -300, -300), (300, 300, 300)), ((-5, -5, -5), (5, 5, 5)), ((50, -20, -30), (90, 20, 30)),
                 ((-100, 40, -2), (-60, 80, 2)), ((0, 0, 0), (200, 200, 200)), ((100, -1, -30), (180, 0.5, 30))]
        for _ in range(30):
            size = rng.choice([4.0, 16.0, 64.0, 256.0])
            corner = tuple(rng.uniform(-200, 200 - size) for _ in range(3))
            boxes.append((corner, tuple(c + size for c in corner)))
        for lo, hi in boxes:
            for limit in (5, 400):
                got = query.galaxy_bright_stars_in_box(conn, lo, hi, edge, limit=limit)
                inside = [
                    (-lum, star_id) for star_id, (p, lum) in zip(ids, stars)
                    if all(lo[a] <= round(p[a] * 1000) / 1000 < hi[a] for a in range(3))
                ]
                assert [s["id"] for s in got] == [star_id for _l, star_id in sorted(inside)[:limit]], (lo, hi, limit)
        top = query.galaxy_bright_stars_in_box(conn, (-300, -300, -300), (300, 300, 300), edge, limit=1)[0]
    finally:
        conn.close()
    best = max(range(len(stars)), key=lambda i: stars[i][1])
    assert top["id"] == ids[best]
    assert top["luminosity_sol"] == pytest.approx(stars[best][1])
    assert (top["x"], top["y"], top["z"]) == pytest.approx(stars[best][0], abs=1e-3)
    assert top["temperature_k"] == 15000.0 and top["yerkes_class"] == "III" and top["system_id"] is None


def test_galaxy_bright_stars_in_box_takes_big_boxes_from_the_galaxy_wide_sample(mysql_config, monkeypatch):
    """A box over the row budget (MAP.47) is answered from the galaxy-wide
    sample: exactly the box's brightest when the sample holds `limit` of
    them there, else the sample's stars inside it, still brightest first."""
    import random

    edge = 3.5
    rng = random.Random(11)
    stars = []
    for _ in range(600):
        r, theta = 150.0 * math.sqrt(rng.random()), rng.uniform(0, 2 * math.pi)
        stars.append(((r * math.cos(theta), r * math.sin(theta), rng.uniform(-10.0, 10.0)), rng.uniform(500, 5e5)))
    monkeypatch.setattr(query, "GALAXY_TILE_BRIGHT_STAR_ROW_BUDGET", 0)
    conn = _db.get_connection(mysql_config)
    try:
        _db.insert_bright_stars(conn, (_bright_row(p, lum, edge) for p, lum in stars))
        conn.commit()
        ids = [row["id"] for row in conn.execute("SELECT id FROM bright_stars ORDER BY id").fetchall()]
        ranked = sorted(zip(ids, stars), key=lambda pair: -pair[1][1])
        reads = []

        def sample_of(count):
            def brightest():
                reads.append(count)
                return query.galaxy_brightest_stars(conn, count)
            return brightest

        assert [s["id"] for s in query.galaxy_brightest_stars(conn, 5)] == [i for i, _s in ranked[:5]]
        for lo, hi in [((-200, -200, -50), (200, 200, 50)), ((0, 0, -50), (200, 200, 50)), ((20, -60, -5), (90, 10, 5))]:
            inside = [i for i, (p, _lum) in ranked if all(lo[a] <= round(p[a] * 1000) / 1000 < hi[a] for a in range(3))]
            whole = query.galaxy_bright_stars_in_box(conn, lo, hi, edge, limit=20, brightest=sample_of(len(stars)))
            assert [s["id"] for s in whole] == inside[:20]
            top = {i for i, _s in ranked[:100]}
            thinned = query.galaxy_bright_stars_in_box(conn, lo, hi, edge, limit=400, brightest=sample_of(100))
            assert [s["id"] for s in thinned] == [i for i in inside if i in top]
        assert reads
    finally:
        conn.close()


def test_galaxy_bright_stars_in_box_is_empty_without_a_scatter(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        assert query.galaxy_bright_stars_in_box(conn, (-10, -10, -10), (10, 10, 10), 3.5) == []
    finally:
        conn.close()


def _generated_sectors(mysql_config, luminosities):
    """One grid-placed sector per entry of `luminosities` along +X, 4 pc
    apart, each holding one single star of that luminosity (solar)."""
    from planetgen.physics import constants

    sector_ids = [
        _place_sector(mysql_config, f"Lit {i}", (2.0 + 4.0 * i, 2.0, 0.0), address=(i, 0, 0))
        for i in range(len(luminosities))
    ]
    conn = _db.get_connection(mysql_config)
    try:
        for sector_id, lum in zip(sector_ids, luminosities):
            conn.execute(
                "UPDATE stars SET luminosity_w = ? WHERE star_system_id IN "
                "(SELECT id FROM star_systems WHERE sector_id = ?)",
                (lum * constants.SOLAR_LUMINOSITY, sector_id),
            )
        conn.commit()
        star_ids = [
            conn.execute(
                "SELECT st.id FROM stars st JOIN star_systems ss ON ss.id = st.star_system_id WHERE ss.sector_id = ?",
                (sector_id,),
            ).fetchone()["id"]
            for sector_id in sector_ids
        ]
    finally:
        conn.close()
    return sector_ids, star_ids


def test_galaxy_generated_stars_in_box_lists_the_brightest_above_the_floor(mysql_config, monkeypatch):
    """MAP.51: a tile lists its generated systems' stars at or above its
    floor, most luminous first, placed at sector center plus system
    offset; a star also pre-placed as a bright star isn't listed twice;
    past the sector budget only every k-th sector is read."""
    import zlib

    lums = [0.002, 0.3, 40.0, 1.0, 0.01]
    sector_ids, star_ids = _generated_sectors(mysql_config, lums)
    lo, hi = (0.0, 0.0, -5.0), (64.0, 64.0, 5.0)
    conn = _db.get_connection(mysql_config)
    try:
        every = query.galaxy_generated_stars_in_box(conn, lo, hi, 0.0, len(lums))
        assert [s["id"] for s in every] == [star_ids[i] for i in (2, 3, 1, 4, 0)]
        top = every[0]
        assert top["luminosity_sol"] == pytest.approx(40.0, rel=1e-3)
        assert (top["x"], top["y"], top["z"]) == pytest.approx((10.0, 2.0, 0.0), abs=0.01)
        assert (top["ring_index"], top["layer_index"], top["ring_slot_index"]) == (2, 0, 0)
        assert top["radius_sol"] > 0 and top["temperature_k"] > 0 and top["system_id"] and top["name"]
        floored = query.galaxy_generated_stars_in_box(conn, lo, hi, 0.25, len(lums))
        assert [s["id"] for s in floored] == [star_ids[i] for i in (2, 3, 1)]
        assert len(query.galaxy_generated_stars_in_box(conn, lo, hi, 0.0, len(lums), limit=2)) == 2
        assert query.galaxy_generated_stars_in_box(conn, (100.0, 0.0, -5.0), (164.0, 64.0, 5.0), 0.0, 1) == []
        assert query.galaxy_generated_stars_in_box(conn, lo, hi, 0.0, 0) == []

        monkeypatch.setattr(query, "GALAXY_TILE_GENERATED_STAR_SECTOR_BUDGET", 2)
        sampled = query.galaxy_generated_stars_in_box(conn, lo, hi, 0.0, len(lums))
        assert {s["id"] for s in sampled} == {
            star for sector, star in zip(sector_ids, star_ids) if zlib.crc32(str(sector).encode()) % 3 == 0
        }
        monkeypatch.undo()

        system_id = top["system_id"]
        _db.insert_bright_stars(conn, [_bright_row((10.0, 2.0, 0.0), 600.0, 4.0)])
        conn.execute("UPDATE bright_stars SET star_system_id = ?", (system_id,))
        conn.commit()
        assert star_ids[2] not in [s["id"] for s in query.galaxy_generated_stars_in_box(conn, lo, hi, 0.0, 5)]
    finally:
        conn.close()


def test_generated_star_floor_drops_fourfold_per_finer_tile_level():
    from planetgen.galaxy.viewport import TILE_MAX_LEVEL

    assert query.generated_star_floor_sol(TILE_MAX_LEVEL) == 0.0
    assert query.generated_star_floor_sol(TILE_MAX_LEVEL - 1) == pytest.approx(query.GENERATED_STAR_FLOOR_SOL_AT_32_PC)
    for level in range(3, TILE_MAX_LEVEL - 1):
        assert query.generated_star_floor_sol(level) == pytest.approx(4 * query.generated_star_floor_sol(level + 1))
    assert query.generated_star_floor_sol(0) is None


def test_galaxy_tiles_lists_generated_stars_by_tile_level(client, mysql_config):
    """The finest tile lists every generated star, a coarser one only the
    brighter ones, and a galaxy-sized one none (the bright stars cover it)."""
    _sector_ids, star_ids = _generated_sectors(mysql_config, [0.001, 5.0])
    finest = tiles_intersecting_sphere(12, (2.0, 2.0, 0.0), 0.0)[0]
    mid = tiles_intersecting_sphere(6, (2.0, 2.0, 0.0), 0.0)[0]
    tiles = client.get(f"/api/galaxy/tiles?tiles={finest},{mid},0/0/0/0").get_json()["tiles"]
    assert [s["id"] for s in tiles[finest]["generated"]] == [star_ids[1], star_ids[0]]
    assert [s["id"] for s in tiles[mid]["generated"]] == [star_ids[1]]
    assert tiles["0/0/0/0"]["generated"] == []


def test_galaxy_tiles_list_every_star_and_point_phenomenon_at_sector_zoom(client, mysql_config, monkeypatch):
    """MAP.80: a finest tile lists every generated star up to its own,
    larger cap; tiles from `POINT_PHENOMENON_MIN_LEVEL` in list the black
    holes, neutron stars and quasars centered in them, brightest first,
    and coarser tiles none."""
    import math

    from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.phenomena.quasar import Quasar

    _sector_ids, star_ids = _generated_sectors(mysql_config, [0.001, 5.0, 0.3])
    sector_id = _sector_ids[0]

    def placement(x):
        return {"center_x_pc": x, "center_y_pc": 2.0, "center_z_pc": 0.0, "galactic_radius_pc": math.hypot(x, 2.0)}

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            hole = _db.insert_black_hole(conn, BlackHole(SystemConfig()), sector_id=sector_id, placement=placement(1.0))
            pulsar = _db.insert_neutron_star(conn, NeutronStar(SystemConfig()), sector_id=sector_id,
                                             placement=placement(3.0))
            quasar = _db.insert_quasar(conn, Quasar(SystemConfig()), sector_id=sector_id, placement=placement(2.5))
            # Bound to a system, so not placed on its own.
            _db.insert_black_hole(conn, BlackHole(SystemConfig()))
    finally:
        conn.close()

    finest = tiles_intersecting_sphere(12, (2.0, 2.0, 0.0), 0.0)[0]
    sector_level = tiles_intersecting_sphere(query.POINT_PHENOMENON_MIN_LEVEL, (2.0, 2.0, 0.0), 0.0)[0]
    coarse = tiles_intersecting_sphere(query.POINT_PHENOMENON_MIN_LEVEL - 1, (2.0, 2.0, 0.0), 0.0)[0]
    monkeypatch.setattr(query, "GALAXY_TILE_MAX_GENERATED_STARS", 1)
    tiles = client.get(f"/api/galaxy/tiles?tiles={finest},{sector_level},{coarse}").get_json()["tiles"]
    assert [s["id"] for s in tiles[finest]["generated"]] == [star_ids[1], star_ids[2], star_ids[0]]
    assert len(tiles[sector_level]["generated"]) == 1
    points = tiles[finest]["points"]
    assert {(p["type"], p["id"]) for p in points} == {("black_hole", hole), ("neutron_star", pulsar), ("quasar", quasar)}
    assert [p["luminosity_sol"] for p in points] == sorted((p["luminosity_sol"] for p in points), reverse=True)
    assert points[0]["type"] == "quasar"
    by_type = {p["type"]: p for p in points}
    assert (by_type["neutron_star"]["x"], by_type["neutron_star"]["y"]) == pytest.approx((3.0, 2.0))
    assert by_type["neutron_star"]["descriptor"] in ("young", "millisecond", "non-pulsing")
    assert by_type["black_hole"]["name"]
    assert tiles[sector_level]["points"] == points
    assert tiles[coarse]["points"] == []


def test_galaxy_shape_reports_the_bright_star_scatter(client, mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        body = client.get("/api/galaxy/shape").get_json()
        assert body["shape"] is None
        assert body["bright_stars"] == {"scattered": False, "min_luminosity_sol": None, "seed": None,
                                        "default_min_luminosity_sol": 1000.0}
        with conn:
            conn.execute("INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc,"
                         " bulge_scale_radius_pc, bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude,"
                         " spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm, edge_pc,"
                         " expected_system_count_at_density_1, outer_ring_index)"
                         " VALUES (1, 1, 1, 1, 1, 2, 0.2, 0.3, 1, 0, 1, 4, 10, 5)")
            _db.record_bright_star_scatter(conn, 500.0, 1234)
    finally:
        conn.close()
    body = client.get("/api/galaxy/shape").get_json()
    assert body["shape"]["edge_pc"] == 4
    assert body["bright_stars"] == {"scattered": True, "min_luminosity_sol": 500.0, "seed": 1234,
                                    "default_min_luminosity_sol": 1000.0}


def test_galaxy_bright_stars_in_one_cell(client, mysql_config):
    """`GET /api/galaxy/bright-stars`: one cell's waiting bright stars,
    brightest first; `all=1` adds the filled ones; bad addresses are 400."""
    assert client.get("/api/galaxy/bright-stars?ring=3&layer=0&slot=1").get_json() == {"items": []}
    row = (3, 0, 1, 12000, 3000, 0, "young", "B2V", "V", 1.4e31, 3.0e6, 22000.0)
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_bright_stars(conn, [row + (800 * 3.828e26, 0.02, 0.03, 7.0, 0.03, 42),
                                           row + (600 * 3.828e26, 0.02, 0.03, 7.0, 0.03, 42)])
    finally:
        conn.close()
    items = client.get("/api/galaxy/bright-stars?ring=3&layer=0&slot=1").get_json()["items"]
    assert [item["luminosity_sol"] for item in items] == pytest.approx([800, 600], rel=1e-2)
    assert client.get("/api/galaxy/bright-stars?ring=3&layer=0&slot=2").get_json() == {"items": []}
    assert len(client.get("/api/galaxy/bright-stars?ring=3&layer=0&slot=1&all=1").get_json()["items"]) == 2
    assert client.get("/api/galaxy/bright-stars?ring=x&layer=0&slot=1").status_code == 400
    assert client.get("/api/galaxy/bright-stars?ring=3").status_code == 400


def test_galaxy_stage_counts_generated_sectors_down_the_ladder(client, mysql_config):
    """`/api/galaxy/stage` counts each child block's generated sectors at
    every level, and at a level-3 block lists the sectors themselves;
    `galaxy_changes` names the stages a new sector makes stale."""
    from planetgen.galaxy.drill import drill_block_sectors, drill_chain_of, format_drill_key
    from planetgen.galaxy.geometry import sector_position_pc

    home = (1705, -20, 3225)
    chain = drill_chain_of(*home)
    level3 = chain[2]
    # Two more sectors in the same level-3 block, one on another layer.
    others = [s for s in drill_block_sectors(level3, -20) if (s.ring, s.slab, s.wedge) != home][:1]
    others += drill_block_sectors(level3, -21)[:1]
    addresses = [home] + [(s.ring, s.slab, s.wedge) for s in others] + [(0, 0, 0)]
    before = query.galaxy_changes(_db.get_connection(mysql_config))["state"]
    ids = [_place_sector(mysql_config, f"Stage {n}", sector_position_pc(*a, 4.0), address=a) for n, a in enumerate(addresses)]

    galaxy = client.get("/api/galaxy/stage").get_json()
    assert galaxy["at"] is None and galaxy["child_m"] == 243
    top = {(c["ring"], c["wedge"], c["slab"]): c["generated"] for c in galaxy["children"]}
    assert top[(chain[0].ring, chain[0].wedge, chain[0].slab)] == 3
    core = drill_chain_of(0, 0, 0)[0]
    assert top[(core.ring, core.wedge, core.slab)] == 1

    def counted(children):
        return [{k: v for k, v in child.items() if k != "look"} for child in children]

    stage = client.get(f"/api/galaxy/stage?at={format_drill_key(chain[0])}").get_json()
    assert stage["child_m"] == 27
    assert counted(stage["children"]) == [
        {"ring": chain[1].ring, "wedge": chain[1].wedge, "slab": chain[1].slab, "generated": 3}]

    stage = client.get(f"/api/galaxy/stage?at={format_drill_key(chain[1])}").get_json()
    assert counted(stage["children"]) == [
        {"ring": level3.ring, "wedge": level3.wedge, "slab": level3.slab, "generated": 3}]

    stage = client.get(f"/api/galaxy/stage?at={format_drill_key(level3)}").get_json()
    assert stage["child_m"] == 1
    assert sorted((s["ring"], s["layer"], s["slot"]) for s in stage["sectors"]) == sorted(addresses[:3])
    assert {s["id"] for s in stage["sectors"]} == set(ids[:3])
    assert all(s["system_count"] == 1 for s in stage["sectors"])
    assert len(stage["children"]) == 3 and all(c["generated"] == 1 for c in stage["children"])

    # A block next door holds none of them.
    empty = client.get(f"/api/galaxy/stage?at=3.{level3.ring}.{level3.wedge + 1}.{level3.slab}").get_json()
    assert empty["children"] == [] and empty["sectors"] == []
    for bad in ("81.0.0.0", "243.0.99.0", "nonsense", "1.0.0.0"):
        assert client.get(f"/api/galaxy/stage?at={bad}").status_code == 400

    changes = query.galaxy_changes(_db.get_connection(mysql_config), before)
    assert not changes["full"]
    assert set(changes["stages"]) == {"galaxy"} | {format_drill_key(b) for b in chain[:3]} | {
        format_drill_key(b) for b in drill_chain_of(0, 0, 0)[:3]}


def test_galaxy_stage_looks_average_their_sectors_stats(client, mysql_config):
    """MAP.86: each stage child, and each sector at a level-3 block, carries
    a `look` from `sector_stats`: the mean fill share of its generated
    sectors, the mean color of those with stars, and how many had one."""
    from planetgen.galaxy.drill import drill_block_sectors, drill_chain_of, format_drill_key
    from planetgen.galaxy.geometry import sector_position_pc

    home = (1705, -20, 3225)
    chain = drill_chain_of(*home)
    others = [s for s in drill_block_sectors(chain[2], -20) if (s.ring, s.slab, s.wedge) != home][:2]
    addresses = [home] + [(s.ring, s.slab, s.wedge) for s in others]
    for n, address in enumerate(addresses):
        _place_sector(mysql_config, f"Look {n}", sector_position_pc(*address, 4.0), address=address)
    looks = [(0.2, (1.0, 0.5, 0.0)), (0.6, (0.0, 0.5, 1.0)), (0.4, None)]
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            for address, (share, color) in zip(addresses, looks):
                conn.execute(
                    "INSERT INTO sector_stats (ring_index, layer_index, ring_slot_index, bright_level_sol, fill_share,"
                    " color_r, color_g, color_b) VALUES (?, ?, ?, 0, ?, ?, ?, ?) ON DUPLICATE KEY UPDATE"
                    " fill_share = ?, color_r = ?, color_g = ?, color_b = ?",
                    (*address, share, *(color or (None,) * 3), share, *(color or (None,) * 3)),
                )
    finally:
        conn.close()

    stage = client.get(f"/api/galaxy/stage?at={format_drill_key(chain[1])}").get_json()
    [child] = stage["children"]
    assert child["generated"] == 3
    assert child["look"] == {"share": pytest.approx(0.4), "color": pytest.approx([0.5, 0.5, 0.5]), "colored": 2}
    galaxy = client.get("/api/galaxy/stage").get_json()
    [top] = [c for c in galaxy["children"] if (c["ring"], c["wedge"], c["slab"]) == (chain[0].ring, chain[0].wedge,
                                                                                     chain[0].slab)]
    assert top["look"] == child["look"]

    stage = client.get(f"/api/galaxy/stage?at={format_drill_key(chain[2])}").get_json()
    by_address = {(s["ring"], s["layer"], s["slot"]): s["look"] for s in stage["sectors"]}
    assert by_address[home] == {"share": pytest.approx(0.2), "color": pytest.approx([1.0, 0.5, 0.0]), "colored": 1}
    assert by_address[addresses[2]] == {"share": pytest.approx(0.4), "color": None, "colored": 0}
    assert {c["look"]["colored"] for c in stage["children"]} == {0, 1}


def test_galaxy_locate_finds_sectors_and_systems_by_name(client, mysql_config):
    """`/api/galaxy/locate` finds sectors and systems whose name contains
    the term, each with its sector address, exact matches first; sectors
    with no address are left out."""
    from planetgen.galaxy.geometry import sector_position_pc

    address = (12, 1, 30)
    sector_id = _place_sector(mysql_config, "Belcana", sector_position_pc(*address, 4.0), address=address)
    other = (13, 0, 40)
    _place_sector(mysql_config, "Belcana Reach", sector_position_pc(*other, 4.0), address=other)
    # No address: never offered, since the map cannot fly to it.
    _place_sector(mysql_config, "Belcana Lost", (50.0, 0.0, 0.0))

    conn = _db.get_connection(mysql_config)
    try:
        conn.execute("UPDATE star_systems SET name = 'Belcana' WHERE sector_id = ?", (sector_id,))
        conn.commit()
    finally:
        conn.close()

    matches = client.get("/api/galaxy/locate?q=belcana").get_json()["matches"]
    assert [(m["kind"], m["name"]) for m in matches[:2]] == [("sector", "Belcana"), ("system", "Belcana")]
    assert all(m["name"] != "Belcana Lost" for m in matches)
    assert {m["name"] for m in matches} == {"Belcana", "Belcana Reach"}
    first = matches[0]
    assert (first["ring"], first["layer"], first["slot"]) == address
    assert first["sector_id"] == sector_id
    system = matches[1]
    assert (system["ring"], system["layer"], system["slot"]) == address
    assert system["sector_name"] == "Belcana"

    assert client.get("/api/galaxy/locate?q=  ").get_json()["matches"] == []
    assert client.get("/api/galaxy/locate?q=nothing-like-this").get_json()["matches"] == []
    # A term with LIKE wildcards in it is a plain substring, not a pattern.
    assert client.get("/api/galaxy/locate?q=%25").get_json()["matches"] == []


def test_galaxy_sectors_in_box_samples_evenly_past_the_cap(mysql_config):
    """A box holding more than `limit` placed sectors returns every Nth by
    id, not the lowest ids -- a neighborhood is generated outward from the core,
    so its lowest ids are only its core-facing half."""
    ids = [_place_sector(mysql_config, f"Row {n}", (float(n), 0.0, 0.0)) for n in range(10)]
    conn = _db.get_connection(mysql_config)
    try:
        sampled = query.galaxy_sectors_in_box(conn, (-1.0, -1.0, -1.0), (20.0, 1.0, 1.0), limit=4)
        everything = query.galaxy_sectors_in_box(conn, (-1.0, -1.0, -1.0), (20.0, 1.0, 1.0), limit=10)
    finally:
        conn.close()
    # ceil(10 / 4) = 3: rows 0, 3, 6, 9 -- spanning the whole run.
    assert [s["id"] for s in sampled] == [ids[0], ids[3], ids[6], ids[9]]
    assert [s["id"] for s in everything] == ids


def test_galaxy_tiles_rejects_bad_requests(client, mysql_config):
    _db.get_connection(mysql_config).close()
    assert client.get("/api/galaxy/tiles?tiles=13/0/0/0").status_code == 400
    assert client.get("/api/galaxy/tiles?tiles=1/2/0/0").status_code == 400
    assert client.get("/api/galaxy/tiles?tiles=nonsense").status_code == 400
    too_many = ",".join(f"12/{i}/0/0" for i in range(129))
    assert client.get(f"/api/galaxy/tiles?tiles={too_many}").status_code == 400
    empty = client.get("/api/galaxy/tiles?tiles=")
    assert empty.status_code == 200
    assert empty.get_json()["tiles"] == {}


def test_galaxy_stamp_changes_when_a_sector_is_placed(client, mysql_config):
    _db.get_connection(mysql_config).close()
    first = client.get("/api/galaxy/stamp").get_json()["stamp"]
    assert re.fullmatch(r"[0-9a-f]{16}", first)
    assert client.get("/api/galaxy/stamp").get_json()["stamp"] == first
    _place_sector(mysql_config, "New Sector", (1.0, 2.0, 3.0))
    assert client.get("/api/galaxy/stamp").get_json()["stamp"] != first

def _galaxy_changes(client, since):
    response = client.get("/api/galaxy/changes", query_string={"since": since})
    assert response.status_code == 200
    return response.get_json()


def test_galaxy_changes_lists_only_the_edited_sectors_tiles(client, mysql_config):
    # A rename bumps `sectors.modified_at` (v27), which names the one tile
    # per level holding that sector -- not the other sector's.
    renamed_id = _place_sector(mysql_config, "Renamed", (5.0, 5.0, 5.0))
    _place_sector(mysql_config, "Untouched", (9000.0, -3000.0, 20.0))
    before = client.get("/api/galaxy/stamp").get_json()

    unchanged = _galaxy_changes(client, before["state"])
    assert unchanged == {"stamp": before["stamp"], "state": before["state"], "full": False, "tiles": [],
                         "stages": []}

    time.sleep(0.01)
    conn = _db.get_connection(mysql_config)
    try:
        conn.execute("UPDATE sectors SET name = ? WHERE id = ?", ("Renamed Again", renamed_id))
        conn.commit()
    finally:
        conn.close()

    changes = _galaxy_changes(client, before["state"])
    assert changes["full"] is False
    assert changes["stamp"] != before["stamp"]
    assert changes["tiles"] == sorted(tile_keys_containing((5.0, 5.0, 5.0)))
    assert _galaxy_changes(client, changes["state"])["tiles"] == []


def test_galaxy_changes_lists_a_new_sectors_tiles(client, mysql_config):
    _place_sector(mysql_config, "Old", (5.0, 5.0, 5.0))
    before = client.get("/api/galaxy/stamp").get_json()
    _place_sector(mysql_config, "New", (-800.0, 40.0, 3.0))
    changes = _galaxy_changes(client, before["state"])
    assert changes["full"] is False
    assert changes["tiles"] == sorted(tile_keys_containing((-800.0, 40.0, 3.0)))


def test_galaxy_changes_lists_the_tiles_of_a_sector_that_lost_a_system(admin_client, mysql_config):
    # Tiles show each sector's system count; deleting a system bumps its
    # sector (`_db.touch_sector`).
    sector_id = _place_sector(mysql_config, "Shrinking", (100.0, 200.0, 300.0))
    conn = _db.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
    finally:
        conn.close()
    before = admin_client.get("/api/galaxy/stamp").get_json()
    time.sleep(0.01)
    assert admin_client.delete(f"/api/systems/{system_id}").status_code == 200
    changes = _galaxy_changes(admin_client, before["state"])
    assert changes["full"] is False
    assert changes["tiles"] == sorted(tile_keys_containing((100.0, 200.0, 300.0)))


def test_galaxy_changes_is_full_after_a_deletion_or_a_bad_since(admin_client, mysql_config):
    # A deleted sector leaves no row to locate its tiles by.
    sector_id = _place_sector(mysql_config, "Doomed", (5.0, 5.0, 5.0))
    _place_sector(mysql_config, "Survivor", (50.0, 5.0, 5.0))
    before = admin_client.get("/api/galaxy/stamp").get_json()
    assert admin_client.delete(f"/api/sectors/{sector_id}").status_code == 200
    changes = _galaxy_changes(admin_client, before["state"])
    assert changes["full"] is True and changes["tiles"] == []

    for since in ("", "nonsense", "0123456789abcdef.1.2.3.4"):
        assert _galaxy_changes(admin_client, since)["full"] is True


def test_search_returns_facets_and_matches_a_class_tag(client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector

    response = client.get("/api/search")
    assert response.status_code == 200
    body = response.get_json()
    assert set(body["facets"]) == {
        "type", "spectral", "luminosity", "class", "body", "life",
        "moon_class", "moon_body", "moon_life", "density", "phenomenon", "phenomenon_class",
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


def test_search_pages_each_result_panel(client, mysql_config):
    for i in range(60):
        _db.save_sector(SpaceSector(f"Pager Sector {i:03d}", edge_ly=10.0), config=mysql_config)

    first = client.get("/api/search?sector_q=Pager&limit=50").get_json()["results"]["sectors"]
    assert (first["total"], first["limit"], first["offset"], first["truncated"]) == (60, 50, 0, True)
    assert [row["name"] for row in first["rows"]][:2] == ["Pager Sector 000", "Pager Sector 001"]
    assert len(first["rows"]) == 50

    second = client.get("/api/search?sector_q=Pager&limit=50&sectors_offset=50").get_json()["results"]["sectors"]
    assert [row["name"] for row in second["rows"]] == [f"Pager Sector {i:03d}" for i in range(50, 60)]
    assert (second["offset"], second["truncated"]) == (50, True)

    # Past the end comes back as the last page, with its real offset.
    past_end = client.get("/api/search?sector_q=Pager&limit=50&sectors_offset=900").get_json()
    assert past_end["results"]["sectors"]["offset"] == 50

    # Without a limit, the old 300-row default still fits every match.
    default = client.get("/api/search?sector_q=Pager").get_json()["results"]["sectors"]
    assert (default["limit"], len(default["rows"]), default["truncated"]) == (300, 60, False)

    assert client.get("/api/search?sector_q=Pager&sectors_offset=-1").status_code == 400
    assert client.get("/api/search?sector_q=Pager&limit=0").status_code == 400


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


def test_login_sets_cookie_scoped_to_root_path_not_api(first_admin_password, client):
    """
    Regression test: the session cookie's `Path` attribute must be `/`,
    not `/api`. The admin pages (`/admin`, `/account`, ...) live at the
    site root -- per RFC 6265 path matching, a cookie scoped to `/api` is
    never attached by the browser to a request for `/admin`, so a
    `Path=/api` cookie would lock a user out of every admin page
    immediately after a successful login. This asserts the cookie's scope
    directly.
    """
    response = client.post("/api/auth/login", json={
        "username": adminAuth.DEFAULT_ADMIN_USERNAME, "password": first_admin_password,
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


def test_change_credentials_rejects_weak_password(default_admin_client, first_admin_password):
    response = default_admin_client.post("/api/auth/change-credentials", json={
        "current_password": first_admin_password, "new_username": "someone", "new_password": "short",
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


def test_seeded_admin_has_no_published_default_password(first_admin_password, client):
    """Security #39: `admin`/`password` no longer logs in on a fresh
    install; only the random first password does."""
    response = client.post("/api/auth/login", json={"username": "admin", "password": "password"})
    assert response.status_code == 401
    response = client.post("/api/auth/login", json={"username": "admin", "password": first_admin_password})
    assert response.status_code == 200
    assert response.get_json()["must_change_credentials"] is True


def test_change_credentials_logs_out_other_sessions_but_not_the_caller(
        mysql_config, first_admin_password, default_admin_client):
    """Security #45: after a credential change, another browser's session
    is gone, the browser that made the change gets a fresh session and
    stays logged in, and the old session cookie no longer works."""
    other = default_admin_client.application.test_client()
    assert other.post("/api/auth/login", json={
        "username": "admin", "password": first_admin_password}).status_code == 200
    assert other.get("/api/auth/me").status_code == 200
    old_cookie = default_admin_client.get_cookie(SESSION_COOKIE_NAME).value

    response = default_admin_client.post("/api/auth/change-credentials", json={
        "current_password": first_admin_password, "new_username": "renamed", "new_password": TEST_ADMIN_PASSWORD,
    })
    assert response.status_code == 200
    assert default_admin_client.get_cookie(SESSION_COOKIE_NAME).value != old_cookie
    assert default_admin_client.get("/api/auth/me").get_json()["username"] == "renamed"
    assert other.get("/api/auth/me").status_code == 401
    stale = default_admin_client.application.test_client()
    stale.set_cookie(SESSION_COOKIE_NAME, old_cookie, domain="localhost")
    assert stale.get("/api/auth/me").status_code == 401


def test_change_credentials_keeps_api_keys(admin_client):
    """Security #45: API keys survive a credential change."""
    key = admin_client.post("/api/auth/api-keys", json={"label": "ci"}).get_json()["key"]
    response = admin_client.post("/api/auth/change-credentials", json={
        "current_password": TEST_ADMIN_PASSWORD, "new_username": "test-admin",
        "new_password": TEST_ADMIN_PASSWORD + "-2",
    })
    assert response.status_code == 200
    bearer = admin_client.application.test_client()
    assert bearer.get("/api/auth/me", headers={"Authorization": f"Bearer {key}"}).status_code == 200


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
    response = admin_client.post("/api/systems", json={"sector_id": "1"})
    assert response.status_code == 400

    response = admin_client.post("/api/systems", json={"position": [0, 0, 0]})
    assert response.status_code == 400  # a position needs a sector

    response = admin_client.post("/api/systems", json={"sector_id": 999999})
    assert response.status_code == 404

    response = admin_client.post("/api/systems", json={"age": "ancient"})
    assert response.status_code == 400

    response = admin_client.post("/api/systems", json={"num_orbits": -1})
    assert response.status_code == 400


def test_create_system_inside_a_sector_keeps_clear_of_its_systems(seeded_sector, admin_client):
    config, sector_id, system_ids = seeded_sector
    response = admin_client.post("/api/systems", json={"sector_id": sector_id, "star_type": "K2V", "planets": False,
                                                        "binary_system": False})
    assert response.status_code == 201, response.get_json()
    body = response.get_json()
    assert body["sector_id"] == sector_id
    placed = body["position"]
    assert all(abs(c) <= 5.0 for c in placed)

    detail = admin_client.get(f"/api/systems/{body['id']}").get_json()
    assert detail["sector_id"] == sector_id
    assert detail["location"]

    conn = _db.get_connection(config)
    try:
        sector = _db.sector_for_placement(conn, sector_id)
    finally:
        conn.close()
    new = next(e for e in sector.entries if e.star_system.name == detail["name"])
    for other in sector.entries:
        if other is not new:
            assert spaceSector.distance_between(new.position, other.position) >= \
                spaceSector.required_separation_ly(new.star_system, other.star_system) - 1e-3


def test_create_system_at_an_explicit_position_inside_the_sector(seeded_sector, admin_client):
    _config, sector_id, _ids = seeded_sector
    response = admin_client.post("/api/systems", json={"sector_id": sector_id, "position": [4.0, -4.0, 2.5],
                                                        "planets": False})
    assert response.status_code == 201
    assert response.get_json()["position"] == pytest.approx([4.0, -4.0, 2.5], abs=1e-3)

    response = admin_client.post("/api/systems", json={"sector_id": sector_id, "position": [40.0, 0, 0]})
    assert response.status_code == 400  # outside the sector


def test_regenerate_keeps_the_system_but_replaces_its_bodies(seeded_sector, admin_client):
    config, sector_id, system_ids = seeded_sector
    system_id = system_ids[0]
    before = admin_client.get(f"/api/systems/{system_id}").get_json()
    old_star_ids = {star["id"] for star in before["stars"]}

    response = admin_client.patch(f"/api/systems/{system_id}", json={
        "regenerate": {"star_type": "K1V", "planets": True, "binary_system": False}})
    assert response.status_code == 200, response.get_json()
    assert response.get_json()["regenerated"] is True

    after = admin_client.get(f"/api/systems/{system_id}").get_json()
    assert (after["name"], after["sector_id"], after["location"]) == (before["name"], before["sector_id"],
                                                                       before["location"])
    assert after["stars"][0]["star_type"].startswith("K1V")
    assert after["stars"][0]["name"] == before["name"]
    assert not old_star_ids & {star["id"] for star in after["stars"]}
    for planet in after["planets"]:
        assert planet["name"].startswith(before["name"])

    conn = _db.get_connection(config)
    try:
        # The stand-in row used while swapping is gone, with its name.
        assert conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE sector_id IS NULL").fetchone()["n"] == 0
    finally:
        conn.close()


def test_regenerate_refuses_to_drop_facilities_unless_told(seeded_sector, admin_client):
    config, sector_id, system_ids = seeded_sector
    system_id = system_ids[1]
    star_id = admin_client.get(f"/api/systems/{system_id}").get_json()["stars"][0]["id"]
    response = admin_client.post("/api/facilities", json={
        "name": "Relay One", "kind": "station", "placement": "orbital", "host_type": "star", "host_id": star_id})
    assert response.status_code == 201, response.get_json()

    response = admin_client.patch(f"/api/systems/{system_id}", json={"regenerate": {"planets": False}})
    assert response.status_code == 409
    response = admin_client.patch(f"/api/systems/{system_id}", json={"regenerate": {"planets": False},
                                                                      "drop_facilities": True})
    assert response.status_code == 200
    assert admin_client.get(f"/api/systems/{system_id}/facilities").get_json()["items"] == []


def test_regenerate_rejects_bad_bodies(admin_client):
    system_id = admin_client.post("/api/systems", json={"planets": False}).get_json()["id"]
    assert admin_client.patch(f"/api/systems/{system_id}", json={}).status_code == 400
    assert admin_client.patch(f"/api/systems/{system_id}", json={"regenerate": []}).status_code == 400
    assert admin_client.patch(f"/api/systems/{system_id}", json={"regenerate": {"age": "ancient"}}).status_code == 400
    assert admin_client.patch(f"/api/systems/{system_id}", json={"regenerate": {"name": "X"}}).status_code == 400
    assert admin_client.patch(f"/api/systems/{system_id}", json={"regenerate": {}, "drop_facilities": 1}
                              ).status_code == 400
    assert admin_client.patch("/api/systems/999999", json={"regenerate": {}}).status_code == 404


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


def _save_wide_binary_with_moons(mysql_config):
    """Saves a wide binary whose primary has a planet with a moon, and
    returns its system id -- retried, since generation is random."""
    for _ in range(60):
        cfg = SystemConfig()
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = True
        system = StarSystem(system_config=cfg)
        if any(p.body_type != "a" and p.moons for p in system.planets):
            return _db.save_system(system, cfg, config=mysql_config)
    raise AssertionError("could not generate a wide binary with a moon")


def _names(mysql_config, table, system_id):
    conn = _db.get_connection(mysql_config)
    try:
        return {row["id"]: row["name"] for row in conn.execute(
            f"SELECT id, name FROM {table} WHERE star_system_id = ?", (system_id,)).fetchall()}
    finally:
        conn.close()


def test_rename_system_carries_its_stars_planets_and_moons(admin_client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    old = admin_client.get(f"/api/systems/{system_id}").get_json()["name"]

    old_first = old.split()[0]
    response = admin_client.patch(f"/api/systems/{system_id}", json={"name": "  Castor   Major "})
    assert response.status_code == 200
    assert response.get_json()["name"] == "Castor Major"
    # A wide pair's stars and its primary's planets carry the system name's
    # first word (GEN.62); the secondary's planets carry its own word.
    stars = _names(mysql_config, "stars", system_id).values()
    assert stars and all(name.split()[0] == "Castor" and len(name.split()) == 2 for name in stars), stars
    conn = _db.get_connection(mysql_config)
    try:
        primary_bodies = [row["name"] for row in conn.execute(
            "SELECT p.name FROM planets p JOIN stars s ON s.id = p.star_id "
            "WHERE p.star_system_id = ? AND s.role = 'primary' "
            "UNION ALL SELECT m.name FROM moons m JOIN planets p ON p.id = m.planet_id "
            "JOIN stars s ON s.id = p.star_id WHERE p.star_system_id = ? AND s.role = 'primary'",
            (system_id, system_id)).fetchall()]
    finally:
        conn.close()
    assert primary_bodies and all(name.startswith("Castor ") for name in primary_bodies), primary_bodies
    for table in ("stars", "planets", "moons"):
        assert not any(name.startswith(old_first + " ") for name in _names(mysql_config, table, system_id).values())


def test_rename_star_planet_and_moon(admin_client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    stars = _names(mysql_config, "stars", system_id)
    conn = _db.get_connection(mysql_config)
    try:
        primary_id = conn.execute(
            "SELECT id FROM stars WHERE star_system_id = ? AND role = 'primary'", (system_id,)).fetchone()["id"]
        planet = conn.execute(
            "SELECT p.id, p.name FROM planets p JOIN moons m ON m.planet_id = p.id "
            "WHERE p.star_id = ? ORDER BY p.orbital_index LIMIT 1", (primary_id,)).fetchone()
        moon_id = conn.execute("SELECT id FROM moons WHERE planet_id = ? LIMIT 1", (planet["id"],)).fetchone()["id"]
    finally:
        conn.close()
    # The primary's planets carry the pair's shared first word (GEN.62).
    numeral = planet["name"][len(stars[primary_id].split()[0]) + 1:]

    # A binary's star is renamed on its own, and its planets and moons follow.
    response = admin_client.patch(f"/api/stars/{primary_id}", json={"name": "Castor Pollux"})
    assert response.status_code == 200
    assert _names(mysql_config, "stars", system_id)[primary_id] == "Castor Pollux"
    assert _names(mysql_config, "planets", system_id)[planet["id"]] == f"Castor {numeral}"
    assert _names(mysql_config, "moons", system_id)[moon_id].startswith(f"Castor {numeral}")

    # A planet renamed by hand keeps its moons' names.
    response = admin_client.patch(f"/api/planets/{planet['id']}", json={"name": "New Terra"})
    assert response.status_code == 200
    assert _names(mysql_config, "planets", system_id)[planet["id"]] == "New Terra"
    assert _names(mysql_config, "moons", system_id)[moon_id].startswith(f"Castor {numeral}")

    response = admin_client.patch(f"/api/moons/{moon_id}", json={"name": "Selene"})
    assert response.status_code == 200
    assert _names(mysql_config, "moons", system_id)[moon_id] == "Selene"

    # Only uniquely named objects (sectors, systems, stars) are checked:
    # another body's name is allowed, a star's is refused.
    response = admin_client.patch(f"/api/moons/{moon_id}", json={"name": "New Terra"})
    assert response.status_code == 200
    response = admin_client.patch(f"/api/moons/{moon_id}", json={"name": "Castor Pollux"})
    assert response.status_code == 409
    assert "star" in response.get_json()["error"]


def test_rename_a_single_star_renames_its_system(admin_client, seeded_sector, mysql_config):
    _config, _sector_id, system_ids = seeded_sector
    star_id = next(iter(_names(mysql_config, "stars", system_ids[0])))

    response = admin_client.patch(f"/api/stars/{star_id}", json={"name": "Sirius"})
    assert response.status_code == 200
    assert admin_client.get(f"/api/systems/{system_ids[0]}").get_json()["name"] == "Sirius"
    assert _names(mysql_config, "stars", system_ids[0])[star_id] == "Sirius"

    # Renaming to its current name isn't a clash with itself...
    assert admin_client.patch(f"/api/systems/{system_ids[0]}", json={"name": "Sirius"}).status_code == 200
    # ...but another system's name is.
    response = admin_client.patch(f"/api/systems/{system_ids[1]}", json={"name": "Sirius"})
    assert response.status_code == 409
    response = admin_client.patch(f"/api/systems/{system_ids[1]}", json={"name": "Test Sector"})
    assert response.status_code == 409  # the sector's name


def test_rename_requires_an_admin(client):
    for kind in ("stars", "planets", "moons"):
        assert client.patch(f"/api/{kind}/1", json={"name": "x"}).status_code == 401


def test_rename_validation(admin_client):
    assert admin_client.patch("/api/planets/999999999", json={"name": "Nope"}).status_code == 404
    assert admin_client.patch("/api/moons/999999999", json={"name": "Nope"}).status_code == 404
    assert admin_client.patch("/api/stars/999999999", json={"name": "Nope"}).status_code == 404
    assert admin_client.patch("/api/planets/1", json={"name": "   "}).status_code == 400
    assert admin_client.patch("/api/planets/1", json={"name": 7}).status_code == 400
    assert admin_client.patch("/api/planets/1", json={"name": "x" * 256}).status_code == 400
    assert admin_client.patch("/api/planets/1", json={"name": "ok", "mass": 1}).status_code == 400


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

    monkeypatch.setattr("planetgen.api.routes.WikiClient", RaisingClient)
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


@pytest.mark.parametrize("wiki_url", [
    "javascript:alert(document.cookie)", "data:text/html,<script>alert(1)</script>", "ftp://wiki.example.com/x",
    "//wiki.example.com/x", "https://", "wiki.example.com/Sector",
])
def test_update_sector_refuses_a_wiki_url_that_is_not_http(admin_client, seeded_sector, wiki_url):
    """The sector page links to `wiki_url`, so only an http(s) URL with a
    host is stored; the /admin form's manual link goes through this same
    PATCH, so it's refused there too."""
    _config, sector_id, _system_ids = seeded_sector
    response = admin_client.patch(f"/api/sectors/{sector_id}", json={"wiki_url": wiki_url})
    assert response.status_code == 400
    assert "wiki_url" in response.get_json()["error"]
    assert admin_client.get(f"/api/sectors/{sector_id}").get_json()["wiki_url"] is None


def test_create_sector_rejects_wiki_url(admin_client):
    # A brand-new sector has never been uploaded anywhere -- wiki_url is
    # an update-only field (SECTOR_UPDATE_FIELDS), not accepted at create.
    response = admin_client.post(
        "/api/sectors", json={"name": "Test", "edge_ly": 10.0, "wiki_url": "https://wiki.example.com/x"}
    )
    assert response.status_code == 400


def test_galaxy_cell_describes_any_address_or_point(client, mysql_config):
    """`GET /api/galaxy/cell` answers for any place in the galaxy, by
    address or by a point inside the cell, with its coordinates and 8
    corners, and names the generated sector there if one exists."""
    from planetgen.galaxy.geometry import sector_position_pc
    from stellarObjects.utils import ly_to_pc

    from planetgen import tuning

    edge_pc = float(tuning.DEFAULT_SECTOR_EDGE_PC)
    center = sector_position_pc(3, -1, 5, edge_pc)
    sector_id = _place_sector(mysql_config, "Cell Sector", center, edge_ly=tuning.DEFAULT_SECTOR_EDGE_LY,
                              address=(3, -1, 5))

    by_address = client.get("/api/galaxy/cell?ring=3&layer=-1&slot=5").get_json()
    assert (by_address["ring_index"], by_address["layer_index"], by_address["ring_slot_index"]) == (3, -1, 5)
    assert by_address["sector_id"] == sector_id
    assert by_address["cartesian_pc"] == pytest.approx(list(center))
    assert len(by_address["vertices_pc"]) == 8
    assert by_address["mean_arc_length_pc"] == pytest.approx(edge_pc, rel=0.15)
    assert by_address["spherical"]["polar_rad"] > math.pi / 2  # below the plane

    by_point = client.get(f"/api/galaxy/cell?x={center[0] + 0.3}&y={center[1]}&z={center[2] - 0.2}").get_json()
    assert by_point["designation"] == by_address["designation"]

    empty = client.get("/api/galaxy/cell?x=9000&y=-120&z=4000").get_json()
    assert empty["sector_id"] is None and len(empty["vertices_pc"]) == 8

    assert client.get("/api/galaxy/cell?ring=0&layer=0&slot=4").status_code == 400
    assert client.get("/api/galaxy/cell?ring=1").status_code == 400
    assert client.get("/api/galaxy/cell?x=nan&y=0&z=0").status_code == 400


def test_facility_routes(admin_client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        moon_id = conn.execute("SELECT id FROM moons WHERE star_system_id = ? LIMIT 1", (system_id,)).fetchone()["id"]
    finally:
        conn.close()

    orbit = admin_client.get(f"/api/facilities/orbit?host_type=moon&host_id={moon_id}").get_json()
    assert orbit["period_years"] > 0 and orbit["orbital_speed_kms"] > 0
    assert 0 < orbit["min_distance_km"] < orbit["distance_km"] <= orbit["max_distance_km"]
    assert admin_client.get("/api/facilities/orbit?host_type=moon&host_id=999999999").status_code == 404

    response = admin_client.post("/api/facilities", json={
        "name": "Moonport", "kind": "station", "placement": "orbital", "host_type": "moon", "host_id": moon_id,
    })
    assert response.status_code == 201
    facility_id = response.get_json()["id"]
    assert admin_client.post("/api/facilities", json={
        "name": "Bad", "kind": "mining-colony", "placement": "orbital", "host_type": "moon", "host_id": moon_id,
    }).status_code == 400

    detail = admin_client.get(f"/api/facilities/{facility_id}").get_json()
    assert (detail["name"], detail["host_type"], detail["host_id"]) == ("Moonport", "moon", moon_id)
    listed = admin_client.get(f"/api/systems/{system_id}/facilities").get_json()["items"]
    assert [f["id"] for f in listed] == [facility_id]

    assert admin_client.delete(f"/api/facilities/{facility_id}").status_code == 200
    assert admin_client.get(f"/api/facilities/{facility_id}").status_code == 404
