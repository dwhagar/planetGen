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

import pytest

from api.app import create_app
from api.config import Config
from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem


@pytest.fixture
def seeded_sector(mysql_config):
    """Saves one small sector (two systems, one G2V, one M5V) and returns
    (mysql_config, sector_id, system_ids) for tests to build requests
    against."""
    sector = SpaceSector("Test Sector", edge_ly=10.0)

    cfg_a = SystemConfig()
    cfg_a.STAR_TYPE = "G2V"
    cfg_a.PLANETS = False
    system_a = StarSystem(system_config=cfg_a)
    sector.add_system(system_a, position=(1.0, 1.0, 1.0), system_config=cfg_a)

    cfg_b = SystemConfig()
    cfg_b.STAR_TYPE = "M5V"
    cfg_b.PLANETS = False
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
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def test_health_ok(client):
    response = client.get("/api/health")
    assert response.status_code == 200
    assert response.get_json() == {"status": "ok"}


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
    assert "markdown_content" in body
    assert {s["id"] for s in body["sector_siblings"]} == set(system_ids)

    response = client.get("/api/systems/999999999")
    assert response.status_code == 404


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


def test_unmatched_route_returns_json_404(client):
    response = client.get("/api/no-such-route")
    assert response.status_code == 404
    assert response.get_json() == {"error": "not found"}


# ---------------------------------------------------------------------
# Write endpoints -- stubs (see routes.py's module docstring). These
# never touch the database, so they don't need `seeded_sector`/a real
# sector or system id to exist first -- every one of them responds 501
# regardless, as long as the request body passes validation.
# ---------------------------------------------------------------------

def test_create_sector_validates_then_returns_not_implemented(client):
    response = client.post("/api/sectors", json={"name": "Test", "edge_ly": 10.0})
    assert response.status_code == 501
    assert "error" in response.get_json()

    response = client.post("/api/sectors", json={"name": "Test"})
    assert response.status_code == 400

    response = client.post("/api/sectors", json={"name": "", "edge_ly": 10.0})
    assert response.status_code == 400

    response = client.post("/api/sectors", json={"name": "Test", "edge_ly": -1})
    assert response.status_code == 400

    response = client.post("/api/sectors", json={"name": "Test", "edge_ly": 10.0, "bogus": 1})
    assert response.status_code == 400

    response = client.post("/api/sectors", data="not json", content_type="text/plain")
    assert response.status_code == 400


def test_update_sector_validates_then_returns_not_implemented(client):
    response = client.patch("/api/sectors/1", json={"name": "Renamed"})
    assert response.status_code == 501

    response = client.patch("/api/sectors/1", json={"edge_ly": 5.0})
    assert response.status_code == 501

    response = client.patch("/api/sectors/1", json={})
    assert response.status_code == 400

    response = client.patch("/api/sectors/1", json={"edge_ly": -1})
    assert response.status_code == 400


def test_delete_sector_returns_not_implemented(client):
    response = client.delete("/api/sectors/1")
    assert response.status_code == 501
    assert "error" in response.get_json()


def test_create_system_requires_json_object_then_returns_not_implemented(client):
    response = client.post("/api/systems", json={"name": "Testworld"})
    assert response.status_code == 501

    response = client.post("/api/systems", data="not json", content_type="text/plain")
    assert response.status_code == 400


def test_update_system_requires_json_object_then_returns_not_implemented(client):
    response = client.patch("/api/systems/1", json={"star_type": "G2V"})
    assert response.status_code == 501


def test_delete_system_returns_not_implemented(client):
    response = client.delete("/api/systems/1")
    assert response.status_code == 501


def test_write_endpoints_are_rate_limited_more_tightly_than_the_default(client):
    # WRITE_RATE_LIMIT is 10/minute -- the 11th write in one minute must
    # be rejected with 429, well before the 50/hour global default would
    # ever kick in.
    for _ in range(10):
        response = client.delete("/api/sectors/1")
        assert response.status_code == 501

    response = client.delete("/api/sectors/1")
    assert response.status_code == 429
    body = response.get_json()
    assert body["error"] == "rate limit exceeded"
