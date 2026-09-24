# tests/test_admin_stats.py

"""
Tests for the admin stats endpoints (`html/api/admin.py`, backed by
`src/adminStats.py`): `GET /api/admin/stats` and `GET
/api/admin/duplicate-names`. Needs a MySQL test server like every other
database-backed test (see `conftest.py`).
"""

import pytest

import adminStats
from api.app import create_app
from api.config import Config
from stellarObjects import _db, adminAuth
from stellarObjects.config import SystemConfig
from stellarObjects.nameUniqueness import strip_decoration
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem


@pytest.fixture
def client(mysql_config):
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def _login(client, change_credentials=True):
    response = client.post("/api/auth/login", json={
        "username": adminAuth.DEFAULT_ADMIN_USERNAME, "password": adminAuth.DEFAULT_ADMIN_PASSWORD,
    })
    assert response.status_code == 200
    if change_credentials:
        response = client.post("/api/auth/change-credentials", json={
            "current_password": adminAuth.DEFAULT_ADMIN_PASSWORD,
            "new_username": "stats-admin",
            "new_password": "a-strong-test-password-123",
        })
        assert response.status_code == 200


@pytest.fixture
def admin_client(mysql_config, client):
    adminAuth.bootstrap_control_schema(mysql_config)
    _db.get_connection(mysql_config).close()
    _login(client)
    return client


def _insert_system(conn, name):
    cfg = SystemConfig()
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    system.star.name = name
    with conn:
        return _db.insert_star_system(conn, system, cfg)


@pytest.fixture
def colliding_names(mysql_config):
    """Two sectors named "Sol", two systems named "Terra", and a system
    named after a sector ("Mars"), so every registry except the body one
    has a row."""
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_sector(conn, SpaceSector(name="Sol"))
        with conn:
            _db.insert_sector(conn, SpaceSector(name="Sol"))
        with conn:
            mars_sector_id = _db.insert_sector(conn, SpaceSector(name="Mars"))
        _insert_system(conn, "Terra")
        _insert_system(conn, "Terra")
        little_mars_id = _insert_system(conn, "Mars")
        _insert_system(conn, "Unique Vega")
    finally:
        conn.close()
    return {"mars_sector_id": mars_sector_id, "little_mars_id": little_mars_id}


def test_admin_endpoints_require_a_fresh_admin(mysql_config, client):
    assert client.get("/api/admin/stats").status_code == 401
    assert client.get("/api/admin/duplicate-names").status_code == 401

    adminAuth.bootstrap_control_schema(mysql_config)
    _login(client, change_credentials=False)
    assert client.get("/api/admin/stats").status_code == 403
    assert client.get("/api/admin/duplicate-names").status_code == 403


def test_stats_reports_health_and_database_numbers(admin_client, colliding_names):
    response = admin_client.get("/api/admin/stats")
    assert response.status_code == 200
    body = response.get_json()

    assert body["api"]["version"]
    assert body["api"]["uptime_seconds"] >= 0
    assert body["mysql"]["version"]

    database = body["database"]
    assert database["reachable"] is True
    assert database["schema_version"] == _db.SCHEMA_VERSION
    assert database["schema_current"] is True
    assert database["counts"] == {"sectors": 3, "star_systems": 4}
    assert {"sectors", "star_systems", "planets", "moons"} <= {t["name"] for t in database["tables"]}

    stamps = {t["table"]: t for t in database["timestamps"]}
    assert set(stamps) == set(_db.TIMESTAMPED_TABLES)
    assert stamps["star_systems"]["newest_created_at"] is not None
    assert stamps["star_systems"]["last_modified_at"] is not None
    assert stamps["black_holes"]["last_modified_at"] is None

    collisions = database["name_collisions"]
    assert collisions["sector"] == 1
    assert collisions["system"] == 2
    assert collisions["body"] == 0
    assert collisions["distinct_base_names"] == 3


def test_duplicate_names_lists_every_row_for_each_base_name(admin_client, colliding_names):
    response = admin_client.get("/api/admin/duplicate-names")
    assert response.status_code == 200
    body = response.get_json()
    assert body["total"] == 3

    items = {item["base_name"]: item for item in body["items"]}
    assert list(items) == ["Mars", "Sol", "Terra"]

    assert items["Sol"]["levels"] == ["sector"]
    assert [(r["kind"], r["name"]) for r in items["Sol"]["rows"]] == [("sector", "Alpha Sol"), ("sector", "Beta Sol")]

    assert [(r["kind"], r["name"]) for r in items["Terra"]["rows"]] == [
        ("system", "Alpha Terra"), ("system", "Beta Terra"),
    ]
    assert all("sector_id" in r for r in items["Terra"]["rows"])

    mars = items["Mars"]
    assert mars["levels"] == ["system"]
    assert [(r["kind"], r["id"], r["name"]) for r in mars["rows"]] == [
        ("sector", colliding_names["mars_sector_id"], "Mars"),
        ("system", colliding_names["little_mars_id"], "Little Mars"),
    ]


def test_duplicate_names_paginates_by_base_name(admin_client, colliding_names):
    body = admin_client.get("/api/admin/duplicate-names?limit=1&offset=1").get_json()
    assert body["total"] == 3
    assert [item["base_name"] for item in body["items"]] == ["Sol"]
    assert admin_client.get("/api/admin/duplicate-names?limit=0").status_code == 400


def test_duplicate_names_includes_planets_with_their_system(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_sector(conn, SpaceSector(name="Kepler"))
        for _ in range(30):
            cfg = SystemConfig()
            cfg.PLANETS = True
            cfg.MAX_PLANETS = True
            system = StarSystem(system_config=cfg)
            planets = [p for p in system.planets if p.body_type != "a"]
            if planets:
                planets[0].name = "Kepler"
                with conn:
                    system_id = _db.insert_star_system(conn, system, cfg)
                break
        else:
            pytest.fail("could not generate a system with a real planet")

        # Other bodies in the generated system can collide on their own,
        # so only the forced one is checked.
        result = adminStats.duplicate_names(conn)
        (item,) = [item for item in result["items"] if item["base_name"] == "Kepler"]
        assert item["levels"] == ["body"]
        planet_rows = [r for r in item["rows"] if r["kind"] == "planet"]
        assert [r["name"] for r in planet_rows] == ["Kepler Kin"]
        assert planet_rows[0]["star_system_id"] == system_id
        assert planet_rows[0]["system_name"] == system.star.name
    finally:
        conn.close()


def test_every_decorated_candidate_strips_back_to_its_base_name():
    for name in adminStats.decorated_name_candidates("Nova Prime"):
        assert strip_decoration(name) == "Nova Prime"
