# tests/test_api_list_tables.py

"""
The Sectors and Systems tables' server side (UX.41): `GET /api/sectors` and
`GET /api/systems` sort, filter and facet counts (`queryDb.list_sectors`,
`list_systems`, `sectors_facets`, `systems_facets`).
"""

from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from tests.test_api import client, seeded_sector  # noqa: F401


def _standalone(mysql_config, name, binary=False):
    cfg = SystemConfig()
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = binary
    system = StarSystem(system_config=cfg)
    system.name = name
    return store.save_system(system, cfg, config=mysql_config)


def _names(response):
    assert response.status_code == 200, response.get_data(as_text=True)
    return [item["name"] for item in response.get_json()["items"]]


def _second_sector(mysql_config):
    sector = SpaceSector("Another Sector", edge_ly=10.0)
    cfg = SystemConfig()
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=cfg), position=(0.0, 0.0, 0.0), system_config=cfg)
    return store.save_sector(sector, config=mysql_config)


# --- sectors -----------------------------------------------------------------------

def test_sectors_sort_by_name_and_system_count(client, seeded_sector):
    mysql_config, _sector_id, _systems = seeded_sector
    _second_sector(mysql_config)
    assert _names(client.get("/api/sectors?sort=name")) == ["Another Sector", "Test Sector"]
    assert _names(client.get("/api/sectors?sort=name&order=desc")) == ["Test Sector", "Another Sector"]
    # Two systems against one.
    assert _names(client.get("/api/sectors?sort=systems&order=desc"))[0] == "Test Sector"
    assert _names(client.get("/api/sectors?sort=systems"))[0] == "Another Sector"


def test_sectors_without_a_sort_keep_the_core_first_order(client, seeded_sector):
    mysql_config, _sector_id, _systems = seeded_sector
    _second_sector(mysql_config)
    body = client.get("/api/sectors").get_json()
    assert [item["name"] for item in body["items"]] == ["Another Sector", "Test Sector"]  # unplaced: by name
    assert "facets" not in body


def test_sectors_quadrant_filter_and_facets(client, seeded_sector):
    mysql_config, _sector_id, _systems = seeded_sector
    _second_sector(mysql_config)
    # Neither saved sector has a galaxy position.
    body = client.get("/api/sectors?quadrant=unplaced&facets=1").get_json()
    assert body["total"] == 2
    assert body["facets"] == {"quadrant": [{"value": "unplaced", "count": 2}]}
    none = client.get("/api/sectors?quadrant=I").get_json()
    assert none["items"] == [] and none["total"] == 0
    # A menu ignores its own filter: the counts stay put.
    assert client.get("/api/sectors?quadrant=I&facets=1").get_json()["facets"]["quadrant"][0]["count"] == 2


def test_sectors_reject_a_bad_sort(client, seeded_sector):
    for query in ("sort=mass", "order=up"):
        assert client.get(f"/api/sectors?{query}").status_code == 400


# --- systems -----------------------------------------------------------------------

def test_systems_sort_by_name_sector_and_binary(client, seeded_sector):
    mysql_config, _sector_id, system_ids = seeded_sector
    _standalone(mysql_config, "Aardvark")
    _standalone(mysql_config, "Zebra", binary=True)
    names = _names(client.get("/api/systems?sort=name&order=desc"))
    assert names == sorted(names, reverse=True) and names[0] == "Zebra"
    by_sector = _names(client.get("/api/systems?sort=sector"))
    assert by_sector[-2:] == ["Aardvark", "Zebra"]  # standalone ones last
    by_sector_desc = _names(client.get("/api/systems?sort=sector&order=desc"))
    assert by_sector_desc[-2:] == ["Aardvark", "Zebra"]
    binary_first = client.get("/api/systems?sort=binary&order=desc").get_json()["items"][0]
    assert binary_first["name"] == "Zebra" and binary_first["is_binary"]


def test_systems_filters(client, seeded_sector):
    mysql_config, sector_id, _systems = seeded_sector
    _standalone(mysql_config, "Aardvark")
    _standalone(mysql_config, "Zebra", binary=True)
    assert client.get("/api/systems?binary=yes").get_json()["total"] == 1
    assert client.get("/api/systems?binary=no").get_json()["total"] == 3
    assert client.get("/api/systems?placement=standalone").get_json()["total"] == 2
    assert client.get("/api/systems?placement=sector").get_json()["total"] == 2
    octant = client.get("/api/systems?placement=sector").get_json()["items"][0]["quadrant"]
    assert client.get(f"/api/systems?octant={octant}").get_json()["total"] >= 1
    paged = client.get("/api/systems?placement=standalone&limit=1&offset=1").get_json()
    assert [item["name"] for item in paged["items"]] == ["Zebra"] and paged["total"] == 2


def test_systems_facets_narrow_each_other(client, seeded_sector):
    mysql_config, _sector_id, _systems = seeded_sector
    _standalone(mysql_config, "Aardvark")
    _standalone(mysql_config, "Zebra", binary=True)
    facets = client.get("/api/systems?facets=1").get_json()["facets"]
    assert {o["value"]: o["count"] for o in facets["placement"]} == {"sector": 2, "standalone": 2}
    assert {o["value"]: o["count"] for o in facets["binary"]} == {"yes": 1, "no": 3}
    assert sum(o["count"] for o in facets["octant"]) == 2
    narrowed = client.get("/api/systems?facets=1&binary=yes").get_json()["facets"]
    # The binary menu still counts both ways; the placement menu counts only binaries.
    assert {o["value"]: o["count"] for o in narrowed["binary"]} == {"yes": 1, "no": 3}
    assert {o["value"]: o["count"] for o in narrowed["placement"]} == {"standalone": 1}
    both = client.get("/api/systems?facets=1&sector_id=none").get_json()["facets"]
    assert {o["value"]: o["count"] for o in both["placement"]} == {"standalone": 2}


def test_systems_reject_a_bad_sort_or_filter(client, seeded_sector):
    for query in ("sort=mass", "order=up", "binary=maybe", "placement=orbit"):
        assert client.get(f"/api/systems?{query}").status_code == 400, query
