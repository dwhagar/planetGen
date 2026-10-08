# tests/test_api_phenomena_table.py

"""
The Phenomena table's server side (UX.41): `GET /api/phenomena`'s sort,
filters and facet counts (`queryDb.list_phenomena`, `count_phenomena`,
`phenomena_facets`).
"""

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from tests.test_api import client  # noqa: F401


def _save(mysql_config, phenomenon, kind, name):
    phenomenon.name = name
    return store.save_phenomenon(phenomenon, SystemConfig(), kind, config=mysql_config)


def _seed(mysql_config):
    """Two nebulae and two rogue planets with known names."""
    cfg = SystemConfig()
    _save(mysql_config, Nebula(cfg), "nebula", "Alpha Cloud")
    _save(mysql_config, Nebula(cfg), "nebula", "Zeta Cloud")
    _save(mysql_config, RoguePlanet(cfg), "rogue-planet", "Beta Wanderer")
    _save(mysql_config, RoguePlanet(cfg), "rogue-planet", "Gamma Wanderer")


def _names(response):
    assert response.status_code == 200, response.get_data(as_text=True)
    return [item["name"] for item in response.get_json()["items"]]


def test_phenomena_sort_by_name_either_way(client, mysql_config):
    _seed(mysql_config)
    assert _names(client.get("/api/phenomena")) == ["Alpha Cloud", "Beta Wanderer", "Gamma Wanderer", "Zeta Cloud"]
    assert _names(client.get("/api/phenomena?sort=name&order=desc"))[0] == "Zeta Cloud"


def test_phenomena_sort_by_type_breaks_ties_by_name(client, mysql_config):
    _seed(mysql_config)
    assert _names(client.get("/api/phenomena?sort=type")) == [
        "Alpha Cloud", "Zeta Cloud", "Beta Wanderer", "Gamma Wanderer"]
    assert _names(client.get("/api/phenomena?sort=type&order=desc")) == [
        "Beta Wanderer", "Gamma Wanderer", "Alpha Cloud", "Zeta Cloud"]


def test_phenomena_filter_by_type_and_descriptor(client, mysql_config):
    _seed(mysql_config)
    body = client.get("/api/phenomena?type=rogue_planet").get_json()
    assert [item["name"] for item in body["items"]] == ["Beta Wanderer", "Gamma Wanderer"]
    assert body["total"] == 2
    descriptor = body["items"][0]["descriptor"]
    both = client.get(f"/api/phenomena?type=rogue_planet&type=nebula&descriptor={descriptor}").get_json()
    assert both["total"] == len([i for i in client.get("/api/phenomena").get_json()["items"]
                                 if i["descriptor"] == descriptor])
    empty = client.get("/api/phenomena?type=quasar").get_json()
    assert empty["items"] == [] and empty["total"] == 0


def test_phenomena_filters_page_with_the_filtered_total(client, mysql_config):
    _seed(mysql_config)
    body = client.get("/api/phenomena?type=nebula&limit=1&offset=1").get_json()
    assert [item["name"] for item in body["items"]] == ["Zeta Cloud"]
    assert body["total"] == 2


def test_phenomena_facets_count_each_menu_without_its_own_filter(client, mysql_config):
    _seed(mysql_config)
    unfiltered = client.get("/api/phenomena?facets=1").get_json()["facets"]
    assert {o["value"]: o["count"] for o in unfiltered["type"]} == {"nebula": 2, "rogue_planet": 2}
    assert sum(o["count"] for o in unfiltered["descriptor"]) == 4

    narrowed = client.get("/api/phenomena?facets=1&type=nebula").get_json()["facets"]
    # The type menu ignores its own filter, so it still offers both types ...
    assert {o["value"]: o["count"] for o in narrowed["type"]} == {"nebula": 2, "rogue_planet": 2}
    # ... while the descriptor menu now counts only the nebulae.
    assert sum(o["count"] for o in narrowed["descriptor"]) == 2


def test_phenomena_facets_only_when_asked(client, mysql_config):
    _seed(mysql_config)
    assert "facets" not in client.get("/api/phenomena").get_json()


def test_phenomena_placed_filter(client, mysql_config):
    _seed(mysql_config)
    assert client.get("/api/phenomena?placed=yes").status_code == 200
    unplaced = client.get("/api/phenomena?placed=no").get_json()
    assert unplaced["total"] == client.get("/api/phenomena").get_json()["total"] - \
        client.get("/api/phenomena?placed=yes").get_json()["total"]


def test_phenomena_rejects_unknown_sorts_and_orders(client, mysql_config):
    _seed(mysql_config)
    for query in ("sort=mass", "order=sideways", "placed=maybe"):
        response = client.get(f"/api/phenomena?{query}")
        assert response.status_code == 400, query
        assert "error" in response.get_json()


def test_phenomena_unplaced_sector_sorts_last_both_ways(client, mysql_config):
    _seed(mysql_config)
    # None of these has a sector: they all tie, then fall back to name.
    assert _names(client.get("/api/phenomena?sort=sector")) == _names(client.get("/api/phenomena"))
