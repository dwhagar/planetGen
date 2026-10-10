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


def test_scattered_phenomena_not_yet_built_are_listed_for_every_kind(client, mysql_config):
    """Boss's empty table: a fresh scatter fills `phenomenon_scatter` and nothing else, so each unbuilt row must be listed."""
    kinds = [("black-hole", "stellar"), ("black-hole", "supermassive"), ("neutron-star", None),
             ("planetary-nebula", None), ("supernova-remnant", None), ("hypervelocity-star", None),
             ("quasar", None)]
    rows = [(0, 0, slot, kind, subtype, 1000 * slot, 0, 0, None, None, None, slot + 1, None)
            for slot, (kind, subtype) in enumerate(kinds)]
    conn = store.get_connection(mysql_config)
    store.insert_phenomenon_scatter(conn, rows)
    store.record_phenomenon_scatter_classes(conn, {(kind, subtype or ""): 1 for kind, subtype in kinds})
    quasar = conn.execute("SELECT id FROM phenomenon_scatter WHERE kind = 'quasar'").fetchone()["id"]
    store.mark_phenomena_built(conn, [quasar])
    conn.commit()
    body = client.get("/api/phenomena?facets=1").get_json()
    assert body["total"] == len(kinds) - 1
    assert sorted({item["type"] for item in body["items"]}) == [
        "black_hole", "hypervelocity_star", "nebula", "neutron_star", "supernova_remnant"]
    assert all(item["scattered"] and item["placed"] for item in body["items"])
    assert {f["value"]: f["count"] for f in body["facets"]["type"]}["black_hole"] == 2
    assert client.get("/api/phenomena?type=neutron_star").get_json()["total"] == 1
    assert "Uncharted neutron star 0.0.2" in _names(client.get("/api/phenomena"))


def test_scattered_phenomena_are_paged_after_the_built_ones_without_scanning(client, mysql_config):
    _seed(mysql_config)
    conn = store.get_connection(mysql_config)
    rows = [(0, 0, slot, "neutron-star", None, slot, 0, 0, None, None, None, slot + 1, None) for slot in range(5)]
    store.insert_phenomenon_scatter(conn, rows)
    store.record_phenomenon_scatter_classes(conn, {("neutron-star", ""): 5})
    conn.commit()
    everything = client.get("/api/phenomena").get_json()
    assert everything["total"] == 9
    assert [i["name"] for i in everything["items"]][4:] == [f"Uncharted neutron star 0.0.{n}" for n in range(5)]
    # A page that straddles the built rows and the scattered ones, and one wholly inside the scattered ones.
    straddle = client.get("/api/phenomena?limit=3&offset=3").get_json()
    assert [i["name"] for i in straddle["items"]] == ["Zeta Cloud", "Uncharted neutron star 0.0.0",
                                                       "Uncharted neutron star 0.0.1"]
    deep = client.get("/api/phenomena?limit=2&offset=7").get_json()
    assert [i["name"] for i in deep["items"]] == ["Uncharted neutron star 0.0.3", "Uncharted neutron star 0.0.4"]
    # Building a scattered phenomenon takes it off the unbuilt count.
    store.mark_phenomena_built(conn, [conn.execute("SELECT MIN(id) AS i FROM phenomenon_scatter").fetchone()["i"]])
    conn.commit()
    assert client.get("/api/phenomena").get_json()["total"] == 8
    assert client.get("/api/phenomena?type=neutron_star").get_json()["total"] == 4
