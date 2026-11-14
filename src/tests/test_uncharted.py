# tests/test_uncharted.py

"""
MAP.400: the objects no sector has stored yet (the scattered black holes, neutron stars, remnants, planetary nebulae, hypervelocity
stars and the nucleus quasar, and the field clouds) are found by the lists, the counts, the search and Nearby, though the object index
(`db/objectindex.py`) holds neither table. Each scattered kind is listed, counted and found by search; a row a sector built is the
stored object's and drops out of the tail.
"""

import pytest

from planetgen.db import near, objectindex, objectlist, query, store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from tests.test_api import client  # noqa: F401

KINDS = [("black-hole", "stellar"), ("black-hole", "supermassive"), ("neutron-star", None), ("planetary-nebula", None),
         ("supernova-remnant", None), ("hypervelocity-star", None), ("quasar", None)]
"""list: One scattered row of each kind, the black hole twice (a stellar and the nucleus), at slot parsecs from the origin."""

CLOUDS = [(0, 0, 0, 0, "M", 40.0, 0.0, 5.0, 0.0), (0, 0, 0, 1, "N", 30.0, 0.0, 6.0, 0.0), (0, 0, 1, 0, "M", 20.0, 0.0, 30.0, 0.0)]
"""list: Three field clouds (`nebula_field` place columns), the third one 30 pc out."""

LISTED = {"black_hole": 2, "neutron_star": 1, "supernova_remnant": 1, "quasar": 1, "hypervelocity_star": 1, "nebula": 4}
"""dict: How many uncharted objects each list kind holds in the seed (the nebulae: three field clouds and the planetary nebula)."""


@pytest.fixture
def scattered(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        rows = [(0, 0, slot, kind, subtype, 1000 * slot, 0, 0, None, None, None, slot + 1, 3.0 if subtype == "stellar" else None, None)
                for slot, (kind, subtype) in enumerate(KINDS)]
        store.insert_phenomenon_scatter(conn, rows)
        store.record_phenomenon_scatter_classes(conn, {(kind, subtype or ""): 1 for kind, subtype in KINDS})
        store.insert_nebula_field(conn, CLOUDS)
        conn.commit()
    finally:
        conn.close()
    return mysql_config


def _walk(config, kind, limit, **options):
    """Every item of the pages from the top, the number of pages and the total."""
    conn = store.get_connection(config)
    try:
        items, after, pages = [], None, 0
        while True:
            page = objectlist.list_objects(conn, kind, limit=limit, after=after, **options)
            items += page["items"]
            pages += 1
            after = page["next"]
            if after is None:
                return items, pages, page["total"]
            assert pages < 100
    finally:
        conn.close()


@pytest.mark.parametrize("kind", sorted(LISTED))
def test_each_scattered_kind_is_listed_and_counted(scattered, kind):
    items, pages, total = _walk(scattered, kind, 50)
    assert total == len(items) == LISTED[kind] and pages == 1
    assert all(item["uncharted"] and item["kind"] == kind and item["placed"] for item in items)
    assert len({item["id"] for item in items}) == len(items), "every object has an ID of its own"


def test_a_scattered_row_links_to_its_scatter_page_and_a_field_cloud_has_none(scattered):
    items, _pages, _total = _walk(scattered, "nebula", 50)
    clouds = [item for item in items if item["id"].startswith("field-")]
    assert sorted(item["id"] for item in clouds) == ["field-0.0.0.0", "field-0.0.0.1", "field-0.0.1.0"]
    assert all(item["scatter_id"] is None and item["name"].startswith("Field cloud ") for item in clouds)
    assert {item["object_class"] for item in clouds} == {"dark"} and all(item["radius_km"] > 0 for item in clouds)
    [planetary] = [item for item in items if not item["id"].startswith("field-")]
    assert planetary["scatter_id"] and planetary["name"].startswith("Uncharted planetary nebula 0.0.")
    hole = _walk(scattered, "black_hole", 50)[0][0]
    assert hole["mass_kg"] is not None and hole["name"].startswith("Uncharted black hole 0.0.")


def test_a_walk_in_small_pages_visits_every_uncharted_object_once_across_the_field_and_the_scatter(scattered):
    items, pages, total = _walk(scattered, "nebula", 3)
    assert pages == 2 and total == len(items) == 4
    assert [item["id"] for item in items[:3]] == ["field-0.0.0.0", "field-0.0.0.1", "field-0.0.1.0"]
    assert not items[3]["id"].startswith("field-")
    by_one, pages, _total = _walk(scattered, "nebula", 1)
    assert pages == 4 and [item["id"] for item in by_one] == [item["id"] for item in items]


def test_the_stored_objects_come_first_and_the_walk_crosses_into_the_uncharted_ones(scattered):
    stored = Nebula(SystemConfig())
    stored.name = "Alpha Cloud"
    store.save_phenomenon(stored, SystemConfig(), "nebula", config=scattered)
    conn = store.get_connection(scattered)
    try:
        objectindex.build(conn)
        conn.commit()
    finally:
        conn.close()
    items, pages, total = _walk(scattered, "nebula", 2)
    assert total == len(items) == 5 and pages == 3
    assert items[0]["name"] == "Alpha Cloud" and not items[0].get("uncharted")
    assert all(item.get("uncharted") for item in items[1:])
    descending, _pages, _total = _walk(scattered, "nebula", 2, descending=True)
    assert descending[0]["name"] == "Alpha Cloud", "the stored objects are the ones the sort orders"


def test_a_built_scatter_row_and_a_stored_field_cloud_leave_the_tail(scattered):
    conn = store.get_connection(scattered)
    try:
        hole = conn.execute("SELECT id FROM phenomenon_scatter WHERE kind = 'black-hole' ORDER BY id LIMIT 1").fetchone()["id"]
        store.mark_phenomena_built(conn, [hole])
        conn.execute("UPDATE nebula_field SET is_stored = 1 WHERE cell_k = 1")
        conn.commit()
        assert objectlist.list_objects(conn, "black_hole")["total"] == 1
        assert objectlist.list_objects(conn, "nebula")["total"] == 3
    finally:
        conn.close()


def test_a_name_a_sector_a_range_or_another_sort_leaves_the_uncharted_objects_out_and_a_class_narrows_them(scattered):
    conn = store.get_connection(scattered)
    try:
        assert objectlist.list_objects(conn, "black_hole", name="Uncharted")["total"] == 0
        assert objectlist.list_objects(conn, "black_hole", sector_id=1)["total"] == 0
        assert objectlist.list_objects(conn, "black_hole", sort="mass")["total"] == 0
        assert objectlist.list_objects(conn, "black_hole", object_class="stellar")["total"] == 1
        assert objectlist.list_objects(conn, "nebula", object_class="dark")["total"] == 3
        assert objectlist.list_objects(conn, "nebula", object_class="emission")["total"] == 0
        assert objectlist.list_objects(conn, "hypervelocity_star", sort="mass")["total"] == 1, "a kind only they have lists in id order"
        assert objectlist.list_objects(conn, "hypervelocity_star", start="Z")["total"] == 1
    finally:
        conn.close()


def test_a_cursor_that_is_not_an_object_is_refused_and_one_that_is_gone_starts_from_the_top(scattered):
    conn = store.get_connection(scattered)
    try:
        with pytest.raises(ValueError):
            objectlist.list_objects(conn, "black_hole", after="not-an-id")
        with pytest.raises(ValueError):
            objectlist.list_objects(conn, "nebula", after="field-x.1.2.3")
        # A field cloud's key is a place in the order, so a key past every cloud (or of one a sector has since stored) goes on with the scattered rows.
        past = objectlist.list_objects(conn, "nebula", after="field-9.9.9.9")
        assert [item["scatter_id"] is not None for item in past["items"]] == [True]
        first = objectlist.list_objects(conn, "black_hole", limit=1)
        second = objectlist.list_objects(conn, "black_hole", limit=1, after=first["next"])
        assert second["next"] is None and second["items"][0]["id"] != first["items"][0]["id"]
    finally:
        conn.close()


def test_the_api_and_the_lists_page_count_what_the_map_draws(scattered, client):  # noqa: F811
    body = client.get("/api/object-lists/black_hole?limit=1").get_json()
    assert body["total"] == 2 and body["next"] and body["items"][0]["uncharted"] is True
    assert client.get("/api/object-lists/hypervelocity_star").get_json()["total"] == 1
    page = client.get("/lists").get_data(as_text=True)
    assert "Hypervelocity stars" in page and "Nebulae</a> (4)" in page and "Black holes</a> (2)" in page
    assert "Supernova remnants</a> (1)" in page
    table = client.get("/lists/black_hole").get_data(as_text=True)
    assert "/uncharted-phenomenon/" in table and "Uncharted black hole" in table
    assert client.get("/lists/nebula").status_code == 200


def test_the_search_counts_and_finds_each_scattered_kind(scattered):
    conn = store.get_connection(scattered)
    try:
        tags = {"phenomenon": set(LISTED)}
        found = query.search(conn, {}, tags)
        counts = {facet["value"]: facet["count"] for facet in found["facets"]["phenomenon"]}
        assert counts == LISTED
        panel = found["results"]["phenomena"]
        assert panel["total"] == sum(LISTED.values()) and all(row["uncharted"] for row in panel["rows"])
        assert {row["type"] for row in panel["rows"]} == set(LISTED)
        for kind, count in LISTED.items():
            rows = query.search(conn, {}, {"phenomenon": {kind}})["results"]["phenomena"]
            assert rows["total"] == count == len(rows["rows"]), kind
            assert all(row["ref"] is None and row["parents"] == [] for row in rows["rows"])
        # A class tag of nebulae narrows them to the field clouds of that class: a planetary nebula has no letter class.
        narrowed = query.search(conn, {}, {"phenomenon_class": {"nebula:M"}})["results"]["phenomena"]
        assert [row["id"] for row in narrowed["rows"]] == ["field-0.0.0.0", "field-0.0.1.0"]
        classes = {facet["value"]: facet["count"] for facet in found["facets"]["phenomenon_class"]}
        assert classes["nebula:M"] == 2 and classes["nebula:N"] == 1
    finally:
        conn.close()


def test_the_search_pages_across_the_stored_and_the_uncharted_rows(scattered):
    stored = Nebula(SystemConfig())
    stored.name = "Alpha Cloud"
    store.save_phenomenon(stored, SystemConfig(), "nebula", config=scattered)
    conn = store.get_connection(scattered)
    try:
        seen = []
        for offset in range(0, 5, 2):
            page = query.search(conn, {}, {"phenomenon": {"nebula"}}, limit=2, offsets={"phenomena": offset})["results"]["phenomena"]
            assert page["total"] == 5
            seen += [(row["name"], bool(row.get("uncharted"))) for row in page["rows"]]
        assert seen[0] == ("Alpha Cloud", False) and len(seen) == 5 and all(flag for _name, flag in seen[1:])
        assert len({name for name, _flag in seen}) == 5
        last = query.search(conn, {}, {"phenomenon": {"nebula"}}, limit=2, offsets={"phenomena": 40})["results"]["phenomena"]
        assert last["offset"] == 4 and len(last["rows"]) == 1, "a stale page link is pulled back to the last page"
    finally:
        conn.close()


def test_nearby_lists_the_uncharted_objects_with_no_reference(scattered):
    conn = store.get_connection(scattered)
    try:
        place = near.place_from_point((0.0, 0.0, 0.0))
        page = near.objects_within(conn, place, 10.0)
        assert page["not_generated"] == {} and page["by_kind"] == {"black_hole": 2, "neutron_star": 1, "nebula": 3, "supernova_remnant": 1,
                                                                    "hypervelocity_star": 1, "quasar": 1}
        assert all(row["uncharted"] and row["ref"] is None for row in page["rows"])
        assert [row["distance_pc"] for row in page["rows"]] == sorted(row["distance_pc"] for row in page["rows"])
        clouds = [row for row in page["rows"] if row["id"].startswith("field-")]
        assert [row["scatter_id"] for row in clouds] == [None, None] and sorted(row["id"] for row in clouds) == ["field-0.0.0.0", "field-0.0.0.1"]
        assert near.objects_within(conn, place, 50.0, kinds=["hypervelocity_star"])["total"] == 1
        assert near.objects_within(conn, place, 50.0, kinds=["nebula"])["total"] == 4
        assert near.objects_within(conn, place, 50.0, kinds=["system"])["total"] == 0
        paged = near.objects_within(conn, place, 50.0, limit=3, offset=3)
        assert len(paged["rows"]) == 3 and paged["total"] == 10
    finally:
        conn.close()


def test_the_index_run_records_how_many_uncharted_objects_each_kind_has(scattered):
    database = "test_map400_uncharted"
    conn = store.get_connection(scattered)
    try:
        objectindex.forget(database)
        objectindex.build(conn, ["nebula", "black_hole"], database=database)
        assert objectindex.recorded_uncharted(database) == {"nebula": 4, "black_hole": 2, "neutron_star": 1, "supernova_remnant": 1,
                                                            "quasar": 1, "hypervelocity_star": 1}
        assert "_uncharted" not in objectindex.recorded(database)
    finally:
        objectindex.forget(database)
        conn.close()
