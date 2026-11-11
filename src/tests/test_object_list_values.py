# tests/test_object_list_values.py

"""UX.228: the star, planet and moon lists filter and sort by the stored physical values of `db.search_values` (UX.226), the way the search page
does: a range keeps the bodies whose value lies in it, a value sort walks the pages in the order of the value (keyset paging, no `OFFSET`), each
item carries the value in the unit asked for, and a kind that takes no value (or a bad key or unit) is refused."""

import pytest

from planetgen.db import objectindex, objectlist, store
from tests.publicids import pid  # noqa: F401
from tests.test_api import client  # noqa: F401
from tests.test_object_index import galaxy  # noqa: F401  (fixture)
from tests.test_object_list import _walk, indexed  # noqa: F401  (fixtures and helper)

pytestmark = pytest.mark.fresh_database


def _stars(config):
    conn = store.get_connection(config)
    try:
        return [dict(row) for row in conn.execute("SELECT id, luminosity_w, age_gy, temperature_k, mass_kg FROM stars").fetchall()]
    finally:
        conn.close()


@pytest.mark.parametrize("key,column", [("luminosity", "luminosity_w"), ("age", "age_gy"), ("temperature", "temperature_k")])
@pytest.mark.parametrize("descending", [False, True])
def test_a_walk_sorted_by_a_value_visits_every_star_that_has_it_once_in_the_order_of_the_value(indexed, key, column, descending):
    config, _sector = indexed
    have = [row for row in _stars(config) if row[column] is not None]
    assert len(have) >= 3
    items, pages, total = _walk(config, "star", 1, sort=key, descending=descending)
    assert pages >= 3 and total == len(have)
    assert len(items) == len({item["id"] for item in items}) == len(have)
    values = [item["values"][key]["value"] for item in items]
    assert values == sorted(values, reverse=descending)
    assert all(item["values"][key]["unit"] for item in items)


def test_a_range_keeps_the_stars_in_it_and_reads_in_the_unit_it_names(indexed):
    config, _sector = indexed
    have = sorted(row["luminosity_w"] for row in _stars(config) if row["luminosity_w"] is not None)
    middle = have[len(have) // 2]
    items, _pages, total = _walk(config, "star", 2, ranges={"luminosity": (None, middle, "W")})
    wanted = [value for value in have if value <= middle]
    assert total == len(items) == len(wanted)
    assert all(item["values"]["luminosity"]["unit"] == "W" for item in items)
    again, _pages, again_total = _walk(config, "star", 2, ranges={"luminosity": (None, middle / 3.828e26, "Lsun")})
    assert again_total == total and {i["id"] for i in again} == {i["id"] for i in items}


def test_two_ranges_both_hold_and_a_range_with_nothing_in_it_lists_nothing(indexed):
    config, _sector = indexed
    rows = [row for row in _stars(config) if row["age_gy"] is not None and row["luminosity_w"] is not None]
    pick = rows[0]
    items, _pages, total = _walk(config, "star", 5, ranges={"age": (pick["age_gy"], pick["age_gy"], "Gyr"),
                                                              "luminosity": (pick["luminosity_w"] * 0.999999, None, "W")})
    assert total >= 1 and pick["id"] in {i["source_id"] for i in items}
    nothing, _pages, none_total = _walk(config, "star", 5, ranges={"age": (1e9, None, "Gyr")})
    assert nothing == [] and none_total == 0


def test_a_range_and_a_sort_page_together_by_keyset(indexed):
    config, _sector = indexed
    items, pages, total = _walk(config, "star", 1, sort="age", descending=True, ranges={"luminosity": (0.0, None, "W")})
    assert pages == total or pages == total + 1 or total == 0
    ages = [item["values"]["age"]["value"] for item in items]
    assert ages == sorted(ages, reverse=True)


@pytest.mark.parametrize("kind,key", [("planet", "gravity"), ("planet", "moon_count"), ("planet", "surface_temperature"), ("moon", "orbit_distance")])
def test_planet_and_moon_lists_take_value_sorts_and_ranges(indexed, kind, key):
    config, _sector = indexed
    items, _pages, total = _walk(config, kind, 3, sort=key)
    values = [item["values"][key]["value"] for item in items]
    assert total == len(items) and values == sorted(values)
    if items:
        least = values[0]
        kept, _pages, kept_total = _walk(config, kind, 3, ranges={key: (least, None, items[0]["values"][key]["unit"])})
        assert kept_total == len([v for v in values if v >= least])


def test_what_takes_no_value_is_refused(indexed):
    config, _sector = indexed
    conn = store.get_connection(config)
    try:
        with pytest.raises(ValueError):
            objectlist.list_objects(conn, "belt", ranges={"age": (1.0, None, "Gyr")})
        with pytest.raises(ValueError):
            objectlist.list_objects(conn, "belt", sort="age")
        with pytest.raises(ValueError):
            objectlist.list_objects(conn, "star", sort="gravity")      # a star has no surface gravity
        assert objectlist.value_sort_keys("belt") == ()
        assert "mass" not in objectlist.value_sort_keys("star") and "luminosity" in objectlist.value_sort_keys("star")
    finally:
        conn.close()


def test_the_api_takes_the_search_ranges_and_refuses_a_bad_one(indexed, client):  # noqa: F811
    config, _sector = indexed
    have = sorted(row["luminosity_w"] for row in _stars(config) if row["luminosity_w"] is not None)
    body = client.get("/api/object-lists/star", query_string={"sort": "luminosity", "order": "desc", "star_max_luminosity": have[1],
                                                              "star_luminosity_unit": "W"}).get_json()
    assert "total" in body, body
    assert body["total"] == 2 and [i["values"]["luminosity"]["unit"] for i in body["items"]] == ["W", "W"]
    assert body["items"][0]["values"]["luminosity"]["value"] >= body["items"][1]["values"]["luminosity"]["value"]
    for bad in ("star_min_luminosity=abc", "star_min_luminosity=2&star_max_luminosity=1", "star_luminosity_unit=kg&star_min_luminosity=1",
                "sort=gravity"):
        assert client.get(f"/api/object-lists/star?{bad}").status_code == 400, bad
    assert client.get("/api/object-lists/belt?star_min_luminosity=1").status_code == 200       # another kind's parameters are not this list's


def test_the_nebula_list_filters_and_sorts_by_the_stored_composition(indexed, client):  # noqa: F811
    """GEN.1276: a nebula's stored composition is searchable like a star's values."""
    config, _sector = indexed
    conn = store.get_connection(config)
    try:
        row = dict(conn.execute("SELECT id, mass_fraction_h, mass_fraction_he, dust_to_gas, gas_oh_dex FROM nebulae").fetchone())
    finally:
        conn.close()
    assert row["mass_fraction_h"] is not None and row["dust_to_gas"] is not None
    items, _pages, total = _walk(config, "nebula", 5, sort="mass_fraction_h")
    assert total == 1 and items[0]["values"]["mass_fraction_h"]["value"] == pytest.approx(row["mass_fraction_h"])
    kept, _pages, kept_total = _walk(config, "nebula", 5, ranges={"mass_fraction_h": (row["mass_fraction_h"] - 0.01, row["mass_fraction_h"] + 0.01, "frac")})
    assert kept_total == 1
    nothing, _pages, none_total = _walk(config, "nebula", 5, ranges={"mass_fraction_h": (row["mass_fraction_h"] + 0.1, None, "frac")})
    assert nothing == [] and none_total == 0
    body = client.get("/api/object-lists/nebula", query_string={"nebula_min_dust_to_gas": 0.0, "sort": "dust_to_gas"}).get_json()
    assert body["total"] == 1, body
