# tests/test_near_values.py

"""
UX.228: What's nearby takes the search's range filters and an order by a stored physical value (`planetgen.db.near.objects_within`
`ranges` and `value_sort`, `GET /api/near` `planet_min_gravity` ..., `planet_sort`): the bodies it lists are exactly those inside the sphere whose
stored value passes, each row carries the value in the unit asked for, an order by a value lists that kind's bodies by it (nearest first among
equals) in pages that add up, and a kind or key that cannot be ranged is refused.
"""

import pytest

from planetgen.db import near, query
from tests.publicids import pid
from tests.test_near import ADDRESSES, EDGE_PC, _conn, _everything, _systems, client, galaxy  # noqa: F401  (fixtures and helpers)
from planetgen.galaxy.spherical import sector_position_pc


def _setup(galaxy_config):
    conn = _conn(galaxy_config)
    center = sector_position_pc(*ADDRESSES[0], EDGE_PC)
    return conn, near.place_from_point(center)


def _stars(conn, ids):
    marks = ",".join(map(str, ids or [0]))
    return {row["id"]: dict(row) for row in conn.execute(
        f"SELECT id, luminosity_w, age_gy, temperature_k, mass_kg FROM stars WHERE star_system_id IN ({marks})").fetchall()}


def _system_ids(rows):
    return {int(row["ref"].split(":")[1]) for row in rows if row["kind"] == "system"}


def test_a_range_on_a_star_value_lists_only_the_stars_that_pass_it_and_says_the_value(galaxy):
    conn, place = _setup(galaxy)
    try:
        _page, everything = _everything(conn, place, 10.0, kinds=["system", "star"])
        stars = _stars(conn, _system_ids(everything))
        assert len(stars) >= 3
        middle = sorted(s["luminosity_w"] for s in stars.values() if s["luminosity_w"] is not None)[len(stars) // 2]
        page, rows = _everything(conn, place, 10.0, ranges={"star": {"luminosity": (None, middle, "W")}})
        assert {row["kind"] for row in rows} == {"star"}                         # a ranged kind alone when `kinds` is not given
        assert {row["id"] for row in rows} == {i for i, s in stars.items() if s["luminosity_w"] is not None and s["luminosity_w"] <= middle}
        assert page["total"] == len(rows) == page["by_kind"]["star"]
        for row in rows:
            assert row["values"]["luminosity"]["unit"] == "W"
            assert row["values"]["luminosity"]["value"] == pytest.approx(stars[row["id"]]["luminosity_w"])
        away = [row["distance_pc"] for row in rows]
        assert away == sorted(away)
        same, again = _everything(conn, place, 10.0, ranges={"star": {"luminosity": (None, middle / 3.828e26, "Lsun")}})
        assert [row["ref"] for row in again] == [row["ref"] for row in rows]       # the unit changes what is read, not what is kept
    finally:
        conn.close()


def test_kinds_asked_for_stay_and_the_ranged_kind_is_added(galaxy):
    conn, place = _setup(galaxy)
    try:
        _page, everything = _everything(conn, place, 10.0, kinds=["system", "star"])
        stars = _stars(conn, _system_ids(everything))
        assert stars
        nothing = {"star": {"age": (1e9, None, "Gyr")}}                             # no star is that old
        page, rows = _everything(conn, place, 10.0, kinds=["system"], ranges=nothing)
        assert {row["kind"] for row in rows} == {"system"} and not any(row["kind"] == "star" for row in rows) and page["by_kind"].get("star", 0) == 0
        assert _system_ids(rows) == _system_ids(everything)                         # every system stays; the stars were filtered to none
    finally:
        conn.close()


@pytest.mark.parametrize("descending", [False, True])
def test_an_order_by_a_value_lists_the_bodies_that_have_it_in_the_order_of_the_value_and_pages_add_up(galaxy, descending):
    conn, place = _setup(galaxy)
    try:
        _page, everything = _everything(conn, place, 10.0, kinds=["system", "star"])
        stars = _stars(conn, _system_ids(everything))
        page, rows = _everything(conn, place, 10.0, value_sort=("star", "age", descending))
        have = {i for i, s in stars.items() if s["age_gy"] is not None}
        assert {row["id"] for row in rows} == have and page["total"] == len(have) == page["by_kind"]["star"]
        ages = [row["values"]["age"]["value"] for row in rows]
        assert ages == sorted(ages, reverse=descending)
        assert all(row["parent"]["kind"] == "system" and row["parent"]["name"] for row in rows)
        paged = []
        for offset in range(0, page["total"], 2):
            paged += near.objects_within(conn, place, 10.0, limit=2, offset=offset, value_sort=("star", "age", descending))["rows"]
        assert [row["ref"] for row in paged] == [row["ref"] for row in rows]
    finally:
        conn.close()


@pytest.mark.parametrize("entity,key", [("planet", "gravity"), ("planet", "moon_count"), ("moon", "orbit_distance")])
def test_planets_and_moons_take_orders_and_ranges_too(galaxy, entity, key):
    conn, place = _setup(galaxy)
    try:
        page, rows = _everything(conn, place, 10.0, value_sort=(entity, key, False))
        values = [row["values"][key]["value"] for row in rows]
        assert values == sorted(values) and page["total"] == len(rows)
        if entity == "moon":
            assert all(row["parent"]["kind"] == "planet" for row in rows)
        if rows:
            least = values[len(values) // 2]
            ranged, kept = _everything(conn, place, 10.0, ranges={entity: {key: (least, None, rows[0]["values"][key]["unit"])}})
            assert len(kept) == len([v for v in values if v >= least])
    finally:
        conn.close()


def test_what_cannot_be_ranged_or_ordered_is_refused(galaxy):
    conn, place = _setup(galaxy)
    try:
        for kwargs in ({"ranges": {"belt": {"age": (1.0, None, "Gyr")}}}, {"ranges": {"star": {"gravity": (1.0, None, "g0")}}},
                       {"value_sort": ("belt", "age", False)}, {"value_sort": ("star", "gravity", False)},
                       {"value_sort": ("star", "age", False), "ranges": {"planet": {"gravity": (1.0, None, "g0")}}}):
            with pytest.raises(near.NearError):
                near.objects_within(conn, place, 5.0, **kwargs)
    finally:
        conn.close()


def test_nearest_takes_the_ranges(galaxy):
    conn, place = _setup(galaxy)
    try:
        result = near.nearest(conn, place, 2, ranges={"star": {"age": (0.0, None, "Gyr")}})
        assert result["rows"] and all(row["kind"] == "star" for row in result["rows"])
    finally:
        conn.close()


def test_the_api_takes_the_search_parameters_and_refuses_what_does_not_fit(client, galaxy):
    conn = _conn(galaxy)
    try:
        system_id = next(iter(_systems(conn)))
    finally:
        conn.close()
    base = f"/api/near?from=system:{pid('system', system_id)}&distance=10"
    body = client.get(base + "&star_max_age=1000&star_age_unit=Gyr&limit=200").get_json()
    assert body["rows"] and all(row["kind"] == "star" and row["values"]["age"]["unit"] == "Gyr" for row in body["rows"])
    ordered = client.get(base + "&star_sort=-age&limit=200").get_json()
    ages = [row["values"]["age"]["value"] for row in ordered["rows"]]
    assert ages == sorted(ages, reverse=True) and {row["kind"] for row in ordered["rows"]} == {"star"}
    assert client.get(base.replace("distance=10", "nearest=2") + "&star_max_age=1000").status_code == 200
    for bad in ("&star_min_age=abc", "&star_min_age=2&star_max_age=1", "&star_sort=gravity", "&star_sort=age&planet_sort=gravity",
                "&nearest=2&star_sort=age", "&star_age_unit=kg&star_min_age=1"):
        assert client.get(base + bad).status_code == 400, bad
