# tests/test_systems_values.py

"""UX.228: the Systems list filters by what a system holds: the range filters of the search (`star_`, `planet_`, `moon_`) keep the systems with at least one
body of that kind whose stored value lies in the ranges, `db.query.list_systems`, `count_systems` and `systems_facets` agree, and `GET /api/systems` takes them."""

import pytest

from planetgen.db import query, store
from tests.publicids import pid  # noqa: F401
from tests.test_api import client  # noqa: F401
from tests.test_object_index import galaxy  # noqa: F401  (fixture)

pytestmark = pytest.mark.fresh_database


def _rows(config, sql, params=()):
    conn = store.get_connection(config)
    try:
        return [dict(row) for row in conn.execute(sql, params).fetchall()]
    finally:
        conn.close()


def _systems_with(config, table, condition, params=()):
    return {row["star_system_id"] for row in _rows(config, f"SELECT DISTINCT star_system_id FROM {table} WHERE {condition}", params)}


def test_a_star_range_keeps_the_systems_with_a_star_in_it_and_the_count_agrees(galaxy):  # noqa: F811
    config, _sector = galaxy
    lums = sorted(row["luminosity_w"] for row in _rows(config, "SELECT luminosity_w FROM stars WHERE luminosity_w IS NOT NULL"))
    assert len(lums) >= 3
    cut = lums[len(lums) // 2]
    expected = _systems_with(config, "stars", "luminosity_w <= ?", (cut,))
    assert expected and len(expected) < len(_rows(config, "SELECT id FROM star_systems"))
    values = {"star": {"luminosity": (None, cut, "W")}}
    conn = query.open_readonly(config)
    try:
        listed = query.list_systems(conn, values=values)
        assert {row["id"] for row in listed} == expected
        assert query.count_systems(conn, values=values) == len(expected)
        facets = query.systems_facets(conn, values=values)
        assert sum(option["count"] for option in facets["placement"]) == len(expected)
    finally:
        conn.close()


def test_planet_and_moon_ranges_count_a_system_once_and_two_kinds_both_have_to_hold(galaxy):  # noqa: F811
    config, _sector = galaxy
    planets = _systems_with(config, "planets", "gravity_g IS NOT NULL AND gravity_g >= 0")
    stars = _systems_with(config, "stars", "age_gy IS NOT NULL AND age_gy >= 0")
    conn = query.open_readonly(config)
    try:
        one = query.list_systems(conn, values={"planet": {"gravity": (0.0, None, "g0")}})
        assert [row["id"] for row in one] == sorted({row["id"] for row in one}, key=lambda i: [r["id"] for r in one].index(i))   # no system twice
        assert {row["id"] for row in one} == planets
        both = query.list_systems(conn, values={"planet": {"gravity": (0.0, None, "g0")}, "star": {"age": (0.0, None, "Gyr")}})
        assert {row["id"] for row in both} == planets & stars
        assert query.count_systems(conn, values={"planet": {"gravity": (0.0, None, "g0")}, "star": {"age": (0.0, None, "Gyr")}}) == len(planets & stars)
        nobody = query.list_systems(conn, values={"star": {"age": (1e9, None, "Gyr")}})
        assert nobody == []
        together = query.list_systems(conn, values={"star": {"age": (0.0, None, "Gyr")}}, stars=(1, 1), sort="stars")
        assert all(row["star_count"] == 1 for row in together)
    finally:
        conn.close()


def test_the_api_takes_the_search_ranges_for_systems_and_refuses_a_bad_one(galaxy, client):  # noqa: F811
    config, _sector = galaxy
    lums = sorted(row["luminosity_w"] for row in _rows(config, "SELECT luminosity_w FROM stars WHERE luminosity_w IS NOT NULL"))
    cut = lums[len(lums) // 2]
    expected = _systems_with(config, "stars", "luminosity_w <= ?", (cut,))
    body = client.get("/api/systems", query_string={"star_max_luminosity": cut, "star_luminosity_unit": "W", "facets": "1", "limit": 200}).get_json()
    assert body["total"] == len(expected) == len(body["items"])
    assert body["facets"]["placement"]
    for bad in ("star_min_luminosity=abc", "star_min_age=2&star_max_age=1", "star_age_unit=kg&star_min_age=1", "planet_min_age=1"):
        response = client.get(f"/api/systems?{bad}")
        if bad.startswith("planet_min_age"):
            assert response.status_code == 200 and response.get_json()["total"] > 0       # a key the kind lacks is not one of its parameters
        else:
            assert response.status_code == 400, bad
