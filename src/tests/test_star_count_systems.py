# tests/test_star_count_systems.py

"""UX.239, the later slices: `GET /api/systems?stars=` and the `stars` sort, and `db.query.list_systems` and `count_systems` under them."""

from planetgen.db import query
from tests.test_api import client, seeded_sector  # noqa: F401
from tests.test_star_count_filter import multiples  # noqa: F401


def _names(response):
    assert response.status_code == 200, response.get_data(as_text=True)
    return [row["name"] for row in response.get_json()["items"]]


def test_the_systems_list_filters_by_stars_and_counts_what_it_filters(multiples):  # noqa: F811
    pairs = _names(multiples.get("/api/systems?stars=2"))
    assert "Pair Kappa" in pairs and "Triple Kappa" not in pairs
    body = multiples.get("/api/systems?stars=3%2B").get_json()
    assert {"Triple Kappa", "Quad Kappa"} <= {row["name"] for row in body["items"]} and all(row["star_count"] >= 3 for row in body["items"])
    assert body["total"] == len(body["items"])
    assert not {"Solo Kappa", "Pair Kappa"} & set(_names(multiples.get("/api/systems?stars=trinary")))
    bad = multiples.get("/api/systems?stars=lots")
    assert bad.status_code == 400 and "stars must be" in bad.get_json()["error"]


def test_the_systems_list_sorts_by_the_number_of_stars(multiples):  # noqa: F811
    ascending = multiples.get("/api/systems?sort=stars&limit=500").get_json()["items"]
    counts = [row["star_count"] for row in ascending]
    assert counts == sorted(counts) and len(set(counts)) > 1
    descending = multiples.get("/api/systems?sort=stars&order=desc&limit=500").get_json()["items"]
    assert [row["star_count"] for row in descending] == sorted(counts, reverse=True)
    assert multiples.get("/api/systems?sort=nonsense").status_code == 400


def test_the_stars_filter_and_the_sort_work_together_with_the_other_filters(multiples):  # noqa: F811
    rows = multiples.get("/api/systems?stars=2%2B&sort=stars&order=desc&limit=500&facets=1").get_json()
    counts = [row["star_count"] for row in rows["items"]]
    assert counts and all(count >= 2 for count in counts) and counts == sorted(counts, reverse=True)
    assert rows["facets"]["binary"], "the facet menus answer with the filter on"


def test_list_systems_in_the_database_layer_agrees_with_the_api(multiples, seeded_sector):  # noqa: F811
    mysql_config, _sector_id, _systems = seeded_sector
    conn = query.open_readonly(mysql_config)
    try:
        listed = {row["name"] for row in query.list_systems(conn, stars=(3, 3))}
        assert "Triple Kappa" in listed and "Pair Kappa" not in listed
        assert query.count_systems(conn, stars=(3, 3)) == len(listed)
        ordered = [row["star_count"] for row in query.list_systems(conn, sort="stars", descending=True)]
        assert ordered == sorted(ordered, reverse=True)
    finally:
        conn.close()
