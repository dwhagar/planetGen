# tests/test_object_list.py

"""UX.141: the lists of every kind read the object index with keyset paging (`db/objectlist.py`, `GET /api/object-lists/<kind>`):
a walk of the pages visits every object once in the order of the sort, in both directions, a filter leaves only what it
names, nothing uses `OFFSET`, and a bad kind, sort or cursor is refused."""

import pytest

from planetgen.db import objectindex, objectlist, store
from tests.publicids import pid  # noqa: F401
from tests.test_api import client  # noqa: F401
from tests.test_object_index import galaxy  # noqa: F401  (fixture)

pytestmark = pytest.mark.fresh_database


@pytest.fixture
def indexed(galaxy):  # noqa: F811
    config, sector_id = galaxy
    conn = store.get_connection(config)
    try:
        objectindex.build(conn)
        conn.commit()
    finally:
        conn.close()
    return config, sector_id


def _walk(config, kind, limit, **options):
    """Every item of the pages from the top, and the number of pages."""
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
            assert pages < 1000
    finally:
        conn.close()


def _all(config, kind):
    conn = store.get_connection(config)
    try:
        return [dict(row) for row in conn.execute("SELECT * FROM object_index WHERE kind = ?", (kind,)).fetchall()]
    finally:
        conn.close()


@pytest.mark.parametrize("sort", sorted(objectlist.SORTS))
@pytest.mark.parametrize("descending", [False, True])
def test_a_walk_of_the_pages_visits_every_object_once_in_the_order_of_the_sort(indexed, sort, descending):
    config, _sector = indexed
    # The fixture's three systems are single stars with a drawn number of planets, so the stars are the kind that is
    # always there (TEST.267: a draw of too few planets made a walk of planets fail one run in four).
    kind = "star"
    everything = _all(config, kind)
    column = objectlist.SORTS[sort][0]
    wanted = [row for row in everything if row[column] is not None]
    assert len(wanted) >= 3, "the galaxy has too few stars for a page walk"
    items, pages, total = _walk(config, kind, 1, sort=sort, descending=descending)
    assert pages >= 3 and total == len(wanted)
    assert len(items) == len({item["id"] for item in items}) == len(wanted)
    key = {"name": "name", "class": "object_class", "mass": "mass_kg", "radius": "radius_km"}[sort]
    values = [item[key].lower() if isinstance(item[key], str) else item[key] for item in items]
    assert values == sorted(values, reverse=descending)


def test_the_filters_leave_only_what_they_name(indexed):
    config, sector_id = indexed
    stars = _all(config, "star")
    prefix = stars[0]["name"][:3]
    items, _pages, total = _walk(config, "star", 3, name=prefix.lower())
    assert items and total == len(items)
    assert all(item["name"].lower().startswith(prefix.lower()) for item in items)
    assert len(items) == len([s for s in stars if s["name"].lower().startswith(prefix.lower())])
    klass = stars[0]["object_class"]
    items, _pages, _total = _walk(config, "star", 3, object_class=klass)
    assert items and all(item["object_class"] == klass for item in items)
    items, _pages, _total = _walk(config, "star", 4, sector_id=sector_id)
    assert items and all(item["sector_id"] == sector_id for item in items)
    assert _walk(config, "star", 4, sector_id=sector_id + 999)[0] == []
    assert _walk(config, "star", 4, name="%")[0] == [], "a percent sign is a character, not a wildcard"


def test_an_item_carries_the_printed_id_its_system_and_its_parent(indexed):
    config, _sector = indexed
    items, _pages, _total = _walk(config, "planet", 50)
    assert items
    item = items[0]
    assert item["kind"] == "planet" and len(item["id"].split("-")) == 3 and item["system_id"] is not None
    assert item["parent_id"] is not None and item["source_id"] > 0


def test_a_planet_or_moon_item_carries_its_stored_habitability(indexed):
    """UX.100: `phi4` and the six `tier_<domain>` ratings come from the kind's own table; other kinds carry none."""
    config, _sector = indexed
    columns = objectlist.HABITABILITY_COLUMNS
    for kind, table in (("planet", "planets"), ("moon", "moons")):
        items, _pages, _total = _walk(config, kind, 50)
        conn = store.get_connection(config)
        try:
            stored = {row["id"]: row for row in conn.execute(f"SELECT id, {', '.join(columns)} FROM {table}").fetchall()}
        finally:
            conn.close()
        for item in items:
            assert {column: item[column] for column in columns} == {column: stored[item["source_id"]][column]
                                                                   for column in columns}
    stars, _pages, _total = _walk(config, "star", 50)
    assert stars and all("phi4" not in star for star in stars)


def test_a_star_planet_or_moon_item_carries_its_stored_physics_values(indexed):
    """UX.187: a star's temperature, luminosity and age and a planet's or moon's orbit size, eccentricity and surface
    temperature come from the kind's own table, as stored; the other kinds carry none."""
    config, _sector = indexed
    for kind, table in (("star", "stars"), ("planet", "planets"), ("moon", "moons")):
        columns = objectlist.PHYSICS_COLUMNS[kind]
        items, _pages, _total = _walk(config, kind, 50)
        conn = store.get_connection(config)
        try:
            stored = {row["id"]: row for row in conn.execute(f"SELECT id, {', '.join(columns)} FROM {table}").fetchall()}
        finally:
            conn.close()
        assert items
        for item in items:
            assert {column: item[column] for column in columns} == {column: stored[item["source_id"]][column]
                                                                   for column in columns}
    belts, _pages, _total = _walk(config, "belt", 50)
    assert all("temperature_k" not in belt and "eccentricity" not in belt for belt in belts)


def test_a_cursor_for_a_gone_object_starts_from_the_top_and_a_bad_one_is_refused(indexed):
    config, _sector = indexed
    conn = store.get_connection(config)
    try:
        top = objectlist.list_objects(conn, "star", limit=2)
        gone = objectlist.list_objects(conn, "star", limit=2, after="00000000000000-0000000-000")
        assert gone["items"] == top["items"]
        with pytest.raises(ValueError):
            objectlist.list_objects(conn, "star", after="not an id")
        with pytest.raises(ValueError):
            objectlist.list_objects(conn, "galaxy")
        with pytest.raises(ValueError):
            objectlist.list_objects(conn, "star", sort="age")
        with pytest.raises(ValueError):
            objectlist.list_objects(conn, "star", limit=0)
        assert len(objectlist.list_objects(conn, "star", limit=10 ** 6)["items"]) <= objectlist.MAX_LIMIT
    finally:
        conn.close()


def test_a_page_reads_the_index_without_offset_or_a_filesort(indexed):
    config, _sector = indexed
    conn = store.get_connection(config)
    try:
        for sort, index in (("name", "idx_object_index_name"), ("mass", "idx_object_index_mass"),
                            ("radius", "idx_object_index_radius"), ("class", "idx_object_index_class")):
            plan = conn.execute(
                f"EXPLAIN SELECT * FROM object_index WHERE kind = 'star' AND {objectlist.SORTS[sort][0]} IS NOT NULL "
                f"ORDER BY {', '.join(objectlist.SORTS[sort])}, uid LIMIT 10").fetchall()
            assert "filesort" not in (plan[0]["Extra"] or "").lower(), (sort, plan)
            assert plan[0]["key"] == index, (sort, plan[0]["key"])
    finally:
        conn.close()
    source = open(objectlist.__file__).read()
    assert "OFFSET" not in source.upper().replace("`OFFSET`", "")


def test_the_api_pages_a_kind_and_refuses_what_it_cannot_serve(client, indexed):  # noqa: F811
    first = client.get("/api/object-lists/star?limit=2&sort=mass&order=desc").get_json()
    assert len(first["items"]) == 2 and first["next"] and first["total"] >= 3 and first["total_capped"] is False
    second = client.get(f"/api/object-lists/star?limit=2&sort=mass&order=desc&after={first['next']}").get_json()
    assert second["items"][0]["mass_kg"] <= first["items"][-1]["mass_kg"]
    assert {i["id"] for i in first["items"]}.isdisjoint(i["id"] for i in second["items"])
    for bad in ("/api/object-lists/galaxy", "/api/object-lists/star?sort=age", "/api/object-lists/star?order=up",
                "/api/object-lists/star?limit=x", "/api/object-lists/star?after=zz", "/api/object-lists/star?sector_id=a"):
        assert client.get(bad).status_code == 400, bad


def test_a_start_begins_the_name_sorted_list_at_the_first_name_from_it_on(indexed):
    config, _sector = indexed
    conn = store.get_connection(config)
    try:
        names = [item["name"] for item in objectlist.list_objects(conn, "planet", limit=200)["items"]]
        assert len(names) >= 4
        cut = names[len(names) // 2]
        page = objectlist.list_objects(conn, "planet", limit=200, start=cut)
        assert [item["name"] for item in page["items"]] == [n for n in names if n.lower() >= cut.lower()]
        assert page["past_end"] is False and page["total"] == len(names)
        down = objectlist.list_objects(conn, "planet", limit=200, descending=True, start=cut)
        assert [item["name"] for item in down["items"]] == [n for n in reversed(names) if n.lower() <= cut.lower()]
        past = objectlist.list_objects(conn, "planet", limit=2, start="zzzzzzzz")
        assert past["past_end"] is True and [item["name"] for item in past["items"]] == names[-2:] and past["next"] is None
        with pytest.raises(ValueError):
            objectlist.list_objects(conn, "planet", sort="mass", start="a")
    finally:
        conn.close()


def test_the_route_takes_from_and_refuses_it_under_another_sort(indexed, client):  # noqa: F811
    assert client.get("/api/object-lists/planet?from=a&limit=2").status_code == 200
    assert client.get("/api/object-lists/planet?from=a&sort=mass").status_code == 400
