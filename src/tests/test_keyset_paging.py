# tests/test_keyset_paging.py

"""
PERF.75: a page that starts after the previous page's last row instead of skipping `offset` rows. For every list
that takes `after=`, the pages read by key are the pages offset paging gives, in both directions, with ties in the
sort column; a row that is gone falls back to `offset`; and a plan test shows the keyset condition uses the index.
Runs against a real throwaway database.
"""

import pytest

from planetgen.db import query as queryDb
from planetgen.db import store
from tests.test_api import client  # noqa: F401
from planetgen.galaxy.spherical import sector_address_at
from tests.test_galactic_motion import EDGE_PC, SCATTER_RADIUS_PC, _seed_scatter_grid

SORTS = list(queryDb.UNCHARTED_SYSTEM_SORTS)


@pytest.fixture
def stars(mysql_config):
    """Twenty-three scattered stars with ties in every sort column."""
    _seed_scatter_grid(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            points = [(SCATTER_RADIUS_PC * (1 + 0.4 * step), 0.0, 40.0 * step) for step in range(4)]
            for index in range(23):
                point = points[index % 4]
                shell, band, wedge = sector_address_at(point, EDGE_PC)
                conn.execute(
                    "INSERT INTO bright_stars (shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc,"
                    " position_z_mpc, population, star_type, yerkes_class, mass_kg, radius_km, temperature_k,"
                    " luminosity_w, age_gy, initial_mass_sol, seed)"
                    " VALUES (?, ?, ?, ?, ?, ?, 'young', ?, 'V', 1e31, 1e6, ?, ?, ?, 5.0, ?)",
                    (shell, band, wedge, *(c * 1000.0 for c in point), "BOA"[index % 3] + "2V",
                     10000 + 1000 * (index % 4), 1e30 * (1 + index % 6), 0.1 * (1 + index % 5), index))
    finally:
        conn.close()
    conn = queryDb.open_readonly(mysql_config)
    try:
        yield conn
    finally:
        conn.close()


def _ids(rows):
    return [row["id"] for row in rows]


@pytest.mark.parametrize("sort", SORTS)
@pytest.mark.parametrize("descending", [False, True])
def test_uncharted_pages_by_key_are_the_pages_by_offset(stars, sort, descending):
    everything = queryDb.list_uncharted_systems(stars, sort=sort, descending=descending)
    assert len(everything) == 23
    for index in range(0, len(everything), 5):
        after = everything[index - 1]["id"] if index else None
        page = queryDb.list_uncharted_systems(stars, limit=5, offset=999, sort=sort, descending=descending,
                                              after=after) if after else \
            queryDb.list_uncharted_systems(stars, limit=5, sort=sort, descending=descending)
        assert _ids(page) == _ids(everything[index:index + 5]), (sort, descending, index)


def test_a_star_that_is_gone_or_no_limit_uses_the_offset(stars):
    by_offset = queryDb.list_uncharted_systems(stars, limit=4, offset=3)
    assert _ids(queryDb.list_uncharted_systems(stars, limit=4, offset=3, after=10 ** 9)) == _ids(by_offset)
    first = queryDb.list_uncharted_systems(stars, limit=1)[0]["id"]
    assert len(queryDb.list_uncharted_systems(stars, after=first)) == 23     # no limit: every star, `after` unused


def test_the_uncharted_route_pages_by_key(client, stars):
    first = client.get("/api/uncharted-systems?limit=5").get_json()["items"]
    second = client.get("/api/uncharted-systems?limit=5&offset=5").get_json()["items"]
    keyed = client.get(f"/api/uncharted-systems?limit=5&offset=999&after={first[-1]['id']}").get_json()["items"]
    assert [item["id"] for item in keyed] == [item["id"] for item in second]
    assert client.get("/api/uncharted-systems?after=not-an-id").status_code == 400


@pytest.fixture
def many_stars(mysql_config):
    """A few thousand stars with statistics, so the optimizer plans as it would on the big table."""
    _seed_scatter_grid(mysql_config)
    point = (SCATTER_RADIUS_PC, 0.0, 0.0)
    shell, band, wedge = sector_address_at(point, EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.executemany(
                "INSERT INTO bright_stars (shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc,"
                " position_z_mpc, population, star_type, yerkes_class, mass_kg, radius_km, temperature_k,"
                " luminosity_w, age_gy, initial_mass_sol, seed) VALUES (?, ?, ?, 1, 2, 3, 'young', 'B2V', 'V',"
                " 1e31, 1e6, 20000, ?, 0.1, 5.0, ?)",
                [(shell, band, (wedge + index % 50) % 500, 1e30 * (1.0 + index / 1000.0), index) for index in range(4000)])
        conn.execute("ANALYZE TABLE bright_stars")
        conn.commit()
    finally:
        conn.close()
    conn = queryDb.open_readonly(mysql_config)
    try:
        yield conn
    finally:
        conn.close()


@pytest.mark.parametrize("descending", [False, True])
@pytest.mark.parametrize("sort", ["luminosity", "sector"])
def test_an_uncharted_page_is_read_in_index_order_with_no_sort(many_stars, sort, descending):
    """The page after a star walks an index and stops at the page (no filesort), first page or later, either way
    round (MySQL here; the same OR form is what the PERF.71 measurements timed on MariaDB 10.11)."""
    columns = queryDb.UNCHARTED_SYSTEM_SORTS[sort].split(",")
    previous = many_stars.execute(
        f"SELECT id, {', '.join(columns)} FROM bright_stars ORDER BY id LIMIT 1 OFFSET 2000").fetchone()
    condition, params = queryDb._keyset_condition(columns, [previous[column] for column in columns], descending,
                                                  "id", previous["id"])
    direction = "DESC" if descending else "ASC"
    order = f"ORDER BY {', '.join(f'{c} {direction}' for c in columns)}, id {direction} LIMIT 5"
    for where, args in ((f"star_system_id IS NULL AND {condition}", params), ("star_system_id IS NULL", [])):
        plan = many_stars.execute(f"EXPLAIN SELECT id FROM bright_stars WHERE {where} {order}", args).fetchall()
        assert plan[0]["key"] is not None, dict(plan[0])
        assert "filesort" not in (plan[0]["Extra"] or "").lower(), dict(plan[0])


def test_the_condition_is_written_with_or_and_never_as_a_row_constructor():
    condition, params = queryDb._keyset_condition(["a", "b"], [1, 2], False, "id", 9)
    assert condition == "((a > ?) OR (a = ? AND b > ?) OR (a = ? AND b = ? AND id > ?))"
    assert params == [1, 1, 2, 1, 2, 9]
    assert "(a, b" not in condition
    assert queryDb._keyset_condition(["a"], [1], True, "id", 9)[0] == "((a < ?) OR (a = ? AND id < ?))"


# --- The scattered tail of the Phenomena list ------------------------------------------------------------------

@pytest.fixture
def phenomena(mysql_config):
    """Three built supernova remnants (some leave a compact remnant of their own) and eleven unbuilt scattered black
    holes, as `list_phenomena` reads them."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
    from tests.test_uncharted_phenomenon import _scatter

    conn = store.get_connection(mysql_config)
    try:
        for _ in range(3):
            store.insert_supernova_remnant(conn, SupernovaRemnant(SystemConfig()))
        conn.commit()
    finally:
        conn.close()
    for index in range(11):
        _scatter(mysql_config, "black-hole", "stellar", seed=index + 1, mass=10.0 + index, record=index == 0)
    conn = store.get_connection(mysql_config)
    try:
        store.record_phenomenon_scatter_classes(conn, {("black-hole", "stellar"): 11})
        conn.commit()
    finally:
        conn.close()
    conn = queryDb.open_readonly(mysql_config)
    try:
        yield conn
    finally:
        conn.close()


@pytest.mark.parametrize("descending", [False, True])
def test_phenomena_pages_inside_the_scattered_tail_by_key_are_the_pages_by_offset(phenomena, descending):
    everything = queryDb.list_phenomena(phenomena, sort="name", descending=descending)
    built = [row for row in everything if not row["scattered"]]
    assert len(built) >= 3 and [row["scattered"] for row in everything] == [False] * len(built) + [True] * 11
    for index in range(len(built) + 1, len(everything), 4):
        after = everything[index - 1]["scatter_id"]
        page = queryDb.list_phenomena(phenomena, limit=4, offset=999, sort="name", descending=descending,
                                      after=after)
        assert [row["id"] for row in page] == [row["id"] for row in everything[index:index + 4]], index
    keyed = queryDb.list_phenomena(phenomena, limit=4, after=everything[len(built)]["scatter_id"])
    assert [row["scatter_id"] for row in keyed] == [row["scatter_id"] for row in everything[len(built) + 1:][:4]]


def test_the_phenomena_route_pages_the_scattered_tail_by_key(client, phenomena):
    everything = client.get("/api/phenomena?limit=40").get_json()["items"]
    first = next(index for index, item in enumerate(everything) if item["scattered"])
    keyed = client.get(f"/api/phenomena?limit=4&offset=999&after={everything[first + 1]['scatter_id']}").get_json()
    assert [item["id"] for item in keyed["items"]] == [item["id"] for item in everything[first + 2:first + 6]]
    assert client.get("/api/phenomena?after=not-an-id").status_code == 400


# --- The /table routes hand the key on and give the next one ----------------------------------------------------

@pytest.fixture
def site(monkeypatch):
    from planetgen.web.app import create_app
    from planetgen.web.lib import apiclient
    from tests.test_web_pages import FakeData, _FakeConfig

    data = FakeData()
    for name in ("get_sectors", "get_systems", "get_uncharted_systems", "get_phenomena", "auth_me"):
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    app = create_app(_FakeConfig)
    app.testing = True
    return app.test_client(), monkeypatch, apiclient


def _star(star_id):
    return {"id": star_id, "star_type": "B2V", "population": "young", "luminosity_sol": 1000.0,
            "temperature_k": 20000.0, "age_gy": 0.1, "designation": "7.2.3", "shell_index": 7, "band_index": 2,
            "wedge_index": 3, "x": 1.0, "y": 2.0, "z": 3.0, "local_x": 0.1, "local_y": 0.2, "local_z": 0.3,
            "member_count": 0, "system_mass_solar": None}


def test_the_uncharted_table_passes_the_key_on_and_gives_the_last_star_as_next(site):
    client, patch, apiclient = site
    seen = {}

    def fake(db, limit=None, offset=None, sort=None, descending=False, after=None, start=None):
        seen.update(after=after, offset=offset)
        return {"items": [_star(41), _star(42)], "total": 99}

    patch.setattr(apiclient, "get_uncharted_systems", fake)
    body = client.get("/table/uncharted_systems?after=40&offset=10&limit=2").get_json()
    assert seen == {"after": "40", "offset": 10} and body["next"] == 42 and len(body["rows"]) == 2
    client.get("/table/uncharted_systems?limit=2")
    assert seen["after"] is None


def _phenomenon(name, scattered, scatter_id=None):
    return {"id": 5, "scatter_id": scatter_id, "type": "black_hole", "name": name, "descriptor": "stellar",
            "radius_ly": 0.0, "sector_id": None, "sector_name": None, "placed": True, "scattered": scattered}


def test_the_phenomena_table_gives_a_next_key_only_when_a_page_ends_in_the_scattered_tail(site):
    client, patch, apiclient = site
    seen = {}

    def fake(db, limit=None, offset=None, sort=None, descending=False, types=(), descriptors=(), placed=None,
             facets=False, after=None, start=None, rows=None):
        seen["after"] = after
        return {"items": fake.rows, "total": 99}

    patch.setattr(apiclient, "get_phenomena", fake)
    fake.rows = [_phenomenon("Built", False), _phenomenon("Uncharted black hole 1", True, 77)]
    assert client.get("/table/phenomena?after=70&limit=2").get_json()["next"] == 77 and seen["after"] == "70"
    fake.rows = [_phenomenon("Built", False)]
    assert client.get("/table/phenomena?limit=2").get_json()["next"] is None
