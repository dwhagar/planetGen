# tests/test_virtual_uncharted_cells.py

"""
OPS.154 (read side): an uncharted cell that holds a waiting object and has no chart-ledger entry needs no sector row of its own. The Sectors list, its counts and the sector
page take such cells from the block light parts by address, after the rows, by distance; a cell that has a row (a rogue planet's, or one the scatters made before the
migration that deletes them) is listed once, as its row, and the bulk rows still there keep the list as it was.
"""

import pytest

from planetgen.api.config import Config
from planetgen.db import chart_ledger, query, store
from planetgen.galaxy import uid as galaxy_uid
from planetgen.galaxy.spherical import provisional_sector_designation
from planetgen.web.app import create_app

pytestmark = pytest.mark.db

EDGE_PC = 4.0
CHARTED = [(9, 2, 5)]
ROGUE_ROW = (10, 3, 4)       # a cell a rogue planet made a row in
WAITING = [(9, 2, 7), (9, 2, 8), (10, 3, 3), (10, 3, 4), (11, 1, 2), (9, 2, 5)]      # cells with a part: one is charted, one has a row


def _part(conn, address):
    conn.execute("INSERT INTO galaxy_block_light_parts (shell_index, band_index, wedge_index, star_count, star_luminosity_sol, star_lpx, star_lpy, star_lpz, star_lt,"
                 " point_count, point_luminosity_sol, point_lpx, point_lpy, point_lpz) VALUES (?, ?, ?, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)", address)


@pytest.fixture
def galaxy(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            for address in CHARTED:
                sector_id = store.uncharted_sector_id(conn, address, EDGE_PC)
                chart_ledger.record(conn, address, store.sector_run_id(conn), 1)
                conn.execute("UPDATE sectors SET name = 'Charted' WHERE id = ?", (sector_id,))
            store.uncharted_sector_id(conn, ROGUE_ROW, EDGE_PC)
            for address in WAITING:
                _part(conn, address)
    finally:
        conn.close()
    return mysql_config


@pytest.fixture
def client(galaxy):
    class TestConfig(Config):
        MYSQL_CONFIG = galaxy
        WRITE_MYSQL_CONFIG = galaxy
        CONTROL_MYSQL_CONFIG = galaxy

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def _uid(address):
    return galaxy_uid.format_sector_uid(galaxy_uid.sector_uid(*address))


def _uncharted(client, query_string=""):
    return client.get(f"/api/sectors?coverage=uncharted&empty=1&{query_string}").get_json()


def test_the_uncharted_list_is_the_rows_then_the_cells_without_a_row_by_distance(client):
    body = _uncharted(client, "limit=50")
    cells = [(i["shell_index"], i["band_index"], i["wedge_index"]) for i in body["items"]]
    # the rogue's row first, then the cells that hold something and have no row, shell then band then wedge; the charted cell is not uncharted
    assert cells == [ROGUE_ROW, (9, 2, 7), (9, 2, 8), (10, 3, 3), (11, 1, 2)]
    virtual = [i for i in body["items"] if i.get("virtual")]
    assert [i["id"] for i in virtual] == [_uid(a) for a in cells[1:]]
    assert all(i["uncharted"] and i["system_count"] == 0 and i["placed"] and i["designation"] == provisional_sector_designation(i["shell_index"], i["band_index"], i["wedge_index"])
               for i in virtual)
    assert body["total"] == 5 and body["uncharted_total"] == 5


def test_the_charted_list_does_not_change(client):
    body = client.get("/api/sectors?limit=50").get_json()
    assert [i["name"] for i in body["items"]] == ["Charted"] and not any(i.get("virtual") for i in body["items"])
    assert body["charted_total"] == 1


def test_pages_follow_one_another_by_the_id_of_the_last_item_or_by_offset(client):
    first = _uncharted(client, "limit=2")
    assert len(first["items"]) == 2 and not first["items"][0].get("virtual") and first["items"][1]["virtual"] is True      # the rogue's row, then the first cell without one
    by_after, by_offset, after = [], [], None
    page = first
    seen = [i["id"] for i in first["items"]]
    while True:
        after = page["items"][-1]["id"]
        page = _uncharted(client, f"limit=2&after={after}")
        if not page["items"]:
            break
        seen.extend(i["id"] for i in page["items"])
    for offset in range(0, 6, 2):
        by_offset.extend(i["id"] for i in _uncharted(client, f"limit=2&offset={offset}")["items"])
    assert len(seen) == 5 == len(set(seen)) and by_offset == seen


def test_a_descending_list_walks_the_cells_backward(client):
    body = _uncharted(client, "limit=50&sort=distance&order=desc")
    assert [i["virtual"] for i in body["items"] if i.get("virtual")] and [(i["shell_index"], i["band_index"], i["wedge_index"]) for i in body["items"] if i.get("virtual")] == [
        (11, 1, 2), (10, 3, 3), (9, 2, 8), (9, 2, 7)]


def test_an_after_that_names_no_row_and_no_waiting_cell_is_a_400(client):
    assert client.get(f"/api/sectors?coverage=uncharted&after={_uid((12, 0, 1))}").status_code == 400


def test_a_cell_without_a_row_answers_by_its_id(client):
    detail = client.get(f"/api/sectors/{_uid((9, 2, 7))}").get_json()
    assert detail["virtual"] and detail["uncharted"] and detail["id"] == _uid((9, 2, 7)) and detail["systems"] == [] and detail["shell_index"] == 9
    assert client.get(f"/api/sectors/{_uid((12, 0, 1))}").status_code == 404      # nothing waiting there, no row: still no sector
    assert client.get(f"/api/sectors/{_uid((9, 2, 5))}").status_code != 404      # a charted cell has its row


def test_the_counts_include_the_cells(galaxy):
    conn = store.get_connection(galaxy)
    try:
        assert query.count_virtual_cells(conn) == 4
        assert query.count_uncharted_sectors(conn) == 5
        assert query.count_sectors(conn, coverage="uncharted", show_empty=True) == 5
        assert query.count_sectors(conn, coverage="uncharted", show_empty=True, virtual=False) == 1
        assert query.count_sectors(conn, coverage="charted") == 1
    finally:
        conn.close()


def test_a_quadrant_filter_keeps_the_cells_in_those_quadrants(galaxy):
    from planetgen.galaxy import spherical
    conn = store.get_connection(galaxy)
    try:
        everything = query.virtual_cells(conn, 50)
        for quadrant in query.SECTOR_QUADRANTS:
            found = query.virtual_cells(conn, 50, quadrants=(quadrant,))
            expected = [i for i in everything if query.SECTOR_QUADRANTS[spherical.sector_quadrant(i["center_x_pc"], i["center_y_pc"]) - 1] == quadrant]
            assert [i["id"] for i in found] == [i["id"] for i in expected]
    finally:
        conn.close()


def test_the_bulk_rows_still_there_keep_the_list_as_it_was(galaxy, monkeypatch):
    """Over the guard the cells are not listed from the parts (they have rows), so the list is the rows alone."""
    monkeypatch.setattr(query, "VIRTUAL_PLACEHOLDER_LIMIT", 0)
    conn = store.get_connection(galaxy)
    try:
        assert not query.virtual_cells_listed(conn)
        assert [i["name"] for i in query.list_sectors(conn, limit=50, coverage="uncharted", show_empty=True)] == [provisional_sector_designation(*ROGUE_ROW)]
        assert query.count_sectors(conn, coverage="uncharted", show_empty=True) == 1
    finally:
        conn.close()


def test_a_cell_that_has_a_row_is_listed_once(galaxy):
    conn = store.get_connection(galaxy)
    try:
        with conn:
            store.uncharted_sector_id(conn, (9, 2, 7), EDGE_PC)      # a rogue lands in a waiting cell: now a row
        items = query.list_sectors(conn, limit=50, coverage="uncharted", show_empty=True)
        cells = [(i["shell_index"], i["band_index"], i["wedge_index"]) for i in items]
        assert len(cells) == len(set(cells)) == 5 - 0 and (9, 2, 7) in cells
        assert [i for i in items if (i["shell_index"], i["band_index"], i["wedge_index"]) == (9, 2, 7)][0].get("virtual") is None
    finally:
        conn.close()


def _search(conn, term, scope=None):
    options = {} if scope is None else {"sector_scope": scope}
    return query.search(conn, {"sector_q": term}, {facet: set() for facet in query.SEARCH_TAG_FACETS}, **options)["results"]["sectors"]


def test_a_designation_typed_whole_finds_the_cell_that_has_no_row(galaxy):
    """OPS.154 part 2: search by name reaches a cell with no row by its designation (with uncharted sectors asked for), and only then."""
    conn = store.get_connection(galaxy)
    try:
        address = WAITING[0]
        found = _search(conn, _uid(address), ("uncharted", True))
        assert [row["id"] for row in found["rows"]] == [_uid(address)] and found["total"] == 1 and found["rows"][0]["virtual"] is True
        assert _search(conn, _uid(address))["rows"] == []      # charted sectors only, as it was
        assert _search(conn, _uid((9, 2, 6)), ("uncharted", True))["rows"] == []      # a cell with no part
        assert _search(conn, _uid(CHARTED[0]), ("uncharted", True))["rows"] == []      # charted
    finally:
        conn.close()
