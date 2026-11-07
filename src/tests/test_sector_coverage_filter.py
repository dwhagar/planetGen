# tests/test_sector_coverage_filter.py

"""
MAP.403: the Sectors list, the API and the search list charted sectors only unless asked for more: `coverage=all` or `uncharted`, and `empty=1` for the
uncharted sectors nothing is known to be in (no rogue planet filed in them, no waiting bright star, no unbuilt scattered phenomenon).
"""

import pytest

from planetgen.api.config import Config
from planetgen.db import chart_ledger, query, store
from planetgen.web.app import create_app

pytestmark = pytest.mark.db

EDGE_PC = 4.0
CHARTED = [(9, 2, 5), (9, 2, 6)]
WITH_SCATTER = (9, 2, 7)       # an unbuilt scattered phenomenon waits in it
WITH_STAR = (10, 3, 4)         # a bright star no system was built around waits in it
EMPTY = (11, 1, 2)             # nothing is known to be in it
BUILT = (12, 1, 3)             # its scattered phenomenon was built: nothing waits any more
NOT_EMPTY = [WITH_SCATTER, WITH_STAR]
EMPTIES = [EMPTY, BUILT]


def _scatter(conn, address, built):
    conn.execute(
        "INSERT INTO phenomenon_scatter (shell_index, band_index, wedge_index, kind, position_x_mpc, position_y_mpc, position_z_mpc, seed, built_at)"
        f" VALUES (?, ?, ?, 'black-hole', 1.0, 2.0, 3.0, 7, {'CURRENT_TIMESTAMP' if built else 'NULL'})", address)


def _bright_star(conn, address):
    conn.execute(
        "INSERT INTO bright_stars (shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc, position_z_mpc, population, star_type, yerkes_class,"
        " mass_kg, radius_km, temperature_k, luminosity_w, age_gy, initial_mass_sol, seed) VALUES (?, ?, ?, 1.0, 2.0, 3.0, 'young', 'B2', 'V', 3e31, 5e6, 2e4,"
        " 1e30, 0.1, 15.0, 9)", address)


@pytest.fixture
def galaxy(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            for address in CHARTED + NOT_EMPTY + EMPTIES:
                sector_id = store.uncharted_sector_id(conn, address, EDGE_PC)
                if address in CHARTED:
                    chart_ledger.record(conn, address, store.sector_run_id(conn), 1)
                    conn.execute("UPDATE sectors SET name = ? WHERE id = ?", (f"Charted {address[2]}", sector_id))
            _scatter(conn, WITH_SCATTER, built=False)
            _scatter(conn, BUILT, built=True)
            _bright_star(conn, WITH_STAR)
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


def _addresses(body):
    return {(item["shell_index"], item["band_index"], item["wedge_index"]) for item in body["items"]}


@pytest.mark.parametrize("query_string, expected, total", [
    ("", set(CHARTED), 2),
    ("coverage=charted&empty=1", set(CHARTED), 2),
    ("coverage=uncharted", set(NOT_EMPTY), 2),
    ("coverage=uncharted&empty=1", set(NOT_EMPTY + EMPTIES), 4),
    ("coverage=all", set(CHARTED + NOT_EMPTY), 4),
    ("coverage=all&empty=1", set(CHARTED + NOT_EMPTY + EMPTIES), 6),
])
def test_the_sectors_endpoint_lists_what_coverage_and_empty_ask_for(client, query_string, expected, total):
    body = client.get(f"/api/sectors?limit=50&{query_string}").get_json()
    assert _addresses(body) == expected and body["total"] == total


def test_the_counts_ignore_what_is_listed(client):
    for query_string in ("", "coverage=all&empty=1", "coverage=uncharted"):
        body = client.get(f"/api/sectors?limit=1&{query_string}").get_json()
        assert (body["charted_total"], body["uncharted_total"], body["empty_total"]) == (2, 4, 2), query_string
    body = client.get("/api/sectors?limit=1&coverage=uncharted&empty=1").get_json()
    assert body["coverage"] == "uncharted" and body["show_empty"] is True
    body = client.get("/api/sectors?limit=1").get_json()
    assert body["coverage"] == "charted" and body["show_empty"] is False


def test_another_coverage_is_a_400(client):
    assert client.get("/api/sectors?coverage=everything").status_code == 400


@pytest.mark.parametrize("sort", [None, "name", "distance", "systems", "density", "position"])
def test_every_sort_keeps_to_the_filter_through_every_page(client, sort):
    sort_part = f"&sort={sort}" if sort else ""
    for query_string, expected in (("coverage=uncharted&empty=1", set(NOT_EMPTY + EMPTIES)), ("", set(CHARTED)), ("coverage=all", set(CHARTED + NOT_EMPTY))):
        seen, offset = [], 0
        while True:
            body = client.get(f"/api/sectors?limit=1&offset={offset}{sort_part}&{query_string}").get_json()
            if not body["items"]:
                break
            seen.append(body["items"][0]["id"])
            offset += 1
            assert offset < 20
        assert len(seen) == len(set(seen)) == len(expected), (sort, query_string, seen)
        listed = _addresses(client.get(f"/api/sectors?limit=50{sort_part}&{query_string}").get_json())
        assert listed == expected


def test_the_database_functions_keep_every_row_unless_asked(galaxy):
    conn = store.get_connection(galaxy)
    try:
        assert len(query.list_sectors(conn)) == 6 and query.count_sectors(conn) == 6
        assert query.count_sectors(conn, coverage="charted") == 2
        assert query.count_sectors(conn, coverage="uncharted", show_empty=False) == 2
        assert query.count_empty_uncharted_sectors(conn) == 2
        with pytest.raises(ValueError):
            query.list_sectors(conn, coverage="nope")
    finally:
        conn.close()


def test_the_search_finds_an_uncharted_sector_only_when_asked(client, galaxy):
    conn = store.get_connection(galaxy)
    try:
        name = conn.execute("SELECT name FROM sectors WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", WITH_SCATTER).fetchone()["name"]
        empty_name = conn.execute("SELECT name FROM sectors WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", EMPTY).fetchone()["name"]
    finally:
        conn.close()

    def found(term, extra=""):
        body = client.get(f"/api/search?sector_q={term}&panels=sectors{extra}").get_json()
        return body["results"]["sectors"]["total"]

    assert found("Charted") == 2
    assert found(name) == 0 and found(name, "&sector_coverage=all") == 1 and found(name, "&sector_coverage=uncharted") == 1
    assert found(name, "&sector_coverage=charted") == 0
    assert found(empty_name, "&sector_coverage=all") == 0 and found(empty_name, "&sector_coverage=all&sector_empty=1") == 1
    assert client.get("/api/search?sector_q=x&sector_coverage=nope").status_code == 400


def test_the_sectors_page_and_its_table_follow_the_filter(client):
    def rows(query_string):
        return client.get(f"/table/sectors?{query_string}").get_json()

    assert rows("")["total"] == 2
    assert rows("coverage=all&empty=1")["total"] == 6
    assert rows("coverage=uncharted")["total"] == 2
    html = client.get("/sectors").get_data(as_text=True)
    assert "2 sectors" in html and "2 uncharted" in html and "2 empty" in html
    assert 'name="coverage"' in html and 'name="empty"' in html
    html = client.get("/sectors?coverage=uncharted&empty=1").get_data(as_text=True)
    assert "coverage=uncharted" in html


def test_the_search_page_keeps_the_sector_scope_in_its_form_and_chips(client):
    html = client.get("/search?sector_q=Charted&sector_coverage=all&sector_empty=1").get_data(as_text=True)
    assert 'name="sector_coverage"' in html and 'name="sector_empty"' in html
    assert "Sectors: all sectors" in html and "Sectors: with the empty ones" in html
    html = client.get("/search?sector_q=Charted").get_data(as_text=True)
    assert "Sectors: all sectors" not in html
