# tests/test_db_web_reads.py

"""
Database-integration tests for the web read path's query work
(PERF.15-17): fewer queries per page, whole-word FULLTEXT name search
with capped counts (schema v46), and the statement time limit on the
web's read-only connections.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import pymysql
import pytest

import queryDb
from api.app import create_app
from api.config import Config
from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem


def _save_sector(config, name, system_names):
    sector = SpaceSector(name, edge_ly=11.5)
    for index, system_name in enumerate(system_names):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "K1V"
        cfg.BINARY_SYSTEM = False
        cfg.MOONS = True
        system = StarSystem(system_config=cfg)
        system.name = system_name
        sector.add_system(system, position=(float(index), 0.0, 0.0), system_config=cfg)
    return _db.save_sector(sector, config=config)


class _Statements:
    """Counts the statements a `Connection` sends."""

    def __init__(self, monkeypatch):
        self.count = 0
        real = _db.Connection._run

        def run(conn, sql, params):
            self.count += 1
            return real(conn, sql, params)

        monkeypatch.setattr(_db.Connection, "_run", run)


def _search(conn, limit=50, **texts):
    return queryDb.search(conn, texts, {facet: set() for facet in queryDb.SEARCH_TAG_FACETS}, limit=limit)


def test_name_search_matches_whole_words_only(mysql_config):
    _save_sector(mysql_config, "Whole Word Sector", ["Kemaral", "Ossiran Mu", "Tobar IV", "Kemaralis"])
    conn = queryDb.open_readonly(mysql_config)
    try:
        def names(term):
            return sorted(row["name"] for row in _search(conn, system_q=term)["results"]["systems"]["rows"])

        assert names("ara") == []
        assert names("kemaral") == ["Kemaral"]
        assert names("Mu") == ["Ossiran Mu"]
        assert names("tobar iv") == ["Tobar IV"]
        assert names("Tobar V") == []
        assert names("  ") == []
    finally:
        conn.close()


def test_galaxy_locate_matches_the_start_of_the_last_word(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            conn.execute(
                "INSERT INTO sectors (name, edge_mpc, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc,"
                " ring_index, layer_index, ring_slot_index) VALUES ('Vashti Corrow', 4000, 1, 1, 1, 1.7, 1, 1, 1)"
            )
    finally:
        conn.close()
    conn = queryDb.open_readonly(mysql_config)
    try:
        assert [m["name"] for m in queryDb.galaxy_locate(conn, "corr")] == ["Vashti Corrow"]
        assert [m["name"] for m in queryDb.galaxy_locate(conn, "Va")] == ["Vashti Corrow"]
        assert queryDb.galaxy_locate(conn, "orrow") == []
    finally:
        conn.close()


def test_result_counts_stop_at_the_cap(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            conn.executemany(
                "INSERT INTO sectors (name, edge_mpc) VALUES (?, 4000)",
                [(f"Crowded Reach {n}",) for n in range(queryDb.SEARCH_COUNT_CAP + 5)],
            )
    finally:
        conn.close()
    conn = queryDb.open_readonly(mysql_config)
    try:
        page = _search(conn, sector_q="crowded")["results"]["sectors"]
        assert page["total"] == queryDb.SEARCH_COUNT_CAP and page["total_capped"]
        small = _search(conn, sector_q="reach 7")["results"]["sectors"]
        assert small["total"] == 1 and not small["total_capped"]
    finally:
        conn.close()


def test_search_facets_are_cached_until_the_content_changes(mysql_config, monkeypatch):
    _save_sector(mysql_config, "Facet Sector", ["Facetstar"])
    calls = []
    real = queryDb._search_facet_type
    monkeypatch.setattr(queryDb, "_search_facet_type", lambda conn: calls.append(1) or real(conn))
    conn = queryDb.open_readonly(mysql_config)
    try:
        first = _search(conn)
        second = _search(conn)
        assert len(calls) == 1
        assert first["facets"] == second["facets"]
    finally:
        conn.close()
    _save_sector(mysql_config, "Second Facet Sector", ["Otherstar"])
    conn = queryDb.open_readonly(mysql_config)
    try:
        _search(conn)
        assert len(calls) == 2
    finally:
        conn.close()


def test_sector_and_system_pages_use_a_fixed_number_of_queries(mysql_config, monkeypatch):
    small_id = _save_sector(mysql_config, "Small Sector", ["Lone"])
    big_id = _save_sector(mysql_config, "Big Sector", [f"Crowd {n}" for n in range(6)])
    counter = _Statements(monkeypatch)
    conn = queryDb.open_readonly(mysql_config)
    try:
        counter.count = 0
        queryDb.sector_detail(conn, small_id)
        small = counter.count
        counter.count = 0
        detail = queryDb.sector_detail(conn, big_id)
        assert counter.count == small
        assert all(system["stars"] for system in detail["systems"])

        system_id = detail["systems"][0]["id"]
        counter.count = 0
        page = queryDb.system_detail(conn, system_id)
        with_planets = counter.count
        assert sum(len(planet["moons"]) for planet in page["planets"]) >= 0
        lone_id = queryDb.sector_detail(conn, small_id)["systems"][0]["id"]
        counter.count = 0
        queryDb.system_detail(conn, lone_id)
        assert counter.count == with_planets

        counter.count = 0
        rows = queryDb.list_systems(conn, limit=50)
        assert len(rows) == 7 and counter.count <= 3
        assert all(row["star_summary"] for row in rows), rows
    finally:
        conn.close()


def test_migration_to_v46_adds_the_fulltext_indexes(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            for table in _db.FULLTEXT_NAME_TABLES:
                conn.execute(f"ALTER TABLE {table} DROP INDEX ft_{table}_name")
            conn.execute("DELETE FROM schema_migrations")
            conn.execute("INSERT INTO schema_migrations (version) VALUES (45)")
    finally:
        conn.close()
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION
    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in _db.FULLTEXT_NAME_TABLES:
            assert _db._has_index(conn, table, f"ft_{table}_name"), table
    finally:
        conn.close()


def test_a_slow_statement_is_stopped(mysql_config):
    conn = queryDb.open_readonly(mysql_config, statement_timeout_s=0.3)
    try:
        with pytest.raises(pymysql.err.OperationalError) as caught:
            conn.execute("SELECT COUNT(*) AS n FROM (SELECT SLEEP(2)) slow").fetchone()
        assert caught.value.args[0] in _db.STATEMENT_TIMEOUT_ERRORS
    finally:
        conn.close()
    unlimited = queryDb.open_readonly(mysql_config)
    try:
        assert unlimited.execute("SELECT SLEEP(0.5) AS s").fetchone()["s"] == 0
    finally:
        unlimited.close()


@pytest.fixture
def client(mysql_config):
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def test_a_timed_out_api_request_is_a_504(client, monkeypatch):
    def slow(*_args, **_kwargs):
        raise pymysql.err.OperationalError(1969, "Query execution was interrupted (max_statement_time exceeded)")

    monkeypatch.setattr("api.routes.list_phenomena", slow)
    response = client.get("/api/phenomena")
    assert response.status_code == 504
    assert "QUERY_TIMEOUT" in response.get_json()["error"]


def test_other_database_errors_stay_500(client, monkeypatch):
    def broken(*_args, **_kwargs):
        raise pymysql.err.OperationalError(2013, "Lost connection to MySQL server during query")

    monkeypatch.setattr("api.routes.list_phenomena", broken)
    response = client.get("/api/phenomena")
    assert response.status_code == 500
