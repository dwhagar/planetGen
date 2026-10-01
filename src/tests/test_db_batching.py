# tests/test_db_batching.py

"""
Database-integration tests for batched sector writes (PERF.12-14):
the schema checked once per process, ids from `id_blocks` (schema v45),
`Connection.batched`'s multi-row INSERTs, the bulk name reservation, the
retry on a deadlock, and several processes saving sectors at once.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import multiprocessing
import pickle

import pymysql
import pytest

from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.roguePlanetData import RoguePlanet
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem


def _sector(name, system_names, rogue_names=()):
    sector = SpaceSector(name, edge_ly=11.5)
    for index, system_name in enumerate(system_names):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        system.name = system_name
        sector.add_system(system, position=(float(index), 0.0, 0.0), system_config=cfg)
    for index, rogue_name in enumerate(rogue_names):
        rogue = RoguePlanet(SystemConfig(), name=rogue_name)
        sector.add_phenomenon(rogue, "rogue-planet", position=(0.0, float(index) + 0.5, 0.0))
    return sector


class _Counter:
    """Counts the statements a `Connection` actually sends."""

    def __init__(self, monkeypatch):
        self.statements = []
        real_run, real_flush = _db.Connection._run, _db.Connection.flush
        counter = self

        def run(conn, sql, params):
            counter.statements.append(sql.split(None, 1)[0].upper())
            return real_run(conn, sql, params)

        def flush(conn):
            if conn._batch:
                counter.statements.extend(
                    ["INSERT"] * sum(-(-len(rows) // _db._BATCH_ROWS) for rows in conn._batch.values()))
            return real_flush(conn)

        monkeypatch.setattr(_db.Connection, "_run", run)
        monkeypatch.setattr(_db.Connection, "flush", flush)


def test_schema_is_applied_once_per_process(mysql_config, monkeypatch):
    _db.get_connection(mysql_config).close()
    calls = []
    monkeypatch.setattr(_db, "_ensure_schema", lambda conn: calls.append(conn))
    for _ in range(3):
        _db.get_connection(mysql_config).close()
    assert calls == []


def test_sector_rows_are_written_in_a_few_multi_row_inserts(mysql_config, monkeypatch):
    sector = _sector("Batch Sector", [f"Batchstar {n}" for n in range(8)], ["Lonely Wanderer"])
    _db.get_connection(mysql_config).close()
    counter = _Counter(monkeypatch)
    sector_id = _db.save_sector(sector, config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        planets = conn.execute("SELECT COUNT(*) AS n FROM planets").fetchone()["n"]
        moons = conn.execute("SELECT COUNT(*) AS n FROM moons").fetchone()["n"]
        systems = conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()
        assert systems["n"] == 8
        assert conn.execute("SELECT COUNT(*) AS n FROM rogue_planets").fetchone()["n"] == 1
        # One INSERT per table/statement shape, not one per row.
        assert counter.statements.count("INSERT") < 40 < planets + moons + 16
        # Every child points at a parent written in the same save.
        orphans = conn.execute(
            "SELECT COUNT(*) AS n FROM moons m LEFT JOIN planets p ON p.id = m.planet_id WHERE p.id IS NULL"
        ).fetchone()["n"]
        assert orphans == 0
        loaded = _db.load_sector(conn, sector_id)
        assert sorted(entry.star_system.name for entry in loaded.entries) == sorted(
            entry.star_system.name for entry in sector.entries)
    finally:
        conn.close()


def test_ids_come_from_id_blocks_above_existing_rows(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("INSERT INTO system_configs (id, markdown) VALUES (5000, 0)")
        _db.forget_id_blocks(mysql_config._key())
        with conn:
            first = _db.insert_system_config(conn, SystemConfig())
            second = _db.insert_system_config(conn, SystemConfig())
        assert first == 5001 and second == 5002
        next_id = conn.execute("SELECT next_id FROM id_blocks WHERE table_name = 'system_configs'").fetchone()
        assert next_id["next_id"] > second
    finally:
        conn.close()


def test_inserts_fall_back_to_auto_increment_without_id_blocks(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("DROP TABLE id_blocks")
        _db.forget_id_blocks(mysql_config._key())
        with conn:
            config_id = _db.insert_system_config(conn, SystemConfig())
        assert conn.execute("SELECT id FROM system_configs").fetchone()["id"] == config_id
    finally:
        conn.close()


def test_migration_to_v45_adds_id_blocks(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("DROP TABLE id_blocks")
            conn.execute("DELETE FROM schema_migrations WHERE version IN (45, 46, 47, 48, 49, 50)")
            conn.execute("INSERT INTO schema_migrations (version) VALUES (44)")
    finally:
        conn.close()
    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION
    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM id_blocks").fetchone()["n"] == 0
    finally:
        conn.close()


def test_equal_names_in_one_sector_get_greek_letters(mysql_config):
    sector = _sector("Twin Sector", ["Kemaral", "Kemaral", "Kemaral"], ["Kemaral"])
    _db.save_sector(sector, config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        systems = sorted(r["name"] for r in conn.execute("SELECT name FROM star_systems").fetchall())
        rogue = conn.execute("SELECT name FROM rogue_planets").fetchone()["name"]
        assert systems == ["Alpha Kemaral", "Beta Kemaral", "Gamma Kemaral"]
        assert rogue == "Delta Kemaral"
        registry = conn.execute(
            "SELECT occurrence_count, first_star_system_id, first_object_table FROM system_name_registry"
            " WHERE base_name = 'Kemaral'"
        ).fetchone()
        first = conn.execute("SELECT id FROM star_systems WHERE name = 'Alpha Kemaral'").fetchone()["id"]
        assert registry["occurrence_count"] == 4
        assert registry["first_star_system_id"] == first
        assert registry["first_object_table"] is None
    finally:
        conn.close()


def test_a_stored_holder_is_renamed_by_a_later_sector(mysql_config):
    _db.save_sector(_sector("First Sector", ["Ossiran"]), config=mysql_config)
    _db.save_sector(_sector("Second Sector", ["Ossiran"]), config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        names = sorted(r["name"] for r in conn.execute("SELECT name FROM star_systems").fetchall())
        assert names == ["Alpha Ossiran", "Beta Ossiran"]
    finally:
        conn.close()


def test_a_deadlocked_save_is_retried_with_the_generated_names(mysql_config, monkeypatch):
    sector = _sector("Retry Sector", ["Retrystar"])
    real_insert = _db.insert_sector
    attempts = []

    def flaky(conn, sector_arg, galaxy_position=None):
        attempts.append(sector_arg.entries[0].star_system.name)
        sector_id = real_insert(conn, sector_arg, galaxy_position=galaxy_position)
        if len(attempts) == 1:
            sector_arg.entries[0].star_system.name = "Mangled"
            raise pymysql.err.OperationalError(1213, "Deadlock found when trying to get lock")
        return sector_id

    monkeypatch.setattr(_db, "insert_sector", flaky)
    monkeypatch.setattr(_db.time, "sleep", lambda seconds: None)
    sector_id = _db.save_sector(sector, config=mysql_config)
    assert attempts == ["Retrystar", "Retrystar"]
    conn = _db.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT sector_id, name FROM star_systems").fetchall()
        assert [(r["sector_id"], r["name"]) for r in rows] == [(sector_id, "Retrystar")]
        assert conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"] == 1
    finally:
        conn.close()


def test_other_errors_are_not_retried(mysql_config, monkeypatch):
    calls = []

    def broken(conn, sector_arg, galaxy_position=None):
        calls.append(1)
        raise pymysql.err.OperationalError(1146, "Table doesn't exist")

    monkeypatch.setattr(_db, "insert_sector", broken)
    with pytest.raises(pymysql.err.OperationalError):
        _db.save_sector(_sector("Broken Sector", []), config=mysql_config)
    assert calls == [1]


def _save_in_worker(args):
    config, payload = args
    sector = pickle.loads(payload)
    return _db.save_sector(sector, config=config)


def test_parallel_writers_with_the_same_names_finish_cleanly(mysql_config):
    """Four processes save sectors whose systems and rogue planets all
    share names, the worst case for the name registry's locks: every
    save succeeds (deadlocks, if any, are retried) and no two rows end
    up with the same name."""
    _db.get_connection(mysql_config).close()
    payloads = [
        pickle.dumps(_sector("Crowded Sector", [f"Common {n}" for n in range(4)], ["Drifter", "Drifter"]))
        for _ in range(4)
    ]
    with multiprocessing.get_context("spawn").Pool(4) as pool:
        sector_ids = pool.map(_save_in_worker, [(mysql_config, payload) for payload in payloads])
    assert len(set(sector_ids)) == 4
    conn = _db.get_connection(mysql_config)
    try:
        for table in ("sectors", "star_systems", "rogue_planets"):
            duplicates = conn.execute(
                f"SELECT name FROM {table} GROUP BY name HAVING COUNT(*) > 1"
            ).fetchall()
            assert duplicates == [], table
        assert conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"] == 16
        assert conn.execute("SELECT COUNT(*) AS n FROM rogue_planets").fetchone()["n"] == 8
    finally:
        conn.close()
