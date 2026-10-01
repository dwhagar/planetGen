# tests/test_db_id_blocks_edges.py

"""
TEST.16 "Id blocks after reset and rollback": the `id_blocks` id
allocator (`_db._allocate_id`/`_reserve_id_block`/`forget_id_blocks`,
schema v45, PERF.13) at its edges -- a `resetDb.py` run in the same
process, a row inserted by hand with an explicit high id between
allocations, a rolled-back transaction, and two processes using up
several blocks of the same table at once.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import multiprocessing
import sys

import pymysql
import pytest

import resetDb
from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem
from tests.bughunt_support import mysql_argv
from tests.conftest import _test_server_kwargs


def _sector(name, system_names):
    sector = SpaceSector(name, edge_ly=11.5)
    for index, system_name in enumerate(system_names):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        system.name = system_name
        sector.add_system(system, position=(float(index), 0.0, 0.0), system_config=cfg)
    return sector


def _insert_configs(conn, count):
    """`count` new `system_configs` rows in one transaction; their ids."""
    with conn:
        return [conn.execute("INSERT INTO system_configs (markdown) VALUES (?)", (0,)).lastrowid
                for _ in range(count)]


def _ids(config, table):
    conn = _db.get_connection(config)
    try:
        return [row["id"] for row in conn.execute(f"SELECT id FROM {table} ORDER BY id").fetchall()]
    finally:
        conn.close()


def test_a_reset_in_the_same_process_restarts_ids_without_collisions(mysql_config, monkeypatch):
    first_sector = _db.save_sector(_sector("Before Reset", ["Haldor", "Imrith"]), config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        before = _insert_configs(conn, 3)
    finally:
        conn.close()
    assert first_sector == 1
    # Cached, partly used blocks for these tables survive in this process.
    assert (mysql_config._key(), "system_configs") in _db._id_blocks

    monkeypatch.setattr(sys, "argv", ["resetDb.py", *mysql_argv(mysql_config), "--yes"])
    resetDb.main()
    assert _ids(mysql_config, "star_systems") == []

    second_sector = _db.save_sector(_sector("After Reset", ["Haldor", "Imrith", "Velos"]), config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        after = _insert_configs(conn, 3)
    finally:
        conn.close()
    assert second_sector == 1
    assert _ids(mysql_config, "star_systems") == [1, 2, 3]
    # The sector's own three configs come first, as they did before the reset.
    assert after == [4, 5, 6] and before == [3, 4, 5]
    assert _ids(mysql_config, "system_configs") == [1, 2, 3, 4, 5, 6]


def test_a_hand_inserted_high_id_is_never_handed_out_again(mysql_config):
    conn = _db.get_connection(mysql_config)
    other = pymysql.connect(database=mysql_config.database, **_test_server_kwargs())
    try:
        first = _insert_configs(conn, 1)
        assert first == [1]
        # One id just past this process's cached block, one further on.
        with other.cursor() as cur:
            cur.execute("INSERT INTO system_configs (id, markdown) VALUES (65, 0), (300, 0)")
        other.commit()
        allocated = first + _insert_configs(conn, 200)
    finally:
        other.close()
        conn.close()

    assert len(set(allocated)) == 201
    assert not {65, 300} & set(allocated)
    assert allocated[:64] == list(range(1, 65))
    assert min(allocated[64:]) > 300
    assert _ids(mysql_config, "system_configs") == sorted(allocated + [65, 300])


def test_a_rolled_back_transaction_leaves_a_gap_never_a_reused_id(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        committed = _insert_configs(conn, 2)
        with pytest.raises(RuntimeError):
            with conn:
                for _ in range(3):
                    conn.execute("INSERT INTO system_configs (markdown) VALUES (?)", (0,))
                raise RuntimeError("roll back")
        # Held back in a batch and never written: still a gap.
        with pytest.raises(RuntimeError):
            with conn:
                with conn.batched():
                    for _ in range(2):
                        conn.execute("INSERT INTO system_configs (markdown) VALUES (?)", (0,))
                    raise RuntimeError("roll back")
        after = _insert_configs(conn, 1)
        # A fresh cache (a restarted process) starts past the whole
        # reserved block: the reservation itself was not rolled back.
        _db.forget_id_blocks(mysql_config._key())
        fresh = _insert_configs(conn, 1)
    finally:
        conn.close()

    assert committed == [1, 2]
    assert after == [8]
    assert fresh == [65]
    assert _ids(mysql_config, "system_configs") == [1, 2, 8, 65]


_ROWS_PER_WORKER = 700  # blocks of 64 + 128 + 256 + 512: four blocks each
_ROWS_PER_COMMIT = 50


def _fill_in_worker(config, barrier, results):
    conn = _db.get_connection(config)
    try:
        barrier.wait()
        ids = []
        while len(ids) < _ROWS_PER_WORKER:
            ids += _insert_configs(conn, min(_ROWS_PER_COMMIT, _ROWS_PER_WORKER - len(ids)))
        results.put(ids)
    finally:
        conn.close()


def test_two_processes_using_up_blocks_of_one_table_never_collide(mysql_config):
    _db.get_connection(mysql_config).close()
    context = multiprocessing.get_context("spawn")
    barrier, results = context.Barrier(2), context.Queue()
    workers = [context.Process(target=_fill_in_worker, args=(mysql_config, barrier, results)) for _ in range(2)]
    for worker in workers:
        worker.start()
    try:
        ids = [results.get(timeout=120) for _ in workers]
    finally:
        for worker in workers:
            worker.join(timeout=60)
    assert [worker.exitcode for worker in workers] == [0, 0]

    assert all(len(worker_ids) == _ROWS_PER_WORKER for worker_ids in ids)
    assert not set(ids[0]) & set(ids[1])
    assert _ids(mysql_config, "system_configs") == sorted(ids[0] + ids[1])
