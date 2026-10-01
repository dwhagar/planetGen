# tests/test_parallel_generation.py

"""
TEST.19 and TEST.22: parallel generation (`--workers`, PERF.7/PERF.8).

- TEST.19: one run seed gives the same tasks' draws, and the same
  galaxy sectors, with one worker, two or several -- the one-worker path
  seeds each task exactly as a worker would.
- TEST.22: every bulk mode (`--shell`, `--block`, `--column`,
  `--center-sector`, random start, `sector --num-sectors N`) on two
  workers counts what it saved and never fills a sector twice.

Tests that take the `mysql_config` fixture (see `conftest.py`) are
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import random
import sys
import uuid

import pymysql
import pytest

import generate
from stellarObjects import _db, workQueue
from stellarObjects._db import MySQLConfig

from tests.conftest import _test_server_kwargs
from tests.test_galaxy_gen import _mysql_argv, _plan_wide_galaxy


def _draws(payload):
    return [random.random() for _ in range(payload)]


def _run_queue(workers, run_seed=99, keys=8):
    drawn = {}
    with workQueue.WorkQueue("seeds", workers=workers, run_seed=run_seed) as queue:
        for n in range(keys):
            queue.submit("draw", f"sector-{n}", _draws, 3,
                         on_done=lambda result, _s, _w, n=n: drawn.__setitem__(n, result))
    return drawn


# ---------------------------------------------------------------------------
# TEST.19: the queue itself
# ---------------------------------------------------------------------------

def test_one_worker_draws_what_two_and_four_workers_draw():
    one = _run_queue(1)
    assert one == _run_queue(2) == _run_queue(4)
    assert len({tuple(values) for values in one.values()}) == len(one)


def test_one_worker_seeds_each_task_from_the_run_seed_and_its_key():
    one = _run_queue(1, run_seed=5)
    for n, values in one.items():
        random.seed(workQueue.task_seed(5, f"sector-{n}"))
        assert values == [random.random() for _ in range(3)]
    assert _run_queue(1, run_seed=6) != one


def test_one_worker_leaves_the_runs_own_random_stream_alone():
    """Tasks don't consume or reseed the calling process's stream, just as
    when they run in a pool, so whatever the run draws after a batch is
    the same at any worker count."""
    def after(workers):
        random.seed(2024)
        queue_seed = random.getrandbits(32)
        _run_queue(workers, run_seed=queue_seed)
        return random.random()

    assert after(1) == after(2)


def test_a_failing_task_on_one_worker_still_restores_the_stream():
    random.seed(7)
    expected = random.random()
    random.seed(7)
    with pytest.raises(ZeroDivisionError):
        with workQueue.WorkQueue("fails", workers=1, run_seed=1) as queue:
            queue.submit("boom", "k", _divide_by_zero, None)
    assert random.random() == expected


def _divide_by_zero(_payload):
    random.random()
    return 1 / 0


# ---------------------------------------------------------------------------
# TEST.19: a whole galaxy run
# ---------------------------------------------------------------------------

@pytest.fixture
def make_database(_mysql_server_available):
    """Like `mysql_config`, but as many fresh databases as a test asks for."""
    kwargs = _test_server_kwargs()
    made = []

    def make():
        config = MySQLConfig(database=f"planetgen_test_{uuid.uuid4().hex[:16]}", **kwargs)
        conn = pymysql.connect(**kwargs)
        try:
            with conn.cursor() as cur:
                cur.execute(f"CREATE DATABASE `{config.database}`")
            conn.commit()
        finally:
            conn.close()
        made.append(config)
        _db.get_control_connection(config, ensure_schema=True).close()
        return config

    yield make
    for config in made:
        _db.close_pool(config)
        conn = pymysql.connect(**kwargs)
        try:
            with conn.cursor() as cur:
                cur.execute(f"DROP DATABASE IF EXISTS `{config.database}`")
            conn.commit()
        finally:
            conn.close()


def _run(command, argv, config, monkeypatch):
    monkeypatch.setenv("PLANETGEN_CONTROL_DATABASE", config.database)
    for counter in generate.RUN_COUNTS:
        generate.RUN_COUNTS[counter] = 0
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", command] + argv + _mysql_argv(config)
        generate.main()
    finally:
        sys.argv = old_argv
    return dict(generate.RUN_COUNTS)


def _systems_per_sector(config):
    conn = _db.get_connection(config)
    try:
        return {(row["ring_index"], row["layer_index"], row["ring_slot_index"]): row["n"]
                for row in conn.execute(
                    "SELECT s.ring_index, s.layer_index, s.ring_slot_index, COUNT(ss.id) AS n"
                    " FROM sectors s LEFT JOIN star_systems ss ON ss.sector_id = s.id"
                    " GROUP BY s.id, s.ring_index, s.layer_index, s.ring_slot_index").fetchall()}
    finally:
        conn.close()


def test_one_two_and_three_workers_fill_the_same_sectors(make_database, monkeypatch):
    """Sector contents come from the OS's random source by design
    (`spaceSector`'s `SystemRandom`, `utils.reseed_rng`), so no two runs
    hold the same stars; what the worker count must not change is which
    sectors a run fills, how many systems each gets, and what it counts."""
    results = {}
    for workers in (1, 2, 3):
        config = make_database()
        _plan_wide_galaxy(config)
        counts = _run("galaxy", ["--ring", "1", "--num-systems", "4", "--workers", str(workers)],
                      config, monkeypatch)
        results[workers] = _systems_per_sector(config)
        assert counts["sectors"] == len(results[workers]) == generate.ring_sector_count(1)
        assert counts["systems"] == sum(results[workers].values())
    assert results[1] == results[2] == results[3]


# ---------------------------------------------------------------------------
# TEST.22: every bulk mode on two workers
# ---------------------------------------------------------------------------

def _sector_addresses(config):
    conn = _db.get_connection(config)
    try:
        return [tuple(row.values()) for row in conn.execute(
            "SELECT ring_index, layer_index, ring_slot_index FROM sectors").fetchall()]
    finally:
        conn.close()


def _table_count(config, table):
    conn = _db.get_connection(config)
    try:
        return conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]
    finally:
        conn.close()


def _check_run(config, counts, before=()):
    """The run's own counts match what it saved, and no address holds two
    sectors."""
    addresses = _sector_addresses(config)
    assert len(addresses) == len(set(addresses)), "a sector was filled twice"
    assert counts["sectors"] == len(addresses) - len(before)
    return addresses


def _workers(n=2):
    return ["--workers", str(n), "--num-systems", "2"]


@pytest.fixture
def galaxy_db(make_database):
    config = make_database()
    _plan_wide_galaxy(config)
    return config


def test_shell_on_two_workers_fills_every_slot_and_layer_once(galaxy_db, monkeypatch):
    counts = _run("galaxy", ["--ring", "0", "--shell", "--limit", "9", "--yes"] + _workers(), galaxy_db, monkeypatch)
    addresses = _check_run(galaxy_db, counts)
    assert len(addresses) == 9
    assert {ring for ring, _layer, _slot in addresses} == {0}
    assert counts["systems"] == _table_count(galaxy_db, "star_systems")
    # Again: only the rest of the shell, never a second copy.
    again = _run("galaxy", ["--ring", "0", "--shell", "--limit", "6", "--yes"] + _workers(), galaxy_db, monkeypatch)
    _check_run(galaxy_db, again, before=addresses)
    assert again["sectors"] == 6


def test_column_on_two_workers_fills_one_slot_on_every_layer_once(galaxy_db, monkeypatch):
    # A column covers every stored layer, so keep the galaxy thin.
    _db.replace_galaxy_layers([(layer, 999) for layer in range(-2, 3)], config=galaxy_db)
    counts = _run("galaxy", ["--ring", "1", "--slot", "2", "--column"] + _workers(), galaxy_db, monkeypatch)
    addresses = _check_run(galaxy_db, counts)
    assert sorted(addresses) == [(1, layer, 2) for layer in range(-2, 3)]
    again = _run("galaxy", ["--ring", "1", "--slot", "2", "--column"] + _workers(), galaxy_db, monkeypatch)
    assert again["sectors"] == 0
    assert sorted(_sector_addresses(galaxy_db)) == sorted(addresses)


def test_block_on_two_workers_fills_each_sector_once(galaxy_db, monkeypatch):
    block = "3.1.0.0"
    expected = list(generate.block_addresses(generate.parse_drill_key(block), 0))
    argv = ["--block", block, "--block-layer", "0", "--limit", "4", "--yes"] + _workers()
    counts = _run("galaxy", argv, galaxy_db, monkeypatch)
    addresses = _check_run(galaxy_db, counts)
    assert len(addresses) == min(4, len(expected))
    again = _run("galaxy", argv, galaxy_db, monkeypatch)
    addresses = _check_run(galaxy_db, again, before=addresses)
    assert set(addresses) <= {tuple(address) for address in expected}


def test_center_sector_on_two_workers_fills_its_neighborhood_once(galaxy_db, monkeypatch):
    _run("galaxy", ["--ring", "2", "--limit", "1"] + _workers(1), galaxy_db, monkeypatch)
    conn = _db.get_connection(galaxy_db)
    try:
        center = conn.execute("SELECT id FROM sectors").fetchone()["id"]
    finally:
        conn.close()
    before = _sector_addresses(galaxy_db)
    counts = _run("galaxy", ["--center-sector", str(center), "--radius-pc", "9"] + _workers(),
                  galaxy_db, monkeypatch)
    addresses = _check_run(galaxy_db, counts, before=before)
    assert len(addresses) > len(before)
    again = _run("galaxy", ["--center-sector", str(center), "--radius-pc", "9"] + _workers(),
                 galaxy_db, monkeypatch)
    assert again["sectors"] == 0
    assert len(_sector_addresses(galaxy_db)) == len(addresses)


def test_random_start_on_two_workers_fills_each_sector_once(galaxy_db, monkeypatch):
    counts = _run("galaxy", ["--max-ring", "3", "--radius-pc", "9"] + _workers(), galaxy_db, monkeypatch)
    addresses = _check_run(galaxy_db, counts)
    assert addresses
    assert counts["systems"] == _table_count(galaxy_db, "star_systems")


def test_unplaced_sectors_on_two_workers_are_each_saved_once(galaxy_db, monkeypatch):
    counts = _run("sector", ["--num-sectors", "5", "--yes", "--workers", "2"], galaxy_db, monkeypatch)
    conn = _db.get_connection(galaxy_db)
    try:
        rows = conn.execute("SELECT id, name, ring_index FROM sectors").fetchall()
    finally:
        conn.close()
    assert len(rows) == counts["sectors"] == 5
    assert all(row["ring_index"] is None for row in rows)
    assert len({row["name"] for row in rows}) == 5
    assert counts["systems"] == _table_count(galaxy_db, "star_systems")
