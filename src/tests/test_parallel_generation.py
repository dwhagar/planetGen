# tests/test_parallel_generation.py

"""
TEST.19 and TEST.22: parallel generation (`--workers`, PERF.7/PERF.8).

- TEST.19: one run seed gives the same tasks' draws with one worker, two
  or several -- the one-worker path seeds each task exactly as a worker
  would.
- GEN.39: one galaxy seed gives the same sectors, contents and all, at
  any worker count, and every one of its 128 bits counts.
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
from planetgen.names import uniqueness
from stellarObjects._db import MySQLConfig

from tests.conftest import _test_server_kwargs
from tests.galaxy_fingerprint import galaxy_rows
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


GALAXY_SEED = bytes.fromhex("0123456789ABCDEF0011223344556677")


def _seeded_galaxy(make_database, seed=GALAXY_SEED):
    config = make_database()
    _plan_wide_galaxy(config, galaxy_seed=seed)
    return config


@pytest.mark.parametrize("workers", [2, 4])
def test_one_galaxy_seed_makes_the_same_sectors_at_any_worker_count(make_database, monkeypatch, workers):
    """GEN.39: each sector draws from its own seed (the galaxy seed and its
    address), so one worker and several fill the same sectors with the
    same systems, stars, planets and phenomena, whatever order they ran."""
    results = {}
    for count in (1, workers):
        config = _seeded_galaxy(make_database)
        counts = _run("galaxy", ["--ring", "1", "--num-systems", "4", "--workers", str(count)],
                      config, monkeypatch)
        assert counts["sectors"] == len(_systems_per_sector(config)) == generate.ring_sector_count(1)
        assert counts["systems"] == sum(_systems_per_sector(config).values())
        results[count] = galaxy_rows(config)
    assert results[1]["star_systems"] and results[1]["planets"]
    assert results[1] == results[workers]


def test_galaxy_seeds_differing_only_in_their_high_64_bits_make_different_sectors(make_database, monkeypatch):
    """GEN.39: every bit of the 128-bit seed counts, not just the low 64."""
    other = bytes([GALAXY_SEED[0] ^ 0x80]) + GALAXY_SEED[1:]
    assert other[8:] == GALAXY_SEED[8:]
    results = []
    for seed in (GALAXY_SEED, other):
        config = _seeded_galaxy(make_database, seed)
        _run("galaxy", ["--ring", "1", "--num-systems", "4", "--workers", "1"], config, monkeypatch)
        results.append(galaxy_rows(config))
    assert results[0]["star_systems"] != results[1]["star_systems"]


def _systems_by_sector(config):
    conn = _db.get_connection(config)
    try:
        rows = conn.execute(
            "SELECT s.ring_index, s.layer_index, s.ring_slot_index, ss.name, ss.position_x_mpc, ss.position_y_mpc,"
            " ss.position_z_mpc FROM star_systems ss JOIN sectors s ON s.id = ss.sector_id").fetchall()
    finally:
        conn.close()
    by_sector = {}
    for row in rows:
        values = tuple(row.values())
        # Which of two colliding names keeps the plain one depends on what
        # else is filled (GEN.57, phase 1): compare without it.
        name = uniqueness.strip_decoration(values[3])
        by_sector.setdefault(values[:3], []).append((name,) + values[4:])
    return {address: sorted(systems) for address, systems in by_sector.items()}


def test_a_sector_comes_out_the_same_whichever_run_fills_it(make_database, monkeypatch):
    """A sector's seed is its address's, not its place in the run: filled
    in a run of three or with all of its ring, a sector holds the same."""
    part = _seeded_galaxy(make_database)
    _run("galaxy", ["--ring", "1", "--limit", "3", "--num-systems", "4", "--workers", "1"], part, monkeypatch)
    whole = _seeded_galaxy(make_database)
    _run("galaxy", ["--ring", "1", "--num-systems", "4", "--workers", "2"], whole, monkeypatch)
    some, every = _systems_by_sector(part), _systems_by_sector(whole)
    assert len(some) == 3
    assert some == {address: every[address] for address in some}


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
