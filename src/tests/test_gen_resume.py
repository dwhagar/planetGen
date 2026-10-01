# tests/test_gen_resume.py

"""
TEST.23 Resume after an interrupted fill and TEST.26 `--force` scatter then
fill (docs/TODO.md).

TEST.23: a `galaxy --ring`, `--ring --shell` or `--block` run stopped
partway (a sector that fails before it is saved, or one whose save fails
inside its transaction) and then run again ends with the same sectors as
one uninterrupted run in a second database: the same set of addresses,
each filled once, with the same number of its own (not bright-star)
systems, and no system or bright-star link left behind by the failed
sector.

TEST.26: the sectors a `plan --bright-stars-only --force` left out (they
were already filled) stay as they were through a later fill of their
ring, get no bright stars from the fill's backfill either, and the
sectors around them build every star the forced scatter placed for them.
"""

import uuid

import pymysql
import pytest

import generate
from stellarObjects import _db, program_constants
from stellarObjects.galaxyGeometry import ring_sector_count, sector_position_pc

from tests.bughunt_support import mysql_argv, run_cli
from tests.conftest import _test_server_kwargs
from tests.test_bright_star_scatter import EDGE_PC, _plan_args, _seed_galaxy
from tests.test_galaxy_gen import _seed_skeleton


@pytest.fixture
def second_mysql_config(mysql_config):
    """A second throwaway database beside `mysql_config`, for the
    uninterrupted run to compare against."""
    kwargs = _test_server_kwargs()
    name = f"planetgen_test_{uuid.uuid4().hex[:16]}"
    conn = pymysql.connect(**kwargs)
    try:
        with conn.cursor() as cur:
            cur.execute(f"CREATE DATABASE `{name}`")
        conn.commit()
    finally:
        conn.close()
    config = _db.MySQLConfig(database=name, **kwargs)
    try:
        yield config
    finally:
        _db.close_pool(config)
        conn = pymysql.connect(**kwargs)
        try:
            with conn.cursor() as cur:
                cur.execute(f"DROP DATABASE IF EXISTS `{name}`")
            conn.commit()
        finally:
            conn.close()


# --- TEST.23 ---------------------------------------------------------------

LAYERS = [(1, 5), (0, 5), (-1, 5)]
"""Three layers out to ring 5 of test_galaxy_gen's toy shape (one system
per sector at density 1, so the backfill draws next to nothing)."""

NUM_SYSTEMS = 2

MODES = {
    "ring": ["--ring", "1", "--layer", "0"],
    "shell": ["--ring", "1", "--shell"],
    "block": ["--block", "3.0.0.0", "--block-layer", "0"],
}


class _Interrupted(RuntimeError):
    pass


def _galaxy(config, mode_argv):
    run_cli("galaxy", [*mode_argv, "--num-systems", str(NUM_SYSTEMS), *mysql_argv(config)])


def _snapshot(config):
    """Per filled address: `(sector rows there, systems in them that
    aren't built around a pre-placed bright star)`, plus the systems
    whose sector is gone and the unbuilt bright stars left in filled
    cells.

    Bright-star systems are left out of the count because, since GEN.30,
    each run backfills once around its own requested sector after its
    sectors are filled: a run cut into `--limit` steps backfills after
    every step, into cells a later step then fills (and builds those
    stars into), while one whole run backfills only once everything is
    filled, so its cells get none. Every star in a filled cell is built
    either way (`unbuilt`)."""
    conn = _db.get_connection(config)
    try:
        sectors = conn.execute(
            "SELECT s.ring_index, s.layer_index, s.ring_slot_index, COUNT(DISTINCT s.id) AS sectors,"
            " COUNT(ss.id) - COUNT(b.id) AS systems FROM sectors s"
            " LEFT JOIN star_systems ss ON ss.sector_id = s.id"
            " LEFT JOIN bright_stars b ON b.star_system_id = ss.id"
            " GROUP BY s.ring_index, s.layer_index, s.ring_slot_index").fetchall()
        orphans = conn.execute(
            "SELECT COUNT(*) AS n FROM star_systems ss LEFT JOIN sectors s ON s.id = ss.sector_id"
            " WHERE s.id IS NULL").fetchone()["n"]
        unbuilt = conn.execute(
            "SELECT COUNT(*) AS n FROM bright_stars b JOIN sectors s ON s.ring_index = b.ring_index"
            " AND s.layer_index = b.layer_index AND s.ring_slot_index = b.ring_slot_index"
            " WHERE b.star_system_id IS NULL").fetchone()["n"]
        dangling = conn.execute(
            "SELECT COUNT(*) AS n FROM bright_stars b LEFT JOIN star_systems ss ON ss.id = b.star_system_id"
            " WHERE b.star_system_id IS NOT NULL AND ss.id IS NULL").fetchone()["n"]
    finally:
        conn.close()
    by_address = {(row["ring_index"], row["layer_index"], row["ring_slot_index"]): (row["sectors"], row["systems"])
                  for row in sectors}
    return by_address, orphans, unbuilt, dangling


def _fail_on_call(monkeypatch, owner, name, call_number):
    """Makes `owner.name` raise `_Interrupted` on its `call_number`-th call
    (1-based); every other call goes through."""
    real = getattr(owner, name)
    calls = [0]

    def wrapper(*args, **kwargs):
        calls[0] += 1
        if calls[0] == call_number:
            raise _Interrupted(f"{name} interrupted at call {call_number}")
        return real(*args, **kwargs)

    monkeypatch.setattr(owner, name, wrapper)
    return calls


def _check_resumed_matches_uninterrupted(mysql_config, second_mysql_config, mode):
    _seed_skeleton(second_mysql_config, layers=LAYERS)
    _galaxy(second_mysql_config, MODES[mode])
    resumed, orphans, unbuilt, dangling = _snapshot(mysql_config)
    whole, *_rest = _snapshot(second_mysql_config)
    assert len(whole) > 4, "the mode should fill more sectors than the interruption point"
    assert set(resumed) == set(whole)
    assert resumed == whole
    assert all(sectors == 1 and systems >= NUM_SYSTEMS for sectors, systems in resumed.values())
    assert (orphans, unbuilt, dangling) == (0, 0, 0)


@pytest.mark.parametrize("mode", sorted(MODES))
def test_a_run_stopped_before_a_save_resumes_to_the_uninterrupted_result(
        mysql_config, second_mysql_config, monkeypatch, mode):
    _seed_skeleton(mysql_config, layers=LAYERS)
    calls = _fail_on_call(monkeypatch, generate, "generate_and_save_sector_at", 4)
    with pytest.raises(_Interrupted):
        _galaxy(mysql_config, MODES[mode])
    assert calls[0] == 4
    partial, *_rest = _snapshot(mysql_config)
    assert len(partial) == 3
    monkeypatch.undo()

    _galaxy(mysql_config, MODES[mode])
    _check_resumed_matches_uninterrupted(mysql_config, second_mysql_config, mode)


@pytest.mark.parametrize("mode", sorted(MODES))
def test_a_save_that_fails_in_its_transaction_leaves_nothing_and_resumes(
        mysql_config, second_mysql_config, monkeypatch, mode):
    # `refresh_containment` runs last inside `insert_sector`'s transaction,
    # after the sector, its systems and their bright-star links are written.
    _seed_skeleton(mysql_config, layers=LAYERS)
    _fail_on_call(monkeypatch, _db, "refresh_containment", 3)
    with pytest.raises(_Interrupted):
        _galaxy(mysql_config, MODES[mode])
    partial, orphans, unbuilt, dangling = _snapshot(mysql_config)
    assert len(partial) == 2
    assert (orphans, unbuilt, dangling) == (0, 0, 0)
    monkeypatch.undo()

    _galaxy(mysql_config, MODES[mode])
    _check_resumed_matches_uninterrupted(mysql_config, second_mysql_config, mode)


def test_resuming_with_limit_steps_fills_each_sector_once(mysql_config, second_mysql_config):
    # A shell run cut into --limit 4 steps (the operator's own way of
    # stopping partway) ends where one whole run does: the same sectors,
    # each filled once with the same systems of its own (`_snapshot`
    # leaves out the bright stars the earlier steps' backfills placed).
    _seed_skeleton(mysql_config, layers=LAYERS)
    total = ring_sector_count(1) * len(LAYERS)
    for _ in range(-(-total // 4)):
        run_cli("galaxy", [*MODES["shell"], "--limit", "4", "--num-systems", str(NUM_SYSTEMS),
                           *mysql_argv(mysql_config)])
    # One more step has nothing left to do.
    before, *_rest = _snapshot(mysql_config)
    run_cli("galaxy", [*MODES["shell"], "--limit", "4", "--num-systems", str(NUM_SYSTEMS),
                       *mysql_argv(mysql_config)])
    assert _snapshot(mysql_config)[0] == before
    assert len(before) == total
    _check_resumed_matches_uninterrupted(mysql_config, second_mysql_config, "shell")


# --- TEST.26 -----------------------------------------------------------------

FORCE_RING = 8
PRE_FILLED = [(FORCE_RING, 0, 1), (FORCE_RING, 0, 3)]
LIMIT = 5


def _cell_counts(config, addresses):
    conn = _db.get_connection(config)
    try:
        return {
            address: (
                conn.execute("SELECT COUNT(*) AS n FROM bright_stars WHERE ring_index = ? AND layer_index = ?"
                             " AND ring_slot_index = ?", address).fetchone()["n"],
                conn.execute("SELECT COUNT(*) AS n FROM bright_stars WHERE ring_index = ? AND layer_index = ?"
                             " AND ring_slot_index = ? AND star_system_id IS NULL", address).fetchone()["n"],
            )
            for address in addresses
        }
    finally:
        conn.close()


def _sector_rows(config):
    conn = _db.get_connection(config)
    try:
        rows = conn.execute(
            "SELECT s.id, s.ring_index, s.layer_index, s.ring_slot_index, COUNT(ss.id) AS systems"
            " FROM sectors s LEFT JOIN star_systems ss ON ss.sector_id = s.id GROUP BY s.id").fetchall()
    finally:
        conn.close()
    return {(row["ring_index"], row["layer_index"], row["ring_slot_index"]): (row["id"], row["systems"])
            for row in rows}


def test_sectors_a_forced_scatter_skipped_fill_correctly_afterwards(mysql_config, monkeypatch):
    # Smaller backfill tiers (GEN.30) keep the run's backfill quick; they
    # still reach the skipped sectors next door.
    monkeypatch.setattr(program_constants, "BRIGHT_STAR_BACKFILL_TIERS", ((10.0, 100.0), (20.0, 250.0)))
    _seed_galaxy(mysql_config)
    args = generate._default_generation_args(config=mysql_config)
    args.num_systems = 1
    for address in PRE_FILLED:
        generate.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    before = _sector_rows(mysql_config)

    forced = generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only", "--force"))
    assert forced["total"] > 0
    assert _cell_counts(mysql_config, PRE_FILLED) == {address: (0, 0) for address in PRE_FILLED}

    # The ring's first slots, around and between the skipped sectors.
    run_cli("galaxy", ["--ring", str(FORCE_RING), "--layer", "0", "--limit", str(LIMIT), "--num-systems", "1",
                       *mysql_argv(mysql_config)])

    after = _sector_rows(mysql_config)
    new = sorted(set(after) - set(before))
    assert new == [(FORCE_RING, 0, slot) for slot in range(LIMIT + len(PRE_FILLED)) if
                   (FORCE_RING, 0, slot) not in PRE_FILLED]
    # The skipped sectors were neither refilled nor touched, and the new
    # fills' backfill put nothing in their cells.
    assert {address: after[address] for address in PRE_FILLED} == {address: before[address]
                                                                     for address in PRE_FILLED}
    assert _cell_counts(mysql_config, PRE_FILLED) == {address: (0, 0) for address in PRE_FILLED}

    conn = _db.get_connection(mysql_config)
    try:
        assert _db.bright_star_scatter_settings(conn)[0] == float(_plan_args(mysql_config).bright_star_min_luminosity)
        for address in new:
            sector_id = after[address][0]
            linked = conn.execute(
                "SELECT b.ring_index, b.layer_index, b.ring_slot_index, ss.sector_id FROM bright_stars b"
                " JOIN star_systems ss ON ss.id = b.star_system_id WHERE ss.sector_id = ?", (sector_id,)).fetchall()
            # Every star placed in this cell (by the forced scatter or the
            # backfill) is built, once, into this sector and no other.
            counts = _cell_counts(mysql_config, [address])[address]
            assert counts[1] == 0
            assert len(linked) == counts[0]
            assert {(row["ring_index"], row["layer_index"], row["ring_slot_index"]) for row in linked} <= {address}
            systems_per_star = conn.execute(
                "SELECT star_system_id, COUNT(*) AS n FROM bright_stars WHERE star_system_id IS NOT NULL"
                " GROUP BY star_system_id HAVING n > 1").fetchall()
            assert systems_per_star == []
        some_scattered = conn.execute(
            "SELECT COUNT(*) AS n FROM bright_stars WHERE ring_index = ? AND layer_index = 0 AND ring_slot_index < ?",
            (FORCE_RING, LIMIT + len(PRE_FILLED))).fetchone()["n"]
    finally:
        conn.close()
    assert some_scattered > 0
