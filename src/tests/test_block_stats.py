"""
DB.25: the stored block totals (`galaxy_block_stats`). Every stored row equals what a drill-down stage computes for
that block from the sectors, a `galaxy` run keeps the table current, a deletion takes the sector out, the migration
and `rebuild` give the same rows, and a window reads its cells in one statement.
"""

import importlib

import pytest
import sqlalchemy

from planetgen.db import alembic_runner, block_stats, edits, query, store
from planetgen.galaxy import seed as galaxySeed
from planetgen.galaxy.drill import DrillBlock, format_drill_key

from tests.test_version_key import SEED, _generate

pytestmark = pytest.mark.db


@pytest.fixture
def galaxy(mysql_config, monkeypatch):
    """A planned galaxy with several generated sectors in ring 1 and 2."""
    _generate(mysql_config, monkeypatch, "plan", "--no-bright-stars", "--max-ring", "40",
              "--seed", galaxySeed.format_seed(SEED))
    for ring in ("1", "2"):
        _generate(mysql_config, monkeypatch, "galaxy", "--ring", ring, "--limit", "4", "--num-systems", "2",
                  "--backfill-from", "none")
    return mysql_config


def _stored(conn):
    """Every stored row as a stage would show its block: `{(m, ring, wedge, slab): (generated, stats)}`."""
    rows = conn.execute("SELECT m, ring_index, wedge, slab FROM galaxy_block_stats").fetchall()
    found = {}
    for level in {r["m"] for r in rows}:
        cells = [(r["ring_index"], r["wedge"], r["slab"]) for r in rows if r["m"] == level]
        for cell, value in query.galaxy_block_stats(conn, level, cells).items():
            found[(level, *cell)] = (value["generated"], value["stats"])
    return found


def _from_stages(conn):
    """The same figures from `galaxy_stage`, walking the ladder down from the galaxy."""
    found = {}

    def walk(at):
        stage = query.galaxy_stage(conn, at)
        if stage["child_m"] == 1:
            return
        for child in stage["children"]:
            found[(stage["child_m"], child["ring"], child["wedge"], child["slab"])] = (child["generated"],
                                                                                    child["stats"])
            if stage["child_m"] > 3:
                walk(format_drill_key(DrillBlock(stage["child_m"], child["ring"], child["wedge"], child["slab"])))

    walk(None)
    return found


def test_a_run_leaves_every_block_total_equal_to_the_stage_it_would_draw(galaxy):
    conn = store.get_connection(galaxy)
    try:
        stored = _stored(conn)
        assert stored and {key[0] for key in stored} == set(block_stats.LEVELS)
        assert stored == _from_stages(conn)
    finally:
        conn.close()


def test_a_rebuild_writes_the_same_rows_a_run_kept(galaxy):
    conn = store.get_connection(galaxy)
    try:
        kept = _stored(conn)
        conn.execute("DELETE FROM galaxy_block_stats")
        assert block_stats.rebuild(conn) == len(kept)
        conn.commit()
        assert _stored(conn) == kept
    finally:
        conn.close()


def test_deleting_a_sector_takes_it_out_of_every_block_above(galaxy):
    conn = store.get_connection(galaxy)
    try:
        total = sum(row[0] for key, row in _stored(conn).items() if key[0] == 243)
        victim = conn.execute("SELECT id FROM sectors WHERE ring_index IS NOT NULL ORDER BY id LIMIT 1").fetchone()["id"]
        edits.delete_sector_with_contents(conn, victim)
        conn.commit()
        assert sum(row[0] for key, row in _stored(conn).items() if key[0] == 243) == total - 1
        assert _stored(conn) == _from_stages(conn)
    finally:
        conn.close()


def test_a_window_reads_its_cells_in_one_statement(galaxy):
    conn = store.get_connection(galaxy)
    try:
        stage = query.galaxy_stage(conn, None)
        wanted = [(c["ring"], c["wedge"], c["slab"]) for c in stage["children"]] + [(9999, 0, 0)]
        statements = []
        original = conn.execute
        conn.execute = lambda sql, params=(): (statements.append(sql), original(sql, params))[1]
        found = query.galaxy_block_stats(conn, 243, wanted)
        assert len(statements) == 1
        assert set(found) == {(c["ring"], c["wedge"], c["slab"]) for c in stage["children"]}
        for child in stage["children"]:
            cell = found[(child["ring"], child["wedge"], child["slab"])]
            assert cell["generated"] == child["generated"] and cell["stats"] == child["stats"]
        assert query.galaxy_block_stats(conn, 243, []) == {}
    finally:
        conn.close()


def test_a_few_blocks_are_refreshed_and_a_missing_table_is_rebuilt(galaxy):
    conn = store.get_connection(galaxy)
    try:
        kept = _stored(conn)
        addresses = [tuple(r.values()) for r in conn.execute(
            "SELECT ring_index, layer_index, ring_slot_index FROM sectors WHERE ring_index IS NOT NULL").fetchall()]
        conn.execute("UPDATE galaxy_block_stats SET sector_count = 0")
        assert block_stats.refresh_sectors(conn, addresses[:2]) > 0
        assert _stored(conn) != kept            # only the blocks above two sectors were put right
        conn.execute("DELETE FROM galaxy_block_stats")
        block_stats.refresh_sectors(conn, addresses[:1])        # an empty table: rebuilt whole
        assert _stored(conn) == kept
        assert block_stats.refresh_chain(conn, (None, 0, 0)) == 0
    finally:
        conn.close()


def test_the_migration_fills_the_table_from_the_sectors_already_there(galaxy, monkeypatch):
    conn = store.get_connection(galaxy)
    try:
        kept = _stored(conn)
    finally:
        conn.close()
    monkeypatch.setattr(alembic_runner, "report_progress", lambda *args: None)     # no Alembic context here
    migration = importlib.import_module("planetgen.db.migrations.versions.0086_galaxy_block_stats")
    engine = alembic_runner._engine(galaxy)
    with engine.begin() as connection:
        connection.execute(sqlalchemy.text("DELETE FROM galaxy_block_stats"))
        migration._fill(connection)
    conn = store.get_connection(galaxy)
    try:
        assert _stored(conn) == kept
    finally:
        conn.close()
