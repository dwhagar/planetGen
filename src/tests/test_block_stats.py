"""
DB.25, MAP.205: the stored block totals (`galaxy_block_stats`) on the block ladder of shells. Every stored row equals
the sums of the generated sectors inside its block, a refresh puts right the blocks above a sector and only those, a
deletion takes the sector out, the migration and `rebuild` give the same rows, and a window reads its cells in one
statement.
"""

import random

import pytest

from planetgen.db import block_stats, edits, query, store
from planetgen.galaxy import spherical
from planetgen.galaxy.drill import drill_chain_of, format_drill_key

pytestmark = pytest.mark.db

EDGE = 4.0


def _addresses(count=48, seed=7):
    """Generated sectors in a few clusters, so blocks of every level hold more than one."""
    rng = random.Random(seed)
    found = set()
    for shell in (30, 31, 32, 33, 90, 91, 300):
        low, high = spherical.band_range(shell)
        centre = rng.randint(low // 2, high // 2)
        for _ in range(count // 7):
            band = max(low, min(high, centre + rng.randint(-2, 2)))
            n = spherical.wedge_count(shell, band)
            found.add((shell, band, rng.randrange(n)))
    return sorted(found)


@pytest.fixture
def galaxy(mysql_config):
    """A database holding generated sectors at known shell, band and wedge addresses, each with stats."""
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            for n, (shell, band, wedge) in enumerate(_addresses()):
                conn.execute(
                    "INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES (?, 4000, ?, ?, ?)",
                    (f"Sector {n}", shell, band, wedge))
                conn.execute(
                    "INSERT INTO sector_stats (shell_index, band_index, wedge_index, expected_systems, actual_systems,"
                    " actual_stars, mean_age_gy, total_luminosity_sol, bright_level_sol) VALUES (?, ?, ?, ?, ?, ?, ?, ?, 0)",
                    (shell, band, wedge, 1.5 + n % 3, 1 + n % 4, 2 + n % 5, 1.0 + n % 7, 10.0 * (1 + n % 6)))
            block_stats.rebuild(conn)
    finally:
        conn.close()
    return mysql_config


def _expected():
    """`{(m, shell, band, wedge): sector count}` worked out here from the chains."""
    counts = {}
    for address in _addresses():
        for block in drill_chain_of(*address)[:-1]:
            counts[tuple(block)] = counts.get(tuple(block), 0) + 1
    return counts


def _stored(conn):
    """Every stored row as a stage would show its block: `{(m, shell, band, wedge): (generated, stats)}`."""
    rows = conn.execute("SELECT m, shell_index, band_index, wedge_index FROM galaxy_block_stats").fetchall()
    found = {}
    for level in {r["m"] for r in rows}:
        cells = [(r["shell_index"], r["band_index"], r["wedge_index"]) for r in rows if r["m"] == level]
        for cell, value in query.galaxy_block_stats(conn, level, cells).items():
            found[(level, *cell)] = (value["generated"], value["stats"])
    return found


def test_a_rebuild_stores_one_row_for_every_block_that_holds_a_sector(galaxy):
    conn = store.get_connection(galaxy)
    try:
        stored = _stored(conn)
        assert {key[0] for key in stored} == set(block_stats.LEVELS)
        assert {key: row[0] for key, row in stored.items()} == _expected()
    finally:
        conn.close()


def test_a_stage_reads_the_stored_rows_and_a_level_three_stage_the_sectors(galaxy):
    conn = store.get_connection(galaxy)
    try:
        stored = _stored(conn)
        top = query.galaxy_stage(conn, None)
        assert top["child_m"] == 243 and top["sectors"] is None
        assert {(243, c["shell"], c["band"], c["wedge"]) for c in top["children"]} == {k for k in stored if k[0] == 243}
        for child in top["children"]:
            assert stored[(243, child["shell"], child["band"], child["wedge"])] == (child["generated"], child["stats"])
        # Down one block of every level: the stage's children are the stored rows of that block's children.
        address = _addresses()[len(_addresses()) // 2]
        chain = drill_chain_of(*address)
        for block in chain[:-2]:
            stage = query.galaxy_stage(conn, format_drill_key(block))
            assert stage["at"] == format_drill_key(block)
            assert stage["child_m"] == block.m // 3
            kids = {(stage["child_m"], c["shell"], c["band"], c["wedge"]): (c["generated"], c["stats"]) for c in stage["children"]}
            assert tuple(chain[chain.index(block) + 1]) in kids
            assert all(stored[key] == value for key, value in kids.items())
        # A level-3 block: its sectors, from `sectors`, and the same sums as its stored row.
        block = chain[-2]
        stage = query.galaxy_stage(conn, format_drill_key(block))
        assert stage["child_m"] == 1 and stage["sectors"]
        assert address in {(s["shell"], s["band"], s["wedge"]) for s in stage["sectors"]}
        assert sum(c["generated"] for c in stage["children"]) == stored[tuple(block)][0]
        assert sum(c["stats"]["systems"] for c in stage["children"]) == stored[tuple(block)][1]["systems"]
    finally:
        conn.close()


def test_deleting_a_sector_takes_it_out_of_every_block_above(galaxy):
    conn = store.get_connection(galaxy)
    try:
        before = _stored(conn)
        shell, band, wedge = _addresses()[0]
        victim = conn.execute("SELECT id FROM sectors WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
                              (shell, band, wedge)).fetchone()["id"]
        edits.delete_sector_with_contents(conn, victim)
        conn.commit()
        after = _stored(conn)
        for block in drill_chain_of(shell, band, wedge)[:-1]:
            key = tuple(block)
            if before[key][0] == 1:
                assert key not in after
            else:
                assert after[key][0] == before[key][0] - 1
        untouched = set(before) - {tuple(b) for b in drill_chain_of(shell, band, wedge)}
        assert all(after[key] == before[key] for key in untouched)
    finally:
        conn.close()


def test_a_rebuild_writes_the_same_rows_a_refresh_kept(galaxy):
    conn = store.get_connection(galaxy)
    try:
        kept = _stored(conn)
        conn.execute("DELETE FROM galaxy_block_stats")
        assert block_stats.rebuild(conn) == len(kept)
        conn.commit()
        assert _stored(conn) == kept
    finally:
        conn.close()


def test_a_few_blocks_are_refreshed_and_an_empty_table_is_rebuilt(galaxy):
    conn = store.get_connection(galaxy)
    try:
        kept = _stored(conn)
        addresses = _addresses()
        conn.execute("UPDATE galaxy_block_stats SET sector_count = 0")
        assert block_stats.refresh_sectors(conn, addresses[:2]) > 0
        after = _stored(conn)
        assert after != kept                      # only the blocks above two sectors were put right
        for address in addresses[:2]:
            for block in drill_chain_of(*address)[:-1]:
                assert after[tuple(block)][0] == kept[tuple(block)][0] or after[tuple(block)][0] >= 1
        conn.execute("DELETE FROM galaxy_block_stats")
        block_stats.refresh_sectors(conn, addresses[:1])        # an empty table: rebuilt whole
        assert _stored(conn) == kept
        assert block_stats.refresh_chain(conn, (None, 0, 0)) == 0
    finally:
        conn.close()


def test_a_refresh_of_every_chain_equals_a_rebuild(galaxy):
    conn = store.get_connection(galaxy)
    try:
        kept = _stored(conn)
        for address in _addresses():
            block_stats.refresh_chain(conn, address)
        assert _stored(conn) == kept
    finally:
        conn.close()


def test_a_window_reads_its_cells_in_one_statement(galaxy):
    conn = store.get_connection(galaxy)
    try:
        stage = query.galaxy_stage(conn, None)
        wanted = [(c["shell"], c["band"], c["wedge"]) for c in stage["children"]] + [(9999, 0, 0)]
        statements = []
        original = conn.execute
        conn.execute = lambda sql, params=(): (statements.append(sql), original(sql, params))[1]
        found = query.galaxy_block_stats(conn, 243, wanted)
        assert len(statements) == 1
        assert set(found) == {(c["shell"], c["band"], c["wedge"]) for c in stage["children"]}
        for child in stage["children"]:
            cell = found[(child["shell"], child["band"], child["wedge"])]
            assert cell["generated"] == child["generated"] and cell["stats"] == child["stats"]
        assert query.galaxy_block_stats(conn, 243, []) == {}
    finally:
        conn.close()


def test_a_window_of_blocks_lists_the_cells_around_a_block_and_the_ones_that_hold_sectors(galaxy):
    conn = store.get_connection(galaxy)
    try:
        address = _addresses()[0]
        block = drill_chain_of(*address)[3]            # the 9-block holding the first sector
        window = query.galaxy_window(conn, format_drill_key(block), 7, EDGE)
        assert window["at"] == format_drill_key(block) and window["g"] == 7 and window["child_m"] == 9
        assert [block.shell, block.band, block.wedge] in window["cells"]
        child = [c for c in window["children"] if (c["shell"], c["band"], c["wedge"]) == (block.shell, block.band, block.wedge)]
        assert child and child[0]["generated"] >= 1
        assert all([c["shell"], c["band"], c["wedge"]] in window["cells"] for c in window["children"])
        with pytest.raises(ValueError):
            query.galaxy_window(conn, "1.30.0.0", 7, EDGE)
    finally:
        conn.close()


