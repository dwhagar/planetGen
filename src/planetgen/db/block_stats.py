# planetgen/db/block_stats.py

"""
Stored aggregates for the Galaxy Map's block ladder (DB.25, `galaxy_block_stats`).

A drill-down stage shows, for every child block of a block, what its generated
sectors hold: how many, their systems, expected systems, stars, summed star
age and luminosity (`query.galaxy_stage`). A map window across the finer
levels of the ladder would need a stage query for each parent block it touches;
with the sums stored per level and cell, `query.galaxy_block_stats` reads the
whole window in one statement.

A row is one block `(m, ring_index, wedge, slab)` of a level above the sector
(`drill.DRILL_LEVELS[:-1]`) holding at least one generated sector; a block
with none has no row. The sums come from `sectors` and `sector_stats` and only
change when a sector is saved or deleted, so:

* `refresh_chain` recomputes the blocks above one sector (a sector visit, a
  deletion, a `galaxy` run saving sectors one at a time);
* `refresh_sectors` does the same for the sectors a run saved, once at its end,
  so the workers never queue on the shared top-level rows;
* `rebuild` computes every row in one pass over the generated sectors
  (`planetgen plan --fill-block-stats`, the migration, and a run on an empty
  table).

When the ladder gains levels (MAP.151b) a `rebuild` fills them.
"""

from planetgen.db import query
from planetgen.galaxy.drill import DRILL_LEVELS, drill_chain_of

LEVELS = DRILL_LEVELS[:-1]
"""tuple: The levels stored: every block size above a sector."""

FULL_REBUILD_ABOVE = 2000
"""int: A refresh that touches more distinct blocks than this is a `rebuild`."""

COLUMNS = ("m", "ring_index", "wedge", "slab", "sector_count", "systems", "expected_systems", "stars", "age_sum",
           "luminosity_sol")
"""tuple: The table's columns in the order `rows` writes them."""

INSERT = (f"INSERT INTO galaxy_block_stats ({', '.join(COLUMNS)}) VALUES ({', '.join('?' for _ in COLUMNS)})")


def chain_blocks(address):
    """The blocks above the sector `(ring, layer, slot)`, one per stored level, coarsest first."""
    ring, layer, slot = address
    return drill_chain_of(ring, layer, slot)[:len(LEVELS)]


def _row(block, sums):
    sectors, systems, expected, stars, age_sum, luminosity = sums
    return (block.m, block.ring, block.wedge, block.slab, int(sectors), int(systems), float(expected), int(stars),
            float(age_sum), float(luminosity))


def clear(conn):
    """Empties the table."""
    conn.execute("DELETE FROM galaxy_block_stats")


def refresh_blocks(conn, blocks):
    """Recomputes the stored row of each block (a row is dropped when its block no longer holds a generated
    sector). Returns how many rows were written."""
    written = 0
    for block in sorted(set(blocks)):
        conn.execute("DELETE FROM galaxy_block_stats WHERE m = ? AND ring_index = ? AND wedge = ? AND slab = ?",
                     (block.m, block.ring, block.wedge, block.slab))
        sums = query.galaxy_block_sums(conn, block)
        if sums is not None:
            conn.execute(INSERT, _row(block, sums))
            written += 1
    return written


def refresh_chain(conn, address):
    """Recomputes the blocks above one sector, `(ring, layer, slot)`; nothing for a sector off the grid."""
    if address is None or any(value is None for value in address):
        return 0
    return refresh_blocks(conn, chain_blocks(tuple(int(value) for value in address)))


def refresh_sectors(conn, addresses, on_progress=None):
    """Recomputes the blocks above each sector in `addresses`, every block once. A set touching more than
    `FULL_REBUILD_ABOVE` blocks, or a table with no rows, is a `rebuild`. Returns the rows written."""
    blocks = {block for address in addresses for block in chain_blocks(tuple(address))}
    has_rows = conn.execute("SELECT 1 FROM galaxy_block_stats LIMIT 1").fetchone() is not None
    if len(blocks) > FULL_REBUILD_ABOVE or (blocks and not has_rows):
        return rebuild(conn, on_progress)
    written = 0
    for done, block in enumerate(sorted(blocks), start=1):
        written += refresh_blocks(conn, [block])
        if on_progress is not None:
            on_progress(done, len(blocks))
    return written


def accumulate(cells, address, systems, expected, stars, mean_age, luminosity):
    """Adds one generated sector's figures to the sums of each block above it in `cells`
    (`{DrillBlock: [sectors, systems, expected, stars, age sum, luminosity]}`), the way
    `query._add_stage_sums` does."""
    for block in chain_blocks(address):
        sums = cells.setdefault(block, [0, 0.0, 0.0, 0.0, 0.0, 0.0])
        sums[0] += 1
        sums[1] += systems or 0
        sums[2] += expected or 0.0
        sums[3] += stars or 0
        sums[4] += (mean_age or 0.0) * (stars or 0)
        sums[5] += luminosity or 0.0


def cell_rows(cells):
    """The table rows for `accumulate`'s result."""
    return [_row(block, sums) for block, sums in sorted(cells.items(), key=lambda item: tuple(item[0]))]


def rebuild(conn, on_progress=None):
    """
    Computes every row from `sectors` and `sector_stats` in one pass, a ring at a time, and replaces the table's
    rows with them in one transaction (readers see the old rows until the commit). Returns the rows written.
    """
    rings = [row["ring_index"] for row in conn.execute(
        "SELECT DISTINCT ring_index FROM sectors WHERE ring_index IS NOT NULL ORDER BY ring_index").fetchall()]
    cells = {}
    for done, ring in enumerate(rings, start=1):
        rows = conn.execute(
            "SELECT s.ring_index, s.layer_index, s.ring_slot_index, st.actual_systems, st.expected_systems,"
            " st.actual_stars, st.mean_age_gy, st.total_luminosity_sol FROM (SELECT * FROM sectors WHERE"
            " ring_index = ?) s LEFT JOIN sector_stats st ON st.ring_index = s.ring_index"
            " AND st.layer_index = s.layer_index AND st.ring_slot_index = s.ring_slot_index", (ring,)).fetchall()
        for r in rows:
            accumulate(cells, (r["ring_index"], r["layer_index"], r["ring_slot_index"]), r["actual_systems"],
                       r["expected_systems"], r["actual_stars"], r["mean_age_gy"], r["total_luminosity_sol"])
        if on_progress is not None:
            on_progress(done, len(rings))
    clear(conn)
    rows = cell_rows(cells)
    for start in range(0, len(rows), 1000):
        conn.executemany(INSERT, rows[start:start + 1000])
    return len(rows)

