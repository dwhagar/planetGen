# planetgen/db/uncharted_index.py

"""
The index of uncharted sectors (MAP.402): a `sectors` row for every cell a scatter left something in.

The scatters (`bright_stars`, `phenomenon_scatter`) keep only the address of the cell each object lies in, so a cell nothing was generated in had no
row, no entry in the Sectors list and nothing for a filter to find, unless a drifting rogue planet landed in it (`store.uncharted_sector_id`).
`index` makes the same placeholder row (a placed sector with no cell in the chart ledger, MAP.399) for every cell that holds a waiting bright star or
an unbuilt scattered phenomenon and has none yet, so the list, the search and the sector pages (`web.sector_page.placed_in_uncharted`) can reach
what the long-range scans detected. It is made at the end of a scatter, where the block light totals are made (`generation.run_plan.refresh_block_lights`
and `refresh_block_lights_unit` call it after `db.block_lights.rebuild`), by `planetgen fill-block-totals` and by the migration that introduced it
(0223, for a galaxy already scattered). It is not part of `block_lights.rebuild`: that runs inside the old migration 0147, on a schema that has no
`sectors.uid` yet.

The rows are made a shell at a time in batches, each batch one transaction run through `db.contention.run_with_retries` so a deadlock with
another writer waits its turn instead of failing the unit. A second run adds nothing (a cell that has a row is left alone, and an insert of an
address that appeared meanwhile is ignored), so a stopped run is started again. A cell whose objects are later built or cleared keeps its row: it is then
an empty uncharted sector (`query.sector_coverage_condition`), which the Sectors list leaves out unless asked.
"""

import math

from planetgen.db import contention, store
from planetgen.galaxy import uid as galaxyUid
from planetgen.galaxy.spherical import provisional_sector_designation, sector_position_pc
from planetgen.util import log

BATCH_CELLS = 500
"""int: Cells one transaction makes rows for."""

_WAITING_CELLS = (
    "SELECT DISTINCT o.band_index, o.wedge_index FROM {table} o WHERE o.shell_index = ? AND {waiting} AND NOT EXISTS "
    "(SELECT 1 FROM sectors s WHERE s.shell_index = o.shell_index AND s.band_index = o.band_index AND s.wedge_index = o.wedge_index)")
"""str: The cells of one shell that hold a waiting object of `{table}` and have no sector row."""

_OBJECTS = (("bright_stars", "o.star_system_id IS NULL"), ("phenomenon_scatter", "o.built_at IS NULL"))
"""tuple: The tables that hold what waits in a cell, and the condition that says an object of it still waits."""

_INSERT_HEAD = ("INSERT IGNORE INTO sectors (name, edge_mpc, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc, shell_index, band_index, wedge_index,"
                " generator_epoch, run_id, uid, parity_dirty) VALUES ")


def shells_with_objects(conn):
    """The shell indexes that hold a waiting bright star or an unbuilt scattered phenomenon, as a sorted list."""
    shells = set()
    for table, waiting in _OBJECTS:
        shells.update(row["shell_index"] for row in conn.execute(
            f"SELECT DISTINCT o.shell_index FROM {table} o WHERE {waiting}").fetchall())
    return sorted(shells)


def cells_without_a_row(conn, shell):
    """The `(band, wedge)` of every cell of `shell` that holds a waiting object and has no sector row, sorted."""
    cells = set()
    for table, waiting in _OBJECTS:
        query = _WAITING_CELLS.format(table=table, waiting=waiting)
        cells.update((row["band_index"], row["wedge_index"]) for row in conn.execute(query, (shell,)).fetchall())
    return sorted(cells)


def _batch(conn, shell, cells, edge_pc, generator_epoch, run_id):
    """Writes the rows of `cells` (`(band, wedge)` of `shell`) in one transaction; returns the rows made."""
    values, params = [], []
    for band, wedge in cells:
        address = (shell, band, wedge)
        x, y, z = sector_position_pc(*address, edge_pc)
        values.append("(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1)")
        params.extend((provisional_sector_designation(*address), edge_pc * 1000.0, x, y, z, math.sqrt(x * x + y * y + z * z), *address,
                       generator_epoch, run_id, galaxyUid.sector_uid(*address)))
    made = conn.execute(_INSERT_HEAD + ", ".join(values), params).rowcount
    conn.commit()
    return made


def index(conn, edge_pc=None, on_progress=None, on_phase=None):
    """
    Makes the placeholder sector row of every cell that holds a waiting bright star or an unbuilt scattered phenomenon and has none (module note).

    Args:
        conn (Connection): A writable connection (not inside a transaction).
        edge_pc (float, optional): The grid's sector edge in parsecs; the galaxy's own by default.
        on_progress (callable, optional): `on_progress(done_shells, total_shells)` after each shell.
        on_phase (callable, optional): `on_phase(text)` with what the run says about itself.

    Returns:
        int: The rows made.
    """
    from planetgen.db import query
    from planetgen.util import generator_epoch as generatorEpoch

    say = on_phase or (lambda text: None)
    say("Looking for the sector shells that hold waiting objects...")
    shells = shells_with_objects(conn)
    made = 0
    if shells:      # nothing is read of the galaxy, and no run is recorded, until there is a cell to make a row for
        edge_pc = query.galaxy_edge_pc(conn) if edge_pc is None else edge_pc
        epoch = generatorEpoch.generator_epoch()
        run_id = store.sector_run_id(conn)
        conn.commit()
    for done, shell in enumerate(shells, start=1):
        cells = cells_without_a_row(conn, shell)
        for start in range(0, len(cells), BATCH_CELLS):
            chunk = cells[start:start + BATCH_CELLS]
            write = lambda c, chunk=chunk: _batch(c, shell, chunk, edge_pc, epoch, run_id)      # noqa: E731
            if getattr(conn, "_config", None) is None:      # a migration's connection: nothing else writes while it runs
                made += write(conn)
            else:
                made += contention.run_with_retries(conn, write, gate="uncharted-sectors", label=f"Uncharted sectors, shell {shell}")
        if on_progress is not None:
            on_progress(done, len(shells))
    say(f"Made {made:,} uncharted sector rows in {len(shells):,} shells." if made else "Every cell that holds a waiting object already has its sector row.")
    log.debug(f"uncharted_index: {made} rows")
    return made


_UNCHARTED = (
    "(sec.shell_index IS NOT NULL AND NOT EXISTS (SELECT 1 FROM chart_ledger cl WHERE cl.shell_index = sec.shell_index "
    "AND cl.band_index = sec.band_index AND cl.wedge_index = sec.wedge_index))")
_UNOCCUPIED = ("NOT EXISTS (SELECT 1 FROM rogue_planets rp WHERE rp.sector_id = sec.id) "
               "AND NOT EXISTS (SELECT 1 FROM star_systems ss WHERE ss.sector_id = sec.id)")
CLEAR_CHUNK = 5000


def clear_unoccupied(conn, on_phase=None):
    """
    Deletes the placeholder rows nothing lives in (no rogue planet filed in them, no system), for a new plan: the scatters they were made for are
    cleared with the old outline, so they would be empty rows for cells that may not exist in the new one. A row something references stays (the
    delete is skipped for it by the database's own keys). Returns the rows deleted.
    """
    deleted, last = 0, 0
    while True:
        ids = [r["id"] for r in conn.execute(
            f"SELECT sec.id FROM sectors sec WHERE sec.id > ? AND {_UNCHARTED} AND {_UNOCCUPIED} ORDER BY sec.id LIMIT ?", (last, CLEAR_CHUNK)).fetchall()]
        if not ids:
            break
        last = ids[-1]
        deleted += conn.execute(f"DELETE IGNORE FROM sectors WHERE id IN ({', '.join('?' for _ in ids)})", ids).rowcount
        conn.commit()
    if on_phase is not None and deleted:
        on_phase(f"Cleared {deleted:,} uncharted sector rows nothing lives in.")
    return deleted
