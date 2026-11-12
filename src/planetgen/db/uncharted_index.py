# planetgen/db/uncharted_index.py

"""
The index of uncharted sectors (MAP.402): a `sectors` row for every cell a scatter left something in.

The scatters (`bright_stars`, `phenomenon_scatter`) keep only the address of the cell each object lies in, so a cell nothing was generated in had no
row, no entry in the Sectors list and nothing for a filter to find, unless a drifting rogue planet landed in it (`store.uncharted_sector_id`).
`index` makes the same placeholder row (a placed sector with no cell in the chart ledger, MAP.399) for every cell that holds a waiting bright star or
an unbuilt scattered phenomenon and has none yet, so the list, the search and the sector pages (`web.sector_page.placed_in_uncharted`) can reach
what the long-range scans detected. It is made at the end of a scatter, where the block light totals are made (`generation.run_plan.refresh_block_lights`
and `refresh_block_lights_unit` call it after `db.block_lights.rebuild`), by `planetgen fill-block-totals` and by the migration that introduced it
(0227, for a galaxy already scattered). It is not part of `block_lights.rebuild`: that runs inside the old migration 0147, on a schema that has no
`sectors.uid` yet.

The rows are made a shell at a time in batches, each batch one transaction run through `db.contention.run_with_retries` so a deadlock with
another writer waits its turn instead of failing the unit. A second run adds nothing (a cell that has a row is left alone, and an insert of an
address that appeared meanwhile is ignored), so a stopped run is started again. A cell whose objects are later built or cleared keeps its row: it is then
an empty uncharted sector (`query.sector_coverage_condition`), which the Sectors list leaves out unless asked.
"""

import math
import time

from planetgen.db import contention, store
from planetgen.galaxy import uid as galaxyUid
from planetgen.galaxy.spherical import provisional_sector_designation, sector_position_pc
from planetgen.util import log
from planetgen.util.ticker import Ticker

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


_PARTS_CELLS = ("SELECT p.band_index, p.wedge_index FROM galaxy_block_light_parts p WHERE p.shell_index = ? AND NOT EXISTS "
                "(SELECT 1 FROM sectors s WHERE s.shell_index = p.shell_index AND s.band_index = p.band_index AND s.wedge_index = p.wedge_index)")
"""str: The cells of one shell that hold a part (a waiting object, GEN.1476) and have no sector row."""


def shells_with_objects(conn, from_parts=False):
    """The shell indexes that hold a waiting bright star or an unbuilt scattered phenomenon, as a sorted list. With `from_parts` they are
    the parts' (`db.block_lights`), which the caller has made current: the two object tables are not read again (GEN.1476)."""
    if from_parts:
        return [int(row["shell_index"]) for row in conn.execute("SELECT DISTINCT shell_index FROM galaxy_block_light_parts ORDER BY shell_index").fetchall()]
    shells = set()
    for table, waiting in _OBJECTS:
        shells.update(row["shell_index"] for row in conn.execute(
            f"SELECT DISTINCT o.shell_index FROM {table} o WHERE {waiting}").fetchall())
    return sorted(shells)


UNIT_CELLS = 300000
"""int: About how many cells with a part one unit of the index works (OPS.153): at about 10,000 rows a second, half a minute of work."""


def shell_part_counts(conn):
    """`[(shell, cells with a part)]`, sorted by shell, for cutting the index into units (OPS.153): one walk of the parts' primary key."""
    return [(int(row["shell_index"]), int(row["n"])) for row in conn.execute(
        "SELECT shell_index, COUNT(*) AS n FROM galaxy_block_light_parts GROUP BY shell_index ORDER BY shell_index").fetchall()]


def shell_runs(counts, max_units, unit_cells=None):
    """
    The shells of `counts` (`shell_part_counts`) cut into runs of neighbors for the units of the index (OPS.153): a run takes shells until it holds `unit_cells`
    cells (more when the cut would make more than `max_units` runs; a shell is never split, as the index works a shell at a time).

    Returns:
        list: `(shells, cells)` for each run, in shell order; the shells in a run are the ones in `counts`.
    """
    unit_cells = UNIT_CELLS if unit_cells is None else unit_cells
    total = sum(n for _shell, n in counts)
    target = max(unit_cells, -(-total // max(max_units, 1)))
    runs, shells, held = [], [], 0
    for shell, n in counts:
        shells.append(shell)
        held += n
        if held >= target:
            runs.append((shells, held))
            shells, held = [], 0
    if shells:
        runs.append((shells, held))
    return runs


def cells_without_a_row(conn, shell, from_parts=False):
    """The `(band, wedge)` of every cell of `shell` that holds a waiting object and has no sector row, sorted. With `from_parts` the cells with a
    part are the cells that hold one (GEN.1476)."""
    cells = set()
    if from_parts:
        return sorted((row["band_index"], row["wedge_index"]) for row in conn.execute(_PARTS_CELLS, (shell,)).fetchall())
    for table, waiting in _OBJECTS:
        query = _WAITING_CELLS.format(table=table, waiting=waiting)
        cells.update((row["band_index"], row["wedge_index"]) for row in conn.execute(query, (shell,)).fetchall())
    return sorted(cells)


def _batch(conn, shell, cells, edge_pc, generator_epoch, run_id):
    """Writes the rows of `cells` (`(band, wedge)` of `shell`) in one transaction; returns the rows made."""
    config = getattr(conn, "_config", None)
    ids = [store._allocate_id(config, "sectors") for _cell in cells] if config is not None else [None] * len(cells)
    with_ids = None not in ids      # the ids come from this process's block, as every other sector row's do (a row made with AUTO_INCREMENT could take an id the block holds)
    values, params = [], []
    for (band, wedge), new_id in zip(cells, ids):
        address = (shell, band, wedge)
        x, y, z = sector_position_pc(*address, edge_pc)
        values.append("(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1)" if with_ids else "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1)")
        if with_ids:
            params.append(new_id)
        params.extend((provisional_sector_designation(*address), edge_pc * 1000.0, x, y, z, math.sqrt(x * x + y * y + z * z), *address,
                       generator_epoch, run_id, galaxyUid.sector_uid(*address)))
    head = _INSERT_HEAD.replace("INTO sectors (", "INTO sectors (id, ", 1) if with_ids else _INSERT_HEAD
    made = conn.execute(head + ", ".join(values), params).rowcount
    conn.commit()
    return made


def index(conn, edge_pc=None, on_progress=None, on_phase=None, run_id=None, shells=None, from_parts=False):
    """
    Makes the placeholder sector row of every cell that holds a waiting bright star or an unbuilt scattered phenomenon and has none (module note).

    Args:
        conn (Connection): A writable connection (not inside a transaction).
        edge_pc (float, optional): The grid's sector edge in parsecs; the galaxy's own by default.
        on_progress (callable, optional): `on_progress(done_shells, total_shells)` after each shell.
        on_phase (callable, optional): `on_phase(text)` with what the run says about itself.
        run_id (int, optional): The run the rows are recorded under; with none, the one run that saves made outside a run share (`store.sector_run_id`).
        shells (list, optional): `shells_with_objects(conn)` when the caller has it already (it sized a bar with it).
        from_parts (bool): Take the shells and cells from the block light parts (GEN.1476), which the caller has just made current, and not from the two scatter tables.

    Returns:
        int: The rows made.
    """
    from planetgen.db import query
    from planetgen.util import generator_epoch as generatorEpoch

    say = on_phase or (lambda text: None)
    say("Looking for the sector shells that hold waiting objects...")
    shells = shells_with_objects(conn, from_parts) if shells is None else shells
    made = 0
    if shells:      # nothing is read of the galaxy, and no run is recorded, until there is a cell to make a row for
        edge_pc = query.galaxy_edge_pc(conn) if edge_pc is None else edge_pc
        epoch = generatorEpoch.generator_epoch()
        run_id = store.sector_run_id(conn, run_id)
        conn.commit()
    ticker = Ticker("Uncharted sectors, shells", total=len(shells), say=say)
    for done, shell in enumerate(shells, start=1):
        cells = cells_without_a_row(conn, shell, from_parts)
        for start in range(0, len(cells), BATCH_CELLS):
            chunk = cells[start:start + BATCH_CELLS]
            write = lambda c, chunk=chunk: _batch(c, shell, chunk, edge_pc, epoch, run_id)      # noqa: E731
            if getattr(conn, "_config", None) is None:      # a migration's connection: nothing else writes while it runs
                made += write(conn)
            else:
                made += contention.run_with_retries(conn, write, gate="uncharted-sectors", label=f"Uncharted sectors, shell {shell}")
        if on_progress is not None:
            on_progress(done, len(shells))
        ticker.tick(done, f"{made:,} rows made")
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
