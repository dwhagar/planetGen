# planetgen/db/block_lights.py

"""
Stored light sums for the Galaxy Map's sector-block nodes (MAP.259, `galaxy_block_lights`).

The single-space map shows a block it has not opened as one node: the brightest objects inside it and, for
the rest, a glow made from their summed light. A block holds sectors that were generated (their sums are
`galaxy_block_stats`) and, far more of them, sectors nothing has been built in, whose only stored objects are
the pre-placed bright stars (`bright_stars` with no system built around them) and the scattered phenomena
(`phenomenon_scatter` not built). This table sums those two kinds per block of every level above a sector
(`drill.DRILL_LEVELS[:-1]`), so a node's glow and its center of light are one row read.

A row is one block `(m, shell_index, band_index, wedge_index)` holding at least one pre-placed object. Columns:

* `star_count`, `star_luminosity_sol` (summed luminosity), `star_lpx`, `star_lpy`, `star_lpz` (luminosity times
  position, so the center of light is the sum over `star_luminosity_sol`; solar luminosities times parsecs) and
  `star_lt` (luminosity times temperature, the light's mean temperature over `star_luminosity_sol`);
* `point_count`, `point_luminosity_sol`, `point_lpx`, `point_lpy`, `point_lpz`: the same for scattered phenomena
  (`map_luminosity_sol`, no temperature).

Only objects not yet built are counted: a star that has become a system, and a phenomenon built into its sector,
are in `galaxy_block_stats` (the generated sector's total luminosity), so the two tables never count one light
twice. The sums change when a sector is saved (its objects become built, or deleted: the unbuilt ones come
back), when a scatter is drawn or cleared, and when luminosities are filled, so:

* `refresh_chain` / `refresh_sectors` recompute the blocks above one sector or a run's sectors, the way
  `block_stats` does (a level-3 block from its sectors' objects, a bigger one from its children's rows);
* `rebuild` makes every row from the two object tables, a unit of shells at a time, and commits each unit: it is
  `planetgen fill-block-totals`, the migration, and the repair and the proof of the parts below;
* `roll` makes every row from the parts (GEN.1476): `galaxy_block_light_parts` holds, for each sector cell that holds an unbuilt
  object, the same sums for that cell, added to by the two scatter writers (`store.insert_bright_stars` and
  `store.insert_phenomenon_scatter` call `note_rows` in the transaction that writes the rows). A scatter, a backfill and a luminosity-free
  rerun therefore never read the two big tables again to total them; a scatter's last stage rolls the parts up. The parts are checked
  before they are trusted (`parts_current`: their object counts equal the tables' unbuilt counts); when they are not, `roll` is a
  `rebuild`, which also writes the parts again, so an old galaxy, a prune or a missed writer heals on the next stage.
"""

import math
import struct

from planetgen.db import query, store
from planetgen.galaxy.drill import DRILL_LEVELS, DRILL_STEP, DrillBlock, drill_children, drill_block_sectors, drill_chain_of, drill_parent
from planetgen.physics import constants

LEVELS = DRILL_LEVELS[:-1]
"""tuple: The levels stored: every block size above a sector."""

FULL_REBUILD_ABOVE = 2000
"""int: A refresh that touches more distinct blocks than this is a `rebuild`."""

UNIT_LEVEL = 729
"""int: A rebuild works in units of the sector shells one block of this level holds (and commits each)."""

SUMS = ("star_count", "star_luminosity_sol", "star_lpx", "star_lpy", "star_lpz", "star_lt", "point_count",
        "point_luminosity_sol", "point_lpx", "point_lpy", "point_lpz")
"""tuple: The summed columns, in the order a row's sums are kept."""

COLUMNS = ("m", "shell_index", "band_index", "wedge_index") + SUMS
"""tuple: The table's columns in the order `rows` writes them."""

INSERT = f"INSERT INTO galaxy_block_lights ({', '.join(COLUMNS)}) VALUES ({', '.join('?' for _ in COLUMNS)})"

_INSERT_HEAD = f"INSERT INTO galaxy_block_lights ({', '.join(COLUMNS)}) VALUES "
_WRITE_ROWS = 1000

_PC_PER_MPC = 1.0e-3
_SOL_W = constants.SOLAR_LUMINOSITY


def chain_blocks(address):
    """The blocks above the sector `(shell, band, wedge)`, one per stored level, coarsest first."""
    shell, band, wedge = address
    return drill_chain_of(shell, band, wedge)[:len(LEVELS)]


def _order(block):
    return (block.m, block.shell, block.band, block.wedge)


def clear(conn):
    """Empties the table."""
    conn.execute("DELETE FROM galaxy_block_lights")


_STAR_SUMS = (
    f"COUNT(*) AS star_count, SUM(luminosity_w) / {_SOL_W!r} AS star_luminosity_sol,"
    f" SUM(luminosity_w * position_x_mpc) / {_SOL_W / _PC_PER_MPC!r} AS star_lpx,"
    f" SUM(luminosity_w * position_y_mpc) / {_SOL_W / _PC_PER_MPC!r} AS star_lpy,"
    f" SUM(luminosity_w * position_z_mpc) / {_SOL_W / _PC_PER_MPC!r} AS star_lpz,"
    f" SUM(luminosity_w * temperature_k) / {_SOL_W!r} AS star_lt")
_POINT_SUMS = (
    "COUNT(*) AS point_count, COALESCE(SUM(map_luminosity_sol), 0) AS point_luminosity_sol,"
    f" COALESCE(SUM(map_luminosity_sol * position_x_mpc), 0) * {_PC_PER_MPC!r} AS point_lpx,"
    f" COALESCE(SUM(map_luminosity_sol * position_y_mpc), 0) * {_PC_PER_MPC!r} AS point_lpy,"
    f" COALESCE(SUM(map_luminosity_sol * position_z_mpc), 0) * {_PC_PER_MPC!r} AS point_lpz")
_STARS_WHERE = "star_system_id IS NULL"
_POINTS_WHERE = "built_at IS NULL"


def _zero():
    return [0] * len(SUMS)


PARTS_COLUMNS = ("shell_index", "band_index", "wedge_index") + SUMS
"""tuple: The columns of `galaxy_block_light_parts`: the sector cell, then the same sums a block row holds."""

_PARTS_HEAD = f"INSERT INTO galaxy_block_light_parts ({', '.join(PARTS_COLUMNS)}) VALUES "
_PARTS_ADD = " ON DUPLICATE KEY UPDATE " + ", ".join(f"{name} = {name} + VALUES({name})" for name in SUMS)


_STAR_VALUES = "luminosity_w, position_x_mpc, position_y_mpc, position_z_mpc, temperature_k"
_POINT_VALUES = "map_luminosity_sol * 1e0 AS map_luminosity_sol, position_x_mpc, position_y_mpc, position_z_mpc"      # `* 1e0`: a FLOAT is sent as a short decimal, a DOUBLE as itself
_FLOAT = struct.Struct("f")


def _float32(value):
    """`value` as a `FLOAT` column keeps it (the phenomena's map luminosity is one), so a part made from the rows being written equals the sum of the rows read back."""
    return _FLOAT.unpack(_FLOAT.pack(value))[0]


def cell_sums(stars, values):
    """
    The sums of one sector cell's objects, the one way every path makes them (GEN.1476): `values` is a list of `(light, x, y, z, temperature)` (`temperature`
    unused and a `light` of `None` left out of the sums, as `SUM` leaves a NULL out, for the phenomena), each sum made with `math.fsum` so it does not depend on
    the order the rows came in, and the same rows give the same bits whether they are summed as they are written, read back by a rebuild or read for a refresh.

    Returns:
        list: The six sums of the stars (the other five zero) or the five of the phenomena (the first six zero), in `SUMS` order.
    """
    lit = [v for v in values if v[0] is not None]
    if stars:
        scale = _SOL_W / _PC_PER_MPC
        return [len(values), math.fsum(v[0] for v in lit) / _SOL_W, math.fsum(v[0] * v[1] for v in lit) / scale, math.fsum(v[0] * v[2] for v in lit) / scale,
                math.fsum(v[0] * v[3] for v in lit) / scale, math.fsum(v[0] * v[4] for v in lit) / _SOL_W, 0, 0.0, 0.0, 0.0, 0.0]
    return [0, 0.0, 0.0, 0.0, 0.0, 0.0, len(values), math.fsum(v[0] for v in lit), math.fsum(v[0] * v[1] for v in lit) * _PC_PER_MPC,
            math.fsum(v[0] * v[2] for v in lit) * _PC_PER_MPC, math.fsum(v[0] * v[3] for v in lit) * _PC_PER_MPC]


def note_rows(conn, table, columns, rows):
    """
    Adds the light of freshly scattered objects to their cells' parts (GEN.1476). Called by `store.insert_bright_stars` and
    `store.insert_phenomenon_scatter` in the transaction that writes `rows`, so a batch that is rolled back leaves no part and a retried one
    adds once. The sums are `cell_sums`, the ones a rebuild makes.

    Args:
        conn (Connection): The connection that wrote `rows`.
        table (str): `bright_stars` or `phenomenon_scatter`.
        columns (tuple): The columns of a row, `store.BRIGHT_STAR_COLUMNS` or `store.PHENOMENON_SCATTER_COLUMNS`.
        rows (list): The rows written (a list, not a generator).

    Returns:
        int: The part rows written.
    """
    if not rows:
        return 0
    stars = table == "bright_stars"
    cell_at = [columns.index(name) for name in ("shell_index", "band_index", "wedge_index")]
    x_at, y_at, z_at = (columns.index(name) for name in ("position_x_mpc", "position_y_mpc", "position_z_mpc"))
    light_at = columns.index("luminosity_w" if stars else "map_luminosity_sol")
    temperature_at = columns.index("temperature_k") if stars else None
    cells = {}
    for row in rows:
        light = row[light_at]
        if light is not None and not stars:
            light = _float32(light)
        cells.setdefault((row[cell_at[0]], row[cell_at[1]], row[cell_at[2]]), []).append(
            (light, row[x_at], row[y_at], row[z_at], row[temperature_at] if stars else None))
    held = set()
    addresses = sorted(cells)
    for start in range(0, len(addresses), 500):
        chunk = addresses[start:start + 500]
        held.update((row["shell_index"], row["band_index"], row["wedge_index"]) for row in conn.execute(
            "SELECT shell_index, band_index, wedge_index FROM galaxy_block_light_parts WHERE (shell_index, band_index, wedge_index) IN ("
            + ", ".join("(?, ?, ?)" for _ in chunk) + ")", [value for cell in chunk for value in cell]).fetchall())
    # A cell with no part yet gets the sums of this batch. A cell that has one holds an earlier batch's rows too, and two sums added are not the sum of all the rows
    # (the last bit differs), so its part is made again from the cell's rows, which include this batch's: one cell, a few rows, and the bits stay those of a rebuild.
    # (`ON DUPLICATE KEY UPDATE` only covers a cell two writers meet at once; it adds, and a refresh or rebuild puts the last bit right.)
    _write_parts(conn, [tuple(address) + tuple(cell_sums(stars, cells[address])) for address in addresses if address not in held], add=True)
    if held:
        sync_cells(conn, sorted(held))
    return len(cells)


def _add(into, row, names):
    for name in names:
        into[SUMS.index(name)] += row[name] or 0


def _object_sums(conn, where, params, by_address=False):
    """The sums of the unbuilt objects `where` picks (`params` its parameters): `{(shell, band, wedge): sums}` with
    `by_address`, else one sums list. By cell, a cell holding one object has the sums the database gives and a cell holding more has `cell_sums` of
    its rows read back, so a cell's sums are the same bits whatever order the database would have added them in (GEN.1476)."""
    group = " GROUP BY shell_index, band_index, wedge_index" if by_address else ""
    keys = "shell_index, band_index, wedge_index, " if by_address else ""
    found = {}
    for table, selected, built, names, values in (
            ("bright_stars", _STAR_SUMS, _STARS_WHERE, SUMS[:6], _STAR_VALUES),
            ("phenomenon_scatter", _POINT_SUMS, _POINTS_WHERE, SUMS[6:], _POINT_VALUES)):
        stars = table == "bright_stars"
        many = []
        for row in conn.execute(f"SELECT {keys}{selected} FROM {table} WHERE {built} AND ({where}){group}",
                                params).fetchall():
            count = row["star_count" if stars else "point_count"]
            if not count:
                continue
            key = (row["shell_index"], row["band_index"], row["wedge_index"]) if by_address else None
            _add(found.setdefault(key, _zero()), row, names)
            if by_address and count > 1:
                many.append(key)
        for start in range(0, len(many), 500):
            chunk = many[start:start + 500]
            rows = conn.execute(
                f"SELECT shell_index, band_index, wedge_index, {values} FROM {table} WHERE {built} AND (shell_index, band_index, wedge_index) IN ("
                + ", ".join("(?, ?, ?)" for _ in chunk) + ") ORDER BY id", [value for cell in chunk for value in cell]).fetchall()
            by_cell = {}
            for row in rows:
                by_cell.setdefault((row["shell_index"], row["band_index"], row["wedge_index"]), []).append(
                    (row["luminosity_w"] if stars else row["map_luminosity_sol"], row["position_x_mpc"], row["position_y_mpc"], row["position_z_mpc"],
                     row["temperature_k"] if stars else None))
            for key, cell in by_cell.items():
                exact = cell_sums(stars, cell)
                into = found[key]
                for index, value in enumerate(exact):
                    if SUMS[index] in names:
                        into[index] = value
    return found if by_address else found.get(None)


def combine(parts):
    """The sums of a block from the sums of the things it holds (its cells, or its child blocks): counts added, every other sum made with `math.fsum`, so the block's
    row is the same bits whichever order they come in and whichever path (a rebuild, a roll or a refresh) made them (GEN.1476)."""
    parts = list(parts)
    return [sum(int(part[index]) for part in parts) if index in _COUNT_AT else math.fsum(part[index] for part in parts) for index in range(len(SUMS))]


_COUNT_AT = (0, 6)


def _row(block, sums):
    return (block.m, block.shell, block.band, block.wedge,
            int(sums[0]), float(sums[1]), float(sums[2]), float(sums[3]), float(sums[4]), float(sums[5]),
            int(sums[6]), float(sums[7]), float(sums[8]), float(sums[9]), float(sums[10]))


def block_sums(conn, block):
    """The sums of the pre-placed objects inside one block, or `None` for a block holding none. A level-3 block
    combines the cells of its sectors; a bigger one combines the stored rows of its children (`combine`)."""
    if block.m == DRILL_STEP:
        sectors = drill_block_sectors(block)
        if not sectors:
            return None
        where, params = query._sectors_where(sectors)
        cells = _object_sums(conn, where, params, by_address=True)
        return combine(cells.values()) if cells else None
    children = [(c.shell, c.band, c.wedge) for _shell, group in drill_children(block) for c in group]
    held = []
    for start in range(0, len(children), 500):
        chunk = children[start:start + 500]
        held.extend([row[name] for name in SUMS] for row in conn.execute(
            f"SELECT {', '.join(SUMS)} FROM galaxy_block_lights WHERE m = ? AND (shell_index, band_index,"
            " wedge_index) IN (" + ", ".join("(?, ?, ?)" for _ in chunk) + ")",
            [block.m // DRILL_STEP] + [value for cell in chunk for value in cell]).fetchall())
    return combine(held) if held else None


def refresh_blocks(conn, blocks):
    """Recomputes the stored row of each block (a row is dropped when its block holds no pre-placed object),
    smallest blocks first because a block above level 3 adds up its children's rows. Returns the rows written."""
    written = 0
    for block in sorted(set(blocks), key=_order):
        conn.execute("DELETE FROM galaxy_block_lights WHERE m = ? AND shell_index = ? AND band_index = ?"
                     " AND wedge_index = ?", (block.m, block.shell, block.band, block.wedge))
        sums = block_sums(conn, block)
        if sums is not None:
            conn.execute(INSERT, _row(block, sums))
            written += 1
    return written


def refresh_chain(conn, address):
    """Recomputes the blocks above one sector, `(shell, band, wedge)`; nothing for a sector off the grid."""
    if address is None or any(value is None for value in address):
        return 0
    sync_cells(conn, [tuple(int(value) for value in address)])
    return refresh_blocks(conn, chain_blocks(tuple(int(value) for value in address)))


def refresh_sectors(conn, addresses, on_progress=None):
    """Recomputes the blocks above each sector in `addresses`, every block once. A set touching more than
    `FULL_REBUILD_ABOVE` blocks, or a table with no rows, is a `rebuild`. Returns the rows written."""
    blocks = {block for address in addresses for block in chain_blocks(tuple(address))}
    has_rows = conn.execute("SELECT 1 FROM galaxy_block_lights LIMIT 1").fetchone() is not None
    if len(blocks) > FULL_REBUILD_ABOVE or (blocks and not has_rows):
        return rebuild(conn, on_progress, parts=True)
    sync_cells(conn, [tuple(address) for address in addresses])
    written = 0
    for done, block in enumerate(sorted(blocks, key=_order), start=1):
        written += refresh_blocks(conn, [block])
        if on_progress is not None:
            on_progress(done, len(blocks))
    return written


def _not_a_number(value):
    raise TypeError(f"a block light row holds only ints and floats, not {value!r}")


def _write_rows(conn, rows):
    """Inserts `rows` (`_row` tuples, plain ints and floats) as multi-row INSERTs of `_WRITE_ROWS` rows. The text is
    made by `store._row_literal`, not by the driver's per-value dispatch, which was a third of a unit's time at
    2.5 million rows per unit (PERF.190)."""
    for start in range(0, len(rows), _WRITE_ROWS):
        conn.execute_literal(_INSERT_HEAD + ",".join(
            store._row_literal(_not_a_number, None, row) for row in rows[start:start + _WRITE_ROWS]))


def _unit_of(shell):
    """The unit a sector shell is rebuilt in: the shell of the level-`UNIT_LEVEL` block holding it."""
    return shell // UNIT_LEVEL


def _part_row(sums):
    """A sums list as the plain ints and floats a parts row holds."""
    return [int(sums[0]), float(sums[1]), float(sums[2]), float(sums[3]), float(sums[4]), float(sums[5]),
            int(sums[6]), float(sums[7]), float(sums[8]), float(sums[9]), float(sums[10])]


def _write_parts(conn, rows, add=False):
    """Writes parts rows `(shell, band, wedge, *sums)`; with `add`, a cell that has a row has the sums added to it (the scatter writers). The
    text is made as `_write_rows` makes it (PERF.190)."""
    tail = _PARTS_ADD if add else ""
    for start in range(0, len(rows), _WRITE_ROWS):
        conn.execute_literal(_PARTS_HEAD + ",".join(
            store._row_literal(_not_a_number, None, row) for row in rows[start:start + _WRITE_ROWS]) + tail)


def clear_parts(conn):
    """Empties the parts."""
    conn.execute("TRUNCATE TABLE galaxy_block_light_parts")


def forget_parts(conn, kind):
    """The scatter of `kind` (`stars` or `points`) was cleared: its sums leave the parts, and a cell left with no object leaves with them."""
    names = SUMS[:6] if kind == "stars" else SUMS[6:]
    conn.execute("UPDATE galaxy_block_light_parts SET " + ", ".join(f"{name} = 0" for name in names))
    conn.execute("DELETE FROM galaxy_block_light_parts WHERE star_count = 0 AND point_count = 0")


def sync_cells(conn, addresses):
    """Makes the parts of each sector cell in `addresses` what the object tables hold for it now: a sector was built, or deleted, or its
    objects pruned, so its unbuilt objects are no longer the ones the scatter wrote (GEN.1476)."""
    for address in addresses:
        shell, band, wedge = (int(value) for value in address)
        conn.execute("DELETE FROM galaxy_block_light_parts WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", (shell, band, wedge))
        found = _object_sums(conn, "shell_index = ? AND band_index = ? AND wedge_index = ?", (shell, band, wedge), by_address=True)
        _write_parts(conn, [tuple(cell) + tuple(_part_row(sums)) for cell, sums in sorted(found.items())])


def _part_shells(conn):
    """Every sector shell that holds a part, ascending."""
    return [int(row["shell_index"]) for row in conn.execute("SELECT DISTINCT shell_index FROM galaxy_block_light_parts ORDER BY shell_index").fetchall()]


def _object_shells(conn):
    """Every sector shell holding a pre-placed object, ascending."""
    shells = set()
    for table, built in (("bright_stars", _STARS_WHERE), ("phenomenon_scatter", _POINTS_WHERE)):
        shells.update(int(row["shell_index"]) for row in conn.execute(
            f"SELECT DISTINCT shell_index FROM {table} WHERE {built}").fetchall())
    return sorted(shells)


def _roll(cells, address, sums, parents):
    """Notes the sums of the sector cell at `address` in its level-`DRILL_STEP` block in `cells` (`{m: {DrillBlock: [sums]}}`); `parents` memoises a block's
    parent. A block combines what was noted in it when it is complete (`_complete`), and that is noted in its parent."""
    level3 = drill_parent((1,) + tuple(address))
    cells[level3.m].setdefault(level3, []).append(sums)


def _complete(cells, level):
    """Takes the rows of the blocks of `level` out of `cells`, sorted, each block the `combine` of what it holds, and notes each block in its parent's list; every
    one is complete because it is called when the shells run past the block shell they all belong to, smallest level first (PERF.187, GEN.1476)."""
    blocks = cells[level]
    rows = []
    for block, held in sorted(blocks.items(), key=lambda item: _order(item[0])):
        sums = combine(held)
        rows.append(_row(block, sums))
        parent = drill_parent(block)
        if parent is not None and parent.m in cells:
            cells[parent.m].setdefault(parent, []).append(sums)
    blocks.clear()
    return rows


def _part_sums(conn, shell):
    """The sums of the cells of one shell as the parts hold them: `{(shell, band, wedge): sums}`, as `_object_sums(by_address=True)` makes them."""
    return {(shell, row["band_index"], row["wedge_index"]): [row[name] for name in SUMS] for row in conn.execute(
        f"SELECT band_index, wedge_index, {', '.join(SUMS)} FROM galaxy_block_light_parts WHERE shell_index = ?", (shell,)).fetchall()}


def _unit_batches(conn, shells, on_shell=None, from_parts=False, on_found=None):
    """The rows of levels `DRILL_STEP` to `UNIT_LEVEL` for the sector `shells` of one unit, a list of rows at a time
    as their blocks are complete. A shell is read on its own, and a block of size `m` holds the sector shells
    `shell * m` to `shell * m + m - 1`, so once the shells pass a block shell its blocks are written and forgotten:
    what is kept is the blocks of the current block shell of each level, not every block of the unit's millions of
    sectors (PERF.187). `on_shell(shell)` is called before each shell is read. With `from_parts` the cells' sums are the parts' (GEN.1476), not
    a grouped read of the two object tables; `on_found(shell, found)` gets each shell's `{address: sums}` (the rebuild writes the parts with it)."""
    cells = {m: {} for m in sorted(LEVELS) if m <= UNIT_LEVEL}
    parents = {}
    previous = None
    for shell in shells:
        if previous is not None:
            for level in cells:
                if shell // level != previous // level and cells[level]:
                    yield _complete(cells, level)
        if on_shell is not None:
            on_shell(shell)
        parents.clear()
        found = _part_sums(conn, shell) if from_parts else _object_sums(conn, "shell_index = ?", (shell,), by_address=True)
        if on_found is not None:
            on_found(shell, found)
        for address, sums in found.items():
            _roll(cells, address, sums, parents)
        previous = shell
    for level in cells:
        if cells[level]:
            yield _complete(cells, level)


def _unit_rows(conn, shells):
    """All the rows of one unit (see `_unit_batches`), for callers that want them together."""
    return [row for batch in _unit_batches(conn, shells) for row in batch]


ABOVE_ROLL_SHARE = 0.8
"""float: The share of `_above_unit`'s progress that adding up the stored rows takes (the rest is writing the rows made)."""

ABOVE_REPORT_ROWS = 20000
"""int: Rows of the unit level added up between two reports of `_above_unit`."""


def _above_unit(conn, on_progress=None):
    """
    Writes the levels above `UNIT_LEVEL` from the stored rows of that level. Returns the rows written.

    Args:
        conn (Connection): A writable connection.
        on_progress (callable, optional): `on_progress(fraction)` from 0 to 1 while the rows are added up and written, so the long
            roll-up of a full galaxy shows on a bar (PERF.251).
    """
    rows_in = conn.execute(f"SELECT m, shell_index, band_index, wedge_index, {', '.join(SUMS)} FROM"
                           " galaxy_block_lights WHERE m = ?", (UNIT_LEVEL,)).fetchall()
    current = {}
    for done, row in enumerate(rows_in, start=1):
        current[DrillBlock(UNIT_LEVEL, row["shell_index"], row["band_index"], row["wedge_index"])] = [row[name] for name in SUMS]
        if on_progress is not None and done % ABOVE_REPORT_ROWS == 0:
            on_progress(ABOVE_ROLL_SHARE * done / len(rows_in))
    rows = []
    while current:      # each level is the `combine` of its children's rows, the order `block_sums` makes a refresh in (GEN.1476)
        held = {}
        for block, sums in current.items():
            parent = drill_parent(block)
            if parent is not None:
                held.setdefault(parent, []).append(sums)
        current = {block: combine(children) for block, children in sorted(held.items(), key=lambda item: _order(item[0]))}
        rows.extend(_row(block, sums) for block, sums in current.items())
    for start in range(0, len(rows), _WRITE_ROWS):
        _write_rows(conn, rows[start:start + _WRITE_ROWS])
        if on_progress is not None:
            on_progress(ABOVE_ROLL_SHARE + (1 - ABOVE_ROLL_SHARE) * min(1.0, (start + _WRITE_ROWS) / len(rows)))
    return len(rows)


def _object_shells_of(conn, from_parts):
    return _part_shells(conn) if from_parts else _object_shells(conn)


def rebuild(conn, on_progress=None, on_phase=None, on_above=None, parts=False, from_parts=False):
    """
    Computes every row from `bright_stars` and `phenomenon_scatter`, replacing the table's rows: the sector shells in
    units (those of one level-`UNIT_LEVEL` block shell), each unit's rows written and committed together so a stopped
    run keeps what it did and a second run replaces it, then the levels above the units added up from the stored rows.
    `on_progress(done, total)` is called after each unit, and `on_phase(text)` before the shell scan, after it and
    before the levels above the units (PERF.189), so a long rebuild never goes quiet. `on_above(fraction)` is called
    from 0 to 1 while the levels above the units are added up and written (`_above_unit`, PERF.251). Returns the rows written.

    With `parts` the cells' sums it finds are also written as the parts (`galaxy_block_light_parts`, GEN.1476), replacing them, so a rebuild
    is the repair of the parts. With `from_parts` the sums are read from the parts and not from the object tables (`roll`); it is the
    caller's to have checked them (`parts_current`). The migration that made this table runs on a schema without the parts: neither is set there.
    """
    if on_phase is not None:
        on_phase("Looking for the sector shells that hold pre-placed objects...")
    shells = _object_shells_of(conn, from_parts)
    units = {}
    for shell in shells:
        units.setdefault(_unit_of(shell), []).append(shell)
    if on_phase is not None:
        on_phase(f"{len(shells):,} sector shells hold pre-placed objects, worked in {len(units):,} units.")
    clear(conn)
    if parts:
        clear_parts(conn)
    conn.commit()
    written = 0

    def note_parts(shell, found):
        _write_parts(conn, [tuple(address) + tuple(_part_row(sums)) for address, sums in sorted(found.items())])

    for done, unit in enumerate(sorted(units), start=1):
        for batch in _unit_batches(conn, units[unit], from_parts=from_parts, on_found=note_parts if parts else None):
            _write_rows(conn, batch)
            written += len(batch)
        conn.commit()
        if on_progress is not None:
            on_progress(done, len(units))
    if on_phase is not None:
        on_phase("Adding up the levels above the units...")
    written += _above_unit(conn, on_above)
    conn.commit()
    return written


def parts_current(conn):
    """Whether the parts hold every unbuilt object once: their counts equal the tables' unbuilt counts. A galaxy scattered before the parts
    existed, a prune or a clear that did not reach them, or a scatter written by a writer that did not note its rows, fails this, and `roll` then
    rebuilds (GEN.1476)."""
    held = conn.execute("SELECT COALESCE(SUM(star_count), 0) AS stars, COALESCE(SUM(point_count), 0) AS points FROM galaxy_block_light_parts").fetchone()
    stars = conn.execute(f"SELECT COUNT(*) AS n FROM bright_stars WHERE {_STARS_WHERE}").fetchone()["n"]
    points = conn.execute(f"SELECT COUNT(*) AS n FROM phenomenon_scatter WHERE {_POINTS_WHERE}").fetchone()["n"]
    return int(held["stars"]) == stars and int(held["points"]) == points


def roll(conn, on_progress=None, on_phase=None, on_above=None):
    """
    The block light rows made from the parts, the stage a scatter ends with (GEN.1476): the same rows `rebuild` makes, without reading the
    two object tables, which at full size is the hour and the memory the stage cost. When `parts_current` fails it is a `rebuild` that also
    writes the parts again. Returns `(rows written, whether the parts were used)`.
    """
    if parts_current(conn):
        return rebuild(conn, on_progress, on_phase, on_above, from_parts=True), True
    if on_phase is not None:
        on_phase("The stored parts do not match the scatter (a galaxy scattered before them, or a prune): making them again from the objects.")
    return rebuild(conn, on_progress, on_phase, on_above, parts=True), False
