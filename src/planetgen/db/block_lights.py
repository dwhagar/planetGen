# planetgen/db/block_lights.py

"""
Stored light sums for the Galaxy Map's sector-block nodes (MAP.259, `galaxy_block_lights`).

The single-space map shows a block it has not opened as one node: the brightest objects inside it and, for
the rest, a glow made from their summed light. A block holds sectors that were generated (their sums are
`galaxy_block_stats`) and, far more of them, sectors nothing has been built in, whose only stored objects are
the pre-placed bright stars (`bright_stars` with no system built around them) and the scattered phenomena
(`phenomenon_scatter` not built). This table sums those two kinds per block of every level above a sector
(`drill.DRILL_LEVELS[:-1]`), so a node's glow and its centre of light are one row read.

A row is one block `(m, shell_index, band_index, wedge_index)` holding at least one pre-placed object. Columns:

* `star_count`, `star_luminosity_sol` (summed luminosity), `star_lpx`, `star_lpy`, `star_lpz` (luminosity times
  position, so the centre of light is the sum over `star_luminosity_sol`; solar luminosities times parsecs) and
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
  `planetgen fill-block-totals`, the migration, and what a scatter does at its end.
"""

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


def _add(into, row, names):
    for name in names:
        into[SUMS.index(name)] += row[name] or 0


def _object_sums(conn, where, params, by_address=False):
    """The sums of the unbuilt objects `where` picks (`params` its parameters): `{(shell, band, wedge): sums}` with
    `by_address`, else one sums list."""
    group = " GROUP BY shell_index, band_index, wedge_index" if by_address else ""
    keys = "shell_index, band_index, wedge_index, " if by_address else ""
    found = {}
    for table, selected, built, names in (
            ("bright_stars", _STAR_SUMS, _STARS_WHERE, SUMS[:6]),
            ("phenomenon_scatter", _POINT_SUMS, _POINTS_WHERE, SUMS[6:])):
        for row in conn.execute(f"SELECT {keys}{selected} FROM {table} WHERE {built} AND ({where}){group}",
                                params).fetchall():
            if not row["star_count" if table == "bright_stars" else "point_count"]:
                continue
            key = (row["shell_index"], row["band_index"], row["wedge_index"]) if by_address else None
            _add(found.setdefault(key, _zero()), row, names)
    return found if by_address else found.get(None)


def _row(block, sums):
    return (block.m, block.shell, block.band, block.wedge,
            int(sums[0]), float(sums[1]), float(sums[2]), float(sums[3]), float(sums[4]), float(sums[5]),
            int(sums[6]), float(sums[7]), float(sums[8]), float(sums[9]), float(sums[10]))


def block_sums(conn, block):
    """The sums of the pre-placed objects inside one block, or `None` for a block holding none. A level-3 block
    adds up the objects of its sectors; a bigger one adds up the stored rows of its children."""
    if block.m == DRILL_STEP:
        sectors = drill_block_sectors(block)
        if not sectors:
            return None
        where, params = query._sectors_where(sectors)
        return _object_sums(conn, where, params)
    children = [(c.shell, c.band, c.wedge) for _shell, group in drill_children(block) for c in group]
    total = None
    for start in range(0, len(children), 500):
        chunk = children[start:start + 500]
        for row in conn.execute(
                f"SELECT {', '.join(SUMS)} FROM galaxy_block_lights WHERE m = ? AND (shell_index, band_index,"
                " wedge_index) IN (" + ", ".join("(?, ?, ?)" for _ in chunk) + ")",
                [block.m // DRILL_STEP] + [value for cell in chunk for value in cell]).fetchall():
            total = _zero() if total is None else total
            _add(total, row, SUMS)
    return total


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
    return refresh_blocks(conn, chain_blocks(tuple(int(value) for value in address)))


def refresh_sectors(conn, addresses, on_progress=None):
    """Recomputes the blocks above each sector in `addresses`, every block once. A set touching more than
    `FULL_REBUILD_ABOVE` blocks, or a table with no rows, is a `rebuild`. Returns the rows written."""
    blocks = {block for address in addresses for block in chain_blocks(tuple(address))}
    has_rows = conn.execute("SELECT 1 FROM galaxy_block_lights LIMIT 1").fetchone() is not None
    if len(blocks) > FULL_REBUILD_ABOVE or (blocks and not has_rows):
        return rebuild(conn, on_progress)
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


def _object_shells(conn):
    """Every sector shell holding a pre-placed object, ascending."""
    shells = set()
    for table, built in (("bright_stars", _STARS_WHERE), ("phenomenon_scatter", _POINTS_WHERE)):
        shells.update(int(row["shell_index"]) for row in conn.execute(
            f"SELECT DISTINCT shell_index FROM {table} WHERE {built}").fetchall())
    return sorted(shells)


def _roll(cells, address, sums, parents):
    """Adds the sums of the sector at `address` to each block above it in `cells` (`{m: {DrillBlock: sums}}`) up to
    the unit's level; `parents` memoises a level-3 block's chain, which most sectors of a row share."""
    level3 = drill_parent((1,) + tuple(address))
    chain = parents.get(level3)
    if chain is None:
        chain = parents[level3] = [block for block in chain_blocks(address) if block.m <= UNIT_LEVEL]
    for block in chain:
        into = cells[block.m].setdefault(block, _zero())
        for index, value in enumerate(sums):
            into[index] += value


def _complete(cells, level):
    """Takes the rows of the blocks of `level` out of `cells`, sorted; every one is complete because it is called
    when the shells run past the block shell they all belong to (PERF.187)."""
    blocks = cells[level]
    rows = [_row(block, sums) for block, sums in sorted(blocks.items(), key=lambda item: _order(item[0]))]
    blocks.clear()
    return rows


def _unit_batches(conn, shells, on_shell=None):
    """The rows of levels `DRILL_STEP` to `UNIT_LEVEL` for the sector `shells` of one unit, a list of rows at a time
    as their blocks are complete. A shell is read on its own, and a block of size `m` holds the sector shells
    `shell * m` to `shell * m + m - 1`, so once the shells pass a block shell its blocks are written and forgotten:
    what is kept is the blocks of the current block shell of each level, not every block of the unit's millions of
    sectors (PERF.187). `on_shell(shell)` is called before each shell is read."""
    cells = {m: {} for m in LEVELS if m <= UNIT_LEVEL}
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
        found = _object_sums(conn, "shell_index = ?", (shell,), by_address=True)
        for address, sums in found.items():
            _roll(cells, address, sums, parents)
        previous = shell
    for level in cells:
        if cells[level]:
            yield _complete(cells, level)


def _unit_rows(conn, shells):
    """All the rows of one unit (see `_unit_batches`), for callers that want them together."""
    return [row for batch in _unit_batches(conn, shells) for row in batch]


def _above_unit(conn):
    """Writes the levels above `UNIT_LEVEL` from the stored rows of that level. Returns the rows written."""
    cells = {}
    for row in conn.execute(f"SELECT m, shell_index, band_index, wedge_index, {', '.join(SUMS)} FROM"
                            " galaxy_block_lights WHERE m = ?", (UNIT_LEVEL,)).fetchall():
        block = DrillBlock(UNIT_LEVEL, row["shell_index"], row["band_index"], row["wedge_index"])
        sums = [row[name] for name in SUMS]
        parent = drill_parent(block)
        while parent is not None:
            into = cells.setdefault(parent, _zero())
            for index, value in enumerate(sums):
                into[index] += value
            parent = drill_parent(parent)
    rows = [_row(block, sums) for block, sums in sorted(cells.items(), key=lambda item: _order(item[0]))]
    _write_rows(conn, rows)
    return len(rows)


def rebuild(conn, on_progress=None, on_phase=None):
    """
    Computes every row from `bright_stars` and `phenomenon_scatter`, replacing the table's rows: the sector shells in
    units (those of one level-`UNIT_LEVEL` block shell), each unit's rows written and committed together so a stopped
    run keeps what it did and a second run replaces it, then the levels above the unit added up from the stored rows.
    `on_progress(done, total)` is called after each unit, and `on_phase(text)` before the shell scan, after it and
    before the levels above the units (PERF.189), so a long rebuild never goes quiet. Returns the rows written.
    """
    if on_phase is not None:
        on_phase("Looking for the sector shells that hold pre-placed objects...")
    shells = _object_shells(conn)
    units = {}
    for shell in shells:
        units.setdefault(_unit_of(shell), []).append(shell)
    if on_phase is not None:
        on_phase(f"{len(shells):,} sector shells hold pre-placed objects, worked in {len(units):,} units.")
    clear(conn)
    conn.commit()
    written = 0
    for done, unit in enumerate(sorted(units), start=1):
        for batch in _unit_batches(conn, units[unit]):
            _write_rows(conn, batch)
            written += len(batch)
        conn.commit()
        if on_progress is not None:
            on_progress(done, len(units))
    if on_phase is not None:
        on_phase("Adding up the levels above the units...")
    written += _above_unit(conn)
    conn.commit()
    return written
