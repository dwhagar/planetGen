# planetgen/db/region_stats.py

"""
Stored aggregates of the named regions (DB.107, `region_stats`; `galaxy/named_regions.py`).

A region is a 9-block of the Galaxy Map's block ladder that holds at least one generated sector. A row is one such cell
`(shell_index, band_index, wedge_index)` (the level, 9, is implied) with:

* `sector_count`: its generated sectors (a row in `sectors`), and `capacity`: the sectors the cell's outline holds, so the
  share charted is `sector_count / capacity` (a standard fill of 109 sectors in a cell of about 780 is about 14 percent);
* the extent the map draws and frames the region with (`named_regions.region_extent`): the centre `extent_x_pc`,
  `extent_y_pc`, `extent_z_pc` and `extent_radius_pc` of the bounding sphere of the largest face-connected clump of its
  generated sectors, never under the framing sphere, and `clump_count`, how many such clumps the cell holds.

A cell with no generated sector has no row, and the region does not exist. The rows only change when a sector is saved or
deleted, so they are kept by the hooks that keep `galaxy_block_stats` (`db/block_stats.py` calls `refresh_chain`,
`refresh_sectors` and `rebuild` here): every row is the same whatever order the sectors were generated or deleted in
(`tests/test_order_independence.py`).

`region_of` and `charted_percent` read the table.
"""

import functools

from planetgen.db import query
from planetgen.galaxy import drill, named_regions, spherical

LEVEL = named_regions.REGION_LEVEL
"""int: The block size of a region."""

COLUMNS = ("shell_index", "band_index", "wedge_index", "sector_count", "capacity", "clump_count", "extent_x_pc",
           "extent_y_pc", "extent_z_pc", "extent_radius_pc")
"""tuple: The table's columns in the order `_row` writes them."""

INSERT = f"INSERT INTO region_stats ({', '.join(COLUMNS)}) VALUES ({', '.join('?' for _ in COLUMNS)})"

CHUNK = 300
"""int: Cells read per statement."""


def region_cell(address):
    """The region (9-block) that holds the sector `(shell, band, wedge)`."""
    shell, band, wedge = (int(value) for value in address)
    return drill.drill_chain_of(shell, band, wedge)[drill.DRILL_LEVELS.index(LEVEL)]


@functools.lru_cache(maxsize=256)
def cell_sectors(cell):
    """Every sector the outline of the region `cell` holds, as `DrillBlock(1, shell, band, wedge)`."""
    return tuple(drill.DrillBlock(1, *address) for address in drill.block_addresses(cell))


def capacity_of(cell):
    """How many sectors the outline of the region `cell` holds."""
    return len(cell_sectors(cell))


def generated_sectors(conn, cell):
    """The `(shell, band, wedge)` of the generated sectors of the region `cell`, sorted."""
    where, params = query._sectors_where(cell_sectors(cell))
    rows = conn.execute(f"SELECT shell_index, band_index, wedge_index FROM sectors WHERE {where}", params).fetchall()
    return sorted((int(r["shell_index"]), int(r["band_index"]), int(r["wedge_index"])) for r in rows)


def _row(cell, addresses, edge_pc):
    extent = named_regions.region_extent(cell, addresses, edge_pc)
    x, y, z = extent.centre_pc
    return (cell.shell, cell.band, cell.wedge, len(addresses), capacity_of(cell), extent.clumps, x, y, z, extent.radius_pc)


def refresh_cells(conn, cells, edge_pc=None):
    """Recomputes the row of each region cell (a row is dropped when its cell no longer holds a generated sector).
    Returns how many rows were written."""
    edge_pc = query.galaxy_edge_pc(conn) if edge_pc is None else edge_pc
    written = 0
    for cell in sorted(set(cells), key=tuple):
        conn.execute("DELETE FROM region_stats WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
                     (cell.shell, cell.band, cell.wedge))
        addresses = generated_sectors(conn, cell)
        if addresses:
            conn.execute(INSERT, _row(cell, addresses, edge_pc))
            written += 1
    return written


def refresh_chain(conn, address):
    """Puts right the region that holds one sector, `(shell, band, wedge)`; nothing for a sector off the grid."""
    if address is None or any(value is None for value in address):
        return 0
    return refresh_cells(conn, [region_cell(address)])


def refresh_sectors(conn, addresses, on_progress=None):
    """Puts right the regions that hold each sector in `addresses`, every region once. Returns the rows written."""
    cells = {region_cell(address) for address in addresses}
    written = 0
    for done, cell in enumerate(sorted(cells, key=tuple), start=1):
        written += refresh_cells(conn, [cell])
        if on_progress is not None:
            on_progress(done, len(cells))
    return written


def rebuild(conn, on_progress=None):
    """
    Makes every row again from `sectors`, a shell at a time: the table emptied, then a row for each region cell that holds a
    generated sector. A cell is written when the shells of its block have all been read (a region's block shell only ever
    grows with the sector shell). Returns the rows written.
    """
    shells = [row["shell_index"] for row in conn.execute(
        "SELECT DISTINCT shell_index FROM sectors WHERE shell_index IS NOT NULL ORDER BY shell_index").fetchall()]
    conn.execute("DELETE FROM region_stats")
    edge_pc = query.galaxy_edge_pc(conn)
    cells, group, written = {}, None, 0

    def flush():
        nonlocal written
        for cell in sorted(cells, key=tuple):
            conn.execute(INSERT, _row(cell, sorted(cells[cell]), edge_pc))
            written += 1
        cells.clear()

    for done, shell in enumerate(shells, start=1):
        for r in conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors WHERE shell_index = ?",
                              (shell,)).fetchall():
            address = (int(r["shell_index"]), int(r["band_index"]), int(r["wedge_index"]))
            cell = region_cell(address)
            if group is not None and cell.shell != group:
                flush()
            group = cell.shell
            cells.setdefault(cell, []).append(address)
        if on_progress is not None:
            on_progress(done, len(shells))
    flush()
    return written


def _region(row):
    cell = drill.DrillBlock(LEVEL, int(row["shell_index"]), int(row["band_index"]), int(row["wedge_index"]))
    count, capacity = int(row["sector_count"]), int(row["capacity"])
    return {"block": cell, "sector_count": count, "capacity": capacity, "charted_percent": 100.0 * count / capacity,
            "clump_count": int(row["clump_count"]),
            "extent_centre_pc": (float(row["extent_x_pc"]), float(row["extent_y_pc"]), float(row["extent_z_pc"])),
            "extent_radius_pc": float(row["extent_radius_pc"])}


def region_row(conn, cell):
    """The stored region of the 9-block `cell` (`block`, `sector_count`, `capacity`, `charted_percent`, `clump_count`,
    `extent_centre_pc`, `extent_radius_pc`), or `None` when it holds no generated sector."""
    row = conn.execute(f"SELECT {', '.join(COLUMNS)} FROM region_stats WHERE shell_index = ? AND band_index = ?"
                       " AND wedge_index = ?", (cell.shell, cell.band, cell.wedge)).fetchone()
    return None if row is None else _region(row)


def region_of(conn, position_pc, edge_pc=None):
    """The region that holds the point `position_pc` (galaxy frame, parsecs): the stored region of the 9-block of the
    sector that contains it, or `None` when that cell holds no generated sector."""
    edge_pc = query.galaxy_edge_pc(conn) if edge_pc is None else edge_pc
    return region_row(conn, region_cell(spherical.sector_address_at(tuple(position_pc), edge_pc)))


def charted_percent(conn, cell):
    """The share of the 9-block `cell` that is generated, in percent: 0 for a cell with no row (no region)."""
    region = region_row(conn, cell)
    return 0.0 if region is None else region["charted_percent"]


def regions(conn, cells=None):
    """The stored regions, every one or just those of `cells` (9-blocks), each as `region_row` gives it, sorted by share
    charted (largest first), ties by the cell."""
    rows = []
    if cells is None:
        rows = conn.execute(f"SELECT {', '.join(COLUMNS)} FROM region_stats").fetchall()
    else:
        cells = list(cells)
        for start in range(0, len(cells), CHUNK):
            chunk = cells[start:start + CHUNK]
            rows += conn.execute(
                f"SELECT {', '.join(COLUMNS)} FROM region_stats WHERE (shell_index, band_index, wedge_index) IN ("
                + ", ".join("(?, ?, ?)" for _ in chunk) + ")", [v for c in chunk for v in (c.shell, c.band, c.wedge)]
            ).fetchall()
    found = [_region(row) for row in rows]
    return sorted(found, key=lambda region: (-region["charted_percent"], tuple(region["block"])))
