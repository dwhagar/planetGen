# planetgen/db/region_stats.py

"""
Stored aggregates of the named regions (DB.107, `region_stats`; `galaxy/named_regions.py`).

A region is a 9-block of the Galaxy Map's block ladder that holds at least one generated sector. A row is one such cell
`(shell_index, band_index, wedge_index)` (the level, 9, is implied) with:

* `sector_count`: its generated sectors (a row in `sectors`), and `capacity`: the sectors the cell's outline holds, so the
  share charted is `sector_count / capacity` (a standard fill of 109 sectors in a cell of about 780 is about 14 percent);
* the extent the map draws and frames the region with (`named_regions.region_extent`): the centre `extent_x_pc`,
  `extent_y_pc`, `extent_z_pc` and `extent_radius_pc` of the bounding sphere of the largest face-connected clump of its
  generated sectors, never under the framing sphere, and `clump_count`, how many such clumps the cell holds;
* `mass_sum_kg` (DB.108): the sum of its generated sectors' masses (`sector_stats.total_mass_kg`, `db/sector_mass.py`), added
  exactly (`math.fsum`) so the same sectors give the same bits whatever order they were generated in; a sector that could not
  be summed counts as nothing.

A cell with no generated sector has no row, and the region does not exist. The rows only change when a sector is saved or
deleted, so they are kept by the hooks that keep `galaxy_block_stats` (`db/block_stats.py` calls `refresh_chain`,
`refresh_sectors` and `rebuild` here): every row is the same whatever order the sectors were generated or deleted in
(`tests/test_order_independence.py`).

`region_of` and `charted_percent` read the table.
"""

import functools
import math

from planetgen.db import interests, query
from planetgen.galaxy import drill, named_regions, spherical
from planetgen.physics import habitability

LEVEL = named_regions.REGION_LEVEL
"""int: The block size of a region."""

COLUMNS = ("shell_index", "band_index", "wedge_index", "sector_count", "capacity", "clump_count", "extent_x_pc",
           "extent_y_pc", "extent_z_pc", "extent_radius_pc", "mass_sum_kg")
"""tuple: The table's columns in the order `_row` writes them."""

INSERT = f"INSERT INTO region_stats ({', '.join(COLUMNS)}) VALUES ({', '.join('?' for _ in COLUMNS)})"

CHUNK = 300
"""int: Cells read per statement."""

LABEL_LEVELS = (27, 81, 243)
"""tuple: The block sizes whose node lists the charted regions inside it (`node_regions`, MAP.314): the levels the Galaxy
Map names regions at from a distance, a family at 81 and 243 and the regions themselves at 27 (a region's own node, 9, carries
its name). Above 243 a node would list too many regions to read, and below 9 the region is the infobox's line (MAP.315)."""


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


def generated_mass_kg(conn, cell):
    """The summed mass, in kilograms, of the generated sectors of the region `cell` (`sector_stats.total_mass_kg`; 0 for a
    sector with no figure)."""
    where, params = query._sectors_where(cell_sectors(cell))
    rows = conn.execute(
        "SELECT st.total_mass_kg AS mass FROM (SELECT shell_index, band_index, wedge_index FROM sectors WHERE " + where + ") s"
        " LEFT JOIN sector_stats st ON st.shell_index = s.shell_index AND st.band_index = s.band_index"
        " AND st.wedge_index = s.wedge_index", params).fetchall()
    return math.fsum(float(r["mass"] or 0.0) for r in rows)


def _row(cell, addresses, edge_pc, mass_kg):
    extent = named_regions.region_extent(cell, addresses, edge_pc)
    x, y, z = extent.centre_pc
    return (cell.shell, cell.band, cell.wedge, len(addresses), capacity_of(cell), extent.clumps, x, y, z, extent.radius_pc,
            mass_kg)


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
            conn.execute(INSERT, _row(cell, addresses, edge_pc, generated_mass_kg(conn, cell)))
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
    cells, masses, group, written = {}, {}, None, 0

    def flush():
        nonlocal written
        for cell in sorted(cells, key=tuple):
            conn.execute(INSERT, _row(cell, sorted(cells[cell]), edge_pc, math.fsum(masses[cell])))
            written += 1
        cells.clear()
        masses.clear()

    for done, shell in enumerate(shells, start=1):
        for r in conn.execute(
                "SELECT s.shell_index, s.band_index, s.wedge_index, st.total_mass_kg FROM (SELECT * FROM sectors WHERE"
                " shell_index = ?) s LEFT JOIN sector_stats st ON st.shell_index = s.shell_index"
                " AND st.band_index = s.band_index AND st.wedge_index = s.wedge_index", (shell,)).fetchall():
            address = (int(r["shell_index"]), int(r["band_index"]), int(r["wedge_index"]))
            cell = region_cell(address)
            if group is not None and cell.shell != group:
                flush()
            group = cell.shell
            cells.setdefault(cell, []).append(address)
            masses.setdefault(cell, []).append(float(r["total_mass_kg"] or 0.0))
        if on_progress is not None:
            on_progress(done, len(shells))
    flush()
    return written


def _region(row):
    cell = drill.DrillBlock(LEVEL, int(row["shell_index"]), int(row["band_index"]), int(row["wedge_index"]))
    count, capacity = int(row["sector_count"]), int(row["capacity"])
    return {"block": cell, "sector_count": count, "capacity": capacity, "charted_percent": 100.0 * count / capacity,
            "clump_count": int(row["clump_count"]), "mass_sum_kg": float(row["mass_sum_kg"]),
            "extent_centre_pc": (float(row["extent_x_pc"]), float(row["extent_y_pc"]), float(row["extent_z_pc"])),
            "extent_radius_pc": float(row["extent_radius_pc"])}


def region_row(conn, cell):
    """The stored region of the 9-block `cell` (`block`, `sector_count`, `capacity`, `charted_percent`, `clump_count`,
    `extent_centre_pc`, `extent_radius_pc`, `mass_sum_kg`), or `None` when it holds no generated sector."""
    row = conn.execute(f"SELECT {', '.join(COLUMNS)} FROM region_stats WHERE shell_index = ? AND band_index = ?"
                       " AND wedge_index = ?", (cell.shell, cell.band, cell.wedge)).fetchone()
    return None if row is None else _region(row)


def extent_of(region):
    """The extent a charted region carries on the wire (MAP.319), from its stored row `region` (`region_row`): `centre`, the
    middle of its bounding sphere as `[x, y, z]` in parsecs, and `radius` in parsecs, which the map frames the region with."""
    return {"centre": list(region["extent_centre_pc"]), "radius": region["extent_radius_pc"]}


def charted_of(region):
    """How much of a charted region is generated, as the region's own node carries it (MAP.314): `sectors` and `percent`, the
    share of its outline's sectors in percent, which the map ranks its names by."""
    return {"sectors": region["sector_count"], "percent": round(region["charted_percent"], 3)}


def labels(cell, key, edge_pc):
    """
    The names a charted region carries on the wire (API.35): `name` ("Fardoo Beez Region") and `parent`, the name of its 81-block
    family ("Fardoo Reach"), under the naming key `key`; `None` when there is no key in force. Names come from the gated phoneme
    codec only. The caller has checked that the region is charted (`region_row`), which is what makes it exist.
    """
    if key is None:
        return None
    family = named_regions.family_of(cell)
    return {"name": named_regions.region_name(cell, key, edge_pc), "parent": named_regions.family_name(family, key, edge_pc)}


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


def cells_in(block):
    """The region cells (9-blocks) a block of a level above 9 holds, in the ladder's order: at most 27 for an 81-block."""
    if block.m < LEVEL:
        raise ValueError(f"a block of {block.m} sectors does not hold whole regions")
    if block.m == LEVEL:
        return [block]
    cells = []
    for _shell, children in drill.drill_children(block):
        for child in children:
            cells.extend(cells_in(child))
    return cells


def node_regions(conn, block, key, edge_pc):
    """
    The charted regions a node of `LABEL_LEVELS` lists on the wire (MAP.314), for the map to name them from a distance: a
    list of `{"key", "name", "parent", "centre", "radius", "sectors", "percent"}` (the region's 9-block key, its name and its
    family's, its bounding sphere as the middle `[x, y, z]` and the `radius`, in parsecs, its generated sectors and its share
    charted, in percent),
    most sectors first, ties by the cell. `None` for a block of another size, a block that holds no charted region, and while
    no naming key is in force (a region has no name then, as on its own node).
    """
    if block.m not in LABEL_LEVELS or key is None:
        return None
    found = regions(conn, cells_in(block))
    if not found:
        return None
    listed = []
    for region in sorted(found, key=lambda region: (-region["sector_count"], tuple(region["block"]))):
        names = labels(region["block"], key, edge_pc)
        listed.append({"key": drill.format_drill_key(region["block"]), "name": names["name"], "parent": names["parent"],
                       "centre": list(region["extent_centre_pc"]), "radius": region["extent_radius_pc"], **charted_of(region)})
    return listed


def inhabited_by_region(conn, edge_pc):
    """
    The inhabited worlds of each region (MAP.321 c): `{cell: {"systems", "worlds", "score"}}` for every region cell that holds a system
    with an `inhabited` interest (`db/interests.py`): how many such systems, how many inhabited worlds in them and the best PHI-6
    score among them (`None` when none is scored). Read from `system_interests` at the time of asking (colonies are few), each system
    put in the cell its position lies in, so it is never out of step with a colony that has just been founded.
    """
    found = {}
    for row in conn.execute(
            "SELECT n, best_score, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc FROM system_interests WHERE reason = ?",
            (interests.INHABITED,)).fetchall():
        cell = region_cell(spherical.sector_address_at(
            (row["galaxy_x_pc"], row["galaxy_y_pc"], row["galaxy_z_pc"]), edge_pc))
        entry = found.setdefault(cell, {"systems": 0, "worlds": 0, "score": None})
        entry["systems"] += 1
        entry["worlds"] += int(row["n"])
        if row["best_score"] is not None and (entry["score"] is None or float(row["best_score"]) > entry["score"]):
            entry["score"] = float(row["best_score"])
    return found


def charted_regions(conn, key, edge_pc):
    """
    The charted regions as the Galaxy Map's list and its galactic-scale nodes read them (MAP.316), each one `{key, m, shell,
    band, wedge, centre, radius, sector_count, capacity, share, clumps}` and, with a naming key in force, `name` and `parent`
    (`labels`): the cell's node key, the middle `[x, y, z]` and `radius` of its extent in parsecs, how many generated sectors it
    holds of how many its outline can, and `share`, that fraction (0 to 1, above 0 for every region). A region with an inhabited
    world also has `inhabited` (MAP.321 c): `systems`, `worlds` and, when one is scored, the best `score` and its PHI-6 `tier`.
    Sorted by share charted (largest first), ties by the cell.
    """
    found = []
    inhabited = inhabited_by_region(conn, edge_pc)
    for region in regions(conn):
        cell = region["block"]
        entry = {"key": drill.format_drill_key(cell), "m": cell.m, "shell": cell.shell, "band": cell.band, "wedge": cell.wedge,
                 "centre": list(region["extent_centre_pc"]), "radius": region["extent_radius_pc"],
                 "sector_count": region["sector_count"], "capacity": region["capacity"],
                 "share": region["sector_count"] / region["capacity"], "clumps": region["clump_count"]}
        entry.update(labels(cell, key, edge_pc) or {})
        if cell in inhabited:
            here = dict(inhabited[cell])
            if here["score"] is None:
                del here["score"]
            else:
                here["tier"] = habitability.tier(here["score"])
                here["score"] = round(here["score"], 3)
            entry["inhabited"] = here
        found.append(entry)
    return found
