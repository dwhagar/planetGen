# planetgen.db.query

"""
Read queries over the planetGen database (src/planetgen/db/schema.sql).

A thin read-only front end over the tables `generation/run_sector.py`/`generation/run_system.py`
populate, for questions like "every G-type system," "everything within 50
light-years of a given system," and "what sectors exist" -- see
docs/TODO.md's Phase 3 ("A way to list/query what's already stored").
Deliberately plain SQL rather than routing through `planetgen.db.store`'s
`load_star_system`/`load_sector` (Phase 2's read path): these are simple,
columnar listings, not full object-graph reconstructions, so a raw query
is the more direct tool for the job -- the read path remains what a
future richer tool (or a re-upload/re-render workflow) would build on.

This tool never writes -- `open_readonly` below connects with the same
`MySQLConfig` every other entry point uses, but the actual enforcement
that the connection can't write is a deployment concern: point
`PLANETGEN_MYSQL_USER`/`PLANETGEN_MYSQL_PASSWORD` at a database account
with `SELECT`-only grants for this tool rather than the read-write account
`planetgen` and the Flask app use (the app writes too: the admin pages
and the API's write endpoints, see `planetgen/api/config.py`) -- MySQL has no per-connection "open this read-only"
flag the way SQLite's `file:...?mode=ro` URI trick gave the old SQLite
version of this function, so the guarantee lives in the account's grants
instead of the connection itself.

The command-line front end is `planetgen.cli.query`.
"""

import datetime
import hashlib
import json
import math
import re
import threading
import time
import zlib

import pymysql

from planetgen.db.store import (
    STATEMENT_TIMEOUT_ERRORS, body_class_scores, body_living, compact_host_of_star, escape_like, get_connection, get_galaxy_bounds, get_galaxy_shape, sector_deletions, surrounding_cloud,
    with_gases,
)
from planetgen.db import compact_pairs, countcache, hostednebula, interests, nebula_blocks, point_columns, search_values
from planetgen.db.scatter_systems import SOURCES as SCATTER_SOURCES
from planetgen.physics import constants, relations
from planetgen.physics.stellar_evolution import phenomenon_map_luminosity_sol
from planetgen.physics import equipment, habitability, habitability_world, life_classes, life_scores, score_filters
from planetgen import tuning
from planetgen.population import facility_ratings, facility_types, tech
from planetgen.generation import cloud as cloudContext
from planetgen.generation.star import compressed_heliosphere_radius
from planetgen.generation.bright_stars import MPC_PER_PC, POPULATIONS as BRIGHT_STAR_POPULATIONS
from planetgen._version import __version__
from planetgen.util import generator_epoch
from planetgen.util.star_count import parse_star_count
from planetgen.galaxy import spherical
from planetgen.galaxy.density import shape_with_terms
from planetgen.galaxy.frame import binary_companion_offset_pc
from planetgen.galaxy.spherical import (
    galaxy_to_local_pc, local_to_galaxy_pc, neighbor_addresses,
    provisional_sector_designation, sector_cell_vertices_pc, sector_position_pc,
)

from planetgen.galaxy import enclosed_mass, keepout, mass_chain, shape_draw
from planetgen.galaxy import objectref as object_ref
from planetgen.galaxy import object_uid
from planetgen.names import naming_key
from planetgen.galaxy.sector import classify_octant
from planetgen.galaxy.window import check_width, parse_window_key, window_cells
from planetgen.galaxy.drill import (
    DRILL_STEP, STAGE_TOP, DrillBlock, drill_block_sectors, drill_chain_of, drill_children, format_drill_key,
    parse_stage_key,
)
from planetgen.db.corridor import (
    BypassBudgetExceeded, UnknownEdges, generated_cells, positions_near_segment, unknown_cells,
)
from planetgen.db import stones as nav_stones
from planetgen.galaxy.nav_graph import build_route_graph, shortest_path
from planetgen.galaxy.navigation import (
    FRAME_GALACTIC, FRAME_SECTOR, FRAME_SYSTEM, course_between, fold_travel_times, route_fold_times, route_warp_times,
    warp_travel_times,
)
from planetgen.physics.constants import SPECTRAL_CLASS_COLORS
from planetgen.tuning import (
    DEFAULT_SECTOR_EDGE_LY, HABITABLE_PLANET_CLASSES, NAV_ADJACENCY_K, NAV_CORRIDOR_FRACTION, NAV_CORRIDOR_MAX_LY,
    NAV_BYPASS_MAX_EDGES, NAV_CORRIDOR_MAX_SYSTEMS, NAV_CORRIDOR_MIN_LY, NAV_CORRIDOR_START_MAX_LY,
    NAV_ISLAND_LINKS, PLANET_CLASSES,
)
from planetgen.physics.units import ly_to_pc, milliparsecs_to_ly, mpc_to_pc, pc_to_ly
from planetgen.util import detmath
from planetgen.physics import unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")
_msun_to_kg = units.converter("Msun", "kg")

_au_to_km = units.converter("au", "km")
_m_to_ly = units.converter("m", "ly")
_km_to_m = units.converter("km", "m")
_s_to_yr = units.converter("s", "yr")
_pc_to_m = units.converter("pc", "m")
_km_to_au = units.converter("km", "au")
_m_to_au = units.converter("m", "au")


def open_readonly(config=None, statement_timeout_s=None):
    """
    Opens a connection for this read-only tool -- see the module
    docstring for why "read-only" is enforced by the configured account's
    grants rather than anything this function does itself.

    Args:
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.
        statement_timeout_s (float, optional): Stop any statement that
            runs longer (PERF.17) -- the web API passes `config.json`'s
            `mysql.statement_timeout_seconds`. `None`: no limit.

    Returns:
        planetgen.db.store.Connection: An open connection.

    Raises:
        SystemExit: If the database can't be reached.
    """
    try:
        return get_connection(config, ensure_schema=False, statement_timeout_s=statement_timeout_s)
    except pymysql.MySQLError as exc:
        raise SystemExit(f"Error: could not open the database ({exc}).")


SECTOR_SORTS = {
    "name": "sec.name", "systems": "c.system_count", "density": "c.system_density",
    "position": "c.quadrant", "distance": "sec.galactic_radius_pc",
}
"""dict: The sort keys `list_sectors` accepts (the Sectors table's column
keys) -> the column they order by (`c` is `sector_system_counts`). A sector
with no value (unplaced) sorts last either way."""

SECTOR_QUADRANTS = ("I", "II", "III", "IV")
"""tuple[str]: The sector Quadrant labels `list_sectors` filters by, plus
`"unplaced"` for a sector with no galaxy position."""

_SECTOR_QUADRANT_SQL = (
    "(CASE WHEN sec.center_x_pc IS NULL THEN 'unplaced' "
    "WHEN sec.center_y_pc >= 0 THEN (CASE WHEN sec.center_x_pc >= 0 THEN 'I' ELSE 'II' END) "
    "ELSE (CASE WHEN sec.center_x_pc < 0 THEN 'III' ELSE 'IV' END) END)"
)
"""str: SQL for a sector's Quadrant label, the same bands as
`planetgen.galaxy.spherical.sector_quadrant`."""


def _sector_quadrant_filter(quadrants):
    """`(where_sql, params)` keeping sectors in any of `quadrants`."""
    if not quadrants:
        return "", []
    return f" WHERE {_SECTOR_QUADRANT_SQL} IN ({', '.join('?' for _ in quadrants)})", list(quadrants)


def refresh_sector_system_counts(conn):
    """
    Makes `sector_system_counts` match the sectors and systems now (PERF.74, PERF.208): one row for each sector with its
    system count, its density (systems per cubic milliparsec, `NULL` for a sector with no galaxy position) and its
    Quadrant (`NULL` likewise), from one `GROUP BY` on the sector index, in one transaction so a reader sees the old
    rows or the new ones. `list_sectors` asks for it through `countcache`, which runs it in the background and not
    more often than every 30 seconds.

    Returns:
        int: The sectors that have a row.
    """
    conn.execute("DELETE FROM sector_system_counts")
    cur = conn.execute(
        "INSERT INTO sector_system_counts (sector_id, system_count, system_density, quadrant) "
        "SELECT sec.id, COALESCE(n.system_count, 0), "
        "CASE WHEN sec.center_x_pc IS NULL THEN NULL ELSE COALESCE(n.system_count, 0) / POW(sec.edge_mpc, 3) END, "
        f"NULLIF({_SECTOR_QUADRANT_SQL}, 'unplaced') "
        "FROM sectors sec LEFT JOIN (SELECT sector_id, COUNT(*) AS system_count FROM star_systems "
        "WHERE sector_id IS NOT NULL GROUP BY sector_id) n ON n.sector_id = sec.id")
    conn.commit()
    return cur.rowcount


def list_sectors(conn, limit=None, offset=None, sort=None, descending=False, quadrants=(), after=None):
    """
    Returns every sector, with its edge length (converted to light-years)
    and how many systems it contains (as of the last refresh of
    `sector_system_counts`, `refresh_sector_system_counts`), nearest the galactic core first
    (`galactic_radius_pc`); sectors never placed in a galaxy have no
    distance and come last, by name.

    The default order and the `distance` and `name` sorts are read by key (PERF.137): the placed sectors by an index
    walk of `galactic_radius_pc` and the unplaced ones, after them, by the name index, so the first page and the last
    cost the same on a galaxy of a million sectors (a sort of the whole table took 1 to 8 s). Sectors at the same
    distance follow the sort's direction by id. The `systems`, `density` and `position` sorts read the same way from
    `sector_system_counts` (PERF.208), each walk an index of it ending in the sector id; a sector made since the last
    refresh has no row there yet and follows the others, by id.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        limit (int, optional): Caps the number of rows returned. `None`
            (the default -- what every existing caller of this function
            still gets) returns every sector.
        offset (int, optional): Skips this many rows first. Ignored
            unless `limit` is also given; meaningless on its own.
        sort (str, optional): A key of `SECTOR_SORTS`; `None` is the
            default order below. Ties fall back to the default order, so
            pages never overlap.
        descending (bool): Reverse `sort`.
        quadrants (iterable[str]): Keep only sectors in these Quadrants
            (`SECTOR_QUADRANTS`, or `"unplaced"`).
        after (int, optional): The row id of the sector the previous page ended on. With a `limit` the page starts
            after that sector in every sort and `offset` is ignored.

    Returns:
        list[dict]: One row per sector, with `id`, `name`,
                           `edge_ly`, `system_count`,
                           `galactic_radius_pc`/`galactic_radius_ly`
                           (`None` if unplaced).
    """
    if sort is not None and sort not in SECTOR_SORTS:
        raise ValueError(f"unknown sector sort {sort!r}")
    _counted(conn, ["sector_system_counts"], refresh_sector_system_counts, lambda c: 0)
    if sort in (None, "name", "distance"):
        rows = _sectors_by_key(conn, limit, offset, sort, descending, quadrants, after)
    else:
        rows = _sectors_by_counts(conn, limit, offset, sort, descending, quadrants, after)
    return [
        {
            "id": r["id"], "name": r["name"], "edge_mpc": r["edge_mpc"], "edge_ly": milliparsecs_to_ly(r["edge_mpc"]),
            "system_count": r["system_count"],
            "center_x_pc": r["center_x_pc"], "center_y_pc": r["center_y_pc"], "center_z_pc": r["center_z_pc"],
            "galactic_radius_pc": r["galactic_radius_pc"],
            "galactic_radius_ly": (
                pc_to_ly(r["galactic_radius_pc"]) if r["galactic_radius_pc"] is not None else None
            ),
            "shell_index": r["shell_index"], "band_index": r["band_index"],
            "wedge_index": r["wedge_index"],
            "placed": r["center_x_pc"] is not None,
        }
        for r in rows
    ]


def _read_sectors(conn, where, params, order, limit, offset):
    """The rows of `list_sectors` for `where` (empty or ` WHERE ...`, with its `params`) in `order`."""
    query = f"""
        SELECT sec.id, sec.name, sec.edge_mpc, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.galactic_radius_pc, sec.shell_index, sec.band_index, sec.wedge_index,
               COALESCE(c.system_count, 0) AS system_count, {_SECTOR_QUADRANT_SQL} AS quadrant
        FROM sectors sec
        LEFT JOIN sector_system_counts c ON c.sector_id = sec.id{where}
        ORDER BY {order}
        """
    params = list(params)
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params.extend([limit, offset])
    return conn.execute(query, params).fetchall()


def _sectors_by_key(conn, limit, offset, sort, descending, quadrants, after):
    """
    The rows of `list_sectors` in the default order or the `name` or `distance` sort, each walk in the order an index
    serves (PERF.137). The `name` sort is one walk of `idx_sectors_name`. The others are two: the placed sectors by
    `idx_sectors_galactic_radius_pc` (backward for a descending `distance`), then the unplaced ones by name. `after`
    starts the page after that sector instead of skipping `offset` rows.
    """
    wanted = set(quadrants)
    backward = descending and sort in ("name", "distance")
    walks = [("all", None, "sec.name", descending)] if sort == "name" else []
    if sort != "name":
        if not wanted or wanted - {"unplaced"}:
            walks.append(("placed", "sec.galactic_radius_pc IS NOT NULL", "sec.galactic_radius_pc", backward))
        if not wanted or "unplaced" in wanted:
            walks.append(("unplaced", "sec.galactic_radius_pc IS NULL", "sec.name", False))
    where_quadrant, quadrant_params = _sector_quadrant_filter(quadrants)
    start = None
    if after is not None and limit is not None:
        start = conn.execute("SELECT name, galactic_radius_pc FROM sectors WHERE id = ?", (after,)).fetchone()
    skip = 0 if start is not None else offset or 0
    rows = []
    for phase, phase_condition, column, walk_backward in walks:
        conditions = [phase_condition] if phase_condition else []
        if where_quadrant:
            conditions.append(where_quadrant.removeprefix(" WHERE "))
        params = list(quadrant_params)
        placed_start = start is not None and start["galactic_radius_pc"] is not None
        if start is not None and not (phase == "unplaced" and placed_start):
            if phase == "placed" and not placed_start:
                continue        # the page began among the unplaced sectors: every placed one is behind it
            key = start["galactic_radius_pc"] if phase == "placed" else start["name"]
            condition, keyset_params = _keyset_condition([column], [key], walk_backward, "sec.id", after)
            conditions.append(condition)
            params += keyset_params
        where = " WHERE " + " AND ".join(conditions) if conditions else ""
        direction = "DESC" if walk_backward else "ASC"
        want = None if limit is None else limit - len(rows)
        page = _read_sectors(conn, where, params, f"{column} {direction}, sec.id {direction}", want, skip)
        rows.extend(page)
        if want is not None and len(page) >= want:
            break
        if skip and not page:       # the offset lies past this walk: what is left of it counts in the next
            counted = conn.execute(f"SELECT COUNT(*) AS n FROM sectors sec{where}", params).fetchone()["n"]
            skip = max(0, skip - counted)
        else:
            skip = 0
    return rows


_SECTOR_COLUMNS = ("sec.id, sec.name, sec.edge_mpc, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc, "
                   "sec.galactic_radius_pc, sec.shell_index, sec.band_index, sec.wedge_index")


_COUNTED_SECTORS = ("sector_system_counts c STRAIGHT_JOIN sectors sec ON sec.id = c.sector_id")


def _sectors_by_counts(conn, limit, offset, sort, descending, quadrants, after):
    """
    The rows of `list_sectors` in the `systems`, `density` and `position` sorts (PERF.208), each walk in the order an
    index of `sector_system_counts` serves, the sector id breaking ties in the sort's direction. `density` and
    `position` are two walks, the sectors with a value and then the ones with none (no galaxy position), which come
    last either way. A last walk reads the sectors made since the last refresh, which have no row yet, by id. `after`
    starts the page after that sector instead of skipping `offset` rows.
    """
    column = SECTOR_SORTS[sort].removeprefix("c.")
    direction = "DESC" if descending else "ASC"
    placed_wanted = [quadrant for quadrant in quadrants if quadrant != "unplaced"]
    placed_condition = ([f"c.quadrant IN ({', '.join('?' for _ in placed_wanted)})"], placed_wanted) if placed_wanted else ([], [])
    newest = conn.execute("SELECT COALESCE(MAX(sector_id), 0) AS n FROM sector_system_counts").fetchone()["n"]
    start = None                # the row the page starts after: its stored values, or {} for a sector with no row yet
    if after is not None and limit is not None:
        start = conn.execute("SELECT system_count, system_density, quadrant FROM sector_system_counts WHERE sector_id = ?",
                             (after,)).fetchone()
        if start is None and after > newest and conn.execute("SELECT 1 FROM sectors WHERE id = ?", (after,)).fetchone():
            start = {}
    skip = 0 if start is not None else offset or 0
    in_tail = start == {}
    started_unplaced = bool(start) and column != "system_count" and start[column] is None

    walks = []      # (name, conditions, params, the column it orders by)
    if column == "system_count":
        parts = placed_condition[0] + (["c.quadrant IS NULL"] if "unplaced" in quadrants else [])
        walks.append(("all", [f"({' OR '.join(parts)})"] if parts else [], placed_condition[1] if parts else [], column))
    else:
        if placed_wanted or not quadrants:
            walks.append(("placed", [f"c.{column} IS NOT NULL"] + placed_condition[0], placed_condition[1], column))
        if not quadrants or "unplaced" in quadrants:
            walks.append(("unplaced", [f"c.{column} IS NULL"], [], None))
    rows = []
    for phase, conditions, params, key_column in walks:
        if in_tail or (started_unplaced and phase != "unplaced"):
            continue
        conditions, params = list(conditions), list(params)
        if start and not (phase == "unplaced" and not started_unplaced):
            if key_column is None:
                conditions.append("c.sector_id < ?" if descending else "c.sector_id > ?")
                params.append(after)
            else:
                condition, keyset_params = _keyset_condition(
                    [f"c.{key_column}"], [start[key_column]], descending, "c.sector_id", after)
                conditions.append(condition)
                params += keyset_params
        where = " WHERE " + " AND ".join(conditions) if conditions else ""
        order = f"c.{key_column} {direction}, c.sector_id {direction}" if key_column else f"c.sector_id {direction}"
        want = None if limit is None else limit - len(rows)
        query = (f"SELECT {_SECTOR_COLUMNS}, c.system_count AS system_count FROM {_COUNTED_SECTORS}{where} ORDER BY {order}")
        page_params = list(params)
        if want is not None:
            query += " LIMIT ? OFFSET ?"
            page_params += [want, skip]
        page = conn.execute(query, page_params).fetchall()
        rows.extend(page)
        if want is not None and len(page) >= want:
            return rows
        if skip and not page:       # the offset lies past this walk: what is left of it counts in the next
            counted = conn.execute(f"SELECT COUNT(*) AS n FROM sector_system_counts c{where}", params).fetchone()["n"]
            skip = max(0, skip - counted)
        else:
            skip = 0
    conditions, params = ["sec.id > ?"], [newest]
    if quadrants:
        conditions.append(_sector_quadrant_filter(quadrants)[0].removeprefix(" WHERE "))
        params += list(quadrants)
    if in_tail:
        conditions.append("sec.id < ?" if descending else "sec.id > ?")
        params.append(after)
    query = (f"SELECT {_SECTOR_COLUMNS}, 0 AS system_count FROM sectors sec WHERE {' AND '.join(conditions)} "
             f"ORDER BY sec.id {direction}")
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params += [limit - len(rows), skip]
    rows.extend(conn.execute(query, params).fetchall())
    return rows


COUNT_FALLBACK_CAP = 10001
"""int: The most a request counts itself when the stored count isn't there yet (PERF.64, PERF.77): `count_*`
counts to this many and stops (30 ms where the full count took 6 s), answering `CappedCount(10000)` when it got
there, or with the table's estimated size when nothing filters it."""


class CappedCount(int):
    """A count that stopped at `COUNT_FALLBACK_CAP`: there are at least this many (shown "10,000+"). The API
    sends it as `total` with `total_capped: true`; the exact figure replaces it once the stored count is made."""

    capped = True


def _capped(n):
    return CappedCount(COUNT_FALLBACK_CAP - 1) if n >= COUNT_FALLBACK_CAP else n


def is_capped(count):
    """Whether `count` (from a `count_*` function) is a `CappedCount`."""
    return bool(getattr(count, "capped", False))


def _estimated_rows(conn, table):
    """The table's row count as the storage engine estimates it (instant, within some tens of percent)."""
    row = conn.execute("SELECT table_rows AS n FROM information_schema.tables "
                       "WHERE table_schema = DATABASE() AND table_name = ?", (table,)).fetchone()
    return int(row["n"] or 0) if row else 0


def _counted(conn, key, compute, fallback, stamp=None):
    """
    `countcache.cached` under the galaxy's content stamp (PERF.64), or the `stamp` given. The key carries the sector
    deletion epoch (`store.note_sector_deleted`), so a count made before a reset or a deletion is never served after it
    (UX.161): a stale answer served while the next is made is fine while the galaxy grows, but the page cache keeps
    what it was given, and after a wipe nothing moves the stamp again to replace it.
    """
    return countcache.cached(conn, [sector_deletions(conn), key], compute, fallback, stamp or galaxy_content_stamp)


def count_sectors(conn, quadrants=()):
    """
    Returns the total number of sectors (that pass the same `quadrants`
    filter as `list_sectors`), ignoring any pagination --
    the denominator `list_sectors(conn, limit=...)` callers (the API's
    `/api/sectors`) need to report how many pages exist. Kept by
    `countcache` (PERF.64): the last count made when the galaxy has changed
    since, an estimate before the first.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.

    Returns:
        int: Total sector count.
    """
    where, params = _sector_quadrant_filter(quadrants)

    def exact(c):
        return c.execute(f"SELECT COUNT(*) AS n FROM sectors sec{where}", params).fetchone()["n"]

    def estimate(c):
        if not where:
            return _estimated_rows(c, "sectors")
        return _capped(c.execute(f"SELECT COUNT(*) AS n FROM (SELECT 1 FROM sectors sec{where} LIMIT ?) c",
                                 list(params) + [COUNT_FALLBACK_CAP]).fetchone()["n"])

    return _counted(conn, ["count_sectors", list(quadrants)], exact, estimate)


def sectors_facets(conn, quadrants=()):
    """
    The option counts for the Sectors table's Quadrant menu: how many sectors
    each Quadrant (and `"unplaced"`) holds. The menu ignores its own filter,
    so choosing one Quadrant still shows the others' counts.

    Returns:
        dict: `{"quadrant": [{"value", "count"}]}` in Quadrant order, only
            those with sectors.
    """
    def exact(c):
        counts = {row["quadrant"]: row["n"] for row in c.execute(
            f"SELECT {_SECTOR_QUADRANT_SQL} AS quadrant, COUNT(*) AS n FROM sectors sec GROUP BY quadrant"
        ).fetchall()}
        return {"quadrant": [{"value": value, "count": counts[value]}
                             for value in (*SECTOR_QUADRANTS, "unplaced") if counts.get(value)]}

    return _counted(conn, ["sectors_facets"], exact, lambda c: {"quadrant": []})


class _NoSector:
    """Sentinel type for `NO_SECTOR` (below) -- a distinct value from
    `None` (`sector_id` unfiltered) meaning "explicitly filter to
    `sector_id IS NULL`" (standalone systems, e.g. `html/browse.py`'s own
    table of systems generated with no sector)."""

    def __repr__(self):
        return "NO_SECTOR"


NO_SECTOR = _NoSector()
"""_NoSector: Pass as `list_systems`/`count_systems`'s `sector_id` to
match only standalone systems (`star_systems.sector_id IS NULL`) --
distinct from the default `None`, which means "don't filter by sector at
all." The API's `/api/systems?sector_id=none` maps onto this."""


STAR_COUNT_SQL = "(SELECT COUNT(*) FROM stars st WHERE st.star_system_id = ss.id)"
"""str: How many stars a system holds, as SQL on a `star_systems ss` row: the rows of `stars` it owns (UX.239)."""

SYSTEM_SORTS = {"name": "ss.name", "sector": "sec.name", "octant": "ss.quadrant", "stars": STAR_COUNT_SQL}
"""dict: The sort keys `list_systems` accepts (the Systems tables' column
keys) -> the SQL they order by. Ties fall back to name, then id. `stars` orders by the number of stars (UX.239)."""


def _star_count_conditions(stars):
    """
    The `WHERE` conditions and their parameters for the stars filter `stars` (`(least, most)` of `parse_star_count`, `most` `None` for "or more") on a
    `star_systems ss` row: `([], [])` for no filter. One definition for the search, the Systems list and its count (UX.239).
    """
    if stars is None:
        return [], []
    least, most = stars
    conditions, params = [f"{STAR_COUNT_SQL} >= ?"], [least]
    if most is not None:
        conditions.append(f"{STAR_COUNT_SQL} <= ?")
        params.append(most)
    return conditions, params


def _after_prefix(text):
    """The smallest text that sorts after every text starting with `text` (its last character stepped up), or `None`
    when there is none: the bound a descending "starting at M" list ends below, so the names that start with M come
    first."""
    characters = list(text)
    while characters:
        code = ord(characters[-1]) + 1
        if code == 0xD800:
            code = 0xE000
        if code <= 0x10FFFF:
            characters[-1] = chr(code)
            return "".join(characters)
        characters.pop()
    return None


def _start_condition(sort, descending, start):
    """
    The condition that keeps a systems list to the rows from `start` on in its `sort` and direction (the jump-by-value
    box, PERF.75): the names from `start` (`sort="name"`), or the systems of the sectors named from `start` on, with
    the systems that have no sector, which sort last, kept (`sort="sector"`). A descending list starts at the names
    that begin with `start` and runs down from them.

    Returns:
        tuple or None: `(condition_sql, params)`, or `None` for an empty `start`, a sort that has no start (octant) or a descending start that nothing sorts above.
    """
    if not start or sort not in ("name", "sector"):
        return None
    if descending:
        bound = _after_prefix(start)
        if bound is None:
            return None
        operator, value = "<", bound
    else:
        operator, value = ">=", start
    if sort == "name":
        return f"ss.name {operator} ?", [value]
    return (f"(ss.sector_id IS NULL OR ss.sector_id IN (SELECT id FROM sectors WHERE name {operator} ?))", [value])


def _systems_filter_clause(star_type_prefix, sector_id, binary=None, octants=(), in_sector=None, start=None,
                           sort="name", descending=False, stars=None):
    """
    Builds the shared `JOIN`/`WHERE`/params fragment `list_systems` and
    `count_systems` both need -- factored out so the count query can't
    silently drift out of sync with what the listing query actually
    matches. `start` keeps the rows from that value on in `sort` (`_start_condition`).

    Returns:
        tuple: `(join_sql, where_sql, params)`, each usable standalone
              (empty strings/list when no filter applies).
    """
    join_sql = ""
    conditions = []
    params = []

    if star_type_prefix is not None:
        join_sql = " JOIN stars s ON s.star_system_id = ss.id"
        conditions.append("s.star_type LIKE ? ESCAPE '\\\\'")
        params.append(f"{escape_like(star_type_prefix)}%")

    if sector_id is NO_SECTOR:
        conditions.append("ss.sector_id IS NULL")
    elif sector_id is not None:
        conditions.append("ss.sector_id = ?")
        params.append(sector_id)

    if binary is not None:
        conditions.append("ss.is_binary = ?")
        params.append(1 if binary else 0)
    if octants:
        conditions.append(f"ss.quadrant IN ({', '.join('?' for _ in octants)})")
        params.extend(octants)
    if in_sector is not None:
        conditions.append("ss.sector_id IS NOT NULL" if in_sector else "ss.sector_id IS NULL")
    star_conditions, star_params = _star_count_conditions(stars)
    conditions.extend(star_conditions)
    params.extend(star_params)
    begun = _start_condition(sort, descending, start)
    if begun:
        conditions.append(begun[0])
        params.extend(begun[1])

    where_sql = (" WHERE " + " AND ".join(conditions)) if conditions else ""
    return join_sql, where_sql, params


def list_systems(conn, star_type_prefix=None, sector_id=None, limit=None, offset=None, sort="name",
                 descending=False, binary=None, octants=(), in_sector=None, after=None, start=None, stars=None):
    """
    Returns systems, optionally filtered by star type and/or sector.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        star_type_prefix (str, optional): Matches any system with at least
            one star (single, primary, or secondary) whose `star_type`
            starts with this text, e.g. `"G"` for every G-type system,
            `"G2V"` for an exact spectral/subclass/luminosity match.
            Case-sensitive, matching the stored spectral class letters.
        sector_id (int or NO_SECTOR, optional): Restricts to one sector's
            systems, or (via the `NO_SECTOR` sentinel) to systems with no
            sector at all. `None` (the default) doesn't filter by sector.
        limit (int, optional): Caps the number of rows returned. `None`
            (the default -- what every existing caller of this function
            still gets) returns every matching system.
        offset (int, optional): Skips this many rows first. Ignored
            unless `limit` is also given; meaningless on its own.
        sort (str): A key of `SYSTEM_SORTS`; a system with no sector sorts
            last by `"sector"` either way.
        descending (bool): Reverse `sort`.
        binary (bool, optional): Keep only binary (or only single-star) systems.
        octants (iterable[str]): Keep only systems in these sector octants.
        in_sector (bool, optional): Keep only systems in a sector (True) or
            standalone ones (False).
        stars (tuple, optional): Keep only systems with that many stars, `(least, most)` of `parse_star_count` (UX.239).
        after (int, optional): A system's row ID: with `sort="name"`, `"sector"` or `"octant"` the page starts
            after that system (a keyset page, PERF.75, PERF.206: the rows before it are not read and `offset` is
            ignored). Other sorts, a system that is gone and a list filtered by star type or sector use `offset`.
        start (str, optional): The jump-by-value box (PERF.75): the list begins at the first name (`sort="name"`) or
            the first sector (`sort="sector"`) from this text on, in the sort's direction, and `offset` counts from
            there. Ignored for the other sorts.

    Returns:
        list[dict]: One row per matching system, with `id`, `name`,
                           `sector_id`, `sector_name`, `quadrant` (its
                           octant in the sector), `is_binary`, `star_count`,
                           `object_count`, `has_habitable_world` (`_system_counts`), and `star_summary`
                           (the single star's `star_type`, or a binary's
                           `binary_type` -- what `html/browse.py`/
                           `html/search.py` show as a system's "Star type"
                           column).
    """
    if sort not in SYSTEM_SORTS:
        raise ValueError(f"unknown system sort {sort!r}")
    # The sector walk begins at its first sector itself (`_systems_by_group`), so its groups carry no start condition.
    walked = (limit is not None and star_type_prefix is None and sector_id is None and in_sector is None
              and sort in _GROUP_WALK_SORTS)
    join_sql, where_sql, params = _systems_filter_clause(
        star_type_prefix, sector_id, binary, octants, in_sector,
        start=None if walked and sort == "sector" else start, sort=sort, descending=descending, stars=stars)
    column = SYSTEM_SORTS[sort]
    direction = "DESC" if descending else "ASC"
    if after is not None and limit is not None and sort == "name":
        previous = conn.execute("SELECT name FROM star_systems WHERE id = ?", (after,)).fetchone()
        if previous is not None:
            condition, keyset_params = _keyset_condition(["ss.name"], [previous["name"]], descending, "ss.id", after)
            where_sql = _and_where(where_sql, condition)
            params = params + keyset_params
            offset = 0
    if walked:
        rows = _systems_by_group(conn, sort, descending, where_sql, params, limit, offset or 0, start, after)
        return _system_rows(conn, rows)
    nulls_last = f"{column} IS NULL, " if sort in ("sector", "octant") else ""
    # The name sort's ties follow its direction, so one walk of the name index serves either way round (PERF.75).
    tie_break = f"ss.id {direction}" if sort == "name" else "ss.name, ss.id"
    # DISTINCT only matters when the star join can repeat a system; without it the sort can stop at the page (PERF.64).
    distinct = "DISTINCT " if join_sql else ""
    query = f"""
        SELECT {distinct}ss.id, ss.name, ss.sector_id, ss.quadrant, ss.is_binary, ss.binary_configuration,
               ss.binary_type, sec.name AS sector_sort
        FROM star_systems ss LEFT JOIN sectors sec ON sec.id = ss.sector_id{join_sql}{where_sql}
        ORDER BY {nulls_last}{column} {direction}, {tie_break}
        """

    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params = params + [limit, offset or 0]

    return _system_rows(conn, conn.execute(query, params).fetchall())


def _system_rows(conn, rows):
    """`list_systems`' dicts for `rows` (`star_systems` rows)."""
    rows = _with_star_types(conn, rows, with_sector_name=True)
    counts = _system_counts(conn, [r["id"] for r in rows])
    return [
        {
            "id": r["id"], "name": r["name"], "sector_id": r["sector_id"], "sector_name": r["sector_name"],
            "quadrant": r["quadrant"], "is_binary": r["is_binary"], "star_summary": _star_summary(r),
            **counts[r["id"]],
        }
        for r in rows
    ]


def _system_counts(conn, ids):
    """
    The figures the Systems tables show beside a system, read for the whole page (PERF.15) from the rows that exist and
    kept nowhere else: `star_count` (every row of `stars`: a star, a brown dwarf, a white dwarf, a neutron star, a
    black hole, in a tree or not), `object_count` (those plus the system's planets, moons and asteroid belts) and
    `has_habitable_world` (a planet or moon scored Yellow or better, `habitability.YELLOW_FLOOR`, the Galaxy Map's
    `habitable` reason).

    Returns:
        dict: system ID -> `{"star_count", "object_count", "has_habitable_world"}`, every ID in `ids` present.
    """
    counts = {i: {"star_count": 0, "object_count": 0, "has_habitable_world": False} for i in ids}
    if not ids:
        return counts
    marks = ", ".join("?" * len(ids))
    tables = ("stars", "planets", "moons", "asteroid_belts")
    # Two queries for the page, not one per table (the page's read count is held by test_db_web_reads).
    for r in conn.execute(" UNION ALL ".join(
            f"SELECT star_system_id AS id, '{table}' AS kind, COUNT(*) AS n FROM {table} WHERE star_system_id IN ({marks}) "
            f"GROUP BY star_system_id" for table in tables), ids * len(tables)).fetchall():
        counts[r["id"]]["object_count"] += r["n"]
        if r["kind"] == "stars":
            counts[r["id"]]["star_count"] = r["n"]
    for r in conn.execute(" UNION ".join(
            f"SELECT star_system_id AS id FROM {table} WHERE star_system_id IN ({marks}) AND phi4 >= ?"
            for table in ("planets", "moons")), [*ids, habitability.YELLOW_FLOOR] * 2).fetchall():
        counts[r["id"]]["has_habitable_world"] = True
    return counts


_GROUP_WALK_SORTS = ("sector", "octant")
"""tuple[str]: The `list_systems` sorts served a group at a time by `_systems_by_group`."""

_GROUP_WALK_COLUMNS = ("ss.id, ss.name, ss.sector_id, ss.quadrant, ss.is_binary, ss.binary_configuration, "
                       "ss.binary_type")


def _and_where(where_sql, condition):
    """`where_sql` (empty or ` WHERE a AND b`) with `condition` added."""
    return f"{where_sql} AND {condition}" if where_sql else f" WHERE {condition}"


def _lexicographic(columns, values, descending, inclusive=False):
    """
    The rows after one row in `ORDER BY <columns>`, all sorted the same way: those whose first column is past its value, or equal to
    it and the second past its value, and so on; with `inclusive` the row itself and the ones that tie with it all through count too.

    Written with `OR` and `AND`, never as a row constructor: `(a, b, id) > (?, ?, ?)` took 5.0 s on MariaDB 10.11 with no index range
    (and a full walk of the index on MySQL 8), the same condition spelled out 2 ms (docs/design/operations/slow-reads-and-timeouts.md).

    Returns:
        tuple: `(condition_sql, params)`; the condition is parenthesised.
    """
    operator = "<" if descending else ">"
    terms, params, equal, equal_params = [], [], [], []
    last = len(columns) - 1
    for number, (column, value) in enumerate(zip(columns, values)):
        compare = operator + "=" if inclusive and number == last else operator
        terms.append("(" + " AND ".join([*equal, f"{column} {compare} ?"]) + ")")
        params.extend([*equal_params, value])
        equal.append(f"{column} = ?")
        equal_params.append(value)
    return "(" + " OR ".join(terms) + ")", params


def _keyset_condition(columns, values, descending, id_column, after_id):
    """
    The rows after one row in `ORDER BY <columns> <direction>, <id_column> <direction>` (a keyset page, PERF.75):
    the rows whose sort key is past `values` (the row's own, in `columns` order) or equal to it with a past id.

    The id sorts the same way as the columns, so one backward or forward walk of an index on the columns (which
    holds the id as its last part) serves the page and the sort needs no filesort. Written as `_lexicographic` is.

    Args:
        columns (list[str]): The sort columns, in order, all sorted the same way.
        values (list): The row's value of each column.
        descending (bool): Whether the columns and the id sort downward.
        id_column (str): The unique tie-break column.
        after_id (int): The row's id.

    Returns:
        tuple: `(condition_sql, params)`; the condition is parenthesised.
    """
    return _lexicographic([*columns, id_column], [*values, after_id], descending)


def _simple_groups(sort, descending):
    """The groups of an octant sort in order: `(condition_sql, params)`, the systems with none last."""
    labels = (("VIII", "VII", "VI", "V", "IV", "III", "II", "I") if descending else
              ("I", "II", "III", "IV", "V", "VI", "VII", "VIII"))
    return [("ss.quadrant = ?", [label]) for label in labels] + [("ss.quadrant IS NULL", [])]


def _read_group(conn, where_sql, params, condition, group_params, need, offset):
    return conn.execute(
        f"SELECT {_GROUP_WALK_COLUMNS} FROM star_systems ss{_and_where(where_sql, condition)} "
        f"ORDER BY ss.name, ss.id LIMIT ? OFFSET ?", list(params) + group_params + [need, offset]).fetchall()


def _systems_by_group(conn, sort, descending, where_sql, params, limit, offset, start=None, after=None):
    """
    One page of `list_systems` sorted by sector or octant, without sorting the table (PERF.70): the
    sort is `<group>, name, id` with the systems that have no sector or octant last, so the page is found by
    counting group after group (an index range each: `idx_star_systems_sector_id`, `idx_star_systems_quadrant`) until `offset` is used up, then reading the rows from the group it lands in
    and the groups after it. The sectors are taken in chunks, counted by one query per chunk. A `start` (the
    jump-by-value box) begins the sector sort at the first sector from that name on.

    With `after` (the id of the system the previous page ended on, PERF.206) the page is read from that system on:
    the groups before its own are not counted, and its own group continues after its name and id (the names
    within a group always run upward), so any page costs one group's index range and the groups after it, however
    deep. A system that is gone, or whose group the filters leave out, falls back to `offset`.

    Returns:
        list: `star_systems` rows, `limit` at most.
    """
    rows = []
    resume = None       # (condition, params) that continues the group the previous page ended in
    if after is not None:
        found = conn.execute("SELECT name, sector_id, quadrant, is_binary FROM star_systems WHERE id = ?",
                             (after,)).fetchone()
        if found is not None:
            resume = _keyset_condition(["ss.name"], [found["name"]], False, "ss.id", after)

    def take(condition, group_params, continuing=False):
        nonlocal offset
        if continuing:
            condition, group_params = f"{condition} AND {resume[0]}", group_params + resume[1]
        rows.extend(_read_group(conn, where_sql, params, condition, group_params, limit - len(rows), offset))
        offset = 0
        return len(rows) >= limit

    if sort != "sector":
        groups = _simple_groups(sort, descending)
        begin = 0
        if resume is not None:
            own = (("ss.quadrant = ?", [found["quadrant"]]) if found["quadrant"] is not None else
                   ("ss.quadrant IS NULL", []))
            if own in groups:
                begin, offset = groups.index(own), 0
            else:
                resume = None
        for number, (condition, group_params) in enumerate(groups):
            if number < begin:
                continue
            if offset:
                n = conn.execute(f"SELECT COUNT(*) AS n FROM star_systems ss{_and_where(where_sql, condition)}",
                                 list(params) + group_params).fetchone()["n"]
                if n <= offset:
                    offset -= n
                    continue
            if take(condition, group_params, continuing=resume is not None and number == begin):
                break
        return rows

    direction = "DESC" if descending else "ASC"
    begun = None if resume is not None else _start_condition("sector", descending, start)
    first = _and_where("", f"name {'<' if descending else '>='} ?") if begun else ""
    first_params = [begun[1][0]] if begun else []
    if resume is not None:
        offset = 0
        if found["sector_id"] is None:
            take("ss.sector_id IS NULL", [], continuing=True)
            return rows
        name = conn.execute("SELECT name FROM sectors WHERE id = ?", (found["sector_id"],)).fetchone()["name"]
        beyond, through = ("<", "<=") if descending else (">", ">=")
        first = _and_where("", f"(name {beyond} ? OR (name = ? AND id {through} ?))")
        first_params = [name, name, found["sector_id"]]
    taken, size = 0, 20
    while True:
        ids = [r["id"] for r in conn.execute(
            f"SELECT id FROM sectors{first} ORDER BY name {direction}, id {direction} LIMIT ? OFFSET ?",
            first_params + [size, taken]).fetchall()]
        if not ids:
            break
        marks = ", ".join("?" for _ in ids)
        counted = {r["sector_id"]: r["n"] for r in conn.execute(
            f"SELECT ss.sector_id, COUNT(*) AS n FROM star_systems ss{_and_where(where_sql, f'ss.sector_id IN ({marks})')} "
            f"GROUP BY ss.sector_id", list(params) + ids).fetchall()}
        for sector in ids:
            n = counted.get(sector, 0)
            continuing = resume is not None and sector == found["sector_id"]
            if n <= offset and not continuing:
                offset -= n
                continue
            if take("ss.sector_id = ?", [sector], continuing=continuing):
                return rows
        taken += len(ids)
        size = min(size * 2, 500)
    take("ss.sector_id IS NULL", [])
    return rows


_STAR_TYPE_COLUMNS = {"single": "single_star_type", "primary": "primary_star_type",
                      "secondary": "secondary_star_type"}


def _with_star_types(conn, rows, with_sector_name=False):
    """
    `rows` (`star_systems` rows) as dicts carrying `single_star_type`/
    `primary_star_type`/`secondary_star_type` (and `sector_name`, if
    asked), read for the whole page in one query each (PERF.15) instead of
    three correlated subqueries per row.
    """
    rows = [dict(row) for row in rows]
    if not rows:
        return rows
    ids = [row["id"] for row in rows]
    marks = ", ".join("?" * len(ids))
    types = {}
    for star in conn.execute(
        f"SELECT star_system_id, role, star_type FROM stars WHERE star_system_id IN ({marks}) ORDER BY id",
        ids,
    ).fetchall():
        types.setdefault((star["star_system_id"], star["role"]), star["star_type"])
    sector_names = {}
    if with_sector_name:
        sector_ids = sorted({row["sector_id"] for row in rows if row["sector_id"] is not None})
        if sector_ids:
            sector_names = {r["id"]: r["name"] for r in conn.execute(
                f"SELECT id, name FROM sectors WHERE id IN ({', '.join('?' * len(sector_ids))})", sector_ids,
            ).fetchall()}
    for row in rows:
        for role, column in _STAR_TYPE_COLUMNS.items():
            row[column] = types.get((row["id"], role))
        if with_sector_name:
            row["sector_name"] = sector_names.get(row["sector_id"])
    return rows


def _star_summary(row):
    """
    Builds the "Star type" summary `list_systems`/`_search_result_systems`
    show, from a row carrying `is_binary`/`binary_configuration`/
    `binary_type`/`single_star_type`/`primary_star_type`/`secondary_star_type`.

    - Single star: that star's own `star_type`.
    - `'close'` (P-type) binary: the merged `BinaryStarProxy`'s `binary_type`
      string (e.g. `"Binary (G/K)"`).
    - `'wide'` (S-type) binary: `binary_type` is NULL (no merged effective
      star exists to summarize -- see `schema.sql`'s "v15" note), so this
      builds an equivalent summary directly from the two stars' own types
      instead, rather than showing a blank cell.

    Args:
        row (dict): A `star_systems` row (or equivalent), joined with the
            per-role star-type subqueries above.

    Returns:
        str or None: The summary string, or `None` for a pre-v15 wide-binary
            row this can't happen for (every `is_binary` row predating v15
            was necessarily `'close'`).
    """
    if not row["is_binary"]:
        return row["single_star_type"]
    if row["binary_configuration"] == "wide":
        return f"{row['primary_star_type']} / {row['secondary_star_type']} (wide binary)"
    return row["binary_type"]


def _start_key(sort, descending, start):
    """The part of a count's stored key that a jump-by-value start adds: none without a usable start."""
    return None if _start_condition(sort, descending, start) is None else [sort, descending, start]


def count_systems(conn, star_type_prefix=None, sector_id=None, binary=None, octants=(), in_sector=None, start=None,
                  sort="name", descending=False, stars=None):
    """
    Returns the total number of systems matching the same filters
    `list_systems` accepts, ignoring any pagination -- the denominator
    `list_systems(conn, limit=...)` callers (the API's `/api/systems`)
    need to report how many pages exist.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        star_type_prefix (str, optional): Same meaning as `list_systems`.
        sector_id (int or NO_SECTOR, optional): Same meaning as `list_systems`.
        start, sort, descending: The jump-by-value box's start (PERF.75), as in `list_systems`: only the systems from
            it on are counted.
        stars (tuple, optional): The stars filter of `list_systems`.

    Returns:
        int: Total matching system count.
    """
    join_sql, where_sql, params = _systems_filter_clause(star_type_prefix, sector_id, binary, octants, in_sector,
                                                         start, sort, descending, stars)
    # A system matches a star type once however many of its stars do; with no star join every row is its own.
    counted = "COUNT(DISTINCT ss.id)" if join_sql else "COUNT(*)"

    def exact(c):
        return c.execute(f"SELECT {counted} AS n FROM star_systems ss{join_sql}{where_sql}", params).fetchone()["n"]

    def estimate(c):
        if not join_sql and not where_sql:
            return _estimated_rows(c, "star_systems")
        return _capped(c.execute(f"SELECT COUNT(*) AS n FROM (SELECT DISTINCT ss.id FROM star_systems ss{join_sql}{where_sql} "
                                 f"LIMIT ?) c", list(params) + [COUNT_FALLBACK_CAP]).fetchone()["n"])

    return _counted(conn, ["count_systems", star_type_prefix, sector_id, binary, list(octants), in_sector,
                           _start_key(sort, descending, start), None if stars is None else list(stars)], exact, estimate)


def systems_facets(conn, sector_id=None, binary=None, octants=(), in_sector=None, stars=None):
    """
    The option counts for the Systems tables' filter menus: where a system is
    (`placement`: `"sector"` or `"standalone"`), whether it is `binary`
    (`"yes"`/`"no"`) and its sector `octant`. Each menu's counts apply every
    filter except its own, so the menus narrow one another. `stars` (UX.239) is a filter of its own with no menu here.

    Returns:
        dict: `{"placement"|"binary"|"octant": [{"value", "count"}]}`,
            options with no systems left out.
    """
    def exact(c):
        def counts(select, **filters):
            base = {"binary": binary, "octants": octants, "in_sector": in_sector, "stars": stars, **filters}
            _join, where, params = _systems_filter_clause(None, sector_id, **base)
            return {row["value"]: row["n"] for row in c.execute(
                f"SELECT {select} AS value, COUNT(*) AS n FROM star_systems ss{where} GROUP BY value",
                params).fetchall()}

        placement = counts("(CASE WHEN ss.sector_id IS NULL THEN 'standalone' ELSE 'sector' END)", in_sector=None)
        binaries = counts("(CASE WHEN ss.is_binary THEN 'yes' ELSE 'no' END)", binary=None)
        found = counts("ss.quadrant", octants=())
        return {
            "placement": [{"value": v, "count": placement[v]} for v in ("sector", "standalone") if placement.get(v)],
            "binary": [{"value": v, "count": binaries[v]} for v in ("yes", "no") if binaries.get(v)],
            "octant": [{"value": v, "count": found[v]} for v in sorted(x for x in found if x is not None)],
        }

    return _counted(conn, ["systems_facets", sector_id, binary, list(octants), in_sector, None if stars is None else list(stars)], exact,
                    lambda c: {"placement": [], "binary": [], "octant": []})


def _body_filter_clause(table_alias, planet_class, min_radius_km, max_radius_km, sector_id, system_id):
    """
    Builds the shared `WHERE`/params fragment `list_planets` and
    `list_moons` both need -- factored out the same way
    `_systems_filter_clause` is, so the two stay consistent with each
    other rather than each hand-rolling the same class/size/sector/system
    branching.

    Args:
        table_alias (str): The `planets`/`moons` table's own alias in the
            calling query (`"p"` or `"m"`) -- prefixes `planet_class`/
            `radius_km` in the clauses this returns.
        planet_class (str, optional): Exact `planet_class` match (e.g.
            `"M"`). `None` means no class filter.
        min_radius_km (float, optional): Lower bound, inclusive.
        max_radius_km (float, optional): Upper bound, inclusive.
        sector_id (int, optional): Restricts to bodies whose system is in
            this sector.
        system_id (int, optional): Restricts to bodies in this one system.

    Returns:
        tuple: `(where_sql, params)`.
    """
    clauses = []
    params = []
    if planet_class is not None:
        clauses.append(f"{table_alias}.planet_class = ?")
        params.append(planet_class)
    size_range = None
    if min_radius_km is not None or max_radius_km is not None:
        size_range = (min_radius_km, max_radius_km)
    _append_size_clause(clauses, params, f"{table_alias}.radius_km", size_range)
    if system_id is not None:
        clauses.append(f"{table_alias}.star_system_id = ?")
        params.append(system_id)
    if sector_id is not None:
        clauses.append("ss.sector_id = ?")
        params.append(sector_id)
    where_sql = (" WHERE " + " AND ".join(clauses)) if clauses else ""
    return where_sql, params


def _body_position(conn, kind, after):
    """The sort values of the planet or moon (`kind` `"planet"` or `"moon"`) with row id `after`, as `(columns, values,
    id_column, system_name, system_id)` for `_keyset_condition`, or `None` when there is no such body."""
    if kind == "planet":
        found = conn.execute(
            "SELECT ss.name AS name, ss.id AS system_id, p.orbital_index AS orbit FROM planets p "
            "JOIN star_systems ss ON ss.id = p.star_system_id WHERE p.id = ?", (after,)).fetchone()
        if found is None:
            return None
        return (["ss.name", "ss.id", "p.orbital_index"], [found["name"], found["system_id"], found["orbit"]], "p.id",
                found["name"], found["system_id"])
    found = conn.execute(
        "SELECT ss.name AS name, ss.id AS system_id, p.orbital_index AS planet_orbit, p.id AS planet_id, "
        "m.orbital_index AS orbit FROM moons m JOIN planets p ON p.id = m.planet_id "
        "JOIN star_systems ss ON ss.id = m.star_system_id WHERE m.id = ?", (after,)).fetchone()
    if found is None:
        return None
    return (["ss.name", "ss.id", "p.orbital_index", "p.id", "m.orbital_index"],
            [found["name"], found["system_id"], found["planet_orbit"], found["planet_id"], found["orbit"]], "m.id",
            found["name"], found["system_id"])


def _bodies_after(conn, kind, where_sql, params, after):
    """
    The keyset condition for `list_planets` and `list_moons` (PERF.207): the bodies after the one with row id `after` in
    their order (`kind` `"planet"` or `"moon"`), as `(where_sql, params)` with it added; unchanged when `after` is `None`
    or names nothing, so the page then starts at the top.

    A planet's order is system name, system, orbit and id; a moon's is system name, system, its planet's orbit and id,
    then its own orbit and id.
    """
    position = None if after is None else _body_position(conn, kind, after)
    if position is None:
        return where_sql, params
    columns, values, id_column, _name, _system_id = position
    condition, keyset_params = _keyset_condition(columns, values, False, id_column, after)
    return _and_where(where_sql, condition), params + keyset_params


BODY_CHUNK_SYSTEMS = 1000
"""int: The most systems `_bodies_by_system` reads bodies of at once (PERF.215)."""


def _bodies_by_system(conn, kind, where_sql, params, query_for, limit, after):
    """
    One page of the unfiltered planets or moons listing in the order of the systems' name index (PERF.215). The order
    spans two tables (the system's name, then the body's orbit), which no index serves in a join: the join sorted every
    body of the galaxy for a page (18 s for the first page with 6,000,000 planets). This walks the systems by name in
    chunks and reads the bodies of each chunk, sorting only those, until the page is full.

    Args:
        query_for (callable): `query_for(where_sql)` is the listing's `SELECT` with that `WHERE` and its `ORDER BY`.
        after (int, optional): The row id the page starts after, as `list_planets` takes it.
    """
    position = None if after is None else _body_position(conn, kind, after)
    first_where, first_params = _bodies_after(conn, kind, where_sql, params, after)
    rows = []
    previous = None             # the (name, id) of the last system read; None before the first chunk
    while len(rows) < limit:
        want = limit - len(rows)
        chunk = min(max(4 * want, 100), BODY_CHUNK_SYSTEMS)
        if previous is not None:
            condition, system_params = "(ss.name > ? OR (ss.name = ? AND ss.id > ?))", [previous[0], previous[0], previous[1]]
        elif position is not None:      # the system the page starts in is read again, from the body after `after`
            condition, system_params = "(ss.name > ? OR (ss.name = ? AND ss.id >= ?))", [position[3], position[3], position[4]]
        else:
            condition, system_params = "1", []
        systems = conn.execute(f"SELECT ss.id, ss.name FROM star_systems ss WHERE {condition} ORDER BY ss.name, ss.id LIMIT ?",
                               system_params + [chunk]).fetchall()
        if not systems:
            break
        chunk_where = _and_where(first_where if previous is None else where_sql,
                                 f"ss.id IN ({', '.join('?' for _ in systems)})")
        chunk_params = (first_params if previous is None else params) + [row["id"] for row in systems]
        rows.extend(conn.execute(query_for(chunk_where) + " LIMIT ?", chunk_params + [want]).fetchall())
        previous = (systems[-1]["name"], systems[-1]["id"])
        if len(systems) < chunk:
            break
    return rows


def list_planets(conn, planet_class=None, min_radius_km=None, max_radius_km=None,
                  sector_id=None, system_id=None, limit=None, offset=None, after=None):
    """
    Returns planets (never asteroid belts -- `body_type` is always `'t'`/
    `'g'` here, `schema.sql`'s `planets` table holds only those), optionally
    filtered by class, radius range, sector, and/or system -- the `planets`
    equivalent of `list_systems` above, closing the gap `docs/TODO.md`'s
    "Open items" > "Search" flagged: the web/API faceted search
    (`search()`/`GET /api/search`) already supports a planet size range and
    class tags, but this CLI's own subcommands never got an equivalent.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        planet_class (str, optional): Exact `planet_class` match (e.g. `"M"`).
        min_radius_km (float, optional): Only planets at least this large.
        max_radius_km (float, optional): Only planets at most this large.
        sector_id (int, optional): Only planets whose system is in this sector.
        system_id (int, optional): Only planets in this one system.
        limit (int, optional): Caps the number of rows returned. `None`
            (the default) returns every matching planet.
        offset (int, optional): Skips this many rows first. Ignored
            unless `limit` is also given, and when `after` is.
        after (int, optional): The row id of the planet the previous page ended on: the page starts after it in the
            order (system name, system, orbit, id) and `offset` is not used, so any page costs the same (PERF.207).

    Returns:
        list[dict]: One row per matching planet, with `id`, `name`,
            `planet_class`, `body_type`, `radius_km`, `star_system_id`,
            `system_name`.
    """
    where_sql, params = _body_filter_clause("p", planet_class, min_radius_km, max_radius_km, sector_id, system_id)

    def query_for(where):
        return f"""
            SELECT p.id, p.name, p.planet_class, p.body_type, p.radius_km, p.star_system_id, ss.name AS system_name
            FROM planets p
            JOIN star_systems ss ON ss.id = p.star_system_id{where}
            ORDER BY ss.name, ss.id, p.orbital_index, p.id
            """

    if limit is not None and not params and (after is not None or not offset):
        return [dict(r) for r in _bodies_by_system(conn, "planet", where_sql, params, query_for, limit, after)]
    where_sql, params = _bodies_after(conn, "planet", where_sql, params, after)
    query = query_for(where_sql)
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params = params + [limit, 0 if after is not None else offset or 0]

    rows = conn.execute(query, params).fetchall()
    return [dict(r) for r in rows]


def list_moons(conn, planet_class=None, min_radius_km=None, max_radius_km=None,
                sector_id=None, system_id=None, limit=None, offset=None, after=None):
    """
    Returns moons, optionally filtered by class, radius range, sector,
    and/or system -- the `moons` counterpart of `list_planets` above (see
    its docstring for the gap this closes).

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        planet_class (str, optional): Exact `planet_class` match (e.g. `"M"`).
        min_radius_km (float, optional): Only moons at least this large.
        max_radius_km (float, optional): Only moons at most this large.
        sector_id (int, optional): Only moons whose system is in this sector.
        system_id (int, optional): Only moons in this one system.
        limit (int, optional): Caps the number of rows returned. `None`
            (the default) returns every matching moon.
        offset (int, optional): Skips this many rows first. Ignored
            unless `limit` is also given, and when `after` is.
        after (int, optional): The row id of the moon the previous page ended on: the page starts after it in the order
            (system name, system, its planet's orbit and id, its orbit, id) and `offset` is not used (PERF.207).

    Returns:
        list[dict]: One row per matching moon, with `id`, `name`,
            `planet_class`, `body_type`, `radius_km`, `star_system_id`,
            `system_name`, `planet_name` (its parent planet's name).
    """
    where_sql, params = _body_filter_clause("m", planet_class, min_radius_km, max_radius_km, sector_id, system_id)

    def query_for(where):
        return f"""
            SELECT m.id, m.name, m.planet_class, m.body_type, m.radius_km, m.star_system_id,
                   ss.name AS system_name, p.name AS planet_name
            FROM moons m
            JOIN planets p ON p.id = m.planet_id
            JOIN star_systems ss ON ss.id = m.star_system_id{where}
            ORDER BY ss.name, ss.id, p.orbital_index, p.id, m.orbital_index, m.id
            """

    if limit is not None and not params and (after is not None or not offset):
        return [dict(r) for r in _bodies_by_system(conn, "moon", where_sql, params, query_for, limit, after)]
    where_sql, params = _bodies_after(conn, "moon", where_sql, params, after)
    query = query_for(where_sql)
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params = params + [limit, 0 if after is not None else offset or 0]

    rows = conn.execute(query, params).fetchall()
    return [dict(r) for r in rows]


class NavUnavailable(Exception):
    """
    Raised by `nav_between` when NAV is unavailable between two systems --
    per the rule set NAV was designed against: either system isn't
    assigned to a sector at all, or the two systems are in different
    sectors and at least one of those sectors has no galaxy placement.
    Distinct from a plain `ValueError` (raised for a system id that
    doesn't exist at all -- see `_load_nav_endpoint`), so the API layer
    can tell "no such system" (404) apart from "these two exist but NAV
    doesn't apply to this pair" (400) without string-matching a message.
    """


def _load_nav_endpoint(conn, system_id):
    """
    Loads the sector placement and position (both sector-local and, when
    the sector itself is galaxy-placed, absolute galaxy-frame) needed to
    resolve one end of a NAV request.

    There is no existing function that combines a sector's galaxy-frame
    center (`sectors.center_x/y/z_pc`, parsecs) with a system's
    sector-local offset (`star_systems.position_x/y/z_mpc`,
    milliparsecs) into one absolute position -- both get converted to
    light-years (`pc_to_ly`/`milliparsecs_to_ly`) and summed componentwise
    here, since light-years is the unit `planetgen.galaxy.navigation`
    already works in for sector-local distances (see
    `spaceSector.distance_between`).

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        system_id (int): The `star_systems.id` to look up.

    Returns:
        dict: `sector_id` (int or None), `position_ly` (sector-local
            `(x, y, z)` tuple, or None if unplaced), `galaxy_position_ly`
            (absolute `(x, y, z)` tuple, or None if `sector_id` is None or
            that sector has no galaxy placement).

    Raises:
        ValueError: If `system_id` doesn't exist.
    """
    row = conn.execute(
        """
        SELECT ss.sector_id, ss.position_x_mpc, ss.position_y_mpc, ss.position_z_mpc,
               sec.center_x_pc, sec.center_y_pc, sec.center_z_pc
        FROM star_systems ss
        LEFT JOIN sectors sec ON sec.id = ss.sector_id
        WHERE ss.id = ?
        """,
        (system_id,),
    ).fetchone()
    if row is None:
        raise ValueError(f"no star_systems row with id {system_id}")

    if row["position_x_mpc"] is None:
        return {"sector_id": row["sector_id"], "position_ly": None, "galaxy_position_ly": None}

    position_ly = (
        milliparsecs_to_ly(row["position_x_mpc"]),
        milliparsecs_to_ly(row["position_y_mpc"]),
        milliparsecs_to_ly(row["position_z_mpc"]),
    )

    galaxy_position_ly = None
    if row["center_x_pc"] is not None:
        galaxy_position_ly = (
            pc_to_ly(row["center_x_pc"]) + position_ly[0],
            pc_to_ly(row["center_y_pc"]) + position_ly[1],
            pc_to_ly(row["center_z_pc"]) + position_ly[2],
        )

    return {"sector_id": row["sector_id"], "position_ly": position_ly, "galaxy_position_ly": galaxy_position_ly}


def _phenomenon_nav_key(phenomenon_type, phenomenon_id):
    """
    The node id a phenomenon endpoint uses in `nav_between`'s position/
    adjacency-graph dicts and `route["path"]` -- a plain string (not a
    tuple) specifically because `route`/`route["positions"]` round-trip
    through `jsonify` (`/api/nav`), and a dict with a tuple key isn't
    JSON-serializable at all, while a `star_systems.id` int key already
    survives that round trip (JSON object keys are always strings, and
    `json.dumps` stringifies an int key for free). The `"phenomenon:"`
    prefix can never collide with a stringified system id.

    Args:
        phenomenon_type (str): One of `_PHENOMENON_TYPE_TO_TABLE`'s keys.
        phenomenon_id (int): The phenomenon's own row id.

    Returns:
        str: e.g. `"phenomenon:nebula:5"`.
    """
    return f"phenomenon:{phenomenon_type}:{phenomenon_id}"


def _load_nav_phenomenon_endpoint(conn, phenomenon_type, phenomenon_id):
    """
    Loads the galaxy-frame position needed to resolve a phenomenon as one
    end of a NAV request -- the phenomenon counterpart to
    `_load_nav_endpoint`, returning the same `sector_id`/`position_ly`/
    `galaxy_position_ly` shape so `nav_between` can treat either kind of
    endpoint identically from that point on.

    A phenomenon's own `sector_id` column (see `_PHENOMENON_TABLES`'s own
    v18/v21 header notes) is only ever a "nearest already-generated
    sector" convenience link, not real containment the way a system's
    `sector_id` foreign key is -- a phenomenon has no sector-LOCAL
    position at all, only a galaxy-frame one, so it can never qualify for
    sector-scope NAV. `sector_id`/`position_ly` are therefore always
    `None` here regardless of that convenience column, which is exactly
    what makes `nav_between`'s own sector-scope eligibility check
    correctly exclude a phenomenon endpoint without needing a separate
    "is this a phenomenon" flag anywhere in that logic.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        phenomenon_type (str): One of `_PHENOMENON_TYPE_TO_TABLE`'s keys.
        phenomenon_id (int): The phenomenon's own row id.

    Returns:
        dict: `sector_id` (always `None`), `position_ly` (always `None`),
            `galaxy_position_ly` (`(x, y, z)` tuple in light-years, or
            `None` if this phenomenon has never been placed in the galaxy).

    Raises:
        ValueError: If `phenomenon_type` is unrecognized, or no such row
            exists.
    """
    table = _PHENOMENON_TYPE_TO_TABLE.get(phenomenon_type)
    if table is None:
        raise ValueError(f"no such phenomenon type: {phenomenon_type!r}")

    row = conn.execute(
        f"SELECT center_x_pc, center_y_pc, center_z_pc FROM {table} WHERE id = ?",
        (phenomenon_id,),
    ).fetchone()
    if row is None:
        raise ValueError(f"no {table} row with id {phenomenon_id}")

    galaxy_position_ly = None
    if row["center_x_pc"] is not None:
        galaxy_position_ly = (
            pc_to_ly(row["center_x_pc"]),
            pc_to_ly(row["center_y_pc"]),
            pc_to_ly(row["center_z_pc"]),
        )
    return {"sector_id": None, "position_ly": None, "galaxy_position_ly": galaxy_position_ly}


def _sector_local_positions(conn, sector_id):
    """
    Returns every placed system's sector-local position (in light-years)
    within one sector -- the position set an in-sector NAV route's
    adjacency graph is built from.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        sector_id (int): The sector to gather positions from.

    Returns:
        dict: `{star_systems.id: (x, y, z)}`, light-years, sector-local.
    """
    rows = conn.execute(
        "SELECT id, position_x_mpc, position_y_mpc, position_z_mpc FROM star_systems "
        "WHERE sector_id = ? AND position_x_mpc IS NOT NULL",
        (sector_id,),
    ).fetchall()
    return {
        row["id"]: (
            milliparsecs_to_ly(row["position_x_mpc"]),
            milliparsecs_to_ly(row["position_y_mpc"]),
            milliparsecs_to_ly(row["position_z_mpc"]),
        )
        for row in rows
    }


def _longest_hop(path, positions):
    return max((detmath.dist(positions[a], positions[b]) for a, b in zip(path, path[1:])), default=0.0)


def _galaxy_route(conn, origin_position, destination_position, from_key, to_key, adjacency_k, unknown_edges=None):
    """
    The cross-sector route between two galaxy-frame positions, searched among
    the systems near the straight line between them instead of every placed
    system (NAV.10): the corridor starts at `NAV_CORRIDOR_FRACTION` of the
    direct distance either side of the line (clamped to
    `NAV_CORRIDOR_MIN_LY`..`NAV_CORRIDOR_START_MAX_LY`) and doubles, up to
    `NAV_CORRIDOR_MAX_LY`, while the route found has a hop longer than half
    the half-width -- a hop that long may be hugging the corridor's edge with
    stepping stones just outside it. A corridor holding more than
    `NAV_CORRIDOR_MAX_SYSTEMS` is not widened. The route is found with A*.

    With `unknown_edges` (`corridor.UnknownEdges`, NAV.48's bypass test) a hop through an uncharted
    sector may not be used, so the answer is `None` when there is no route through known space.

    Returns:
        tuple: `(positions, found)`: the positions the graph was built from
            (`{node_id: (x, y, z)}`, light-years) and `shortest_path`'s answer.
    """
    direct_ly = detmath.dist(origin_position, destination_position)
    half_width = min(max(direct_ly * NAV_CORRIDOR_FRACTION, NAV_CORRIDOR_MIN_LY), NAV_CORRIDOR_START_MAX_LY)
    while True:
        positions = positions_near_segment(conn, origin_position, destination_position, half_width)
        # The ends are nodes whatever they are: a phenomenon is no star_systems
        # row, and a system's own position is the one NAV measured.
        positions[from_key] = origin_position
        positions[to_key] = destination_position
        graph = build_route_graph(positions, adjacency_k, NAV_ISLAND_LINKS)
        found = shortest_path(
            graph, from_key, to_key, positions,
            blocked=None if unknown_edges is None else (lambda a, b: unknown_edges(positions[a], positions[b])))
        settled = (found is not None and _longest_hop(found[0], positions) <= half_width / 2.0)
        if settled or half_width >= NAV_CORRIDOR_MAX_LY or len(positions) > NAV_CORRIDOR_MAX_SYSTEMS:
            return positions, found
        half_width = min(half_width * 2.0, NAV_CORRIDOR_MAX_LY)


def _with_stones(conn, path, positions, distance_ly):
    """
    NAV.47: `path` with every hop through unknown space broken up by the stepping stones
    `stones.chain_for_jump` finds for it. Each stone's position is added to `positions`.

    Returns:
        tuple: `(path, distance_ly, {key: stone})`; the arguments' own values when nothing changed.
    """
    cells = unknown_cells(conn, [positions[node] for node in path])
    if not any(cells):
        return path, distance_ly, {}
    new_path, found = [path[0]], {}
    for (a, b), crossed in zip(zip(path, path[1:]), cells):
        if crossed:
            for key, stone in nav_stones.chain_for_jump(conn, positions[a], positions[b], crossed):
                positions[key] = stone["position_ly"]
                found[key] = stone
                new_path.append(key)
        new_path.append(b)
    if not found:
        return path, distance_ly, {}
    return new_path, sum(detmath.dist(positions[a], positions[b]) for a, b in zip(new_path, new_path[1:])), found


def _route_places(conn, path):
    """
    NAV.42: where each stop of a route sits in its sector. `{node: {"sector_id", "sector_position_ly"}}`, the
    position the system's sector-local `(x, y, z)` in light-years; both `None` for a phenomenon end and for a system
    with no sector or no position.
    """
    places = {node: {"sector_id": None, "sector_position_ly": None} for node in path}
    system_ids = [node for node in path if not isinstance(node, str)]
    if system_ids:
        marks = ",".join("?" * len(system_ids))
        for row in conn.execute(
                "SELECT id, sector_id, position_x_mpc, position_y_mpc, position_z_mpc FROM star_systems "
                f"WHERE id IN ({marks})", system_ids).fetchall():
            places[row["id"]]["sector_id"] = row["sector_id"]
            if row["position_x_mpc"] is not None:
                places[row["id"]]["sector_position_ly"] = (
                    milliparsecs_to_ly(row["position_x_mpc"]), milliparsecs_to_ly(row["position_y_mpc"]),
                    milliparsecs_to_ly(row["position_z_mpc"]))
    return places


def _hop_course(hop, places, galaxy_frame, local_positions):
    """
    NAV.42: the course from one stop of a route to the next, in the frame NAV uses for that pair -- the Sector
    Local Frame when both stops are systems of one sector, the Galactic Frame otherwise
    (docs/design/orbits/navigation-frames.md).

    Returns:
        dict: `bearing_deg`, `mark_deg`, `elevation_deg` and `frame`.
    """
    a, b = places[hop["from"]], places[hop["to"]]
    if (a["sector_id"] is not None and a["sector_id"] == b["sector_id"]
            and a["sector_position_ly"] is not None and b["sector_position_ly"] is not None):
        course = course_between(a["sector_position_ly"], b["sector_position_ly"], frame=FRAME_SECTOR)
    elif galaxy_frame is not None:
        course = course_between(galaxy_frame[hop["from"]], galaxy_frame[hop["to"]], frame=FRAME_GALACTIC)
    else:  # an in-sector search over sector-local positions
        course = course_between(local_positions[hop["from"]], local_positions[hop["to"]], frame=FRAME_SECTOR)
    return {"bearing_deg": course.bearing_deg, "mark_deg": course.mark_deg, "elevation_deg": course.elevation_deg,
            "frame": course.frame}


def _bypass_route(conn, origin_position, destination_position, from_key, to_key, adjacency_k):
    """
    NAV.48's bypass test: whether a route exists that crosses no uncharted sector, and how long it is.

    Returns:
        dict: `found` (a route through known space exists), `distance_ly` (its length, else `None`) and
            `checked` (false when the search gave up after `tuning.NAV_BYPASS_MAX_EDGES` hops: the answer is
            then unknown and `found` is false).
    """
    try:
        _positions, found = _galaxy_route(conn, origin_position, destination_position, from_key, to_key,
                                          adjacency_k, unknown_edges=UnknownEdges(conn, NAV_BYPASS_MAX_EDGES))
    except BypassBudgetExceeded:
        return {"found": False, "distance_ly": None, "checked": False}
    return {"found": found is not None, "distance_ly": found[1] if found is not None else None, "checked": True}


def nav_between(conn, from_id, to_id, adjacency_k=NAV_ADJACENCY_K,
                 from_kind="system", to_kind="system", from_type=None, to_type=None, stay_minutes=0.0,
                 chart=False, stones=True):
    """
    Resolves full NAV information between two endpoints -- each either a
    star system or a standalone phenomenon (nebula/asteroid field/black
    hole/neutron star) -- a direct course (distance/bearing/mark/warp and fold
    travel times, from `planetgen.galaxy.navigation`) plus an optimal route
    via adjacent systems (`planetgen.galaxy.nav_graph`), or raises if NAV
    doesn't apply to this pair.

    NAV availability rules (see docs/reference/api.md's NAV section for the
    user-facing statement of these):
        - Both endpoints are systems in the SAME sector -> available,
          scoped to that sector's own systems (sector-local positions) --
          checked first/preferred, same as before this function supported
          phenomenon endpoints at all.
        - Otherwise, both endpoints have a galaxy-frame position (a system
          in a galaxy-placed sector, or a galaxy-placed phenomenon) ->
          available at galaxy scope. A phenomenon endpoint can only ever
          reach this branch: it has no sector-LOCAL position at all (see
          `_load_nav_phenomenon_endpoint`), so it never qualifies for
          sector scope regardless of its own "nearest sector" convenience
          link.
        - Anything else (an unplaced system, a phenomenon never placed in
          the galaxy, or two systems in different sectors where either
          lacks a galaxy placement) -> unavailable.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        from_id (int): The origin's own row id -- `star_systems.id` when
            `from_kind == "system"`, else the phenomenon's own table id.
        to_id (int): Same, for the destination.
        adjacency_k (int): Passed through to
            `navGraph.build_knn_adjacency` as `k`; the islands that
            graph leaves are then joined (`navGraph.join_islands`,
            `NAV_ISLAND_LINKS` nearest islands each).
        from_kind (str): `"system"` (default) or `"phenomenon"`.
        to_kind (str): Same, for the destination.
        from_type (str, optional): Required when `from_kind ==
            "phenomenon"` -- one of `_PHENOMENON_TYPE_TO_TABLE`'s keys.
        to_type (str, optional): Same, for the destination.
        stay_minutes (float): NAV.11: the time spent at each system the
            route stops at between the two ends (default 0); only the
            route's total times use it.
        chart (bool): NAV.48: each hop through unknown space also lists its
            uncharted `unknown_cells`, and a route with such a hop gets a
            `bypass` entry (`_bypass_route`: whether a route through known
            space exists). `nav_chart_plan` asks for these.
        stones (bool): NAV.47: a hop through unknown space is broken up by stops
            at scattered stars, black holes, neutron stars, quasars and
            hypervelocity stars in the uncharted cells when a chain of them
            exists (`stones.chain_for_jump`); each such stop's `stop_places`
            entry names its `kind` and `name`. `nav_chart_plan` turns this off,
            since it asks which cells lie on the direct jump.

    Returns:
        dict: `scope` (`"sector"` or `"galaxy"`), `direct` (a
            `navigation.Course`), `warp_times` (a list of
            `navigation.WarpLeg`, for `direct.distance_ly`), `fold_times`
            (a list of `navigation.FoldLeg`, same distance),
            `origin_position`/`destination_position` (the `(x, y, z)`
            light-year positions `direct` was computed from, in `scope`'s
            frame -- sector-local for `"sector"`, galaxy-frame for
            `"galaxy"`), and `route`: `None` if the two endpoints resolve
            to the same node (the joined graph always has a path
            otherwise), else `{"path": [...node ids...], "distance_ly": float,
            "positions": {node_id: (x, y, z), ...}, "hops": [{"from", "to",
            "distance_ly", "unknown_space", "bearing_deg", "mark_deg", "elevation_deg", "frame" (NAV.42: the
            hop's course), "warp_times", "fold_times"}], "stop_places": [{"node", "sector_id",
            "sector_position_ly"}] (NAV.42: one per id in `path`),
            "longest_hop_ly": float, "stops": int, "stay_minutes": float,
            "warp_times", "fold_times"}` (NAV.11: a hop's lists are
            `WarpLeg`/`FoldLeg` dicts for that hop alone; the route's are
            `WarpLeg`/`FoldLeg` for the whole route, a stop at every system
            between the ends for `stay_minutes` each) (one entry per id in
            `path`, same frame as `origin_position`/`destination_position`
            -- for rendering the route, e.g. `planetgen/web/maps/navmap.py`, without
            a second position lookup). A node id is a plain `star_systems.
            id` int for a system hop (every INTERMEDIATE hop always is,
            regardless of either endpoint's own kind -- only `path[0]`/
            `path[-1]` can ever be a phenomenon), or a
            `_phenomenon_nav_key`-shaped string for a phenomenon endpoint.

    Raises:
        ValueError: If an endpoint id (or `from_type`/`to_type`) doesn't
            resolve to a real row.
        NavUnavailable: If NAV doesn't apply to this pair, per the rules
            above. `str(exc)` explains why.
    """
    def resolve(ref_id, kind, phenomenon_type):
        if kind == "phenomenon":
            return _load_nav_phenomenon_endpoint(conn, phenomenon_type, ref_id)
        return _load_nav_endpoint(conn, ref_id)

    def node_key(ref_id, kind, phenomenon_type):
        return ref_id if kind == "system" else _phenomenon_nav_key(phenomenon_type, ref_id)

    origin = resolve(from_id, from_kind, from_type)
    destination = resolve(to_id, to_kind, to_type)
    from_key = node_key(from_id, from_kind, from_type)
    to_key = node_key(to_id, to_kind, to_type)

    sector_scope_ok = (
        origin["sector_id"] is not None
        and origin["sector_id"] == destination["sector_id"]
        and origin["position_ly"] is not None
        and destination["position_ly"] is not None
    )
    galaxy_scope_ok = (
        origin["galaxy_position_ly"] is not None
        and destination["galaxy_position_ly"] is not None
    )

    sector_offset = None
    if sector_scope_ok:
        scope = "sector"
        origin_position, destination_position = origin["position_ly"], destination["position_ly"]
        if galaxy_scope_ok:
            # NAV.12: the route may leave the sector, so it is searched in the
            # galaxy frame and shifted back into this sector's local frame.
            sector_offset = tuple(g - l for g, l in zip(origin["galaxy_position_ly"], origin["position_ly"]))
            positions = None
        else:
            positions = _sector_local_positions(conn, origin["sector_id"])
    elif galaxy_scope_ok:
        scope = "galaxy"
        origin_position, destination_position = origin["galaxy_position_ly"], destination["galaxy_position_ly"]
        positions = None  # read below, from the corridor round the direct line (NAV.10)
    else:
        raise NavUnavailable(
            "NAV requires both endpoints to share a sector, or both to have a galaxy placement"
        )

    # Both frames are centered on their coordinate origin: sector-local
    # positions are offsets from the sector's center, galaxy-frame ones
    # from the galactic core. Bearing 000 points at that center.
    frame = FRAME_SECTOR if scope == "sector" else FRAME_GALACTIC
    direct = course_between(origin_position, destination_position, frame=frame)

    route = None
    if from_key != to_key:
        # The graph's islands are joined (NAV.34), so two placed
        # endpoints always have a route.
        if positions is not None:
            graph = build_route_graph(positions, adjacency_k, NAV_ISLAND_LINKS)
            found = shortest_path(graph, from_key, to_key, positions)
            galaxy_frame = None
        else:
            search_from, search_to = ((origin["galaxy_position_ly"], destination["galaxy_position_ly"])
                                      if sector_offset is not None else (origin_position, destination_position))
            positions, found = _galaxy_route(conn, search_from, search_to, from_key, to_key, adjacency_k)
            galaxy_frame = positions
            search_ends = (search_from, search_to)
        if found is not None:
            path, distance_ly = found
            stops_found = {}
            if stones and galaxy_frame is not None:
                path, distance_ly, stops_found = _with_stones(conn, path, positions, distance_ly)
            route_positions = {node_id: positions[node_id] for node_id in path}
            hops = [{"from": a, "to": b, "distance_ly": detmath.dist(positions[a], positions[b]),
                     "unknown_space": False} for a, b in zip(path, path[1:])]
            places = _route_places(conn, path)
            for hop in hops:
                hop.update(_hop_course(hop, places, galaxy_frame, positions))
                hop["warp_times"] = [leg._asdict() for leg in warp_travel_times(hop["distance_ly"])]
                hop["fold_times"] = [leg._asdict() for leg in fold_travel_times(hop["distance_ly"])]
            bypass = None
            if galaxy_frame is not None:
                cells = unknown_cells(conn, [galaxy_frame[node_id] for node_id in path])
                for hop, hop_cells in zip(hops, cells):
                    # A hop to or from a stop in unknown space (NAV.47) is in it whatever cells it crosses.
                    hop["unknown_space"] = bool(hop_cells) or nav_stones.is_stone(hop["from"]) \
                        or nav_stones.is_stone(hop["to"])
                    if chart:
                        hop["unknown_cells"] = hop_cells
                if chart and any(cells):
                    bypass = _bypass_route(conn, *search_ends, from_key, to_key, adjacency_k)
                if sector_offset is not None:
                    route_positions = {node_id: tuple(c - o for c, o in zip(position, sector_offset))
                                       for node_id, position in route_positions.items()}
            stays = max(len(hops) - 1 - len(stops_found), 0)  # not at a stop in unknown space (NAV.47)
            route = {
                "path": path,
                "distance_ly": distance_ly,
                "positions": route_positions,
                "stop_places": [
                    {"node": node_id, **places[node_id],
                     **({"kind": stops_found[node_id]["kind"], "name": nav_stones.describe(stops_found[node_id])}
                        if node_id in stops_found else {})}
                    for node_id in path],
                "hops": hops,
                "longest_hop_ly": max((hop["distance_ly"] for hop in hops), default=0.0),
                "stops": max(len(hops) - 1, 0),
                "stay_minutes": stay_minutes,
                "warp_times": route_warp_times([hop["distance_ly"] for hop in hops], stay_minutes, stays=stays),
                "fold_times": route_fold_times([hop["distance_ly"] for hop in hops], stay_minutes, stays=stays),
            }
            if bypass is not None:
                route["bypass"] = bypass

    return {
        "scope": scope,
        "direct": direct,
        "warp_times": warp_travel_times(direct.distance_ly),
        "fold_times": fold_travel_times(direct.distance_ly),
        "origin_position": origin_position,
        "destination_position": destination_position,
        "route": route,
    }


NAV_DEFAULT_HELIOPAUSE_AU = 120.0
"""float: The heliopause used for a system that stores none (Boss's text
in `docs/design/orbits/navigation-frames.md`)."""

NAV_IN_SYSTEM_NOTE = (
    "Travel times use the same warp and fold tables as courses between systems; "
    "they are not tuned for distances inside a system."
)
"""str: Shown with a course that has a leg inside a system (NAV.16)."""


def _system_heliopause_km(conn, system_id):
    """A system's heliopause radius in km, as `system_detail` reports it
    (pressed in by a surrounding cloud), else `NAV_DEFAULT_HELIOPAUSE_AU`."""
    row = conn.execute("SELECT * FROM star_systems WHERE id = ?", (system_id,)).fetchone()
    _open, heliopause_au = _heliopause_au(conn, row, containing_cloud(conn, row))
    return _au_to_km(heliopause_au if heliopause_au is not None else NAV_DEFAULT_HELIOPAUSE_AU)


def _km_to_ly(km):
    return _m_to_ly(_km_to_m(km))


def _nav_endpoint(conn, ref):
    """
    Resolves one `from`/`to` reference for `nav_course`.

    Returns:
        dict: `ref`, `kind`, `id`, `name`, `anchor` (`(kind, id, type)` as
            `nav_between` takes it: a system or a phenomenon) and
            `system_km` (a body's system-local position in km, `None` for a
            system, a phenomenon, or anything NAV treats as the system itself).

    Raises:
        ValueError: For a bad reference or a missing row.
        NavUnavailable: For a sector (not a place a ship can go).
    """
    kind, object_id = object_ref.parse(ref)
    if kind == "sector":
        raise NavUnavailable("a sector is not an endpoint; pick a system or something in one")
    if kind in object_ref.CONTAINER_KINDS:
        raise NavUnavailable("a cluster or a superbubble is not an endpoint; pick a system or something in one")
    obj = resolve_object(conn, kind, object_id)
    if kind in object_ref.PHENOMENON_KINDS:
        anchor = ("phenomenon", object_id, kind)
    elif kind == "system":
        anchor = ("system", object_id, None)
    else:
        system = next(parent for parent in obj["parents"] if parent["kind"] == "system")
        anchor = ("system", object_ref.parse(system["ref"])[1], None)
    system_km = None
    if kind in object_ref.BODY_KINDS:
        system_km = obj["positions"]["system_km"]
        if system_km is None:  # an asteroid belt is a ring: take its mean radius
            row = conn.execute("SELECT distance_km FROM asteroid_belts WHERE id = ?", (object_id,)).fetchone()
            system_km = [row["distance_km"], 0.0, 0.0]
    return {"ref": obj["ref"], "kind": kind, "id": object_id, "name": obj["name"],
            "anchor": anchor, "system_km": system_km}


def _nav_leg(kind, origin_ref, destination_ref, course):
    return {
        "kind": kind, "from": origin_ref, "to": destination_ref, "direct": course,
        "warp_times": warp_travel_times(course.distance_ly), "fold_times": fold_travel_times(course.distance_ly),
    }


def nav_chart_plan(conn, from_ref, to_ref, border=False, adjacency_k=NAV_ADJACENCY_K):
    """
    NAV.48: the uncharted sectors that block a course, to be generated so it can be plotted again.

    The route is the one `nav_course` finds; its hops through unknown space name the cells to chart
    (the cells of those hops only, not every cell on the straight line), left out where they lie
    outside the galaxy's outline. `border` adds the cells sharing a face with those that are not
    charted either.

    Returns:
        dict: `unknown_hops` (how many), `cells` (`(shell, band, wedge)` to chart, in the order the
            route enters them), `outside_galaxy` (uncharted cells past the outline, never offered),
            `route_distance_ly` and `bypass` (`_bypass_route`'s answer, `None` when no hop is unknown).

    Raises:
        ValueError: For a bad reference or a missing row.
        NavUnavailable: As `nav_between`, and for a sector endpoint.
    """
    origin, destination = _nav_endpoint(conn, from_ref), _nav_endpoint(conn, to_ref)
    if origin["anchor"] == destination["anchor"] and origin["anchor"][0] == "system":
        return {"unknown_hops": 0, "cells": [], "outside_galaxy": 0, "route_distance_ly": 0.0, "bypass": None}
    from_kind, from_id, from_type = origin["anchor"]
    to_kind, to_id, to_type = destination["anchor"]
    result = nav_between(conn, from_id, to_id, adjacency_k, from_kind=from_kind, to_kind=to_kind,
                         from_type=from_type, to_type=to_type, chart=True, stones=False)
    route = result["route"]
    hops = [hop for hop in (route["hops"] if route else []) if hop.get("unknown_space")]
    bounds = get_galaxy_bounds(conn)
    cells = list(dict.fromkeys(cell for hop in hops for cell in hop["unknown_cells"]))
    inside = [cell for cell in cells if bounds.qualifies(cell[0], cell[1])]
    outside = len(cells) - len(inside)
    if border and inside:
        listed = set(inside)
        around = dict.fromkeys(
            near for cell in inside for near in neighbor_addresses(*cell)
            if bounds.qualifies(near[0], near[1]) and near not in listed)
        charted = generated_cells(conn, set(around))
        inside += [cell for cell in around if cell not in charted]
    return {"unknown_hops": len(hops), "cells": inside, "outside_galaxy": outside,
            "route_distance_ly": route["distance_ly"] if route else 0.0,
            "bypass": route.get("bypass") if route else None}


def nav_course(conn, from_ref, to_ref, adjacency_k=NAV_ADJACENCY_K, stay_minutes=0.0):
    """
    NAV between any two objects (NAV.16): `nav_between`'s result plus
    `origin`/`destination` (`{ref, kind, name}`), `legs` and `note`.

    Two systems or phenomena are `nav_between`'s course, one `"between"`
    leg. A body (star, planet, moon, belt, comet) adds an `"out"` leg from
    it to its system's heliopause on the side facing the destination (a
    destination body, an `"into"` leg from the heliopause on the side
    facing the origin to the body), each in the System Local Frame. Two
    objects in one system are a single `"within"` leg (`scope`
    `"system"`; positions in light-years from the system's origin, no
    route). Each leg is `{kind, from, to, direct, warp_times, fold_times}`;
    `direct` is a `navigation.Course`.

    `stay_minutes` (NAV.11) is the stay at each stop of the route.

    Raises:
        ValueError: For a bad reference or a missing row.
        NavUnavailable: As `nav_between`, and for a sector endpoint.
    """
    origin, destination = _nav_endpoint(conn, from_ref), _nav_endpoint(conn, to_ref)
    ends = {"origin": {k: origin[k] for k in ("ref", "kind", "name")},
            "destination": {k: destination[k] for k in ("ref", "kind", "name")}}

    if origin["anchor"] == destination["anchor"] and origin["anchor"][0] == "system":
        # Inside one system (a system itself sits at the origin).
        a = [_km_to_ly(v) for v in origin["system_km"] or (0.0, 0.0, 0.0)]
        b = [_km_to_ly(v) for v in destination["system_km"] or (0.0, 0.0, 0.0)]
        course = course_between(a, b, frame=FRAME_SYSTEM)
        return {
            "scope": "system", "direct": course,
            "warp_times": warp_travel_times(course.distance_ly), "fold_times": fold_travel_times(course.distance_ly),
            "origin_position": tuple(a), "destination_position": tuple(b), "route": None,
            "legs": [_nav_leg("within", origin["ref"], destination["ref"], course)],
            "note": NAV_IN_SYSTEM_NOTE, **ends,
        }

    def anchor_args(point):
        kind, anchor_id, anchor_type = point["anchor"]
        return anchor_id, kind, anchor_type

    from_id, from_kind, from_type = anchor_args(origin)
    to_id, to_kind, to_type = anchor_args(destination)
    result = nav_between(conn, from_id, to_id, adjacency_k, from_kind=from_kind, to_kind=to_kind,
                         from_type=from_type, to_type=to_type, stay_minutes=stay_minutes)

    toward = [d - o for o, d in zip(result["origin_position"], result["destination_position"])]
    length = math.sqrt(sum(v * v for v in toward))
    heading = [v / length for v in toward] if length else [1.0, 0.0, 0.0]

    legs = []
    if origin["system_km"] is not None:
        edge = [v * _system_heliopause_km(conn, from_id) for v in heading]
        legs.append(_nav_leg("out", origin["ref"], object_ref.format("system", from_id), course_between(
            [_km_to_ly(v) for v in origin["system_km"]], [_km_to_ly(v) for v in edge], frame=FRAME_SYSTEM)))
    legs.append(_nav_leg("between", object_ref.format(from_kind if from_kind == "system" else from_type, from_id),
                         object_ref.format(to_kind if to_kind == "system" else to_type, to_id), result["direct"]))
    if destination["system_km"] is not None:
        edge = [-v * _system_heliopause_km(conn, to_id) for v in heading]
        legs.append(_nav_leg("into", object_ref.format("system", to_id), destination["ref"], course_between(
            [_km_to_ly(v) for v in edge], [_km_to_ly(v) for v in destination["system_km"]], frame=FRAME_SYSTEM)))
    return {**result, "legs": legs, "note": NAV_IN_SYSTEM_NOTE if len(legs) > 1 else None, **ends}


def sector_address(row):
    """A `sectors` row's `(shell_index, band_index, wedge_index)`, or
    `None` when it has no grid address (never placed, or hand-placed)."""
    address = (row["shell_index"], row["band_index"], row["wedge_index"])
    return None if any(value is None for value in address) else address


def sector_neighbors(conn, sector):
    """
    Every cell sharing a face with a galaxy-placed sector
    (`spherical.neighbor_addresses`: the wedges either side, the
    bands above and below, and the overlapping wedges in the shells inside
    and outside), each tagged with whether a real `sectors` row already
    exists there. Drives the Sector Map's (`planetgen/web/maps/starmap.py`)
    neighboring-sector indicators -- an existing neighbor links straight
    to it; a not-yet-generated one shows its address so it can be fed to
    `planetgen galaxy --shell I --band J --wedge K`.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        sector (dict or Row): This sector's own row -- needs
            `shell_index`, `band_index`, `wedge_index`, and `edge_mpc`.

    Returns:
        list[dict]: `shell_index`, `band_index`, `wedge_index`,
            `direction_pc` (`[x, y, z]`, this neighbor's galaxy-frame
            center minus this sector's own), `exists` (bool),
            `sector_id`/`sector_name` (both `None` when `exists` is
            `False`), and `designation` (`provisional_sector_designation`).
            `[]` for a sector with no grid address.
    """
    address = sector_address(sector)
    if address is None or not sector["edge_mpc"]:
        return []

    edge_pc = mpc_to_pc(sector["edge_mpc"])
    this_position = sector_position_pc(*address, edge_pc)
    addresses = neighbor_addresses(*address)

    clauses = " OR ".join(["(shell_index = ? AND band_index = ? AND wedge_index = ?)"] * len(addresses))
    params = [value for neighbor in addresses for value in neighbor]
    rows = conn.execute(
        f"SELECT id, name, shell_index, band_index, wedge_index FROM sectors WHERE {clauses}", tuple(params),
    ).fetchall()
    existing = {(r["shell_index"], r["band_index"], r["wedge_index"]): r for r in rows}

    results = []
    for neighbor in addresses:
        position = sector_position_pc(*neighbor, edge_pc)
        match = existing.get(neighbor)
        results.append({
            "shell_index": neighbor[0], "band_index": neighbor[1], "wedge_index": neighbor[2],
            "direction_pc": [position[i] - this_position[i] for i in range(3)],
            "exists": match is not None,
            "sector_id": match["id"] if match else None,
            "sector_name": match["name"] if match else None,
            "designation": provisional_sector_designation(*neighbor),
        })
    return results


def _center_distance_ly(row):
    """A `star_systems` row's distance from its sector's center, in
    light-years (its `position_*_mpc` is already center-relative), or
    `None` if it was never placed."""
    if row["position_x_mpc"] is None:
        return None
    return milliparsecs_to_ly(
        math.sqrt(row["position_x_mpc"] ** 2 + row["position_y_mpc"] ** 2 + row["position_z_mpc"] ** 2)
    )


_FACILITY_SELECT = """
    SELECT f.*,
           COALESCE(st.name, p.name, m.name, af.name, sec.name) AS host_name,
           ft.name AS type_name, ft.icon AS type_icon
    FROM facilities f
    LEFT JOIN facility_types ft ON ft.id = f.type_id
    LEFT JOIN stars st ON st.id = f.star_id
    LEFT JOIN planets p ON p.id = f.planet_id
    LEFT JOIN moons m ON m.id = f.moon_id
    LEFT JOIN asteroid_fields af ON af.id = f.asteroid_field_id
    LEFT JOIN sectors sec ON sec.id = f.sector_id
"""


def _facility_dict(row):
    """One `facilities` row (plus `host_name`, `None` for an asteroid
    belt, which has no name) as the API returns it."""
    facility = {key: row[key] for key in (
        "id", "name", "kind", "placement", "host_type", "star_system_id", "star_id", "planet_id", "moon_id",
        "asteroid_belt_id", "asteroid_field_id", "sector_id", "center_x_pc", "center_y_pc", "center_z_pc",
        "orbit_distance_km", "orbit_period_years", "orbital_speed_kms", "orbit_phase_deg", "description",
        "velocity_x_kms", "velocity_y_kms", "velocity_z_kms", "host_name",
    )}
    host_columns = {"star": "star_id", "planet": "planet_id", "moon": "moon_id", "asteroid_belt": "asteroid_belt_id",
                    "asteroid_field": "asteroid_field_id", "space": "sector_id"}
    facility["host_id"] = row[host_columns[row["host_type"]]]
    for key in ("type_id", "type_name", "type_icon", "affiliation_polity_id", "affiliation_name"):
        facility[key] = row[key]
    facility.update(facility_ratings.describe(row))
    facility["fields"] = {}
    return facility


def _attach_facility_fields(conn, facilities):
    """Fills each facility's `fields` (`{key: value}` for the values it holds) in one query. Returns `facilities`."""
    if not facilities:
        return facilities
    by_id = {facility["id"]: facility for facility in facilities}
    marks = ", ".join(["?"] * len(by_id))
    for row in conn.execute(
            f"SELECT v.facility_id, d.`key`, d.datatype, v.value_int, v.value_dec, v.value_text, v.value_bool "
            f"FROM facility_field_values v JOIN facility_type_fields d ON d.id = v.field_id "
            f"WHERE v.facility_id IN ({marks}) ORDER BY d.sort_order, d.id", tuple(by_id)).fetchall():
        value = row[facility_types.value_column(row["datatype"])]
        by_id[row["facility_id"]]["fields"][row["key"]] = bool(value) if row["datatype"] == "bool" else value
    return facilities


def facility_detail(conn, facility_id):
    """One facility (schema v42), or `None`."""
    row = conn.execute(_FACILITY_SELECT + " WHERE f.id = ?", (facility_id,)).fetchone()
    return None if row is None else _attach_facility_fields(conn, [_facility_dict(row)])[0]


def facilities_for_system(conn, system_id):
    """Every facility on a host inside one star system (its stars,
    planets, moons and belts), by name."""
    rows = conn.execute(_FACILITY_SELECT + " WHERE f.star_system_id = ? ORDER BY f.name, f.id",
                        (system_id,)).fetchall()
    return _attach_facility_fields(conn, [_facility_dict(row) for row in rows])


def facilities_in_sector(conn, sector_id):
    """The facilities a sector holds outside its systems: stand-alone ones
    parked in it, and those on its asteroid fields."""
    rows = conn.execute(
        _FACILITY_SELECT + " WHERE f.sector_id = ? OR af.sector_id = ? ORDER BY f.name, f.id",
        (sector_id, sector_id),
    ).fetchall()
    return _attach_facility_fields(conn, [_facility_dict(row) for row in rows])


FACILITY_SORTS = ("name", "worst", "kind", "type")
"""tuple: What `list_facilities` sorts by; `worst` puts the worst overall rating first."""

_WORST = ("GREATEST(COALESCE(f.rating_crime, 0), COALESCE(f.rating_housing, 0), COALESCE(f.rating_resources, 0), "
          "COALESCE(f.rating_maintenance, 0), COALESCE(f.rating_health, 0))")


def list_facilities(conn, kind=None, type_id=None, polity_id=None, affiliation_name=None, ratings=None,
                    overall=None, sort="name", descending=False, limit=50, offset=0):
    """
    A page of facilities, optionally filtered, and how many pass the filters.

    Args:
        kind (str, optional): A facility kind.
        type_id (int, optional): A facility type.
        polity_id (int, optional): Facilities affiliated with this polity.
        affiliation_name (str, optional): Facilities with exactly this affiliation text.
        ratings (dict, optional): `{dimension: colour word}` -- each named rating equals that colour.
        overall (str, optional): A colour word -- the worst of the facility's assessed ratings.
        sort (str): One of `FACILITY_SORTS`. `worst` lists the worst overall rating first (ascending order is
            worst first; unassessed facilities come last either way).

    Returns:
        tuple: `(items, total)`.
    """
    where, params = [], []
    for column, value in (("f.kind", kind), ("f.type_id", type_id), ("f.affiliation_polity_id", polity_id),
                          ("f.affiliation_name", affiliation_name)):
        if value is not None:
            where.append(f"{column} = ?")
            params.append(value)
    for name, color in (ratings or {}).items():
        where.append(f"f.rating_{name} = ?")
        params.append(facility_ratings.value_of(color))
    if overall is not None:
        where.append(f"{_WORST} = ?")
        params.append(facility_ratings.value_of(overall))
    clause = f" WHERE {' AND '.join(where)}" if where else ""
    total = conn.execute(f"SELECT COUNT(*) AS n FROM facilities f{clause}", tuple(params)).fetchone()["n"]
    direction = "DESC" if descending else "ASC"
    order = {
        "name": f"f.name {direction}, f.id",
        "kind": f"f.kind {direction}, f.name, f.id",
        "type": f"ft.name IS NULL, ft.name {direction}, f.name, f.id",
        "worst": f"{_WORST} = 0, {_WORST} {'ASC' if descending else 'DESC'}, f.name, f.id",
    }[sort]
    rows = conn.execute(_FACILITY_SELECT + clause + f" ORDER BY {order} LIMIT ? OFFSET ?",
                        (*params, limit, offset)).fetchall()
    return _attach_facility_fields(conn, [_facility_dict(row) for row in rows]), total


def colonized_body_ids(conn, system_id):
    """
    The planets and moons in one system with a colony on them -- a colony
    makes its world inhabited (schema v42).

    Returns:
        dict: `{"planets": set of ids, "moons": set of ids}`.
    """
    found = {"planets": set(), "moons": set()}
    for row in conn.execute(
        "SELECT planet_id, moon_id FROM facilities WHERE star_system_id = ? AND kind = 'colony'"
        " AND placement = 'terrestrial'",
        (system_id,),
    ).fetchall():
        if row["planet_id"] is not None:
            found["planets"].add(row["planet_id"])
        if row["moon_id"] is not None:
            found["moons"].add(row["moon_id"])
    return found


def containing_cloud(conn, row):
    """
    The nebula or supernova remnant a row sits inside (schema v39), from
    its `inside_nebula_id`/`inside_remnant_id`: `store.surrounding_cloud`.

    Args:
        conn (Connection): An open connection.
        row (Mapping): Any row carrying those two columns (a star system,
            rogue planet, interstellar comet, black hole, neutron star,
            asteroid field or nebula).

    Returns:
        dict or None: `{"type": "nebula" | "supernova_remnant", "id",
            "name", "class", "descriptor", "density_cm3", "temperature_k"}`,
            or `None` in open space.
    """
    cloud = surrounding_cloud(conn, row)
    return None if cloud is None else {key: cloud[key] for key in cloudContext.LISTED_KEYS}


def sector_detail(conn, sector_id):
    """
    Returns one sector's full web-display detail: name, size, galaxy
    placement, and every system placed in it (each with its own star
    roster) -- everything `html/sector.py`'s systems table and Sector Map
    (`planetgen/web/maps/starmap.py`) need, in one function.

    Distinct from `planetgen.db.store.load_sector`, which reconstructs
    the *generation* object graph (config/provenance, no database ids) --
    this is a flat, ids-and-display-fields read, the same relationship
    `list_sectors`/`list_systems` above already have to
    `planetgen.db.store.load_star_system`.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        sector_id (int): The `sectors.id` to look up.

    Returns:
        dict: `id`, `name`, `edge_mpc`, `edge_ly`, `center_x_pc`/
            `center_y_pc`/`center_z_pc`, `shell_index`, `band_index`,
            `wedge_index`, `placed`, `system_count`, and `systems` (one entry per system
            placed in this sector, nearest the sector's center first,
            then any with no position by name: `id`, `name`, `quadrant`, `location`,
            `is_binary`, `binary_type`, `binary_separation_km`,
            `position_x_mpc`/`position_y_mpc`/
            `position_z_mpc`, `center_distance_ly` (`None` with no
            position), and `stars` -- 1 entry (single), 2
            (primary then secondary) or more (then each `member`), each `role`/`star_type`/
            `temperature_k`/`radius_km`/`luminosity_w`), `phenomena`
            (see `phenomena_near_sector`), `neighbors` (see
            `sector_neighbors`, `[]` for a sector with no galaxy
            placement), and `stats` (`sector_stats`, PERF.11: `None` off the grid).

    Raises:
        ValueError: If no such sector exists.
    """
    sector = conn.execute("SELECT * FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    if sector is None:
        raise ValueError(f"no sectors row with id {sector_id}")

    system_rows = conn.execute(
        """
        SELECT id, name, is_binary, quadrant, location, binary_type, binary_separation_km,
               position_x_mpc, position_y_mpc, position_z_mpc, runaway_class, runaway_speed_kms,
               inside_nebula_id, inside_remnant_id
        FROM star_systems
        WHERE sector_id = ?
        ORDER BY position_x_mpc IS NULL,
                 position_x_mpc * position_x_mpc + position_y_mpc * position_y_mpc
                     + position_z_mpc * position_z_mpc,
                 name, id
        """,
        (sector_id,),
    ).fetchall()

    nearest = nearest_systems(conn, "star_systems", [row["id"] for row in system_rows])
    stars_by_system = {}
    for star in conn.execute(
        "SELECT st.star_system_id, st.role, st.name, st.star_type, st.temperature_k, st.radius_km, st.luminosity_w,"
        " st.metallicity_fe_h"
        " FROM stars st JOIN star_systems ss ON ss.id = st.star_system_id WHERE ss.sector_id = ?"
        " ORDER BY st.star_system_id, CASE st.role WHEN 'secondary' THEN 1 WHEN 'member' THEN 2 ELSE 0 END, st.id",
        (sector_id,),
    ).fetchall():
        star = dict(star)
        stars_by_system.setdefault(star.pop("star_system_id"), []).append(star)
    systems = []
    for row in system_rows:
        star_rows = stars_by_system.get(row["id"], [])
        systems.append({
            "id": row["id"], "name": row["name"], "quadrant": row["quadrant"], "location": row["location"],
            "is_binary": row["is_binary"], "binary_type": row["binary_type"],
            "binary_separation_km": row["binary_separation_km"],
            "position_x_mpc": row["position_x_mpc"], "position_y_mpc": row["position_y_mpc"],
            "position_z_mpc": row["position_z_mpc"],
            "center_distance_ly": _center_distance_ly(row),
            "runaway_class": row["runaway_class"], "runaway_speed_kms": row["runaway_speed_kms"],
            "inside": containing_cloud(conn, row),
            "nearest": nearest.get(row["id"], []),
            "stars": [dict(star_row) for star_row in star_rows],
        })
    star_count = sum(max(1, len(system["stars"])) for system in systems)

    return {
        "id": sector["id"], "name": sector["name"], "edge_mpc": sector["edge_mpc"],
        "edge_ly": milliparsecs_to_ly(sector["edge_mpc"]),
        "center_x_pc": sector["center_x_pc"], "center_y_pc": sector["center_y_pc"],
        "center_z_pc": sector["center_z_pc"], "galactic_radius_pc": sector["galactic_radius_pc"],
        "shell_index": sector["shell_index"], "band_index": sector["band_index"],
        "wedge_index": sector["wedge_index"],
        "placed": sector["center_x_pc"] is not None,
        "system_count": len(systems),
        "star_count": star_count,
        "interstellar_debris_count": interstellar_debris_count(star_count),
        "systems": systems,
        "phenomena": phenomena_near_sector(conn, sector_id),
        "neighbors": sector_neighbors(conn, sector),
        "stats": sector_stats(conn, sector["shell_index"], sector["band_index"], sector["wedge_index"]),
    }


SECTOR_STATS_FIELDS = (
    "bright_level_sol", "relative_density", "expected_systems", "actual_systems", "actual_stars",
    "mean_temperature_k", "mean_luminosity_sol", "mean_age_gy", "total_luminosity_sol", "habitable_worlds",
)
"""tuple: The `sector_stats` columns `sector_stats` reads (PERF.11)."""


def sector_stats(conn, shell_index, band_index, wedge_index):
    """
    One grid sector's stored stats (PERF.11, GEN.44): its bright-star
    level, its expected density and, once filled, what it holds
    (`SECTOR_STATS_FIELDS`). `None` for a sector off the grid or with no
    row yet.
    """
    if shell_index is None or band_index is None or wedge_index is None:
        return None
    row = conn.execute(
        f"SELECT {', '.join(SECTOR_STATS_FIELDS)} FROM sector_stats"
        " WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
        (shell_index, band_index, wedge_index),
    ).fetchone()
    return None if row is None else {field: row[field] for field in SECTOR_STATS_FIELDS}


def interstellar_debris_count(star_count):
    """
    About how many interstellar comets and planetesimals (meters to
    kilometers across) drift through a sector holding `star_count` stars:
    `INTERSTELLAR_DEBRIS_DENSITY_PC3` scales with stellar density like
    every other interstellar object, so it's a fixed number per star
    (~7e15). Far too many to be rows -- a figure for the sector page
    ("about 10^16 interstellar comets and planetesimals"); the
    `interstellar_comets` rows are only the notable ones.

    Returns:
        float: The estimated count (0 for an empty sector).
    """
    return (tuning.INTERSTELLAR_DEBRIS_DENSITY_PC3
            / tuning.REFERENCE_STELLAR_DENSITY_PC3 * star_count)


_PHENOMENON_TABLES = (
    ("nebulae", "nebula", "nebula_type", "radius_ly"),
    ("asteroid_fields", "asteroid_field", "density", "radius_ly"),
    # v21: black_holes/neutron_stars gained the same galaxy-frame placement
    # shape (schema.sql's "v21" header note) -- included here the same way
    # nebulae/asteroid_fields are, so a sector-generated compact remnant
    # shows up on the Sector Map/Galaxy Map too. Neither table has a
    # radius_ly column of its own (a compact remnant's real physical size
    # -- an event horizon a few km across, a neutron star ~10-20 km -- is
    # utterly negligible at sector/galaxy map scale), so it's queried as a
    # literal 0: a point object for the sphere-vs-cube overlap check below,
    # which is exactly correct (only shows up where its own center falls
    # within a sector's bounding sphere, never "spills into" a neighboring
    # sector's cube the way a nebula's real radius can).
    ("black_holes", "black_hole",
     "(CASE WHEN has_accretion_disk THEN 'accreting' ELSE 'quiescent' END)", "0"),
    ("neutron_stars", "neutron_star", "pulsar_type", "0"),
    # v28: the last three types gained the same placement columns
    # (schema.sql's "v28" header note). A supernova remnant has a real
    # radius_ly of its own; a rogue planet's radius_km and a comet's
    # nucleus are negligible at light-year scale, so both are points like
    # the compact remnants above. The comet is `interstellar_comet`, not
    # bare `comet`, to stay unambiguous next to the unrelated `comets`
    # table (a star system's own planet-orbiting comets).
    ("supernova_remnants", "supernova_remnant", "morphology", "radius_ly"),
    ("rogue_planets", "rogue_planet",
     "(CASE WHEN mass_bin = 'brown-dwarf' THEN 'brown dwarf' "
     "WHEN planet_type = 'g' THEN 'gas giant' ELSE 'terrestrial' END)", "0"),
    ("interstellar_comets", "interstellar_comet",
     "(CASE WHEN is_active THEN 'active' ELSE 'dormant' END)", "0"),
    # v31: a galaxy's active nucleus, always at the galactic center
    # (schema.sql's "v31" header note). Its jets can reach far past the
    # galaxy, but the engine itself is light-days across, so it's a point.
    ("quasars", "quasar",
     "(CASE WHEN is_radio_loud THEN 'radio-loud' ELSE 'radio-quiet' END)", "0"),
)
"""tuple: `(table_name, type_label, descriptor_expr, radius_expr)` for
every standalone phenomenon table -- all seven have galaxy-frame placement
columns since v28 (see `schema.sql`'s "v18"/"v21"/"v28" header notes),
though any single row may still be unplaced. `descriptor_expr` is a SQL
expression for each table's own one-line flavor field (a nebula's
`nebula_type`, a black hole's accretion state), normalized to a common
`descriptor` key so callers don't need to know which table a given `type`
came from; `radius_expr` is likewise a SQL expression (a plain column
where the object has a light-year-scale size, a literal `0` for the
point-like types)."""

_RADIUS_PHENOMENON_TABLES = tuple(
    table for table, _type_label, _descriptor_expr, radius_expr in _PHENOMENON_TABLES if radius_expr != "0"
)
"""tuple: The `_PHENOMENON_TABLES` tables with a real `radius_ly` column
(nebulae, asteroid fields, supernova remnants) --
`_widest_placed_phenomenon_radius_ly` checks just these."""


_PHENOMENON_MASS_COLUMNS = {"black_holes": "mass_solar", "neutron_stars": "mass_solar", "quasars": "black_hole_mass_solar"}
"""dict: The solar-mass column of each `_PHENOMENON_TABLES` table that has one (MAP.208: the maps
draw these as bright as a star of the mass); other types' rows carry `mass_solar` as `None`."""


_PHENOMENON_CLASS_COLUMNS = {
    "nebulae": "nebula_class",
    "supernova_remnants": "remnant_class",
    "asteroid_fields": "field_class",
    "rogue_planets": "planet_class",
}
"""dict: The letter-class column (v38; a rogue planet's `PLANET_CLASSES`
letter since v47) of each `_PHENOMENON_TABLES` table that has one; every
other type's rows carry `class` as `None`."""


def rogue_planet_gravity_g(mass_kg, radius_km):
    """A rogue planet's surface gravity in g (UX.99), from its mass and radius; `None` when either is missing."""
    if not mass_kg or not radius_km:
        return None
    return relations.surface_gravity(mass_kg, radius_km * 1000.0) / constants.EARTH_GRAVITY


def _placed_phenomenon_rows(conn, bbox=None, sector_id=None):
    """
    Reads every galaxy-placed row (non-NULL `center_x_pc`) from all four
    `_PHENOMENON_TABLES`, normalized to one
    common shape -- shared by `phenomena_near_sector` (which then filters
    by distance) and `galaxy_placed_phenomena` (which doesn't need to).

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        bbox (tuple, optional): `(center_x_pc, center_y_pc, center_z_pc,
            margin_pc)` -- when given, adds a `center_x_pc BETWEEN ...`
            (and y/z) SQL `WHERE` clause to each table's own query, the
            same bounding-box prefilter `galaxy_sectors_in_view` already
            uses for `sectors` (see that function's own docstring and
            `schema.sql`'s "v25"/"v26" header notes). `None` (the default,
            used by `galaxy_placed_phenomena`'s own whole-galaxy listing,
            which genuinely needs every row) runs the old unconditional
            scan. `phenomena_near_sector` is the only caller that passes
            this, since it's the only one that only ever needs a
            neighborhood, not the whole galaxy -- without it, this was a
            confirmed-in-production full-table-scan-across-four-tables on
            every single `sector_detail` call (see `schema.sql`'s "v26"
            header note), the same failure mode "v25"'s note already
            documented for `sectors` before that migration.
        sector_id (int, optional): When given (instead of `bbox`), reads
            only the placed rows whose own `sector_id` is this sector --
            an indexed lookup `phenomena_near_sector` uses so a sector
            always lists the phenomena generated as part of it.

    Returns:
        list[dict]: `id`, `type` (one of `_PHENOMENON_TABLES`' type
            labels), `name`, `descriptor`, `radius_ly` (0 for the
            point-like types), `class` (see `_PHENOMENON_CLASS_COLUMNS`), `gravity_g` (rogue planets only: the
            surface gravity in g, UX.99), `mass_solar` (a black hole, neutron star or quasar's mass, MAP.208;
            else `None`), `x`/`y`/`z` (`center_x/y/z_pc`), `galactic_radius_pc`, and for a supernova remnant `age_years` (MAP.286).
    """
    where_bbox = ""
    bbox_params = ()
    if sector_id is not None:
        where_bbox = " AND sector_id = ?"
        bbox_params = (sector_id,)
    elif bbox is not None:
        cx, cy, cz, margin_pc = bbox
        where_bbox = (
            " AND center_x_pc BETWEEN ? AND ?"
            " AND center_y_pc BETWEEN ? AND ?"
            " AND center_z_pc BETWEEN ? AND ?"
        )
        bbox_params = (
            cx - margin_pc, cx + margin_pc,
            cy - margin_pc, cy + margin_pc,
            cz - margin_pc, cz + margin_pc,
        )

    rows = []
    for table, type_label, descriptor_expr, radius_expr in _PHENOMENON_TABLES:
        query_rows = conn.execute(
            f"""
            SELECT id, name, {descriptor_expr} AS descriptor, {radius_expr} AS radius_ly,
                   {_PHENOMENON_CLASS_COLUMNS.get(table, "NULL")} AS class_code,
                   center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc,
                   {"mass_kg, radius_km" if table == "rogue_planets" else "NULL AS mass_kg, NULL AS radius_km"},
                   {_PHENOMENON_MASS_COLUMNS.get(table, "NULL")} AS mass_solar,
                   {"age_years" if table == "supernova_remnants" else "NULL"} AS age_years
            FROM {table}
            WHERE center_x_pc IS NOT NULL
            {where_bbox}
            """,
            bbox_params,
        ).fetchall()
        for row in query_rows:
            rows.append({
                "id": row["id"], "type": type_label, "name": row["name"],
                "descriptor": row["descriptor"], "radius_ly": row["radius_ly"],
                "class": row["class_code"],
                "x": row["center_x_pc"], "y": row["center_y_pc"], "z": row["center_z_pc"],
                "galactic_radius_pc": row["galactic_radius_pc"], "mass_solar": row["mass_solar"],
            })
            if table == "rogue_planets":
                rows[-1]["gravity_g"] = rogue_planet_gravity_g(row["mass_kg"], row["radius_km"])
            elif table == "supernova_remnants":
                rows[-1]["age_years"] = row["age_years"]
    return rows


def _widest_placed_phenomenon_radius_ly(conn):
    """
    The largest `radius_ly` among currently placed rows of the phenomenon
    tables with a real, non-zero radius (`_RADIUS_PHENOMENON_TABLES`), or
    `0.0` if none of them has any placed row at all.

    `phenomena_near_sector` uses this to size its own bounding-box margin
    (see that function's docstring): a bare `MAX(radius_ly)` aggregate
    query, not a full row fetch, so this stays cheap even on a large
    galaxy (MySQL can satisfy it from an index or a single scan of one
    narrow column, never the whole row for every placed phenomenon the
    way the old unconditional `_placed_phenomenon_rows` scan did) while
    keeping the bounding box exactly as correct as that old unconditional
    scan for a phenomenon of any size, rather than assuming a fixed
    ceiling that a future (or deliberately test-constructed) oversized row
    could silently fall outside of.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.

    Returns:
        float: The widest placed radius, light-years.
    """
    widest = 0.0
    for table in _RADIUS_PHENOMENON_TABLES:
        row = conn.execute(
            f"SELECT MAX(radius_ly) AS widest FROM {table} WHERE center_x_pc IS NOT NULL"
        ).fetchone()
        if row and row["widest"] is not None:
            widest = max(widest, row["widest"])
    return widest


PHENOMENON_SORTS = {
    "name": "name", "type": "type", "descriptor": "descriptor", "radius": "radius_ly",
    "sector": "sector_name", "placed": "placed",
}
"""dict: The sort keys `list_phenomena` accepts (the Phenomena table's column
keys) -> the column of the union they order by."""


_SCATTER_TYPES = {
    "black-hole": "black_hole", "neutron-star": "neutron_star", "planetary-nebula": "nebula",
    "supernova-remnant": "supernova_remnant", "hypervelocity-star": "hypervelocity_star",
}
"""dict: A `phenomenon_scatter` row's `kind` -> the Phenomena table's `type`."""


def _phenomena_union():
    """The `SELECT` that stacks every `_PHENOMENON_TABLES` table into one
    shape (`type`, `id`, `name`, `descriptor`, `radius_ly`, `sector_id`,
    `sector_name`, `placed`, `scattered`). A `black_holes`/`neutron_stars` row with
    `star_id` set is a normal system's own compact-remnant star, already on
    that system's page, so it is left out; every other table is always
    standalone. The scatter's unbuilt rows are not here: they number in the
    hundreds of millions, so `_scattered_classes` counts them from the
    stored class totals and `_scattered_page` reads just the page asked for."""
    return " UNION ALL ".join(
        f"""
        SELECT '{type_label}' AS type, t.id AS id, t.name AS name,
               {descriptor_expr} AS descriptor, {radius_expr} AS radius_ly,
               t.sector_id AS sector_id, sec.name AS sector_name,
               (t.center_x_pc IS NOT NULL) AS placed, 0 AS scattered
        FROM {table} t
        LEFT JOIN sectors sec ON sec.id = t.sector_id
        {"WHERE t.star_id IS NULL" if table in ("black_holes", "neutron_stars") else ""}
        """
        for table, type_label, descriptor_expr, radius_expr in _PHENOMENON_TABLES
    )


def _scattered_classes(conn, types=(), descriptors=(), placed=None, ignore=()):
    """
    The scattered classes that still have unbuilt phenomena and pass the
    Phenomena table's filters: `[{"kind", "subtype", "type", "descriptor",
    "unbuilt"}]`. A scattered phenomenon is placed on the galaxy, so
    `placed=False` matches none. `ignore` names a filter (`"types"` or
    `"descriptors"`) to leave out, for the filter menus' own counts.
    """
    if placed is False:
        return []
    classes = []
    for row in conn.execute(
            "SELECT kind, subtype, placed - built AS unbuilt FROM phenomenon_scatter_classes"
            " WHERE placed > built ORDER BY kind, subtype").fetchall():
        kind, subtype = row["kind"], row["subtype"]
        entry = {"kind": kind, "subtype": subtype, "type": _SCATTER_TYPES.get(kind, kind),
                 "descriptor": subtype or "scattered", "unbuilt": int(row["unbuilt"])}
        if types and "types" not in ignore and entry["type"] not in types:
            continue
        if descriptors and "descriptors" not in ignore and entry["descriptor"] not in descriptors:
            continue
        classes.append(entry)
    return classes


def parse_scatter_start(text):
    """The jump-by-value start of the Phenomena table (PERF.135) as a scatter id: a whole number of at least 1.

    Returns:
        int or None: The id, or `None` for an empty text.

    Raises:
        ValueError: If the text is not a whole number from 1 up.
    """
    text = (text or "").strip()
    if not text:
        return None
    try:
        number = int(text)
    except ValueError:
        raise ValueError(f"{text!r} is not a scatter id")
    if number < 1:
        raise ValueError(f"{text!r} is not a scatter id")
    return number


def _scattered_page(conn, classes, limit, offset, after=None, first=None, last=False):
    """`limit` unbuilt scatter rows from `offset` among the `classes`, in id order, shaped like the Phenomena rows
    (`scattered` True, named "Uncharted <kind> <shell>.<band>.<wedge>"). With `after` (a scatter row id) the page
    starts after that row and `offset` is not used (PERF.75). `first` (the jump-by-value box, PERF.135) keeps only the
    rows from that id on, `offset` counting from there; `last` takes the last `limit` rows instead."""
    if not classes:
        return []
    by_class = {(entry["kind"], entry["subtype"]): entry for entry in classes}
    clauses, params = [], []
    for entry in classes:
        clauses.append("(kind = ? AND COALESCE(subtype, '') = ?)")
        params.extend([entry["kind"], entry["subtype"]])
    later = ""
    if first is not None:
        later, params = "id >= ? AND ", [first] + params
    if after is not None:
        later, params, offset = "id > ? AND " + later, [after] + params, 0
    query = ("SELECT id, kind, COALESCE(subtype, '') AS subtype, shell_index, band_index, wedge_index,"
             f" member_count, system_mass_solar FROM phenomenon_scatter WHERE built_at IS NULL AND {later}({' OR '.join(clauses)}) "
             f"ORDER BY id{' DESC' if last else ''}")
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params.extend([limit, 0 if last else offset or 0])
    rows = conn.execute(query, params).fetchall()
    return [
        {"id": row["id"], "scatter_id": row["id"], "type": by_class[(row["kind"], row["subtype"])]["type"],
         "name": f"Uncharted {row['kind'].replace('-', ' ')} {row['shell_index']}.{row['band_index']}."
                 f"{row['wedge_index']}",
         "descriptor": by_class[(row["kind"], row["subtype"])]["descriptor"], "radius_ly": 0.0,
         "sector_id": None, "sector_name": None, "placed": True, "scattered": True, **_member_fields(row)}
        for row in (reversed(rows) if last else rows)
    ]


def _phenomena_where(types=(), descriptors=(), placed=None):
    """`(where_sql, params)` for the Phenomena table's filters: any of the
    `types`, any of the `descriptors` (an empty list filters nothing), and
    `placed` (True/False, or None for both)."""
    clauses, params = [], []
    for column, values in (("type", types), ("descriptor", descriptors)):
        if values:
            clauses.append(f"{column} IN ({', '.join('?' for _ in values)})")
            params.extend(values)
    if placed is not None:
        clauses.append("placed = ?")
        params.append(1 if placed else 0)
    return (" WHERE " + " AND ".join(clauses)) if clauses else "", params


def parse_built_cursor(text):
    """`(type, row id)` from the `<type>:<id>` text a page of built phenomena hands to the next (PERF.205), or `None`
    when `text` is a scatter id or is no cursor. The id is the row's database id, as `list_phenomena` returns it."""
    type_label, _colon, row_id = str(text).partition(":")
    if type_label in _PHENOMENON_TYPE_TO_TABLE and row_id.isdigit():
        return type_label, int(row_id)
    return None


SECTOR_WALK_MIN_ROWS = 50_000
"""int: Under this many built phenomena (the sum of the tables' highest ids) the Phenomena table's sector sort sorts them whole, which a
page of a few thousand rows costs next to nothing; from it up the sort walks the sectors by name (`_built_by_sector`, PERF.211)."""

SECTOR_CHUNK = 200
"""int: The sectors the sector sort reads at first for a page; each further read doubles it, to `SECTOR_CHUNK_MAX`."""
SECTOR_CHUNK_MAX = 5000

_COLLATION = " COLLATE utf8mb4_unicode_ci"


def _descriptor_key(descriptor_expr):
    """The expression of a table's descriptor as its descriptor index holds it (`idx_<table>_descriptor_name`): a CASE of text
    carries the collation of the columns, so that the sort of one table and the merge of all agree."""
    return descriptor_expr + _COLLATION if descriptor_expr.startswith("(") else descriptor_expr


_PLACED_KEY = "(t.center_x_pc IS NOT NULL)"
_BUILT_OUTER_KEYS = {"type": "type", "descriptor": "descriptor", "radius": "radius_ly", "placed": "placed"}
"""dict: The Phenomena sorts that order by a key before the name -> the column of the merge they order by."""


def _built_sort_key(sort, type_label, descriptor_expr, radius_expr):
    """`(sql, constant)` for how a phenomenon table orders by `sort`: the expression its index serves (`sql`), or the one value
    (`constant`) every row of the table has (the type, a radius of 0), or neither for the name."""
    if sort == "type":
        return None, type_label
    if sort == "placed":
        return _PLACED_KEY, None
    if sort == "descriptor":
        return _descriptor_key(descriptor_expr), None
    if sort == "radius":
        return (None, 0) if radius_expr == "0" else (f"t.{radius_expr}", None)
    return None, None


def _built_start(conn, sort, after):
    """What a page after the built row `after` (`(type, row id)`) continues from: the row's name, type, id, its key for `sort`
    and its sector, or `None` when the row is gone (the page then starts at the top)."""
    table = _PHENOMENON_TYPE_TO_TABLE[after[0]]
    entry = next(e for e in _PHENOMENON_TABLES if e[0] == table)
    key_sql, _constant = _built_sort_key(sort, entry[1], entry[2], entry[3])
    row = conn.execute(
        f"SELECT t.name AS name, t.sector_id AS sector_id, sec.name AS sector_name"
        f"{', ' + key_sql + ' AS sort_key' if key_sql else ''} FROM {table} t "
        f"LEFT JOIN sectors sec ON sec.id = t.sector_id WHERE t.id = ?", (after[1],)).fetchone()
    if row is None:
        return None
    return {"name": row["name"], "type": after[0], "id": after[1], "sector_id": row["sector_id"],
            "sector_name": row["sector_name"], "key": row["sort_key"] if key_sql else (entry[1] if sort == "type" else
                                                                                      0 if sort == "radius" else None)}


def _after_condition(key_sql, constant, type_label, start, descending):
    """`(sql, params)` for a table's rows that come after `start` in the order `key, name, type, id` (all in one direction), `None`
    when none of them do. Written as `_lexicographic` is, so the index serves the range; the type is not a column,
    so it decides whether a row equal to `start` in key and name comes after it (the type does, or the id does)."""
    ascending = not descending
    if constant is not None and key_sql is None and start["key"] is not None:
        if constant != start["key"]:
            return ("", []) if (constant > start["key"]) == ascending else None
    if (type_label > start["type"]) == ascending and type_label != start["type"]:
        closed = True
    elif type_label == start["type"]:
        closed = None
    else:
        closed = False
    columns = ([key_sql] if key_sql else []) + ["t.name"]
    values = ([start["key"]] if key_sql else []) + [start["name"]]
    if closed is None:
        columns.append("t.id")
        values.append(start["id"])
    return _lexicographic(columns, values, descending, inclusive=bool(closed))


def _built_by_key(conn, limit, sort, descending, types, descriptors, placed, start=None, extra=None):
    """
    A page of the built phenomena by `sort` read in index order (PERF.205, PERF.211): from each phenomenon table the `limit`
    rows that come next in `ORDER BY key, name, id`, found by walking the table's index on that key and name, and the best `limit`
    of those together in `ORDER BY key, name, type, id`. The page costs the same at any depth, since no table is read or sorted
    past the rows it returns, where the union of every table sorted whole for each page. Ties follow the direction of the sort.

    Args:
        limit (int): Rows to return.
        sort (str): `name`, `type`, `descriptor`, `radius` or `placed` (the key each table's index serves).
        descending (bool): Whether the keys, then the names, types and ids, sort downward.
        types, descriptors, placed: The Phenomena table's filters (`_phenomena_where`).
        start (dict, optional): The row the previous page ended on (`_built_start`).
        extra (str, optional): One more condition on a table `t` (`t.sector_id IS NULL`).

    Returns:
        list[dict]: The rows, as `list_phenomena` returns them.
    """
    direction = "DESC" if descending else "ASC"
    parts, params = [], []
    for table, type_label, descriptor_expr, radius_expr in _PHENOMENON_TABLES:
        if types and type_label not in types:
            continue
        key_sql, constant = _built_sort_key(sort, type_label, descriptor_expr, radius_expr)
        clauses, clause_params = [], []
        if table in ("black_holes", "neutron_stars"):
            clauses.append("t.star_id IS NULL")
        if extra:
            clauses.append(extra)
        if descriptors:
            clauses.append(f"{_descriptor_key(descriptor_expr)} IN ({', '.join('?' for _ in descriptors)})")
            clause_params.extend(descriptors)
        if placed is not None:
            clauses.append(f"{_PLACED_KEY} = ?")
            clause_params.append(1 if placed else 0)
        if start is not None:
            after = _after_condition(key_sql, constant, type_label, start, descending)
            if after is None:
                continue
            if after[0]:
                clauses.append(after[0])
                clause_params.extend(after[1])
        order = (f"{key_sql} {direction}, " if key_sql else "") + f"t.name {direction}, t.id {direction}"
        parts.append(
            f"(SELECT '{type_label}' AS type, t.id AS id, t.name AS name, {_descriptor_key(descriptor_expr)} AS descriptor, "
            f"{radius_expr} AS radius_ly, t.sector_id AS sector_id, sec.name AS sector_name, "
            f"{_PLACED_KEY} AS placed, 0 AS scattered FROM {table} t "
            f"LEFT JOIN sectors sec ON sec.id = t.sector_id "
            f"{'WHERE ' + ' AND '.join(clauses) + ' ' if clauses else ''}"
            f"ORDER BY {order} LIMIT ?)")
        params.extend(clause_params + [limit])
    if not parts:
        return []
    outer = f"{_BUILT_OUTER_KEYS[sort]} {direction}, " if sort in _BUILT_OUTER_KEYS else ""
    rows = conn.execute(
        f"SELECT * FROM ({' UNION ALL '.join(parts)}) AS phenomena "
        f"ORDER BY {outer}name {direction}, type {direction}, id {direction} LIMIT ?", params + [limit]).fetchall()
    return [_built_item(row) for row in rows]


def _built_row_bound(conn):
    """About how many built phenomena there are, from above: the sum of the highest ids of the phenomenon tables (index reads)."""
    reads = " + ".join(f"(SELECT COALESCE(MAX(id), 0) FROM {table})" for table, *_rest in _PHENOMENON_TABLES)
    return int(conn.execute(f"SELECT {reads} AS n").fetchone()["n"])


def _sector_rows(conn, limit, descending, types, descriptors, placed, sector_ids=None, start=None):
    """Up to `limit` built phenomena that are in a sector, in the order `sector name, sector id, name, type, id`: those of the
    `sector_ids`, or of every sector, and, with `start`, only those after that row."""
    direction = "DESC" if descending else "ASC"
    parts, params = [], []
    for table, type_label, descriptor_expr, radius_expr in _PHENOMENON_TABLES:
        if types and type_label not in types:
            continue
        clauses = ["t.sector_id IN (" + ", ".join("?" for _ in sector_ids) + ")"] if sector_ids is not None \
            else ["t.sector_id IS NOT NULL"]
        clause_params = list(sector_ids or [])
        if table in ("black_holes", "neutron_stars"):
            clauses.append("t.star_id IS NULL")
        if descriptors:
            clauses.append(f"{_descriptor_key(descriptor_expr)} IN ({', '.join('?' for _ in descriptors)})")
            clause_params.extend(descriptors)
        if placed is not None:
            clauses.append(f"{_PLACED_KEY} = ?")
            clause_params.append(1 if placed else 0)
        parts.append(
            f"SELECT '{type_label}' AS type, t.id AS id, t.name AS name, {_descriptor_key(descriptor_expr)} AS descriptor, "
            f"{radius_expr} AS radius_ly, t.sector_id AS sector_id, sec.name AS sector_name, "
            f"{_PLACED_KEY} AS placed, 0 AS scattered FROM {table} t "
            f"JOIN sectors sec ON sec.id = t.sector_id WHERE {' AND '.join(clauses)}")
        params.extend(clause_params)
    if not parts:
        return []
    after = ""
    if start is not None:
        condition, found = _lexicographic(["sector_name", "sector_id", "name", "type", "id"],
                                          [start["sector_name"], start["sector_id"], start["name"], start["type"], start["id"]], descending)
        after = f"WHERE {condition} "
        params.extend(found)
    rows = conn.execute(
        f"SELECT * FROM ({' UNION ALL '.join(parts)}) AS phenomena {after}"
        f"ORDER BY sector_name {direction}, sector_id {direction}, name {direction}, type {direction}, id {direction} LIMIT ?",
        params + [limit]).fetchall()
    return [_built_item(row) for row in rows]


def _built_by_sector(conn, limit, descending, types, descriptors, placed, after=None):
    """
    A page of the built phenomena by sector name (PERF.211), a phenomenon with no sector after those that have one (in either
    direction): the sectors are read by name in index order, a few hundred at a time (more each time while the page is not full), and
    only the phenomena of those are sorted, so a page costs about what it returns whatever its depth. No index could hold a name
    that lives in another table. With fewer than `SECTOR_WALK_MIN_ROWS` phenomena built they are sorted whole instead.
    Ties follow the direction of the sort, a sector's name then its id.

    Args:
        limit, descending, types, descriptors, placed: As `_built_by_key`.
        after (tuple, optional): `(type, row id)` of the row the previous page ended on.
    """
    start = None if after is None else _built_start(conn, "sector", after)
    rows = []
    if start is None or start["sector_id"] is not None:
        if _built_row_bound(conn) < SECTOR_WALK_MIN_ROWS:
            rows = _sector_rows(conn, limit, descending, types, descriptors, placed, start=start)
        else:
            direction = "DESC" if descending else "ASC"
            chunk, beyond = SECTOR_CHUNK, None
            while len(rows) < limit:
                if beyond is not None:
                    condition, where_params = _lexicographic(["name", "id"], list(beyond), descending)
                    where = f"WHERE {condition} "
                elif start is not None:
                    condition, where_params = _lexicographic(["name", "id"], [start["sector_name"], start["sector_id"]], descending,
                                                             inclusive=True)
                    where = f"WHERE {condition} "
                else:
                    where, where_params = "", []
                sectors = conn.execute(
                    f"SELECT id, name FROM sectors {where}ORDER BY name {direction}, id {direction} LIMIT ?",
                    where_params + [chunk]).fetchall()
                if not sectors:
                    break
                rows += _sector_rows(conn, limit - len(rows), descending, types, descriptors, placed,
                                     sector_ids=[sector["id"] for sector in sectors],
                                     start=start if beyond is None else None)
                if len(sectors) < chunk:
                    break
                beyond = (sectors[-1]["name"], sectors[-1]["id"])
                chunk = min(chunk * 2, SECTOR_CHUNK_MAX)
    if len(rows) < limit:
        rows += _built_by_key(conn, limit - len(rows), "name", descending, types, descriptors, placed,
                              start=start if start is not None and start["sector_id"] is None else None,
                              extra="t.sector_id IS NULL")
    return rows


def _built_page(conn, limit, sort, descending, types, descriptors, placed, after=None):
    """A page of the built phenomena by `sort`, from the top or after the `(type, row id)` in `after`."""
    if sort == "sector":
        return _built_by_sector(conn, limit, descending, types, descriptors, placed, after)
    start = None if after is None else _built_start(conn, sort, after)
    return _built_by_key(conn, limit, sort, descending, types, descriptors, placed, start)


def _built_item(row):
    """One built phenomenon of the listing, from a row of the union."""
    return {
        "id": row["id"], "type": row["type"], "name": row["name"],
        "descriptor": row["descriptor"], "radius_ly": row["radius_ly"],
        "sector_id": row["sector_id"], "sector_name": row["sector_name"],
        "placed": bool(row["placed"]), "scattered": bool(row["scattered"]),
    }


def list_phenomena(conn, limit=None, offset=None, sort="name", descending=False, types=(), descriptors=(),
                   placed=None, after=None, start=None, last=False):
    """
    Returns every exotic phenomenon (nebula/asteroid field/black hole/
    neutron star/supernova remnant/rogue planet/interstellar comet/quasar --
    every table in `_PHENOMENON_TABLES`), across every sector and
    regardless of galaxy placement -- `GET /api/phenomena`'s own flat
    listing (`html/phenomena.py`), unlike `galaxy_placed_phenomena` (which
    only returns the galaxy-placed subset, for the Galaxy Map) or
    `phenomena_near_sector` (one sector's own neighborhood).

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        limit (int, optional): Caps the number of rows returned.
        offset (int, optional): Skips this many rows first. Ignored unless
            `limit` is also given.
        sort (str): A key of `PHENOMENON_SORTS`; ties fall back to name,
            type and id (in the direction of the sort) so pages never overlap. The first page and the pages after a
            row (`after`) of any sort are read by key; a page by `offset` past the first sorts the whole union.
        descending (bool): Reverse the sort.
        types (iterable[str]): Keep only these types (any of).
        descriptors (iterable[str]): Keep only these descriptors (any of).
        placed (bool, optional): Keep only placed or only unplaced ones.
        after (int or tuple, optional): A scattered row's `scatter_id`: the page starts after it in the tail of unbuilt
            scattered phenomena (they follow every built one, in id order), and the built ones and `offset` are
            not read (a keyset page, PERF.75). Or the `(type, row id)` of a built row (`parse_built_cursor`): the page
            starts after it in the order of `sort`, then the tail (PERF.205, PERF.211: `offset` is not used, and every sort
            reads its built rows in index order, so a page costs the same at any depth). Without `limit` it is not used.
        start (str, optional): The jump-by-value box (PERF.135), a scatter id: the list is the tail of unbuilt
            scattered phenomena from that id on (the built ones are not listed), `offset` counting from there.
            Needs `limit`; text that is not an id raises `ValueError`.
        last (bool): With `start`, the last `limit` scattered rows instead (where the id lies past the end).

    Excludes a `black_holes`/`neutron_stars` row with `star_id` set -- that
    shape is a normal star system's own compact-remnant star (already
    shown on that system's own `system.py` page), not a standalone exotic
    phenomenon; every other table here has no `star_id` at all (always
    standalone, see their own table comments) and needs no such filter.

    Returns:
        list[dict]: One row per phenomenon, ordered by `sort`: `id`, `type`
            (`"nebula"`, `"asteroid_field"`, `"black_hole"`,
            `"neutron_star"`, `"supernova_remnant"`, `"rogue_planet"`,
            `"interstellar_comet"` or `"quasar"`), `name`, `descriptor`,
            `radius_ly`, `sector_id`/`sector_name` (both `None` if this
            phenomenon has never been linked to a sector -- see
            `schema.sql`'s "v18" header note), `placed` (bool --
            whether it has a galaxy position at all) and `scattered` (bool
            -- placed by the scatter but not yet built, so `id` is its
            `phenomenon_scatter` row, with a page of its own; its `type` can also be
            `"hypervelocity_star"`). A scattered row also has `member_count` and `system_mass_solar`
            (`_member_fields`, UX.114).
    """
    if sort not in PHENOMENON_SORTS:
        raise ValueError(f"unknown phenomenon sort {sort!r}")
    first = parse_scatter_start(start)
    if first is not None and limit is not None:
        return _scattered_page(conn, _scattered_classes(conn, types, descriptors, placed), limit, offset or 0,
                               after=after, first=None if last else first, last=last)
    if isinstance(after, tuple) and limit is not None:
        return _built_then_scattered(conn, limit, sort, descending, types, descriptors, placed, after)
    if after is not None and limit is not None:
        return _scattered_page(conn, _scattered_classes(conn, types, descriptors, placed), limit, 0, after=after)
    begin = offset or 0
    if begin == 0 and limit is not None:
        return _built_then_scattered(conn, limit, sort, descending, types, descriptors, placed, None)
    where, filter_params = _phenomena_where(types, descriptors, placed)
    column = PHENOMENON_SORTS[sort]
    direction = "DESC" if descending else "ASC"
    tie = "sector_id " + direction + ", " if sort == "sector" else ""
    query = (f"SELECT * FROM ({_phenomena_union()}) AS phenomena{where} "
             f"ORDER BY {column} IS NULL, {column} {direction}, {tie}name {direction}, type {direction}, id {direction}")
    params = list(filter_params)
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params.extend([limit, begin])
    items = [_built_item(row) for row in conn.execute(query, params).fetchall()]
    # The unbuilt scattered phenomena come after the built ones, in id order. The built ones are counted only when the
    # page starts past some of them, to know where in the tail it starts (PERF.205: not for the first page).
    room = None if limit is None else limit - len(items)
    if room is None or room > 0:
        into_tail = 0
        if begin:
            built_total = conn.execute(
                f"SELECT COUNT(*) AS n FROM ({_phenomena_union()}) AS phenomena{where}", filter_params).fetchone()["n"]
            into_tail = max(0, begin - built_total)
        items += _scattered_page(conn, _scattered_classes(conn, types, descriptors, placed), room, into_tail)
    return items


def _built_then_scattered(conn, limit, sort, descending, types, descriptors, placed, after):
    """A page by `sort` of the built phenomena after `after` (`(type, row id)`, or the top), and when the built ones end
    inside it the start of the scattered tail that follows them (PERF.205, PERF.211)."""
    items = _built_page(conn, limit, sort, descending, types, descriptors, placed, after)
    room = limit - len(items)
    if room > 0:
        items += _scattered_page(conn, _scattered_classes(conn, types, descriptors, placed), room, 0)
    return items


def _phenomena_stamp(conn):
    """
    The stamp the stored Phenomena counts are kept under (UX.110): the galaxy's content stamp plus what the
    phenomena themselves hold. The content stamp moves with sectors and systems only, so a scatter (the stored class
    totals), a sector fill that builds scattered phenomena, or a standalone phenomenon saved or edited left the
    stored count at its old figure for hours (`countcache.MAX_AGE_SECONDS`). One primary-key and one `modified_at`
    index read per phenomenon table, and the class totals.
    """
    reads = ", ".join(
        f"(SELECT CONCAT(COALESCE(MAX(id), 0), '/', COALESCE(MAX(modified_at), '')) FROM {table}) AS t{n}"
        for n, (table, *_rest) in enumerate(_PHENOMENON_TABLES))
    row = conn.execute(f"SELECT {reads}, (SELECT COALESCE(SUM(placed), 0) FROM phenomenon_scatter_classes) AS placed, "
                       "(SELECT COALESCE(SUM(built), 0) FROM phenomenon_scatter_classes) AS built").fetchone()
    token = "|".join(str(row[key]) for key in row.keys())
    return galaxy_content_stamp(conn) + hashlib.sha256(token.encode("utf-8")).hexdigest()[:16]


def count_phenomena(conn, types=(), descriptors=(), placed=None, start=None):
    """
    Returns the number of exotic phenomena across every type in
    `_PHENOMENON_TABLES` that pass the same filters as `list_phenomena`,
    ignoring any pagination -- the denominator `list_phenomena(conn,
    limit=...)` callers (the API's `/api/phenomena`) need to report how many
    pages exist.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        start (str, optional): The jump-by-value box's scatter id (PERF.135): only the unbuilt scattered phenomena
            from that id on are counted, as `list_phenomena` lists them.

    Returns:
        int: Phenomenon count. Kept by `countcache` (PERF.69): the last count made when the galaxy has changed
            since, a capped count (`CappedCount`) before the first.
    """
    first = parse_scatter_start(start)
    if first is not None:
        def exact_tail(c):
            classes = _scattered_classes(c, types, descriptors, placed)
            if not classes:
                return 0
            clauses = " OR ".join("(kind = ? AND COALESCE(subtype, '') = ?)" for _ in classes)
            marks = [value for entry in classes for value in (entry["kind"], entry["subtype"])]
            return c.execute(f"SELECT COUNT(*) AS n FROM phenomenon_scatter WHERE built_at IS NULL AND id >= ? AND ({clauses})",
                             [first] + marks).fetchone()["n"]

        return _counted(conn, ["count_phenomena", list(types), list(descriptors), placed, first], exact_tail,
                        lambda c: _capped(min(exact_tail(c), COUNT_FALLBACK_CAP)), _phenomena_stamp)
    where, params = _phenomena_where(types, descriptors, placed)

    def exact(c):
        built = c.execute(
            f"SELECT COUNT(*) AS n FROM ({_phenomena_union()}) AS phenomena{where}", params).fetchone()["n"]
        return built + sum(entry["unbuilt"] for entry in _scattered_classes(c, types, descriptors, placed))

    def estimate(c):
        built = c.execute(f"SELECT COUNT(*) AS n FROM (SELECT 1 FROM ({_phenomena_union()}) AS phenomena{where} "
                          f"LIMIT ?) c", list(params) + [COUNT_FALLBACK_CAP]).fetchone()["n"]
        return _capped(built + sum(entry["unbuilt"] for entry in _scattered_classes(c, types, descriptors, placed)))

    return _counted(conn, ["count_phenomena", list(types), list(descriptors), placed], exact, estimate,
                    _phenomena_stamp)


def phenomena_facets(conn, types=(), descriptors=(), placed=None):
    """
    The option counts the Phenomena table's filter menus show: how many
    phenomena each `type` has and each `descriptor` (a class, kind or state:
    nebula class, remnant morphology, rogue planet size ...). Each count
    applies every filter except its own menu's, so picking a type narrows the
    descriptor options and the other way round.

    Returns:
        dict: `{"type": [{"value", "count"}], "descriptor": [{"value",
            "count"}]}`, options ordered by value.
    """
    def exact(c):
        facets = {}
        for column, own in (("type", "types"), ("descriptor", "descriptors")):
            filters = {"types": types, "descriptors": descriptors, own: ()}
            where, params = _phenomena_where(filters["types"], filters["descriptors"], placed)
            rows = c.execute(
                f"SELECT {column} AS value, COUNT(*) AS n FROM ({_phenomena_union()}) AS phenomena{where} "
                f"GROUP BY {column} ORDER BY {column}", params).fetchall()
            merged = {row["value"]: row["n"] for row in rows}
            for entry in _scattered_classes(c, types, descriptors, placed, ignore=(own,)):
                merged[entry[column]] = merged.get(entry[column], 0) + entry["unbuilt"]
            facets[column] = [{"value": value, "count": merged[value]} for value in sorted(merged)]
        return facets

    return _counted(conn, ["phenomena_facets", list(types), list(descriptors), placed], exact,
                    lambda c: {"type": [], "descriptor": []}, _phenomena_stamp)


_PHENOMENON_TYPE_TO_TABLE = {
    type_label: table for table, type_label, _de, _re in _PHENOMENON_TABLES
}
"""dict: `type` value (as returned by `list_phenomena`/`galaxy_placed_phenomena`)
-> its backing table name, e.g. `"nebula"` -> `"nebulae"` -- the reverse of
`_PHENOMENON_TABLES`'s own `(table, type_label, ...)` order, used by
`phenomenon_detail` to find the one table a `(type, id)` pair actually
means without hand-listing the mapping a second time."""


OBJECT_SIBLING_LIMIT = 200
"""int: The most sibling references `resolve_object` returns."""


def _sibling_refs(conn, kind, table, where, params, object_id):
    """The references of the other rows of `table` matching `where` (the
    object's parent), by id, up to `OBJECT_SIBLING_LIMIT`; `(refs, truncated)`."""
    rows = conn.execute(
        f"SELECT id FROM {table} WHERE {where} AND id <> ? ORDER BY id LIMIT {OBJECT_SIBLING_LIMIT + 1}",
        (*params, object_id),
    ).fetchall()
    refs = [object_ref.format(kind, row["id"]) for row in rows[:OBJECT_SIBLING_LIMIT]]
    return refs, len(rows) > OBJECT_SIBLING_LIMIT


def _enclosed_mass_kg(conn, radius_pc):
    """The mass, kg, the stored galaxy holds inside `radius_pc` of its centre (`galaxy.enclosed_mass`, GEN.370): from
    its stored shape and seed, or the shape the fallback seed draws when no skeleton is stored."""
    skeleton = get_galaxy_shape(conn)
    if skeleton is None:
        return enclosed_mass.enclosed_mass_kg(
            shape_draw.default_shape(None), mass_chain.FALLBACK_SEED, radius_pc)
    return enclosed_mass.enclosed_mass_kg(skeleton.shape, skeleton.galaxy_seed, radius_pc)


def keep_out_radius(conn, kind, object_id):
    """
    The keep-out radius of one object (NAV.24, `planetgen.galaxy.keepout`).

    Returns:
        keepout.KeepOut: `radius_km`, `basis` and `note`.

    Raises:
        ValueError: For a kind with no such rule (a sector) or a missing row.
    """
    def row(table, columns, row_id=object_id):
        found = conn.execute(f"SELECT {columns} FROM {table} WHERE id = ?", (row_id,)).fetchone()
        if found is None:
            raise ValueError(f"no {table} row with id {row_id}")
        return found

    def system_perimeter(system_id):
        system = row("star_systems", "binary_system_perimeter_km", system_id)
        stars = conn.execute(
            "SELECT system_perimeter_km FROM stars WHERE star_system_id = ?", (system_id,)).fetchall()
        sizes = [s["system_perimeter_km"] for s in stars] + [system["binary_system_perimeter_km"]]
        sizes = [v for v in sizes if v is not None]
        return keepout.KeepOut(max(sizes), "perimeter", None) if sizes else keepout.NO_KEEP_OUT

    if kind in ("planet", "moon"):
        body = row(object_ref.TABLES[kind], "hill_radius_km, radius_km")
        if body["hill_radius_km"] is None:
            return keepout.KeepOut(body["radius_km"], "radius", None)
        return keepout.KeepOut(body["hill_radius_km"], "hill", None)
    if kind == "star":
        return system_perimeter(row("stars", "star_system_id")["star_system_id"])
    if kind == "system":
        row("star_systems", "id")
        return system_perimeter(object_id)
    if kind in ("black_hole", "neutron_star"):
        table = _PHENOMENON_TYPE_TO_TABLE[kind]
        body = row(table, "star_id, mass_solar, radius_km, galactic_radius_pc"
                   if kind == "neutron_star" else "star_id, mass_solar, event_horizon_radius_km AS radius_km, "
                   "galactic_radius_pc")
        if body["galactic_radius_pc"] is None and body["star_id"] is not None:
            # An anchored remnant is a star of its own system.
            return system_perimeter(row("stars", "star_system_id", body["star_id"])["star_system_id"])
        if body["galactic_radius_pc"] is None:
            return keepout.NO_KEEP_OUT  # never placed in the galaxy: nowhere to measure from
        return keepout.compact_keep_out(
            _msun_to_kg(body["mass_solar"]), body["galactic_radius_pc"], body["radius_km"],
            _enclosed_mass_kg(conn, body["galactic_radius_pc"]))
    if kind in ("quasar", "rogue_planet"):
        table = _PHENOMENON_TYPE_TO_TABLE[kind]
        body = row(table, "black_hole_mass_solar * " + str(constants.SOLAR_MASS_TO_KG)
                   + " AS mass_kg, event_horizon_radius_km AS radius_km, galactic_radius_pc"
                   if kind == "quasar" else "mass_kg, radius_km, galactic_radius_pc")
        if body["galactic_radius_pc"] is None:
            return keepout.NO_KEEP_OUT
        return keepout.compact_keep_out(
            body["mass_kg"], body["galactic_radius_pc"], body["radius_km"],
            _enclosed_mass_kg(conn, body["galactic_radius_pc"]))
    if kind in ("nebula", "supernova_remnant", "asteroid_field"):
        row(_PHENOMENON_TYPE_TO_TABLE[kind], "id")
        return keepout.KeepOut(None, "none", keepout.PASS_THROUGH_NOTE)
    if kind in ("belt", "comet", "interstellar_comet", "facility"):
        row(object_ref.TABLES.get(kind) or _PHENOMENON_TYPE_TO_TABLE[kind], "id")
        return keepout.NO_KEEP_OUT
    if kind in object_ref.CONTAINER_KINDS:     # a place in the galaxy, not a body: nothing to keep out of
        row("galaxy_containers", "id")
        return keepout.NO_KEEP_OUT
    raise ValueError(f"no keep-out radius for a {kind}")


def keep_out_radius_dict(conn, kind, object_id):
    """`keep_out_radius` as the plain dict `resolve_object` returns; a
    sector, which has none, gets the empty one."""
    found = keepout.NO_KEEP_OUT if kind == "sector" else keep_out_radius(conn, kind, object_id)
    return found._asdict()


def resolve_object(conn, kind, object_id):
    """
    Resolves one object reference (NAV.7, `planetgen.galaxy.objectref`).

    Returns:
        dict: `ref`, `kind`, `id`, `name`; `parents` (the chain from the
            galaxy down to the direct parent, each `{ref, kind, name}`);
            `siblings` (the references of the other objects of the same
            kind under the same parent, by id, at most
            `OBJECT_SIBLING_LIMIT`) and `siblings_truncated`; and
            `positions`: `galaxy_pc` (galaxy frame, parsecs), `sector_ly`
            (sector-local, light-years) and `system_km` (system-local,
            kilometers from the system's origin, the barycenter of its
            stars), each `[x, y, z]` or `None` where the object has no
            position in that frame. A body inside a system takes the
            system's galaxy and sector positions (its offset is far below
            the precision of either); an asteroid belt is a ring and has no
            `system_km`.

    Raises:
        ValueError: For an unknown kind or a missing row.
    """
    if kind not in object_ref.KINDS:
        raise ValueError(f"unknown object kind: {kind!r}")
    galaxy = {"ref": object_ref.GALAXY, "kind": "galaxy", "name": "Galaxy"}
    positions = {"galaxy_pc": None, "sector_ly": None, "system_km": None}

    def one(table, row_id, columns):
        row = conn.execute(f"SELECT {columns} FROM {table} WHERE id = ?", (row_id,)).fetchone()
        if row is None:
            raise ValueError(f"no {table} row with id {row_id}")
        return row

    def system_parents_and_positions(system_id):
        """Parents (galaxy, sector, system) and the system's positions."""
        row = conn.execute(
            """
            SELECT ss.name, ss.sector_id, sec.name AS sector_name,
                   ss.position_x_mpc, ss.position_y_mpc, ss.position_z_mpc,
                   sec.center_x_pc, sec.center_y_pc, sec.center_z_pc
            FROM star_systems ss LEFT JOIN sectors sec ON sec.id = ss.sector_id
            WHERE ss.id = ?
            """, (system_id,)).fetchone()
        if row is None:
            raise ValueError(f"no star_systems row with id {system_id}")
        chain = [galaxy]
        if row["sector_id"] is not None:
            chain.append({"ref": object_ref.format("sector", row["sector_id"]), "kind": "sector",
                          "name": row["sector_name"]})
        chain.append({"ref": object_ref.format("system", system_id), "kind": "system", "name": row["name"]})
        pos = {"galaxy_pc": None, "sector_ly": None, "system_km": None}
        if row["position_x_mpc"] is not None:
            mpc = (row["position_x_mpc"], row["position_y_mpc"], row["position_z_mpc"])
            pos["sector_ly"] = [milliparsecs_to_ly(v) for v in mpc]
            if row["center_x_pc"] is not None:
                center = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
                pos["galaxy_pc"] = list(local_to_galaxy_pc(center, tuple(mpc_to_pc(v) for v in mpc)))   # the sector frame is oriented (TEST.453)
        return chain, pos

    def star_offset(star_id):
        """A star's system-local offset (km) from its system's origin."""
        if star_id is None:
            return [0.0, 0.0, 0.0]
        row = conn.execute(
            """
            SELECT s.role, ss.binary_primary_position_x_km AS px, ss.binary_primary_position_y_km AS py,
                   ss.binary_primary_position_z_km AS pz, ss.binary_secondary_position_x_km AS sx,
                   ss.binary_secondary_position_y_km AS sy, ss.binary_secondary_position_z_km AS sz,
                   s.offset_x_au, s.offset_y_au, s.offset_z_au
            FROM stars s JOIN star_systems ss ON ss.id = s.star_system_id WHERE s.id = ?
            """, (star_id,)).fetchone()
        if row is None:
            return [0.0, 0.0, 0.0]
        if row["role"] == "member":    # GEN.630: a star beyond the pair has no column of the system's; its stored offset
            return [_au_to_km(row[f"offset_{axis}_au"] or 0.0) for axis in "xyz"]
        values = (row["sx"], row["sy"], row["sz"]) if row["role"] == "secondary" else (row["px"], row["py"], row["pz"])
        return [v or 0.0 for v in values]

    def add(a, b):
        return [x + y for x, y in zip(a, b)]

    def offset(row):
        return [row["position_x_km"], row["position_y_km"], row["position_z_km"]]

    if kind == "sector":
        row = one("sectors", object_id, "name, center_x_pc, center_y_pc, center_z_pc")
        name, parents = row["name"], [galaxy]
        if row["center_x_pc"] is not None:
            positions["galaxy_pc"] = [row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]]
        positions["sector_ly"] = [0.0, 0.0, 0.0]
        siblings, truncated = _sibling_refs(conn, kind, "sectors", "1 = 1", (), object_id)
    elif kind == "system":
        sector_id = one("star_systems", object_id, "sector_id")["sector_id"]
        parents, positions = system_parents_and_positions(object_id)
        name = parents.pop()["name"]
        positions["system_km"] = [0.0, 0.0, 0.0]
        siblings, truncated = _sibling_refs(conn, kind, "star_systems", "sector_id <=> ?", (sector_id,), object_id)
    elif kind in object_ref.PHENOMENON_KINDS:
        row = one(_PHENOMENON_TYPE_TO_TABLE[kind], object_id, "name, center_x_pc, center_y_pc, center_z_pc")
        name, parents = row["name"], [galaxy]
        if row["center_x_pc"] is not None:
            positions["galaxy_pc"] = [row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]]
        siblings, truncated = [], False
    elif kind in object_ref.CONTAINER_KINDS:    # UX.184
        row = one("galaxy_containers", object_id, "name, kind, center_x_pc, center_y_pc, center_z_pc")
        if row["kind"] != kind:
            raise ValueError(f"no {kind} with id {object_id}")
        name, parents = row["name"], [galaxy]
        positions["galaxy_pc"] = [row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]]
        siblings, truncated = _sibling_refs(conn, kind, "galaxy_containers", "kind = ?", (kind,), object_id)
    elif kind == "facility":    # UX.183
        row = one("facilities", object_id,
                  "name, host_type, star_system_id, star_id, planet_id, moon_id, asteroid_belt_id, asteroid_field_id, "
                  "sector_id, center_x_pc, center_y_pc, center_z_pc")
        name = row["name"]
        host_column = {"star": "star_id", "planet": "planet_id", "moon": "moon_id", "asteroid_belt": "asteroid_belt_id",
                       "asteroid_field": "asteroid_field_id", "space": "sector_id"}[row["host_type"]]
        if row["star_system_id"] is not None:
            parents, positions = system_parents_and_positions(row["star_system_id"])
            host_kind, host_table = {"star": ("star", "stars"), "planet": ("planet", "planets"), "moon": ("moon", "moons"),
                                     "asteroid_belt": ("belt", "asteroid_belts")}.get(row["host_type"], (None, None))
            host_id = row[host_column]
            if host_kind == "moon":     # its planet is above it in the chain
                planet = one("moons", host_id, "planet_id")["planet_id"]
                if planet is not None:
                    parents.append({"ref": object_ref.format("planet", planet), "kind": "planet",
                                    "name": one("planets", planet, "name")["name"]})
            if host_kind is not None:
                host_name = (f"Asteroid belt {one(host_table, host_id, 'orbital_index')['orbital_index'] + 1}"
                             if host_kind == "belt" else one(host_table, host_id, "name")["name"])
                parents.append({"ref": object_ref.format(host_kind, host_id), "kind": host_kind, "name": host_name})
        else:
            parents = [galaxy]
            field_sector = None
            if row["asteroid_field_id"] is not None:
                field = one("asteroid_fields", row["asteroid_field_id"], "name, sector_id")
                field_sector = field["sector_id"]
            sector_id = row["sector_id"] if row["sector_id"] is not None else field_sector
            if sector_id is not None:
                parents.append({"ref": object_ref.format("sector", sector_id), "kind": "sector",
                                "name": one("sectors", sector_id, "name")["name"]})
            if row["asteroid_field_id"] is not None:
                parents.append({"ref": object_ref.format("asteroid_field", row["asteroid_field_id"]),
                                "kind": "asteroid_field", "name": field["name"]})
            if row["center_x_pc"] is not None:
                positions["galaxy_pc"] = [row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]]
        if row["star_system_id"] is not None:
            siblings, truncated = _sibling_refs(conn, kind, "facilities", f"{host_column} = ?", (row[host_column],), object_id)
        else:
            siblings, truncated = _sibling_refs(conn, kind, "facilities", "star_system_id IS NULL AND sector_id <=> ?",
                                                (row["sector_id"],), object_id)
    else:
        table = object_ref.TABLES[kind]
        columns = {
            "star": "name, star_system_id",
            "planet": f"name, star_system_id, star_id, position_x_km, position_y_km, position_z_km, {_tiers('')}",
            "moon": "name, star_system_id, planet_id, rogue_planet_id, star_id, position_x_km, position_y_km, "
                    f"position_z_km, {_tiers('')}",
            "belt": "star_system_id, orbital_index",
            "comet": "name, star_system_id, star_id, position_x_km, position_y_km, position_z_km",
        }[kind]
        row = one(table, object_id, columns)
        free_moon = kind == "moon" and row["rogue_planet_id"] is not None
        if free_moon:    # DB.94: a moon of a free planet has no system; its parent is the rogue planet
            rogue = one("rogue_planets", row["rogue_planet_id"], "name, center_x_pc, center_y_pc, center_z_pc")
            parents = [galaxy, {"ref": object_ref.format("rogue_planet", row["rogue_planet_id"]), "kind": "rogue_planet",
                                "name": rogue["name"]}]
            if rogue["center_x_pc"] is not None:
                positions["galaxy_pc"] = [rogue["center_x_pc"], rogue["center_y_pc"], rogue["center_z_pc"]]
        else:
            parents, positions = system_parents_and_positions(row["star_system_id"])
        name = f"Asteroid belt {row['orbital_index'] + 1}" if kind == "belt" else row["name"]
        if kind == "star":
            positions["system_km"] = star_offset(object_id)
        elif kind in ("planet", "comet"):
            positions["system_km"] = add(star_offset(row["star_id"]), offset(row))
        elif kind == "moon" and not free_moon:
            planet = one("planets", row["planet_id"], "name, position_x_km, position_y_km, position_z_km")
            parents.append({"ref": object_ref.format("planet", row["planet_id"]), "kind": "planet",
                            "name": planet["name"]})
            positions["system_km"] = add(add(star_offset(row["star_id"]), offset(planet)), offset(row))
        if free_moon:
            siblings, truncated = _sibling_refs(conn, kind, table, "rogue_planet_id = ?", (row["rogue_planet_id"],),
                                                object_id)
        elif kind == "moon":
            siblings, truncated = _sibling_refs(conn, kind, table, "planet_id = ?", (row["planet_id"],), object_id)
        else:
            siblings, truncated = _sibling_refs(
                conn, kind, table, "star_system_id = ?", (row["star_system_id"],), object_id)
    found = {
        "ref": object_ref.format(kind, object_id), "kind": kind, "id": object_id, "name": name,
        "parents": parents, "siblings": siblings, "siblings_truncated": truncated,
        "positions": positions,
        "keep_out": keep_out_radius_dict(conn, kind, object_id),
    }
    if kind in ("planet", "moon"):
        found.update({column: row[column] for column in _tier_columns()})
    return found


def nebula_shape(conn, nebula_id):
    """
    One nebula's shape (GEN.75, `planetgen.galaxy.nebula_shape`): the one
    stored with it (`nebulae.shape_*` and `nebula_shape_balls`), or, for a
    nebula saved before v56, the one a fresh save would have stored,
    drawn from its own properties.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        nebula_id (int): `nebulae.id`.

    Returns:
        tuple: `(row, shape)`: the nebula's row (a dict) and its `NebulaShape`.

    Raises:
        ValueError: If no such nebula exists.
    """
    from planetgen.galaxy import nebula_shape as shapes

    row = conn.execute("SELECT * FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
    if row is None:
        raise ValueError(f"no such nebula: {nebula_id!r}")
    row = _with_printed_uid(row)
    if row.get("shape_scale") is None:
        return row, shapes.shape_for_nebula(row["nebula_class"], row["radius_ly"], row["density_cm3"],
                                            row["temperature_k"], row["extinction_av"], row["dominant_species"])
    balls = conn.execute(
        "SELECT ball_index, center_x, center_y, center_z, radius FROM nebula_shape_balls"
        " WHERE nebula_id = ? ORDER BY ball_index", (nebula_id,)).fetchall()
    return row, shapes.from_columns(row, [tuple(b[k] for k in ("ball_index", "center_x", "center_y", "center_z", "radius"))
                                          for b in balls])


NEBULA_SURROUNDINGS_STARS = 300
"""int: Most bright stars `nebula_surroundings` lists around a nebula."""

NEBULA_SURROUNDINGS_REACH = 2.5
"""float: How far round a nebula `nebula_surroundings` looks for stars, in
nebula radii (the half-width of its box)."""

NEBULA_SURROUNDINGS_MIN_PC = 30.0
"""float: The least half-width of that box, parsecs, so a small nebula still
shows some stars round it."""


def nebula_surroundings(conn, nebula_id):
    """
    The brightest stars round one nebula, for its page's 3D view (MAP.105):
    the pre-placed bright stars (`galaxy_bright_stars_in_box`) in a box
    `NEBULA_SURROUNDINGS_REACH` radii each side of its centre, at most
    `NEBULA_SURROUNDINGS_STARS`.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        nebula_id (int): `nebulae.id`.

    Returns:
        dict: `radius_pc`, `half_width_pc` (of the box) and `stars`: the
            most luminous first, each `x`/`y`/`z` (parsecs from the
            nebula's centre), `luminosity_sol` and `temperature_k`. No
            stars for a nebula never placed in the galaxy.

    Raises:
        ValueError: If no such nebula exists.
    """
    row = conn.execute(
        "SELECT center_x_pc, center_y_pc, center_z_pc, radius_ly FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
    if row is None:
        raise ValueError(f"no such nebula: {nebula_id!r}")
    radius_pc = ly_to_pc(row["radius_ly"])
    half = max(NEBULA_SURROUNDINGS_MIN_PC, NEBULA_SURROUNDINGS_REACH * radius_pc)
    result = {"radius_pc": radius_pc, "half_width_pc": half, "stars": []}
    if row["center_x_pc"] is None:
        return result
    center = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
    skeleton = get_galaxy_shape(conn)
    edge_pc = skeleton.edge_pc if skeleton is not None else ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    stars = galaxy_bright_stars_in_box(
        conn, tuple(c - half for c in center), tuple(c + half for c in center), edge_pc,
        limit=NEBULA_SURROUNDINGS_STARS)
    result["stars"] = [
        {"x": round(star["x"] - center[0], 3), "y": round(star["y"] - center[1], 3),
         "z": round(star["z"] - center[2], 3), "luminosity_sol": star["luminosity_sol"],
         "temperature_k": star["temperature_k"]}
        for star in stars]
    return result


def _with_printed_uid(row):
    """A copy of a table row as a dict with its `uid` (a BINARY(12) column, GEN.170) as
    the text the site prints, so the row can be written as JSON."""
    out = dict(row)
    if isinstance(out.get("uid"), (bytes, bytearray)):
        out["uid"] = object_uid.format_id(object_uid.from_bytes(out["uid"]))
    return out


def phenomenon_detail(conn, phenomenon_type, phenomenon_id):
    """
    Returns one phenomenon's full row -- every column its own table has,
    plus its sector's name (see `list_phenomena`'s identical `sector_id`/
    `sector_name` convention) -- for `GET /api/phenomena/<type>/<id>`
    (`html/phenomenon.py`'s detail page). Unlike `list_phenomena`'s
    normalized `descriptor`/`radius_ly` (a common shape across all seven
    types, for a flat list), this returns the row as-is: each type has its
    own genuinely different set of fields (a nebula's `nebula_type`/
    `composition`/`formation_cause` vs. a black hole's
    `has_accretion_disk`/`hawking_temperature_k`/...), and a detail page
    showing just one phenomenon has no need to normalize them into a
    shared shape the way a mixed-type list does.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        phenomenon_type (str): One of `_PHENOMENON_TYPE_TO_TABLE`'s keys
            (`"nebula"`, `"asteroid_field"`, `"black_hole"`,
            `"neutron_star"`, `"supernova_remnant"`, `"rogue_planet"`, or
            `"interstellar_comet"`).
        phenomenon_id (int): The row's own `id` in its backing table.

    Returns:
        dict: Every column of the phenomenon's own row, plus `type`,
            `sector_name` (`None` if it has no `sector_id`), `nearest`
            (its stored nearest star systems, see `nearest_systems`) and, for a black hole, neutron star or rogue
            planet in a charted galaxy, `reference_radius_pc` (the observer its Einstein radius is worked out for). An
            asteroid field also gets `composition` (`[{"component",
            "concentration"}]`, from `asteroid_field_composition`) and an
            interstellar comet `composition` (its components, from
            `interstellar_comet_composition`), in their saved order.

    Raises:
        ValueError: If `phenomenon_type` isn't a recognized type, or no
            such row exists.
    """
    table = _PHENOMENON_TYPE_TO_TABLE.get(phenomenon_type)
    if table is None:
        raise ValueError(f"no such phenomenon type: {phenomenon_type!r}")

    row = conn.execute(
        f"""
        SELECT t.*, sec.name AS sector_name
        FROM {table} t
        LEFT JOIN sectors sec ON sec.id = t.sector_id
        WHERE t.id = ?
        """,
        (phenomenon_id,),
    ).fetchone()
    if row is None:
        raise ValueError(f"no {table} row with id {phenomenon_id}")

    detail = _with_printed_uid(row)
    detail["type"] = phenomenon_type
    if table == "rogue_planets":
        detail["gravity_g"] = rogue_planet_gravity_g(detail.get("mass_kg"), detail.get("radius_km"))
    detail["nearest"] = nearest_systems(conn, table, [phenomenon_id]).get(phenomenon_id, [])
    free = " AND star_id IS NULL" if "star_id" in detail else ""    # a black hole or neutron star of a system is not a free one
    for key, op, order in (("previous", "<", "DESC"), ("next", ">", "ASC")):
        peer = conn.execute(f"SELECT uid, name FROM {table} WHERE id {op} ?{free} ORDER BY id {order} LIMIT 1",
                            (phenomenon_id,)).fetchone()
        detail[key] = None if peer is None else {"id": _with_printed_uid(peer)["uid"], "name": peer["name"]}
    if table == "supernova_remnants":
        centre = (("black_holes", detail.get("compact_remnant_black_hole_id")),
                  ("neutron_stars", detail.get("compact_remnant_neutron_star_id")))
        for centre_table, centre_id in centre:
            if centre_id is not None:
                found = conn.execute(f"SELECT name FROM {centre_table} WHERE id = ?", (centre_id,)).fetchone()
                detail["compact_remnant_name"] = found["name"] if found else None
    if phenomenon_type in ("black_hole", "neutron_star", "rogue_planet"):
        # The Einstein radius (GEN.154) is seen from the galaxy's reference radius: no home position is stored (GEN.404).
        skeleton = get_galaxy_shape(conn)
        if skeleton is not None:
            detail["reference_radius_pc"] = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * skeleton.shape.disk_scale_length_pc
    # DB.2: the composition rows saved beside an asteroid field or an
    # interstellar comet, in their saved order (`composition_summary` is
    # the same list written out as a phrase).
    if table == "asteroid_fields":
        detail["composition"] = [
            {"component": comp["component"], "concentration": comp["concentration"]}
            for comp in conn.execute(
                "SELECT component, concentration FROM asteroid_field_composition"
                " WHERE field_id = ? ORDER BY position", (phenomenon_id,)).fetchall()
        ]
    elif table == "interstellar_comets":
        detail["composition"] = [
            comp["component"] for comp in conn.execute(
                "SELECT component FROM interstellar_comet_composition"
                " WHERE comet_id = ? ORDER BY position", (phenomenon_id,)).fetchall()
        ]
    return detail


_SCATTER_BUILT_TABLES = {
    "black-hole": ("black_holes", "black_hole", "AND star_id IS NULL"),
    "neutron-star": ("neutron_stars", "neutron_star", "AND star_id IS NULL"),
    "supernova-remnant": ("supernova_remnants", "supernova_remnant", ""),
    "planetary-nebula": ("nebulae", "nebula", ""),
    "quasar": ("quasars", "quasar", ""),
}
"""dict: A `phenomenon_scatter` kind -> `(table, phenomenon type, extra WHERE)` of the object a sector's save makes
from a row of that kind."""

SCATTER_BUILT_MATCH_PC = 0.05
"""float: How far (parsecs) a built object's centre may be from its scatter row's point and still be that row's."""


def _scattered_built_target(conn, row):
    """
    The object a built scatter row became, as `{"kind", "id", "name"}` (`kind` is a phenomenon type, or `"system"`
    for a hypervelocity star), or `None` when it cannot be found. The save does not link the two, so the object
    is found by what it kept: a hypervelocity star is the sector's hypervelocity system with the row's speed
    (the speed is a continuous draw, so two never match), every other kind the nearest object of its table in the
    sector, within `SCATTER_BUILT_MATCH_PC` of the row's point.
    """
    sector = conn.execute(
        "SELECT id FROM sectors WHERE shell_index = ? AND band_index = ? AND wedge_index = ? ORDER BY id LIMIT 1",
        (row["shell_index"], row["band_index"], row["wedge_index"])).fetchone()
    if sector is None:
        return None
    if row["kind"] == "hypervelocity-star":
        speed = math.sqrt(sum((row[f"velocity_{axis}_kms"] or 0.0) ** 2 for axis in "xyz"))
        found = conn.execute(
            "SELECT id, name FROM star_systems WHERE sector_id = ? AND runaway_class = 'hypervelocity'"
            " AND ABS(runaway_speed_kms - ?) <= ? ORDER BY ABS(runaway_speed_kms - ?), id LIMIT 1",
            (sector["id"], speed, 1e-6 * max(speed, 1.0), speed)).fetchone()
        return None if found is None else {"kind": "system", "id": found["id"], "name": found["name"]}
    entry = _SCATTER_BUILT_TABLES.get(row["kind"])
    if entry is None:
        return None
    table, phenomenon_type, extra = entry
    point = [row[f"position_{axis}_mpc"] / MPC_PER_PC for axis in "xyz"]
    candidates = conn.execute(
        f"SELECT id, name, center_x_pc, center_y_pc, center_z_pc FROM {table} WHERE sector_id = ? {extra}"
        " AND center_x_pc IS NOT NULL", (sector["id"],)).fetchall()
    best = min(
        ((detmath.dist(point, (c["center_x_pc"], c["center_y_pc"], c["center_z_pc"])), c) for c in candidates),
        key=lambda pair: (pair[0], pair[1]["id"]), default=None)
    if best is None or best[0] > SCATTER_BUILT_MATCH_PC:
        return None
    return {"kind": phenomenon_type, "id": best[1]["id"], "name": best[1]["name"]}


def scattered_phenomenon_detail(conn, scatter_id):
    """
    One `phenomenon_scatter` row as a page needs it (MAP.211): an uncharted black hole, neutron star, remnant,
    planetary nebula, quasar or hypervelocity star has no table row yet, so this is its page.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        scatter_id (int): The row's id.

    Returns:
        dict: `id`, `kind`, `subtype`, `type` (as `_SCATTER_TYPES`), `built` (bool), `address` (`[shell, band,
            wedge]`), `x`/`y`/`z` (galaxy-frame parsecs), `distance_from_core_pc`, `mass_solar` (or `None`),
            `member_count` and `system_mass_solar` (UX.114; an unbuilt row with members also has `members` and
            `nodes`, as `scatter_member_tree`), and
            for a row with a velocity `speed_kms`, `direction` (a unit vector) and `flight_years` (the time to
            fly its distance from the core at that speed, the age bound of a star ejected from the centre). A
            built row also has `target`: `{"kind", "id", "name"}` of what it became (`"system"` for a
            hypervelocity star), or `None` when that cannot be found. `sector_id`/`sector_name` are those of its
            cell when that is charted.

    Raises:
        ValueError: If no such row exists.
    """
    row = conn.execute("SELECT * FROM phenomenon_scatter WHERE id = ?", (int(scatter_id),)).fetchone()
    if row is None:
        raise ValueError(f"no phenomenon_scatter row with id {scatter_id}")
    point = [row[f"position_{axis}_mpc"] / MPC_PER_PC for axis in "xyz"]
    distance_pc = math.sqrt(sum(value * value for value in point))
    detail = {
        "id": row["id"], "kind": row["kind"], "subtype": row["subtype"],
        "type": _SCATTER_TYPES.get(row["kind"], row["kind"]), "built": row["built_at"] is not None,
        "address": [row["shell_index"], row["band_index"], row["wedge_index"]],
        "x": point[0], "y": point[1], "z": point[2], "distance_from_core_pc": distance_pc,
        "mass_solar": row["mass_solar"], **_member_fields(row),
    }
    if detail["member_count"] and not detail["built"]:
        detail.update(scatter_member_tree(conn, "phenomenon_scatter", row["id"]))
    velocity = [row[f"velocity_{axis}_kms"] for axis in "xyz"]
    if all(value is not None for value in velocity):
        speed = math.sqrt(sum(value * value for value in velocity))
        detail["speed_kms"] = speed
        detail["direction"] = [value / speed for value in velocity] if speed > 0.0 else None
        if speed > 0.0:
            detail["flight_years"] = _s_to_yr(_pc_to_m(distance_pc) / _km_to_m(speed))
    sector = conn.execute(
        "SELECT id, name FROM sectors WHERE shell_index = ? AND band_index = ? AND wedge_index = ? ORDER BY id LIMIT 1",
        (row["shell_index"], row["band_index"], row["wedge_index"])).fetchone()
    detail["sector_id"] = None if sector is None else sector["id"]
    detail["sector_name"] = None if sector is None else sector["name"]
    if detail["built"]:
        detail["target"] = _scattered_built_target(conn, row)
    return detail


def galaxy_placed_phenomena(conn):
    """
    Every galaxy-placed standalone phenomenon (all `_PHENOMENON_TABLES`
    types) -- the phenomenon counterpart to `galaxy_placed_sectors`, plotted as small
    dots on the same Galaxy Map (`planetgen/web/maps/galaxymap.py`).

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.

    Returns:
        list[dict]: See `_placed_phenomenon_rows`.
    """
    return _placed_phenomenon_rows(conn)


def nearest_systems(conn, object_table, object_ids):
    """
    The stored nearest star systems (`nearest_systems`, schema v41) of
    several objects of one kind.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        object_table (str): `'star_systems'` or a phenomenon table.
        object_ids (iterable): The objects' ids.

    Returns:
        dict: `{object_id: [{"id", "name", "distance_ly"}, ...]}`, nearest
            first; an object with no stored neighbors is left out.
    """
    object_ids = sorted(set(object_ids))
    found = {}
    for start in range(0, len(object_ids), 500):
        batch = object_ids[start:start + 500]
        marks = ", ".join("?" * len(batch))
        for row in conn.execute(
            f"SELECT n.object_id, n.neighbor_system_id, n.distance_pc, ss.name"
            f" FROM nearest_systems n JOIN star_systems ss ON ss.id = n.neighbor_system_id"
            f" WHERE n.object_table = ? AND n.object_id IN ({marks})"
            f" ORDER BY n.object_id, n.neighbor_rank",
            (object_table, *batch),
        ).fetchall():
            found.setdefault(row["object_id"], []).append({
                "id": row["neighbor_system_id"], "name": row["name"], "distance_ly": pc_to_ly(row["distance_pc"]),
            })
    return found


_PHENOMENON_TYPE_TABLES = {type_label: table for table, type_label, _descriptor, _radius in _PHENOMENON_TABLES}
"""dict: `_PHENOMENON_TABLES`' type label -> table."""


def _sector_reach_pc(sector, edge_pc):
    """
    The radius of a sphere around `sector`'s center that holds its whole
    cell: the farthest cell corner, or the cube's half diagonal for a
    sector with no grid address (or when that is bigger). A polar cap's
    cell reaches past the cube's half diagonal (its rim sits well
    over an edge from its center), so sizing this from the
    cube alone could miss a cloud that reaches into the real cell.
    """
    reach = edge_pc * math.sqrt(3) / 2
    if sector["shell_index"] is None or sector["band_index"] is None or sector["wedge_index"] is None:
        return reach
    center = (sector["center_x_pc"], sector["center_y_pc"], sector["center_z_pc"])
    for vertex in sector_cell_vertices_pc(
            sector["shell_index"], sector["band_index"], sector["wedge_index"], edge_pc):
        reach = max(reach, detmath.dist(vertex, center))
    return reach


def phenomena_near_sector(conn, sector_id):
    """
    Every galaxy-placed standalone phenomenon (any `_PHENOMENON_TABLES`
    type) generated as part of `sector_id` (its own `sector_id`), wherever
    it sits, plus every cloud from elsewhere (a phenomenon with a nonzero
    `radius_ly`) whose sphere could plausibly reach into this sector -- the
    data `planetgen/web/maps/starmap.py`'s Sector Map draws (translucent clouds for
    nebulae/asteroid fields/supernova remnants, point markers for the
    point-like types, whose own `radius_ly` is always 0 -- see
    `_PHENOMENON_TABLES`) and `planetgen/web/sector_page.py` lists alongside
    the sector's systems.

    A point-like object from another sector (a rogue planet, comet, black
    hole or neutron star) is never included: it sits in its own sector's
    cell, so on this sector's map it could only ever be drawn outside the
    wireframe (MAP.45, where neighbors' rogue planets were a third of what
    a sector drew).

    An exact cell-vs-sphere overlap test isn't worth the complexity here,
    so a cloud from elsewhere is compared against a sphere that holds the
    whole cell (`_sector_reach_pc`) instead: a safe upper bound that can
    only ever include a few extra clouds whose sphere clips that sphere
    but not the cell itself (out near a corner), never silently miss a
    real overlap. Each entry's `home` says which kind it is, so the map
    can draw a neighbor's cloud as one.

    Reads via `_placed_phenomenon_rows`'s own `bbox` argument -- a SQL
    bounding-box prefilter, not the whole-galaxy scan that function's
    default (`bbox=None`) runs -- sized to that reach plus
    whatever the widest currently-placed `radius_ly` actually is
    (`_widest_placed_phenomenon_radius_ly`), so it stays exactly as
    correct as scanning every row (a phenomenon of any size, however
    large, that could plausibly overlap is still included) while letting
    MySQL range-scan `idx_{table}_center` instead of examining every row
    in every phenomenon table on every call -- see `schema.sql`'s "v26" header
    note for the production failure this fixes (a genuine full-table-scan
    -times-four on every `GET /api/sectors/<id>`, "the exact same failure
    mode `schema.sql`'s "v25" note already documented for `sectors`).

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        sector_id (int): The `sectors.id` to check against.

    Returns:
        list[dict]: One entry per candidate phenomenon, nearest first:
            `id`, `type` (one of `_PHENOMENON_TABLES`' type labels),
            `name`, `descriptor`, `class` (letter class or `None`),
            `radius_ly`, `mass_solar` (black holes, neutron stars and quasars), `age_years` (a supernova remnant's, MAP.286; else
            `None`), `distance_ly` (sector center to phenomenon
            center), and `offset_x_ly`/`offset_y_ly`/`offset_z_ly` (the
            phenomenon's center relative to the sector's own center, in
            light-years -- the same frame `starmap.py` already places
            stars in), `octant` (the sector octant its center sits in,
            `star_systems.quadrant`'s labels), `home` (generated as part
            of this sector, rather than a neighbor's cloud reaching in) and
            `nearest` (its stored nearest star systems,
            `nearest_systems`). Empty if this sector has no galaxy
            placement of its own.

    Raises:
        ValueError: If no such sector exists.
    """
    sector = conn.execute(
        "SELECT center_x_pc, center_y_pc, center_z_pc, edge_mpc, shell_index, band_index, wedge_index"
        " FROM sectors WHERE id = ?",
        (sector_id,),
    ).fetchone()
    if sector is None:
        raise ValueError(f"no sectors row with id {sector_id}")
    if sector["center_x_pc"] is None:
        return []

    reach_pc = _sector_reach_pc(sector, mpc_to_pc(sector["edge_mpc"]))
    margin_pc = reach_pc + ly_to_pc(_widest_placed_phenomenon_radius_ly(conn))
    bbox = (sector["center_x_pc"], sector["center_y_pc"], sector["center_z_pc"], margin_pc)

    # A phenomenon generated as part of this sector always belongs in its
    # listing, even if an older placement put its center outside the cell.
    home_rows = _placed_phenomenon_rows(conn, sector_id=sector_id)
    home_keys = {(row["type"], row["id"]) for row in home_rows}
    nearby_rows = [
        row for row in _placed_phenomenon_rows(conn, bbox=bbox)
        if (row["type"], row["id"]) not in home_keys and row["radius_ly"]
    ]

    matches = []
    for phenomenon in home_rows + nearby_rows:
        dx = phenomenon["x"] - sector["center_x_pc"]
        dy = phenomenon["y"] - sector["center_y_pc"]
        dz = phenomenon["z"] - sector["center_z_pc"]
        distance_pc = math.sqrt(dx * dx + dy * dy + dz * dz)
        is_home = (phenomenon["type"], phenomenon["id"]) in home_keys
        if not is_home and distance_pc > reach_pc + ly_to_pc(phenomenon["radius_ly"]):
            continue
        center = (sector["center_x_pc"], sector["center_y_pc"], sector["center_z_pc"])
        octant, _magnitudes = classify_octant(
            galaxy_to_local_pc(center, (phenomenon["x"], phenomenon["y"], phenomenon["z"])))
        matches.append({
            "id": phenomenon["id"], "type": phenomenon["type"], "name": phenomenon["name"],
            "descriptor": phenomenon["descriptor"], "radius_ly": phenomenon["radius_ly"],
            "class": phenomenon["class"], "gravity_g": phenomenon.get("gravity_g"),
            "mass_solar": phenomenon.get("mass_solar"), "age_years": phenomenon.get("age_years"),
            "distance_ly": pc_to_ly(distance_pc),
            "offset_x_ly": pc_to_ly(dx), "offset_y_ly": pc_to_ly(dy), "offset_z_ly": pc_to_ly(dz),
            "octant": octant,
            "home": is_home,
        })
    by_type = {}
    for match in matches:
        by_type.setdefault(match["type"], []).append(match["id"])
    for kind, ids in by_type.items():
        stored = nearest_systems(conn, _PHENOMENON_TYPE_TABLES[kind], ids)
        for match in matches:
            if match["type"] == kind:
                match["nearest"] = stored.get(match["id"], [])
    matches.sort(key=lambda match: match["distance_ly"])
    return matches


def _heliopause_au(conn, system, cloud):
    """
    `(open-space heliopause, heliopause inside `cloud`)` in AU for a
    `star_systems` row: a close pair's shared bubble, otherwise the
    primary (or single) star's. `(None, None)` when nothing stores one.
    """
    if system["binary_configuration"] == "close" and system["binary_heliosphere_radius_km"] is not None:
        radius_km = system["binary_heliosphere_radius_km"]
    else:
        star = conn.execute(
            "SELECT heliosphere_radius_km FROM stars WHERE star_system_id = ? AND role IN ('primary', 'single')"
            " ORDER BY id LIMIT 1", (system["id"],)).fetchone()
        radius_km = star["heliosphere_radius_km"] if star else None
    if radius_km is None:
        return None, None
    open_space_au = _km_to_au(radius_km)
    if cloud is None:
        return open_space_au, open_space_au
    return open_space_au, compressed_heliosphere_radius(open_space_au, cloud["density_cm3"], cloud["temperature_k"])


def system_detail(conn, system_id):
    """
    Returns one system's full web-display detail: stars, planets (each
    with its own moons), asteroid belts, comets, description text, and
    enough sector context to render `html/system.py` -- in one function, the
    same DB-row-shaped read `sector_detail` gives sectors (see that
    function's docstring for why this is distinct from
    `planetgen.db.store.load_star_system`'s generation object graph).

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        system_id (int): The `star_systems.id` to look up.

    Returns:
        dict: `id`, `name`, `sector_id`, `quadrant`, `location`,
            `is_binary`, `binary_type`, `binary_configuration` (`'close'`,
            `'wide'`, or `None` -- see `schema.sql`'s "v15" note),
            `binary_separation_km` and `binary_heliosphere_radius_km` (a
            close pair's, for orbits around it),
            `binary_mutual_position_x/y/z_km` (the secondary's position
            relative to the primary -- NULL for a single star; see
            `schema.sql`'s "v14"/"v15" notes), `stars` (id/role/name/
            star_type/mass_kg/radius_km/temperature_k/luminosity_w/
            heliosphere_radius_km/offset_x_au/offset_y_au/offset_z_au, the
            offsets set for a `member` star only, and `parent_node_no`/`slot`, the node of the orbit tree the star hangs in and which of its
            two places it holds, `None` for a single star --
            `id` matches a `'wide'` binary's `planets`/`belts` rows' own
            `star_id`, disambiguating which star each orbits), `orbits` (the system's `orbit_nodes` rows in node order, MAP.213:
            `node_no`, `parent_node_no`, `slot`, `kind`, `mass_solar`, `separation_km`, `eccentricity`, `inclination_deg`,
            `period_years`, `interaction` and `flags`; empty for a single star), `planets`
            (each a `planets` row, including its own `star_id`, plus its
            own `moons` list; every planet and moon also carries
            `habitable`, `life_stage` and `inhabited` (a colony counts --
            see `_with_life_fields`), `belts` (`asteroid_belts` rows, including
            `star_id`, plus `composition`: its `asteroid_belt_composition`
            rows as `{component, concentration}`, largest share first), `comets` (`comets` rows, including `star_id` --
            no `orbital_index`, see `insert_comet`'s docstring), and
            `sector_siblings` (`{id, name}` for every other system in the
            same sector, empty if standalone -- for linkifying
            `location`'s "nearest: ..." names without a second round
            trip), and `nearest_neighbors` (`{id, name, distance_ly}` for
            the closest same-sector systems, nearest first -- see
            `_nearest_sector_siblings`), `inside` (`containing_cloud`), and
            `heliopause_au` (the system's heliopause, pressed in by the
            cloud it sits inside, if any -- the System Local Frame's edge for
            navigation) with `heliopause_open_space_au` (the same without
            the cloud); both `None` when no star row has one.

    Raises:
        ValueError: If no such system exists.
    """
    system = conn.execute("SELECT * FROM star_systems WHERE id = ?", (system_id,)).fetchone()
    if system is None:
        raise ValueError(f"no star_systems row with id {system_id}")

    stars = conn.execute(
        "SELECT id, role, name, star_type, yerkes_class, mass_kg, radius_km, temperature_k, luminosity_w, heliosphere_radius_km,"
        " offset_x_au, offset_y_au, offset_z_au, parent_node_id AS parent_node_no, slot"
        " FROM stars WHERE star_system_id = ?"
        " ORDER BY CASE role WHEN 'primary' THEN 0 WHEN 'single' THEN 0 WHEN 'secondary' THEN 1 ELSE 2 END, id",
        (system_id,),
    ).fetchall()
    orbits = [dict(r) for r in conn.execute(
        "SELECT node_no, parent_node_no, slot, kind, mass_solar, separation_km, eccentricity, inclination_deg, period_years,"
        " interaction, flags FROM orbit_nodes WHERE star_system_id = ? ORDER BY node_no", (system_id,)).fetchall()]
    # mass_kg (already selected above) is what lets planetgen/web/maps/systemmap.py
    # split binary_mutual_position_x/y/z_km below into each star's own
    # mass-weighted offset from the system's barycenter, rather than
    # (incorrectly) anchoring the whole system on the primary alone.

    planet_rows = with_gases(conn, conn.execute(
        "SELECT * FROM planets WHERE star_system_id = ? ORDER BY orbital_index", (system_id,)
    ).fetchall(), False)
    colonized = colonized_body_ids(conn, system_id)
    moon_rows = with_gases(conn, conn.execute(
        "SELECT * FROM moons WHERE star_system_id = ? ORDER BY planet_id, orbital_index", (system_id,)
    ).fetchall(), True)
    planet_scores = body_class_scores(conn, [planet["id"] for planet in planet_rows], False)
    moon_scores = body_class_scores(conn, [moon["id"] for moon in moon_rows], True)
    planet_living = body_living(conn, [planet["id"] for planet in planet_rows], False)
    moon_living = body_living(conn, [moon["id"] for moon in moon_rows], True)
    moons_by_planet = {}
    for moon in moon_rows:
        moons_by_planet.setdefault(moon["planet_id"], []).append(moon)
    planets = []
    hosts = {}  # one lookup per host star, not per planet
    for planet in planet_rows:
        if planet["star_id"] not in hosts:
            hosts[planet["star_id"]] = compact_host_of_star(conn, planet["star_id"])
        host = hosts[planet["star_id"]]
        planet_dict = _with_life_fields(_with_printed_uid(planet), colonized["planets"], host, planet_scores.get(planet["id"]),
                                       planet_living.get(planet["id"]))
        planet_dict["moons"] = [_with_life_fields(_with_printed_uid(m), colonized["moons"], host, moon_scores.get(m["id"]), moon_living.get(m["id"]))
                                for m in moons_by_planet.get(planet["id"], [])]
        planets.append(planet_dict)

    belts = [_with_printed_uid(b) for b in conn.execute(
        "SELECT * FROM asteroid_belts WHERE star_system_id = ? ORDER BY orbital_index", (system_id,)
    ).fetchall()]
    belt_composition = {}
    for row in conn.execute(
        "SELECT c.belt_id, c.component, c.concentration FROM asteroid_belt_composition c"
        " JOIN asteroid_belts b ON b.id = c.belt_id WHERE b.star_system_id = ?"
        " ORDER BY c.belt_id, c.position",
        (system_id,),
    ).fetchall():
        belt_composition.setdefault(row["belt_id"], []).append(
            {"component": row["component"], "concentration": row["concentration"]})
    for belt in belts:
        belt["composition"] = belt_composition.get(belt["id"], [])

    comets = conn.execute(
        "SELECT * FROM comets WHERE star_system_id = ? ORDER BY id", (system_id,)
    ).fetchall()

    sector_siblings = []
    nearest_neighbors = []
    if system["sector_id"] is not None:
        sibling_rows = conn.execute(
            "SELECT id, name, position_x_mpc, position_y_mpc, position_z_mpc"
            " FROM star_systems WHERE sector_id = ?",
            (system["sector_id"],),
        ).fetchall()
        sector_siblings = [{"id": r["id"], "name": r["name"]} for r in sibling_rows]
        # The stored neighbors (schema v41) are searched across sector
        # boundaries; a system generated before them falls back to its
        # own sector's.
        nearest_neighbors = (nearest_systems(conn, "star_systems", [system_id]).get(system_id)
                             or _nearest_sector_siblings(system, sibling_rows))

    cloud = containing_cloud(conn, system)
    open_space_au, heliopause_au = _heliopause_au(conn, system, cloud)

    return {
        "id": system["id"], "name": system["name"], "sector_id": system["sector_id"],
        "quadrant": system["quadrant"], "location": system["location"],
        "is_binary": system["is_binary"], "binary_type": system["binary_type"],
        "binary_configuration": system["binary_configuration"],
        "binary_separation_km": system["binary_separation_km"],
        "binary_heliosphere_radius_km": system["binary_heliosphere_radius_km"],
        "binary_mutual_position_x_km": system["binary_mutual_position_x_km"],
        "binary_mutual_position_y_km": system["binary_mutual_position_y_km"],
        "binary_mutual_position_z_km": system["binary_mutual_position_z_km"],
        "runaway_class": system["runaway_class"], "runaway_speed_kms": system["runaway_speed_kms"],
        "stars": [dict(s) for s in stars],
        "orbits": orbits,
        "planets": planets,
        "belts": belts,
        "comets": [_with_printed_uid(c) for c in comets],
        "sector_siblings": sector_siblings,
        "nearest_neighbors": nearest_neighbors,
        "inside": cloud,
        "heliopause_au": heliopause_au,
        "heliopause_open_space_au": open_space_au,
    }


def _life_classes_of(body, living):
    """The `life_classes` list of a planet or moon row (GEN.1209): one `{life_class_id, share, life_stage, animal_grade, kinds}` for each class that lives on it,
    the dominant one first, from its `planet_life` / `moon_life` `living` rows (`store.body_living`). A body stored before they existed has its one class from its own
    columns, with the whole share; one with no living class has an empty list."""
    class_id = life_classes.living_class_id(body.get("life_class"), body.get("life_class_variant"), body.get("life_stage"))
    if class_id is None:
        return []
    if not living:
        return [{"life_class_id": class_id, "share": 1.0, "life_stage": body["life_stage"], "animal_grade": body.get("life_animal_grade"),
                 "kinds": {"microbial": body.get("life_microbial"), "flora": body.get("life_flora"), "fauna": body.get("life_fauna"),
                           "sentient": body.get("life_sentient")}}]
    ordered = sorted(living, key=lambda entry: entry["life_class"] != class_id)
    return [{"life_class_id": entry["life_class"], "share": entry["share"], "life_stage": entry["stage"], "animal_grade": entry["animal_grade"],
             "kinds": dict(zip(("microbial", "flora", "fauna", "sentient"), entry["kinds"]))} for entry in ordered if entry["stage"] is not None]


def _with_life_fields(body, colonized=(), compact_host=None, class_scores=None, living=None):
    """
    Adds the system page's per-body life summary to one `planets`/`moons`
    row dict: `habitable` (its overall PHI-6 colour is Blue or Green, `habitability.is_habitable`, the same test
    `StarSystem.count_habitable` uses), `habitable_class` (its class is one of `tuning.HABITABLE_PLANET_CLASSES`),
    `life_stage` (the most advanced evolutionary milestone its timeline reached, or `None`: the stored column, kept
    for a habitable class only) and
    `inhabited` (reached a technological civilization --
    the page only ever describes a timeline for a habitable class, see
    `Planet._generate_life_and_flavor_paragraphs` -- or a colony stands on
    it: `colonized` is that table's ids from `colonized_body_ids`, schema
    v42), `life_class_id` (the class that lives on it as one id such as `"1a"`, or `None`; none for a class the body suits but has no life stage yet, GEN.1001) and `life_class_scores`
    (`{class id: score 1 to 9}` for every class it supports, from `class_scores`, the body's
    `planet_life_class_scores` or `moon_life_class_scores` rows; empty when it supports none, API.37) and `life_classes` (the list of the classes that live on it with their
    shares, stages, animal grades and kind scores, the dominant one first, from `living`, its `planet_life` or `moon_life` rows, GEN.1209).
    """
    body["life_class_id"] = life_classes.living_class_id(body.get("life_class"), body.get("life_class_variant"), body.get("life_stage"))
    body["life_class_scores"] = dict(sorted((class_scores or {}).items()))
    body["life_classes"] = _life_classes_of(body, living)
    body["habitable_class"] = body["planet_class"] in HABITABLE_PLANET_CLASSES
    body["habitable"] = habitability.is_habitable(body.get("phi4"))
    body["life_stage"] = body["life_stage"] if body["habitable_class"] else None
    body["compact_host"] = compact_host
    body["equipment"] = equipment.summary(body)
    body["inhabited"] = body["life_stage"] == "technological_civilization" or body["id"] in colonized
    return body


NEAREST_NEIGHBOR_COUNT = 3
"""int: How many nearest same-sector systems `system_detail` lists -- the
same count the stored `location` string was generated with
(`store._location_for_entry`)."""


def _nearest_sector_siblings(system, sibling_rows, count=NEAREST_NEIGHBOR_COUNT):
    """
    The `count` systems in `system`'s own sector nearest to it, computed
    from their current rows rather than parsed back out of the stored
    `location` string. That string is written at generation time, before
    each neighbor's name has been made unique (`store.insert_star_system`'s
    `reserve_system_name` can still add a suffix afterwards) and never
    updated when a system is renamed, so the names it lists often no
    longer match any row and can't be linked.

    Args:
        system (dict-like): The `star_systems` row being shown.
        sibling_rows (list): `id`/`name`/`position_*_mpc` rows for every
                             system in the same sector (including
                             `system` itself, which is skipped).
        count (int): See `NEAREST_NEIGHBOR_COUNT`.

    Returns:
        list[dict]: Nearest first, each with `id`, `name`, `distance_ly`.
            Empty when `system` has no stored position.
    """
    origin = (system["position_x_mpc"], system["position_y_mpc"], system["position_z_mpc"])
    if any(v is None for v in origin):
        return []
    neighbors = []
    for row in sibling_rows:
        position = (row["position_x_mpc"], row["position_y_mpc"], row["position_z_mpc"])
        if row["id"] == system["id"] or any(v is None for v in position):
            continue
        distance_ly = milliparsecs_to_ly(detmath.dist(origin, position))
        neighbors.append({"id": row["id"], "name": row["name"], "distance_ly": distance_ly})
    neighbors.sort(key=lambda entry: entry["distance_ly"])
    return neighbors[:count]


def galaxy_placed_sectors(conn):
    """
    Every sector with a galaxy position, plus its live system count -- the
    data the `/galaxy` Galaxy Map (`planetgen/web/maps/galaxymap.py`) plots.
    Unplaced sectors (`center_x_pc IS NULL`) have nothing to plot and are
    excluded at the query itself.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.

    Returns:
        list[dict]: `id`, `name`, `x`/`y`/`z` (`center_x/y/z_pc`),
            `galactic_radius_pc`, `shell_index`, `system_count`.
    """
    rows = conn.execute(
        """
        SELECT sec.id, sec.name, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.galactic_radius_pc, sec.shell_index, COALESCE(counts.n, 0) AS system_count
        FROM sectors sec
        LEFT JOIN (SELECT sector_id, COUNT(*) AS n FROM star_systems
                   WHERE sector_id IS NOT NULL GROUP BY sector_id) counts ON counts.sector_id = sec.id
        WHERE sec.center_x_pc IS NOT NULL
        ORDER BY sec.galactic_radius_pc
        """,
    ).fetchall()
    return [
        {
            "id": r["id"], "name": r["name"],
            "x": r["center_x_pc"], "y": r["center_y_pc"], "z": r["center_z_pc"],
            "galactic_radius_pc": r["galactic_radius_pc"], "shell_index": r["shell_index"],
            "system_count": r["system_count"],
        }
        for r in rows
    ]


MADE_SECTOR_LIMIT = 2000
"""int: Most sectors `sectors_made` returns by address: a bigger run is shown
by its first sectors (the response still says the real total)."""





def sectors_made(conn, since, until=None, limit=MADE_SECTOR_LIMIT):
    """
    The galaxy-placed sectors created from `since` to `until` (Unix
    seconds; `until=None` means now) -- what a generate job made, for the
    Galaxy Map's "Show on Galaxy Map" (ADM.31).

    Returns:
        dict: `total` (how many), `items` (up to `limit`, lowest id first:
            `id`, `name`, `x`/`y`/`z` pc, `shell_index`, `band_index`,
            `wedge_index`).
    """
    where = "center_x_pc IS NOT NULL AND UNIX_TIMESTAMP(created_at) >= ?"
    params = [since]
    if until is not None:
        where += " AND UNIX_TIMESTAMP(created_at) <= ?"
        params.append(until)
    total = conn.execute(f"SELECT COUNT(*) AS n FROM sectors WHERE {where}", params).fetchone()["n"]
    rows = conn.execute(
        "SELECT id, name, center_x_pc, center_y_pc, center_z_pc, shell_index, band_index, wedge_index "
        f"FROM sectors WHERE {where} ORDER BY id LIMIT ?", [*params, int(limit)],
    ).fetchall()
    return {
        "total": total,
        "items": [
            {"id": r["id"], "name": r["name"], "x": r["center_x_pc"], "y": r["center_y_pc"], "z": r["center_z_pc"],
             "shell_index": r["shell_index"], "band_index": r["band_index"], "wedge_index": r["wedge_index"]}
            for r in rows
        ],
    }


def galaxy_density_shape(conn):
    """
    The galaxy's stored density-skeleton shape (the singleton
    `galaxy_shape` row `planetgen plan` writes -- see
    `planetgen.galaxy.density.GalaxyShape` and `store.get_galaxy_shape`),
    serialized to a plain JSON-able dict. This is the real
    exponential-disk-plus-bulge-plus-spiral-arm model already used to gate
    and weight actual sector generation (`planetgen`'s `_BatchDensity`);
    exposing it here lets the Galaxy Map (`planetgen/web/maps/galaxymap.py`) shade
    its "expected density" cloud from this same model instead of a
    generic illustrative gradient, so un-generated space still reads as
    the spiral it's predicted to be.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.

    Returns:
        dict or None: Every `GalaxyShape` field and
            `galaxy.density.model_terms` and `halo_terms` entry, plus `edge_pc`,
            `outer_shell_index`, `expected_system_count_at_density_1`,
            `pinned_stellar_mass_solar` and `pins` (DB.81: the stellar mass the plan pinned, or `None`, and the
            levers it was told, `{name: value}`) and
            `half_height_pc` (how far above and below the plane the
            layers reach; left out before the layers are stored).
            `None` if `planetgen plan` has never been run against this
            database (no `galaxy_shape` row yet).
    """
    skeleton = get_galaxy_shape(conn)
    if skeleton is None:
        return None
    shape = shape_with_terms(skeleton.shape)
    shape["edge_pc"] = skeleton.edge_pc
    shape["outer_shell_index"] = skeleton.outer_shell_index
    shape["expected_system_count_at_density_1"] = skeleton.expected_system_count_at_density_1
    shape["pinned_stellar_mass_solar"] = skeleton.pinned_stellar_mass_solar
    shape["pins"] = skeleton.pins
    reach = conn.execute("SELECT shell_index, high_band_index FROM galaxy_shell").fetchall()
    if reach:
        # How far above the plane the outline goes: the northern edge of each shell's highest band, on its outer face.
        shape["half_height_pc"] = max(
            spherical.shell_bounds_pc(row["shell_index"], skeleton.edge_pc)[1]
            * abs(math.cos(spherical.band_colatitude_bounds(row["shell_index"], row["high_band_index"])[0]))
            for row in reach)
    return shape


# ---------------------------------------------------------------------
# Interactive 3D Galaxy Map viewport queries -- the placed-sector shape
# the node readers below reuse. Unlike `galaxy_placed_sectors`/`galaxy_density_shape` above (each
# called once per page load for the flat, whole-galaxy overview map),
# these are scoped to a moving viewport.
# ---------------------------------------------------------------------

GALAXY_VIEW_MAX_PLACED = 2000
"""int: Cap on how many placed (already-generated) sectors
`galaxy_sectors_in_view` returns, closest-first -- a view centered on a
heavily-generated region could otherwise return an unbounded response."""


def galaxy_sectors_in_view(conn, center_x_pc, center_y_pc, center_z_pc, radius_pc, limit=GALAXY_VIEW_MAX_PLACED):
    """
    Every galaxy-placed sector within `radius_pc` of `(center_x_pc,
    center_y_pc, center_z_pc)`, closest first -- the "placed" tier of the
    interactive 3D Galaxy Map's live viewport (see `planetgen.galaxy.
    viewport`'s module docstring), unlike `galaxy_placed_sectors`
    (the whole galaxy, once, for the flat overview map's own quadrant/ring
    tables).

    `sectors.center_x/y/z_pc` are covered by a composite index
    (`idx_sectors_center` -- schema v25; see `schema.sql`'s "v25" header
    note) MySQL can range-scan on the leading `center_x_pc` column for the
    bounding-box `WHERE` clause below, rather than a full table scan --
    this was a genuine, confirmed-in-production full-table-scan-per-call
    before that index existed, the same failure mode v22's own note
    documents for the pre-v22 `/api/search`. The exact Euclidean-distance
    filter and closest-`limit` truncation are still done in Python once
    that's narrowed the row-fetch volume, the same division of labor
    `phenomena_near_sector` already uses for its own bounding-sphere
    prefilter.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        center_x_pc, center_y_pc, center_z_pc (float): The view center,
            galaxy-frame parsecs.
        radius_pc (float): The view radius, parsecs.
        limit (int): See `GALAXY_VIEW_MAX_PLACED`.

    Returns:
        list[dict]: `id`, `name`, `x`/`y`/`z`, `galactic_radius_pc`,
            `shell_index`, `band_index`, `wedge_index`, `designation`
            (`provisional_sector_designation`, `None` if this sector has
            no grid address), `system_count`, `edge_ly` (this
            sector's own real edge length -- lets a client compute its
            true stellar density, `system_count / edge_ly ** 3`, rather
            than just its raw system count), `distance_pc` (from the
            given center).
    """
    rows = conn.execute(
        """
        SELECT sec.id, sec.name, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.galactic_radius_pc, sec.shell_index, sec.band_index, sec.wedge_index,
               sec.edge_mpc,
               (SELECT COUNT(*) FROM star_systems ss WHERE ss.sector_id = sec.id) AS system_count
        FROM sectors sec
        WHERE sec.center_x_pc IS NOT NULL
          AND sec.center_x_pc BETWEEN ? AND ?
          AND sec.center_y_pc BETWEEN ? AND ?
          AND sec.center_z_pc BETWEEN ? AND ?
        """,
        (
            center_x_pc - radius_pc, center_x_pc + radius_pc,
            center_y_pc - radius_pc, center_y_pc + radius_pc,
            center_z_pc - radius_pc, center_z_pc + radius_pc,
        ),
    ).fetchall()

    radius_sq = radius_pc * radius_pc
    candidates = []
    for r in rows:
        dx = r["center_x_pc"] - center_x_pc
        dy = r["center_y_pc"] - center_y_pc
        dz = r["center_z_pc"] - center_z_pc
        distance_sq = dx * dx + dy * dy + dz * dz
        if distance_sq > radius_sq:
            continue
        candidates.append((distance_sq, r))
    candidates.sort(key=lambda pair: pair[0])
    candidates = candidates[:limit]

    results = []
    for distance_sq, r in candidates:
        entry = _placed_sector_entry(r, r["system_count"])
        entry["distance_pc"] = math.sqrt(distance_sq)
        results.append(entry)
    return results


# ---------------------------------------------------------------------
# Box readers -- what `db.nodes.galaxy_nodes` reads a node's box with (GET /api/galaxy/nodes), and GET
# /api/galaxy/stamp. The `GALAXY_TILE_*` names are the readers' caps, kept from the cube tiles: each cap now applies to
# one node's box.
# ---------------------------------------------------------------------

GALAXY_TILE_MAX_PLACED = 250
"""int: Most placed sectors one box returns. Past this, a box returns an
evenly spaced sample (every Nth sector by id, so the same sample every
time) -- a zoomed-out view can't show thousands of sub-pixel sectors
anyway, and zooming in switches to smaller nodes that hold the rest.

The sample has to be spread out, not just the first ids: a local
neighborhood is generated in order outward from the core, so its lowest
ids are its core-facing side. Taking the lowest 250 used to draw a
generated 100 ly sphere as a core-facing bowl, cut off flat where the cap
ran out. A fully generated node can hold thousands of sectors, so this
isn't only a far-zoomed-out case."""

GALAXY_TILE_MAX_CLOUDS = 200
"""int: Most nebulae and supernova remnants one node lists
(`galaxy_clouds_in_box`), the largest first. A node wide enough to hold
more shows the map at a scale where the smaller ones are under a pixel
anyway."""

def galaxy_sectors_in_box(conn, lo, hi, limit=GALAXY_TILE_MAX_PLACED):
    """
    Placed sectors whose center lies in the half-open box `[lo, hi)`,
    ordered by id, at most `limit` -- the "placed" half of one node. When
    the box holds more than `limit`, it's every Nth of them by id, `N =
    ceil(count / limit)` (see `GALAXY_TILE_MAX_PLACED` for why not the
    lowest ids).

    Two queries rather than `galaxy_sectors_in_view`'s one: the first
    reads only the (indexed, `idx_sectors_center`) sector columns, and the
    per-sector system count runs only for the sectors actually returned,
    so a box covering the whole galaxy costs one read of the sector
    columns instead of a count for every placed sector.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        lo (tuple): `(x, y, z)` inclusive lower corner, parsecs.
        hi (tuple): `(x, y, z)` exclusive upper corner, parsecs.
        limit (int): See `GALAXY_TILE_MAX_PLACED`.

    Returns:
        list[dict]: The same entries `galaxy_sectors_in_view` returns,
            minus `distance_pc`.
    """
    rows = conn.execute(
        """
        SELECT id, name, center_x_pc, center_y_pc, center_z_pc,
               galactic_radius_pc, shell_index, band_index, wedge_index, edge_mpc
        FROM (
            SELECT sec.id, sec.name, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
                   sec.galactic_radius_pc, sec.shell_index, sec.band_index, sec.wedge_index,
                   sec.edge_mpc,
                   ROW_NUMBER() OVER (ORDER BY sec.id) AS row_num,
                   COUNT(*) OVER () AS row_total
            FROM sectors sec
            WHERE sec.center_x_pc >= ? AND sec.center_x_pc < ?
              AND sec.center_y_pc >= ? AND sec.center_y_pc < ?
              AND sec.center_z_pc >= ? AND sec.center_z_pc < ?
        ) ranked
        WHERE MOD(row_num - 1, CEILING(row_total / ?)) = 0
        ORDER BY id
        LIMIT ?
        """,
        (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], int(limit), int(limit)),
    ).fetchall()
    if not rows:
        return []

    ids = [r["id"] for r in rows]
    placeholders = ", ".join("?" for _ in ids)
    counts = {
        c["sector_id"]: c["system_count"]
        for c in conn.execute(
            f"SELECT sector_id, COUNT(*) AS system_count FROM star_systems "
            f"WHERE sector_id IN ({placeholders}) GROUP BY sector_id",
            ids,
        ).fetchall()
    }
    return [_placed_sector_entry(r, counts.get(r["id"], 0)) for r in rows]


def _placed_sector_entry(r, system_count):
    """One placed-tier dict (see `galaxy_sectors_in_view`'s Returns) from a
    `sectors` row, without `distance_pc`."""
    address = sector_address(r)
    edge_pc = mpc_to_pc(r["edge_mpc"]) if r["edge_mpc"] else None
    return {
        "id": r["id"], "name": r["name"],
        "x": r["center_x_pc"], "y": r["center_y_pc"], "z": r["center_z_pc"],
        "galactic_radius_pc": r["galactic_radius_pc"],
        "shell_index": r["shell_index"], "band_index": r["band_index"], "wedge_index": r["wedge_index"],
        "designation": provisional_sector_designation(*address) if address is not None else None,
        "system_count": system_count,
        "edge_ly": pc_to_ly(edge_pc) if edge_pc else None,
    }


def galaxy_sector_count_in_box(conn, lo, hi):
    """
    How many placed sectors have their center in the half-open box
    `[lo, hi)`: unlike a node's `placed` list, nothing is left out, so the
    generated-star budget a node shares out (`generated_star_budget`) goes
    by every sector the node holds.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        lo (tuple): `(x, y, z)` inclusive lower corner, parsecs.
        hi (tuple): `(x, y, z)` exclusive upper corner, parsecs.

    Returns:
        int: The count.
    """
    return int(conn.execute(
        "SELECT COUNT(*) AS n FROM sectors WHERE center_x_pc >= ? AND center_x_pc < ? AND center_y_pc >= ? "
        "AND center_y_pc < ? AND center_z_pc >= ? AND center_z_pc < ? AND shell_index IS NOT NULL",
        (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]),
    ).fetchone()["n"])


_CLOUD_TABLES = (
    ("nebulae", "nebula", "nebula_type"),
    ("supernova_remnants", "supernova_remnant", "morphology"),
)
"""tuple: `(table, type_label, descriptor_column)` for the phenomena the
Galaxy Map draws as clouds -- the ones with a real extent worth seeing
at galaxy scale (an asteroid field is far smaller)."""


def _cloud_margins_pc(conn):
    """`{table: widest placed radius in pc}` for each `_CLOUD_TABLES`
    table with a placed row -- `galaxy_clouds_in_box`'s prefilter margin,
    read once per node request."""
    margins = {}
    for table, _type_label, _descriptor_column in _CLOUD_TABLES:
        widest = conn.execute(
            f"SELECT MAX(radius_ly) AS r FROM {table} WHERE center_x_pc IS NOT NULL"
        ).fetchone()
        if widest is not None and widest["r"] is not None:
            margins[table] = ly_to_pc(float(widest["r"]))
    return margins


def galaxy_clouds_in_box(conn, lo, hi, max_clouds=GALAXY_TILE_MAX_CLOUDS, margins=None):
    """
    Every placed nebula and supernova remnant whose sphere reaches into the
    box `[lo, hi)` -- not just the ones centered there, so a big cloud
    shows from every node it covers (the map drops the repeats by type and
    id). The Galaxy Map (`static/galaxymap3d.js`) draws each as a
    translucent cloud.

    Each table is read with a bounding-box prefilter widened by its own
    widest radius (so `idx_<table>_center` can range-scan it), then the
    sphere is tested against the box exactly.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        lo (tuple): `(x, y, z)` lower corner, parsecs.
        hi (tuple): `(x, y, z)` upper corner, parsecs.
        max_clouds (int): See `GALAXY_TILE_MAX_CLOUDS`.
        margins (dict, optional): An already-read `_cloud_margins_pc`.

    Returns:
        list[dict]: Largest first: `type` (`"nebula"` or
            `"supernova_remnant"`), `id`, `name`, `descriptor` (a nebula's
            type, a remnant's morphology), `class` (letter class or
            `None`), `radius_pc`, and `x`/`y`/`z` (center, parsecs).
    """
    if margins is None:
        margins = _cloud_margins_pc(conn)
    clouds = []
    for table, type_label, descriptor_column in _CLOUD_TABLES:
        if table not in margins:
            continue
        margin = margins[table]
        rows = conn.execute(
            f"""
            SELECT id, name, {descriptor_column} AS descriptor, radius_ly,
                   {_PHENOMENON_CLASS_COLUMNS[table]} AS class_code,
                   center_x_pc, center_y_pc, center_z_pc
            FROM {table}
            WHERE center_x_pc BETWEEN ? AND ? AND center_y_pc BETWEEN ? AND ?
              AND center_z_pc BETWEEN ? AND ?
            """,
            (lo[0] - margin, hi[0] + margin, lo[1] - margin, hi[1] + margin, lo[2] - margin, hi[2] + margin),
        ).fetchall()
        for row in rows:
            radius_pc = ly_to_pc(float(row["radius_ly"] or 0.0))
            center = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
            gap = math.sqrt(sum(max(lo[i] - center[i], 0.0, center[i] - hi[i]) ** 2 for i in range(3)))
            if radius_pc <= 0 or gap > radius_pc:
                continue
            clouds.append({
                "type": type_label, "id": row["id"], "name": row["name"],
                "descriptor": row["descriptor"], "class": row["class_code"], "radius_pc": radius_pc,
                "x": center[0], "y": center[1], "z": center[2],
            })
    clouds.sort(key=lambda cloud: (-cloud["radius_pc"], cloud["type"], cloud["id"]))
    return clouds[:max_clouds]


GALAXY_TILE_BRIGHTEST_SAMPLE = 100000
"""int: How many of the galaxy's most luminous bright stars `db.nodes.galaxy_nodes`
reads in one go (once per request, only when a node needs it) to answer
nodes too big to query on their own -- see
`GALAXY_TILE_BRIGHT_STAR_ROW_BUDGET`. A node of a few kiloparsecs holds a
few hundred of them, so zoomed out every node still fills its
`GALAXY_TILE_BRIGHT_BUDGET`."""

GALAXY_TILE_BRIGHT_STAR_ROW_BUDGET = 50000
"""int: Most `bright_stars` rows one node's own query may expect to read
(the cheaper of `galaxy_bright_stars_in_box`'s two ways). A node that
would read more -- a big node with millions of stars, or one mostly past
the galaxy's edge, where the luminosity walk finds few stars inside --
takes its stars from the galaxy-wide sample instead: the brightest
stars inside it, thinned by luminosity when the sample holds fewer than
the node's cap there. Before this a zoomed-out request could run past
the web app's 30 s API timeout and the stars vanished (MAP.47)."""

BRIGHT_STAR_MAX_EXACT_RANGES = 2000
"""int: Most `(shell, band)` index ranges one node's bright-star query
lists (each with its wedges) before it reads each shell's bands in one range
instead."""


def _box_angle_intervals(lo, hi):
    """
    The galactic longitudes (counterclockwise from +X, in `[0, 2*pi)`) the
    box `[lo, hi)` covers seen from the axis, as a list of `(start, end)`
    radians -- one interval, or two when it wraps past +X, or the whole
    circle when the box holds (or touches) the axis.
    """
    full = [(0.0, 2 * math.pi)]
    if lo[0] <= 0 <= hi[0] and lo[1] <= 0 <= hi[1]:
        return full
    # A convex region clear of the axis spans less than half a turn, and
    # its extreme longitudes are at its corners.
    center = math.atan2((lo[1] + hi[1]) / 2, (lo[0] + hi[0]) / 2)
    offsets = [
        (math.atan2(y, x) - center + math.pi) % (2 * math.pi) - math.pi
        for x in (lo[0], hi[0]) for y in (lo[1], hi[1])
    ]
    start = (center + min(offsets)) % (2 * math.pi)
    end = start + (max(offsets) - min(offsets))
    if end <= 2 * math.pi:
        return [(start, end)]
    return [(start, 2 * math.pi), (0.0, end - 2 * math.pi)]


def _band_wedge_ranges(shell, band, intervals):
    """`intervals`' wedges of one band, as inclusive `(first, last)` pairs
    (rounded outward, so every wedge a longitude touches is in)."""
    n = spherical.wedge_count(shell, band)
    if n == 1:
        return [(0, 0)]
    ranges = []
    for start, end in intervals:
        first = max(0, int(math.floor(start * n / (2 * math.pi))))
        last = min(n - 1, int(math.floor(end * n / (2 * math.pi))))
        ranges.append((first, last))
    return ranges


def _box_colatitude_range(lo, hi):
    """`(north, south)` colatitudes, radians from galactic north, the box
    `[lo, hi)` spans: the extremes sit at the corners of the box's
    `(distance from the axis, z)` rectangle."""
    near = detmath.hypot(min(max(0.0, lo[0]), hi[0]), min(max(0.0, lo[1]), hi[1]))
    far = detmath.hypot(max(abs(lo[0]), abs(hi[0])), max(abs(lo[1]), abs(hi[1])))
    angles = [math.atan2(rho, z) for rho in (near, far) for z in (lo[2], hi[2])]
    return min(angles), max(angles)


def _bright_star_bands(conn, lo, hi, edge_pc):
    """
    `([(shell, band_min, band_max), ...], share)`: the shells the box
    `[lo, hi)` reaches, each with the bands it reaches there (the box's
    colatitudes at that shell's radii), cut to the stored outline
    (`galaxy_shell`) when there is one -- no bright star sits outside it --
    and the share of the outline's sectors those bands hold (about; it
    counts every wedge of a band, the longitudes are cut elsewhere), `None`
    without an outline.
    """
    corner_far = math.sqrt(max(abs(lo[0]), abs(hi[0])) ** 2 + max(abs(lo[1]), abs(hi[1])) ** 2
                           + max(abs(lo[2]), abs(hi[2])) ** 2)
    near = math.sqrt(sum(min(max(0.0, lo[a]), hi[a]) ** 2 for a in range(3)))
    shell_lo = spherical.shell_index_at(near, edge_pc)
    shell_hi = spherical.shell_index_at(corner_far, edge_pc)
    theta_north, theta_south = _box_colatitude_range(lo, hi)
    extents = {
        row["shell_index"]: (row["high_band_index"], row["cell_count"])
        for row in conn.execute(
            "SELECT shell_index, high_band_index, cell_count FROM galaxy_shell WHERE shell_index BETWEEN ? AND ?",
            (shell_lo, shell_hi),
        ).fetchall()
    }
    has_outline = bool(extents) or conn.execute("SELECT 1 FROM galaxy_shell LIMIT 1").fetchone() is not None
    total_cells = int(conn.execute("SELECT COALESCE(SUM(cell_count), 0) AS n FROM galaxy_shell").fetchone()["n"])
    bands, cells = [], 0.0
    for shell in range(shell_lo, shell_hi + 1):
        if has_outline and shell not in extents:
            continue
        if shell == 0:
            bands.append((0, 0, 0))
            cells += 1
            continue
        inner, outer = spherical.shell_bounds_pc(shell, edge_pc)
        cos_north = 1.0 if hi[2] > 0 and inner <= 0 else (min(1.0, hi[2] / inner) if hi[2] > 0 else hi[2] / outer)
        cos_south = -1.0 if lo[2] < 0 and inner <= 0 else (max(-1.0, lo[2] / inner) if lo[2] < 0 else lo[2] / outer)
        # A block's box can lie wholly above or below a thin shell's reach in z (the ratio then passes 1).
        north = max(theta_north, math.acos(max(-1.0, min(1.0, cos_north)))) - 1e-9
        south = min(theta_south, math.acos(max(-1.0, min(1.0, cos_south)))) + 1e-9
        if north > south:
            continue
        low = spherical.band_at(shell, south)
        high = spherical.band_at(shell, north)
        if has_outline:
            reach = extents[shell][0]
            low, high = max(low, -reach), min(high, reach)
            if low > high:
                continue
        bands.append((shell, low, high))
        if has_outline:
            cells += extents[shell][1] * (math.cos(north) - math.cos(south)) / 2
    return bands, (min(1.0, cells / total_cells) if has_outline and total_cells else None)


def bright_star_scatter_status(conn):
    """
    Whether the galaxy's bright-star scatter (`planetgen plan`) has run,
    and with what threshold and seed -- so the Generate page can say so
    and offer the right next step.

    Returns:
        dict: `scattered` (bool), `min_luminosity_sol` and `seed` (both
            `None` until a scatter runs), `default_min_luminosity_sol`
            (`BRIGHT_STAR_MIN_LUMINOSITY_SOL`, what a plain plan uses) and
            `limits` (UX.116): the settings stored with the plan, `None`
            without one: `core_min_luminosity_sol` and `core_mass_limit_sol`
            (what the core takes: the one set, else the galaxy-wide one times
            `CORE_TIER_FACTOR`), `core_follows` (neither was set),
            `halo_min_luminosity_sol`, `cell_star_ceiling` and `fill_luminosity_cap_sol`.
    """
    row = conn.execute(
        "SELECT bright_star_min_luminosity_sol, bright_star_seed, bright_star_mass_limit_sol, core_min_luminosity_sol,"
        " core_mass_limit_sol, halo_min_luminosity_sol, cell_star_ceiling, fill_luminosity_cap_sol"
        " FROM galaxy_shape WHERE id = 1").fetchone()
    scattered = row is not None and row["bright_star_min_luminosity_sol"] is not None
    limits = None
    if row is not None:
        level = float(row["bright_star_min_luminosity_sol"]) if scattered else tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL
        mass = (float(row["bright_star_mass_limit_sol"]) if row["bright_star_mass_limit_sol"] is not None
                else tuning.PHENOMENON_MIN_MASS_SOLAR)
        core_floor, core_mass = row["core_min_luminosity_sol"], row["core_mass_limit_sol"]
        limits = {
            "core_min_luminosity_sol": float(core_floor) if core_floor is not None else level * tuning.CORE_TIER_FACTOR,
            "core_mass_limit_sol": float(core_mass) if core_mass is not None else mass * tuning.CORE_TIER_FACTOR,
            "core_follows": core_floor is None and core_mass is None,
            "halo_min_luminosity_sol": float(row["halo_min_luminosity_sol"]),
            "cell_star_ceiling": float(row["cell_star_ceiling"]),
            "fill_luminosity_cap_sol": float(row["fill_luminosity_cap_sol"]),
        }
    return {
        "scattered": scattered,
        "min_luminosity_sol": float(row["bright_star_min_luminosity_sol"]) if scattered else None,
        "seed": row["bright_star_seed"] if scattered else None,
        "default_min_luminosity_sol": tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL,
        "limits": limits,
    }


def _radius_sol(radius_km):
    """A star's radius in solar radii, from kilometres (`None` stays `None`)."""
    return None if radius_km is None else radius_km * 1000.0 / constants.SOLAR_RADIUS_M


def _bright_star_entry(row):
    """One `bright_stars` row as the dict `galaxy_bright_stars_in_box` and
    `bright_stars_in_sector` return."""
    return {
        "id": row["id"],
        "x": row["position_x_mpc"] / MPC_PER_PC, "y": row["position_y_mpc"] / MPC_PER_PC,
        "z": row["position_z_mpc"] / MPC_PER_PC,
        "luminosity_sol": row["luminosity_w"] / constants.SOLAR_LUMINOSITY,
        "temperature_k": row["temperature_k"], "radius_sol": _radius_sol(row["radius_km"]),
        "star_type": row["star_type"], "population": row["population"],
        "tier": LOWEST_PRIORITY_TIER.get(row["population"], 0),
        "yerkes_class": row["yerkes_class"], "shell_index": row["shell_index"],
        "band_index": row["band_index"], "wedge_index": row["wedge_index"],
        "system_id": row["star_system_id"],
        **_member_fields(row),
    }


def scatter_member_tree(conn, table, row_id):
    """
    The stored members of an unbuilt scatter row, as a page lists them (UX.114): `members` (`member_no`, `kind`,
    `state`, `mass_solar`, `parent_node_no`, `flags`, `decay_age_yr`, which only a member of a young group has, GEN.327, and `radius_au`, GEN.952) in member
    order and `nodes` (`node_no`, `parent_node_no`, `kind`
    (`close`, `wide`, `compact` or `group`), `mass_solar`, `separation_au`, `period_years`, `eccentricity`, `interaction`, `flags`),
    both empty for a single object, and `compact_pairs` (UX.186: the X-ray binaries and spider pulsars among them,
    `db.compact_pairs`). Reads the two child tables by their primary key, so it is cheap at any galaxy size.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        table (str): `"bright_stars"` or `"phenomenon_scatter"`.
        row_id (int): The row's id.
    """
    source = SCATTER_SOURCES[table]
    members = [
        {"member_no": r["member_no"], "kind": r["kind"], "state": r["state"], "mass_solar": r["mass_solar"],
         "parent_node_no": r["parent_node_no"], "decay_age_yr": r["decay_age_yr"], "flags": r["flags"],
         "radius_au": r["radius_au"]}
        for r in conn.execute(
            "SELECT member_no, kind, state, mass_solar, parent_node_no, decay_age_yr, flags, radius_au FROM scatter_members"
            " WHERE source = ? AND source_id = ? ORDER BY member_no", (source, int(row_id))).fetchall()]
    nodes = [
        {"node_no": r["node_no"], "parent_node_no": r["parent_node_no"], "kind": r["kind"],
         "mass_solar": r["mass_solar"], "separation_au": _m_to_au(_km_to_m(r["separation_km"])),
         "period_years": r["period_years"], "eccentricity": r["eccentricity"], "interaction": r["interaction"],
         "flags": r["flags"]}
        for r in conn.execute(
            "SELECT node_no, parent_node_no, kind, mass_solar, separation_km, period_years, eccentricity, interaction, flags"
            " FROM scatter_nodes WHERE source = ? AND source_id = ? ORDER BY node_no",
            (source, int(row_id))).fetchall()]
    return {"members": members, "nodes": nodes, "compact_pairs": compact_pairs.compact_pairs(members, nodes)}


def _member_fields(row):
    """
    `member_count` (0 for a single object) and `system_mass_solar` (the mass the point carries, `None` for a single
    object) of a scatter row (UX.114). Only the readers that ask for the two columns carry them: the map nodes, which
    read thousands of rows, leave them out and the keys are absent.
    """
    if "member_count" not in row.keys():
        return {}
    return {"member_count": int(row["member_count"] or 0), "system_mass_solar": row["system_mass_solar"]}


def bright_stars_in_sector(conn, shell_index, band_index, wedge_index, unfilled_only=True):
    """
    The pre-placed bright stars (`bright_stars`) in one sector cell, most
    luminous first -- the stars a sector page lists for a cell that hasn't
    been filled yet (filling builds each into a system, see
    `planetgen fill`). Reads `idx_bright_stars_address`, so it's cheap
    at any galaxy size.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        shell_index, band_index, wedge_index (int): The cell's address.
        unfilled_only (bool): Leave out stars already built into a system
            (the default); `False` lists every star placed in the cell,
            each filled one with its `system_id`.

    Returns:
        list[dict]: As `galaxy_bright_stars_in_box`; empty when no scatter
            has run or the cell holds none.
    """
    unfilled = " AND star_system_id IS NULL" if unfilled_only else ""
    rows = conn.execute(
        f"""
        SELECT id, position_x_mpc, position_y_mpc, position_z_mpc, luminosity_w, temperature_k, radius_km,
               star_type, population, yerkes_class, shell_index, band_index, wedge_index, star_system_id,
               member_count, system_mass_solar
        FROM bright_stars
        WHERE shell_index = ? AND band_index = ? AND wedge_index = ?{unfilled}
        ORDER BY luminosity_w DESC, id
        """,
        (int(shell_index), int(band_index), int(wedge_index)),
    ).fetchall()
    return [_bright_star_entry(row) for row in rows]


UNCHARTED_SYSTEM_SORTS = {
    "luminosity": "luminosity_w", "temperature": "temperature_k", "type": "star_type",
    "sector": "shell_index,band_index,wedge_index", "age": "age_gy",
}
"""dict: The sort keys `list_uncharted_systems` accepts -> the `bright_stars` column they order by."""


_UNCHARTED_START_COLUMNS = {"luminosity": "luminosity_w", "temperature": "temperature_k", "age": "age_gy",
                            "type": "star_type", "sector": "shell_index"}
"""dict: The Uncharted stars sorts that have a jump-by-value start (PERF.135) -> the column it compares: the first
column the sort orders by."""


def parse_uncharted_start(sort, text):
    """
    The jump-by-value start of the Uncharted stars table as the value `bright_stars` compares (PERF.135): a luminosity
    in solar luminosities (as the column shows it, converted to watts), a temperature in kelvin, an age in Gyr, a
    shell index (a whole number) or the start of a star type.

    Returns:
        float, int, str or None: The value, or `None` for an empty text or a sort with no start.

    Raises:
        ValueError: If the text is not a finite number of at least zero (a whole one for the sector sort).
    """
    text = (text or "").strip()
    if not text or sort not in _UNCHARTED_START_COLUMNS:
        return None
    if sort == "type":
        return text
    try:
        number = float(text)
    except ValueError:
        raise ValueError(f"{text!r} is not a number")
    if sort == "luminosity":
        number *= constants.SOLAR_LUMINOSITY
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"{text!r} is not a number from 0 up")
    if sort == "sector":
        if number != int(number):
            raise ValueError(f"{text!r} is not a whole shell number")
        return int(number)
    return number


def _uncharted_start_condition(sort, descending, start):
    """`(condition_sql, params)` keeping the Uncharted stars from `start` on in `sort` and direction (the first row at or
    after the value; a type start keeps the types that begin with it first), or `None` for no start."""
    value = parse_uncharted_start(sort, start)
    if value is None:
        return None
    column = _UNCHARTED_START_COLUMNS[sort]
    if sort == "type" and descending:
        bound = _after_prefix(value)
        return None if bound is None else (f"{column} < ?", [bound])
    return f"{column} {'<=' if descending else '>='} ?", [value]


def count_uncharted_systems(conn, start=None, sort="luminosity", descending=True):
    """
    How many scattered stars are not yet built into a system (UX.87): the `bright_stars` rows with no
    `star_system_id`; with `start` (the jump-by-value box, PERF.135) only those from that value on in `sort`.

    Returns:
        int: The count. Kept by `countcache` (PERF.112): the scan of the whole scatter table is made once under a
            galaxy content stamp and read from the stored answer after that (a star built into a system, a
            rescatter and a reset each move the stamp, so the next read makes a new one in the background); a
            capped count (`CappedCount`) before the first.
    """
    begun = _uncharted_start_condition(sort, descending, start)
    where = "star_system_id IS NULL" + (f" AND {begun[0]}" if begun else "")
    params = begun[1] if begun else []

    def exact(c):
        return c.execute(f"SELECT COUNT(*) AS n FROM bright_stars WHERE {where}", params).fetchone()["n"]

    def estimate(c):
        return _capped(c.execute(f"SELECT COUNT(*) AS n FROM (SELECT 1 FROM bright_stars WHERE {where} LIMIT ?) c",
                                 list(params) + [COUNT_FALLBACK_CAP]).fetchone()["n"])

    return _counted(conn, ["count_uncharted_systems"] + ([sort, descending, start.strip()] if begun else []),
                    exact, estimate)


def list_uncharted_systems(conn, limit=None, offset=None, sort="luminosity", descending=True, after=None, start=None,
                           last=False):
    """
    The scattered stars no system was built around yet (UX.87), a page of
    them: each is an "uncharted" system that can be generated by itself.
    Needs no sector: it reads `bright_stars`, which holds every star the
    mass-limit and luminosity scatters placed.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        limit, offset (int, optional): The page (`offset` ignored without `limit`).
        sort (str): A key of `UNCHARTED_SYSTEM_SORTS`; ties fall back to the id, in the same direction.
        descending (bool): Reverse the sort (brightest first by default).
        after (int, optional): A star's id: the page starts after that star in this sort (a keyset page, PERF.75:
            the rows before it are not read and `offset` is ignored). A star that is gone uses `offset`.
        start (str, optional): The jump-by-value box (PERF.135): the list begins at the first star at or after this
            value in the sort and its direction (`parse_uncharted_start` says what the value is); `offset` and
            `after` count from there. A text that is not a value raises `ValueError`.
        last (bool): The last `limit` stars of the sort instead of the first (where a start lies past the end).

    Returns:
        list[dict]: One row per star: `id` (the `bright_stars` id), the
            star (`star_type`, `yerkes_class`, `population`, `mass_solar`,
            `radius_solar`, `temperature_k`, `luminosity_sol`, `age_gy`),
            its sector (`shell_index`, `band_index`, `wedge_index`,
            `designation`), `x`/`y`/`z` (galaxy frame, parsecs),
            `local_x`/`local_y`/`local_z` (parsecs from the sector's
            center, where the star's sector is in the grid),
            `galactic_radius_pc`, and `member_count` and
            `system_mass_solar` (`_member_fields`, UX.114).
    """
    if sort not in UNCHARTED_SYSTEM_SORTS:
        raise ValueError(f"unknown uncharted system sort {sort!r}")
    if last and limit is not None:
        return list(reversed(list_uncharted_systems(conn, limit=limit, sort=sort, descending=not descending)))
    direction = "DESC" if descending else "ASC"
    columns = UNCHARTED_SYSTEM_SORTS[sort].split(",")
    order = ", ".join(f"{column} {direction}" for column in columns)
    where, params = "star_system_id IS NULL", []
    begun = _uncharted_start_condition(sort, descending, start)
    if begun:
        where += f" AND {begun[0]}"
        params = list(begun[1])
    if after is not None and limit is not None:
        previous = conn.execute(f"SELECT {', '.join(columns)} FROM bright_stars WHERE id = ?", (after,)).fetchone()
        if previous is not None:
            condition, keyset_params = _keyset_condition(columns, [previous[column] for column in columns],
                                                         descending, "id", after)
            where += f" AND {condition}"
            params = params + keyset_params
            offset = 0
    query = (f"SELECT id, shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc, position_z_mpc,"
             f" population, star_type, yerkes_class, mass_kg, radius_km, temperature_k, luminosity_w, age_gy,"
             f" member_count, system_mass_solar FROM bright_stars WHERE {where} ORDER BY {order}, id {direction}")
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params = params + [int(limit), int(offset or 0)]
    skeleton = get_galaxy_shape(conn)
    edge_pc = skeleton.edge_pc if skeleton else float(tuning.DEFAULT_SECTOR_EDGE_PC)
    return [_uncharted_item(row, edge_pc) for row in conn.execute(query, params).fetchall()]


def _uncharted_item(row, edge_pc):
    """One `bright_stars` row as `list_uncharted_systems` returns it."""
    address = (row["shell_index"], row["band_index"], row["wedge_index"])
    x, y, z = (row[f"position_{axis}_mpc"] / MPC_PER_PC for axis in "xyz")
    local = galaxy_to_local_pc(sector_position_pc(*address, edge_pc), (x, y, z))
    return {
        "id": row["id"], "star_type": row["star_type"], "yerkes_class": row["yerkes_class"],
        "population": row["population"], "mass_solar": _kg_to_msun(row["mass_kg"]),
        "radius_solar": _radius_sol(row["radius_km"]), "temperature_k": row["temperature_k"],
        "luminosity_sol": row["luminosity_w"] / constants.SOLAR_LUMINOSITY, "age_gy": row["age_gy"],
        "shell_index": address[0], "band_index": address[1], "wedge_index": address[2],
        "designation": provisional_sector_designation(*address),
        "x": x, "y": y, "z": z, "local_x": local[0], "local_y": local[1], "local_z": local[2],
        "galactic_radius_pc": math.sqrt(x * x + y * y + z * z), **_member_fields(row),
    }


def uncharted_system_detail(conn, bright_star_id):
    """
    One scattered star as its page needs it (UX.114): `list_uncharted_systems`'s item, plus `system_id` (the system
    built around it, or `None` while it is uncharted), `rotation_period_hours` (GEN.846; `None` on a row made before v185), `hosted_nebula` (the nebula it lights, as
    `hostednebula.summary`, or `None`: the one its sector builds, GEN.530) and, for a star with company, its `members` and `nodes`
    (`scatter_member_tree`). A star already built into a system lists no members: its members are the system's stars.

    Raises:
        ValueError: If no such `bright_stars` row exists.
    """
    row = conn.execute(
        "SELECT id, shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc, position_z_mpc,"
        " population, star_type, yerkes_class, mass_kg, radius_km, temperature_k, luminosity_w, age_gy,"
        " member_count, system_mass_solar, star_system_id, rotation_period_hours FROM bright_stars WHERE id = ?",
        (int(bright_star_id),)).fetchone()
    if row is None:
        raise ValueError(f"no bright_stars row with id {bright_star_id}")
    skeleton = get_galaxy_shape(conn)
    detail = _uncharted_item(row, skeleton.edge_pc if skeleton else float(tuning.DEFAULT_SECTOR_EDGE_PC))
    detail["system_id"] = row["star_system_id"]
    detail["rotation_period_hours"] = row["rotation_period_hours"]
    detail["hosted_nebula"] = None if row["star_system_id"] is not None else hostednebula.hosted_by_rows(
        conn, [row]).get(row["id"])
    if detail["member_count"] and row["star_system_id"] is None:
        detail.update(scatter_member_tree(conn, "bright_stars", row["id"]))
    return detail


UNCHARTED_SCATTER_LIMIT = 500
"""int: Most scattered phenomena one uncharted sector lists (MAP.162)."""


def uncharted_sector_contents(conn, shell_index, band_index, wedge_index):
    """
    What the scatters left in one sector cell that nothing has built yet
    (MAP.162): its waiting bright stars (`bright_stars_in_sector`) and its
    unbuilt scattered phenomena (`phenomenon_scatter`), so a cell with no
    generated contents can still be opened and its objects looked at.

    Returns:
        dict: `stars` (as `bright_stars_in_sector`), `scattered`, one
            `{"id", "kind", "subtype", "type", "x", "y", "z"}` per unbuilt
            scatter row (`type` as `_SCATTER_TYPES`; positions in
            galaxy-frame parsecs), at most `UNCHARTED_SCATTER_LIMIT`, and `hosted_nebulae`, the nebulae the
            waiting stars light (`hostednebula.hosted_in_cell`, GEN.530).
    """
    address = (int(shell_index), int(band_index), int(wedge_index))
    try:
        stars = bright_stars_in_sector(conn, *address)
    except pymysql.err.ProgrammingError:  # no bright_stars table yet (before v43)
        stars = []
    rows = conn.execute(
        f"""
        SELECT id, kind, subtype, mass_solar, {point_columns.LAYOUTS["phenomenon_scatter"].select("position", "velocity")},
               member_count, system_mass_solar
        FROM phenomenon_scatter
        WHERE shell_index = ? AND band_index = ? AND wedge_index = ? AND built_at IS NULL
        ORDER BY id
        LIMIT ?
        """,
        (*address, UNCHARTED_SCATTER_LIMIT),
    ).fetchall()
    scattered = [
        {"id": row["id"], "kind": row["kind"], "subtype": row["subtype"],
         "type": _SCATTER_TYPES.get(row["kind"], row["kind"]), "mass_solar": row["mass_solar"],
         "x": row["position_x_mpc"] / MPC_PER_PC, "y": row["position_y_mpc"] / MPC_PER_PC,
         "z": row["position_z_mpc"] / MPC_PER_PC,
         "speed_kms": None if row["velocity_x_kms"] is None else math.sqrt(
             row["velocity_x_kms"] ** 2 + row["velocity_y_kms"] ** 2 + row["velocity_z_kms"] ** 2),
         **_member_fields(row)}
        for row in rows
    ]
    try:
        hosted = hostednebula.hosted_in_cell(conn, address)
    except pymysql.err.ProgrammingError:
        hosted = []
    return {"stars": stars, "scattered": scattered, "hosted_nebulae": hosted}


BRIGHT_STAR_POPULATION_WEIGHT = {"young": 1.0, "intermediate": 10.0, "old": 400.0, "bulge": 400.0,
                                 "nuclear_cluster": 400.0, "nuclear_disc": 400.0}
"""dict: The map key of a bright star is its luminosity times its population's
weight (MAP.154). Ranking by luminosity alone puts the young blue stars first
and leaves no room for the old giants of the bulge and the thick disk, which
top out near 2,500 Lsun while a young star reaches a million; the weights
bring the top of each population to about the same key (400 x 2,500 =
1,000,000). One key for every node level is what makes the lists nested: a
child node lists every star its parent listed inside the child's box."""


LOWEST_PRIORITY_TIER = {"halo": 1, "satellite": 2, "satellite_young": 2, "satellite_intermediate": 2}
"""dict: The tier of the populations every map lists last (MAP.250); a population not here is tier 0. Stars of the
stellar halo (tier 1) and of the satellite galaxies (tier 2) are the lowest priority of any star budget: a list that
cannot hold every star drops the satellite galaxies' first, then the halo's, and only then a star of the disc, the
bulge, the nuclei or a globular cluster, whatever the brightness. The tier comes before `_map_key`, so the order is
still one for every node level and the lists stay nested."""


def _map_tier(star):
    """A bright star's tier (`LOWEST_PRIORITY_TIER`): 0 listed first, the larger last."""
    return LOWEST_PRIORITY_TIER.get(star["population"], 0)


def _map_key(star):
    """A bright star's rank key (`BRIGHT_STAR_POPULATION_WEIGHT`), larger is listed first."""
    return star["luminosity_sol"] * BRIGHT_STAR_POPULATION_WEIGHT.get(star["population"], 1.0)


EVEN_DROP_BAND_DEX = 0.5
"""float: The width, in powers of ten of the map key (or of the size share of a point), of the bands a list is cut
by (MAP.251). A list that cannot hold everything keeps whole bands from the top (the most massive and the brightest
are the last to go) and thins the band the cut falls in by a fixed number every object carries (the best tenth of the list, at least one object, is kept whatever the band, `EVEN_DROP_ANCHOR_SHARE`), so that band loses
the same share of its objects in every region instead of its dimmest end being truncated wherever that lies."""

EVEN_DROP_ANCHOR_SHARE = 0.1
"""float: The share of a full list, at least one object, that is the plain best of it by rank whatever the band, so
the single most massive or brightest object is never the one a thinned band leaves out."""

EVEN_DROP_POOL = 4
"""int: How many times the list's limit the database reads of each ranking before `_even_pick` cuts it, so the band
the cut falls in has more candidates than it gets slots and the cut has something to thin."""


def _thin_number(token):
    """A fixed number in [0, 1) for an object (from its id and kind): the same object gets the same number in every
    node at every level, so a list that grows keeps every object a shorter one had (the lists stay nested)."""
    mixed = (zlib.crc32(str(token).encode("utf-8")) * 2654435761) & 0xFFFFFFFF
    return mixed / 4294967296.0


def _even_pick(items, limit, order, band, token, tier=lambda item: 0):
    """
    At most `limit` of `items`, cut evenly (MAP.251): the lowest `tier` first (the halo and the satellite galaxies
    are dropped before anything else), inside a tier the highest `band` first (the most massive and the brightest
    go last), and inside the band the cut falls in the lowest `_thin_number(token(item))` first, so that band loses
    the same share of its objects everywhere; the first `EVEN_DROP_ANCHOR_SHARE` of the slots, at least one, go to
    the best by `order` whatever the band. Returned in (`tier`, `order`) order, `order` ascending. A pure
    function of `items`, and one order for every node level, so the lists stay nested wherever the database hands
    over the same candidates.
    """
    if len(items) <= limit:
        return sorted(items, key=lambda item: (tier(item), order(item)))
    ranked = sorted(items, key=lambda item: (tier(item), order(item)))
    anchors = max(1, int(math.ceil(limit * EVEN_DROP_ANCHOR_SHARE))) if limit > 0 else 0
    kept = ranked[:anchors]
    kept += sorted(ranked[anchors:], key=lambda item: (tier(item), -band(item), _thin_number(token(item))))[:limit - anchors]
    return sorted(kept, key=lambda item: (tier(item), order(item)))


def _star_order(star):
    return (-_map_key(star), -star["id"])


def _decade_band(value):
    return int(math.floor(math.log10(max(value, 1e-30)) / EVEN_DROP_BAND_DEX))


def _top_by_key(stars, limit, even=True):
    """
    At most `limit` of `stars`: the lowest `_map_tier` first (the halo and the satellite galaxies come last), then
    the highest `_map_key` first (ties by descending id). With `even` (the nodes') the band of the key the cut falls
    in is thinned evenly (`_even_pick`, MAP.251); without it the list is the plain best `limit` (the galaxy-wide
    sample's). The order of the list is the same at every node level.
    """
    if not even:
        return sorted(stars, key=lambda star: (_map_tier(star), _star_order(star)))[:limit]
    return _even_pick(stars, limit, _star_order, lambda star: _decade_band(_map_key(star)),
                      lambda star: star["id"], _map_tier)


def galaxy_brightest_stars(conn, count=GALAXY_TILE_BRIGHTEST_SAMPLE):
    """
    The galaxy's bright stars with the highest map key
    (`BRIGHT_STAR_POPULATION_WEIGHT`), highest first (ties by descending
    id) -- `db.nodes.galaxy_nodes`' sample for nodes too big to query on their own.
    One walk down `idx_bright_stars_population` per population and height:
    the `count` most luminous of each, ranked together by the weighted key and
    cut to `count` (the old giants never reach the young stars' luminosity, so
    unweighted the sample held only young stars).

    Returns:
        list[dict]: As `galaxy_bright_stars_in_box`.
    """
    picks = []
    for population in BRIGHT_STAR_POPULATIONS:
        for off_plane in (0, 1):
            picks += [
                _bright_star_entry(row) for row in conn.execute(
                    """
                    SELECT id, position_x_mpc, position_y_mpc, position_z_mpc, luminosity_w, temperature_k,
                           radius_km, star_type, population, yerkes_class, shell_index, band_index,
                           wedge_index, star_system_id
                    FROM bright_stars FORCE INDEX (idx_bright_stars_population)
                    WHERE population = ? AND off_plane = ?
                    ORDER BY luminosity_w DESC, id DESC
                    LIMIT ?
                    """,
                    (population, off_plane, int(count)),
                ).fetchall()
            ]
    return _top_by_key(picks, int(count), even=False)


def _luminosity_walk_is_cheaper(total, cell_share, angle_share, limit):
    """
    Whether `galaxy_bright_stars_in_box` reads the box's stars by walking the luminosity index (else by address) -- the
    two answer the same, so this only chooses the cheaper (PERF.157). The address path reads the box's share of the
    stars; the walk reads, down each of its `2 x populations` arms, until `limit x EVEN_DROP_POOL` stars of the
    arm are inside the box, so about that many over the box's share. The share is the box's share of the outline's
    cells and longitudes, which counts its height above the plane, where the share of the disc's area (which sets
    whether the box is answered from the sample) does not: a 1 kpc box a kiloparsec up holds a sliver of the stars.
    """
    share = max((1.0 if cell_share is None else cell_share) * angle_share, 1e-12)
    walks = 2 * len(BRIGHT_STAR_POPULATIONS)
    return limit * EVEN_DROP_POOL * walks / share < share * total


def galaxy_bright_stars_in_box(conn, lo, hi, edge_pc, limit=None, unfilled_only=False,
                               brightest=None):
    """
    The pre-placed bright stars (`bright_stars`) in the box `[lo, hi)` with
    the highest map key (`BRIGHT_STAR_POPULATION_WEIGHT`), at most `limit` --
    the stars the Galaxy Map draws before (and after) their sectors are
    filled. The key is the same for every box, so the list of a box inside a
    bigger one holds all of the bigger one's stars inside it whenever its
    `limit` is no smaller (MAP.154).

    `bright_stars` is indexed by address and by luminosity, not by
    position, so the box is turned into address ranges: one per ring and
    layer it reaches (each with the ring's slots in its longitudes) for a
    small box, one per ring when that would be too many ranges. A box
    holding a good share of the galaxy instead walks the luminosity index
    from the top, which finds `limit` stars inside it quickly; the choice
    is by the estimated rows each way reads.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        lo (tuple): `(x, y, z)` inclusive lower corner, parsecs.
        hi (tuple): `(x, y, z)` exclusive upper corner, parsecs.
        edge_pc (float): The sector edge, parsecs.
        limit (int or None): The most stars; `None` is `GALAXY_TILE_MAX_DETAIL_BRIGHT_STARS`.
        unfilled_only (bool): Only the stars whose sector hasn't been
            filled yet (`system_id` is `None`) -- what a page listing the
            stars still waiting in a box shows.
        brightest (callable or None): Returns the galaxy-wide sample
            (`galaxy_brightest_stars`). When given and the box would read
            more than `GALAXY_TILE_BRIGHT_STAR_ROW_BUDGET` rows either way,
            the answer is the sample's stars inside the box instead: the
            same stars whenever the sample holds `limit` of them there,
            else only those.

    Returns:
        list[dict]: Highest key first (ties by descending id): `id`, `x`/`y`/`z` (parsecs),
            `luminosity_sol`, `temperature_k`, `radius_sol`, `star_type`,
            `yerkes_class`, `shell_index`, `band_index`,
            `wedge_index` and `system_id` (the star's system once its
            sector is filled, else `None`).
    """
    if limit is None:
        limit = GALAXY_TILE_MAX_DETAIL_BRIGHT_STARS
    total = int(conn.execute("SELECT COALESCE(MAX(id), 0) AS n FROM bright_stars").fetchone()["n"])
    if not total:
        return []
    bands, cell_share = _bright_star_bands(conn, lo, hi, edge_pc)
    if not bands:
        return []
    box = (
        "position_x_mpc >= ? AND position_x_mpc < ? AND position_y_mpc >= ? AND position_y_mpc < ? "
        "AND position_z_mpc >= ? AND position_z_mpc < ?"
    )
    if unfilled_only:
        box += " AND star_system_id IS NULL"
    box_params = [int(math.ceil(v * MPC_PER_PC)) for pair in zip(lo, hi) for v in pair]

    # Rows read each way, assuming stars spread evenly over the sectors: the
    # address path reads the box's longitudes of the bands it reaches, the
    # luminosity path about `limit / share of the galaxy`.
    intervals = _box_angle_intervals(lo, hi)
    angle_share = min(1.0, sum(end - start for start, end in intervals) / (2 * math.pi))
    outer = int(conn.execute("SELECT COALESCE(MAX(shell_index), 0) AS r FROM galaxy_shell").fetchone()["r"]) or bands[-1][0]
    disk = math.pi * ((outer + 1) * edge_pc) ** 2
    box_share = min(1.0, (hi[0] - lo[0]) * (hi[1] - lo[1]) / disk)
    # The address index is (shell, band, wedge), so only the box's own
    # longitudes of each band are read.
    by_address = (1.0 if cell_share is None else cell_share) * angle_share * total
    by_luminosity = limit / max(box_share, 1e-12)

    if brightest is not None and min(by_address, by_luminosity) > GALAXY_TILE_BRIGHT_STAR_ROW_BUDGET:
        found = []
        for star in brightest():
            if unfilled_only and star["system_id"] is not None:
                continue
            # The stored positions are whole micro-parsecs; test them as
            # the query does.
            if all(math.ceil(lo[a] * MPC_PER_PC) <= round(star[k] * MPC_PER_PC) < math.ceil(hi[a] * MPC_PER_PC)
                   for a, k in enumerate(("x", "y", "z"))):
                found.append(star)
        return _top_by_key(found, limit)

    columns = ("id, position_x_mpc, position_y_mpc, position_z_mpc, luminosity_w, temperature_k, radius_km, "
               "star_type, population, yerkes_class, shell_index, band_index, wedge_index, star_system_id")
    # Both keys descending, so the luminosity path walks its index
    # backwards and stops at `limit`: with `id` ascending MariaDB and MySQL
    # sort every star in the box first, millions in a zoomed-out node.
    order = "luminosity_w DESC, id DESC"
    if _luminosity_walk_is_cheaper(total, cell_share, angle_share, limit):
        # One walk per population and height down (population, off_plane,
        # luminosity), each stopping at `limit`: a single walk down the
        # luminosity index would meet only young stars (the weights of
        # `BRIGHT_STAR_POPULATION_WEIGHT` apply once they are merged).
        selects, params = [], []
        for population in BRIGHT_STAR_POPULATIONS:
            for off_plane in (0, 1):
                selects.append(
                    f"(SELECT {columns} FROM bright_stars FORCE INDEX (idx_bright_stars_population) "
                    f"WHERE population = ? AND off_plane = ? AND {box} ORDER BY {order} LIMIT ?)")
                params += [population, off_plane] + box_params + [int(limit) * EVEN_DROP_POOL]
        rows = conn.execute(" UNION ALL ".join(selects), params).fetchall()
        return _top_by_key([_bright_star_entry(row) for row in rows], limit)
    clauses, params = [], []
    exact = sum(t - b + 1 for _s, b, t in bands) <= BRIGHT_STAR_MAX_EXACT_RANGES
    for shell, bottom, top in bands:
        if exact:
            for band in range(bottom, top + 1):
                wedges = _band_wedge_ranges(shell, band, intervals)
                wedge_test = " OR ".join("wedge_index BETWEEN ? AND ?" for _ in wedges)
                clauses.append(f"(shell_index = ? AND band_index = ? AND ({wedge_test}))")
                params += [shell, band] + [v for pair in wedges for v in pair]
        else:
            clauses.append("(shell_index = ? AND band_index BETWEEN ? AND ?)")
            params += [shell, bottom, top]
    # The address path reads the box's stars once and keeps the `limit`
    # most luminous of each population, then the `limit` of highest key.
    rows = conn.execute(
        f"""
        SELECT * FROM (
            SELECT {columns},
                   ROW_NUMBER() OVER (PARTITION BY population ORDER BY {order}) AS rank_in_population
            FROM bright_stars FORCE INDEX (idx_bright_stars_address)
            WHERE ({" OR ".join(clauses)}) AND {box}
        ) ranked
        WHERE rank_in_population <= ?
        """,
        params + box_params + [int(limit) * EVEN_DROP_POOL],
    ).fetchall()
    return _top_by_key([_bright_star_entry(row) for row in rows], limit)


GALAXY_TILE_MAX_DETAIL_STARS = 1200
"""int: Most stars of generated systems a finest node (`TILE_MAX_LEVEL`,
16 pc) lists (MAP.80): the map fetches these around a sector it is
zoomed to, so the sector shows every star it holds (a 16 pc node of a
filled region holds a few hundred). Only a node deep in the bulge holds
more, and then its faintest are left out."""

GALAXY_TILE_MAX_DETAIL_BRIGHT_STARS = 300
"""int: Most pre-placed bright stars a finest node lists."""

GALAXY_TILE_GENERATED_STAR_SECTOR_BUDGET = 1500
"""int: Most generated sectors one node reads stars from. A coarse node
over a large filled region reads an even sample of its sectors (every
k-th, by a hash of the id) instead, so the request stays bounded; zooming in reaches
nodes small enough to read all of them."""


SECTOR_ALLOWANCE_FUDGE = 3
"""int: A node's budget is shared out by sector: each generated sector may
put `SECTOR_ALLOWANCE_FUDGE * budget / sectors` stars (at least one) on
the node, its brightest, so one dense sector can't take the whole budget
while a sparse sector, which has fewer stars than that, keeps all it has
and the leftover room goes to the brighter stars of the dense ones."""

def galaxy_generated_stars_in_box(conn, lo, hi, min_luminosity_sol, sector_count,
                                  limit=GALAXY_TILE_MAX_DETAIL_STARS, per_sector=None):
    """
    The most luminous stars of generated systems whose sector's center
    lies in the box `[lo, hi)`, at least `min_luminosity_sol` each, at
    most `limit` -- what the Galaxy Map draws as points of light once
    sectors are filled (MAP.51). A star also listed as a pre-placed
    bright star (`bright_stars.star_system_id`) is left out, so it isn't
    drawn twice; a bright star's companion is still listed.

    The query walks the sectors in the box by `idx_sectors_center`, then
    their systems and stars by their `sector_id`/`star_system_id`
    indexes, so it reads at most `GALAXY_TILE_GENERATED_STAR_SECTOR_BUDGET`
    sectors' stars: when the box holds more, about one in k of them,
    picked by a hash of the id so the sample has no stripes.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        lo (tuple): `(x, y, z)` inclusive lower corner, parsecs.
        hi (tuple): `(x, y, z)` exclusive upper corner, parsecs.
        min_luminosity_sol (float): The faintest star listed.
        sector_count (int): How many generated (grid-placed) sectors the
            box holds, as `galaxy_sector_count_in_box` counts them.
        limit (int): See `GALAXY_TILE_STAR_BUDGET`.
        per_sector (int or None): Most stars one sector contributes, its
            brightest (`generated_star_budget`); `None` is no limit.

    Returns:
        list[dict]: Most luminous first (ties by descending id): `id`
            (the `stars` row), `name`, `x`/`y`/`z` (galaxy-frame
            parsecs), `luminosity_sol`, `temperature_k`, `radius_sol`,
            `star_type`, `shell_index`, `band_index`, `wedge_index`
            and `system_id`.
    """
    if sector_count <= 0:
        return []
    stride = max(1, int(math.ceil(sector_count / float(GALAXY_TILE_GENERATED_STAR_SECTOR_BUDGET))))
    select = """
        SELECT st.id, st.name, st.role, st.luminosity_w, st.temperature_k, st.radius_km, st.star_type,
               ss.id AS system_id, ss.binary_separation_km, ss.position_x_mpc, ss.position_y_mpc, ss.position_z_mpc,
               sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.shell_index, sec.band_index, sec.wedge_index"""
    rank = (", ROW_NUMBER() OVER (PARTITION BY sec.id ORDER BY st.luminosity_w DESC, st.id DESC) AS sector_rank"
            if per_sector is not None else "")
    source = """
        FROM sectors sec FORCE INDEX (idx_sectors_center)
        JOIN star_systems ss ON ss.sector_id = sec.id
        JOIN stars st ON st.star_system_id = ss.id
        LEFT JOIN bright_stars b ON b.star_system_id = ss.id AND st.role NOT IN ('secondary', 'member')
        WHERE sec.center_x_pc >= ? AND sec.center_x_pc < ?
          AND sec.center_y_pc >= ? AND sec.center_y_pc < ?
          AND sec.center_z_pc >= ? AND sec.center_z_pc < ?
          AND sec.shell_index IS NOT NULL AND MOD(CRC32(sec.id), ?) = 0
          AND ss.position_x_mpc IS NOT NULL
          AND st.luminosity_w >= ? AND b.id IS NULL"""
    params = [lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], stride, float(min_luminosity_sol) * constants.SOLAR_LUMINOSITY]
    if per_sector is None:
        sql = select + source + "\n        ORDER BY st.luminosity_w DESC, st.id DESC\n        LIMIT ?"
    else:
        sql = ("SELECT * FROM (" + select + rank + source + ") ranked WHERE sector_rank <= ?"
               "\n        ORDER BY luminosity_w DESC, id DESC\n        LIMIT ?")
        params.append(int(per_sector))
    rows = conn.execute(sql, params + [int(limit)]).fetchall()
    return [_generated_star_entry(row) for row in rows]


def _generated_star_entry(row):
    """One generated star's row as the dict `galaxy_generated_stars_in_box` returns."""
    # A system's stored offset is along its sector's own axes (MAP.186), so it is
    # rotated into the galaxy frame like the Sector Map does; a binary's
    # secondary sits its real separation from the primary (`binary_companion_offset_pc`).
    center = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
    offset = [row[f"position_{axis}_mpc"] / MPC_PER_PC for axis in "xyz"]
    point = local_to_galaxy_pc(center, offset)
    if row["role"] in ("secondary", "member"):
        shift = binary_companion_offset_pc(row["binary_separation_km"])
        point = tuple(point[i] + shift[i] for i in range(3))
    return {
        "id": row["id"], "name": row["name"],
        "x": round(point[0], 3), "y": round(point[1], 3), "z": round(point[2], 3),
        "luminosity_sol": float("%.4g" % (row["luminosity_w"] / constants.SOLAR_LUMINOSITY)),
        "temperature_k": round(row["temperature_k"]),
        "radius_sol": float("%.3g" % _radius_sol(row["radius_km"])),
        "star_type": row["star_type"], "shell_index": row["shell_index"], "band_index": row["band_index"],
        "wedge_index": row["wedge_index"], "system_id": row["system_id"],
    }


GALAXY_TILE_MAX_POINTS = 200
"""int: Most point-like phenomena one node lists, the most luminous
first. A 64 pc node holds a handful."""

_POINT_PHENOMENON_TABLES = (
    ("black_holes", "black_hole", "(CASE WHEN has_accretion_disk THEN 'accreting' ELSE 'quiescent' END)", "mass_solar"),
    ("neutron_stars", "neutron_star", "pulsar_type", "mass_solar"),
    ("quasars", "quasar", "(CASE WHEN is_radio_loud THEN 'radio-loud' ELSE 'radio-quiet' END)", "black_hole_mass_solar"),
)
"""tuple: `(table, type_label, descriptor_expr, mass_column)` for the phenomena the
Galaxy Map draws as points (the same labels and descriptors as
`_PHENOMENON_TABLES`; a quasar's mass is its black hole's)."""


def galaxy_point_phenomena_in_box(conn, lo, hi, limit=GALAXY_TILE_MAX_POINTS):
    """
    The placed black holes, neutron stars and quasars whose center lies
    in the box `[lo, hi)` (MAP.80), at most `limit`, read by each table's
    `idx_<table>_center`. One bound to a star system (`star_id`) has no
    placement of its own and is drawn as that system's star instead.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        lo (tuple): `(x, y, z)` inclusive lower corner, parsecs.
        hi (tuple): `(x, y, z)` exclusive upper corner, parsecs.
        limit (int): See `GALAXY_TILE_MAX_POINTS`.

    Returns:
        list[dict]: Heaviest first (by `map_luminosity_sol`; ties by type, then id), the heaviest half of a full list
            and the rest spread over the box (`_even_pick`, MAP.251): `type`
            (`"black_hole"`, `"neutron_star"` or `"quasar"`), `id`,
            `name`, `descriptor`, `mass_solar`, `luminosity_sol` (its own light), `map_luminosity_sol` (the
            luminosity of a star of its mass, `phenomenon_map_luminosity_sol`:
            what the maps draw it at, MAP.208) and `x`/`y`/`z` (parsecs).
    """
    found = []
    for table, type_label, descriptor_expr, mass_column in _POINT_PHENOMENON_TABLES:
        rows = conn.execute(
            f"""
            SELECT id, name, {descriptor_expr} AS descriptor, {mass_column} AS mass_solar, luminosity_w,
                   center_x_pc, center_y_pc, center_z_pc
            FROM {table} FORCE INDEX (idx_{table}_center)
            WHERE center_x_pc >= ? AND center_x_pc < ? AND center_y_pc >= ? AND center_y_pc < ?
              AND center_z_pc >= ? AND center_z_pc < ?
            ORDER BY {mass_column} DESC, id
            LIMIT ?
            """,
            (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], int(limit) * EVEN_DROP_POOL),
        ).fetchall()
        found += [{
            "type": type_label, "id": row["id"], "name": row["name"], "descriptor": row["descriptor"],
            "mass_solar": None if row["mass_solar"] is None else float("%.4g" % row["mass_solar"]),
            "luminosity_sol": float("%.4g" % (row["luminosity_w"] / constants.SOLAR_LUMINOSITY)),
            "map_luminosity_sol": float("%.4g" % phenomenon_map_luminosity_sol(type_label, row["mass_solar"])),
            "x": round(row["center_x_pc"], 3), "y": round(row["center_y_pc"], 3), "z": round(row["center_z_pc"], 3),
        } for row in rows]
    return _even_pick(found, limit, lambda point: (-point["map_luminosity_sol"], point["type"], point["id"]),
                      lambda point: _decade_band(point["map_luminosity_sol"]),
                      lambda point: (point["type"], point["id"]))


def galaxy_system_interests_in_box(conn, lo, hi, limit=GALAXY_TILE_MAX_POINTS):
    """
    The star systems worth a mark on the map (a habitable world, a gas giant: `interests.REASONS`) whose position lies in the box
    `[lo, hi)` (MAP.305), at most `limit`, those with the most bodies behind their reasons first (ties by id), read from
    `system_interests` by its position key. The planet and life class reasons (`interests.PER_CLASS_PREFIXES`) are not read here, `galaxy_class_interests_in_box` does.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        lo (tuple): `(x, y, z)` inclusive lower corner, parsecs.
        hi (tuple): `(x, y, z)` exclusive upper corner, parsecs.
        limit (int): See `GALAXY_TILE_MAX_POINTS`.

    Returns:
        list[dict]: `id`, `name`, `x`/`y`/`z` (parsecs) and `reasons`, a dict by reason id of `n` (how many bodies) and, for a
        reason with a score, `score` (the largest) and `tier` (its PHI-6 colour); the `civilization` reason has instead `tech` (the technology
        level of its most advanced species) and `band` (its name), when the system's species are drawn (MAP.344).
    """
    rows = conn.execute(
        """
        SELECT i.star_system_id, i.reason, i.n, i.best_score, i.galaxy_x_pc, i.galaxy_y_pc, i.galaxy_z_pc
        FROM system_interests i FORCE INDEX (idx_system_interests_position)
        WHERE i.galaxy_x_pc >= ? AND i.galaxy_x_pc < ? AND i.galaxy_y_pc >= ? AND i.galaxy_y_pc < ?
          AND i.galaxy_z_pc >= ? AND i.galaxy_z_pc < ?
          AND i.reason NOT LIKE ? AND i.reason NOT LIKE ?
        ORDER BY i.n DESC, i.star_system_id, i.reason
        LIMIT ?
        """,
        (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], *(prefix + "%" for prefix in interests.PER_CLASS_PREFIXES),
         max(0, int(limit)) * len(interests.BASE_REASONS)),
    ).fetchall()
    return _interest_systems(conn, rows, limit)


def galaxy_class_interests_in_box(conn, lo, hi, per_class):
    """
    The star systems with a planet or moon of each planet class and with each life class living on one (the `class-` and `lifeclass-` reasons of
    `interests.REASONS`, MAP.220, MAP.351) whose position lies in the box `[lo, hi)`: for each class at most `per_class` systems, those with
    the most bodies of it first (ties by id). The base
    reasons are `galaxy_system_interests_in_box`'s; a system is listed here with only its class reasons.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        lo (tuple): `(x, y, z)` inclusive lower corner, parsecs.
        hi (tuple): `(x, y, z)` exclusive upper corner, parsecs.
        per_class (int): The most systems listed for one class.

    Returns:
        list[dict]: As `galaxy_system_interests_in_box`, with only class reasons (`class-` and `lifeclass-`) in `reasons`.
    """
    if int(per_class) <= 0:
        return []
    rows = conn.execute(
        """
        SELECT star_system_id, reason, n, best_score, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc FROM (
            SELECT i.star_system_id, i.reason, i.n, i.best_score, i.galaxy_x_pc, i.galaxy_y_pc, i.galaxy_z_pc,
                   ROW_NUMBER() OVER (PARTITION BY i.reason ORDER BY i.n DESC, i.star_system_id) AS rk
            FROM system_interests i FORCE INDEX (idx_system_interests_position)
            WHERE i.galaxy_x_pc >= ? AND i.galaxy_x_pc < ? AND i.galaxy_y_pc >= ? AND i.galaxy_y_pc < ?
              AND i.galaxy_z_pc >= ? AND i.galaxy_z_pc < ? AND (i.reason LIKE ? OR i.reason LIKE ?)
        ) ranked WHERE rk <= ?
        ORDER BY n DESC, star_system_id, reason
        """,
        (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], *(prefix + "%" for prefix in interests.PER_CLASS_PREFIXES), int(per_class)),
    ).fetchall()
    return _interest_systems(conn, rows, None)


def _interest_systems(conn, rows, limit):
    """The `system_interests` `rows` as the systems they belong to (`galaxy_system_interests_in_box`' shape), the `limit` with the most
    bodies kept (all of them for a `limit` of None), named."""
    systems = {}
    for row in rows:
        system = systems.setdefault(row["star_system_id"], {
            "id": row["star_system_id"], "x": round(row["galaxy_x_pc"], 3), "y": round(row["galaxy_y_pc"], 3),
            "z": round(row["galaxy_z_pc"], 3), "reasons": {}, "bodies": 0})
        reason = {"n": int(row["n"])}
        if row["best_score"] is not None and row["reason"] in interests.TECH_REASONS:
            reason["tech"] = round(float(row["best_score"]), 2)
            reason["band"] = tech.band_name(reason["tech"])
        elif row["best_score"] is not None:
            reason["score"] = round(float(row["best_score"]), 3)
            reason["tier"] = habitability.tier(float(row["best_score"]))
        system["reasons"][row["reason"]] = reason
        system["bodies"] += int(row["n"])
    chosen = sorted(systems.values(), key=lambda system: (-system["bodies"], system["id"]))
    if limit is not None:
        chosen = chosen[:max(0, int(limit))]
    names = {}
    for start in range(0, len(chosen), 500):
        ids = [system["id"] for system in chosen[start:start + 500]]
        for row in conn.execute(f"SELECT id, name FROM star_systems WHERE id IN ({', '.join('?' * len(ids))})", ids).fetchall():
            names[row["id"]] = row["name"]
    out = []
    for system in chosen:
        del system["bodies"]
        system["name"] = names.get(system["id"])
        out.append(system)
    return out


SCATTERED_POINT_CLASSES = (
    ("quasar", None, "quasar", 1.0),
    ("black-hole", "supermassive", "black_hole", 1.0),
    ("black-hole", "intermediate", "black_hole", 0.8),
    ("black-hole", "stellar", "black_hole", 0.55),
    ("neutron-star", None, "neutron_star", 0.35),
    ("neutron-star", "magnetar", "neutron_star", 0.35),
)
"""tuple: `(scatter kind, subtype, map type, size share)` for the scattered
objects the Galaxy Map draws as points before their sector is filled
(MAP.164), biggest first. The share is the size of a row with no stored
mass (one scattered before v80), from its mass class: 1 for the nucleus
black hole down to 0.35 for a neutron star. A row with a mass is sized by
`scattered_point_size`."""

SCATTERED_SIZE_MASS_DECADES = 8.0
"""float: The mass, as powers of ten of a solar mass, at which a scattered
point reaches full size (MAP.165); a solar mass is the smallest."""

SCATTERED_SIZE_MIN = 0.2
"""float: The size share of a solar mass."""


def scattered_point_size(mass_solar):
    """The 0..1 size share of a scattered black hole or neutron star of
    `mass_solar` (MAP.165): a straight line in the logarithm of the mass from
    `SCATTERED_SIZE_MIN` at one solar mass to 1 at 10^`SCATTERED_SIZE_MASS_DECADES`,
    so a heavier object is never drawn smaller and every class keeps its place
    (neutron stars, stellar, intermediate and supermassive black holes)."""
    decades = max(0.0, min(SCATTERED_SIZE_MASS_DECADES, math.log10(max(mass_solar, 1.0))))
    return round(SCATTERED_SIZE_MIN + (1.0 - SCATTERED_SIZE_MIN) * decades / SCATTERED_SIZE_MASS_DECADES, 4)

SCATTERED_SIZE_BAND = 0.1
"""float: The width of the size share bands a list of scattered points is cut by (`_even_pick`, MAP.251)."""

SCATTERED_POINT_COARSE_CLASSES = 3
"""int: Nodes coarser than `POINT_PHENOMENON_MIN_LEVEL` list only this many
of the biggest `SCATTERED_POINT_CLASSES` (the nucleus, quasar or black hole, and
intermediate-mass black holes), so the largest are visible from the whole-galaxy view."""

GALAXY_TILE_MAX_SCATTERED_POINTS = 60
"""int: Most scattered points one coarse node lists."""


SCATTERED_MAX_SHELL_RANGES = 1000
"""int: Most `(shell, bands)` index ranges one tile's scattered-point query lists before it reads the whole run of shells
in one range instead (PERF.157). The one range reads every row of those shells, bands and wedges the box does not reach
included, and the optimizer then often scans the table (3.6 million rows for a 256 pc tile at the centre); a range per
shell, with the bands the box reaches there, reads only those. It was 64."""


def galaxy_scattered_points_in_box(conn, lo, hi, edge_pc, limit, coarse):
    """The scattered black holes, neutron stars and quasars (`phenomenon_scatter`, not yet built into a sector) whose position lies in the box `[lo, hi)`, biggest class first, at most `limit` -- drawn as points before their sector is filled (MAP.164). Each point has `charted`: true once its sector is built (a coarse node lists those too), so "Charted only" dims the others (MAP.111)."""
    classes = SCATTERED_POINT_CLASSES[:SCATTERED_POINT_COARSE_CLASSES] if coarse else SCATTERED_POINT_CLASSES
    bands, _share = _bright_star_bands(conn, lo, hi, edge_pc)
    if not bands:
        return []
    if len(bands) > SCATTERED_MAX_SHELL_RANGES:
        where_address = "shell_index BETWEEN ? AND ?"
        address_params = [bands[0][0], bands[-1][0]]
    else:
        where_address = "(" + " OR ".join("(shell_index = ? AND band_index BETWEEN ? AND ?)" for _ in bands) + ")"
        address_params = [v for band in bands for v in band]
    kinds = " OR ".join("(kind = ? AND subtype <=> ?)" for _ in classes)
    kind_params = [v for cls in classes for v in cls[:2]]
    box_params = [int(math.ceil(v * MPC_PER_PC)) for pair in zip(lo, hi) for v in pair]
    rows = conn.execute(
        f"""
        SELECT id, kind, subtype, built_at, mass_solar, position_x_mpc, position_y_mpc, position_z_mpc
        FROM phenomenon_scatter
        WHERE {"" if coarse else "built_at IS NULL AND "}{where_address} AND ({kinds})
          AND position_x_mpc >= ? AND position_x_mpc < ? AND position_y_mpc >= ? AND position_y_mpc < ?
          AND position_z_mpc >= ? AND position_z_mpc < ?
        ORDER BY (kind = 'quasar') DESC, (subtype <=> 'supermassive') DESC, (subtype <=> 'intermediate') DESC, (subtype <=> 'stellar') DESC, mass_solar DESC, id
        LIMIT ?
        """,
        (*address_params, *kind_params, *box_params, int(limit) * EVEN_DROP_POOL),
    ).fetchall()
    share = {(kind, subtype): (map_type, size) for kind, subtype, map_type, size in classes}
    points = []
    for row in rows:
        map_type, size = share[(row["kind"], row["subtype"])]
        if row["mass_solar"] is not None:
            size = scattered_point_size(row["mass_solar"])
        label = ("Magnetar" if row["subtype"] == "magnetar"
                 else {"black_hole": "Black hole", "quasar": "Quasar"}.get(map_type, "Neutron star"))
        points.append({
            "type": map_type, "id": f"s{row['id']}", "name": f"{label} (uncharted)" if row["built_at"] is None else label,
            "descriptor": row["subtype"] or "scattered", "luminosity_sol": 0.0, "scattered": True, "size": size,
            "charted": row["built_at"] is not None,
            "map_luminosity_sol": float("%.4g" % phenomenon_map_luminosity_sol(map_type, row["mass_solar"])),
            "mass_solar": None if row["mass_solar"] is None else float("%.4g" % row["mass_solar"]),
            "x": round(row["position_x_mpc"] / MPC_PER_PC, 3), "y": round(row["position_y_mpc"] / MPC_PER_PC, 3),
            "z": round(row["position_z_mpc"] / MPC_PER_PC, 3),
        })
    return _even_pick(points, limit, lambda point: (-point["size"], point["id"]),
                      lambda point: int(point["size"] / SCATTERED_SIZE_BAND), lambda point: point["id"])


TILE_PIECE_SECONDS = 3.0
"""float: The statement time limit of each piece of a Galaxy Map node (PERF.76): placed sectors, sector count,
clouds, bright stars, generated stars, points and scattered points. A piece that runs past it is left out of
the node, which is served marked `incomplete`."""


def _tile_piece(conn, incomplete, name, compute, default):
    """One piece of a node under `TILE_PIECE_SECONDS`: its value, or `default` (and `name` added to
    `incomplete`) when the database stopped it for taking too long."""
    try:
        with conn.statement_limit(TILE_PIECE_SECONDS):
            return compute()
    except pymysql.err.OperationalError as exc:
        if exc.args and exc.args[0] in STATEMENT_TIMEOUT_ERRORS:
            incomplete.append(name)
            return default
        raise


def galaxy_edge_pc(conn):
    """The galaxy's sector edge in parsecs: the stored skeleton's, else the default."""
    skeleton = get_galaxy_shape(conn)
    return skeleton.edge_pc if skeleton else float(tuning.DEFAULT_SECTOR_EDGE_PC)


def galaxy_stage(conn, at=None):
    """
    One drill-down stage's generated counts (the Galaxy Map drill-down,
    docs/design/galaxy/spherical-block-ladder.md section 4): how many generated
    sectors each child block of `at` holds. Totals (allowed sectors) are the
    page's own math, so only generated counts come from here, and children
    with none are left out.

    Without `at`, the galaxy: every level-243 block holding a generated
    sector. A block above level 3 reads its children from the stored block
    totals (`galaxy_block_stats`). A level-3 `at` has sectors for children,
    read from `sectors` and `sector_stats` by their addresses, and `sectors`
    lists each generated one.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        at (str or None): A block key, `m.shell.band.wedge` (see
            `drill.parse_stage_key`), or `None` for the galaxy.

    Returns:
        dict: `at` (the canonical key, or `None`), `edge_pc` (the sector edge,
            parsecs, which a window of sectors needs to find its cells), `child_m`, `children`
            (`[{shell, band, wedge, generated, stats}]`, at `child_m = 1` a
            sector's own address), and `sectors` (`[{shell, band, wedge, id,
            name, system_count, stats}]` at `child_m = 1`, else `None`).
            `stats` is what the generated sectors hold, from `sector_stats`
            (DB.14, `_stage_stats`): `{systems, expected_systems, stars,
            mean_age_gy, luminosity_sol, habitable_worlds}`.

    Raises:
        ValueError: On a malformed or impossible `at`.
    """
    edge_pc = galaxy_edge_pc(conn)
    if at is None or at == "":
        rows = _block_rows(conn, STAGE_TOP)
        counts = {}
        for r in rows:
            _add_stage_sums(counts, DrillBlock(STAGE_TOP, int(r["shell_index"]), int(r["band_index"]), int(r["wedge_index"])), r)
        return {"at": None, "child_m": STAGE_TOP, "children": _stage_children(counts), "sectors": None, "edge_pc": edge_pc}

    block = parse_stage_key(at)
    child_m = block.m // DRILL_STEP
    children = [c for _shell, group in drill_children(block) for c in group]

    if child_m == 1:
        if not children:
            return {"at": format_drill_key(block), "child_m": 1, "children": [], "sectors": [], "edge_pc": edge_pc}
        where, params = _sectors_where(children)
        rows = conn.execute(
            f"SELECT * FROM ({_stage_sectors_with_stats(where, 'id, name, ')}) placed "
            f"ORDER BY shell_index, band_index DESC, wedge_index",
            params,
        ).fetchall()
        system_counts = {}
        if rows:
            ids = [r["id"] for r in rows]
            placeholders = ", ".join("?" for _ in ids)
            system_counts = {
                c["sector_id"]: int(c["n"]) for c in conn.execute(
                    f"SELECT sector_id, COUNT(*) AS n FROM star_systems WHERE sector_id IN ({placeholders}) "
                    f"GROUP BY sector_id",
                    ids,
                ).fetchall()
            }
        sectors = [{
            "shell": r["shell_index"], "band": r["band_index"], "wedge": r["wedge_index"],
            "id": r["id"], "name": r["name"], "system_count": system_counts.get(r["id"], 0),
            "stats": _stage_stats(_sector_sums(r)),
        } for r in rows]
        counts = {}
        for r in rows:
            _add_stage_sums(counts, DrillBlock(1, r["shell_index"], r["band_index"], r["wedge_index"]), {
                "n": 1, "systems": r["actual_systems"], "expected": r["expected_systems"], "stars": r["actual_stars"],
                "age_sum": (r["mean_age_gy"] or 0.0) * (r["actual_stars"] or 0), "luminosity": r["total_luminosity_sol"],
                "habitable": r["habitable_worlds"], "mass": r["total_mass_kg"],
            })
        return {"at": format_drill_key(block), "child_m": 1, "children": _stage_children(counts), "sectors": sectors, "edge_pc": edge_pc}

    counts = {}
    for r in _block_rows(conn, child_m, [(c.shell, c.band, c.wedge) for c in children]):
        _add_stage_sums(counts, DrillBlock(child_m, int(r["shell_index"]), int(r["band_index"]), int(r["wedge_index"])), r)
    stage = {"at": format_drill_key(block), "child_m": child_m, "children": _stage_children(counts), "sectors": None, "edge_pc": edge_pc}
    return _with_region_names(conn, stage, block, edge_pc)


def _with_region_names(conn, stage, block, edge_pc):
    """`stage` with the region names of the wire (API.35): `name` and `parent` on each child that is a charted 9-block (a child
    is listed only when it holds a generated sector), and on the stage itself when `block` is a charted 9-block. Nothing is added
    without a naming key in force."""
    from planetgen.db import region_stats

    key = naming_key.active_key()
    if key is None:
        return stage
    if stage["child_m"] == region_stats.LEVEL:
        for child in stage["children"]:
            child.update(region_stats.labels(DrillBlock(region_stats.LEVEL, child["shell"], child["band"], child["wedge"]), key, edge_pc))
    if block.m == region_stats.LEVEL and region_stats.region_row(conn, block) is not None:
        stage.update(region_stats.labels(block, key, edge_pc))
    return stage


def _sectors_where(sectors):
    """The `WHERE` on `sectors`' own columns (and its parameters) for an explicit list of sector addresses, one
    clause per run of wedges of a (shell, band): a level-3 block holds a few dozen sectors, a handful of runs."""
    runs = {}
    for sector in sorted(sectors):
        runs.setdefault((sector.shell, sector.band), []).append(sector.wedge)
    clauses, params = [], []
    for (shell, band), wedges in runs.items():
        start = previous = wedges[0]
        for wedge in wedges[1:] + [None]:
            if wedge is not None and wedge == previous + 1:
                previous = wedge
                continue
            clauses.append("(shell_index = ? AND band_index = ? AND wedge_index BETWEEN ? AND ?)")
            params += [shell, band, start, previous]
            if wedge is not None:
                start = previous = wedge
    return " OR ".join(clauses), params


_BLOCK_STATS_COLUMNS = (
    "shell_index, band_index, wedge_index, sector_count AS n, systems, expected_systems AS expected, stars, age_sum,"
    " luminosity_sol AS luminosity, habitable_worlds AS habitable, mass_sum_kg AS mass"
)


def _block_rows(conn, m, cells=None):
    """The stored `galaxy_block_stats` rows of level `m`: every one, or just those of `cells`, a list of `(shell,
    band, wedge)`."""
    if cells is None:
        return conn.execute(f"SELECT {_BLOCK_STATS_COLUMNS} FROM galaxy_block_stats WHERE m = ?", [m]).fetchall()
    cells = list(cells)
    if not cells:
        return []
    return conn.execute(
        f"SELECT {_BLOCK_STATS_COLUMNS} FROM galaxy_block_stats WHERE m = ? AND (shell_index, band_index, wedge_index) IN ("
        + ", ".join("(?, ?, ?)" for _ in cells) + ")",
        [m] + [value for cell in cells for value in cell]).fetchall()


def galaxy_block_sums(conn, block):
    """
    The sums of the generated sectors inside one block (DB.25): `[sectors, systems, expected systems, stars, stars
    times mean age, luminosity, habitable worlds, mass in kilograms]` (`_add_stage_sums`' order), or `None` for a block with no generated sector. A
    level-3 block adds up its sectors from `sectors` and `sector_stats`; a bigger one adds up the stored rows of its
    children, so a refresh climbs the ladder one level at a time and never reads a block's sectors in one go.
    """
    counts = {}
    if block.m == DRILL_STEP:
        sectors = drill_block_sectors(block)
        if not sectors:
            return None
        where, params = _sectors_where(sectors)
        row = conn.execute(f"SELECT {_STAGE_LOOK_SUMS} FROM ({_stage_sectors_with_stats(where)}) placed", params).fetchone()
        if not row or not row["n"]:
            return None
        _add_stage_sums(counts, block, row)
        return counts[block]
    child_m = block.m // DRILL_STEP
    children = [c for _shell, group in drill_children(block) for c in group]
    for r in _block_rows(conn, child_m, [(c.shell, c.band, c.wedge) for c in children]):
        _add_stage_sums(counts, block, r)
    return counts.get(block)


def galaxy_block_stats(conn, m, cells):
    """The stored aggregates (DB.25, `galaxy_block_stats`) of the cells of one ladder level a map window shows, in one statement.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        m (int): The level, one of `drill.DRILL_LEVELS` above 1.
        cells (iterable): `(shell, band, wedge)` of each cell wanted.

    Returns:
        dict: `(shell, band, wedge)` -> `{"generated", "stats"}` (`_stage_stats`) for each cell holding a generated
            sector; the others are left out, as `galaxy_stage` leaves them out.
    """
    found = {}
    for r in _block_rows(conn, m, cells):
        counts = {}
        _add_stage_sums(counts, None, r)
        found[(int(r["shell_index"]), int(r["band_index"]), int(r["wedge_index"]))] = {
            "generated": counts[None][0], "stats": _stage_stats(counts[None])}
    return found


def galaxy_window(conn, at, g, edge_pc):
    """
    A Galaxy Map data window of cells bigger than a sector (MAP.151b,
    `galaxy.window`): `g` cells across around the cell `at`, read in one
    statement from the stored block totals (`galaxy_block_stats`).

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        at (str): A cell key `m.shell.band.wedge` of a level above a sector.
        g (int): Cells across, one of `window.WINDOW_WIDTHS`.
        edge_pc (float): The sector edge in parsecs.

    Returns:
        dict: The shape `window.crop_stages` makes: `at`, `g`, `child_m`,
            `cells`, `children` (the cells that hold a generated sector, as
            `galaxy_stage` lists them) and `sectors` (`None`).

    Raises:
        ValueError: On a malformed key or width, or a sector cell (a window
            of sectors lists names and comes from the parents' stages).
    """
    centre = parse_window_key(at)
    if centre.m == 1:
        raise ValueError("a window of sectors is read from the stages of its parents")
    width = check_width(g)
    cells = window_cells(centre, edge_pc, width)
    found = galaxy_block_stats(conn, centre.m, [(c.shell, c.band, c.wedge) for c in cells])
    children = [{"shell": shell, "band": band, "wedge": wedge, **value}
                for (shell, band, wedge), value in found.items()]
    children.sort(key=lambda c: (c["shell"], -c["band"], c["wedge"]))
    return {
        "at": format_drill_key(centre), "g": width, "child_m": centre.m,
        "cells": [[c.shell, c.band, c.wedge] for c in cells], "children": children, "sectors": None,
    }


def galaxy_nebula_cover(conn, at=None, g=None, edge_pc=None):
    """
    The nebula cover of the cells a Galaxy Map stage or window shows (MAP.174): the share of each cell inside a nebula,
    for the cells that hold a generated sector and the cells that do not, from the stored aggregates
    (`db.nebula_blocks`, MAP.172), which are made from the seed's nebula field and the stored nebulae and so need no
    sector. It is read beside the stage and the window, not in them, because those are cached by the sectors that
    changed and the aggregates change without one (a re-drawn field, a nebula stored after a run's sectors).

    Cells of a sector have no stored cover (the aggregates are of the blocks above a sector), so their list is empty.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        at (str or None): A block key `m.shell.band.wedge`. Without `g`, the stage of `at` (none: the galaxy); with
            `g`, the middle cell of a window.
        g (int or str, optional): The window's width in cells, as `galaxy_window` takes it.
        edge_pc (float, optional): The sector edge in parsecs (a window needs it); read when left out.

    Returns:
        dict: `at` (the canonical key or `None`), `child_m` (the cells' size) and `cells` (`[[shell, band, wedge,
            cover], ...]`, sorted, only cells with a cover above zero).

    Raises:
        ValueError: On a malformed or impossible `at` or `g`.
    """
    if g is not None:
        centre = parse_window_key(at or "")
        width = check_width(g)
        child_m = centre.m
        if child_m == 1:
            return {"at": format_drill_key(centre), "child_m": 1, "cells": []}
        edge_pc = galaxy_edge_pc(conn) if edge_pc is None else edge_pc
        cells = [(c.shell, c.band, c.wedge) for c in window_cells(centre, edge_pc, width)]
        key = format_drill_key(centre)
    elif at is None or at == "":
        found = nebula_blocks.covers(conn, STAGE_TOP)
        return {"at": None, "child_m": STAGE_TOP,
                "cells": [[*cell, cover] for cell, cover in sorted(found.items())]}
    else:
        block = parse_stage_key(at)
        child_m = block.m // DRILL_STEP
        key = format_drill_key(block)
        if child_m == 1:
            return {"at": key, "child_m": 1, "cells": []}
        cells = [(c.shell, c.band, c.wedge) for _shell, group in drill_children(block) for c in group]
    found = nebula_blocks.covers(conn, child_m, cells)
    return {"at": key, "child_m": child_m, "cells": [[*cell, cover] for cell, cover in sorted(found.items())]}


def _sector_sums(row):
    """One sector's sums, in `_add_stage_sums`' order, for its own `stats`."""
    stars = row["actual_stars"] or 0
    return [1, row["actual_systems"] or 0, row["expected_systems"] or 0.0, stars,
            (row["mean_age_gy"] or 0.0) * stars, row["total_luminosity_sol"] or 0.0, row["habitable_worlds"] or 0,
            row["total_mass_kg"] or 0.0]


def _stage_sectors_with_stats(where, extra=""):
    """The placed sectors matching `where` (on `sectors`' own columns),
    each with the `sector_stats` the map colors from (DB.14):
    `actual_systems`, `expected_systems`, `actual_stars`, `mean_age_gy`,
    `total_luminosity_sol`, `habitable_worlds` and `total_mass_kg` (NULL without a row or a fill)."""
    return (
        f"SELECT {extra}s.shell_index, s.band_index, s.wedge_index, st.actual_systems, st.expected_systems,"
        " st.actual_stars, st.mean_age_gy, st.total_luminosity_sol, st.habitable_worlds, st.total_mass_kg"
        f" FROM (SELECT * FROM sectors WHERE {where}) s LEFT JOIN sector_stats st"
        " ON st.shell_index = s.shell_index AND st.band_index = s.band_index"
        " AND st.wedge_index = s.wedge_index"
    )


_STAGE_LOOK_SUMS = (
    "COUNT(*) AS n, SUM(COALESCE(actual_systems, 0)) AS systems, SUM(COALESCE(expected_systems, 0)) AS expected,"
    " SUM(COALESCE(actual_stars, 0)) AS stars, SUM(COALESCE(mean_age_gy, 0) * COALESCE(actual_stars, 0)) AS age_sum,"
    " SUM(COALESCE(total_luminosity_sol, 0)) AS luminosity, SUM(COALESCE(habitable_worlds, 0)) AS habitable,"
    " SUM(CAST(COALESCE(total_mass_kg, 0) AS DECIMAL(65, 0))) AS mass"
)
"""The per-group sums `galaxy_stage` folds into each child's stats."""


def _add_stage_sums(counts, child, row):
    """Adds one group's sums (`_STAGE_LOOK_SUMS`) to `counts[child]`:
    `[sectors, systems, expected systems, stars, stars times mean age,
    luminosity, habitable worlds, mass in kilograms]`."""
    sums = counts.setdefault(child, [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    sums[0] += int(row["n"])
    sums[1] += float(row["systems"] or 0.0)
    sums[2] += float(row["expected"] or 0.0)
    sums[3] += float(row["stars"] or 0.0)
    sums[4] += float(row["age_sum"] or 0.0)
    sums[5] += float(row["luminosity"] or 0.0)
    sums[6] += float(row["habitable"] or 0.0)
    sums[7] += float(row["mass"] or 0.0)


def _stage_stats(sums):
    """A child's `stats` from its sums (DB.14): the `systems` and `stars`
    its generated sectors hold, the `expected_systems` the density model
    gave them, `mean_age_gy` (the stars' mean age, each sector's mean
    weighted by its star count; `None` without stars),
    `luminosity_sol` (the sum of their luminosity), `habitable_worlds`
    (MAP.143: the Blue or Green planets and moons they hold) and `mass_kg` (MAP.322: the total mass of what they hold, in
    kilograms, 0 for a sector with no figure), which the map ranks a region's sector names by."""
    _n, systems, expected, stars, age_sum, luminosity, habitable, mass = sums
    return {
        "systems": int(systems),
        "expected_systems": round(expected, 3),
        "stars": int(stars),
        "mean_age_gy": round(age_sum / stars, 4) if stars else None,
        "luminosity_sol": luminosity,
        "habitable_worlds": int(habitable),
        "mass_kg": float(mass or 0.0),
    }


def _stage_children(counts):
    """`galaxy_stage`'s `children` list from `{DrillBlock: sums}`
    (`_add_stage_sums`)."""
    return [
        {"shell": b.shell, "band": b.band, "wedge": b.wedge, "generated": sums[0], "stats": _stage_stats(sums)}
        for b, sums in sorted(counts.items(), key=lambda item: (item[0].shell, -item[0].band, item[0].wedge))
    ]


GALAXY_LOCATE_LIMIT = 8
"""int: Most matches `galaxy_locate` returns."""


def _locate_match(conn, column, term):
    """`galaxy_locate`'s name test: whole words, the last one as typed so
    far (`_name_match`). A lone word too short for the FULLTEXT index
    matches the start of the name through its ordinary index instead of
    scanning every row with REGEXP."""
    words = re.findall(r"\w+", term)
    if len(words) == 1 and len(words[0]) < _fulltext_min_word(conn):
        return f"{column} LIKE ? ESCAPE '\\\\'", [_search_like_pattern(term)[1:]]
    return _name_match(conn, column, term, prefix_last=True)


def galaxy_locate(conn, term, limit=GALAXY_LOCATE_LIMIT):
    """
    Sectors, star systems, planets and moons whose name contains `term`
    (asteroid belts have no names), with each one's sector address,
    reference and parent chain, for the Galaxy Map's address bar (the
    drill-down's section 9.3): picking a match flies to that sector, and a
    planet or moon links through `/object/<ref>` (NAV.9). Exact names come
    first, then names that start with `term`, then the rest, by name.
    Sectors without an address (placed before the grid) and
    systems outside a sector are left out.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        term (str): Part of a name; blank finds nothing.
        limit (int): Most matches to return.

    Returns:
        list[dict]: `{kind ("sector", "system", "planet" or "moon"), ref,
            id, name, parents ([{ref, kind, name}], sector down to the
            direct parent), sector_id, sector_name, shell_index, band_index, wedge_index}`; a
            planet or moon also has `system_id` and `system_name`.
    """
    term = (term or "").strip()
    if not term:
        return []
    order = "CASE WHEN {col} = ? THEN 0 WHEN {col} LIKE ? ESCAPE '\\\\' THEN 1 ELSE 2 END, {col}, id"
    prefix = _search_like_pattern(term)[1:]
    sector_match, sector_params = _locate_match(conn, "name", term)
    sectors = conn.execute(
        "SELECT id, name, shell_index, band_index, wedge_index FROM sectors "
        f"WHERE shell_index IS NOT NULL AND {sector_match} "
        f"ORDER BY {order.format(col='name')} LIMIT ?",
        (*sector_params, term, prefix, limit),
    ).fetchall()
    system_match, system_params = _locate_match(conn, "ss.name", term)
    systems = conn.execute(
        "SELECT ss.id, ss.name, s.id AS sector_id, s.name AS sector_name, "
        "s.shell_index, s.band_index, s.wedge_index "
        "FROM star_systems ss JOIN sectors s ON s.id = ss.sector_id "
        f"WHERE s.shell_index IS NOT NULL AND {system_match} "
        f"ORDER BY {order.format(col='ss.name').replace(', id', ', ss.id')} LIMIT ?",
        (*system_params, term, prefix, limit),
    ).fetchall()
    body_order = order.format(col="b.name").replace(", id", ", b.id")
    planets = conn.execute(
        "SELECT b.id, b.name, b.star_system_id, ss.name AS system_name, s.id AS sector_id, s.name AS sector_name, "
        "s.shell_index, s.band_index, s.wedge_index "
        "FROM planets b JOIN star_systems ss ON ss.id = b.star_system_id JOIN sectors s ON s.id = ss.sector_id "
        f"WHERE s.shell_index IS NOT NULL AND {_locate_match(conn, 'b.name', term)[0]} "
        f"ORDER BY {body_order} LIMIT ?",
        (*_locate_match(conn, "b.name", term)[1], term, prefix, limit),
    ).fetchall()
    moons = conn.execute(
        "SELECT b.id, b.name, b.star_system_id, b.planet_id, p.name AS planet_name, ss.name AS system_name, "
        "s.id AS sector_id, s.name AS sector_name, s.shell_index, s.band_index, s.wedge_index "
        "FROM moons b JOIN planets p ON p.id = b.planet_id JOIN star_systems ss ON ss.id = b.star_system_id "
        "JOIN sectors s ON s.id = ss.sector_id "
        f"WHERE s.shell_index IS NOT NULL AND {_locate_match(conn, 'b.name', term)[0]} "
        f"ORDER BY {body_order} LIMIT ?",
        (*_locate_match(conn, "b.name", term)[1], term, prefix, limit),
    ).fetchall()

    def where(r):
        return {"sector_id": r["sector_id"], "sector_name": r["sector_name"],
                "shell_index": int(r["shell_index"]), "band_index": int(r["band_index"]), "wedge_index": int(r["wedge_index"])}

    def sector_parent(r):
        return {"ref": object_ref.format("sector", r["sector_id"]), "kind": "sector", "name": r["sector_name"]}

    def system_parent(r, system_id):
        return {"ref": object_ref.format("system", system_id), "kind": "system", "name": r["system_name"]}

    found = [{
        "kind": "sector", "ref": object_ref.format("sector", r["id"]), "id": r["id"], "name": r["name"],
        "sector_id": r["id"], "sector_name": r["name"], "parents": [],
        "shell_index": int(r["shell_index"]), "band_index": int(r["band_index"]), "wedge_index": int(r["wedge_index"]),
    } for r in sectors] + [{
        "kind": "system", "ref": object_ref.format("system", r["id"]), "id": r["id"], "name": r["name"],
        "sector_id": r["sector_id"], "sector_name": r["sector_name"], "parents": [sector_parent(r)],
        "shell_index": int(r["shell_index"]), "band_index": int(r["band_index"]), "wedge_index": int(r["wedge_index"]),
    } for r in systems] + [{
        "kind": "planet", "ref": object_ref.format("planet", r["id"]), "id": r["id"], "name": r["name"],
        "system_id": r["star_system_id"], "system_name": r["system_name"],
        "parents": [sector_parent(r), system_parent(r, r["star_system_id"])], **where(r),
    } for r in planets] + [{
        "kind": "moon", "ref": object_ref.format("moon", r["id"]), "id": r["id"], "name": r["name"],
        "system_id": r["star_system_id"], "system_name": r["system_name"],
        "parents": [sector_parent(r), system_parent(r, r["star_system_id"]),
                    {"ref": object_ref.format("planet", r["planet_id"]), "kind": "planet", "name": r["planet_name"]}],
        **where(r),
    } for r in moons]
    lowered = term.lower()

    def rank(match):
        name = (match["name"] or "").lower()
        return (0 if name == lowered else 1 if name.startswith(lowered) else 2, name, match["kind"], match["id"])

    return sorted(found, key=rank)[:limit]


def galaxy_stage_keys(shell, band, wedge):
    """The stage keys a change to sector `(shell, band, wedge)` makes stale:
    `"galaxy"` and its blocks above level 1 up to the stages' top (the level-243 blocks, `drill.STAGE_TOP`)."""
    return ["galaxy"] + [format_drill_key(b) for b in drill_chain_of(shell, band, wedge)[:-1] if b.m <= STAGE_TOP]


def galaxy_node_keys(shell, band, wedge):
    """The node keys a change to sector `(shell, band, wedge)` makes stale (MAP.260): every block of the ladder over
    it, from the root down, and the sector's own node (`drill_chain_of`)."""
    return [format_drill_key(b) for b in drill_chain_of(shell, band, wedge)]


TILE_FORMAT = 9
"""int: The version of what the Galaxy Map's stored nodes and drill-down stages hold and how they are packed
(`web/lib/tilecache.trim_node`, `tilepack`, `db.nodes.galaxy_nodes`, `galaxy_stage`). A change that alters those bytes for
the same rows goes up by one (`tests/test_tile_format.py` fails until it does), and the stored files are rebuilt once;
a release that does not alter them leaves every stored node valid (PERF.133). 2: the block ladder's levels above 243
(MAP.258) moved the rows and wedges of every level, so a stored stage names other blocks. 3: nodes replace the cube tiles
(MAP.262), whose files are dropped with the generation. 5: a node of a charted region carries its `name` and `parent`
(API.35), and so does the stage of a 9-block and each charted 9-block child of a stage. 6: a 9-block's node carries its
`extent` (MAP.319). 7: a 9-block's node carries `charted` (its sectors and share), and a node of 27, 81 or 243 sectors lists the charted regions inside it, `regions` (MAP.314). 8: the node's `worlds` section is `interests` (MAP.305): the systems worth a mark, each with its `reasons` (a habitable world, a gas giant). 9: a node of 729, 243, 81, 27 or 9 sectors lists its nebula regions, `nebula_regions` (MAP.173)."""

GALAXY_CHANGES_MAX_SECTORS = 1000
"""int: Most changed sectors `galaxy_changes` lists nodes for. More than
that (a big generation run, say) reports `full` instead -- refetching
everything is cheaper than invalidating tens of thousands of nodes one
by one."""

_STATE_TOKEN_RE = re.compile(r"^([0-9a-f]{16})\.([0-9a-f]{16})\.(\d+)\.(\d+)\.(\d{17}|0)\.(\d+)\.(\d+)$")


def _state_token(state):
    """`galaxy_content_state`'s dict as the opaque string the API hands
    out and `galaxy_changes` reads back."""
    return ("{base}.{tile_base}.{deletions}.{sector_max_id}.{sector_modified}.{system_max_id}.{bright_max_id}"
            .format(**state))


def _parse_state_token(token):
    """The dict `_state_token` encoded, or `None` for anything else."""
    match = _STATE_TOKEN_RE.match(str(token or ""))
    if not match:
        return None
    base, tile_base, deletions, sector_max_id, sector_modified, system_max_id, bright_max_id = match.groups()
    return {
        "base": base, "tile_base": tile_base, "deletions": int(deletions), "sector_max_id": int(sector_max_id),
        "sector_modified": sector_modified, "system_max_id": int(system_max_id),
        "bright_max_id": int(bright_max_id),
    }


def _timestamp_digits(value):
    """A `TIMESTAMP(3)` value as 17 digits (`YYYYMMDDHHMMSSmmm`), or `"0"`
    for `NULL` -- compact, and orders the same way the timestamp does."""
    if value is None:
        return "0"
    if isinstance(value, str):
        value = datetime.datetime.fromisoformat(value)
    return value.strftime("%Y%m%d%H%M%S") + f"{value.microsecond // 1000:03d}"


def _timestamp_literal(digits):
    """`_timestamp_digits`' output back as a literal MySQL compares
    against a `TIMESTAMP(3)` column."""
    if digits == "0":
        return "1970-01-01 00:00:01.000"
    d = digits
    return f"{d[0:4]}-{d[4:6]}-{d[6:8]} {d[8:10]}:{d[10:12]}:{d[12:14]}.{d[14:17]}"


def galaxy_shape_inputs(conn):
    """
    What the stored galaxy shape is made from, as plain values: the `galaxy_shape` row, the levers it was told
    (`galaxy_shape_pins`) and a signature of the `galaxy_components` rows (their count, first and last id and sums). The
    shape itself (`galaxy_density_shape`) is a function of these and of the code, and building it works out a mean system
    mass for every stored cluster, which took minutes on a galaxy of ten thousands of clusters: `galaxy_content_state` runs
    on a page view's freshness check, so it reads these instead. `None` when the galaxy was never planned.
    """
    row = conn.execute("SELECT * FROM galaxy_shape WHERE id = 1").fetchone()
    if row is None:
        return None
    components = None
    if conn.execute("SELECT 1 FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE() "
                    "AND TABLE_NAME = 'galaxy_components'").fetchone():
        components = dict(conn.execute(
            "SELECT COUNT(*) AS n, MIN(id) AS first_id, MAX(id) AS last_id, SUM(mass_solar) AS mass, SUM(age_gy) AS age, "
            "SUM(metallicity_fe_h) AS fe_h FROM galaxy_components").fetchone())
    return {"row": {key: (value.hex() if isinstance(value, (bytes, bytearray)) else value) for key, value in dict(row).items()},
            "pins": {pin["name"]: pin["value"] for pin in conn.execute("SELECT name, value FROM galaxy_shape_pins").fetchall()},
            "components": components}


def galaxy_content_state(conn):
    """
    Everything `db.nodes.galaxy_nodes`' output depends on, summarized as a handful
    of numbers:

    - `base`: a hash of the stored galaxy shape, the bright-star scatter
      (`bright_stars`' highest id and the scatter's seed) and this code's
      version. Planned slots, density clouds, `edge_pc` and `has_shape`
      depend on the shape, every node's `stars` on the scatter, so a new
      `base` means every page and count cached by the stamp is stale.
    - `tile_base`: the same without the release: the shape, the scatter's
      seed, the naming key, `TILE_FORMAT` and the generator epoch. The
      stored nodes follow this one, so a release that changes neither the
      node format nor what the generator makes keeps every stored node
      (PERF.133).
    - `deletions`/`sector_max_id`: the deletion epoch (`store.
      note_sector_deleted`, PERF.38: it goes up with every sector deleted,
      where a count of the placed sectors cost 0.13 s per million of them
      on every check), and the highest sector id. Together they tell a new
      sector (higher id) from a deleted one (a new epoch).
    - `sector_modified`: the newest `sectors.modified_at` (v27). A
      sector's own edit (a rename, say) bumps it, and so does deleting one
      of its systems (`store.touch_sector`).
    - `bright_max_id`: the highest `bright_stars` id; it grows while a
      backfill runs, so it is kept out of `base` and `galaxy_changes`
      reports it as `busy` rather than `full` (PERF.34).
    - `system_max_id`: the highest star-system id. A new system changes
      its sector's system count. System edits don't touch a node (nodes
      only show the count), so `star_systems.modified_at` isn't used.

    Cheap by design -- one primary-key read and four index-only maxima --
    since the web layer checks it about once a minute per database.

    Returns:
        dict: The keys above; `sector_modified` as `_timestamp_digits`.
    """
    maxima = conn.execute(
        "SELECT (SELECT COALESCE(MAX(id), 0) FROM sectors) AS sector_max_id, "
        "(SELECT MAX(modified_at) FROM sectors) AS sector_modified, "
        "(SELECT COALESCE(MAX(id), 0) FROM star_systems) AS system_max_id"
    ).fetchone()
    bright = conn.execute(
        "SELECT (SELECT COALESCE(MAX(id), 0) FROM bright_stars) AS max_id, "
        "(SELECT bright_star_seed FROM galaxy_shape WHERE id = 1) AS seed"
    ).fetchone()
    shape_inputs = galaxy_shape_inputs(conn)
    base = hashlib.sha256(json.dumps(
        {"shape": shape_inputs, "version": __version__,
         "bright_stars": bright["seed"],
         # GEN.70: the naming key renames the codec-named objects in every node and page.
         "naming_key": naming_key.active_key()}, sort_keys=True, default=str,
    ).encode("utf-8")).hexdigest()[:16]
    tile_base = hashlib.sha256(json.dumps(
        {"shape": shape_inputs, "tile_format": TILE_FORMAT, "epoch": generator_epoch.generator_epoch(),
         "bright_stars": bright["seed"], "naming_key": naming_key.active_key()}, sort_keys=True, default=str,
    ).encode("utf-8")).hexdigest()[:16]
    return {
        "base": base,
        "tile_base": tile_base,
        "deletions": sector_deletions(conn),
        "sector_max_id": int(maxima["sector_max_id"]),
        "sector_modified": _timestamp_digits(maxima["sector_modified"]),
        "system_max_id": int(maxima["system_max_id"]),
        "bright_max_id": int(bright["max_id"]),
    }


def galaxy_content_stamp(conn, state=None):
    """
    A short token that changes whenever anything `db.nodes.galaxy_nodes` returns
    could change -- a hash of `galaxy_content_state`. The web layer and
    the browser key their node caches by it, so a stale node is never
    served after sectors are generated, renamed or deleted.

    Args:
        state (dict, optional): An already-read `galaxy_content_state`.

    Returns:
        str: 16 hex characters.
    """
    if state is None:
        state = galaxy_content_state(conn)
    return hashlib.sha256(_state_token(state).encode("utf-8")).hexdigest()[:16]


def galaxy_tile_stamp(state):
    """
    `galaxy_content_stamp`'s counterpart for the stored nodes: it changes when a node could change, and not with the
    release alone (see `galaxy_content_state`'s `tile_base`).

    Args:
        state (dict): A `galaxy_content_state`.

    Returns:
        str: 16 hex characters.
    """
    token = "{tile_base}.{deletions}.{sector_max_id}.{sector_modified}.{system_max_id}.{bright_max_id}".format(**state)
    return hashlib.sha256(token.encode("utf-8")).hexdigest()[:16]


def galaxy_changes(conn, since=None):
    """
    What changed in the galaxy's nodes and stages since an earlier state -- how the
    web layer's node cache refreshes only the blocks an edit touched
    instead of throwing every cached node away.

    Uses v27's `sectors.modified_at` (see `galaxy_content_state`): each
    sector edited since `since`, each new sector, and each sector that
    gained a system is "changed", and so is every block over it
    (`galaxy_node_keys`, `galaxy_stage_keys`). A sector
    never moves once placed, so its address is where it was before too.

    Deletions leave no row behind, so a deleted sector can't be located;
    the deletion epoch moves, and the answer is `full` instead. So is a new
    shape, node format or generator epoch (`tile_base`), an unreadable `since`, or more than
    `GALAXY_CHANGES_MAX_SECTORS` changed sectors. A new release alone (`base`) makes `full` true for what is cached by
    `stamp` (pages, counts) but not `tile_full`, which is what the stored nodes follow (PERF.133).

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        since (str or None): A `state` from an earlier call (or
            `GET /api/galaxy/stamp`), or `None`.

    Returns:
        dict: `stamp` (`galaxy_content_stamp` at the call), `tile_stamp`
            (`galaxy_tile_stamp`: what the stored nodes are keyed by),
            `state` (the token to pass as `since` next time), `full`
            (everything cached by `stamp` may have changed), `tile_full`
            (every node may have changed), `stages` (sorted
            keys of the drill-down stages whose counts changed,
            `galaxy_stage_keys`; empty when `tile_full`), `nodes` (sorted keys of the nodes that changed,
            `galaxy_node_keys`: every block over a changed sector; empty when `tile_full`), and, only when true, `busy`: `full` because so much
            changed at once (a generation run, a backfill), not because
            something was deleted, so a cache may keep its nodes a while
            longer instead of refetching all of them under that load
            (PERF.34). Unless `full`, also `sector_ids` (the row ids of the
            sectors edited, added or given a system since `since`, for the
            page cache's tags, PERF.104) and `grew` (a sector or a system
            was added).
    """
    state = galaxy_content_state(conn)
    result = {
        "stamp": galaxy_content_stamp(conn, state), "tile_stamp": galaxy_tile_stamp(state),
        "state": _state_token(state), "full": False, "tile_full": False, "stages": [], "nodes": [],
    }

    def everything(busy=False):
        result["full"] = result["tile_full"] = True
        if busy:
            result["busy"] = True
        return result

    previous = _parse_state_token(since)
    if previous is None or previous["tile_base"] != state["tile_base"]:
        return everything()
    if previous["base"] != state["base"]:
        # A release that changes neither the node format nor the generator epoch: whatever else is cached by
        # `stamp` (pages, counts) is stale, the stored nodes are not (PERF.133). The nodes that changed are
        # still listed below.
        result["full"] = True
    if previous == state:
        return result
    if previous["bright_max_id"] != state["bright_max_id"]:
        # A backfill is adding stars to cells all over (PERF.34): which
        # nodes they touch isn't worked out, so everything counts as changed.
        return everything(busy=True)

    if previous["deletions"] != state["deletions"]:
        return everything()

    rows = conn.execute(
        """
        SELECT center_x_pc, center_y_pc, center_z_pc, shell_index, band_index, wedge_index
        FROM sectors
        WHERE center_x_pc IS NOT NULL
          AND (id > ? OR modified_at > ?
               OR id IN (SELECT sector_id FROM star_systems WHERE id > ?))
        LIMIT ?
        """,
        (
            previous["sector_max_id"], _timestamp_literal(previous["sector_modified"]),
            previous["system_max_id"], GALAXY_CHANGES_MAX_SECTORS + 1,
        ),
    ).fetchall()
    if len(rows) > GALAXY_CHANGES_MAX_SECTORS:
        return everything(busy=True)

    stages = set()
    node_keys = set()
    # The page cache tags every sector, placed or not (PERF.104).
    result["sector_ids"] = sorted(row["id"] for row in conn.execute(
        """
        SELECT id FROM sectors
        WHERE id > ? OR modified_at > ? OR id IN (SELECT sector_id FROM star_systems WHERE id > ?)
        LIMIT ?
        """,
        (previous["sector_max_id"], _timestamp_literal(previous["sector_modified"]),
         previous["system_max_id"], GALAXY_CHANGES_MAX_SECTORS + 1)).fetchall())
    if len(result["sector_ids"]) > GALAXY_CHANGES_MAX_SECTORS:
        result["sector_ids"] = []
        return everything(busy=True)
    result["grew"] = (previous["sector_max_id"] != state["sector_max_id"]
                      or previous["system_max_id"] != state["system_max_id"])
    for r in rows:
        if sector_address(r) is not None:
            stages.update(galaxy_stage_keys(*sector_address(r)))
            node_keys.update(galaxy_node_keys(*sector_address(r)))
    result["stages"] = sorted(stages)
    result["nodes"] = sorted(node_keys)
    return result

# ---------------------------------------------------------------------
# Faceted search -- backs GET /api/search and html/search.py. Ported
# from html/search.py's own query layer (same SQL, same "which result
# panels actually have a reason to run" logic) so the CGI page and the
# API build on the exact same functions rather than the CGI page querying
# the database directly -- see this module's own docstring.
# ---------------------------------------------------------------------

SEARCH_RESULT_LIMIT = 300
SEARCH_AUTOCOMPLETE_LIMIT = 500

SEARCH_TAG_FACETS = (
    "type", "spectral", "luminosity",
    "class", "body", "life", "equipment",
    "moon_class", "moon_body", "moon_life", "moon_equipment",
    *score_filters.FACETS, *(f"moon_{name}" for name in score_filters.FACETS),
    "density",
    "phenomenon", "phenomenon_class",
)

_SEARCH_PHENOMENON_LABELS = {
    "nebula": "Nebula", "asteroid_field": "Asteroid Field", "black_hole": "Black Hole",
    "neutron_star": "Neutron Star", "supernova_remnant": "Supernova Remnant",
    "rogue_planet": "Rogue Planet", "interstellar_comet": "Interstellar Comet", "quasar": "Quasar",
}

_SEARCH_PHENOMENON_CLASS_COLUMNS = {
    "nebula": "nebula_class", "supernova_remnant": "remnant_class", "asteroid_field": "field_class",
}
"""The phenomenon types with a class (schema v38) and its column. A
`phenomenon_class` tag is `"<type>:<class>"`, e.g. `"nebula:D"`, since an
asteroid field's letters overlap a nebula's."""

# generation/star.py's own Yerkes-class-to-descriptive-label mapping (Star.__init__),
# duplicated here (not imported) since it's a plain literal there, not an
# importable constant. "D" is an alternate white-dwarf code seen elsewhere
# in generation/star.py alongside "VII".
_SEARCH_YERKES_LABELS = {
    "0": "Hypergiant",
    "IA": "Supergiant",
    "IAB": "Intermediate-size Luminous Supergiant",
    "IB": "Less Luminous Supergiant",
    "II": "Bright Giant",
    "III": "Giant",
    "IV": "Subgiant",
    "V": "Main Sequence",
    "VII": "White Dwarf",
    "D": "White Dwarf",
}
_SEARCH_YERKES_ORDER = ["0", "IA", "IAB", "IB", "II", "III", "IV", "V", "VII", "D"]

_SEARCH_BODY_LABELS = {"t": "Terrestrial", "g": "Gas Giant"}


FULLTEXT_STOPWORDS = frozenset((
    "a", "about", "an", "are", "as", "at", "be", "by", "com", "de", "en", "for", "from", "how", "i", "in",
    "is", "it", "la", "of", "on", "or", "that", "the", "this", "to", "was", "what", "when", "where", "who",
    "will", "with", "und", "www",
))
"""frozenset: InnoDB's default full-text stopwords, which its index
leaves out -- `_name_match` checks these words with REGEXP instead."""

_token_sizes = {}


def _fulltext_token_sizes(conn):
    """The server's `innodb_ft_min_token_size` and `_max_token_size` (3
    and 84 unless changed), cached per database: shorter and longer words
    aren't in a FULLTEXT index."""
    config = getattr(conn, "_config", None)
    key = config._key() if config is not None else None
    if key not in _token_sizes:
        try:
            row = conn.execute(
                "SELECT @@innodb_ft_min_token_size AS n, @@innodb_ft_max_token_size AS x").fetchone()
            _token_sizes[key] = (int(row["n"]), int(row["x"]))
        except pymysql.MySQLError:
            _token_sizes[key] = (3, 84)
    return _token_sizes[key]


def _fulltext_min_word(conn):
    """The server's `innodb_ft_min_token_size`: shorter words aren't in a
    FULLTEXT index."""
    return _fulltext_token_sizes(conn)[0]


def _word_pattern(word, prefix=False):
    """A REGEXP matching `word` as a whole word (or, with `prefix`, the
    start of one) -- `word` is `\\w` characters only, so needs no escaping."""
    end = "" if prefix else "([^[:alnum:]_]|$)"
    return f"(^|[^[:alnum:]_]){word}{end}"


def _name_match(conn, column, term, prefix_last=False):
    """
    A WHERE fragment matching names that contain every word of `term` as
    a whole word (PERF.16, Boss 2026-10-01: "full text index, to match
    whole words"), so "ara" no longer finds "Kemaral". Words long enough
    for the FULLTEXT index (v46) go through `MATCH ... AGAINST` in boolean
    mode; shorter words, stopwords (Greek letters like "Mu", numerals
    like "IV") and words too long for the index are checked with a
    word-boundary REGEXP on the rows the index found -- or on every row,
    when no word fits the index.

    Args:
        conn (planetgen.db.store.Connection): An open connection.
        column (str): The (aliased) `name` column, a literal.
        term (str): What was typed.
        prefix_last (bool): Let the last word match the start of a word,
            for the Galaxy Map's address bar, which searches as you type.

    Returns:
        tuple: `(sql, params)`; `("1 = 0", [])` when `term` has no words.
    """
    words = re.findall(r"\w+", term or "")
    if not words:
        return "1 = 0", []
    min_word, max_word = _fulltext_token_sizes(conn)
    against, clauses, params = [], [], []
    for index, word in enumerate(words):
        prefix = prefix_last and index == len(words) - 1
        if min_word <= len(word) <= max_word and word.lower() not in FULLTEXT_STOPWORDS:
            against.append(f"+{word}*" if prefix else f"+{word}")
        else:
            clauses.append(f"{column} REGEXP ?")
            params.append(_word_pattern(word, prefix))
    if against:
        clauses.insert(0, f"MATCH({column}) AGAINST (? IN BOOLEAN MODE)")
        params.insert(0, " ".join(against))
    return " AND ".join(clauses), params


def _search_match(conn, column, term):
    """
    The name test of a search result panel (UX.97): `_name_match`, except that a lone word too short for the
    FULLTEXT index matches the start of the name through the ordinary name index, instead of a REGEXP over every
    row (a 2-letter word took 8.4 s on 2,000,000 systems). "ar" therefore finds "Arcturus", not "Alpha Ar".
    """
    words = re.findall(r"\w+", term or "")
    if len(words) == 1 and len(words[0]) < _fulltext_min_word(conn):
        return f"{column} LIKE ? ESCAPE '\\\\'", [_search_like_pattern(words[0])[1:]]
    return _name_match(conn, column, term)


def _search_like_pattern(term):
    """Escapes `%`/`_`/`\\` in a user-supplied substring so it's safe to
    use as a SQL LIKE pattern (paired with `ESCAPE '\\\\'` in the query)."""
    escaped = term.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")
    return f"%{escaped}%"


def _append_size_clause(clauses, params, column, size_range):
    """
    Appends a `column BETWEEN`/`>=`/`<=` clause for a min/max size filter
    (in km, over any of `stars`/`planets`/`moons`' own `radius_km`) to an
    in-progress `clauses`/`params` pair, shared by
    `_search_result_stars`/`_search_result_planets`/`_search_result_moons`
    rather than duplicating the same three-way "both bounds, min only,
    max only" branching in each.

    Args:
        clauses (list[str]): SQL WHERE fragments to append to, in place.
        params (list): Query parameters to append to, in place, matching
                       `clauses`' own `?` placeholders in order.
        column (str): The (already-aliased, e.g. `"p.radius_km"`) column
                      to filter on.
        size_range (tuple[float or None, float or None] or None): `(min_km,
            max_km)` -- either may be `None` for "no lower/upper bound".
            `None` itself (not a tuple) means no size filter at all.
    """
    if size_range is None:
        return
    min_km, max_km = size_range
    if min_km is not None and max_km is not None:
        clauses.append(f"{column} BETWEEN ? AND ?")
        params.extend((min_km, max_km))
    elif min_km is not None:
        clauses.append(f"{column} >= ?")
        params.append(min_km)
    elif max_km is not None:
        clauses.append(f"{column} <= ?")
        params.append(max_km)


# --- Facet option discovery -- each returns a list of {"value", "label",
# "count", "tooltip"} dicts, one per distinct value actually present in
# the database (never a fixed/static enumeration). ---

def _search_facet_type(conn):
    star_c = conn.execute("SELECT COUNT(*) AS c FROM stars").fetchone()["c"]
    planet_c = conn.execute("SELECT COUNT(*) AS c FROM planets").fetchone()["c"]
    moon_c = conn.execute("SELECT COUNT(*) AS c FROM moons").fetchone()["c"]
    belt_c = conn.execute("SELECT COUNT(*) AS c FROM asteroid_belts").fetchone()["c"]
    opts = []
    for value, label, count in (
        ("star", "Stars", star_c),
        ("planet", "Planets", planet_c),
        ("moon", "Moons", moon_c),
        ("belt", "Asteroid Belts", belt_c),
    ):
        if count:
            opts.append({"value": value, "label": label, "count": count, "tooltip": None})
    return opts


def _search_facet_spectral(conn):
    """Every spectral class the generator produces (`SPECTRAL_CLASS_COLORS`), with its star count (0 when none yet),
    then any other letter the stars hold."""
    counts = {row["v"]: row["c"] for row in conn.execute(
        "SELECT SUBSTR(star_type, 1, 1) AS v, COUNT(*) AS c FROM stars GROUP BY v ORDER BY v").fetchall()}
    letters = [*SPECTRAL_CLASS_COLORS, *sorted(letter for letter in counts if letter not in SPECTRAL_CLASS_COLORS)]
    opts = []
    for letter in letters:
        color = SPECTRAL_CLASS_COLORS.get(letter)
        tip = f"{color} star" if color else None
        opts.append({"value": letter, "label": f"{letter}-Type Star", "count": counts.get(letter, 0), "tooltip": tip})
    return opts


def _search_facet_luminosity(conn):
    """Every Yerkes class (`_SEARCH_YERKES_ORDER`) with its star count (0 when none yet), then any other class the
    stars hold."""
    counts = {row["v"]: row["c"] for row in conn.execute(
        "SELECT yerkes_class AS v, COUNT(*) AS c FROM stars WHERE yerkes_class IS NOT NULL GROUP BY v").fetchall()}
    classes = [*_SEARCH_YERKES_ORDER, *sorted(value for value in counts if value not in _SEARCH_YERKES_ORDER)]
    return [
        {"value": value, "label": _SEARCH_YERKES_LABELS.get(value, value), "count": counts.get(value, 0),
         "tooltip": f"Yerkes class {value}"}
        for value in classes
    ]


def _search_facet_class(conn):
    rows = conn.execute(
        "SELECT planet_class AS v, COUNT(*) AS c FROM planets WHERE planet_class IS NOT NULL GROUP BY v ORDER BY v"
    ).fetchall()
    opts = []
    for row in rows:
        info = PLANET_CLASSES.get(row["v"], {})
        opts.append({"value": row["v"], "label": f"Class {row['v']} Planet", "count": row["c"],
                     "tooltip": info.get("description")})
    return opts


def _search_facet_body(conn):
    rows = conn.execute("SELECT body_type AS v, COUNT(*) AS c FROM planets GROUP BY v ORDER BY v").fetchall()
    return [
        {"value": row["v"], "label": _SEARCH_BODY_LABELS.get(row["v"], row["v"]), "count": row["c"], "tooltip": None}
        for row in rows
    ]


def _life_class_label(class_id):
    """The facet label of a life class id: `"Class 1a: oxygenic light-eaters ..."`."""
    return f"Class {class_id}: {life_classes.describe(class_id).name}" if class_id in life_classes.BY_ID else str(class_id)


def _life_class_facet(conn, table):
    """The life classes living on the bodies of `table` (`planets` or `moons`), with how many bodies each lives on: any class of a body counts, not only the dominant
    one (GEN.1209). A body stored before the `planet_life` / `moon_life` rows existed counts its one class from its own columns."""
    child, key = ("moon_life", "moon_id") if table == "moons" else ("planet_life", "planet_id")
    rows = conn.execute(
        f"SELECT v, COUNT(*) AS c FROM ("
        f"SELECT CONCAT(b.life_class, COALESCE(b.life_class_variant, '')) AS v FROM {table} b WHERE b.life_class IS NOT NULL AND b.life_stage IS NOT NULL "
        f"AND NOT EXISTS (SELECT 1 FROM {child} x WHERE x.{key} = b.id) "
        f"UNION ALL SELECT x.life_class AS v FROM {child} x WHERE x.life_stage IS NOT NULL) t GROUP BY v ORDER BY v"
    ).fetchall()
    return [{"value": row["v"], "label": _life_class_label(row["v"]), "count": row["c"], "tooltip": None} for row in rows]


def _life_class_clause(alias, class_ids, clauses, params):
    """Adds the filter on the life classes `class_ids` (ids such as `"1a"`, `"6"`) of the table `alias`: a body matches when any class that lives on it is one of them
    (GEN.1209), the dominant one in its own columns or another in its `planet_life` / `moon_life` rows."""
    if not class_ids:
        return
    known = [class_id for class_id in class_ids if class_id in life_classes.BY_ID]   # a stray value matches nothing
    if not known:
        clauses.append("1 = 0")
        return
    child, key = ("moon_life", "moon_id") if alias.startswith("m") else ("planet_life", "planet_id")
    pairs = sorted((int(class_id[:1]), class_id[1:]) for class_id in known)
    ids = sorted(known)
    clauses.append(
        f"((({alias}.life_class, COALESCE({alias}.life_class_variant, '')) IN ({','.join('(?, ?)' for _ in pairs)}) AND {alias}.life_stage IS NOT NULL) "
        f"OR EXISTS (SELECT 1 FROM {child} lc WHERE lc.{key} = {alias}.id AND lc.life_stage IS NOT NULL AND lc.life_class IN ({','.join('?' for _ in ids)})))")
    for number, variant in pairs:
        params.extend((number, variant))
    params.extend(ids)


def _search_facet_life(conn):
    return _life_class_facet(conn, "planets")


def _search_equipment_options(conn, table):
    """GEN.655: the items a human needs on a body (`equipment.TICKS`, the required ones), with how many bodies of
    `table` need each; the ones none need are left out."""
    row = conn.execute(f"SELECT {equipment.tick_counts_select('t')} FROM {table} t WHERE t.eq_body IS NOT NULL").fetchone()
    return [{"value": value, "label": label, "count": int(row[f"t{i}"]), "tooltip": f"A human needs {label.lower()} here"}
            for i, (value, label, _sql) in enumerate(equipment.TICKS) if row[f"t{i}"]]


def _search_facet_equipment(conn):
    return _search_equipment_options(conn, "planets")


def _search_facet_moon_equipment(conn):
    return _search_equipment_options(conn, "moons")


def _score_options(conn, table, facet):
    """UX.214: the ticks of one score filter facet (`score_filters.FACETS`) with how many bodies of `table` each holds
    for; the ones none hold for are left out."""
    ticks = score_filters.FACETS[facet][1]
    kind = table[:-1]
    if facet == "supports":
        found = {(row["life_class"], least): int(row[f"t{least}"]) for row in conn.execute(
            score_filters.supports_counts_sql(kind)).fetchall() for least, _word in score_filters.SUPPORTS_LEVELS}
        return [{"value": value, "label": label, "count": found[tuple_], "tooltip": label}
                for value, label, _sql in ticks
                for tuple_ in [(value.split(":")[0], int(value.split(":")[1]))] if found.get(tuple_)]
    row = conn.execute(
        f"SELECT {score_filters.counts_select(facet, 't', kind)} FROM {table} t WHERE t.phi4 IS NOT NULL").fetchone()
    return [{"value": value, "label": label, "count": int(row[f"t{i}"]), "tooltip": label}
            for i, (value, label, _sql) in enumerate(ticks) if row[f"t{i}"]]


def _score_facet_defs(conn):
    """The facet definitions of every score filter, for the planet and the moon panel."""
    return tuple((f"{prefix}{facet}", _score_options(conn, table, facet))
                 for prefix, table in (("", "planets"), ("moon_", "moons")) for facet in score_filters.FACETS)


def _score_clauses(alias, tags, prefix, clauses):
    """Adds one condition for each score filter facet ticked (`score_filters.condition`; a stray value is dropped, which
    leaves a facet of only strays matching nothing). `prefix` is `""` for planets and `"moon_"` for moons."""
    kind = "moon" if prefix else "planet"
    for facet in score_filters.FACETS:
        selected = tags.get(f"{prefix}{facet}", set())
        if selected:
            clauses.append(score_filters.condition(facet, alias, selected, kind) or "1 = 0")


def _score_order(alias, sort, default):
    """The ORDER BY of a body panel. `sort` is a `score_filters.SORT_COLUMNS` key: the worst bodies first, or with `-` in
    front the best first (a domain colour runs Red to Blue, so the stored tier number the other way). Any other `sort`
    (`system`, the page's own order) keeps `default`. Bodies with no score come last either way."""
    key = (sort or "").lstrip("-")
    column = score_filters.SORT_COLUMNS.get(key)
    if column is None:
        return default
    best_first = (sort or "").startswith("-")
    descending = best_first == (key in score_filters.HIGH_FIRST)
    return f"{alias}.{column} IS NULL, {alias}.{column} {'DESC' if descending else 'ASC'}, {default}"


def _equipment_clause(alias, tags, clauses, params):
    """Adds the OR of the ticked equipment items (`equipment.TICKS`; a stray value is dropped) as one condition."""
    if tags:
        clauses.append(equipment.tick_condition(alias, set(tags)) or "1 = 0")


def _search_facet_moon_class(conn):
    rows = conn.execute(
        "SELECT planet_class AS v, COUNT(*) AS c FROM moons WHERE planet_class IS NOT NULL GROUP BY v ORDER BY v"
    ).fetchall()
    opts = []
    for row in rows:
        info = PLANET_CLASSES.get(row["v"], {})
        opts.append({"value": row["v"], "label": f"Class {row['v']} Moon", "count": row["c"],
                     "tooltip": info.get("description")})
    return opts


def _search_facet_moon_body(conn):
    rows = conn.execute("SELECT body_type AS v, COUNT(*) AS c FROM moons GROUP BY v ORDER BY v").fetchall()
    return [
        {"value": row["v"], "label": _SEARCH_BODY_LABELS.get(row["v"], row["v"]), "count": row["c"], "tooltip": None}
        for row in rows
    ]


def _search_facet_moon_life(conn):
    return _life_class_facet(conn, "moons")


def _search_facet_density(conn):
    rows = conn.execute("SELECT density AS v, COUNT(*) AS c FROM asteroid_belts GROUP BY v ORDER BY v").fetchall()
    return [{"value": row["v"], "label": row["v"].capitalize(), "count": row["c"], "tooltip": None} for row in rows]


def _search_phenomenon_tables():
    """`(table, type)` for every standalone phenomenon table."""
    return [(table, type_label) for table, type_label, _descriptor, _radius in _PHENOMENON_TABLES]


def _search_facet_phenomenon(conn):
    options = []
    for table, type_label in _search_phenomenon_tables():
        count = conn.execute(f"SELECT COUNT(*) AS c FROM {table}").fetchone()["c"]
        if count:
            options.append({"value": type_label, "label": _SEARCH_PHENOMENON_LABELS.get(type_label, type_label),
                            "count": count, "tooltip": None})
    return options


def _search_phenomenon_class_label(type_label, value):
    if type_label == "asteroid_field":
        return f"{value} asteroid field"
    entry = tuning.NEBULA_CLASSES.get(value)
    return f"{value}: {entry['name']}" if entry else value


def _search_facet_phenomenon_class(conn):
    options = []
    tables = dict((type_label, table) for table, type_label in _search_phenomenon_tables())
    for type_label, column in _SEARCH_PHENOMENON_CLASS_COLUMNS.items():
        rows = conn.execute(
            f"SELECT {column} AS v, COUNT(*) AS c FROM {tables[type_label]} GROUP BY v ORDER BY v"
        ).fetchall()
        options.extend(
            {"value": f"{type_label}:{row['v']}", "label": _search_phenomenon_class_label(type_label, row["v"]),
             "count": row["c"], "tooltip": None}
            for row in rows if row["v"]
        )
    return options


def _search_result_phenomena(conn, type_tags, class_tags, limit, offset):
    """Every standalone phenomenon matching the Phenomenon and Phenomenon
    Class tags. Tags in one group combine with OR, the two groups with AND,
    and a class tag only narrows its own type (UX.44): Black Hole + Nebula
    + nebula class D is every black hole and the class D nebulae. With no
    Phenomenon tag, the class tags pick their own types."""
    classes = {}
    for tag in class_tags:
        type_label, _sep, value = tag.partition(":")
        if type_label in _SEARCH_PHENOMENON_CLASS_COLUMNS and value:
            classes.setdefault(type_label, set()).add(value)
    parts, params = [], []
    for table, type_label in _search_phenomenon_tables():
        if type_tags and type_label not in type_tags:
            continue
        if class_tags and not type_tags and type_label not in classes:
            continue
        column = _SEARCH_PHENOMENON_CLASS_COLUMNS.get(type_label)
        class_sql = column if column else "NULL"
        where = ""
        if type_label in classes:
            where = f" WHERE {column} IN ({','.join('?' * len(classes[type_label]))})"
            params.extend(sorted(classes[type_label]))
        parts.append(f"SELECT '{type_label}' AS type, id, name, {class_sql} AS phenomenon_class, sector_id "
                     f"FROM {table}{where}")
    if not parts:
        return {"rows": [], "total": 0, "limit": limit, "offset": 0, "truncated": False}
    rows, page = _search_page(
        conn, "SELECT type, id, name, phenomenon_class, sector_id",
        f"FROM ({' UNION ALL '.join(parts)}) AS p", "name, type, id", params, limit, offset,
    )
    return {"rows": _with_refs(conn, [dict(r) for r in rows], None), **page}


def _search_name_list(conn, table, limit=SEARCH_AUTOCOMPLETE_LIMIT):
    rows = conn.execute(f"SELECT DISTINCT name FROM {table} ORDER BY name LIMIT ?", (limit,)).fetchall()
    return [row["name"] for row in rows]


# --- Result panels ---

SEARCH_RESULT_PANELS = ("sectors", "systems", "stars", "planets", "moons", "belts", "phenomena")


SEARCH_COUNT_CAP = 300
"""int: Most matches a result panel counts (PERF.16): past it the panel
says "300+" (`total_capped`) instead of counting every row of a large
galaxy. Paging further raises the cap to two pages past the current one."""


def _search_page(conn, select_sql, from_sql, order_sql, params, limit, offset):
    """
    Runs one result panel's query a page at a time: a count over
    `from_sql` (its FROM/JOIN/WHERE, sharing `params`) that stops at
    `SEARCH_COUNT_CAP` (or two pages past `offset`, if further), then
    `limit` rows from `offset` in `order_sql` order. An `offset` past the
    last match (a stale page link) is pulled back to the last page's
    first row.

    Returns:
        tuple[list, dict]: `(rows, page)` -- `page` is the panel's
            `{"total", "total_capped", "limit", "offset", "truncated"}`
            (`total_capped`: there are more than `total` matches;
            `truncated`: `rows` holds fewer than `total`, i.e. there are
            more pages).
    """
    cap = max(SEARCH_COUNT_CAP, offset + 2 * limit)
    total = conn.execute(
        f"SELECT COUNT(*) AS n FROM (SELECT 1 AS one {from_sql} LIMIT ?) capped", list(params) + [cap + 1]
    ).fetchone()["n"]
    capped = total > cap
    if capped:
        total = cap
    if total and offset >= total:
        offset = ((total - 1) // limit) * limit
    rows = conn.execute(
        f"{select_sql} {from_sql} ORDER BY {order_sql} LIMIT ? OFFSET ?", list(params) + [limit, offset]
    ).fetchall()
    return rows, {"total": total, "total_capped": capped, "limit": limit, "offset": offset,
                  "truncated": len(rows) < total}


def _with_refs(conn, rows, kind):
    """
    Adds `ref` (the hit's object reference) and `parents` (its chain below the galaxy: sector, system, and for
    a moon its planet, each `{ref, kind, name}`) to search or locate rows, from the names the panel's own JOIN
    already selected plus one batched lookup of the sector names (NAV.9). `kind` is a reference kind, or `None`
    for the phenomenon panel, whose rows carry their own `type`.
    """
    sector_ids = {r["sector_id"] for r in rows if r.get("sector_id") is not None and kind != "sector"}
    names = {}
    if sector_ids:
        marks = ",".join("?" * len(sector_ids))
        names = {r["id"]: r["name"] for r in conn.execute(
            f"SELECT id, name FROM sectors WHERE id IN ({marks})", sorted(sector_ids)).fetchall()}
    for row in rows:
        row_kind = kind or row["type"]
        parents = []
        if kind != "sector" and row.get("sector_id") is not None:
            parents.append({"ref": object_ref.format("sector", row["sector_id"]), "kind": "sector",
                            "name": names.get(row["sector_id"])})
        system_id = row["id"] if kind == "system" else row.get("star_system_id")
        if kind in ("star", "planet", "moon", "belt") and system_id is not None:
            parents.append({"ref": object_ref.format("system", system_id), "kind": "system",
                            "name": row.get("system_name")})
        if kind == "moon" and row.get("rogue_planet_id") is not None:    # DB.94: a moon of a free planet
            parents.append({"ref": object_ref.format("rogue_planet", row["rogue_planet_id"]), "kind": "rogue_planet",
                            "name": row.get("planet_name")})
        elif kind == "moon":
            parents.append({"ref": object_ref.format("planet", row["planet_id"]), "kind": "planet",
                            "name": row.get("planet_name")})
        row["ref"] = object_ref.format(row_kind, row["id"])
        row["parents"] = parents
    return rows


def _search_result_sectors(conn, term, limit, offset):
    match, params = _search_match(conn, "name", term)
    rows, page = _search_page(
        conn, "SELECT id, name, edge_mpc", f"FROM sectors WHERE {match}", "name, id", params, limit, offset,
    )
    return {"rows": _with_refs(conn, [dict(r) for r in rows], "sector"), **page}


def _search_result_systems(conn, term, limit, offset, stars=None):
    clauses, params = [], []
    if term:
        match, params = _search_match(conn, "ss.name", term)
        clauses.append(match)
    star_clauses, star_params = _star_count_conditions(stars)
    clauses.extend(star_clauses)
    params = [*params, *star_params]
    rows, page = _search_page(
        conn,
        "SELECT ss.id, ss.name, ss.sector_id, ss.is_binary, ss.binary_configuration, ss.binary_type",
        f"FROM star_systems ss WHERE {' AND '.join(clauses)}",
        "ss.name, ss.id",
        params, limit, offset,
    )
    rows = _with_star_types(conn, rows)
    counts = _system_counts(conn, [r["id"] for r in rows])
    return {
        "rows": _with_refs(conn, [
            {
                "id": r["id"], "name": r["name"], "sector_id": r["sector_id"], "is_binary": r["is_binary"],
                "star_summary": _star_summary(r), **counts[r["id"]],
            }
            for r in rows
        ], "system"),
        **page,
    }


def _search_result_stars(conn, spectral_tags, luminosity_tags, term, limit, offset, size_range=None, ranges=None, sort=None):
    clauses, params = [], []
    if spectral_tags:
        clauses.append(f"SUBSTR(s.star_type, 1, 1) IN ({','.join('?' * len(spectral_tags))})")
        params.extend(sorted(spectral_tags))
    if luminosity_tags:
        clauses.append(f"s.yerkes_class IN ({','.join('?' * len(luminosity_tags))})")
        params.extend(sorted(luminosity_tags))
    _append_size_clause(clauses, params, "s.radius_km", size_range)
    search_values.clauses("star", "s", ranges, clauses, params)
    if term:
        match, match_params = _search_match(conn, "s.name", term)
        clauses.append(match)
        params.extend(match_params)
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    keys = search_values.used("star", ranges, sort)
    rows, page = _search_page(
        conn,
        "SELECT s.id, s.name, s.role, s.star_type, s.radius_km, s.star_system_id, ss.name AS system_name, ss.sector_id"
        + search_values.select_columns("star", "s", keys),
        f"FROM stars s JOIN star_systems ss ON ss.id = s.star_system_id WHERE 1=1{where}",
        search_values.order_by("star", "s", sort, "ss.name, s.name, s.id"),
        params, limit, offset,
    )
    return {"rows": _with_refs(conn, _with_values("star", rows, keys), "star"), **page}


def _tier_columns():
    """The stored habitability of a planet or moon (UX.100, UX.152): `phi4`, the six `tier_<domain>` columns and the four
    life scores."""
    return ("phi4", *(f"tier_{domain}" for domain in habitability_world.DOMAINS), *life_scores.FIELDS)


def _tiers(alias):
    """`_tier_columns` as a select list of table alias `alias` (none: bare column names)."""
    prefix = f"{alias}." if alias else ""
    return ", ".join(f"{prefix}{column}" for column in _tier_columns())


def _equipment_columns(alias):
    """The select list of a body's survival equipment levels (`equipment.FIELDS`) on the table `alias`."""
    return ", ".join(f"{alias}.{column}" for column in equipment.FIELDS)


def _with_equipment(row):
    """A search row with its equipment line (`equipment.line`, empty for a body that needs nothing), its items
    (`equipment_items`: `item`, `severity` and the PHI-6 `domain` of each, which the chips are drawn from, UX.154) and the
    levels it was worked out from taken out."""
    row["equipment"] = equipment.line(row)
    row["equipment_items"] = [{"item": item.name, "severity": item.severity, "domain": item.domain}
                              for item in equipment.items(row)]
    for column in equipment.FIELDS:
        row.pop(column, None)
    return row


def _with_values(entity, rows, keys):
    """The rows of a numeric-key panel as dicts, each with the `values` of the keys in use (`search_values.row_values`)."""
    out = []
    for row in rows:
        row = dict(row)
        if keys:
            row["values"] = search_values.row_values(entity, row, keys)
        out.append(row)
    return out


def _search_result_planets(conn, class_tags, body_tags, life_tags, term, limit, offset, size_range=None,
                           equipment_tags=(), tags=None, sort=None, ranges=None):
    clauses, params = [], []
    _equipment_clause("p", equipment_tags, clauses, params)
    _score_clauses("p", tags or {}, "", clauses)
    if class_tags:
        clauses.append(f"p.planet_class IN ({','.join('?' * len(class_tags))})")
        params.extend(sorted(class_tags))
    if body_tags:
        clauses.append(f"p.body_type IN ({','.join('?' * len(body_tags))})")
        params.extend(sorted(body_tags))
    _life_class_clause("p", life_tags, clauses, params)
    _append_size_clause(clauses, params, "p.radius_km", size_range)
    search_values.clauses("planet", "p", ranges, clauses, params)
    if term:
        match, match_params = _search_match(conn, "p.name", term)
        clauses.append(match)
        params.extend(match_params)
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    keys = search_values.used("planet", ranges, sort)
    rows, page = _search_page(
        conn,
        f"""
        SELECT p.id, p.name, p.planet_class, p.body_type, CASE WHEN p.life_stage IS NULL THEN NULL ELSE CONCAT(p.life_class, COALESCE(p.life_class_variant, '')) END AS life_class_id, {_equipment_columns('p')}, {_tiers("p")}, p.radius_km,
               p.star_system_id, ss.name AS system_name, ss.sector_id{search_values.select_columns("planet", "p", keys)}
        """,
        f"FROM planets p JOIN star_systems ss ON ss.id = p.star_system_id WHERE 1=1{where}",
        search_values.order_by("planet", "p", sort, _score_order("p", sort, "ss.name, ss.id, p.orbital_index, p.id")),
        params, limit, offset,
    )
    return {"rows": _with_refs(conn, [_with_equipment(r) for r in _with_values("planet", rows, keys)], "planet"), **page}


def _search_result_moons(conn, class_tags, body_tags, life_tags, term, limit, offset, size_range=None,
                         equipment_tags=(), tags=None, sort=None, ranges=None):
    clauses, params = [], []
    _equipment_clause("m", equipment_tags, clauses, params)
    _score_clauses("m", tags or {}, "moon_", clauses)
    if class_tags:
        clauses.append(f"m.planet_class IN ({','.join('?' * len(class_tags))})")
        params.extend(sorted(class_tags))
    if body_tags:
        clauses.append(f"m.body_type IN ({','.join('?' * len(body_tags))})")
        params.extend(sorted(body_tags))
    _life_class_clause("m", life_tags, clauses, params)
    _append_size_clause(clauses, params, "m.radius_km", size_range)
    search_values.clauses("moon", "m", ranges, clauses, params)
    if term:
        match, match_params = _search_match(conn, "m.name", term)
        clauses.append(match)
        params.extend(match_params)
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    keys = search_values.used("moon", ranges, sort)
    rows, page = _search_page(
        conn,
        f"""
        SELECT m.id, m.name, m.planet_class, m.body_type, CASE WHEN m.life_stage IS NULL THEN NULL ELSE CONCAT(m.life_class, COALESCE(m.life_class_variant, '')) END AS life_class_id, {_equipment_columns('m')}, {_tiers("m")}, m.radius_km,
               m.planet_id, m.rogue_planet_id, COALESCE(p.name, rp.name) AS planet_name, m.star_system_id,
               ss.name AS system_name, COALESCE(ss.sector_id, rp.sector_id) AS sector_id{search_values.select_columns("moon", "m", keys)}
        """,
        f"""
        FROM moons m
        LEFT JOIN planets p ON p.id = m.planet_id
        LEFT JOIN rogue_planets rp ON rp.id = m.rogue_planet_id
        LEFT JOIN star_systems ss ON ss.id = m.star_system_id
        WHERE 1=1{where}
        """,
        search_values.order_by("moon", "m", sort, _score_order("m", sort, "ss.name, ss.id, rp.name, rp.id, p.orbital_index, m.orbital_index, m.id")),
        params, limit, offset,
    )
    return {"rows": _with_refs(conn, [_with_equipment(r) for r in _with_values("moon", rows, keys)], "moon"), **page}


def _search_result_belts(conn, density_tags, limit, offset):
    clauses, params = [], []
    if density_tags:
        clauses.append(f"ab.density IN ({','.join('?' * len(density_tags))})")
        params.extend(sorted(density_tags))
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    rows, page = _search_page(
        conn,
        "SELECT ab.id, ab.density, ab.composition_summary, ab.star_system_id, ss.name AS system_name, ss.sector_id",
        f"FROM asteroid_belts ab JOIN star_systems ss ON ss.id = ab.star_system_id WHERE 1=1{where}",
        "ss.name, ss.id, ab.orbital_index, ab.id",
        params, limit, offset,
    )
    return {"rows": _with_refs(conn, [dict(r) for r in rows], "belt"), **page}


SEARCH_FACET_CACHE_SECONDS = 600
"""int: Longest the search page's facet counts and name lists are reused
(PERF.15) even when nothing they depend on looks changed -- a backstop
for edits `_search_content_key` can't see."""

_facet_cache = {}
_facet_cache_lock = threading.Lock()


def _search_content_key(conn):
    """One cheap row (index-only maxima and the deletion epoch) that changes
    whenever a sector or system is added, edited or deleted -- what the
    facet counts depend on."""
    row = conn.execute(
        "SELECT (SELECT COALESCE(MAX(id), 0) FROM sectors) AS sector_max_id, "
        "(SELECT MAX(modified_at) FROM sectors) AS sector_modified, "
        "(SELECT COALESCE(MAX(id), 0) FROM star_systems) AS system_max_id, "
        "(SELECT MAX(modified_at) FROM star_systems) AS system_modified"
    ).fetchone()
    return (str(sector_deletions(conn)),) + tuple(str(row[column]) for column in (
        "sector_max_id", "sector_modified", "system_max_id", "system_modified"))


def _search_facets_cached(conn):
    """
    `(facet_defs, autocomplete)` for `search` -- the 12 facet queries
    (counts over every star, planet, moon, belt and phenomenon) and 5
    name lists, cached per database until `_search_content_key` changes
    or `SEARCH_FACET_CACHE_SECONDS` pass (PERF.15), so a search request
    runs one cheap query for them instead of 17 scans.
    """
    config = getattr(conn, "_config", None)
    database = config._key() if config is not None else None
    key = _search_content_key(conn)
    now = time.monotonic()
    with _facet_cache_lock:
        cached = _facet_cache.get(database)
    if cached is not None and cached[0] == key and now - cached[1] < SEARCH_FACET_CACHE_SECONDS:
        return cached[2]
    facet_defs = (
        ("type", _search_facet_type(conn)),
        ("spectral", _search_facet_spectral(conn)),
        ("luminosity", _search_facet_luminosity(conn)),
        ("class", _search_facet_class(conn)),
        ("body", _search_facet_body(conn)),
        ("life", _search_facet_life(conn)),
        ("equipment", _search_facet_equipment(conn)),
        ("moon_class", _search_facet_moon_class(conn)),
        ("moon_body", _search_facet_moon_body(conn)),
        ("moon_life", _search_facet_moon_life(conn)),
        ("moon_equipment", _search_facet_moon_equipment(conn)),
        *_score_facet_defs(conn),
        ("density", _search_facet_density(conn)),
        ("phenomenon", _search_facet_phenomenon(conn)),
        ("phenomenon_class", _search_facet_phenomenon_class(conn)),
    )
    autocomplete = {
        "sectors": _search_name_list(conn, "sectors"),
        "systems": _search_name_list(conn, "star_systems"),
        "stars": _search_name_list(conn, "stars"),
        "planets": _search_name_list(conn, "planets"),
        "moons": _search_name_list(conn, "moons"),
    }

    value = (facet_defs, autocomplete)
    if database is not None:
        with _facet_cache_lock:
            _facet_cache[database] = (key, now, value)
    return value


def search(conn, texts, tags, sizes=None, limit=SEARCH_RESULT_LIMIT, offsets=None, only=None, sorts=None, stars=None, ranges=None):
    """
    Runs the faceted search behind `GET /api/search`/`html/search.py`:
    the same click-to-filter attribute tags (object type; star spectral/
    luminosity class; planet/moon class, body type, supported life
    chemistry; asteroid belt density) plus a per-entity name search this
    project's search page has always offered -- ported from
    `html/search.py`'s previous direct-SQL implementation (same queries,
    same "which result panels actually have a reason to run" logic: a
    panel only appears when one of its own tags/name field/size range is
    active, or its object type is explicitly selected), plus a min/max
    size (`radius_km`) range filter for stars/planets/moons each
    (`sizes`) -- there's no discrete set of values to offer as a facet
    for a continuous quantity like size, so it's a range, not a tag.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        texts (dict): `{"sector_q", "system_q", "star_q", "planet_q",
            "moon_q"}` -> search term (`""`/absent for "not searched").
        tags (dict): `{facet: set(value, ...)}`, one entry per name in
            `SEARCH_TAG_FACETS` -- an absent/empty set means no active
            filter for that facet.
        sizes (dict, optional): `{"star", "planet", "moon"} -> (min_km,
            max_km)`, each bound `None` for "unbounded" -- an absent key
            (or `None` altogether) means no size filter for that entity.
        limit (int): Rows per result panel page.
        offsets (dict, optional): `{panel: offset}` for any of
            `SEARCH_RESULT_PANELS` -- each panel pages independently; an
            absent panel starts at 0.
        only (iterable, optional): Run just these panels (a table scrolling
            one panel asks for that panel alone); the others come back `None`.
        stars (tuple, optional): The stars filter `(least, most)` of `parse_star_count` (UX.239): the Systems panel runs for it (with or without a
            name) and holds only the systems with that many stars in the `stars` table.
        ranges (dict, optional): `{"star"/"planet"/"moon": {key: (least, most, unit)}}` of `search_values.ranges_from_args` (UX.226): each
            bound `None` for "unbounded", `unit` the unit key both are in. A panel with a range runs for it, and its rows carry the `values`
            of the keys in use. `sorts` also takes a key of `search_values.KEYS` (`gravity` smallest first, `-gravity` largest first; no
            value last) for `star`, `planet` and `moon`.

    Returns:
        dict: `facets` (`{facet: [{"value","label","count","tooltip"}, ...]}`,
            one list per `SEARCH_TAG_FACETS` name), `autocomplete`
            (`sectors`/`systems`/`stars`/`planets`/`moons` -> list of
            distinct names), `facet_labels` (`{"facet:value": label}`,
            for rendering an active-filter chip without a second lookup),
            `values` (the numeric keys of `search_values.KEYS` with their units, UX.226), `active_ranges` (the ranges of
            this request, one `{entity, key, label, least, most, unit, symbol}` each), and `results` (`sectors`/`systems`/`stars`/`planets`/`moons`/
            `belts` -> `{"rows": [...], "total", "limit", "offset",
            "truncated"}`, one page of that panel's matches -- `truncated`
            meaning `rows` isn't every match -- or `None` for a panel
            with no reason to run).
    """
    sizes = sizes or {}
    star_size, planet_size, moon_size = sizes.get("star"), sizes.get("planet"), sizes.get("moon")
    ranges = ranges or {}

    spectral_tags, luminosity_tags = tags.get("spectral", set()), tags.get("luminosity", set())
    class_tags, body_tags, life_tags = tags.get("class", set()), tags.get("body", set()), tags.get("life", set())
    moon_class_tags, moon_body_tags = tags.get("moon_class", set()), tags.get("moon_body", set())
    moon_life_tags = tags.get("moon_life", set())
    equipment_tags, moon_equipment_tags = tags.get("equipment", set()), tags.get("moon_equipment", set())
    density_tags = tags.get("density", set())
    type_tags = tags.get("type", set())

    facet_defs, autocomplete = _search_facets_cached(conn)
    facets = {name: options for name, options in facet_defs}
    facet_labels = {
        f"{name}:{opt['value']}": opt["label"]
        for name, options in facet_defs for opt in options
    }

    star_has_reason = bool(spectral_tags or luminosity_tags or texts.get("star_q") or star_size or ranges.get("star"))
    planet_score = any(tags.get(name) for name in score_filters.FACETS)
    moon_score = any(tags.get(f"moon_{name}") for name in score_filters.FACETS)
    planet_has_reason = bool(class_tags or body_tags or life_tags or equipment_tags or planet_score or texts.get("planet_q") or planet_size or ranges.get("planet"))
    moon_has_reason = bool(moon_class_tags or moon_body_tags or moon_life_tags or moon_equipment_tags or moon_score or texts.get("moon_q") or moon_size or ranges.get("moon"))
    belt_has_reason = bool(density_tags)

    if type_tags:
        stars_included = "star" in type_tags
        planets_included = "planet" in type_tags
        moons_included = "moon" in type_tags
        belts_included = "belt" in type_tags
    else:
        stars_included = star_has_reason
        planets_included = planet_has_reason
        moons_included = moon_has_reason
        belts_included = belt_has_reason

    offsets = offsets or {}

    def _page(panel):
        return limit, offsets.get(panel, 0)

    results = {panel: None for panel in SEARCH_RESULT_PANELS}
    wanted = set(SEARCH_RESULT_PANELS if only is None else only)
    if texts.get("sector_q") and "sectors" in wanted:
        results["sectors"] = _search_result_sectors(conn, texts["sector_q"], *_page("sectors"))
    if (texts.get("system_q") or stars is not None) and "systems" in wanted:
        results["systems"] = _search_result_systems(conn, texts.get("system_q", ""), *_page("systems"), stars=stars)
    if stars_included and "stars" in wanted:
        results["stars"] = _search_result_stars(
            conn, spectral_tags, luminosity_tags, texts.get("star_q", ""), *_page("stars"), size_range=star_size,
            ranges=ranges.get("star"), sort=(sorts or {}).get("star"),
        )
    if planets_included and "planets" in wanted:
        results["planets"] = _search_result_planets(
            conn, class_tags, body_tags, life_tags, texts.get("planet_q", ""), *_page("planets"),
            size_range=planet_size, equipment_tags=equipment_tags,
            tags=tags, sort=(sorts or {}).get("planet"), ranges=ranges.get("planet"),
        )
    if moons_included and "moons" in wanted:
        results["moons"] = _search_result_moons(
            conn, moon_class_tags, moon_body_tags, moon_life_tags, texts.get("moon_q", ""), *_page("moons"),
            size_range=moon_size, equipment_tags=moon_equipment_tags,
            tags=tags, sort=(sorts or {}).get("moon"), ranges=ranges.get("moon"),
        )
    if belts_included and "belts" in wanted:
        results["belts"] = _search_result_belts(conn, density_tags, *_page("belts"))
    phenomenon_tags, phenomenon_class_tags = tags.get("phenomenon", set()), tags.get("phenomenon_class", set())
    if (phenomenon_tags or phenomenon_class_tags) and "phenomena" in wanted:
        results["phenomena"] = _search_result_phenomena(
            conn, phenomenon_tags, phenomenon_class_tags, *_page("phenomena"))

    return {"facets": facets, "autocomplete": autocomplete, "facet_labels": facet_labels, "results": results,
            "values": search_values.describe(), "active_ranges": search_values.active(ranges)}
