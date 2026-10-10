# planetgen.db.query

"""
Read queries over the planetGen database (src/planetgen/db/schema.sql).

A thin read-only front end over the tables `sectorGen.py`/`systemGen.py`
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

import pymysql

from planetgen.db.store import escape_like, get_connection, get_galaxy_bounds, get_galaxy_shape, surrounding_cloud
from planetgen.physics import constants
from planetgen.physics.habitability_world import EQUIPMENT_LABELS, EQUIPMENT_NAMES
from planetgen import tuning
from planetgen.generation.star import compressed_heliosphere_radius
from planetgen.generation.bright_stars import MPC_PER_PC, POPULATIONS as BRIGHT_STAR_POPULATIONS
from planetgen._version import __version__
from planetgen.galaxy.density import shape_with_terms
from planetgen.galaxy.geometry import (
    galaxy_to_local_pc, layer_index_at, neighbor_addresses, provisional_sector_designation, ring_index_at,
    ring_sector_count, sector_cell_vertices_pc, sector_position_pc,
)
from planetgen.galaxy import keepout
from planetgen.galaxy import objectref as object_ref
from planetgen.galaxy import object_uid
from planetgen.names import naming_key
from planetgen.galaxy.sector import classify_octant
from planetgen.galaxy.viewport import (
    TILE_MAX_LEVEL,
    TILE_ROOT_EDGE_PC,
    parse_tile_key,
    planned_slots_in_tile,
    tile_bounds_pc,
    tile_keys_containing,
)
from planetgen.galaxy.drill import (
    DRILL_TOP, DrillBlock, drill_chain_of, drill_wedge_count, format_drill_key, parse_drill_key,
)
from planetgen.db.corridor import positions_near_segment, unknown_space_flags
from planetgen.galaxy.nav_graph import build_route_graph, shortest_path
from planetgen.galaxy.navigation import (
    FRAME_GALACTIC, FRAME_SECTOR, FRAME_SYSTEM, course_between, fold_travel_times, warp_travel_times,
)
from planetgen.physics.constants import SPECTRAL_CLASS_COLORS
from planetgen.generation.evolution import life_stage_from_paragraphs
from planetgen.tuning import (
    DEFAULT_SECTOR_EDGE_LY, HABITABLE_PLANET_CLASSES, NAV_ADJACENCY_K, NAV_CORRIDOR_FRACTION, NAV_CORRIDOR_MAX_LY,
    NAV_CORRIDOR_MAX_SYSTEMS, NAV_CORRIDOR_MIN_LY, NAV_CORRIDOR_START_MAX_LY, NAV_ISLAND_LINKS, PLANET_CLASSES,
)
from planetgen.physics.units import ly_to_pc, milliparsecs_to_ly, mpc_to_pc, pc_to_ly


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
    "name": "sec.name", "systems": "system_count", "density": "(COUNT(ss.id) / POW(sec.edge_mpc, 3))",
    "position": "quadrant", "distance": "sec.galactic_radius_pc",
}
"""dict: The sort keys `list_sectors` accepts (the Sectors table's column
keys) -> the SQL they order by. A sector with no value (unplaced) sorts last
either way."""

SECTOR_QUADRANTS = ("I", "II", "III", "IV")
"""tuple[str]: The sector Quadrant labels `list_sectors` filters by, plus
`"unplaced"` for a sector with no galaxy position."""

_SECTOR_QUADRANT_SQL = (
    "(CASE WHEN sec.center_x_pc IS NULL THEN 'unplaced' "
    "WHEN sec.center_y_pc >= 0 THEN (CASE WHEN sec.center_x_pc >= 0 THEN 'I' ELSE 'II' END) "
    "ELSE (CASE WHEN sec.center_x_pc < 0 THEN 'III' ELSE 'IV' END) END)"
)
"""str: SQL for a sector's Quadrant label, the same bands as
`planetgen.galaxy.geometry.sector_quadrant`."""


def _sector_quadrant_filter(quadrants):
    """`(where_sql, params)` keeping sectors in any of `quadrants`."""
    if not quadrants:
        return "", []
    return f" WHERE {_SECTOR_QUADRANT_SQL} IN ({', '.join('?' for _ in quadrants)})", list(quadrants)


def list_sectors(conn, limit=None, offset=None, sort=None, descending=False, quadrants=()):
    """
    Returns every sector, with its edge length (converted to light-years)
    and how many systems it contains, nearest the galactic core first
    (`galactic_radius_pc`); sectors never placed in a galaxy have no
    distance and come last, by name.

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

    Returns:
        list[dict]: One row per sector, with `id`, `name`,
                           `edge_ly`, `system_count`,
                           `galactic_radius_pc`/`galactic_radius_ly`
                           (`None` if unplaced).
    """
    if sort is not None and sort not in SECTOR_SORTS:
        raise ValueError(f"unknown sector sort {sort!r}")
    where, params = _sector_quadrant_filter(quadrants)
    order = "sec.galactic_radius_pc IS NULL, sec.galactic_radius_pc, sec.name, sec.id"
    if sort is not None:
        column = SECTOR_SORTS[sort]
        unplaced_last = "sec.center_x_pc IS NULL, " if sort in ("density", "distance", "position") else ""
        order = f"{unplaced_last}{column} {'DESC' if descending else 'ASC'}, " + order
    query = f"""
        SELECT sec.id, sec.name, sec.edge_mpc, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.galactic_radius_pc, sec.ring_index, sec.layer_index, sec.ring_slot_index,
               COUNT(ss.id) AS system_count, {_SECTOR_QUADRANT_SQL} AS quadrant
        FROM sectors sec
        LEFT JOIN star_systems ss ON ss.sector_id = sec.id{where}
        -- Every selected column, not just the key: MariaDB's
        -- ONLY_FULL_GROUP_BY doesn't see columns that depend on sec.id.
        GROUP BY sec.id, sec.name, sec.edge_mpc, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
                 sec.galactic_radius_pc, sec.ring_index, sec.layer_index, sec.ring_slot_index
        ORDER BY {order}
        """
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params.extend([limit, offset or 0])

    rows = conn.execute(query, params).fetchall()
    return [
        {
            "id": r["id"], "name": r["name"], "edge_mpc": r["edge_mpc"], "edge_ly": milliparsecs_to_ly(r["edge_mpc"]),
            "system_count": r["system_count"],
            "center_x_pc": r["center_x_pc"], "center_y_pc": r["center_y_pc"], "center_z_pc": r["center_z_pc"],
            "galactic_radius_pc": r["galactic_radius_pc"],
            "galactic_radius_ly": (
                pc_to_ly(r["galactic_radius_pc"]) if r["galactic_radius_pc"] is not None else None
            ),
            "ring_index": r["ring_index"], "layer_index": r["layer_index"],
            "ring_slot_index": r["ring_slot_index"],
            "placed": r["center_x_pc"] is not None,
        }
        for r in rows
    ]


def count_sectors(conn, quadrants=()):
    """
    Returns the total number of sectors (that pass the same `quadrants`
    filter as `list_sectors`), ignoring any pagination --
    the denominator `list_sectors(conn, limit=...)` callers (the API's
    `/api/sectors`) need to report how many pages exist.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.

    Returns:
        int: Total sector count.
    """
    where, params = _sector_quadrant_filter(quadrants)
    return conn.execute(f"SELECT COUNT(*) AS n FROM sectors sec{where}", params).fetchone()["n"]


def sectors_facets(conn, quadrants=()):
    """
    The option counts for the Sectors table's Quadrant menu: how many sectors
    each Quadrant (and `"unplaced"`) holds. The menu ignores its own filter,
    so choosing one Quadrant still shows the others' counts.

    Returns:
        dict: `{"quadrant": [{"value", "count"}]}` in Quadrant order, only
            those with sectors.
    """
    counts = {row["quadrant"]: row["n"] for row in conn.execute(
        f"SELECT {_SECTOR_QUADRANT_SQL} AS quadrant, COUNT(*) AS n FROM sectors sec GROUP BY quadrant").fetchall()}
    return {"quadrant": [{"value": value, "count": counts[value]}
                         for value in (*SECTOR_QUADRANTS, "unplaced") if counts.get(value)]}


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


SYSTEM_SORTS = {"name": "ss.name", "sector": "sec.name", "octant": "ss.quadrant", "binary": "ss.is_binary"}
"""dict: The sort keys `list_systems` accepts (the Systems tables' column
keys) -> the SQL they order by. Ties fall back to name, then id."""


def _systems_filter_clause(star_type_prefix, sector_id, binary=None, octants=(), in_sector=None):
    """
    Builds the shared `JOIN`/`WHERE`/params fragment `list_systems` and
    `count_systems` both need -- factored out so the count query can't
    silently drift out of sync with what the listing query actually
    matches.

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

    where_sql = (" WHERE " + " AND ".join(conditions)) if conditions else ""
    return join_sql, where_sql, params


def list_systems(conn, star_type_prefix=None, sector_id=None, limit=None, offset=None, sort="name",
                 descending=False, binary=None, octants=(), in_sector=None):
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

    Returns:
        list[dict]: One row per matching system, with `id`, `name`,
                           `sector_id`, `sector_name`, `quadrant` (its
                           octant in the sector), `is_binary`, and `star_summary`
                           (the single star's `star_type`, or a binary's
                           `binary_type` -- what `html/browse.py`/
                           `html/search.py` show as a system's "Star type"
                           column).
    """
    if sort not in SYSTEM_SORTS:
        raise ValueError(f"unknown system sort {sort!r}")
    join_sql, where_sql, params = _systems_filter_clause(star_type_prefix, sector_id, binary, octants, in_sector)
    column = SYSTEM_SORTS[sort]
    direction = "DESC" if descending else "ASC"
    nulls_last = f"{column} IS NULL, " if sort in ("sector", "octant") else ""
    query = f"""
        SELECT DISTINCT ss.id, ss.name, ss.sector_id, ss.quadrant, ss.is_binary, ss.binary_configuration,
               ss.binary_type, sec.name AS sector_sort
        FROM star_systems ss LEFT JOIN sectors sec ON sec.id = ss.sector_id{join_sql}{where_sql}
        ORDER BY {nulls_last}{column} {direction}, ss.name, ss.id
        """

    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params = params + [limit, offset or 0]

    rows = _with_star_types(conn, conn.execute(query, params).fetchall(), with_sector_name=True)
    return [
        {
            "id": r["id"], "name": r["name"], "sector_id": r["sector_id"], "sector_name": r["sector_name"],
            "quadrant": r["quadrant"], "is_binary": r["is_binary"], "star_summary": _star_summary(r),
        }
        for r in rows
    ]


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


def count_systems(conn, star_type_prefix=None, sector_id=None, binary=None, octants=(), in_sector=None):
    """
    Returns the total number of systems matching the same filters
    `list_systems` accepts, ignoring any pagination -- the denominator
    `list_systems(conn, limit=...)` callers (the API's `/api/systems`)
    need to report how many pages exist.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        star_type_prefix (str, optional): Same meaning as `list_systems`.
        sector_id (int or NO_SECTOR, optional): Same meaning as `list_systems`.

    Returns:
        int: Total matching system count.
    """
    join_sql, where_sql, params = _systems_filter_clause(star_type_prefix, sector_id, binary, octants, in_sector)
    query = f"SELECT COUNT(DISTINCT ss.id) AS n FROM star_systems ss{join_sql}{where_sql}"
    return conn.execute(query, params).fetchone()["n"]


def systems_facets(conn, sector_id=None, binary=None, octants=(), in_sector=None):
    """
    The option counts for the Systems tables' filter menus: where a system is
    (`placement`: `"sector"` or `"standalone"`), whether it is `binary`
    (`"yes"`/`"no"`) and its sector `octant`. Each menu's counts apply every
    filter except its own, so the menus narrow one another.

    Returns:
        dict: `{"placement"|"binary"|"octant": [{"value", "count"}]}`,
            options with no systems left out.
    """
    def counts(select, **filters):
        base = {"binary": binary, "octants": octants, "in_sector": in_sector, **filters}
        _join, where, params = _systems_filter_clause(None, sector_id, **base)
        return {row["value"]: row["n"] for row in conn.execute(
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


def list_planets(conn, planet_class=None, min_radius_km=None, max_radius_km=None,
                  sector_id=None, system_id=None, limit=None, offset=None):
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
            unless `limit` is also given.

    Returns:
        list[dict]: One row per matching planet, with `id`, `name`,
            `planet_class`, `body_type`, `radius_km`, `star_system_id`,
            `system_name`.
    """
    where_sql, params = _body_filter_clause("p", planet_class, min_radius_km, max_radius_km, sector_id, system_id)
    query = f"""
        SELECT p.id, p.name, p.planet_class, p.body_type, p.radius_km, p.star_system_id, ss.name AS system_name
        FROM planets p
        JOIN star_systems ss ON ss.id = p.star_system_id{where_sql}
        ORDER BY ss.name, p.orbital_index
        """
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params = params + [limit, offset or 0]

    rows = conn.execute(query, params).fetchall()
    return [dict(r) for r in rows]


def list_moons(conn, planet_class=None, min_radius_km=None, max_radius_km=None,
                sector_id=None, system_id=None, limit=None, offset=None):
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
            unless `limit` is also given.

    Returns:
        list[dict]: One row per matching moon, with `id`, `name`,
            `planet_class`, `body_type`, `radius_km`, `star_system_id`,
            `system_name`, `planet_name` (its parent planet's name).
    """
    where_sql, params = _body_filter_clause("m", planet_class, min_radius_km, max_radius_km, sector_id, system_id)
    query = f"""
        SELECT m.id, m.name, m.planet_class, m.body_type, m.radius_km, m.star_system_id,
               ss.name AS system_name, p.name AS planet_name
        FROM moons m
        JOIN planets p ON p.id = m.planet_id
        JOIN star_systems ss ON ss.id = m.star_system_id{where_sql}
        ORDER BY ss.name, p.orbital_index, m.orbital_index
        """
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params = params + [limit, offset or 0]

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
    return max((math.dist(positions[a], positions[b]) for a, b in zip(path, path[1:])), default=0.0)


def _galaxy_route(conn, origin_position, destination_position, from_key, to_key, adjacency_k):
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

    Returns:
        tuple: `(positions, found)`: the positions the graph was built from
            (`{node_id: (x, y, z)}`, light-years) and `shortest_path`'s answer.
    """
    direct_ly = math.dist(origin_position, destination_position)
    half_width = min(max(direct_ly * NAV_CORRIDOR_FRACTION, NAV_CORRIDOR_MIN_LY), NAV_CORRIDOR_START_MAX_LY)
    while True:
        positions = positions_near_segment(conn, origin_position, destination_position, half_width)
        # The ends are nodes whatever they are: a phenomenon is no star_systems
        # row, and a system's own position is the one NAV measured.
        positions[from_key] = origin_position
        positions[to_key] = destination_position
        graph = build_route_graph(positions, adjacency_k, NAV_ISLAND_LINKS)
        found = shortest_path(graph, from_key, to_key, positions)
        settled = (found is not None and _longest_hop(found[0], positions) <= half_width / 2.0)
        if settled or half_width >= NAV_CORRIDOR_MAX_LY or len(positions) > NAV_CORRIDOR_MAX_SYSTEMS:
            return positions, found
        half_width = min(half_width * 2.0, NAV_CORRIDOR_MAX_LY)


def nav_between(conn, from_id, to_id, adjacency_k=NAV_ADJACENCY_K,
                 from_kind="system", to_kind="system", from_type=None, to_type=None):
    """
    Resolves full NAV information between two endpoints -- each either a
    star system or a standalone phenomenon (nebula/asteroid field/black
    hole/neutron star) -- a direct course (distance/bearing/mark/warp and fold
    travel times, from `planetgen.galaxy.navigation`) plus an optimal route
    via adjacent systems (`planetgen.galaxy.nav_graph`), or raises if NAV
    doesn't apply to this pair.

    NAV availability rules (see docs/api.md's NAV section for the
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
            "distance_ly", "unknown_space"}], "longest_hop_ly": float}` (one entry per id in
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
        if found is not None:
            path, distance_ly = found
            route_positions = {node_id: positions[node_id] for node_id in path}
            hops = [{"from": a, "to": b, "distance_ly": math.dist(positions[a], positions[b]),
                     "unknown_space": False} for a, b in zip(path, path[1:])]
            if galaxy_frame is not None:
                flags = unknown_space_flags(conn, [galaxy_frame[node_id] for node_id in path])
                for hop, flag in zip(hops, flags):
                    hop["unknown_space"] = flag
                if sector_offset is not None:
                    route_positions = {node_id: tuple(c - o for c, o in zip(position, sector_offset))
                                       for node_id, position in route_positions.items()}
            route = {
                "path": path,
                "distance_ly": distance_ly,
                "positions": route_positions,
                "hops": hops,
                "longest_hop_ly": max((hop["distance_ly"] for hop in hops), default=0.0),
            }

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
in `docs/design/navigation-frames.md`)."""

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
    return (heliopause_au if heliopause_au is not None else NAV_DEFAULT_HELIOPAUSE_AU) * constants.AU_TO_KM


def _km_to_ly(km):
    return km * 1000.0 / constants.LY_TO_M


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


def nav_course(conn, from_ref, to_ref, adjacency_k=NAV_ADJACENCY_K):
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
                         from_type=from_type, to_type=to_type)

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
    """A `sectors` row's `(ring_index, layer_index, ring_slot_index)`, or
    `None` when it has no grid address (never placed, or hand-placed)."""
    address = (row["ring_index"], row["layer_index"], row["ring_slot_index"])
    return None if any(value is None for value in address) else address


def sector_neighbors(conn, sector):
    """
    Every cell sharing a face with a galaxy-placed sector
    (`galaxyGeometry.neighbor_addresses`: the slots either side, the
    layers above and below, and the overlapping slots in the rings inside
    and outside), each tagged with whether a real `sectors` row already
    exists there. Drives the Sector Map's (`planetgen/web/maps/starmap.py`)
    neighboring-sector indicators -- an existing neighbor links straight
    to it; a not-yet-generated one shows its address so it can be fed to
    `planetgen galaxy --ring I --layer J --slot K`.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        sector (dict or Row): This sector's own row -- needs
            `ring_index`, `layer_index`, `ring_slot_index`, and `edge_mpc`.

    Returns:
        list[dict]: `ring_index`, `layer_index`, `ring_slot_index`,
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

    clauses = " OR ".join(["(ring_index = ? AND layer_index = ? AND ring_slot_index = ?)"] * len(addresses))
    params = [value for neighbor in addresses for value in neighbor]
    rows = conn.execute(
        f"SELECT id, name, ring_index, layer_index, ring_slot_index FROM sectors WHERE {clauses}", tuple(params),
    ).fetchall()
    existing = {(r["ring_index"], r["layer_index"], r["ring_slot_index"]): r for r in rows}

    results = []
    for neighbor in addresses:
        position = sector_position_pc(*neighbor, edge_pc)
        match = existing.get(neighbor)
        results.append({
            "ring_index": neighbor[0], "layer_index": neighbor[1], "ring_slot_index": neighbor[2],
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
           COALESCE(st.name, p.name, m.name, af.name, sec.name) AS host_name
    FROM facilities f
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
    return facility


def facility_detail(conn, facility_id):
    """One facility (schema v42), or `None`."""
    row = conn.execute(_FACILITY_SELECT + " WHERE f.id = ?", (facility_id,)).fetchone()
    return None if row is None else _facility_dict(row)


def facilities_for_system(conn, system_id):
    """Every facility on a host inside one star system (its stars,
    planets, moons and belts), by name."""
    rows = conn.execute(_FACILITY_SELECT + " WHERE f.star_system_id = ? ORDER BY f.name, f.id",
                        (system_id,)).fetchall()
    return [_facility_dict(row) for row in rows]


def facilities_in_sector(conn, sector_id):
    """The facilities a sector holds outside its systems: stand-alone ones
    parked in it, and those on its asteroid fields."""
    rows = conn.execute(
        _FACILITY_SELECT + " WHERE f.sector_id = ? OR af.sector_id = ? ORDER BY f.name, f.id",
        (sector_id, sector_id),
    ).fetchall()
    return [_facility_dict(row) for row in rows]


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
    its `inside_nebula_id`/`inside_remnant_id`: `_db.surrounding_cloud`.

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
    return surrounding_cloud(conn, row)


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
            `center_y_pc`/`center_z_pc`, `ring_index`, `layer_index`,
            `ring_slot_index`, `placed`, `system_count`, and `systems` (one entry per system
            placed in this sector, nearest the sector's center first,
            then any with no position by name: `id`, `name`, `quadrant`, `location`,
            `is_binary`, `binary_type`, `position_x_mpc`/`position_y_mpc`/
            `position_z_mpc`, `center_distance_ly` (`None` with no
            position), and `stars` -- 1 entry (single) or 2
            (primary then secondary), each `role`/`star_type`/
            `temperature_k`/`radius_km`/`luminosity_w`), `phenomena`
            (see `phenomena_near_sector`), `neighbors` (see
            `sector_neighbors`, `[]` for a sector with no galaxy
            placement), and `wiki_url` (`sectors.wiki_url` -- `None` if
            this sector has never been uploaded to, or manually linked
            to, a wiki page; see `schema.sql`'s "v22" header note), and
            `stats` (`sector_stats`, PERF.11: `None` off the grid).

    Raises:
        ValueError: If no such sector exists.
    """
    sector = conn.execute("SELECT * FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    if sector is None:
        raise ValueError(f"no sectors row with id {sector_id}")

    system_rows = conn.execute(
        """
        SELECT id, name, is_binary, quadrant, location, binary_type,
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
        "SELECT st.star_system_id, st.role, st.name, st.star_type, st.temperature_k, st.radius_km, st.luminosity_w"
        " FROM stars st JOIN star_systems ss ON ss.id = st.star_system_id WHERE ss.sector_id = ?"
        " ORDER BY st.star_system_id, CASE st.role WHEN 'secondary' THEN 1 ELSE 0 END, st.id",
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
        "ring_index": sector["ring_index"], "layer_index": sector["layer_index"],
        "ring_slot_index": sector["ring_slot_index"],
        "placed": sector["center_x_pc"] is not None,
        "system_count": len(systems),
        "star_count": star_count,
        "interstellar_debris_count": interstellar_debris_count(star_count),
        "systems": systems,
        "phenomena": phenomena_near_sector(conn, sector_id),
        "neighbors": sector_neighbors(conn, sector),
        "wiki_url": sector["wiki_url"],
        "stats": sector_stats(conn, sector["ring_index"], sector["layer_index"], sector["ring_slot_index"]),
    }


SECTOR_STATS_FIELDS = (
    "bright_level_sol", "relative_density", "expected_systems", "actual_systems", "actual_stars",
    "mean_temperature_k", "mean_luminosity_sol", "mean_age_gy", "total_luminosity_sol",
)
"""tuple: The `sector_stats` columns `sector_stats` reads (PERF.11)."""


def sector_stats(conn, ring_index, layer_index, ring_slot_index):
    """
    One grid sector's stored stats (PERF.11, GEN.44): its bright-star
    level, its expected density and, once filled, what it holds
    (`SECTOR_STATS_FIELDS`). `None` for a sector off the grid or with no
    row yet.
    """
    if ring_index is None or layer_index is None or ring_slot_index is None:
        return None
    row = conn.execute(
        f"SELECT {', '.join(SECTOR_STATS_FIELDS)} FROM sector_stats"
        " WHERE ring_index = ? AND layer_index = ? AND ring_slot_index = ?",
        (ring_index, layer_index, ring_slot_index),
    ).fetchone()
    return None if row is None else {field: row[field] for field in SECTOR_STATS_FIELDS}


def interstellar_debris_count(star_count):
    """
    About how many interstellar comets and planetesimals (meters to
    kilometers across) drift through a sector holding `star_count` stars:
    `INTERSTELLAR_DEBRIS_DENSITY_PC3` scales with stellar density like
    every other interstellar object, so it's a fixed number per star
    (~7e12). Far too many to be rows -- a figure for the sector page
    ("about 10^13 interstellar comets and planetesimals"); the
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


_PHENOMENON_CLASS_COLUMNS = {
    "nebulae": "nebula_class",
    "supernova_remnants": "remnant_class",
    "asteroid_fields": "field_class",
    "rogue_planets": "planet_class",
}
"""dict: The letter-class column (v38; a rogue planet's `PLANET_CLASSES`
letter since v47) of each `_PHENOMENON_TABLES` table that has one; every
other type's rows carry `class` as `None`."""


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
            point-like types), `class` (see `_PHENOMENON_CLASS_COLUMNS`),
            `x`/`y`/`z` (`center_x/y/z_pc`),
            `galactic_radius_pc`.
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
                   center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc
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
                "galactic_radius_pc": row["galactic_radius_pc"],
            })
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


def _scattered_page(conn, classes, limit, offset):
    """`limit` unbuilt scatter rows from `offset` among the `classes`, in id order, shaped like the Phenomena rows
    (`scattered` True, named "Uncharted <kind> <ring>.<layer>.<slot>")."""
    if not classes:
        return []
    by_class = {(entry["kind"], entry["subtype"]): entry for entry in classes}
    clauses, params = [], []
    for entry in classes:
        clauses.append("(kind = ? AND COALESCE(subtype, '') = ?)")
        params.extend([entry["kind"], entry["subtype"]])
    query = ("SELECT id, kind, COALESCE(subtype, '') AS subtype, ring_index, layer_index, ring_slot_index"
             f" FROM phenomenon_scatter WHERE built_at IS NULL AND ({' OR '.join(clauses)}) ORDER BY id")
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params.extend([limit, offset or 0])
    return [
        {"id": row["id"], "type": by_class[(row["kind"], row["subtype"])]["type"],
         "name": f"Uncharted {row['kind'].replace('-', ' ')} {row['ring_index']}.{row['layer_index']}."
                 f"{row['ring_slot_index']}",
         "descriptor": by_class[(row["kind"], row["subtype"])]["descriptor"], "radius_ly": 0.0,
         "sector_id": None, "sector_name": None, "placed": True, "scattered": True}
        for row in conn.execute(query, params).fetchall()
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


def list_phenomena(conn, limit=None, offset=None, sort="name", descending=False, types=(), descriptors=(),
                   placed=None):
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
            type and id so pages never overlap.
        descending (bool): Reverse the sort.
        types (iterable[str]): Keep only these types (any of).
        descriptors (iterable[str]): Keep only these descriptors (any of).
        placed (bool, optional): Keep only placed or only unplaced ones.

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
            `phenomenon_scatter` row, with no page; its `type` can also be
            `"hypervelocity_star"`).
    """
    if sort not in PHENOMENON_SORTS:
        raise ValueError(f"unknown phenomenon sort {sort!r}")
    where, params = _phenomena_where(types, descriptors, placed)
    column = PHENOMENON_SORTS[sort]
    built_total = conn.execute(
        f"SELECT COUNT(*) AS n FROM ({_phenomena_union()}) AS phenomena{where}", params).fetchone()["n"]
    query = (f"SELECT * FROM ({_phenomena_union()}) AS phenomena{where} "
             f"ORDER BY {column} IS NULL, {column} {'DESC' if descending else 'ASC'}, name, type, id")
    start = offset or 0
    built_limit = limit
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params.extend([limit, start])
    items = [
        {
            "id": row["id"], "type": row["type"], "name": row["name"],
            "descriptor": row["descriptor"], "radius_ly": row["radius_ly"],
            "sector_id": row["sector_id"], "sector_name": row["sector_name"],
            "placed": bool(row["placed"]), "scattered": bool(row["scattered"]),
        }
        for row in conn.execute(query, params).fetchall()
    ] if start < built_total or limit is None else []
    # The unbuilt scattered phenomena come after the built ones, in id order.
    room = None if limit is None else limit - len(items)
    if room is None or room > 0:
        items += _scattered_page(conn, _scattered_classes(conn, types, descriptors, placed), room,
                                 max(0, start - built_total))
    return items


def count_phenomena(conn, types=(), descriptors=(), placed=None):
    """
    Returns the number of exotic phenomena across every type in
    `_PHENOMENON_TABLES` that pass the same filters as `list_phenomena`,
    ignoring any pagination -- the denominator `list_phenomena(conn,
    limit=...)` callers (the API's `/api/phenomena`) need to report how many
    pages exist.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.

    Returns:
        int: Phenomenon count.
    """
    where, params = _phenomena_where(types, descriptors, placed)
    built = conn.execute(
        f"SELECT COUNT(*) AS n FROM ({_phenomena_union()}) AS phenomena{where}", params).fetchone()["n"]
    return built + sum(entry["unbuilt"] for entry in _scattered_classes(conn, types, descriptors, placed))


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
    facets = {}
    for column, own in (("type", "types"), ("descriptor", "descriptors")):
        filters = {"types": types, "descriptors": descriptors, own: ()}
        where, params = _phenomena_where(filters["types"], filters["descriptors"], placed)
        rows = conn.execute(
            f"SELECT {column} AS value, COUNT(*) AS n FROM ({_phenomena_union()}) AS phenomena{where} "
            f"GROUP BY {column} ORDER BY {column}", params).fetchall()
        merged = {row["value"]: row["n"] for row in rows}
        for entry in _scattered_classes(conn, types, descriptors, placed, ignore=(own,)):
            merged[entry[column]] = merged.get(entry[column], 0) + entry["unbuilt"]
        facets[column] = [{"value": value, "count": merged[value]} for value in sorted(merged)]
    return facets


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
            body["mass_solar"] * constants.SOLAR_MASS_TO_KG, body["galactic_radius_pc"], body["radius_km"])
    if kind in ("quasar", "rogue_planet"):
        table = _PHENOMENON_TYPE_TO_TABLE[kind]
        body = row(table, "black_hole_mass_solar * " + str(constants.SOLAR_MASS_TO_KG)
                   + " AS mass_kg, event_horizon_radius_km AS radius_km, galactic_radius_pc"
                   if kind == "quasar" else "mass_kg, radius_km, galactic_radius_pc")
        if body["galactic_radius_pc"] is None:
            return keepout.NO_KEEP_OUT
        return keepout.compact_keep_out(body["mass_kg"], body["galactic_radius_pc"], body["radius_km"])
    if kind in ("nebula", "supernova_remnant", "asteroid_field"):
        row(_PHENOMENON_TYPE_TO_TABLE[kind], "id")
        return keepout.KeepOut(None, "none", keepout.PASS_THROUGH_NOTE)
    if kind in ("belt", "comet", "interstellar_comet"):
        row(object_ref.TABLES.get(kind) or _PHENOMENON_TYPE_TO_TABLE[kind], "id")
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
                pos["galaxy_pc"] = [c + mpc_to_pc(v) for c, v in zip(center, mpc)]
        return chain, pos

    def star_offset(star_id):
        """A star's system-local offset (km) from its system's origin."""
        if star_id is None:
            return [0.0, 0.0, 0.0]
        row = conn.execute(
            """
            SELECT s.role, ss.binary_primary_position_x_km AS px, ss.binary_primary_position_y_km AS py,
                   ss.binary_primary_position_z_km AS pz, ss.binary_secondary_position_x_km AS sx,
                   ss.binary_secondary_position_y_km AS sy, ss.binary_secondary_position_z_km AS sz
            FROM stars s JOIN star_systems ss ON ss.id = s.star_system_id WHERE s.id = ?
            """, (star_id,)).fetchone()
        if row is None:
            return [0.0, 0.0, 0.0]
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
    else:
        table = object_ref.TABLES[kind]
        columns = {
            "star": "name, star_system_id",
            "planet": "name, star_system_id, star_id, position_x_km, position_y_km, position_z_km",
            "moon": "name, star_system_id, planet_id, star_id, position_x_km, position_y_km, position_z_km",
            "belt": "star_system_id, orbital_index",
            "comet": "name, star_system_id, star_id, position_x_km, position_y_km, position_z_km",
        }[kind]
        row = one(table, object_id, columns)
        parents, positions = system_parents_and_positions(row["star_system_id"])
        name = f"Asteroid belt {row['orbital_index'] + 1}" if kind == "belt" else row["name"]
        if kind == "star":
            positions["system_km"] = star_offset(object_id)
        elif kind in ("planet", "comet"):
            positions["system_km"] = add(star_offset(row["star_id"]), offset(row))
        elif kind == "moon":
            planet = one("planets", row["planet_id"], "name, position_x_km, position_y_km, position_z_km")
            parents.append({"ref": object_ref.format("planet", row["planet_id"]), "kind": "planet",
                            "name": planet["name"]})
            positions["system_km"] = add(add(star_offset(row["star_id"]), offset(planet)), offset(row))
        if kind == "moon":
            siblings, truncated = _sibling_refs(conn, kind, table, "planet_id = ?", (row["planet_id"],), object_id)
        else:
            siblings, truncated = _sibling_refs(
                conn, kind, table, "star_system_id = ?", (row["star_system_id"],), object_id)
    return {
        "ref": object_ref.format(kind, object_id), "kind": kind, "id": object_id, "name": name,
        "parents": parents, "siblings": siblings, "siblings_truncated": truncated,
        "positions": positions,
        "keep_out": keep_out_radius_dict(conn, kind, object_id),
    }


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
    """A copy of a table row as a dict with its `uid` (a BINARY(10) column, GEN.170) as
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
            `sector_name` (`None` if it has no `sector_id`) and `nearest`
            (its stored nearest star systems, see `nearest_systems`). An
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
    detail["nearest"] = nearest_systems(conn, table, [phenomenon_id]).get(phenomenon_id, [])
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
    sector with no grid address (or when that is bigger). A small ring's
    cell reaches past the cube's half diagonal (a ring-0 pie wedge's outer
    corners sit a whole edge from its center), so sizing this from the
    cube alone could miss a cloud that reaches into the real cell.
    """
    reach = edge_pc * math.sqrt(3) / 2
    if sector["ring_index"] is None or sector["layer_index"] is None or sector["ring_slot_index"] is None:
        return reach
    center = (sector["center_x_pc"], sector["center_y_pc"], sector["center_z_pc"])
    for vertex in sector_cell_vertices_pc(
            sector["ring_index"], sector["layer_index"], sector["ring_slot_index"], edge_pc):
        reach = max(reach, math.dist(vertex, center))
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
            `radius_ly`, `distance_ly` (sector center to phenomenon
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
        "SELECT center_x_pc, center_y_pc, center_z_pc, edge_mpc, ring_index, layer_index, ring_slot_index"
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
            "class": phenomenon["class"],
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
    open_space_au = radius_km / constants.AU_TO_KM
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
            `schema.sql`'s "v14"/"v15" notes), `wikijs_url`/`mediawiki_url`
            (`star_systems.wikijs_url`/`mediawiki_url` -- `None` for
            whichever wiki (or both) this system hasn't been uploaded to
            yet; see `schema.sql`'s "v22" header note), `stars` (id/role/name/
            star_type/mass_kg/radius_km/temperature_k/luminosity_w/
            heliosphere_radius_km --
            `id` matches a `'wide'` binary's `planets`/`belts` rows' own
            `star_id`, disambiguating which star each orbits), `planets`
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
        "SELECT id, role, name, star_type, mass_kg, radius_km, temperature_k, luminosity_w, heliosphere_radius_km"
        " FROM stars WHERE star_system_id = ?"
        " ORDER BY CASE role WHEN 'primary' THEN 0 WHEN 'single' THEN 0 ELSE 1 END",
        (system_id,),
    ).fetchall()
    # mass_kg (already selected above) is what lets planetgen/web/maps/systemmap.py
    # split binary_mutual_position_x/y/z_km below into each star's own
    # mass-weighted offset from the system's barycenter, rather than
    # (incorrectly) anchoring the whole system on the primary alone.

    planet_rows = conn.execute(
        "SELECT * FROM planets WHERE star_system_id = ? ORDER BY orbital_index", (system_id,)
    ).fetchall()
    planet_stages = _life_stages(conn, "planet_evolutionary_paragraphs", "planet_id", "planets", system_id)
    moon_stages = _life_stages(conn, "moon_evolutionary_paragraphs", "moon_id", "moons", system_id)
    colonized = colonized_body_ids(conn, system_id)
    moons_by_planet = {}
    for moon in conn.execute(
        "SELECT * FROM moons WHERE star_system_id = ? ORDER BY planet_id, orbital_index", (system_id,)
    ).fetchall():
        moons_by_planet.setdefault(moon["planet_id"], []).append(moon)
    planets = []
    for planet in planet_rows:
        planet_dict = _with_life_fields(_with_printed_uid(planet), planet_stages, colonized["planets"])
        planet_dict["moons"] = [_with_life_fields(_with_printed_uid(m), moon_stages, colonized["moons"])
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
        "wikijs_url": system["wikijs_url"], "mediawiki_url": system["mediawiki_url"],
        "runaway_class": system["runaway_class"], "runaway_speed_kms": system["runaway_speed_kms"],
        "stars": [dict(s) for s in stars],
        "planets": planets,
        "belts": belts,
        "comets": [_with_printed_uid(c) for c in comets],
        "sector_siblings": sector_siblings,
        "nearest_neighbors": nearest_neighbors,
        "inside": cloud,
        "heliopause_au": heliopause_au,
        "heliopause_open_space_au": open_space_au,
    }


def _life_stages(conn, paragraph_table, id_column, body_table, system_id):
    """
    `{body id: life stage}` for every planet (or moon) in one system that
    has an evolutionary timeline -- see `evolution.life_stage_from_paragraphs`.
    One query for the whole system rather than one per body.
    """
    rows = conn.execute(
        f"SELECT p.{id_column} AS body_id, p.paragraph FROM {paragraph_table} p"
        f" JOIN {body_table} b ON b.id = p.{id_column}"
        f" WHERE b.star_system_id = ? ORDER BY p.{id_column}, p.position",
        (system_id,),
    ).fetchall()
    paragraphs = {}
    for row in rows:
        paragraphs.setdefault(row["body_id"], []).append(row["paragraph"])
    return {body_id: life_stage_from_paragraphs(texts) for body_id, texts in paragraphs.items()}


def _with_life_fields(body, stages, colonized=()):
    """
    Adds the system page's per-body life summary to one `planets`/`moons`
    row dict: `habitable` (its class is one of
    `tuning.HABITABLE_PLANET_CLASSES`, the same test
    `StarSystem.count_habitable` uses), `life_stage` (the most advanced
    evolutionary milestone its timeline reached, or `None`) and
    `inhabited` (habitable and reached a technological civilization --
    the page only ever describes a timeline for a habitable class, see
    `Planet._generate_life_and_flavor_paragraphs` -- or a colony stands on
    it: `colonized` is that table's ids from `colonized_body_ids`, schema
    v42).
    """
    body["habitable"] = body["planet_class"] in HABITABLE_PLANET_CLASSES
    body["life_stage"] = stages.get(body["id"]) if body["habitable"] else None
    body["inhabited"] = body["life_stage"] == "technological_civilization" or body["id"] in colonized
    return body


NEAREST_NEIGHBOR_COUNT = 3
"""int: How many nearest same-sector systems `system_detail` lists -- the
same count the stored `location` string was generated with
(`_db._location_for_entry`)."""


def _nearest_sector_siblings(system, sibling_rows, count=NEAREST_NEIGHBOR_COUNT):
    """
    The `count` systems in `system`'s own sector nearest to it, computed
    from their current rows rather than parsed back out of the stored
    `location` string. That string is written at generation time, before
    each neighbor's name has been made unique (`_db.insert_star_system`'s
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
        distance_ly = milliparsecs_to_ly(math.dist(origin, position))
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
            `galactic_radius_pc`, `ring_index`, `system_count`.
    """
    rows = conn.execute(
        """
        SELECT sec.id, sec.name, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.galactic_radius_pc, sec.ring_index, COALESCE(counts.n, 0) AS system_count
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
            "galactic_radius_pc": r["galactic_radius_pc"], "ring_index": r["ring_index"],
            "system_count": r["system_count"],
        }
        for r in rows
    ]


MADE_SECTOR_LIMIT = 2000
"""int: Most sectors `sectors_made` returns by address: a bigger run is shown
by its first sectors (the response still says the real total)."""


LAYER_SPEC_ROWS = 41
"""int: At most this many layers get a row in `galaxy_layer_specs`; a galaxy
with more shows an even sample that always includes the two ends and layer 0."""


def galaxy_layer_specs(conn, edge_pc):
    """
    The galaxy's layers as the Generate page shows them (ADM.28).

    Returns:
        dict | None: `None` before any layer is stored; else `count`,
            `lowest`, `highest` (layer indexes), `height_pc` (one layer's),
            `thickness_pc` (all layers), `radius_pc` (the widest layer's
            outer edge), `charted` (generated galaxy-placed sectors, all
            layers), `sampled` (whether `rows` is a sample) and `rows`,
            highest layer first: `layer`, `bottom_pc`, `top_pc`,
            `radius_pc` (outer edge) and `charted`.
    """
    from planetgen.galaxy.geometry import layer_bounds_pc, ring_bounds_pc
    layers = conn.execute("SELECT layer_index, outer_ring_index FROM galaxy_layer ORDER BY layer_index DESC").fetchall()
    if not layers:
        return None
    charted = {row["layer_index"]: row["n"] for row in conn.execute(
        "SELECT layer_index, COUNT(*) AS n FROM sectors WHERE center_x_pc IS NOT NULL GROUP BY layer_index").fetchall()}
    keep = layers
    if len(layers) > LAYER_SPEC_ROWS:
        picks = {round(i * (len(layers) - 1) / (LAYER_SPEC_ROWS - 1)) for i in range(LAYER_SPEC_ROWS)}
        picks.update(i for i, row in enumerate(layers) if row["layer_index"] == 0)
        keep = [layers[i] for i in sorted(picks)]
    rows = []
    for row in keep:
        bottom, top = layer_bounds_pc(row["layer_index"], edge_pc)
        rows.append({"layer": row["layer_index"], "bottom_pc": bottom, "top_pc": top,
                     "radius_pc": ring_bounds_pc(row["outer_ring_index"], edge_pc)[1],
                     "charted": charted.get(row["layer_index"], 0)})
    widest = max(row["outer_ring_index"] for row in layers)
    return {"count": len(layers), "lowest": layers[-1]["layer_index"], "highest": layers[0]["layer_index"],
            "height_pc": edge_pc, "thickness_pc": edge_pc * len(layers),
            "radius_pc": ring_bounds_pc(widest, edge_pc)[1], "charted": sum(charted.values()),
            "sampled": len(rows) < len(layers), "rows": rows}


def sectors_made(conn, since, until=None, limit=MADE_SECTOR_LIMIT):
    """
    The galaxy-placed sectors created from `since` to `until` (Unix
    seconds; `until=None` means now) -- what a generate job made, for the
    Galaxy Map's "Show on Galaxy Map" (ADM.31).

    Returns:
        dict: `total` (how many), `items` (up to `limit`, lowest id first:
            `id`, `name`, `x`/`y`/`z` pc, `ring_index`, `layer_index`,
            `ring_slot_index`).
    """
    where = "center_x_pc IS NOT NULL AND UNIX_TIMESTAMP(created_at) >= ?"
    params = [since]
    if until is not None:
        where += " AND UNIX_TIMESTAMP(created_at) <= ?"
        params.append(until)
    total = conn.execute(f"SELECT COUNT(*) AS n FROM sectors WHERE {where}", params).fetchone()["n"]
    rows = conn.execute(
        "SELECT id, name, center_x_pc, center_y_pc, center_z_pc, ring_index, layer_index, ring_slot_index "
        f"FROM sectors WHERE {where} ORDER BY id LIMIT ?", [*params, int(limit)],
    ).fetchall()
    return {
        "total": total,
        "items": [
            {"id": r["id"], "name": r["name"], "x": r["center_x_pc"], "y": r["center_y_pc"], "z": r["center_z_pc"],
             "ring_index": r["ring_index"], "layer_index": r["layer_index"], "ring_slot_index": r["ring_slot_index"]}
            for r in rows
        ],
    }


def galaxy_density_shape(conn):
    """
    The galaxy's stored density-skeleton shape (the singleton
    `galaxy_shape` row `planetgen plan` writes -- see
    `planetgen.galaxy.density.GalaxyShape` and `_db.get_galaxy_shape`),
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
            `galaxyDensity.model_terms` entry, plus `edge_pc`,
            `outer_ring_index`, `expected_system_count_at_density_1`.
            `None` if `planetgen plan` has never been run against this
            database (no `galaxy_shape` row yet).
    """
    skeleton = get_galaxy_shape(conn)
    if skeleton is None:
        return None
    shape = shape_with_terms(skeleton.shape)
    shape["edge_pc"] = skeleton.edge_pc
    shape["outer_ring_index"] = skeleton.outer_ring_index
    shape["expected_system_count_at_density_1"] = skeleton.expected_system_count_at_density_1
    return shape


# ---------------------------------------------------------------------
# Interactive 3D Galaxy Map viewport queries -- the placed-sector shape
# the cube tiles below reuse. Unlike `galaxy_placed_sectors`/`galaxy_density_shape` above (each
# called once per page load for the flat, whole-galaxy overview map),
# these are scoped to a moving viewport -- see `planetgen.galaxy.
# viewport`'s own module docstring for the three content tiers
# (placed/planned/density) this combines.
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
            `ring_index`, `layer_index`, `ring_slot_index`, `designation`
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
               sec.galactic_radius_pc, sec.ring_index, sec.layer_index, sec.ring_slot_index,
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
# Cube tiles -- backs GET /api/galaxy/tiles and GET /api/galaxy/stamp. See
# `planetgen.galaxy.viewport`'s "Cube tiles" section for the model.
# ---------------------------------------------------------------------

GALAXY_TILE_MAX_PLACED = 250
"""int: Most placed sectors one tile returns. Past this, a tile returns an
evenly spaced sample (every Nth sector by id, so the same sample every
time) -- a zoomed-out view can't show thousands of sub-pixel sectors
anyway, and zooming in switches to smaller tiles that hold the rest.

The sample has to be spread out, not just the first ids: a local
neighborhood is generated in order outward from the core, so its lowest
ids are its core-facing side. Taking the lowest 250 used to draw a
generated 100 ly sphere as a core-facing bowl, cut off flat where the cap
ran out. A fully generated tile can hold thousands of sectors, so this
isn't only a far-zoomed-out case."""

GALAXY_TILE_FILLED_SCALE = 2048
"""int: A tile's filled-sector summary (`galaxy_filled_in_box`) groups
sectors into cells at most `tile edge / GALAXY_TILE_FILLED_SCALE` across.
The map fetches tiles one to two view radii across (1.6 orbit radii) and
never draws blocks narrower than 4 pixels, so on any screen up to about
2,000 pixels tall a cell is never bigger than the blocks the map draws
with that tile, and nests inside them."""

GALAXY_TILE_MAX_FILLED_CELLS = 5000
"""int: Most cells one tile's filled summary lists. Past this the cells
grow three times bigger until they fit."""

GALAXY_TILE_MAX_CLOUDS = 200
"""int: Most nebulae and supernova remnants one tile lists
(`galaxy_clouds_in_box`), the largest first. A tile wide enough to hold
more shows the map at a scale where the smaller ones are under a pixel
anyway."""

MAX_TILES_PER_REQUEST = 128
"""int: Most tiles one `/api/galaxy/tiles` request may ask for. The map
needs at most 27 view tiles plus about 64 planned tiles at once."""


def galaxy_sectors_in_box(conn, lo, hi, limit=GALAXY_TILE_MAX_PLACED):
    """
    Placed sectors whose center lies in the half-open box `[lo, hi)`,
    ordered by id, at most `limit` -- the "placed" half of one tile. When
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
               galactic_radius_pc, ring_index, layer_index, ring_slot_index, edge_mpc
        FROM (
            SELECT sec.id, sec.name, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
                   sec.galactic_radius_pc, sec.ring_index, sec.layer_index, sec.ring_slot_index,
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
        "ring_index": r["ring_index"], "layer_index": r["layer_index"], "ring_slot_index": r["ring_slot_index"],
        "designation": provisional_sector_designation(*address) if address is not None else None,
        "system_count": system_count,
        "edge_ly": pc_to_ly(edge_pc) if edge_pc else None,
    }


def galaxy_filled_in_box(conn, lo, hi, tile_edge_pc, edge_pc, max_cells=GALAXY_TILE_MAX_FILLED_CELLS):
    """
    Every placed sector whose center lies in the half-open box `[lo, hi)`,
    counted into cells the Galaxy Map sums into its blocks
    (`static/galaxymap3d.js`): unlike a tile's `placed` list, nothing is
    left out, so a block's filled count is right at every zoom.

    Cells are `g` sectors a side, `g` a power of 3 (see
    `GALAXY_TILE_FILLED_SCALE`), laid out like the map's blocks: cell ring
    `I` is sector rings `I*g .. I*g + g - 1`, cell layer `S` is sector
    layers `S*g - (g-1)/2 .. S*g + (g-1)/2`, and a cell ring has
    `max(3, round(2 pi (I + 1/2)))` equal wedges counterclockwise from +X,
    each holding the sectors whose center angle falls in it. At `g = 1` a
    cell is one sector, listed with its id, name and system count.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        lo (tuple): `(x, y, z)` inclusive lower corner, parsecs.
        hi (tuple): `(x, y, z)` exclusive upper corner, parsecs.
        tile_edge_pc (float): The tile's edge (sets the cell size).
        edge_pc (float): The sector edge.
        max_cells (int): See `GALAXY_TILE_MAX_FILLED_CELLS`.

    Returns:
        dict: `g` and `cells`: at `g = 1`, `[ring, layer, slot, id,
            system_count, name]` per sector; otherwise `[ring, layer,
            wedge, count]` per cell with any placed sector, in cell units.
    """
    g = 1
    while g * 3 * edge_pc <= tile_edge_pc / GALAXY_TILE_FILLED_SCALE:
        g *= 3
    box = (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2])
    where = (
        "center_x_pc >= ? AND center_x_pc < ? AND center_y_pc >= ? AND center_y_pc < ? "
        "AND center_z_pc >= ? AND center_z_pc < ? AND ring_index IS NOT NULL"
    )
    if g == 1:
        rows = conn.execute(
            f"SELECT id, name, ring_index, layer_index, ring_slot_index FROM sectors WHERE {where} "
            f"ORDER BY id LIMIT ?",
            box + (max_cells + 1,),
        ).fetchall()
        if len(rows) <= max_cells:
            counts = {}
            if rows:
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
            return {"g": 1, "cells": [
                [r["ring_index"], r["layer_index"], r["ring_slot_index"], r["id"], counts.get(r["id"], 0), r["name"]]
                for r in rows
            ]}
        g = 3
    while True:
        half = (g - 1) // 2
        rows = conn.execute(
            f"""
            SELECT cell_ring, cell_layer,
                   FLOOR(MOD(ATAN2(center_y_pc, center_x_pc) + 2 * PI(), 2 * PI()) / (2 * PI())
                         * GREATEST(3, ROUND(2 * PI() * (cell_ring + 0.5)))) AS cell_wedge,
                   COUNT(*) AS n
            FROM (
                SELECT FLOOR(ring_index / ?) AS cell_ring, FLOOR((layer_index + ?) / ?) AS cell_layer,
                       center_x_pc, center_y_pc
                FROM sectors WHERE {where}
            ) placed
            GROUP BY cell_ring, cell_layer, cell_wedge
            LIMIT ?
            """,
            (g, half, g) + box + (max_cells + 1,),
        ).fetchall()
        if len(rows) <= max_cells:
            wedges = lambda ring: max(3, round(2 * math.pi * (ring + 0.5)))  # noqa: E731
            return {"g": g, "cells": [
                [int(r["cell_ring"]), int(r["cell_layer"]),
                 min(int(r["cell_wedge"]), wedges(int(r["cell_ring"])) - 1), int(r["n"])]
                for r in rows
            ]}
        g *= 3


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
    read once per tile request."""
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
    shows from every tile it covers (the map drops the repeats by type and
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


GALAXY_TILE_MAX_BRIGHT_STARS = 400
"""int: Most pre-placed bright stars (`bright_stars`, v43) one tile lists,
the most luminous first. The map draws each as a point of light the same
few pixels across at every zoom, so a zoomed-out view of 27 tiles shows
about ten thousand, tracing the spiral arms; zooming in switches to
smaller tiles that hold the dimmer ones."""

GALAXY_TILE_BRIGHTEST_SAMPLE = 100000
"""int: How many of the galaxy's most luminous bright stars `galaxy_tiles`
reads in one go (once per request, only when a tile needs it) to answer
tiles too big to query on their own -- see
`GALAXY_TILE_BRIGHT_STAR_ROW_BUDGET`. A tile of a few kiloparsecs holds a
few hundred of them, so zoomed out every tile still fills its
`GALAXY_TILE_MAX_BRIGHT_STARS`."""

GALAXY_TILE_BRIGHT_STAR_ROW_BUDGET = 50000
"""int: Most `bright_stars` rows one tile's own query may expect to read
(the cheaper of `galaxy_bright_stars_in_box`'s two ways). A tile that
would read more -- a big tile with millions of stars, or one mostly past
the galaxy's edge, where the luminosity walk finds few stars inside --
takes its stars from the galaxy-wide sample instead: the brightest
stars inside it, thinned by luminosity when the sample holds fewer than
the tile's cap there. Before this a zoomed-out request could run past
the web app's 30 s API timeout and the stars vanished (MAP.47)."""

BRIGHT_STAR_MAX_EXACT_RANGES = 2000
"""int: Most `(ring, layer)` index ranges one tile's bright-star query
lists before it reads each ring's layers in one range instead."""


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


def _ring_slot_ranges(ring, intervals):
    """`intervals`' slots of `ring`, as inclusive `(first, last)` pairs
    (rounded outward, so every slot a longitude touches is in)."""
    n = ring_sector_count(ring)
    ranges = []
    for start, end in intervals:
        first = max(0, int(math.floor(start * n / (2 * math.pi))))
        last = min(n - 1, int(math.floor(end * n / (2 * math.pi))))
        ranges.append((first, last))
    return ranges


def _bright_star_bands(conn, lo, hi, edge_pc):
    """
    `[(ring, layer_min, layer_max), ...]`: the rings the box `[lo, hi)`
    reaches, each with the layers it reaches, cut to the stored outline
    (`galaxy_column`) when there is one -- no bright star sits outside it.
    """
    corners_r = [math.hypot(x, y) for x in (lo[0], hi[0]) for y in (lo[1], hi[1])]
    near_x = min(max(0.0, lo[0]), hi[0])
    near_y = min(max(0.0, lo[1]), hi[1])
    ring_lo = ring_index_at(math.hypot(near_x, near_y), edge_pc)
    ring_hi = ring_index_at(max(corners_r), edge_pc)
    layer_lo = layer_index_at(lo[2], edge_pc)
    layer_hi = layer_index_at(hi[2], edge_pc)
    columns = {
        row["ring_index"]: (row["layer_index_min"], row["layer_index_max"])
        for row in conn.execute(
            "SELECT ring_index, layer_index_min, layer_index_max FROM galaxy_column WHERE ring_index BETWEEN ? AND ?",
            (ring_lo, ring_hi),
        ).fetchall()
    }
    has_outline = bool(columns) or conn.execute("SELECT 1 FROM galaxy_column LIMIT 1").fetchone() is not None
    bands = []
    for ring in range(ring_lo, ring_hi + 1):
        bottom, top = layer_lo, layer_hi
        if has_outline:
            if ring not in columns:
                continue
            bottom, top = max(bottom, columns[ring][0]), min(top, columns[ring][1])
        if bottom <= top:
            bands.append((ring, bottom, top))
    return bands


def bright_star_scatter_status(conn):
    """
    Whether the galaxy's bright-star scatter (`planetgen plan`) has run,
    and with what threshold and seed -- so the Generate page can say so
    and offer the right next step.

    Returns:
        dict: `scattered` (bool), `min_luminosity_sol` and `seed` (both
            `None` until a scatter runs), and `default_min_luminosity_sol`
            (`BRIGHT_STAR_MIN_LUMINOSITY_SOL`, what a plain plan uses).
    """
    row = conn.execute(
        "SELECT bright_star_min_luminosity_sol, bright_star_seed FROM galaxy_shape WHERE id = 1").fetchone()
    scattered = row is not None and row["bright_star_min_luminosity_sol"] is not None
    return {
        "scattered": scattered,
        "min_luminosity_sol": float(row["bright_star_min_luminosity_sol"]) if scattered else None,
        "seed": row["bright_star_seed"] if scattered else None,
        "default_min_luminosity_sol": tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL,
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
        "yerkes_class": row["yerkes_class"], "ring_index": row["ring_index"],
        "layer_index": row["layer_index"], "ring_slot_index": row["ring_slot_index"],
        "system_id": row["star_system_id"],
    }


def bright_stars_in_sector(conn, ring_index, layer_index, ring_slot_index, unfilled_only=True):
    """
    The pre-placed bright stars (`bright_stars`) in one sector cell, most
    luminous first -- the stars a sector page lists for a cell that hasn't
    been filled yet (filling builds each into a system, see
    `planetgen fill`). Reads `idx_bright_stars_address`, so it's cheap
    at any galaxy size.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        ring_index, layer_index, ring_slot_index (int): The cell's address.
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
               star_type, population, yerkes_class, ring_index, layer_index, ring_slot_index, star_system_id
        FROM bright_stars
        WHERE ring_index = ? AND layer_index = ? AND ring_slot_index = ?{unfilled}
        ORDER BY luminosity_w DESC, id
        """,
        (int(ring_index), int(layer_index), int(ring_slot_index)),
    ).fetchall()
    return [_bright_star_entry(row) for row in rows]


UNCHARTED_SCATTER_LIMIT = 500
"""int: Most scattered phenomena one uncharted sector lists (MAP.162)."""


def uncharted_sector_contents(conn, ring_index, layer_index, ring_slot_index):
    """
    What the scatters left in one sector cell that nothing has built yet
    (MAP.162): its waiting bright stars (`bright_stars_in_sector`) and its
    unbuilt scattered phenomena (`phenomenon_scatter`), so a cell with no
    generated contents can still be opened and its objects looked at.

    Returns:
        dict: `stars` (as `bright_stars_in_sector`) and `scattered`, one
            `{"id", "kind", "subtype", "type", "x", "y", "z"}` per unbuilt
            scatter row (`type` as `_SCATTER_TYPES`; positions in
            galaxy-frame parsecs), at most `UNCHARTED_SCATTER_LIMIT`.
    """
    address = (int(ring_index), int(layer_index), int(ring_slot_index))
    try:
        stars = bright_stars_in_sector(conn, *address)
    except pymysql.err.ProgrammingError:  # no bright_stars table yet (before v43)
        stars = []
    rows = conn.execute(
        """
        SELECT id, kind, subtype, position_x_mpc, position_y_mpc, position_z_mpc
        FROM phenomenon_scatter
        WHERE ring_index = ? AND layer_index = ? AND ring_slot_index = ? AND built_at IS NULL
        ORDER BY id
        LIMIT ?
        """,
        (*address, UNCHARTED_SCATTER_LIMIT),
    ).fetchall()
    scattered = [
        {"id": row["id"], "kind": row["kind"], "subtype": row["subtype"],
         "type": _SCATTER_TYPES.get(row["kind"], row["kind"]),
         "x": row["position_x_mpc"] / MPC_PER_PC, "y": row["position_y_mpc"] / MPC_PER_PC,
         "z": row["position_z_mpc"] / MPC_PER_PC}
        for row in rows
    ]
    return {"stars": stars, "scattered": scattered}


def _stratified(stars, limit):
    """
    At most `limit` of `stars`, an equal share from each population: the
    most luminous of each, a population with fewer than its share leaving
    the rest to the others. Ranking all of them by luminosity alone puts
    the young blue stars first and leaves no room for the old giants of
    the bulge and the thick disk, which top out near 2,500 Lsun while a
    young star reaches a million.

    Returns:
        list[dict]: Most luminous first (ties by descending id).
    """
    groups = {}
    for star in sorted(stars, key=lambda star: (-star["luminosity_sol"], -star["id"])):
        groups.setdefault(star["population"], []).append(star)
    picks, remaining = [], limit
    # Smallest group first, so what it leaves unused is shared by the rest.
    pending = sorted(groups.values(), key=len)
    for index, group in enumerate(pending):
        share = -(-remaining // (len(pending) - index))
        picks += group[:share]
        remaining -= min(len(group), share)
    picks.sort(key=lambda star: (-star["luminosity_sol"], -star["id"]))
    return picks


def galaxy_brightest_stars(conn, count=GALAXY_TILE_BRIGHTEST_SAMPLE):
    """
    The galaxy's most luminous bright stars of each population, most
    luminous first (ties by descending id) -- `galaxy_tiles`' sample for
    tiles too big to query on their own. One walk down
    `idx_bright_stars_population` per population and height: the `count`
    stars on the plane, split equally among the four populations, and half
    as many 250 pc or more off it (GEN.117). Ranked together by
    luminosity, the sample held only young stars: the old giants of the
    bulge and the thick disk never reach the luminosity of the young
    ones, so the bulge and everything off the plane's thin young layer
    were missing from every zoomed-out tile.

    Returns:
        list[dict]: As `galaxy_bright_stars_in_box`.
    """
    def top(population, off_plane, limit):
        return [
            _bright_star_entry(row) for row in conn.execute(
                """
                SELECT id, position_x_mpc, position_y_mpc, position_z_mpc, luminosity_w, temperature_k,
                       radius_km, star_type, population, yerkes_class, ring_index, layer_index,
                       ring_slot_index, star_system_id
                FROM bright_stars FORCE INDEX (idx_bright_stars_population)
                WHERE population = ? AND off_plane = ?
                ORDER BY luminosity_w DESC, id DESC
                LIMIT ?
                """,
                (population, off_plane, limit),
            ).fetchall()
        ]

    picks = []
    for off_plane, total in ((0, int(count)), (1, int(count) // 2)):
        limits = {population: -(-total // len(BRIGHT_STAR_POPULATIONS)) for population in BRIGHT_STAR_POPULATIONS}
        found = {population: top(population, off_plane, limit) for population, limit in limits.items()}
        # A population short of its share leaves the rest to those that
        # filled theirs, once.
        spare = total - sum(len(stars) for stars in found.values())
        full = [population for population, stars in found.items() if len(stars) >= limits[population]]
        if spare >= len(full) > 0:
            for population in full:
                found[population] = top(population, off_plane, limits[population] + spare // len(full))
        for stars in found.values():
            picks += stars
    picks.sort(key=lambda star: (-star["luminosity_sol"], -star["id"]))
    return picks


BRIGHT_STAR_PLANE_HALF_THICKNESS_PC = 250.0
"""float: A tile's bright stars are picked in three height bands (GEN.117):
the plane, `|z|` below this, and the two sides above and below it. The
Galaxy Map's old picks, the most luminous first, kept only the young blue
stars on the plane and none of the old giants above it."""

BRIGHT_STAR_OFF_PLANE_SHARE = 0.5
"""float: The share of a tile's bright-star budget (`limit`) split between
the bands off the plane, when the tile reaches them: each gets half of
it, and the plane the rest. A band holding fewer stars than its share
leaves the rest to the plane."""


def _height_bands(lo, hi):
    """The box `[lo, hi)` cut at the plane band's edges: `[(lo, hi, on_plane)]`,
    only the bands the box reaches."""
    edge = BRIGHT_STAR_PLANE_HALF_THICKNESS_PC
    cuts = [-edge, edge]
    bounds = [lo[2]] + [z for z in cuts if lo[2] < z < hi[2]] + [hi[2]]
    bands = []
    for bottom, top in zip(bounds, bounds[1:]):
        bands.append(((lo[0], lo[1], bottom), (hi[0], hi[1], top), -edge <= bottom and top <= edge))
    return bands


def galaxy_bright_stars_in_box(conn, lo, hi, edge_pc, limit=GALAXY_TILE_MAX_BRIGHT_STARS, unfilled_only=False,
                               brightest=None):
    """
    The most luminous pre-placed bright stars (`bright_stars`) in the box
    `[lo, hi)`, at most `limit` -- the stars the Galaxy Map draws before
    (and after) their sectors are filled.

    A box reaching beyond the galactic plane's band picks within each
    height band (`BRIGHT_STAR_PLANE_HALF_THICKNESS_PC`) so the old giants
    above and below the plane are not all crowded out by the young stars
    on it (GEN.117): each band off the plane gets its share of `limit`
    (`BRIGHT_STAR_OFF_PLANE_SHARE`), and what a band leaves unused goes
    to the plane. See `_galaxy_bright_stars_in_band` for the rest.
    """
    bands = _height_bands(lo, hi)
    if len(bands) == 1:
        return _galaxy_bright_stars_in_band(conn, lo, hi, edge_pc, limit, unfilled_only, brightest,
                                            off_plane=0 if bands[0][2] else 1)
    off = [band for band in bands if not band[2]]
    share = int(limit * BRIGHT_STAR_OFF_PLANE_SHARE) // max(1, len(off))
    picks = []
    for band_lo, band_hi, on_plane in bands:
        if not on_plane:
            picks += _galaxy_bright_stars_in_band(conn, band_lo, band_hi, edge_pc, share, unfilled_only, brightest,
                                                  off_plane=1)
    plane = [band for band in bands if band[2]]
    room = limit - len(picks)
    for band_lo, band_hi, _on in plane:
        picks += _galaxy_bright_stars_in_band(conn, band_lo, band_hi, edge_pc, room, unfilled_only, brightest,
                                              off_plane=0)
    picks.sort(key=lambda star: (-star["luminosity_sol"], -star["id"]))
    return picks[:limit]


def _galaxy_bright_stars_in_band(conn, lo, hi, edge_pc, limit, unfilled_only, brightest, off_plane):
    """
    The most luminous pre-placed bright stars in the box `[lo, hi)`, at
    most `limit`, with no height banding (`galaxy_bright_stars_in_box`):
    `off_plane` says which band the box is in (1 for 250 pc or more off
    the plane). The picks are an equal share from each population, the
    most luminous of each (`_stratified`).

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
        limit (int): See `GALAXY_TILE_MAX_BRIGHT_STARS`.
        unfilled_only (bool): Only the stars whose sector hasn't been
            filled yet (`system_id` is `None`) -- what a page listing the
            stars still waiting in a box shows.
        brightest (callable or None): Returns the galaxy-wide sample
            (`galaxy_brightest_stars`). When given and the box would read
            more than `GALAXY_TILE_BRIGHT_STAR_ROW_BUDGET` rows either way,
            the answer is the sample's stars inside the box instead: the
            same stars whenever the sample holds `limit` of them there,
            else only those (thinned by luminosity).

    Returns:
        list[dict]: Most luminous first (ties by descending id): `id`, `x`/`y`/`z` (parsecs),
            `luminosity_sol`, `temperature_k`, `radius_sol`, `star_type`,
            `yerkes_class`, `ring_index`, `layer_index`,
            `ring_slot_index` and `system_id` (the star's system once its
            sector is filled, else `None`).
    """
    total = int(conn.execute("SELECT COALESCE(MAX(id), 0) AS n FROM bright_stars").fetchone()["n"])
    if not total:
        return []
    bands = _bright_star_bands(conn, lo, hi, edge_pc)
    if not bands:
        return []
    box = (
        "position_x_mpc >= ? AND position_x_mpc < ? AND position_y_mpc >= ? AND position_y_mpc < ? "
        "AND position_z_mpc >= ? AND position_z_mpc < ?"
    )
    if unfilled_only:
        box += " AND star_system_id IS NULL"
    box_params = [int(math.ceil(v * MPC_PER_PC)) for pair in zip(lo, hi) for v in pair]

    # Rows read each way, assuming stars spread evenly over the disk: the
    # address path reads its rings' full layer bands (the slot test only
    # filters), the luminosity path about `limit / share of the galaxy`.
    outer_ring = int(conn.execute("SELECT COALESCE(MAX(ring_index), 0) AS r FROM galaxy_column").fetchone()["r"])
    outer_ring = outer_ring or bands[-1][0]
    disk = math.pi * ((outer_ring + 1) * edge_pc) ** 2
    intervals = _box_angle_intervals(lo, hi)
    ring_share = sum(2 * ring + 1 for ring, _b, _t in bands) / float((outer_ring + 1) ** 2)
    angle_share = min(1.0, sum(end - start for start, end in intervals) / (2 * math.pi))
    box_share = min(1.0, (hi[0] - lo[0]) * (hi[1] - lo[1]) / disk)
    # The address index is (ring, layer, slot), so only the box's own
    # longitudes of each ring are read.
    by_address = ring_share * angle_share * total
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
        return _stratified(found, limit)

    columns = ("id, position_x_mpc, position_y_mpc, position_z_mpc, luminosity_w, temperature_k, radius_km, "
               "star_type, population, yerkes_class, ring_index, layer_index, ring_slot_index, star_system_id")
    # Both keys descending, so the luminosity path walks its index
    # backwards and stops at `limit`: with `id` ascending MariaDB and MySQL
    # sort every star in the box first, millions in a zoomed-out tile.
    order = "luminosity_w DESC, id DESC"
    if by_luminosity < by_address:
        # One walk per population down (population, off_plane, luminosity),
        # each stopping at `limit`: a single walk down the luminosity index
        # would meet only young stars (`_stratified`).
        selects, params = [], []
        for population in BRIGHT_STAR_POPULATIONS:
            selects.append(
                f"(SELECT {columns} FROM bright_stars FORCE INDEX (idx_bright_stars_population) "
                f"WHERE population = ? AND off_plane = ? AND {box} ORDER BY {order} LIMIT ?)")
            params += [population, off_plane] + box_params + [int(limit)]
        rows = conn.execute(" UNION ALL ".join(selects), params).fetchall()
        return _stratified([_bright_star_entry(row) for row in rows], limit)
    clauses, params = [], []
    exact = sum(t - b + 1 for _r, b, t in bands) <= BRIGHT_STAR_MAX_EXACT_RANGES
    for ring, bottom, top in bands:
        slots = _ring_slot_ranges(ring, intervals)
        slot_test = " OR ".join("ring_slot_index BETWEEN ? AND ?" for _ in slots)
        slot_params = [v for pair in slots for v in pair]
        layers = [(layer, layer) for layer in range(bottom, top + 1)] if exact else [(bottom, top)]
        for first, last in layers:
            clauses.append(f"(ring_index = ? AND layer_index BETWEEN ? AND ? AND ({slot_test}))")
            params += [ring, first, last] + slot_params
    # The address path reads the box's stars once and keeps the `limit`
    # most luminous of each population.
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
        params + box_params + [int(limit)],
    ).fetchall()
    return _stratified([_bright_star_entry(row) for row in rows], limit)


GALAXY_TILE_MAX_GENERATED_STARS = 1000
"""int: Most stars of generated systems one tile lists (MAP.51), the most
luminous first. A finest tile (16 pc, 64 sectors) of a filled region
holds a few hundred, so up close every star shows, red dwarfs included."""

GENERATED_STAR_FLOOR_SOL_AT_32_PC = 0.004
"""float: The faintest generated star (solar luminosities) a 32 pc tile
lists. Each coarser tile level's floor is four times the last (the floor
grows with the square of the tile's edge, as a star's apparent
brightness falls with the square of its distance), so zooming in shows
fainter and fainter stars: about 1 L_sun at 512 pc tiles, 260 L_sun at
8 kpc. The finest tiles (`TILE_MAX_LEVEL`) list every star."""

GENERATED_STAR_MAX_FLOOR_SOL = 1000.0
"""float: Tiles whose floor is past this list no generated stars at all:
they span a good share of the galaxy, where the pre-placed bright stars
(`galaxy_bright_stars_in_box`) already show everything that bright."""

GALAXY_TILE_MAX_DETAIL_STARS = 4000
"""int: Most stars of generated systems a finest tile (`TILE_MAX_LEVEL`,
16 pc) lists (MAP.80): the map fetches these around a sector it is
zoomed to, so the sector shows every star it holds. Only a tile deep in
the bulge holds more, and then its faintest are left out."""

GALAXY_TILE_GENERATED_STAR_SECTOR_BUDGET = 1500
"""int: Most generated sectors one tile reads stars from. A coarse tile
over a large filled region reads an even sample of its sectors (every
k-th, by a hash of the id) instead, so the request stays bounded; zooming in reaches
tiles small enough to read all of them."""


def generated_star_floor_sol(level):
    """
    The faintest generated star a tile of `level` lists, in solar
    luminosities (0 = every star), or `None` when it lists none -- see
    `GENERATED_STAR_FLOOR_SOL_AT_32_PC`.
    """
    if level >= TILE_MAX_LEVEL:
        return 0.0
    floor = GENERATED_STAR_FLOOR_SOL_AT_32_PC * (TILE_ROOT_EDGE_PC / 2 ** level / 32.0) ** 2
    return None if floor > GENERATED_STAR_MAX_FLOOR_SOL else floor


GALAXY_TILE_STAR_BUDGET = {3: 150, 4: 200, 5: 250, 6: 300, 7: 350, 8: 400, 9: 500, 10: 650, 11: 850}
"""dict: MAP.116's one table of how many generated stars a Galaxy Map tile
lists, by tile level (the levels in `generated_star_floor_sol`'s reach and
coarser than `TILE_MAX_LEVEL`). The map shows a roughly constant number of
tiles at any zoom, so a coarser tile gets fewer stars: a dense filled
region stays readable two or three zoom levels out from a sector, where
1000 per tile used to blur into a smear. The brightest stars win; the
brightness floor (`generated_star_floor_sol`) is the other half of the rule.
A finest tile lists up to `GALAXY_TILE_MAX_DETAIL_STARS` instead."""

SECTOR_ALLOWANCE_FUDGE = 3
"""int: A tile's budget is shared out by sector: each generated sector may
put `SECTOR_ALLOWANCE_FUDGE * budget / sectors` stars (at least one) on
the tile, its brightest, so one dense sector can't take the whole budget
while a sparse sector, which has fewer stars than that, keeps all it has
and the leftover room goes to the brighter stars of the dense ones."""

GALAXY_TILE_POINT_BUDGET = {10: 40, 11: 100}
"""dict: Most point phenomena (black holes, neutron stars, quasars) a tile
of that level lists (MAP.116); a finest tile lists up to
`GALAXY_TILE_MAX_POINTS`. The most luminous come first, so a remnant shows
where the budget reaches it."""


GALAXY_VIEW_MAX_TILES = 27
"""int: Most tiles a view's own level needs: the view sphere is at most as
wide as a tile (`viewport.tile_level_for_view_radius`), so it touches at
most 3 tiles along each axis."""

GALAXY_VIEW_MAX_DETAIL_TILES = 8
"""int: Most finest tiles the map adds around the target once zoomed to a
sector (`DETAIL_RADIUS_PC` is 8 pc, half a finest tile's edge, so the
sphere touches at most 2 along each axis)."""

GALAXY_VIEW_MAX_STARS = 70000
"""int: MAP.109's stated cap on the stars (pre-placed and generated) one
view's tiles can carry, about 8 MB of JSON before compression. It holds by
construction -- every tile's lists are capped (`GALAXY_TILE_MAX_BRIGHT_STARS`,
`GALAXY_TILE_STAR_BUDGET`, `GALAXY_TILE_MAX_DETAIL_STARS`) -- and
`galaxy_view_star_cap` adds those caps up, so a test fails when a budget is
raised past it. A real view carries far less: the caps are for a tile packed
with stars, and most of the galaxy's tiles are not."""


def galaxy_view_star_cap():
    """
    The most stars one view can fetch: `GALAXY_VIEW_MAX_TILES` tiles of its
    own level (each with the bright-star cap and the biggest level budget)
    plus `GALAXY_VIEW_MAX_DETAIL_TILES` finest tiles.
    """
    coarse = GALAXY_VIEW_MAX_TILES * (GALAXY_TILE_MAX_BRIGHT_STARS + max(GALAXY_TILE_STAR_BUDGET.values()))
    detail = GALAXY_VIEW_MAX_DETAIL_TILES * (GALAXY_TILE_MAX_BRIGHT_STARS + GALAXY_TILE_MAX_DETAIL_STARS)
    return coarse + detail


def generated_star_budget(level, sector_count):
    """
    `(tile budget, per-sector allowance)` for the generated stars of a tile of
    `level` holding `sector_count` generated sectors -- `GALAXY_TILE_STAR_BUDGET`
    shared by `SECTOR_ALLOWANCE_FUDGE`. A finest tile lists up to
    `GALAXY_TILE_MAX_DETAIL_STARS` with no per-sector allowance.
    """
    if level >= TILE_MAX_LEVEL:
        return GALAXY_TILE_MAX_DETAIL_STARS, None
    budget = GALAXY_TILE_STAR_BUDGET.get(level, GALAXY_TILE_MAX_GENERATED_STARS)
    allowance = max(1, int(math.ceil(SECTOR_ALLOWANCE_FUDGE * budget / float(max(1, sector_count)))))
    return budget, allowance


def galaxy_generated_stars_in_box(conn, lo, hi, min_luminosity_sol, sector_count,
                                  limit=GALAXY_TILE_MAX_GENERATED_STARS, per_sector=None):
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
            box holds, as its `galaxy_filled_in_box` summary counts them.
        limit (int): See `GALAXY_TILE_MAX_GENERATED_STARS`.
        per_sector (int or None): Most stars one sector contributes, its
            brightest (`generated_star_budget`); `None` is no limit.

    Returns:
        list[dict]: Most luminous first (ties by descending id): `id`
            (the `stars` row), `name`, `x`/`y`/`z` (galaxy-frame
            parsecs), `luminosity_sol`, `temperature_k`, `radius_sol`,
            `star_type`, `ring_index`, `layer_index`, `ring_slot_index`
            and `system_id`.
    """
    if sector_count <= 0:
        return []
    stride = max(1, int(math.ceil(sector_count / float(GALAXY_TILE_GENERATED_STAR_SECTOR_BUDGET))))
    select = """
        SELECT st.id, st.name, st.luminosity_w, st.temperature_k, st.radius_km, st.star_type,
               ss.id AS system_id, ss.position_x_mpc, ss.position_y_mpc, ss.position_z_mpc,
               sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.ring_index, sec.layer_index, sec.ring_slot_index"""
    rank = (", ROW_NUMBER() OVER (PARTITION BY sec.id ORDER BY st.luminosity_w DESC, st.id DESC) AS sector_rank"
            if per_sector is not None else "")
    source = """
        FROM sectors sec FORCE INDEX (idx_sectors_center)
        JOIN star_systems ss ON ss.sector_id = sec.id
        JOIN stars st ON st.star_system_id = ss.id
        LEFT JOIN bright_stars b ON b.star_system_id = ss.id AND st.role <> 'secondary'
        WHERE sec.center_x_pc >= ? AND sec.center_x_pc < ?
          AND sec.center_y_pc >= ? AND sec.center_y_pc < ?
          AND sec.center_z_pc >= ? AND sec.center_z_pc < ?
          AND sec.ring_index IS NOT NULL AND MOD(CRC32(sec.id), ?) = 0
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
    return [{
        "id": row["id"], "name": row["name"],
        "x": round(row["center_x_pc"] + row["position_x_mpc"] / MPC_PER_PC, 3),
        "y": round(row["center_y_pc"] + row["position_y_mpc"] / MPC_PER_PC, 3),
        "z": round(row["center_z_pc"] + row["position_z_mpc"] / MPC_PER_PC, 3),
        "luminosity_sol": float("%.4g" % (row["luminosity_w"] / constants.SOLAR_LUMINOSITY)),
        "temperature_k": round(row["temperature_k"]),
        "radius_sol": float("%.3g" % _radius_sol(row["radius_km"])),
        "star_type": row["star_type"], "ring_index": row["ring_index"], "layer_index": row["layer_index"],
        "ring_slot_index": row["ring_slot_index"], "system_id": row["system_id"],
    } for row in rows]


POINT_PHENOMENON_MIN_LEVEL = 10
"""int: The coarsest tile level (64 pc tiles) that lists the point-like
phenomena -- black holes, neutron stars (pulsars) and quasars (MAP.80).
That is the level the map fetches once it is zoomed to a sector, so a
sector shows every one of them; farther out they would be lost among the
stars anyway."""

GALAXY_TILE_MAX_POINTS = 200
"""int: Most point-like phenomena one tile lists, the most luminous
first. A 64 pc tile holds a handful."""

_POINT_PHENOMENON_TABLES = (
    ("black_holes", "black_hole", "(CASE WHEN has_accretion_disk THEN 'accreting' ELSE 'quiescent' END)"),
    ("neutron_stars", "neutron_star", "pulsar_type"),
    ("quasars", "quasar", "(CASE WHEN is_radio_loud THEN 'radio-loud' ELSE 'radio-quiet' END)"),
)
"""tuple: `(table, type_label, descriptor_expr)` for the phenomena the
Galaxy Map draws as points (the same labels and descriptors as
`_PHENOMENON_TABLES`)."""


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
        list[dict]: Most luminous first (ties by type, then id): `type`
            (`"black_hole"`, `"neutron_star"` or `"quasar"`), `id`,
            `name`, `descriptor`, `luminosity_sol` and `x`/`y`/`z`
            (parsecs).
    """
    found = []
    for table, type_label, descriptor_expr in _POINT_PHENOMENON_TABLES:
        rows = conn.execute(
            f"""
            SELECT id, name, {descriptor_expr} AS descriptor, luminosity_w,
                   center_x_pc, center_y_pc, center_z_pc
            FROM {table} FORCE INDEX (idx_{table}_center)
            WHERE center_x_pc >= ? AND center_x_pc < ? AND center_y_pc >= ? AND center_y_pc < ?
              AND center_z_pc >= ? AND center_z_pc < ?
            ORDER BY luminosity_w DESC, id
            LIMIT ?
            """,
            (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], int(limit)),
        ).fetchall()
        found += [{
            "type": type_label, "id": row["id"], "name": row["name"], "descriptor": row["descriptor"],
            "luminosity_sol": float("%.4g" % (row["luminosity_w"] / constants.SOLAR_LUMINOSITY)),
            "x": round(row["center_x_pc"], 3), "y": round(row["center_y_pc"], 3), "z": round(row["center_z_pc"], 3),
        } for row in rows]
    found.sort(key=lambda point: (-point["luminosity_sol"], point["type"], point["id"]))
    return found[:limit]


SCATTERED_POINT_CLASSES = (
    ("quasar", None, "quasar", 1.0),
    ("black-hole", "supermassive", "black_hole", 1.0),
    ("black-hole", "intermediate", "black_hole", 0.8),
    ("black-hole", "stellar", "black_hole", 0.55),
    ("neutron-star", None, "neutron_star", 0.35),
)
"""tuple: `(scatter kind, subtype, map type, size share)` for the scattered
objects the Galaxy Map draws as points before their sector is filled
(MAP.164), biggest first. A scatter row stores no mass (the object is built
from its seed when its sector is made), so the share comes from the mass
class the row does carry: 1 for the nucleus black hole down to 0.35 for a
neutron star."""

SCATTERED_POINT_COARSE_CLASSES = 3
"""int: Tiles coarser than `POINT_PHENOMENON_MIN_LEVEL` list only this many
of the biggest `SCATTERED_POINT_CLASSES` (the nucleus, quasar or black hole, and
intermediate-mass black holes), so the largest are visible from the whole-galaxy view."""

GALAXY_TILE_MAX_SCATTERED_POINTS = 60
"""int: Most scattered points one coarse tile lists."""


def galaxy_scattered_points_in_box(conn, lo, hi, edge_pc, limit, coarse):
    """
    The scattered black holes, neutron stars and quasars
    (`phenomenon_scatter`, not yet built into a sector) whose position lies in the box `[lo, hi)`,
    biggest class first, at most `limit` -- drawn as points like the placed
    ones (MAP.164), so the map shows them before their sector is filled.
    `coarse` keeps to the `SCATTERED_POINT_COARSE_CLASSES` biggest classes
    and lists them built or not, since a coarse tile has no placed points
    of its own and the nucleus must not vanish from the galaxy view when
    its sector is filled.

    Returns:
        list[dict]: Same keys as `galaxy_point_phenomena_in_box`, with
            `scattered` (True), `size` (0..1, the class's size share) and
            an `id` that is the scatter row's id as text.
    """
    classes = SCATTERED_POINT_CLASSES[:SCATTERED_POINT_COARSE_CLASSES] if coarse else SCATTERED_POINT_CLASSES
    bands = _bright_star_bands(conn, lo, hi, edge_pc)
    if not bands:
        return []
    if len(bands) > 64:
        where_address = "ring_index BETWEEN ? AND ?"
        address_params = [bands[0][0], bands[-1][0]]
    else:
        where_address = "(" + " OR ".join("(ring_index = ? AND layer_index BETWEEN ? AND ?)" for _ in bands) + ")"
        address_params = [v for band in bands for v in band]
    kinds = " OR ".join("(kind = ? AND subtype <=> ?)" for _ in classes)
    kind_params = [v for cls in classes for v in cls[:2]]
    box_params = [int(math.ceil(v * MPC_PER_PC)) for pair in zip(lo, hi) for v in pair]
    rows = conn.execute(
        f"""
        SELECT id, kind, subtype, built_at, position_x_mpc, position_y_mpc, position_z_mpc
        FROM phenomenon_scatter
        WHERE {"" if coarse else "built_at IS NULL AND "}{where_address} AND ({kinds})
          AND position_x_mpc >= ? AND position_x_mpc < ? AND position_y_mpc >= ? AND position_y_mpc < ?
          AND position_z_mpc >= ? AND position_z_mpc < ?
        ORDER BY (kind = 'quasar') DESC, (subtype <=> 'supermassive') DESC, (subtype <=> 'intermediate') DESC, (subtype <=> 'stellar') DESC, id
        LIMIT ?
        """,
        (*address_params, *kind_params, *box_params, int(limit)),
    ).fetchall()
    share = {(kind, subtype): (map_type, size) for kind, subtype, map_type, size in classes}
    points = []
    for row in rows:
        map_type, size = share[(row["kind"], row["subtype"])]
        label = {"black_hole": "Black hole", "quasar": "Quasar"}.get(map_type, "Neutron star")
        points.append({
            "type": map_type, "id": f"s{row['id']}", "name": f"{label} (uncharted)" if row["built_at"] is None else label,
            "descriptor": row["subtype"] or "scattered", "luminosity_sol": 0.0, "scattered": True, "size": size,
            "x": round(row["position_x_mpc"] / MPC_PER_PC, 3), "y": round(row["position_y_mpc"] / MPC_PER_PC, 3),
            "z": round(row["position_z_mpc"] / MPC_PER_PC, 3),
        })
    points.sort(key=lambda point: (-point["size"], point["id"]))
    return points[:limit]


def galaxy_tiles(conn, tile_keys):
    """
    The contents of each requested cube tile -- the interactive 3D Galaxy
    Map's data source. Every part of
    the result depends only on its tile key and the database's contents
    (see `galaxy_content_stamp`), so callers can cache each part by key.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        tile_keys (list[str]): `"level/ix/iy/iz"` keys (see
            `galaxyViewport.parse_tile_key`), at most
            `MAX_TILES_PER_REQUEST`.

    Returns:
        dict: `tiles` (`{key: {"placed": [...], "planned": [...],
            "filled": {...}, "clouds": [...], "stars": [...],
            "generated": [...], "points": [...]}}`, see `galaxy_sectors_in_box`,
            `galaxyViewport.planned_slots_in_tile`, `galaxy_filled_in_box`,
            `galaxy_clouds_in_box`, `galaxy_bright_stars_in_box` and
            `galaxy_generated_stars_in_box` (its floor from
            `generated_star_floor_sol`, all of them in a finest tile up
            to `GALAXY_TILE_MAX_DETAIL_STARS`) and
            `galaxy_point_phenomena_in_box` (tiles of
            `POINT_PHENOMENON_MIN_LEVEL` and finer)),
            `edge_pc`, `has_shape`. Predicted density isn't served: the
            page evaluates the shape itself (`static/galaxyprisms.js`).

    Raises:
        ValueError: On a malformed key or too many keys.
    """
    parsed = [(key, parse_tile_key(key)) for key in dict.fromkeys(tile_keys)]
    if len(parsed) > MAX_TILES_PER_REQUEST:
        raise ValueError(f"at most {MAX_TILES_PER_REQUEST} tiles per request, got {len(parsed)}")

    skeleton = get_galaxy_shape(conn)
    if skeleton is not None:
        edge_pc = skeleton.edge_pc
        shape = skeleton.shape
        expected_system_count = skeleton.expected_system_count_at_density_1
        bounds = get_galaxy_bounds(conn) if parsed else None
    else:
        edge_pc = ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
        shape = None
        expected_system_count = None
        bounds = None

    cloud_margins = _cloud_margins_pc(conn) if parsed else {}
    sample = []

    def brightest():
        if not sample:
            sample.append(galaxy_brightest_stars(conn))
        return sample[0]

    tiles = {}
    for key, (level, ix, iy, iz) in parsed:
        lo, hi = tile_bounds_pc(level, ix, iy, iz)
        placed = galaxy_sectors_in_box(conn, lo, hi)
        exclude_addresses = {
            address for address in (sector_address(sector) for sector in placed) if address is not None
        }
        planned = planned_slots_in_tile(
            level, ix, iy, iz, edge_pc, shape, expected_system_count, exclude_addresses, bounds,
        )
        filled = galaxy_filled_in_box(conn, lo, hi, hi[0] - lo[0], edge_pc)
        clouds = galaxy_clouds_in_box(conn, lo, hi, margins=cloud_margins)
        stars = galaxy_bright_stars_in_box(conn, lo, hi, edge_pc, brightest=brightest)
        floor = generated_star_floor_sol(level)
        generated = []
        if floor is not None:
            sector_count = len(filled["cells"]) if filled["g"] == 1 else sum(cell[3] for cell in filled["cells"])
            limit, per_sector = generated_star_budget(level, sector_count)
            generated = galaxy_generated_stars_in_box(conn, lo, hi, floor, sector_count, limit=limit,
                                                      per_sector=per_sector)
        points = (galaxy_point_phenomena_in_box(
            conn, lo, hi, limit=GALAXY_TILE_POINT_BUDGET.get(level, GALAXY_TILE_MAX_POINTS))
            if level >= POINT_PHENOMENON_MIN_LEVEL else [])
        # MAP.164: scattered ones not yet built into a sector, the biggest
        # classes at every level so the largest show from the whole galaxy.
        coarse_scatter = level < POINT_PHENOMENON_MIN_LEVEL
        points = points + galaxy_scattered_points_in_box(
            conn, lo, hi, edge_pc,
            GALAXY_TILE_MAX_SCATTERED_POINTS if coarse_scatter else GALAXY_TILE_POINT_BUDGET.get(level, GALAXY_TILE_MAX_POINTS),
            coarse_scatter)
        tiles[key] = {
            "placed": placed, "planned": planned, "filled": filled, "clouds": clouds, "stars": stars,
            "generated": generated, "points": points,
        }

    return {"tiles": tiles, "edge_pc": edge_pc, "has_shape": shape is not None}


def _block_slot_range(block, sector_ring):
    """The slots of `sector_ring` whose centers fall in `block`'s wedge,
    `(first, last)` inclusive -- every drill level's wedges nest, so this
    is exactly the slots whose chain passes through `block`."""
    wedges = drill_wedge_count(block.m, block.ring)
    n = ring_sector_count(sector_ring)
    first = -((wedges - 2 * block.wedge * n) // (2 * wedges))
    last = -((wedges - 2 * (block.wedge + 1) * n) // (2 * wedges)) - 1
    return max(0, first), min(n - 1, last)


def galaxy_stage(conn, at=None):
    """
    One drill-down stage's generated counts (the Galaxy Map drill-down,
    docs/design/galaxy-drilldown-navigation.md section 7): how many
    generated sectors each child block of `at` holds. Totals (allowed
    sectors) are the page's own math, so only generated counts come from
    here, and children with none are left out.

    Without `at`, the galaxy: every level-243 block holding a generated
    sector. With a level-3 `at`, the children are sectors, and `sectors`
    lists each generated one.

    The query reads only the block's own rows: member rings and layers by
    the address index, and in each member ring just the slot range its
    wedge covers (wedges nest at every level, so that range is exact).

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        at (str or None): A block key, `m.ring.wedge.slab` (see
            `galaxyDrill.parse_drill_key`), or `None` for the galaxy.

    Returns:
        dict: `at` (the canonical key, or `None`), `child_m`, `children`
            (`[{ring, wedge, slab, generated, stats}]`, at `child_m = 1` a
            sector's `wedge` is its slot and `slab` its layer), and
            `sectors` (`[{ring, layer, slot, id, name, system_count,
            stats}]` at `child_m = 1`, else `None`). `stats` is what the
            generated sectors hold, from `sector_stats` (DB.14,
            `_stage_stats`): `{systems, expected_systems, stars,
            mean_age_gy, luminosity_sol}`.

    Raises:
        ValueError: On a malformed or impossible `at`.
    """
    if at is None or at == "":
        block = None
        child_m = DRILL_TOP
        rows = conn.execute(
            f"SELECT ring_index, FLOOR((layer_index + ?) / ?) AS slab, ring_slot_index, {_STAGE_LOOK_SUMS} "
            f"FROM ({_stage_sectors_with_stats('ring_index IS NOT NULL')}) placed "
            "GROUP BY ring_index, slab, ring_slot_index",
            ((DRILL_TOP - 1) // 2, DRILL_TOP),
        ).fetchall()
        counts = {}
        for r in rows:
            top = drill_chain_of(int(r["ring_index"]), int(r["slab"]) * DRILL_TOP, int(r["ring_slot_index"]))[0]
            _add_stage_sums(counts, top, r)
        return {"at": None, "child_m": child_m, "children": _stage_children(counts), "sectors": None}

    block = parse_drill_key(at)
    child_m = block.m // (3 if block.m == 3 else 9)
    half = (block.m - 1) // 2
    ring_clauses = []
    params = []
    for i in range(block.ring * block.m, block.ring * block.m + block.m):
        first, last = _block_slot_range(block, i)
        if first <= last:
            ring_clauses.append("(ring_index = ? AND ring_slot_index BETWEEN ? AND ?)")
            params += [i, first, last]
    if not ring_clauses:
        raise ValueError(f"no such block: {at!r}")
    where = f"layer_index BETWEEN ? AND ? AND ({' OR '.join(ring_clauses)})"
    layer_params = [block.slab * block.m - half, block.slab * block.m + half]
    level = (243, 27, 3, 1).index(child_m)

    if child_m == 1:
        rows = conn.execute(
            f"SELECT * FROM ({_stage_sectors_with_stats(where, 'id, name, ')}) placed "
            f"ORDER BY ring_index, layer_index, ring_slot_index",
            layer_params + params,
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
            "ring": r["ring_index"], "layer": r["layer_index"], "slot": r["ring_slot_index"],
            "id": r["id"], "name": r["name"], "system_count": system_counts.get(r["id"], 0),
            "stats": _stage_stats(_sector_sums(r)),
        } for r in rows]
        counts = {}
        for r, sector in zip(rows, sectors):
            _add_stage_sums(counts, DrillBlock(1, sector["ring"], sector["slot"], sector["layer"]), {
                "n": 1, "systems": r["actual_systems"], "expected": r["expected_systems"], "stars": r["actual_stars"],
                "age_sum": (r["mean_age_gy"] or 0.0) * (r["actual_stars"] or 0), "luminosity": r["total_luminosity_sol"],
            })
        return {"at": format_drill_key(block), "child_m": 1, "children": _stage_children(counts), "sectors": sectors}

    child_half = (child_m - 1) // 2
    rows = conn.execute(
        f"SELECT ring_index, FLOOR((layer_index + ?) / ?) AS slab, ring_slot_index, {_STAGE_LOOK_SUMS} "
        f"FROM ({_stage_sectors_with_stats(where)}) placed GROUP BY ring_index, slab, ring_slot_index",
        [child_half, child_m] + layer_params + params,
    ).fetchall()
    counts = {}
    for r in rows:
        child = drill_chain_of(int(r["ring_index"]), int(r["slab"]) * child_m, int(r["ring_slot_index"]))[level]
        _add_stage_sums(counts, child, r)
    return {"at": format_drill_key(block), "child_m": child_m, "children": _stage_children(counts), "sectors": None}


def _sector_sums(row):
    """One sector's sums, in `_add_stage_sums`' order, for its own `stats`."""
    stars = row["actual_stars"] or 0
    return [1, row["actual_systems"] or 0, row["expected_systems"] or 0.0, stars,
            (row["mean_age_gy"] or 0.0) * stars, row["total_luminosity_sol"] or 0.0]


def _stage_sectors_with_stats(where, extra=""):
    """The placed sectors matching `where` (on `sectors`' own columns),
    each with the `sector_stats` the map colors from (DB.14):
    `actual_systems`, `expected_systems`, `actual_stars`, `mean_age_gy`
    and `total_luminosity_sol` (NULL without a row or a fill)."""
    return (
        f"SELECT {extra}s.ring_index, s.layer_index, s.ring_slot_index, st.actual_systems, st.expected_systems,"
        " st.actual_stars, st.mean_age_gy, st.total_luminosity_sol"
        f" FROM (SELECT * FROM sectors WHERE {where}) s LEFT JOIN sector_stats st"
        " ON st.ring_index = s.ring_index AND st.layer_index = s.layer_index"
        " AND st.ring_slot_index = s.ring_slot_index"
    )


_STAGE_LOOK_SUMS = (
    "COUNT(*) AS n, SUM(COALESCE(actual_systems, 0)) AS systems, SUM(COALESCE(expected_systems, 0)) AS expected,"
    " SUM(COALESCE(actual_stars, 0)) AS stars, SUM(COALESCE(mean_age_gy, 0) * COALESCE(actual_stars, 0)) AS age_sum,"
    " SUM(COALESCE(total_luminosity_sol, 0)) AS luminosity"
)
"""The per-group sums `galaxy_stage` folds into each child's stats."""


def _add_stage_sums(counts, child, row):
    """Adds one group's sums (`_STAGE_LOOK_SUMS`) to `counts[child]`:
    `[sectors, systems, expected systems, stars, stars times mean age,
    luminosity]`."""
    sums = counts.setdefault(child, [0, 0.0, 0.0, 0.0, 0.0, 0.0])
    sums[0] += int(row["n"])
    sums[1] += float(row["systems"] or 0.0)
    sums[2] += float(row["expected"] or 0.0)
    sums[3] += float(row["stars"] or 0.0)
    sums[4] += float(row["age_sum"] or 0.0)
    sums[5] += float(row["luminosity"] or 0.0)


def _stage_stats(sums):
    """A child's `stats` from its sums (DB.14): the `systems` and `stars`
    its generated sectors hold, the `expected_systems` the density model
    gave them, `mean_age_gy` (the stars' mean age, each sector's mean
    weighted by its star count; `None` without stars) and
    `luminosity_sol` (the sum of their luminosity)."""
    _n, systems, expected, stars, age_sum, luminosity = sums
    return {
        "systems": int(systems),
        "expected_systems": round(expected, 3),
        "stars": int(stars),
        "mean_age_gy": round(age_sum / stars, 4) if stars else None,
        "luminosity_sol": luminosity,
    }


def _stage_children(counts):
    """`galaxy_stage`'s `children` list from `{DrillBlock: sums}`
    (`_add_stage_sums`)."""
    return [
        {"ring": b.ring, "wedge": b.wedge, "slab": b.slab, "generated": sums[0], "stats": _stage_stats(sums)}
        for b, sums in sorted(counts.items(), key=lambda item: (item[0].slab, item[0].ring, item[0].wedge))
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
    Sectors and star systems whose name contains `term`, with each one's
    sector address, for the Galaxy Map's address bar (the drill-down's
    section 9.3): picking a match flies to that sector. Exact names come
    first, then names that start with `term`, then the rest, by name.
    Sectors without an address (placed before the cylindrical grid) and
    systems outside a sector are left out.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        term (str): Part of a name; blank finds nothing.
        limit (int): Most matches to return.

    Returns:
        list[dict]: `{kind ("sector" or "system"), id, name, sector_id,
            sector_name, ring, layer, slot}`.
    """
    term = (term or "").strip()
    if not term:
        return []
    order = "CASE WHEN {col} = ? THEN 0 WHEN {col} LIKE ? ESCAPE '\\\\' THEN 1 ELSE 2 END, {col}, id"
    prefix = _search_like_pattern(term)[1:]
    sector_match, sector_params = _locate_match(conn, "name", term)
    sectors = conn.execute(
        "SELECT id, name, ring_index, layer_index, ring_slot_index FROM sectors "
        f"WHERE ring_index IS NOT NULL AND {sector_match} "
        f"ORDER BY {order.format(col='name')} LIMIT ?",
        (*sector_params, term, prefix, limit),
    ).fetchall()
    system_match, system_params = _locate_match(conn, "ss.name", term)
    systems = conn.execute(
        "SELECT ss.id, ss.name, s.id AS sector_id, s.name AS sector_name, "
        "s.ring_index, s.layer_index, s.ring_slot_index "
        "FROM star_systems ss JOIN sectors s ON s.id = ss.sector_id "
        f"WHERE s.ring_index IS NOT NULL AND {system_match} "
        f"ORDER BY {order.format(col='ss.name').replace(', id', ', ss.id')} LIMIT ?",
        (*system_params, term, prefix, limit),
    ).fetchall()
    found = [{
        "kind": "sector", "id": r["id"], "name": r["name"], "sector_id": r["id"], "sector_name": r["name"],
        "ring": int(r["ring_index"]), "layer": int(r["layer_index"]), "slot": int(r["ring_slot_index"]),
    } for r in sectors] + [{
        "kind": "system", "id": r["id"], "name": r["name"], "sector_id": r["sector_id"],
        "sector_name": r["sector_name"],
        "ring": int(r["ring_index"]), "layer": int(r["layer_index"]), "slot": int(r["ring_slot_index"]),
    } for r in systems]
    lowered = term.lower()

    def rank(match):
        name = (match["name"] or "").lower()
        return (0 if name == lowered else 1 if name.startswith(lowered) else 2, name, match["kind"], match["id"])

    return sorted(found, key=rank)[:limit]


def galaxy_stage_keys(ring, layer, slot):
    """The stage keys a change to sector `(ring, layer, slot)` makes stale:
    `"galaxy"` and its level-243, 27 and 3 blocks."""
    return ["galaxy"] + [format_drill_key(b) for b in drill_chain_of(ring, layer, slot)[:3]]


GALAXY_CHANGES_MAX_SECTORS = 1000
"""int: Most changed sectors `galaxy_changes` lists tiles for. More than
that (a big generation run, say) reports `full` instead -- refetching
everything is cheaper than invalidating tens of thousands of tiles one
by one."""

_STATE_TOKEN_RE = re.compile(r"^([0-9a-f]{16})\.(\d+)\.(\d+)\.(\d{17}|0)\.(\d+)\.(\d+)$")


def _state_token(state):
    """`galaxy_content_state`'s dict as the opaque string the API hands
    out and `galaxy_changes` reads back."""
    return "{base}.{sectors}.{sector_max_id}.{sector_modified}.{system_max_id}.{bright_max_id}".format(**state)


def _parse_state_token(token):
    """The dict `_state_token` encoded, or `None` for anything else."""
    match = _STATE_TOKEN_RE.match(str(token or ""))
    if not match:
        return None
    base, sectors, sector_max_id, sector_modified, system_max_id, bright_max_id = match.groups()
    return {
        "base": base, "sectors": int(sectors), "sector_max_id": int(sector_max_id),
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


def galaxy_content_state(conn):
    """
    Everything `galaxy_tiles`' output depends on, summarized as a handful
    of numbers:

    - `base`: a hash of the stored galaxy shape, the bright-star scatter
      (`bright_stars`' highest id and the scatter's seed) and this code's
      version. Planned slots, density clouds, `edge_pc` and `has_shape`
      depend on the shape, every tile's `stars` on the scatter, and a
      release may change the tile format, so a new `base` means every
      tile is stale.
    - `sectors`/`sector_max_id`: how many sectors are placed, and the
      highest sector id. Together they tell a new sector (higher id) from
      a deleted one (the count drops).
    - `sector_modified`: the newest `sectors.modified_at` (v27). A
      sector's own edit (a rename, say) bumps it, and so does deleting one
      of its systems (`_db.touch_sector`).
    - `bright_max_id`: the highest `bright_stars` id; it grows while a
      backfill runs, so it is kept out of `base` and `galaxy_changes`
      reports it as `busy` rather than `full` (PERF.34).
    - `system_max_id`: the highest star-system id. A new system changes
      its sector's system count. System edits don't touch a tile (tiles
      only show the count), so `star_systems.modified_at` isn't used.

    Cheap by design -- one indexed count and four index-only maxima --
    since the web layer checks it about once a minute per database.

    Returns:
        dict: The keys above; `sector_modified` as `_timestamp_digits`.
    """
    sector_row = conn.execute(
        "SELECT COUNT(*) AS n FROM sectors WHERE center_x_pc IS NOT NULL"
    ).fetchone()
    maxima = conn.execute(
        "SELECT (SELECT COALESCE(MAX(id), 0) FROM sectors) AS sector_max_id, "
        "(SELECT MAX(modified_at) FROM sectors) AS sector_modified, "
        "(SELECT COALESCE(MAX(id), 0) FROM star_systems) AS system_max_id"
    ).fetchone()
    bright = conn.execute(
        "SELECT (SELECT COALESCE(MAX(id), 0) FROM bright_stars) AS max_id, "
        "(SELECT bright_star_seed FROM galaxy_shape WHERE id = 1) AS seed"
    ).fetchone()
    base = hashlib.sha256(json.dumps(
        {"shape": galaxy_density_shape(conn), "version": __version__,
         "bright_stars": bright["seed"],
         # GEN.70: the naming key renames the codec-named objects in every tile and page.
         "naming_key": naming_key.active_key()}, sort_keys=True, default=str,
    ).encode("utf-8")).hexdigest()[:16]
    return {
        "base": base,
        "sectors": int(sector_row["n"]),
        "sector_max_id": int(maxima["sector_max_id"]),
        "sector_modified": _timestamp_digits(maxima["sector_modified"]),
        "system_max_id": int(maxima["system_max_id"]),
        "bright_max_id": int(bright["max_id"]),
    }


def galaxy_content_stamp(conn, state=None):
    """
    A short token that changes whenever anything `galaxy_tiles` returns
    could change -- a hash of `galaxy_content_state`. The web layer and
    the browser key their tile caches by it, so a stale tile is never
    served after sectors are generated, renamed or deleted.

    Args:
        state (dict, optional): An already-read `galaxy_content_state`.

    Returns:
        str: 16 hex characters.
    """
    if state is None:
        state = galaxy_content_state(conn)
    return hashlib.sha256(_state_token(state).encode("utf-8")).hexdigest()[:16]


def galaxy_changes(conn, since=None):
    """
    What changed in the galaxy's tiles since an earlier state -- how the
    web layer's tile cache refreshes only the cubes an edit touched
    instead of throwing every cached tile away.

    Uses v27's `sectors.modified_at` (see `galaxy_content_state`): each
    sector edited since `since`, each new sector, and each sector that
    gained a system is "changed", and so is the one tile per level that
    holds its center (`galaxyViewport.tile_keys_containing`). A sector
    never moves once placed, so its center is where it was before too.

    Deletions leave no row behind, so a deleted sector can't be located;
    the placed count drops, and the answer is `full` instead. So is a new
    shape or release (`base`), an unreadable `since`, or more than
    `GALAXY_CHANGES_MAX_SECTORS` changed sectors.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        since (str or None): A `state` from an earlier call (or
            `GET /api/galaxy/stamp`), or `None`.

    Returns:
        dict: `stamp` (`galaxy_content_stamp` now), `state` (the token to
            pass as `since` next time), `full` (every tile may have
            changed), `tiles` (sorted keys of the changed tiles; empty
            when `full`), and `stages` (sorted keys of the drill-down
            stages whose counts changed, `galaxy_stage_keys`; empty when
            `full`), and, only when true, `busy`: `full` because so much
            changed at once (a generation run, a backfill), not because
            something was deleted, so a cache may keep its tiles a while
            longer instead of refetching all of them under that load
            (PERF.34).
    """
    state = galaxy_content_state(conn)
    result = {
        "stamp": galaxy_content_stamp(conn, state), "state": _state_token(state), "full": False, "tiles": [],
        "stages": [],
    }
    previous = _parse_state_token(since)
    if previous is None or previous["base"] != state["base"]:
        result["full"] = True
        return result
    if previous == state:
        return result
    if previous["bright_max_id"] != state["bright_max_id"]:
        # A backfill is adding stars to cells all over (PERF.34): which
        # tiles they touch isn't worked out, so everything counts as changed.
        result["full"] = result["busy"] = True
        return result

    new_placed = conn.execute(
        "SELECT COUNT(*) AS n FROM sectors WHERE id > ? AND center_x_pc IS NOT NULL",
        (previous["sector_max_id"],),
    ).fetchone()["n"]
    if previous["sectors"] + int(new_placed) != state["sectors"]:
        result["full"] = True
        return result

    rows = conn.execute(
        """
        SELECT center_x_pc, center_y_pc, center_z_pc, ring_index, layer_index, ring_slot_index
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
        result["full"] = result["busy"] = True
        return result

    keys = set()
    stages = set()
    for r in rows:
        keys.update(tile_keys_containing((r["center_x_pc"], r["center_y_pc"], r["center_z_pc"])))
        if sector_address(r) is not None:
            stages.update(galaxy_stage_keys(*sector_address(r)))
    result["tiles"] = sorted(keys)
    result["stages"] = sorted(stages)
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

# starData.py's own Yerkes-class-to-descriptive-label mapping (Star.__init__),
# duplicated here (not imported) since it's a plain literal there, not an
# importable constant. "D" is an alternate white-dwarf code seen elsewhere
# in starData.py alongside "VII".
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
    rows = conn.execute(
        "SELECT SUBSTR(star_type, 1, 1) AS v, COUNT(*) AS c FROM stars GROUP BY v ORDER BY v"
    ).fetchall()
    opts = []
    for row in rows:
        letter = row["v"]
        color = SPECTRAL_CLASS_COLORS.get(letter)
        tip = f"{color} star" if color else None
        opts.append({"value": letter, "label": f"{letter}-Type Star", "count": row["c"], "tooltip": tip})
    return opts


def _search_facet_luminosity(conn):
    rows = conn.execute(
        "SELECT yerkes_class AS v, COUNT(*) AS c FROM stars WHERE yerkes_class IS NOT NULL GROUP BY v"
    ).fetchall()
    rows = sorted(
        rows, key=lambda r: _SEARCH_YERKES_ORDER.index(r["v"]) if r["v"] in _SEARCH_YERKES_ORDER else len(_SEARCH_YERKES_ORDER)
    )
    return [
        {"value": row["v"], "label": _SEARCH_YERKES_LABELS.get(row["v"], row["v"]), "count": row["c"],
         "tooltip": f"Yerkes class {row['v']}"}
        for row in rows
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


def _search_facet_life(conn):
    rows = conn.execute(
        "SELECT life_chemical AS v, COUNT(*) AS c FROM planets WHERE life_chemical IS NOT NULL GROUP BY v ORDER BY v"
    ).fetchall()
    return [{"value": row["v"], "label": row["v"], "count": row["c"], "tooltip": None} for row in rows]


def _search_equipment_options(conn, table):
    """GEN.89: the five equipment tiers a human needs (value "0" to "4"),
    with how many bodies of `table` need each; the ones none need are left out."""
    rows = conn.execute(
        f"SELECT equipment_tier AS v, COUNT(*) AS c FROM {table} WHERE equipment_tier IS NOT NULL"
        " GROUP BY v ORDER BY v").fetchall()
    return [{"value": str(row["v"]), "label": EQUIPMENT_LABELS[row["v"]], "count": row["c"],
             "tooltip": ("Conditions are ideal for a human here" if row["v"] == 0
                         else f"A human needs {EQUIPMENT_NAMES[row['v']]} here")} for row in rows]


def _search_facet_equipment(conn):
    return _search_equipment_options(conn, "planets")


def _search_facet_moon_equipment(conn):
    return _search_equipment_options(conn, "moons")


def _equipment_clause(alias, tags, clauses, params):
    """Adds `equipment_tier IN (...)` for the tags that are tiers (a stray value is dropped)."""
    tiers = sorted({int(tag) for tag in tags if tag in {str(i) for i in range(len(EQUIPMENT_NAMES))}})
    if tags:
        clauses.append(f"{alias}.equipment_tier IN ({','.join('?' * len(tiers))})" if tiers else "1 = 0")
        params.extend(tiers)


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
    rows = conn.execute(
        "SELECT life_chemical AS v, COUNT(*) AS c FROM moons WHERE life_chemical IS NOT NULL GROUP BY v ORDER BY v"
    ).fetchall()
    return [{"value": row["v"], "label": row["v"], "count": row["c"], "tooltip": None} for row in rows]


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
    return {"rows": [dict(r) for r in rows], **page}


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


def _search_result_sectors(conn, term, limit, offset):
    match, params = _name_match(conn, "name", term)
    rows, page = _search_page(
        conn, "SELECT id, name, edge_mpc", f"FROM sectors WHERE {match}", "name, id", params, limit, offset,
    )
    return {"rows": [dict(r) for r in rows], **page}


def _search_result_systems(conn, term, limit, offset):
    match, params = _name_match(conn, "ss.name", term)
    rows, page = _search_page(
        conn,
        "SELECT ss.id, ss.name, ss.sector_id, ss.is_binary, ss.binary_configuration, ss.binary_type",
        f"FROM star_systems ss WHERE {match}",
        "ss.name, ss.id",
        params, limit, offset,
    )
    return {
        "rows": [
            {
                "id": r["id"], "name": r["name"], "sector_id": r["sector_id"], "is_binary": r["is_binary"],
                "star_summary": _star_summary(r),
            }
            for r in _with_star_types(conn, rows)
        ],
        **page,
    }


def _search_result_stars(conn, spectral_tags, luminosity_tags, term, limit, offset, size_range=None):
    clauses, params = [], []
    if spectral_tags:
        clauses.append(f"SUBSTR(s.star_type, 1, 1) IN ({','.join('?' * len(spectral_tags))})")
        params.extend(sorted(spectral_tags))
    if luminosity_tags:
        clauses.append(f"s.yerkes_class IN ({','.join('?' * len(luminosity_tags))})")
        params.extend(sorted(luminosity_tags))
    _append_size_clause(clauses, params, "s.radius_km", size_range)
    if term:
        match, match_params = _name_match(conn, "s.name", term)
        clauses.append(match)
        params.extend(match_params)
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    rows, page = _search_page(
        conn,
        "SELECT s.name, s.role, s.star_type, s.radius_km, s.star_system_id, ss.name AS system_name, ss.sector_id",
        f"FROM stars s JOIN star_systems ss ON ss.id = s.star_system_id WHERE 1=1{where}",
        "ss.name, s.name, s.id",
        params, limit, offset,
    )
    return {"rows": [dict(r) for r in rows], **page}


def _search_result_planets(conn, class_tags, body_tags, life_tags, term, limit, offset, size_range=None,
                           equipment_tags=()):
    clauses, params = [], []
    _equipment_clause("p", equipment_tags, clauses, params)
    if class_tags:
        clauses.append(f"p.planet_class IN ({','.join('?' * len(class_tags))})")
        params.extend(sorted(class_tags))
    if body_tags:
        clauses.append(f"p.body_type IN ({','.join('?' * len(body_tags))})")
        params.extend(sorted(body_tags))
    if life_tags:
        clauses.append(f"p.life_chemical IN ({','.join('?' * len(life_tags))})")
        params.extend(sorted(life_tags))
    _append_size_clause(clauses, params, "p.radius_km", size_range)
    if term:
        match, match_params = _name_match(conn, "p.name", term)
        clauses.append(match)
        params.extend(match_params)
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    rows, page = _search_page(
        conn,
        """
        SELECT p.name, p.planet_class, p.body_type, p.life_chemical, p.equipment_tier, p.radius_km,
               p.star_system_id, ss.name AS system_name, ss.sector_id
        """,
        f"FROM planets p JOIN star_systems ss ON ss.id = p.star_system_id WHERE 1=1{where}",
        "ss.name, ss.id, p.orbital_index, p.id",
        params, limit, offset,
    )
    return {"rows": [dict(r) for r in rows], **page}


def _search_result_moons(conn, class_tags, body_tags, life_tags, term, limit, offset, size_range=None,
                         equipment_tags=()):
    clauses, params = [], []
    _equipment_clause("m", equipment_tags, clauses, params)
    if class_tags:
        clauses.append(f"m.planet_class IN ({','.join('?' * len(class_tags))})")
        params.extend(sorted(class_tags))
    if body_tags:
        clauses.append(f"m.body_type IN ({','.join('?' * len(body_tags))})")
        params.extend(sorted(body_tags))
    if life_tags:
        clauses.append(f"m.life_chemical IN ({','.join('?' * len(life_tags))})")
        params.extend(sorted(life_tags))
    _append_size_clause(clauses, params, "m.radius_km", size_range)
    if term:
        match, match_params = _name_match(conn, "m.name", term)
        clauses.append(match)
        params.extend(match_params)
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    rows, page = _search_page(
        conn,
        """
        SELECT m.name, m.planet_class, m.body_type, m.life_chemical, m.equipment_tier, m.radius_km, p.name AS planet_name,
               m.star_system_id, ss.name AS system_name, ss.sector_id
        """,
        f"""
        FROM moons m
        JOIN planets p ON p.id = m.planet_id
        JOIN star_systems ss ON ss.id = m.star_system_id
        WHERE 1=1{where}
        """,
        "ss.name, ss.id, p.orbital_index, m.orbital_index, m.id",
        params, limit, offset,
    )
    return {"rows": [dict(r) for r in rows], **page}


def _search_result_belts(conn, density_tags, limit, offset):
    clauses, params = [], []
    if density_tags:
        clauses.append(f"ab.density IN ({','.join('?' * len(density_tags))})")
        params.extend(sorted(density_tags))
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    rows, page = _search_page(
        conn,
        "SELECT ab.density, ab.composition_summary, ab.star_system_id, ss.name AS system_name, ss.sector_id",
        f"FROM asteroid_belts ab JOIN star_systems ss ON ss.id = ab.star_system_id WHERE 1=1{where}",
        "ss.name, ss.id, ab.orbital_index, ab.id",
        params, limit, offset,
    )
    return {"rows": [dict(r) for r in rows], **page}


SEARCH_FACET_CACHE_SECONDS = 600
"""int: Longest the search page's facet counts and name lists are reused
(PERF.15) even when nothing they depend on looks changed -- a backstop
for edits `_search_content_key` can't see."""

_facet_cache = {}
_facet_cache_lock = threading.Lock()


def _search_content_key(conn):
    """One cheap row (index-only maxima and one count) that changes
    whenever a sector or system is added, edited or deleted -- what the
    facet counts depend on."""
    row = conn.execute(
        "SELECT (SELECT COUNT(*) FROM sectors) AS sectors, "
        "(SELECT COALESCE(MAX(id), 0) FROM sectors) AS sector_max_id, "
        "(SELECT MAX(modified_at) FROM sectors) AS sector_modified, "
        "(SELECT COALESCE(MAX(id), 0) FROM star_systems) AS system_max_id, "
        "(SELECT MAX(modified_at) FROM star_systems) AS system_modified"
    ).fetchone()
    return tuple(str(row[column]) for column in
                 ("sectors", "sector_max_id", "sector_modified", "system_max_id", "system_modified"))


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


def search(conn, texts, tags, sizes=None, limit=SEARCH_RESULT_LIMIT, offsets=None, only=None):
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

    Returns:
        dict: `facets` (`{facet: [{"value","label","count","tooltip"}, ...]}`,
            one list per `SEARCH_TAG_FACETS` name), `autocomplete`
            (`sectors`/`systems`/`stars`/`planets`/`moons` -> list of
            distinct names), `facet_labels` (`{"facet:value": label}`,
            for rendering an active-filter chip without a second lookup),
            and `results` (`sectors`/`systems`/`stars`/`planets`/`moons`/
            `belts` -> `{"rows": [...], "total", "limit", "offset",
            "truncated"}`, one page of that panel's matches -- `truncated`
            meaning `rows` isn't every match -- or `None` for a panel
            with no reason to run).
    """
    sizes = sizes or {}
    star_size, planet_size, moon_size = sizes.get("star"), sizes.get("planet"), sizes.get("moon")

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

    star_has_reason = bool(spectral_tags or luminosity_tags or texts.get("star_q") or star_size)
    planet_has_reason = bool(class_tags or body_tags or life_tags or equipment_tags or texts.get("planet_q") or planet_size)
    moon_has_reason = bool(moon_class_tags or moon_body_tags or moon_life_tags or moon_equipment_tags or texts.get("moon_q") or moon_size)
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
    if texts.get("system_q") and "systems" in wanted:
        results["systems"] = _search_result_systems(conn, texts["system_q"], *_page("systems"))
    if stars_included and "stars" in wanted:
        results["stars"] = _search_result_stars(
            conn, spectral_tags, luminosity_tags, texts.get("star_q", ""), *_page("stars"), size_range=star_size
        )
    if planets_included and "planets" in wanted:
        results["planets"] = _search_result_planets(
            conn, class_tags, body_tags, life_tags, texts.get("planet_q", ""), *_page("planets"),
            size_range=planet_size, equipment_tags=equipment_tags,
        )
    if moons_included and "moons" in wanted:
        results["moons"] = _search_result_moons(
            conn, moon_class_tags, moon_body_tags, moon_life_tags, texts.get("moon_q", ""), *_page("moons"),
            size_range=moon_size, equipment_tags=moon_equipment_tags,
        )
    if belts_included and "belts" in wanted:
        results["belts"] = _search_result_belts(conn, density_tags, *_page("belts"))
    phenomenon_tags, phenomenon_class_tags = tags.get("phenomenon", set()), tags.get("phenomenon_class", set())
    if (phenomenon_tags or phenomenon_class_tags) and "phenomena" in wanted:
        results["phenomena"] = _search_result_phenomena(
            conn, phenomenon_tags, phenomenon_class_tags, *_page("phenomena"))

    return {"facets": facets, "autocomplete": autocomplete, "facet_labels": facet_labels, "results": results}
