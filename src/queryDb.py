# src/queryDb.py

"""
List/query CLI for the planetGen database (src/stellarObjects/schema.sql).

A thin read-only front end over the tables `sectorGen.py`/`systemGen.py`
populate, for questions like "every G-type system," "everything within 50
light-years of a given system," and "what sectors exist" -- see
docs/TODO.md's Phase 3 ("A way to list/query what's already stored").
Deliberately plain SQL rather than routing through `stellarObjects._db`'s
`load_star_system`/`load_sector` (Phase 2's read path): these are simple,
columnar listings, not full object-graph reconstructions, so a raw query
is the more direct tool for the job -- the read path remains what a
future richer tool (or a re-upload/re-render workflow) would build on.

This tool never writes -- `open_readonly` below connects with the same
`MySQLConfig` every other entry point uses, but the actual enforcement
that the connection can't write is a deployment concern: point
`PLANETGEN_MYSQL_USER`/`PLANETGEN_MYSQL_PASSWORD` at a database account
with `SELECT`-only grants for this tool rather than the read-write account
`generate.py` and the Flask app use (the app writes too: the admin pages
and the API's write endpoints, see `html/api/config.py`) -- MySQL has no per-connection "open this read-only"
flag the way SQLite's `file:...?mode=ro` URI trick gave the old SQLite
version of this function, so the guarantee lives in the account's grants
instead of the connection itself.

Run directly as `python src/queryDb.py`: this file lives alongside
`stellarObjects/` under `src/`, so Python's own sys.path[0] (the running
script's directory) already makes `stellarObjects` importable -- no
sys.path shim needed, unlike the root-level entry scripts
(`sectorGen.py`/`systemGen.py`) that stay one directory further away.
"""

import argparse
import datetime
import hashlib
import json
import math
import re
import threading
import time

import pymysql

from stellarObjects._db import (add_mysql_connection_args, escape_like, get_connection, get_galaxy_shape,
                                mysql_config_from_args, surrounding_cloud)
from stellarObjects import physical_constants, program_constants
from stellarObjects.starData import compressed_heliosphere_radius
from stellarObjects.brightStars import MPC_PER_PC
from stellarObjects._version import VersionAction, __version__, version_banner
from stellarObjects.galaxyGeometry import (
    galaxy_to_local_pc, neighbor_addresses, provisional_sector_designation, ring_sector_count, sector_cell_vertices_pc,
    sector_position_pc,
)
from stellarObjects.spaceSector import classify_octant
from stellarObjects.galaxyViewport import (
    TILE_MAX_LEVEL,
    TILE_ROOT_EDGE_PC,
    parse_tile_key,
    planned_slots_in_tile,
    tile_bounds_pc,
    tile_keys_containing,
)
from stellarObjects.galaxyDrill import (
    DRILL_TOP, DrillBlock, drill_chain_of, drill_wedge_count, format_drill_key, parse_drill_key,
)
from stellarObjects.navGraph import build_knn_adjacency, shortest_path
from stellarObjects.navigation import (
    FRAME_GALACTIC, FRAME_SECTOR, course_between, fold_travel_times, warp_travel_times,
)
from stellarObjects.physical_constants import SPECTRAL_CLASS_COLORS
from stellarObjects.evolution import life_stage_from_paragraphs
from stellarObjects.program_constants import (
    DEFAULT_SECTOR_EDGE_LY, HABITABLE_PLANET_CLASSES, NAV_ADJACENCY_K, PLANET_CLASSES,
)
from stellarObjects.utils import ly_to_pc, milliparsecs_to_ly, mpc_to_pc, pc_to_ly


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
        stellarObjects._db.Connection: An open connection.

    Raises:
        SystemExit: If the database can't be reached.
    """
    try:
        return get_connection(config, ensure_schema=False, statement_timeout_s=statement_timeout_s)
    except pymysql.MySQLError as exc:
        raise SystemExit(f"Error: could not open the database ({exc}).")


def list_sectors(conn, limit=None, offset=None):
    """
    Returns every sector, with its edge length (converted to light-years)
    and how many systems it contains, nearest the galactic core first
    (`galactic_radius_pc`); sectors never placed in a galaxy have no
    distance and come last, by name.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        limit (int, optional): Caps the number of rows returned. `None`
            (the default -- what every existing caller of this function
            still gets) returns every sector.
        offset (int, optional): Skips this many rows first. Ignored
            unless `limit` is also given; meaningless on its own.

    Returns:
        list[dict]: One row per sector, with `id`, `name`,
                           `edge_ly`, `system_count`,
                           `galactic_radius_pc`/`galactic_radius_ly`
                           (`None` if unplaced).
    """
    query = """
        SELECT sec.id, sec.name, sec.edge_mpc, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.galactic_radius_pc, sec.ring_index, sec.layer_index, sec.ring_slot_index, COUNT(ss.id) AS system_count
        FROM sectors sec
        LEFT JOIN star_systems ss ON ss.sector_id = sec.id
        -- Every selected column, not just the key: MariaDB's
        -- ONLY_FULL_GROUP_BY doesn't see columns that depend on sec.id.
        GROUP BY sec.id, sec.name, sec.edge_mpc, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
                 sec.galactic_radius_pc, sec.ring_index, sec.layer_index, sec.ring_slot_index
        ORDER BY sec.galactic_radius_pc IS NULL, sec.galactic_radius_pc, sec.name, sec.id
        """
    params = []
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


def count_sectors(conn):
    """
    Returns the total number of sectors, ignoring any pagination --
    the denominator `list_sectors(conn, limit=...)` callers (the API's
    `/api/sectors`) need to report how many pages exist.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.

    Returns:
        int: Total sector count.
    """
    return conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"]


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


def _systems_filter_clause(star_type_prefix, sector_id):
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

    where_sql = (" WHERE " + " AND ".join(conditions)) if conditions else ""
    return join_sql, where_sql, params


def list_systems(conn, star_type_prefix=None, sector_id=None, limit=None, offset=None):
    """
    Returns systems, optionally filtered by star type and/or sector.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
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

    Returns:
        list[dict]: One row per matching system, with `id`, `name`,
                           `sector_id`, `sector_name`, `quadrant` (its
                           octant in the sector), `is_binary`, and `star_summary`
                           (the single star's `star_type`, or a binary's
                           `binary_type` -- what `html/browse.py`/
                           `html/search.py` show as a system's "Star type"
                           column).
    """
    join_sql, where_sql, params = _systems_filter_clause(star_type_prefix, sector_id)
    query = f"""
        SELECT DISTINCT ss.id, ss.name, ss.sector_id, ss.quadrant, ss.is_binary, ss.binary_configuration,
               ss.binary_type
        FROM star_systems ss{join_sql}{where_sql} ORDER BY ss.name, ss.id
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


def count_systems(conn, star_type_prefix=None, sector_id=None):
    """
    Returns the total number of systems matching the same filters
    `list_systems` accepts, ignoring any pagination -- the denominator
    `list_systems(conn, limit=...)` callers (the API's `/api/systems`)
    need to report how many pages exist.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        star_type_prefix (str, optional): Same meaning as `list_systems`.
        sector_id (int or NO_SECTOR, optional): Same meaning as `list_systems`.

    Returns:
        int: Total matching system count.
    """
    join_sql, where_sql, params = _systems_filter_clause(star_type_prefix, sector_id)
    query = f"SELECT COUNT(DISTINCT ss.id) AS n FROM star_systems ss{join_sql}{where_sql}"
    return conn.execute(query, params).fetchone()["n"]


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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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


def systems_within_radius(conn, system_id, radius_ly):
    """
    Finds every other system in the same sector as `system_id`, within
    `radius_ly` light-years, nearest first.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        system_id (int): The `star_systems.id` to measure distances from.
        radius_ly (float): The search radius, in light-years.

    Returns:
        list[dict]: One entry per match, nearest first, each with `id`,
                   `name`, `distance_ly`.

    Raises:
        SystemExit: If `system_id` doesn't exist or isn't placed in a
                   sector (no position to measure from).
    """
    origin = conn.execute(
        "SELECT sector_id, position_x_mpc, position_y_mpc, position_z_mpc "
        "FROM star_systems WHERE id = ?",
        (system_id,),
    ).fetchone()
    if origin is None:
        raise SystemExit(f"Error: no star_systems row with id {system_id}.")
    if origin["sector_id"] is None or origin["position_x_mpc"] is None:
        raise SystemExit(f"Error: system {system_id} isn't placed in a sector (no position to measure from).")

    origin_ly = (
        milliparsecs_to_ly(origin["position_x_mpc"]),
        milliparsecs_to_ly(origin["position_y_mpc"]),
        milliparsecs_to_ly(origin["position_z_mpc"]),
    )

    candidates = conn.execute(
        "SELECT id, name, position_x_mpc, position_y_mpc, position_z_mpc "
        "FROM star_systems WHERE sector_id = ? AND id != ?",
        (origin["sector_id"], system_id),
    ).fetchall()

    results = []
    for row in candidates:
        candidate_ly = (
            milliparsecs_to_ly(row["position_x_mpc"]),
            milliparsecs_to_ly(row["position_y_mpc"]),
            milliparsecs_to_ly(row["position_z_mpc"]),
        )
        distance_ly = math.dist(origin_ly, candidate_ly)
        if distance_ly <= radius_ly:
            results.append({"id": row["id"], "name": row["name"], "distance_ly": distance_ly})

    results.sort(key=lambda entry: entry["distance_ly"])
    return results


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
    here, since light-years is the unit `stellarObjects.navigation`
    already works in for sector-local distances (see
    `spaceSector.distance_between`).

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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


def _galaxy_frame_positions(conn):
    """
    Returns every placed system's absolute galaxy-frame position (in
    light-years) across every galaxy-placed sector -- the position set a
    cross-sector NAV route's adjacency graph is built from. Systems in a
    sector with no galaxy placement (`sectors.center_x_pc IS NULL`, e.g. a
    standalone sector in a database with no galaxy at all) are excluded,
    same as an unplaced system within a sector -- neither has an absolute
    position to route through.

    This necessarily only sees sectors that have actually been generated
    and stored (see `galaxyGen.ensure_sector_generated`'s lazy generation),
    not every sector a galaxy's skeleton says *could* exist -- there is no
    position to route through for a sector nothing has visited yet either.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.

    Returns:
        dict: `{star_systems.id: (x, y, z)}`, light-years, galaxy-frame.
    """
    rows = conn.execute(
        """
        SELECT ss.id, ss.position_x_mpc, ss.position_y_mpc, ss.position_z_mpc,
               sec.center_x_pc, sec.center_y_pc, sec.center_z_pc
        FROM star_systems ss
        JOIN sectors sec ON sec.id = ss.sector_id
        WHERE ss.position_x_mpc IS NOT NULL AND sec.center_x_pc IS NOT NULL
        """
    ).fetchall()
    positions = {}
    for row in rows:
        positions[row["id"]] = (
            pc_to_ly(row["center_x_pc"]) + milliparsecs_to_ly(row["position_x_mpc"]),
            pc_to_ly(row["center_y_pc"]) + milliparsecs_to_ly(row["position_y_mpc"]),
            pc_to_ly(row["center_z_pc"]) + milliparsecs_to_ly(row["position_z_mpc"]),
        )
    return positions


def nav_between(conn, from_id, to_id, adjacency_k=NAV_ADJACENCY_K,
                 from_kind="system", to_kind="system", from_type=None, to_type=None):
    """
    Resolves full NAV information between two endpoints -- each either a
    star system or a standalone phenomenon (nebula/asteroid field/black
    hole/neutron star) -- a direct course (distance/bearing/mark/warp and fold
    travel times, from `stellarObjects.navigation`) plus an optimal route
    via adjacent systems (`stellarObjects.navGraph`), or raises if NAV
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
        from_id (int): The origin's own row id -- `star_systems.id` when
            `from_kind == "system"`, else the phenomenon's own table id.
        to_id (int): Same, for the destination.
        adjacency_k (int): Passed through to
            `navGraph.build_knn_adjacency` as `k`.
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
            to the same node or no path exists through the adjacency
            graph, else `{"path": [...node ids...], "distance_ly": float,
            "positions": {node_id: (x, y, z), ...}}` (one entry per id in
            `path`, same frame as `origin_position`/`destination_position`
            -- for rendering the route, e.g. `html/lib/navmap.py`, without
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

    if sector_scope_ok:
        scope = "sector"
        positions = _sector_local_positions(conn, origin["sector_id"])
        origin_position, destination_position = origin["position_ly"], destination["position_ly"]
    elif galaxy_scope_ok:
        scope = "galaxy"
        positions = _galaxy_frame_positions(conn)
        # A phenomenon endpoint is never itself a row _galaxy_frame_positions
        # reads (it isn't a star_systems row at all) -- added as this one-
        # off query's own extra graph node instead, under its own
        # _phenomenon_nav_key so it can't collide with any real system id.
        if from_kind == "phenomenon":
            positions[from_key] = origin["galaxy_position_ly"]
        if to_kind == "phenomenon":
            positions[to_key] = destination["galaxy_position_ly"]
        origin_position, destination_position = origin["galaxy_position_ly"], destination["galaxy_position_ly"]
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
        graph = build_knn_adjacency(positions, adjacency_k)
        found = shortest_path(graph, from_key, to_key)
        if found is not None:
            path, distance_ly = found
            route = {
                "path": path,
                "distance_ly": distance_ly,
                "positions": {node_id: positions[node_id] for node_id in path},
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
    exists there. Drives the Sector Map's (`html/lib/starmap.py`)
    neighboring-sector indicators -- an existing neighbor links straight
    to it; a not-yet-generated one shows its address so it can be fed to
    `generate.py galaxy --ring I --layer J --slot K`.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
        "host_name",
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
            "name", "class", "density_cm3", "temperature_k"}`, or `None`
            in open space.
    """
    return surrounding_cloud(conn, row)


def sector_detail(conn, sector_id):
    """
    Returns one sector's full web-display detail: name, size, galaxy
    placement, and every system placed in it (each with its own star
    roster) -- everything `html/sector.py`'s systems table and Sector Map
    (`html/lib/starmap.py`) need, in one function.

    Distinct from `stellarObjects._db.load_sector`, which reconstructs
    the *generation* object graph (config/provenance, no database ids) --
    this is a flat, ids-and-display-fields read, the same relationship
    `list_sectors`/`list_systems` above already have to
    `stellarObjects._db.load_star_system`.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
            to, a wiki page; see `schema.sql`'s "v22" header note).

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
    }


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
    return (program_constants.INTERSTELLAR_DEBRIS_DENSITY_PC3
            / program_constants.REFERENCE_STELLAR_DENSITY_PC3 * star_count)


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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
        conn (stellarObjects._db.Connection): An open, read-only connection.

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


def list_phenomena(conn, limit=None, offset=None):
    """
    Returns every exotic phenomenon (nebula/asteroid field/black hole/
    neutron star/supernova remnant/rogue planet/interstellar comet --
    every table in `_PHENOMENON_TABLES`), across every sector and
    regardless of galaxy placement -- `GET /api/phenomena`'s own flat
    listing (`html/phenomena.py`), unlike `galaxy_placed_phenomena` (which
    only returns the galaxy-placed subset, for the Galaxy Map) or
    `phenomena_near_sector` (one sector's own neighborhood).

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        limit (int, optional): Caps the number of rows returned.
        offset (int, optional): Skips this many rows first. Ignored unless
            `limit` is also given.

    Excludes a `black_holes`/`neutron_stars` row with `star_id` set -- that
    shape is a normal star system's own compact-remnant star (already
    shown on that system's own `system.py` page), not a standalone exotic
    phenomenon; every other table here has no `star_id` at all (always
    standalone, see their own table comments) and needs no such filter.

    Returns:
        list[dict]: One row per phenomenon, ordered by name: `id`, `type`
            (`"nebula"`, `"asteroid_field"`, `"black_hole"`,
            `"neutron_star"`, `"supernova_remnant"`, `"rogue_planet"`, or
            `"interstellar_comet"`), `name`, `descriptor`, `radius_ly`,
            `sector_id`/`sector_name` (both `None` if this phenomenon has
            never been linked to a sector -- see `schema.sql`'s "v18"
            header note), and `placed` (bool -- whether it has a galaxy
            position at all, `center_x_pc IS NOT NULL`).
    """
    union_parts = [
        f"""
        SELECT '{type_label}' AS type, t.id AS id, t.name AS name,
               {descriptor_expr} AS descriptor, {radius_expr} AS radius_ly,
               t.sector_id AS sector_id, sec.name AS sector_name,
               t.center_x_pc AS center_x_pc
        FROM {table} t
        LEFT JOIN sectors sec ON sec.id = t.sector_id
        {"WHERE t.star_id IS NULL" if table in ("black_holes", "neutron_stars") else ""}
        """
        for table, type_label, descriptor_expr, radius_expr in _PHENOMENON_TABLES
    ]
    query = "SELECT * FROM (" + " UNION ALL ".join(union_parts) + ") AS phenomena ORDER BY name"
    params = []
    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params.extend([limit, offset or 0])

    rows = conn.execute(query, params).fetchall()
    return [
        {
            "id": row["id"], "type": row["type"], "name": row["name"],
            "descriptor": row["descriptor"], "radius_ly": row["radius_ly"],
            "sector_id": row["sector_id"], "sector_name": row["sector_name"],
            "placed": row["center_x_pc"] is not None,
        }
        for row in rows
    ]


def count_phenomena(conn):
    """
    Returns the total number of exotic phenomena across every type in
    `_PHENOMENON_TABLES`, ignoring any
    pagination -- the denominator `list_phenomena(conn, limit=...)`
    callers (the API's `/api/phenomena`) need to report how many pages
    exist.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.

    Returns:
        int: Total phenomenon count.
    """
    return sum(
        conn.execute(
            f"SELECT COUNT(*) AS n FROM {table}"
            + (" WHERE star_id IS NULL" if table in ("black_holes", "neutron_stars") else "")
        ).fetchone()["n"]
        for table, _type_label, _descriptor_expr, _radius_expr in _PHENOMENON_TABLES
    )


_PHENOMENON_TYPE_TO_TABLE = {
    type_label: table for table, type_label, _de, _re in _PHENOMENON_TABLES
}
"""dict: `type` value (as returned by `list_phenomena`/`galaxy_placed_phenomena`)
-> its backing table name, e.g. `"nebula"` -> `"nebulae"` -- the reverse of
`_PHENOMENON_TABLES`'s own `(table, type_label, ...)` order, used by
`phenomenon_detail` to find the one table a `(type, id)` pair actually
means without hand-listing the mapping a second time."""


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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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

    detail = dict(row)
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
    dots on the same Galaxy Map (`html/lib/galaxymap.py`).

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.

    Returns:
        list[dict]: See `_placed_phenomenon_rows`.
    """
    return _placed_phenomenon_rows(conn)


def nearest_systems(conn, object_table, object_ids):
    """
    The stored nearest star systems (`nearest_systems`, schema v41) of
    several objects of one kind.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
    data `html/lib/starmap.py`'s Sector Map draws (translucent clouds for
    nebulae/asteroid fields/supernova remnants, point markers for the
    point-like types, whose own `radius_ly` is always 0 -- see
    `_PHENOMENON_TABLES`) and `html/web/sector_page.py` lists alongside
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
    open_space_au = radius_km / physical_constants.AU_TO_KM
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
    `stellarObjects._db.load_star_system`'s generation object graph).

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
    # mass_kg (already selected above) is what lets html/lib/systemmap.py
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
        planet_dict = _with_life_fields(dict(planet), planet_stages, colonized["planets"])
        planet_dict["moons"] = [_with_life_fields(dict(m), moon_stages, colonized["moons"])
                                for m in moons_by_planet.get(planet["id"], [])]
        planets.append(planet_dict)

    belts = [dict(b) for b in conn.execute(
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
        "comets": [dict(c) for c in comets],
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
    `program_constants.HABITABLE_PLANET_CLASSES`, the same test
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
    data the `/galaxy` Galaxy Map (`html/lib/galaxymap.py`) plots.
    Unplaced sectors (`center_x_pc IS NULL`) have nothing to plot and are
    excluded at the query itself.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.

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


def galaxy_density_shape(conn):
    """
    The galaxy's stored density-skeleton shape (the singleton
    `galaxy_shape` row `generate.py plan` writes -- see
    `stellarObjects.galaxyDensity.GalaxyShape` and `_db.get_galaxy_shape`),
    serialized to a plain JSON-able dict. This is the real
    exponential-disk-plus-bulge-plus-spiral-arm model already used to gate
    and weight actual sector generation (`generate.py`'s `_BatchDensity`);
    exposing it here lets the Galaxy Map (`html/lib/galaxymap.py`) shade
    its "expected density" cloud from this same model instead of a
    generic illustrative gradient, so un-generated space still reads as
    the spiral it's predicted to be.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.

    Returns:
        dict or None: Every `GalaxyShape` field plus `edge_pc`,
            `outer_ring_index`, `expected_system_count_at_density_1`.
            `None` if `generate.py plan` has never been run against this
            database (no `galaxy_shape` row yet).
    """
    skeleton = get_galaxy_shape(conn)
    if skeleton is None:
        return None
    shape = dict(skeleton.shape._asdict())
    shape["edge_pc"] = skeleton.edge_pc
    shape["outer_ring_index"] = skeleton.outer_ring_index
    shape["expected_system_count_at_density_1"] = skeleton.expected_system_count_at_density_1
    return shape


# ---------------------------------------------------------------------
# Interactive 3D Galaxy Map viewport queries -- the placed-sector shape
# the cube tiles below reuse. Unlike `galaxy_placed_sectors`/`galaxy_density_shape` above (each
# called once per page load for the flat, whole-galaxy overview map),
# these are scoped to a moving viewport -- see `stellarObjects.
# galaxyViewport`'s own module docstring for the three content tiers
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
    interactive 3D Galaxy Map's live viewport (see `stellarObjects.
    galaxyViewport`'s module docstring), unlike `galaxy_placed_sectors`
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
# `stellarObjects.galaxyViewport`'s "Cube tiles" section for the model.
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
    ring_lo = int(math.floor(math.hypot(near_x, near_y) / edge_pc))
    ring_hi = int(math.floor(max(corners_r) / edge_pc))
    layer_lo = int(math.floor(lo[2] / edge_pc + 0.5))
    layer_hi = int(math.floor(hi[2] / edge_pc + 0.5))
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
    Whether the galaxy's bright-star scatter (`generate.py plan`) has run,
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
        "default_min_luminosity_sol": program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL,
    }


def _radius_sol(radius_km):
    """A star's radius in solar radii, from kilometres (`None` stays `None`)."""
    return None if radius_km is None else radius_km * 1000.0 / physical_constants.SOLAR_RADIUS_M


def _bright_star_entry(row):
    """One `bright_stars` row as the dict `galaxy_bright_stars_in_box` and
    `bright_stars_in_sector` return."""
    return {
        "id": row["id"],
        "x": row["position_x_mpc"] / MPC_PER_PC, "y": row["position_y_mpc"] / MPC_PER_PC,
        "z": row["position_z_mpc"] / MPC_PER_PC,
        "luminosity_sol": row["luminosity_w"] / physical_constants.SOLAR_LUMINOSITY,
        "temperature_k": row["temperature_k"], "radius_sol": _radius_sol(row["radius_km"]),
        "star_type": row["star_type"],
        "yerkes_class": row["yerkes_class"], "ring_index": row["ring_index"],
        "layer_index": row["layer_index"], "ring_slot_index": row["ring_slot_index"],
        "system_id": row["star_system_id"],
    }


def bright_stars_in_sector(conn, ring_index, layer_index, ring_slot_index, unfilled_only=True):
    """
    The pre-placed bright stars (`bright_stars`) in one sector cell, most
    luminous first -- the stars a sector page lists for a cell that hasn't
    been filled yet (filling builds each into a system, see
    `generate.py fill`). Reads `idx_bright_stars_address`, so it's cheap
    at any galaxy size.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
               star_type, yerkes_class, ring_index, layer_index, ring_slot_index, star_system_id
        FROM bright_stars
        WHERE ring_index = ? AND layer_index = ? AND ring_slot_index = ?{unfilled}
        ORDER BY luminosity_w DESC, id
        """,
        (int(ring_index), int(layer_index), int(ring_slot_index)),
    ).fetchall()
    return [_bright_star_entry(row) for row in rows]


def galaxy_brightest_stars(conn, count=GALAXY_TILE_BRIGHTEST_SAMPLE):
    """
    The galaxy's `count` most luminous bright stars, most luminous first
    (ties by descending id) -- `galaxy_tiles`' sample for tiles too big to
    query on their own. One walk down `idx_bright_stars_luminosity`.

    Returns:
        list[dict]: As `galaxy_bright_stars_in_box`.
    """
    rows = conn.execute(
        """
        SELECT id, position_x_mpc, position_y_mpc, position_z_mpc, luminosity_w, temperature_k, radius_km,
               star_type, yerkes_class, ring_index, layer_index, ring_slot_index, star_system_id
        FROM bright_stars FORCE INDEX (idx_bright_stars_luminosity)
        ORDER BY luminosity_w DESC, id DESC
        LIMIT ?
        """,
        (int(count),),
    ).fetchall()
    return [_bright_star_entry(row) for row in rows]


def galaxy_bright_stars_in_box(conn, lo, hi, edge_pc, limit=GALAXY_TILE_MAX_BRIGHT_STARS, unfilled_only=False,
                               brightest=None):
    """
    The most luminous pre-placed bright stars (`bright_stars`) in the box
    `[lo, hi)`, at most `limit` -- the stars the Galaxy Map draws before
    (and after) their sectors are filled.

    `bright_stars` is indexed by address and by luminosity, not by
    position, so the box is turned into address ranges: one per ring and
    layer it reaches (each with the ring's slots in its longitudes) for a
    small box, one per ring when that would be too many ranges. A box
    holding a good share of the galaxy instead walks the luminosity index
    from the top, which finds `limit` stars inside it quickly; the choice
    is by the estimated rows each way reads.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
                if len(found) >= limit:
                    break
        return found

    clauses, params = [], []
    # Both keys descending, so the luminosity path walks its index
    # backwards and stops at `limit`: with `id` ascending MariaDB and MySQL
    # sort every star in the box first, millions in a zoomed-out tile.
    order = "luminosity_w DESC, id DESC"
    if by_luminosity < by_address:
        where, index = box, "idx_bright_stars_luminosity"
    else:
        exact = sum(t - b + 1 for _r, b, t in bands) <= BRIGHT_STAR_MAX_EXACT_RANGES
        for ring, bottom, top in bands:
            slots = _ring_slot_ranges(ring, intervals)
            slot_test = " OR ".join("ring_slot_index BETWEEN ? AND ?" for _ in slots)
            slot_params = [v for pair in slots for v in pair]
            layers = [(layer, layer) for layer in range(bottom, top + 1)] if exact else [(bottom, top)]
            for first, last in layers:
                clauses.append(f"(ring_index = ? AND layer_index BETWEEN ? AND ? AND ({slot_test}))")
                params += [ring, first, last] + slot_params
        where, index = "(" + " OR ".join(clauses) + ") AND " + box, "idx_bright_stars_address"
    rows = conn.execute(
        f"""
        SELECT id, position_x_mpc, position_y_mpc, position_z_mpc, luminosity_w, temperature_k, radius_km,
               star_type, yerkes_class, ring_index, layer_index, ring_slot_index, star_system_id
        FROM bright_stars FORCE INDEX ({index})
        WHERE {where}
        ORDER BY {order}
        LIMIT ?
        """,
        params + box_params + [int(limit)],
    ).fetchall()
    return [_bright_star_entry(row) for row in rows]


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


def galaxy_generated_stars_in_box(conn, lo, hi, min_luminosity_sol, sector_count,
                                  limit=GALAXY_TILE_MAX_GENERATED_STARS):
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
        lo (tuple): `(x, y, z)` inclusive lower corner, parsecs.
        hi (tuple): `(x, y, z)` exclusive upper corner, parsecs.
        min_luminosity_sol (float): The faintest star listed.
        sector_count (int): How many generated (grid-placed) sectors the
            box holds, as its `galaxy_filled_in_box` summary counts them.
        limit (int): See `GALAXY_TILE_MAX_GENERATED_STARS`.

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
    rows = conn.execute(
        """
        SELECT st.id, st.name, st.luminosity_w, st.temperature_k, st.radius_km, st.star_type,
               ss.id AS system_id, ss.position_x_mpc, ss.position_y_mpc, ss.position_z_mpc,
               sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.ring_index, sec.layer_index, sec.ring_slot_index
        FROM sectors sec FORCE INDEX (idx_sectors_center)
        JOIN star_systems ss ON ss.sector_id = sec.id
        JOIN stars st ON st.star_system_id = ss.id
        LEFT JOIN bright_stars b ON b.star_system_id = ss.id AND st.role <> 'secondary'
        WHERE sec.center_x_pc >= ? AND sec.center_x_pc < ?
          AND sec.center_y_pc >= ? AND sec.center_y_pc < ?
          AND sec.center_z_pc >= ? AND sec.center_z_pc < ?
          AND sec.ring_index IS NOT NULL AND MOD(CRC32(sec.id), ?) = 0
          AND ss.position_x_mpc IS NOT NULL
          AND st.luminosity_w >= ? AND b.id IS NULL
        ORDER BY st.luminosity_w DESC, st.id DESC
        LIMIT ?
        """,
        (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], stride,
         float(min_luminosity_sol) * physical_constants.SOLAR_LUMINOSITY, int(limit)),
    ).fetchall()
    return [{
        "id": row["id"], "name": row["name"],
        "x": round(row["center_x_pc"] + row["position_x_mpc"] / MPC_PER_PC, 3),
        "y": round(row["center_y_pc"] + row["position_y_mpc"] / MPC_PER_PC, 3),
        "z": round(row["center_z_pc"] + row["position_z_mpc"] / MPC_PER_PC, 3),
        "luminosity_sol": float("%.4g" % (row["luminosity_w"] / physical_constants.SOLAR_LUMINOSITY)),
        "temperature_k": round(row["temperature_k"]),
        "radius_sol": float("%.3g" % _radius_sol(row["radius_km"])),
        "star_type": row["star_type"], "ring_index": row["ring_index"], "layer_index": row["layer_index"],
        "ring_slot_index": row["ring_slot_index"], "system_id": row["system_id"],
    } for row in rows]


def galaxy_tiles(conn, tile_keys):
    """
    The contents of each requested cube tile -- the interactive 3D Galaxy
    Map's data source. Every part of
    the result depends only on its tile key and the database's contents
    (see `galaxy_content_stamp`), so callers can cache each part by key.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        tile_keys (list[str]): `"level/ix/iy/iz"` keys (see
            `galaxyViewport.parse_tile_key`), at most
            `MAX_TILES_PER_REQUEST`.

    Returns:
        dict: `tiles` (`{key: {"placed": [...], "planned": [...],
            "filled": {...}, "clouds": [...], "stars": [...],
            "generated": [...]}}`, see `galaxy_sectors_in_box`,
            `galaxyViewport.planned_slots_in_tile`, `galaxy_filled_in_box`,
            `galaxy_clouds_in_box`, `galaxy_bright_stars_in_box` and
            `galaxy_generated_stars_in_box` (its floor from
            `generated_star_floor_sol`)),
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
    else:
        edge_pc = ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
        shape = None
        expected_system_count = None

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
            level, ix, iy, iz, edge_pc, shape, expected_system_count, exclude_addresses,
        )
        filled = galaxy_filled_in_box(conn, lo, hi, hi[0] - lo[0], edge_pc)
        clouds = galaxy_clouds_in_box(conn, lo, hi, margins=cloud_margins)
        stars = galaxy_bright_stars_in_box(conn, lo, hi, edge_pc, brightest=brightest)
        floor = generated_star_floor_sol(level)
        generated = []
        if floor is not None:
            sector_count = len(filled["cells"]) if filled["g"] == 1 else sum(cell[3] for cell in filled["cells"])
            generated = galaxy_generated_stars_in_box(conn, lo, hi, floor, sector_count)
        tiles[key] = {
            "placed": placed, "planned": planned, "filled": filled, "clouds": clouds, "stars": stars,
            "generated": generated,
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
        at (str or None): A block key, `m.ring.wedge.slab` (see
            `galaxyDrill.parse_drill_key`), or `None` for the galaxy.

    Returns:
        dict: `at` (the canonical key, or `None`), `child_m`, `children`
            (`[{ring, wedge, slab, generated}]`, at `child_m = 1` a
            sector's `wedge` is its slot and `slab` its layer), and
            `sectors` (`[{ring, layer, slot, id, name, system_count}]` at
            `child_m = 1`, else `None`).

    Raises:
        ValueError: On a malformed or impossible `at`.
    """
    if at is None or at == "":
        block = None
        child_m = DRILL_TOP
        rows = conn.execute(
            "SELECT ring_index, FLOOR((layer_index + ?) / ?) AS slab, ring_slot_index, COUNT(*) AS n "
            "FROM sectors WHERE ring_index IS NOT NULL GROUP BY ring_index, slab, ring_slot_index",
            ((DRILL_TOP - 1) // 2, DRILL_TOP),
        ).fetchall()
        counts = {}
        for r in rows:
            top = drill_chain_of(int(r["ring_index"]), int(r["slab"]) * DRILL_TOP, int(r["ring_slot_index"]))[0]
            counts[top] = counts.get(top, 0) + int(r["n"])
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
            f"SELECT id, name, ring_index, layer_index, ring_slot_index FROM sectors WHERE {where} "
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
        } for r in rows]
        counts = {DrillBlock(1, s["ring"], s["slot"], s["layer"]): 1 for s in sectors}
        return {"at": format_drill_key(block), "child_m": 1, "children": _stage_children(counts), "sectors": sectors}

    child_half = (child_m - 1) // 2
    rows = conn.execute(
        f"SELECT ring_index, FLOOR((layer_index + ?) / ?) AS slab, ring_slot_index, COUNT(*) AS n "
        f"FROM sectors WHERE {where} GROUP BY ring_index, slab, ring_slot_index",
        [child_half, child_m] + layer_params + params,
    ).fetchall()
    counts = {}
    for r in rows:
        child = drill_chain_of(int(r["ring_index"]), int(r["slab"]) * child_m, int(r["ring_slot_index"]))[level]
        counts[child] = counts.get(child, 0) + int(r["n"])
    return {"at": format_drill_key(block), "child_m": child_m, "children": _stage_children(counts), "sectors": None}


def _stage_children(counts):
    """`galaxy_stage`'s `children` list from `{DrillBlock: count}`."""
    return [
        {"ring": b.ring, "wedge": b.wedge, "slab": b.slab, "generated": n}
        for b, n in sorted(counts.items(), key=lambda item: (item[0].slab, item[0].ring, item[0].wedge))
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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

_STATE_TOKEN_RE = re.compile(r"^([0-9a-f]{16})\.(\d+)\.(\d+)\.(\d{17}|0)\.(\d+)$")


def _state_token(state):
    """`galaxy_content_state`'s dict as the opaque string the API hands
    out and `galaxy_changes` reads back."""
    return "{base}.{sectors}.{sector_max_id}.{sector_modified}.{system_max_id}".format(**state)


def _parse_state_token(token):
    """The dict `_state_token` encoded, or `None` for anything else."""
    match = _STATE_TOKEN_RE.match(str(token or ""))
    if not match:
        return None
    base, sectors, sector_max_id, sector_modified, system_max_id = match.groups()
    return {
        "base": base, "sectors": int(sectors), "sector_max_id": int(sector_max_id),
        "sector_modified": sector_modified, "system_max_id": int(system_max_id),
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
         "bright_stars": [bright["max_id"], bright["seed"]]}, sort_keys=True, default=str,
    ).encode("utf-8")).hexdigest()[:16]
    return {
        "base": base,
        "sectors": int(sector_row["n"]),
        "sector_max_id": int(maxima["sector_max_id"]),
        "sector_modified": _timestamp_digits(maxima["sector_modified"]),
        "system_max_id": int(maxima["system_max_id"]),
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
        since (str or None): A `state` from an earlier call (or
            `GET /api/galaxy/stamp`), or `None`.

    Returns:
        dict: `stamp` (`galaxy_content_stamp` now), `state` (the token to
            pass as `since` next time), `full` (every tile may have
            changed), `tiles` (sorted keys of the changed tiles; empty
            when `full`), and `stages` (sorted keys of the drill-down
            stages whose counts changed, `galaxy_stage_keys`; empty when
            `full`).
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
        result["full"] = True
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
    "class", "body", "life",
    "moon_class", "moon_body", "moon_life",
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
        conn (stellarObjects._db.Connection): An open connection.
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
    entry = program_constants.NEBULA_CLASSES.get(value)
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
    Class tags (both narrow: a class tag keeps only its own type's rows)."""
    classes = {}
    for tag in class_tags:
        type_label, _sep, value = tag.partition(":")
        if type_label in _SEARCH_PHENOMENON_CLASS_COLUMNS and value:
            classes.setdefault(type_label, set()).add(value)
    parts, params = [], []
    for table, type_label in _search_phenomenon_tables():
        if type_tags and type_label not in type_tags:
            continue
        if class_tags and type_label not in classes:
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


def _search_result_planets(conn, class_tags, body_tags, life_tags, term, limit, offset, size_range=None):
    clauses, params = [], []
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
        SELECT p.name, p.planet_class, p.body_type, p.life_chemical, p.radius_km,
               p.star_system_id, ss.name AS system_name, ss.sector_id
        """,
        f"FROM planets p JOIN star_systems ss ON ss.id = p.star_system_id WHERE 1=1{where}",
        "ss.name, ss.id, p.orbital_index, p.id",
        params, limit, offset,
    )
    return {"rows": [dict(r) for r in rows], **page}


def _search_result_moons(conn, class_tags, body_tags, life_tags, term, limit, offset, size_range=None):
    clauses, params = [], []
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
        SELECT m.name, m.planet_class, m.body_type, m.life_chemical, m.radius_km, p.name AS planet_name,
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
        ("moon_class", _search_facet_moon_class(conn)),
        ("moon_body", _search_facet_moon_body(conn)),
        ("moon_life", _search_facet_moon_life(conn)),
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


def search(conn, texts, tags, sizes=None, limit=SEARCH_RESULT_LIMIT, offsets=None):
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
        conn (stellarObjects._db.Connection): An open, read-only connection.
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
    density_tags = tags.get("density", set())
    type_tags = tags.get("type", set())

    facet_defs, autocomplete = _search_facets_cached(conn)
    facets = {name: options for name, options in facet_defs}
    facet_labels = {
        f"{name}:{opt['value']}": opt["label"]
        for name, options in facet_defs for opt in options
    }

    star_has_reason = bool(spectral_tags or luminosity_tags or texts.get("star_q") or star_size)
    planet_has_reason = bool(class_tags or body_tags or life_tags or texts.get("planet_q") or planet_size)
    moon_has_reason = bool(moon_class_tags or moon_body_tags or moon_life_tags or texts.get("moon_q") or moon_size)
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
    if texts.get("sector_q"):
        results["sectors"] = _search_result_sectors(conn, texts["sector_q"], *_page("sectors"))
    if texts.get("system_q"):
        results["systems"] = _search_result_systems(conn, texts["system_q"], *_page("systems"))
    if stars_included:
        results["stars"] = _search_result_stars(
            conn, spectral_tags, luminosity_tags, texts.get("star_q", ""), *_page("stars"), size_range=star_size
        )
    if planets_included:
        results["planets"] = _search_result_planets(
            conn, class_tags, body_tags, life_tags, texts.get("planet_q", ""), *_page("planets"),
            size_range=planet_size,
        )
    if moons_included:
        results["moons"] = _search_result_moons(
            conn, moon_class_tags, moon_body_tags, moon_life_tags, texts.get("moon_q", ""), *_page("moons"),
            size_range=moon_size,
        )
    if belts_included:
        results["belts"] = _search_result_belts(conn, density_tags, *_page("belts"))
    phenomenon_tags, phenomenon_class_tags = tags.get("phenomenon", set()), tags.get("phenomenon_class", set())
    if phenomenon_tags or phenomenon_class_tags:
        results["phenomena"] = _search_result_phenomena(
            conn, phenomenon_tags, phenomenon_class_tags, *_page("phenomena"))

    return {"facets": facets, "autocomplete": autocomplete, "facet_labels": facet_labels, "results": results}


def process_args():
    """
    Parses command-line arguments for the three subcommands: `sectors`,
    `systems`, and `near`.

    Returns:
        argparse.Namespace: The parsed arguments, including `command`
                            (which subcommand was invoked).
    """
    parser = argparse.ArgumentParser(
        description="List/query what's already stored in the planetGen database.",
    )
    parser.add_argument('--version', action=VersionAction, banner=version_banner('queryDb.py'))
    add_mysql_connection_args(parser)

    subparsers = parser.add_subparsers(dest='command', required=True)

    subparsers.add_parser('sectors', help="List every sector, with its size and system count.")

    systems_parser = subparsers.add_parser('systems', help="List systems, optionally filtered.")
    systems_parser.add_argument('--star-type', type=str,
                                help="Only systems with a star whose type starts with this "
                                     "(e.g. 'G' for every G-type system, 'G2V' for an exact match).")
    systems_parser.add_argument('--sector-id', type=int, help="Only systems in this sector.")

    near_parser = subparsers.add_parser(
        'near', help="Find systems within a radius of another system, in the same sector.",
    )
    near_parser.add_argument('system_id', type=int, help="The star_systems.id to measure distances from.")
    near_parser.add_argument('--radius', type=float, required=True,
                             help="Search radius in light-years (e.g. 50 for 'everything within 50 ly').")

    planets_parser = subparsers.add_parser(
        'planets', help="List planets, optionally filtered by class, radius, sector, or system.",
    )
    planets_parser.add_argument('--class', dest='planet_class', type=str,
                                help="Only planets of this exact class (e.g. 'M').")
    planets_parser.add_argument('--min-radius-km', type=float, help="Only planets at least this large.")
    planets_parser.add_argument('--max-radius-km', type=float, help="Only planets at most this large.")
    planets_parser.add_argument('--sector-id', type=int, help="Only planets whose system is in this sector.")
    planets_parser.add_argument('--system-id', type=int, help="Only planets in this one system.")

    moons_parser = subparsers.add_parser(
        'moons', help="List moons, optionally filtered by class, radius, sector, or system.",
    )
    moons_parser.add_argument('--class', dest='planet_class', type=str,
                              help="Only moons of this exact class (e.g. 'M').")
    moons_parser.add_argument('--min-radius-km', type=float, help="Only moons at least this large.")
    moons_parser.add_argument('--max-radius-km', type=float, help="Only moons at most this large.")
    moons_parser.add_argument('--sector-id', type=int, help="Only moons whose system is in this sector.")
    moons_parser.add_argument('--system-id', type=int, help="Only moons in this one system.")

    return parser.parse_args()


def main():
    """
    The main entry point: dispatches to the requested subcommand and prints
    a plain-text listing of the results.
    """
    args = process_args()
    conn = open_readonly(mysql_config_from_args(args))
    try:
        if args.command == 'sectors':
            sectors = list_sectors(conn)
            if not sectors:
                print("No sectors stored.")
                return
            for sector in sectors:
                print(f"[{sector['id']}] {sector['name']} "
                      f"(edge {sector['edge_ly']:.2f} ly, {sector['system_count']} systems)")

        elif args.command == 'systems':
            systems = list_systems(conn, star_type_prefix=args.star_type, sector_id=args.sector_id)
            if not systems:
                print("No matching systems.")
                return
            for system in systems:
                kind = "binary" if system["is_binary"] else "single"
                sector_note = f"sector {system['sector_id']}" if system["sector_id"] is not None else "standalone"
                print(f"[{system['id']}] {system['name']} ({kind}, {sector_note})")

        elif args.command == 'near':
            matches = systems_within_radius(conn, args.system_id, args.radius)
            if not matches:
                print(f"No other systems within {args.radius} ly.")
                return
            for match in matches:
                print(f"[{match['id']}] {match['name']} -- {match['distance_ly']:.2f} ly")

        elif args.command == 'planets':
            planets = list_planets(
                conn, planet_class=args.planet_class, min_radius_km=args.min_radius_km,
                max_radius_km=args.max_radius_km, sector_id=args.sector_id, system_id=args.system_id,
            )
            if not planets:
                print("No matching planets.")
                return
            for planet in planets:
                print(f"[{planet['id']}] {planet['name']} (Class {planet['planet_class']}, "
                      f"{planet['radius_km']:.0f} km) -- {planet['system_name']}")

        elif args.command == 'moons':
            moons = list_moons(
                conn, planet_class=args.planet_class, min_radius_km=args.min_radius_km,
                max_radius_km=args.max_radius_km, sector_id=args.sector_id, system_id=args.system_id,
            )
            if not moons:
                print("No matching moons.")
                return
            for moon in moons:
                print(f"[{moon['id']}] {moon['name']} (Class {moon['planet_class']}, "
                      f"{moon['radius_km']:.0f} km) -- {moon['system_name']} / {moon['planet_name']}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
