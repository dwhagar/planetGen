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
with `SELECT`-only grants for this tool (and the read-only Flask API,
`html/api/config.py`) rather than the read-write account `sectorGen.py`/
`systemGen.py` use -- MySQL has no per-connection "open this read-only"
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
import hashlib
import json
import math

import pymysql

from stellarObjects._db import add_mysql_connection_args, get_connection, get_galaxy_shape, mysql_config_from_args
from stellarObjects._version import VersionAction, __version__, version_banner
from stellarObjects.galaxyGeometry import provisional_sector_designation, sector_position_pc
from stellarObjects.galaxyViewport import (
    PLANNED_RADIUS_CAP_PC,
    density_points_for_tile,
    density_sample_points,
    parse_tile_key,
    planned_slots_in_tile,
    planned_slots_in_view,
    tile_bounds_pc,
)
from stellarObjects.navGraph import build_knn_adjacency, shortest_path
from stellarObjects.navigation import course_between, warp_travel_times
from stellarObjects.physical_constants import SPECTRAL_CLASS_COLORS
from stellarObjects.sectorGeometry import lateral_neighbor_slots, radial_neighbor_slot
from stellarObjects.program_constants import DEFAULT_SECTOR_EDGE_LY, NAV_ADJACENCY_K, PLANET_CLASSES
from stellarObjects.utils import ly_to_pc, milliparsecs_to_ly, mpc_to_pc, pc_to_ly


def open_readonly(config=None):
    """
    Opens a connection for this read-only tool -- see the module
    docstring for why "read-only" is enforced by the configured account's
    grants rather than anything this function does itself.

    Args:
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.

    Returns:
        stellarObjects._db.Connection: An open connection.

    Raises:
        SystemExit: If the database can't be reached.
    """
    try:
        return get_connection(config, ensure_schema=False)
    except pymysql.MySQLError as exc:
        raise SystemExit(f"Error: could not open the database ({exc}).")


def list_sectors(conn, limit=None, offset=None):
    """
    Returns every sector, with its edge length (converted to light-years)
    and how many systems it contains.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        limit (int, optional): Caps the number of rows returned. `None`
            (the default -- what every existing caller of this function
            still gets) returns every sector.
        offset (int, optional): Skips this many rows first. Ignored
            unless `limit` is also given; meaningless on its own.

    Returns:
        list[dict]: One row per sector, with `id`, `name`,
                           `edge_ly`, `system_count`.
    """
    query = """
        SELECT sec.id, sec.name, sec.edge_mpc, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.shell_index, sec.shell_slot_index, COUNT(ss.id) AS system_count
        FROM sectors sec
        LEFT JOIN star_systems ss ON ss.sector_id = sec.id
        GROUP BY sec.id
        ORDER BY sec.name
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
            "shell_index": r["shell_index"], "shell_slot_index": r["shell_slot_index"],
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
        conditions.append("s.star_type LIKE ?")
        params.append(f"{star_type_prefix}%")

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
                           `sector_id`, `is_binary`, and `star_summary`
                           (the single star's `star_type`, or a binary's
                           `binary_type` -- what `html/browse.py`/
                           `html/search.py` show as a system's "Star type"
                           column).
    """
    join_sql, where_sql, params = _systems_filter_clause(star_type_prefix, sector_id)
    query = f"""
        SELECT DISTINCT ss.id, ss.name, ss.sector_id, ss.is_binary, ss.binary_configuration, ss.binary_type,
               (SELECT s.star_type FROM stars s WHERE s.star_system_id = ss.id AND s.role = 'single' LIMIT 1)
                   AS single_star_type,
               (SELECT s.star_type FROM stars s WHERE s.star_system_id = ss.id AND s.role = 'primary' LIMIT 1)
                   AS primary_star_type,
               (SELECT s.star_type FROM stars s WHERE s.star_system_id = ss.id AND s.role = 'secondary' LIMIT 1)
                   AS secondary_star_type
        FROM star_systems ss{join_sql}{where_sql} ORDER BY ss.name
        """

    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params = params + [limit, offset or 0]

    rows = conn.execute(query, params).fetchall()
    return [
        {
            "id": r["id"], "name": r["name"], "sector_id": r["sector_id"], "is_binary": r["is_binary"],
            "star_summary": _star_summary(r),
        }
        for r in rows
    ]


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
        phenomenon_type (str): One of `_PHENOMENON_TYPE_TO_TABLE`'s keys
            (`"nebula"`, `"asteroid_field"`, `"black_hole"`, `"neutron_star"`).
        phenomenon_id (int): The phenomenon's own row id.

    Returns:
        dict: `sector_id` (always `None`), `position_ly` (always `None`),
            `galaxy_position_ly` (`(x, y, z)` tuple in light-years, or
            `None` if this phenomenon has never been placed in the galaxy).

    Raises:
        ValueError: If `phenomenon_type` is unrecognized, is a type with no
            galaxy-frame placement columns at all (`_UNPLACED_PHENOMENON_
            TABLES`'s own three types -- see `_SUPERNOVA_REMNANT_TABLE`'s
            own docstring), or no such row exists.
    """
    if phenomenon_type not in _PLACEABLE_PHENOMENON_TYPES:
        raise ValueError(
            f"{phenomenon_type!r} has no galaxy-frame placement and can never be a NAV endpoint"
        )
    table = _PHENOMENON_TYPE_TO_TABLE[phenomenon_type]

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
    hole/neutron star) -- a direct course (distance/azimuth/altitude/warp
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
            `navigation.WarpLeg`, for `direct.distance_ly`),
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

    direct = course_between(origin_position, destination_position)

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
        "origin_position": origin_position,
        "destination_position": destination_position,
        "route": route,
    }


def sector_neighbors(conn, sector):
    """
    Every immediately-surrounding sector address for a galaxy-placed
    sector -- its exact same-shell (lateral) Voronoi neighbors
    (`sectorGeometry.lateral_neighbor_slots`) plus its nearest inward and
    outward radial neighbor (`sectorGeometry.radial_neighbor_slot`), each
    tagged with whether a real `sectors` row already exists there. Drives
    the Sector Map's (`html/lib/starmap.py`) neighboring-sector
    indicators -- an existing neighbor's indicator links straight to it
    (`sector.py`); a not-yet-generated one shows its address so it can be
    fed to `generate.py galaxy --shell K --slot N`, the same convention
    the Galaxy Map's own "planned" tier already uses
    (`galaxyViewport.planned_slots_in_view`).

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        sector (dict or Row): This sector's own row -- needs
            `shell_index`, `shell_slot_index`, and `edge_mpc`.

    Returns:
        list[dict]: `shell_index`, `shell_slot_index`, `direction_pc`
            (`[x, y, z]`, this neighbor's galaxy-frame center minus this
            sector's own -- a direction, not clamped to the sector's own
            cube), `exists` (bool), `sector_id`/`sector_name` (both `None`
            when `exists` is `False`), and `designation`
            (`provisional_sector_designation`, always present so a
            not-yet-generated neighbor has a human-readable label even
            with no name of its own). `[]` for a sector with no galaxy
            placement.
    """
    shell_index, shell_slot_index, edge_mpc = sector["shell_index"], sector["shell_slot_index"], sector["edge_mpc"]
    if shell_index is None or shell_slot_index is None or not edge_mpc:
        return []

    edge_pc = mpc_to_pc(edge_mpc)
    edge_ly = pc_to_ly(edge_pc)
    this_position = sector_position_pc(shell_index, shell_slot_index, edge_pc)

    addresses = [(shell_index, slot) for slot in lateral_neighbor_slots(shell_index, shell_slot_index, edge_pc)]
    for direction in (-1, 1):
        radial_slot = radial_neighbor_slot(shell_index, shell_slot_index, edge_pc, direction)
        if radial_slot is not None:
            addresses.append((shell_index + direction, radial_slot))
    if not addresses:
        return []

    clauses = " OR ".join(["(shell_index = ? AND shell_slot_index = ?)"] * len(addresses))
    params = [value for address in addresses for value in address]
    rows = conn.execute(
        f"SELECT id, name, shell_index, shell_slot_index FROM sectors WHERE {clauses}", tuple(params),
    ).fetchall()
    existing = {(r["shell_index"], r["shell_slot_index"]): r for r in rows}

    results = []
    for addr_shell, addr_slot in addresses:
        position = sector_position_pc(addr_shell, addr_slot, edge_pc)
        direction_pc = [
            position[0] - this_position[0], position[1] - this_position[1], position[2] - this_position[2],
        ]
        match = existing.get((addr_shell, addr_slot))
        results.append({
            "shell_index": addr_shell, "shell_slot_index": addr_slot,
            "direction_pc": direction_pc,
            "exists": match is not None,
            "sector_id": match["id"] if match else None,
            "sector_name": match["name"] if match else None,
            "designation": provisional_sector_designation(addr_shell, addr_slot, edge_pc, edge_ly),
        })
    return results


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
            `center_y_pc`/`center_z_pc`, `shell_index`, `shell_slot_index`,
            `placed`, `system_count`, and `systems` (one entry per system
            placed in this sector: `id`, `name`, `quadrant`, `location`,
            `is_binary`, `binary_type`, `position_x_mpc`/`position_y_mpc`/
            `position_z_mpc`, and `stars` -- 1 entry (single) or 2
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
               position_x_mpc, position_y_mpc, position_z_mpc
        FROM star_systems
        WHERE sector_id = ?
        ORDER BY name
        """,
        (sector_id,),
    ).fetchall()

    systems = []
    for row in system_rows:
        star_rows = conn.execute(
            "SELECT role, star_type, temperature_k, radius_km, luminosity_w"
            " FROM stars WHERE star_system_id = ?"
            " ORDER BY CASE role WHEN 'secondary' THEN 1 ELSE 0 END",
            (row["id"],),
        ).fetchall()
        systems.append({
            "id": row["id"], "name": row["name"], "quadrant": row["quadrant"], "location": row["location"],
            "is_binary": row["is_binary"], "binary_type": row["binary_type"],
            "position_x_mpc": row["position_x_mpc"], "position_y_mpc": row["position_y_mpc"],
            "position_z_mpc": row["position_z_mpc"],
            "stars": [dict(star_row) for star_row in star_rows],
        })

    return {
        "id": sector["id"], "name": sector["name"], "edge_mpc": sector["edge_mpc"],
        "edge_ly": milliparsecs_to_ly(sector["edge_mpc"]),
        "center_x_pc": sector["center_x_pc"], "center_y_pc": sector["center_y_pc"],
        "center_z_pc": sector["center_z_pc"],
        "shell_index": sector["shell_index"], "shell_slot_index": sector["shell_slot_index"],
        "placed": sector["center_x_pc"] is not None,
        "system_count": len(systems),
        "systems": systems,
        "phenomena": phenomena_near_sector(conn, sector_id),
        "neighbors": sector_neighbors(conn, sector),
        "wiki_url": sector["wiki_url"],
    }


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
)
"""tuple: `(table_name, type_label, descriptor_expr, radius_expr)` for each
galaxy-placeable standalone phenomenon `phenomena_near_sector`/
`galaxy_placed_phenomena` read from -- see `schema.sql`'s "v18"/"v21"
header notes. `descriptor_expr` is a SQL expression for each table's own
one-line flavor field (a nebula's `nebula_type`, a black hole's accretion
state), normalized to a common `descriptor` key so callers don't need to
know which table a given `type` came from; `radius_expr` is likewise a SQL
expression (a plain column for nebulae/asteroid fields, a literal `0` for
the two point-like compact-remnant types)."""


_SUPERNOVA_REMNANT_TABLE = ("supernova_remnants", "supernova_remnant", "morphology", "radius_ly")
"""tuple: The `supernova_remnant` counterpart to one `_PHENOMENON_TABLES`
entry -- kept OUT of that tuple deliberately, since `supernova_remnants`
never gained the v21 galaxy-frame placement columns
(`center_x_pc`/`center_y_pc`/`center_z_pc`/`galactic_radius_pc`) the other
four phenomenon tables did (see `schema.sql`'s own "v16"/"v21" header
notes) -- so it can never appear on the Galaxy Map or as a NAV endpoint,
and can't share `_placed_phenomenon_rows`/`galaxy_placed_phenomena`/
`phenomena_near_sector`'s common query shape, which all select
`center_x_pc`. It still has its own real `radius_ly` (unlike the two
point-like compact-remnant types), so it's fully usable everywhere that
only needs the phenomenon's own physical size -- `list_phenomena`/
`count_phenomena`'s flat listing, `phenomenon_detail`'s page, and
`lib/phenomenonmap.py`'s AU-scale diagram."""

_ROGUE_PLANET_TABLE = (
    "rogue_planets", "rogue_planet",
    "(CASE WHEN planet_type = 'g' THEN 'gas giant' ELSE 'terrestrial' END)", "0",
)
"""tuple: The `rogue_planet` counterpart to one `_PHENOMENON_TABLES` entry
-- same "no galaxy-frame placement columns at all" reasoning as
`_SUPERNOVA_REMNANT_TABLE` (see that constant's own docstring; `schema.sql`
never gave `rogue_planets` `center_x/y/z_pc` either). Its `radius_expr` is
a literal `0`, not a real column, for the same reason `black_holes`/
`neutron_stars` use one in `_PHENOMENON_TABLES`: a rogue planet's own
`radius_km` is planet-scale, utterly negligible next to the light-year
scale `list_phenomena`'s shared `radius_ly` column otherwise means.

Confirmed missing end-to-end before this was added: `generate.py`
(`generate_sector_phenomena`) has always generated and saved these at a
non-trivial rate (`program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM`'s
own `"rogue-planet": 0.1` -- roughly one per ten star systems, far more
common than a nebula), but no query function anywhere ever read the
`rogue_planets` table, so a generated rogue planet was completely
invisible in every listing/page despite existing in the database the
whole time."""

_INTERSTELLAR_COMET_TABLE = (
    "interstellar_comets", "interstellar_comet",
    "(CASE WHEN is_active THEN 'active' ELSE 'dormant' END)", "0",
)
"""tuple: The `interstellar_comet` counterpart to one `_PHENOMENON_TABLES`
entry -- same "no galaxy-frame placement columns at all" reasoning as
`_ROGUE_PLANET_TABLE` immediately above (and the same "confirmed missing
end-to-end" history: generated via `program_constants.
PHENOMENON_RATE_PER_STAR_SYSTEM`'s `"comet"` rate, never once queried).
`radius_expr` is a literal `0` for the same reason -- `nucleus_diameter_km`
is negligible at this shared column's light-year scale. Named
`interstellar_comet`, not bare `comet`, to stay unambiguous next to the
unrelated `comets` table (a star system's own planet-orbiting comets,
`queryDb.system_detail`'s own `comets` key -- a completely different
table this constant has nothing to do with)."""

_UNPLACED_PHENOMENON_TABLES = (_SUPERNOVA_REMNANT_TABLE, _ROGUE_PLANET_TABLE, _INTERSTELLAR_COMET_TABLE)
"""tuple: Every phenomenon type with no galaxy-frame placement columns at
all (as opposed to `_PHENOMENON_TABLES`' own four, which simply may or may
not be placed yet) -- `_PHENOMENON_TABLES + _UNPLACED_PHENOMENON_TABLES`
is `list_phenomena`/`count_phenomena`/`_PHENOMENON_TYPE_TO_TABLE`'s own
"every type" tuple, shared here so a future phenomenon type only ever
needs adding in one place."""

_UNPLACED_PHENOMENON_TABLE_NAMES = frozenset(table for table, *_rest in _UNPLACED_PHENOMENON_TABLES)
"""frozenset: Just the table names out of `_UNPLACED_PHENOMENON_TABLES` --
`list_phenomena`'s own union query selects a literal `NULL AS center_x_pc`
for any of these (they have no such column at all to select), unlike
`_PHENOMENON_TABLES`' own four, which select their real column."""


def _placed_phenomenon_rows(conn, bbox=None):
    """
    Reads every galaxy-placed row (non-NULL `center_x_pc`) from all four
    v18/v21-placeable standalone-phenomenon tables, normalized to one
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

    Returns:
        list[dict]: `id`, `type` (`"nebula"`, `"asteroid_field"`,
            `"black_hole"`, or `"neutron_star"`), `name`, `descriptor`,
            `radius_ly` (0 for the two compact-remnant types), `x`/`y`/`z`
            (`center_x/y/z_pc`), `galactic_radius_pc`.
    """
    where_bbox = ""
    bbox_params = ()
    if bbox is not None:
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
                "x": row["center_x_pc"], "y": row["center_y_pc"], "z": row["center_z_pc"],
                "galactic_radius_pc": row["galactic_radius_pc"],
            })
    return rows


def _widest_placed_phenomenon_radius_ly(conn):
    """
    The largest `radius_ly` among currently placed `nebulae`/
    `asteroid_fields` rows (the only two phenomenon tables with a real,
    non-zero radius -- see `_PHENOMENON_TABLES`), or `0.0` if neither
    table has any placed row at all.

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
    for table in ("nebulae", "asteroid_fields"):
        row = conn.execute(
            f"SELECT MAX(radius_ly) AS widest FROM {table} WHERE center_x_pc IS NOT NULL"
        ).fetchone()
        if row and row["widest"] is not None:
            widest = max(widest, row["widest"])
    return widest


def list_phenomena(conn, limit=None, offset=None):
    """
    Returns every exotic phenomenon (nebula/asteroid field/black hole/
    neutron star/supernova remnant/rogue planet/interstellar comet -- the
    four in `_PHENOMENON_TABLES` plus `_UNPLACED_PHENOMENON_TABLES`),
    across every sector and regardless of galaxy placement -- `GET
    /api/phenomena`'s own flat listing (`html/phenomena.py`), unlike
    `galaxy_placed_phenomena` (which only returns the galaxy-placed
    subset, for the Galaxy Map) or `phenomena_near_sector` (one sector's
    own neighborhood) -- neither of which a supernova remnant/rogue
    planet/interstellar comet can ever appear in, since none of those
    three tables have galaxy-frame placement columns at all (see
    `_SUPERNOVA_REMNANT_TABLE`'s own docstring); their `placed` is
    therefore always `False` here.

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
            position at all, `center_x_pc IS NOT NULL`; always `False` for
            any of `_UNPLACED_PHENOMENON_TABLES`'s own three types).
    """
    union_parts = [
        f"""
        SELECT '{type_label}' AS type, t.id AS id, t.name AS name,
               {descriptor_expr} AS descriptor, {radius_expr} AS radius_ly,
               t.sector_id AS sector_id, sec.name AS sector_name,
               {"NULL" if table in _UNPLACED_PHENOMENON_TABLE_NAMES else "t.center_x_pc"} AS center_x_pc
        FROM {table} t
        LEFT JOIN sectors sec ON sec.id = t.sector_id
        {"WHERE t.star_id IS NULL" if table in ("black_holes", "neutron_stars") else ""}
        """
        for table, type_label, descriptor_expr, radius_expr in _PHENOMENON_TABLES + _UNPLACED_PHENOMENON_TABLES
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
    `_PHENOMENON_TABLES` plus `_UNPLACED_PHENOMENON_TABLES`, ignoring any
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
        for table, _type_label, _descriptor_expr, _radius_expr in _PHENOMENON_TABLES + _UNPLACED_PHENOMENON_TABLES
    )


_PHENOMENON_TYPE_TO_TABLE = {
    type_label: table for table, type_label, _de, _re in _PHENOMENON_TABLES + _UNPLACED_PHENOMENON_TABLES
}
"""dict: `type` value (as returned by `list_phenomena`/`galaxy_placed_phenomena`)
-> its backing table name, e.g. `"nebula"` -> `"nebulae"` -- the reverse of
`_PHENOMENON_TABLES`'s own `(table, type_label, ...)` order, used by
`phenomenon_detail` to find the one table a `(type, id)` pair actually
means without hand-listing the mapping a second time. Includes
`_UNPLACED_PHENOMENON_TABLES` too, since `phenomenon_detail`'s plain
`SELECT t.*` works for those exactly the same as for the other four types
even though they can't participate in the galaxy-placement-only helpers
below."""

_PLACEABLE_PHENOMENON_TYPES = frozenset(type_label for _t, type_label, _de, _re in _PHENOMENON_TABLES)
"""frozenset: The `type` values that DO have galaxy-frame placement
columns (everything in `_PHENOMENON_TABLES`) -- deliberately excludes
every one of `_UNPLACED_PHENOMENON_TABLES`'s own types, unlike
`_PHENOMENON_TYPE_TO_TABLE` above. Guards `_load_nav_phenomenon_endpoint`
against ever running its `center_x_pc` SELECT against one of those tables,
none of which has any such column, and would otherwise raise a raw SQL
error instead of the clean `ValueError` a NAV
request for an inherently unplaceable phenomenon type should get."""


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
        dict: Every column of the phenomenon's own row, plus `type` and
            `sector_name` (`None` if it has no `sector_id`).

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
    return detail


def galaxy_placed_phenomena(conn):
    """
    Every galaxy-placed nebula/asteroid field/black hole/neutron star --
    the phenomenon counterpart to `galaxy_placed_sectors`, plotted as small
    dots on the same Galaxy Map (`html/lib/galaxymap.py`).

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.

    Returns:
        list[dict]: See `_placed_phenomenon_rows`.
    """
    return _placed_phenomenon_rows(conn)


def phenomena_near_sector(conn, sector_id):
    """
    Every galaxy-placed nebula/asteroid field/black hole/neutron star
    whose sphere could plausibly reach into `sector_id`'s own cube -- the
    data `html/lib/starmap.py`'s Sector Map draws (translucent clouds for
    nebulae/asteroid fields, point markers for the two compact-remnant
    types, whose own `radius_ly` is always 0 -- see `_PHENOMENON_TABLES`).

    An exact cube-vs-sphere overlap test isn't worth the complexity here,
    so this compares against each cube's own *bounding* sphere (radius =
    half its space diagonal, `edge_pc * sqrt(3) / 2`) instead: a safe,
    exact upper bound that can only ever include a few extra phenomena
    whose sphere clips the bounding sphere but not the cube itself (out
    near a corner), never silently miss a real overlap -- consistent with
    this project's existing "sector cubes already accept small real-world
    gaps/overlaps" tolerance (`docs/design/galaxy-coordinate-system.md`
    section 3) rather than a false-negative risk.

    Reads via `_placed_phenomenon_rows`'s own `bbox` argument -- a SQL
    bounding-box prefilter, not the whole-galaxy scan that function's
    default (`bbox=None`) runs -- sized to `half_diagonal_pc` plus
    whatever the widest currently-placed `radius_ly` actually is
    (`_widest_placed_phenomenon_radius_ly`), so it stays exactly as
    correct as scanning every row (a phenomenon of any size, however
    large, that could plausibly overlap is still included) while letting
    MySQL range-scan `idx_{table}_center` instead of examining every row
    in all four tables on every call -- see `schema.sql`'s "v26" header
    note for the production failure this fixes (a genuine full-table-scan
    -times-four on every `GET /api/sectors/<id>`, "the exact same failure
    mode `schema.sql`'s "v25" note already documented for `sectors`).

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        sector_id (int): The `sectors.id` to check against.

    Returns:
        list[dict]: One entry per candidate phenomenon: `id`, `type`
            (`"nebula"`, `"asteroid_field"`, `"black_hole"`, or
            `"neutron_star"`), `name`, `descriptor`,
            `radius_ly`, `distance_ly` (sector center to phenomenon
            center), and `offset_x_ly`/`offset_y_ly`/`offset_z_ly` (the
            phenomenon's center relative to the sector's own center, in
            light-years -- the same frame `starmap.py` already places
            stars in). Empty if this sector has no galaxy placement of its
            own.

    Raises:
        ValueError: If no such sector exists.
    """
    sector = conn.execute(
        "SELECT center_x_pc, center_y_pc, center_z_pc, edge_mpc FROM sectors WHERE id = ?",
        (sector_id,),
    ).fetchone()
    if sector is None:
        raise ValueError(f"no sectors row with id {sector_id}")
    if sector["center_x_pc"] is None:
        return []

    half_diagonal_pc = mpc_to_pc(sector["edge_mpc"]) * math.sqrt(3) / 2
    margin_pc = half_diagonal_pc + ly_to_pc(_widest_placed_phenomenon_radius_ly(conn))
    bbox = (sector["center_x_pc"], sector["center_y_pc"], sector["center_z_pc"], margin_pc)

    matches = []
    for phenomenon in _placed_phenomenon_rows(conn, bbox=bbox):
        dx = phenomenon["x"] - sector["center_x_pc"]
        dy = phenomenon["y"] - sector["center_y_pc"]
        dz = phenomenon["z"] - sector["center_z_pc"]
        distance_pc = math.sqrt(dx * dx + dy * dy + dz * dz)
        if distance_pc > half_diagonal_pc + ly_to_pc(phenomenon["radius_ly"]):
            continue
        matches.append({
            "id": phenomenon["id"], "type": phenomenon["type"], "name": phenomenon["name"],
            "descriptor": phenomenon["descriptor"], "radius_ly": phenomenon["radius_ly"],
            "distance_ly": pc_to_ly(distance_pc),
            "offset_x_ly": pc_to_ly(dx), "offset_y_ly": pc_to_ly(dy), "offset_z_ly": pc_to_ly(dz),
        })
    return matches


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
            `binary_mutual_position_x/y/z_km` (the secondary's position
            relative to the primary -- NULL for a single star; see
            `schema.sql`'s "v14"/"v15" notes), `markdown_content`,
            `wikitext_content`, `wikijs_url`/`mediawiki_url`
            (`star_systems.wikijs_url`/`mediawiki_url` -- `None` for
            whichever wiki (or both) this system hasn't been uploaded to
            yet; see `schema.sql`'s "v22" header note), `stars` (id/role/name/
            star_type/mass_kg/radius_km/temperature_k/luminosity_w --
            `id` matches a `'wide'` binary's `planets`/`belts` rows' own
            `star_id`, disambiguating which star each orbits), `planets`
            (each a `planets` row, including its own `star_id`, plus its
            own `moons` list), `belts` (`asteroid_belts` rows, including
            `star_id`), `comets` (`comets` rows, including `star_id` --
            no `orbital_index`, see `insert_comet`'s docstring), and
            `sector_siblings` (`{id, name}` for every other system in the
            same sector, empty if standalone -- for linkifying
            `location`'s "nearest: ..." names without a second round
            trip).

    Raises:
        ValueError: If no such system exists.
    """
    system = conn.execute("SELECT * FROM star_systems WHERE id = ?", (system_id,)).fetchone()
    if system is None:
        raise ValueError(f"no star_systems row with id {system_id}")

    stars = conn.execute(
        "SELECT id, role, name, star_type, mass_kg, radius_km, temperature_k, luminosity_w"
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
    planets = []
    for planet in planet_rows:
        moon_rows = conn.execute(
            "SELECT * FROM moons WHERE planet_id = ? ORDER BY orbital_index", (planet["id"],)
        ).fetchall()
        planet_dict = dict(planet)
        planet_dict["moons"] = [dict(m) for m in moon_rows]
        planets.append(planet_dict)

    belts = conn.execute(
        "SELECT * FROM asteroid_belts WHERE star_system_id = ? ORDER BY orbital_index", (system_id,)
    ).fetchall()

    comets = conn.execute(
        "SELECT * FROM comets WHERE star_system_id = ? ORDER BY id", (system_id,)
    ).fetchall()

    sector_siblings = []
    if system["sector_id"] is not None:
        sibling_rows = conn.execute(
            "SELECT id, name FROM star_systems WHERE sector_id = ?", (system["sector_id"],)
        ).fetchall()
        sector_siblings = [{"id": r["id"], "name": r["name"]} for r in sibling_rows]

    return {
        "id": system["id"], "name": system["name"], "sector_id": system["sector_id"],
        "quadrant": system["quadrant"], "location": system["location"],
        "is_binary": system["is_binary"], "binary_type": system["binary_type"],
        "binary_configuration": system["binary_configuration"],
        "binary_mutual_position_x_km": system["binary_mutual_position_x_km"],
        "binary_mutual_position_y_km": system["binary_mutual_position_y_km"],
        "binary_mutual_position_z_km": system["binary_mutual_position_z_km"],
        "markdown_content": system["markdown_content"], "wikitext_content": system["wikitext_content"],
        "wikijs_url": system["wikijs_url"], "mediawiki_url": system["mediawiki_url"],
        "stars": [dict(s) for s in stars],
        "planets": planets,
        "belts": [dict(b) for b in belts],
        "comets": [dict(c) for c in comets],
        "sector_siblings": sector_siblings,
    }


def galaxy_placed_sectors(conn):
    """
    Every sector with a galaxy position, plus its live system count -- the
    data `html/galaxy.py`'s Galaxy Map (`html/lib/galaxymap.py`) plots.
    Unplaced sectors (`center_x_pc IS NULL`) have nothing to plot and are
    excluded at the query itself.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.

    Returns:
        list[dict]: `id`, `name`, `x`/`y`/`z` (`center_x/y/z_pc`),
            `galactic_radius_pc`, `shell_index`, `system_count`.
    """
    rows = conn.execute(
        """
        SELECT sec.id, sec.name, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.galactic_radius_pc, sec.shell_index,
               (SELECT COUNT(*) FROM star_systems ss WHERE ss.sector_id = sec.id) AS system_count
        FROM sectors sec
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
            `outer_shell_index`, `expected_system_count_at_density_1`.
            `None` if `generate.py plan` has never been run against this
            database (no `galaxy_shape` row yet).
    """
    skeleton = get_galaxy_shape(conn)
    if skeleton is None:
        return None
    shape = dict(skeleton.shape._asdict())
    shape["edge_pc"] = skeleton.edge_pc
    shape["outer_shell_index"] = skeleton.outer_shell_index
    shape["expected_system_count_at_density_1"] = skeleton.expected_system_count_at_density_1
    return shape


# ---------------------------------------------------------------------
# Interactive 3D Galaxy Map viewport queries -- backs GET /api/galaxy/view
# and `html/galaxy_view.py`, called fresh every time the map's own camera
# moves. Unlike `galaxy_placed_sectors`/`galaxy_density_shape` above (each
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
            `shell_index`, `shell_slot_index`, `designation`
            (`provisional_sector_designation`, `None` if this sector has
            no `shell_slot_index` -- a placement predating the v8 schema's
            per-slot addressing), `system_count`, `edge_ly` (this
            sector's own real edge length -- lets a client compute its
            true stellar density, `system_count / edge_ly ** 3`, rather
            than just its raw system count), `distance_pc` (from the
            given center).
    """
    rows = conn.execute(
        """
        SELECT sec.id, sec.name, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.galactic_radius_pc, sec.shell_index, sec.shell_slot_index, sec.edge_mpc,
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
        shell_index, shell_slot_index = r["shell_index"], r["shell_slot_index"]
        if shell_index is not None and shell_slot_index is not None and r["edge_mpc"]:
            designation = provisional_sector_designation(
                shell_index, shell_slot_index, mpc_to_pc(r["edge_mpc"]), pc_to_ly(mpc_to_pc(r["edge_mpc"])),
            )
        else:
            designation = None
        results.append({
            "id": r["id"], "name": r["name"],
            "x": r["center_x_pc"], "y": r["center_y_pc"], "z": r["center_z_pc"],
            "galactic_radius_pc": r["galactic_radius_pc"],
            "shell_index": shell_index, "shell_slot_index": shell_slot_index,
            "designation": designation,
            "system_count": r["system_count"],
            "edge_ly": pc_to_ly(mpc_to_pc(r["edge_mpc"])) if r["edge_mpc"] else None,
            "distance_pc": math.sqrt(distance_sq),
        })
    return results


def galaxy_view(conn, center_x_pc, center_y_pc, center_z_pc, radius_pc):
    """
    The interactive 3D Galaxy Map's full live-viewport payload: every
    real, already-generated sector nearby (`galaxy_sectors_in_view`), every
    real, not-yet-generated sector address this galaxy's own density model
    predicts would qualify (`galaxyViewport.planned_slots_in_view`, only
    enumerated up to its own `PLANNED_RADIUS_CAP_PC` -- see that module's
    docstring), and -- for whatever's left of `radius_pc` beyond that cap
    -- a coarse illustrative density point cloud
    (`galaxyViewport.density_sample_points`).

    Works even before `generate.py plan` has ever been run: with no stored
    `galaxy_shape`, planned slots are returned unfiltered (every enumerated
    address, not just "qualifying" ones -- there's no density model yet to
    qualify them against) and the density tier is simply empty (nothing to
    sample), so the address scheme itself (shells/slots -- independent of
    any density/population model, see `docs/design/galaxy-coordinate-
    system.md` section 3) is still fully explorable.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        center_x_pc, center_y_pc, center_z_pc (float): The view center,
            galaxy-frame parsecs.
        radius_pc (float): The view radius, parsecs.

    Returns:
        dict: `placed` (`galaxy_sectors_in_view`'s own list), `planned`
            (`planned_slots_in_view`'s own list -- covers up to
            `PLANNED_RADIUS_CAP_PC` of the requested `radius_pc`, whatever
            that is; the address scheme has no reason to leave this
            empty), `density` (`density_sample_points`'s own list, `[]`
            whenever `radius_pc` doesn't exceed `PLANNED_RADIUS_CAP_PC` --
            `planned` already covers the whole view exactly in that case,
            so there's nothing left for an illustrative tier to add),
            `edge_pc`, `has_shape` (bool -- whether a real density model
            gates `planned`'s own qualification, i.e. whether `generate.py
            plan` has been run).
    """
    skeleton = get_galaxy_shape(conn)
    if skeleton is not None:
        edge_pc = skeleton.edge_pc
        shape = skeleton.shape
        expected_system_count = skeleton.expected_system_count_at_density_1
    else:
        edge_pc = ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
        shape = None
        expected_system_count = None
    edge_ly = pc_to_ly(edge_pc)

    center_pc = (center_x_pc, center_y_pc, center_z_pc)

    placed = galaxy_sectors_in_view(conn, center_x_pc, center_y_pc, center_z_pc, radius_pc)
    exclude_addresses = {
        (sector["shell_index"], sector["shell_slot_index"])
        for sector in placed
        if sector["shell_index"] is not None and sector["shell_slot_index"] is not None
    }

    planned = planned_slots_in_view(
        center_pc, radius_pc, edge_pc, edge_ly, shape, expected_system_count, exclude_addresses,
    )

    # The density cloud is only worth sampling for the part of the view
    # planned_slots_in_view's own radius cap couldn't cover with exact
    # addresses -- a view already entirely within that cap has nothing
    # left for an illustrative tier to add.
    density = density_sample_points(center_pc, radius_pc, shape) if radius_pc > PLANNED_RADIUS_CAP_PC else []

    return {
        "placed": placed, "planned": planned, "density": density,
        "edge_pc": edge_pc, "has_shape": shape is not None,
    }



# ---------------------------------------------------------------------
# Cube tiles -- backs GET /api/galaxy/tiles and GET /api/galaxy/stamp. See
# `stellarObjects.galaxyViewport`'s "Cube tiles" section for the model.
# ---------------------------------------------------------------------

GALAXY_TILE_MAX_PLACED = 250
"""int: Most placed sectors one tile returns. Small tiles never hold this
many; for big, zoomed-out tiles it's a sample (lowest ids first, so the
same sample every time) -- a zoomed-out view can't show thousands of
sub-pixel sectors anyway, and zooming in switches to smaller tiles that
hold the rest."""

MAX_TILES_PER_REQUEST = 128
"""int: Most tiles one `/api/galaxy/tiles` request may ask for. The map
needs at most 27 view tiles plus about 64 planned tiles at once."""


def galaxy_sectors_in_box(conn, lo, hi, limit=GALAXY_TILE_MAX_PLACED):
    """
    Placed sectors whose center lies in the half-open box `[lo, hi)`,
    lowest id first, at most `limit` -- the "placed" half of one tile.

    Two queries rather than `galaxy_sectors_in_view`'s one: the first
    reads only the (indexed, `idx_sectors_center`) sector columns with a
    `LIMIT`, and the per-sector system count runs only for the sectors
    actually returned, so a box covering the whole galaxy costs one
    bounded read instead of a count for every placed sector.

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
        SELECT sec.id, sec.name, sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               sec.galactic_radius_pc, sec.shell_index, sec.shell_slot_index, sec.edge_mpc
        FROM sectors sec
        WHERE sec.center_x_pc >= ? AND sec.center_x_pc < ?
          AND sec.center_y_pc >= ? AND sec.center_y_pc < ?
          AND sec.center_z_pc >= ? AND sec.center_z_pc < ?
        ORDER BY sec.id
        LIMIT ?
        """,
        (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], int(limit)),
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
    shell_index, shell_slot_index = r["shell_index"], r["shell_slot_index"]
    edge_pc = mpc_to_pc(r["edge_mpc"]) if r["edge_mpc"] else None
    if shell_index is not None and shell_slot_index is not None and edge_pc:
        designation = provisional_sector_designation(shell_index, shell_slot_index, edge_pc, pc_to_ly(edge_pc))
    else:
        designation = None
    return {
        "id": r["id"], "name": r["name"],
        "x": r["center_x_pc"], "y": r["center_y_pc"], "z": r["center_z_pc"],
        "galactic_radius_pc": r["galactic_radius_pc"],
        "shell_index": shell_index, "shell_slot_index": shell_slot_index,
        "designation": designation,
        "system_count": system_count,
        "edge_ly": pc_to_ly(edge_pc) if edge_pc else None,
    }


def galaxy_tiles(conn, tile_keys, density_key=None):
    """
    The contents of each requested cube tile, plus optionally one density
    cloud -- the interactive 3D Galaxy Map's data source. Every part of
    the result depends only on its tile key and the database's contents
    (see `galaxy_content_stamp`), so callers can cache each part by key.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        tile_keys (list[str]): `"level/ix/iy/iz"` keys (see
            `galaxyViewport.parse_tile_key`), at most
            `MAX_TILES_PER_REQUEST`.
        density_key (str or None): A tile key to anchor a density cloud on
            (`galaxyViewport.density_points_for_tile`), or `None` for none.

    Returns:
        dict: `tiles` (`{key: {"placed": [...], "planned": [...]}}`, see
            `galaxy_sectors_in_box`/`galaxyViewport.planned_slots_in_tile`),
            `density` (`{"key": density_key, "points": [...]}`, or `None`
            when no `density_key` was given; `points` is empty when the
            galaxy has no shape yet), `edge_pc`, `has_shape`.

    Raises:
        ValueError: On a malformed key or too many keys.
    """
    parsed = [(key, parse_tile_key(key)) for key in dict.fromkeys(tile_keys)]
    if len(parsed) > MAX_TILES_PER_REQUEST:
        raise ValueError(f"at most {MAX_TILES_PER_REQUEST} tiles per request, got {len(parsed)}")
    density_tile = parse_tile_key(density_key) if density_key else None

    skeleton = get_galaxy_shape(conn)
    if skeleton is not None:
        edge_pc = skeleton.edge_pc
        shape = skeleton.shape
        expected_system_count = skeleton.expected_system_count_at_density_1
    else:
        edge_pc = ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
        shape = None
        expected_system_count = None
    edge_ly = pc_to_ly(edge_pc)

    tiles = {}
    for key, (level, ix, iy, iz) in parsed:
        lo, hi = tile_bounds_pc(level, ix, iy, iz)
        placed = galaxy_sectors_in_box(conn, lo, hi)
        exclude_addresses = {
            (sector["shell_index"], sector["shell_slot_index"])
            for sector in placed
            if sector["shell_index"] is not None and sector["shell_slot_index"] is not None
        }
        planned = planned_slots_in_tile(
            level, ix, iy, iz, edge_pc, edge_ly, shape, expected_system_count, exclude_addresses,
        )
        tiles[key] = {"placed": placed, "planned": planned}

    density = None
    if density_tile is not None:
        points = density_points_for_tile(*density_tile, shape) if shape is not None else []
        density = {"key": density_key, "points": points}

    return {"tiles": tiles, "density": density, "edge_pc": edge_pc, "has_shape": shape is not None}


def galaxy_content_stamp(conn):
    """
    A short token that changes whenever anything `galaxy_tiles` returns
    could change: the placed-sector set (count and highest id), the
    highest star-system id (new systems change a sector's system count),
    the stored galaxy shape, and this code's own version (so a release
    that changes the tile format never reuses old cached tiles). The web
    layer and the browser key their tile caches by it, so a stale tile is
    never served after new sectors are generated.

    Cheap by design -- one indexed aggregate and two primary-key maxima --
    since it runs once per Galaxy Map page load.

    Returns:
        str: 16 hex characters.
    """
    sector_row = conn.execute(
        "SELECT COUNT(*) AS n, COALESCE(MAX(id), 0) AS max_id FROM sectors WHERE center_x_pc IS NOT NULL"
    ).fetchone()
    system_row = conn.execute("SELECT COALESCE(MAX(id), 0) AS max_id FROM star_systems").fetchone()
    parts = {
        "sectors": [int(sector_row["n"]), int(sector_row["max_id"])],
        "systems": int(system_row["max_id"]),
        "shape": galaxy_density_shape(conn),
        "version": __version__,
    }
    digest = hashlib.sha256(json.dumps(parts, sort_keys=True, default=str).encode("utf-8"))
    return digest.hexdigest()[:16]

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
)

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


def _search_name_list(conn, table, limit=SEARCH_AUTOCOMPLETE_LIMIT):
    rows = conn.execute(f"SELECT DISTINCT name FROM {table} ORDER BY name LIMIT ?", (limit,)).fetchall()
    return [row["name"] for row in rows]


# --- Result panels ---

def _search_result_sectors(conn, term):
    rows = conn.execute(
        "SELECT id, name, edge_mpc FROM sectors WHERE name LIKE ? ESCAPE '\\\\' ORDER BY name LIMIT ?",
        (_search_like_pattern(term), SEARCH_RESULT_LIMIT + 1),
    ).fetchall()
    truncated = len(rows) > SEARCH_RESULT_LIMIT
    return {"rows": [dict(r) for r in rows[:SEARCH_RESULT_LIMIT]], "truncated": truncated}


def _search_result_systems(conn, term):
    rows = conn.execute(
        """
        SELECT ss.id, ss.name, ss.sector_id, ss.is_binary, ss.binary_configuration, ss.binary_type,
               (SELECT s.star_type FROM stars s WHERE s.star_system_id = ss.id AND s.role = 'single' LIMIT 1)
                   AS single_star_type,
               (SELECT s.star_type FROM stars s WHERE s.star_system_id = ss.id AND s.role = 'primary' LIMIT 1)
                   AS primary_star_type,
               (SELECT s.star_type FROM stars s WHERE s.star_system_id = ss.id AND s.role = 'secondary' LIMIT 1)
                   AS secondary_star_type
        FROM star_systems ss
        WHERE ss.name LIKE ? ESCAPE '\\\\'
        ORDER BY ss.name
        LIMIT ?
        """,
        (_search_like_pattern(term), SEARCH_RESULT_LIMIT + 1),
    ).fetchall()
    truncated = len(rows) > SEARCH_RESULT_LIMIT
    rows = rows[:SEARCH_RESULT_LIMIT]
    return {
        "rows": [
            {
                "id": r["id"], "name": r["name"], "sector_id": r["sector_id"], "is_binary": r["is_binary"],
                "star_summary": _star_summary(r),
            }
            for r in rows
        ],
        "truncated": truncated,
    }


def _search_result_stars(conn, spectral_tags, luminosity_tags, term, size_range=None):
    clauses, params = [], []
    if spectral_tags:
        clauses.append(f"SUBSTR(s.star_type, 1, 1) IN ({','.join('?' * len(spectral_tags))})")
        params.extend(sorted(spectral_tags))
    if luminosity_tags:
        clauses.append(f"s.yerkes_class IN ({','.join('?' * len(luminosity_tags))})")
        params.extend(sorted(luminosity_tags))
    _append_size_clause(clauses, params, "s.radius_km", size_range)
    if term:
        clauses.append("s.name LIKE ? ESCAPE '\\\\'")
        params.append(_search_like_pattern(term))
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    params.append(SEARCH_RESULT_LIMIT + 1)
    rows = conn.execute(
        f"""
        SELECT s.name, s.role, s.star_type, s.radius_km, s.star_system_id, ss.name AS system_name, ss.sector_id
        FROM stars s
        JOIN star_systems ss ON ss.id = s.star_system_id
        WHERE 1=1{where}
        ORDER BY ss.name, s.name
        LIMIT ?
        """,
        params,
    ).fetchall()
    truncated = len(rows) > SEARCH_RESULT_LIMIT
    return {"rows": [dict(r) for r in rows[:SEARCH_RESULT_LIMIT]], "truncated": truncated}


def _search_result_planets(conn, class_tags, body_tags, life_tags, term, size_range=None):
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
        clauses.append("p.name LIKE ? ESCAPE '\\\\'")
        params.append(_search_like_pattern(term))
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    params.append(SEARCH_RESULT_LIMIT + 1)
    rows = conn.execute(
        f"""
        SELECT p.name, p.planet_class, p.body_type, p.life_chemical, p.radius_km,
               p.star_system_id, ss.name AS system_name, ss.sector_id
        FROM planets p
        JOIN star_systems ss ON ss.id = p.star_system_id
        WHERE 1=1{where}
        ORDER BY ss.name, p.orbital_index
        LIMIT ?
        """,
        params,
    ).fetchall()
    truncated = len(rows) > SEARCH_RESULT_LIMIT
    return {"rows": [dict(r) for r in rows[:SEARCH_RESULT_LIMIT]], "truncated": truncated}


def _search_result_moons(conn, class_tags, body_tags, life_tags, term, size_range=None):
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
        clauses.append("m.name LIKE ? ESCAPE '\\\\'")
        params.append(_search_like_pattern(term))
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    params.append(SEARCH_RESULT_LIMIT + 1)
    rows = conn.execute(
        f"""
        SELECT m.name, m.planet_class, m.body_type, m.life_chemical, m.radius_km, p.name AS planet_name,
               m.star_system_id, ss.name AS system_name, ss.sector_id
        FROM moons m
        JOIN planets p ON p.id = m.planet_id
        JOIN star_systems ss ON ss.id = m.star_system_id
        WHERE 1=1{where}
        ORDER BY ss.name, p.orbital_index, m.orbital_index
        LIMIT ?
        """,
        params,
    ).fetchall()
    truncated = len(rows) > SEARCH_RESULT_LIMIT
    return {"rows": [dict(r) for r in rows[:SEARCH_RESULT_LIMIT]], "truncated": truncated}


def _search_result_belts(conn, density_tags):
    clauses, params = [], []
    if density_tags:
        clauses.append(f"ab.density IN ({','.join('?' * len(density_tags))})")
        params.extend(sorted(density_tags))
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    params.append(SEARCH_RESULT_LIMIT + 1)
    rows = conn.execute(
        f"""
        SELECT ab.density, ab.composition_summary, ab.star_system_id, ss.name AS system_name, ss.sector_id
        FROM asteroid_belts ab
        JOIN star_systems ss ON ss.id = ab.star_system_id
        WHERE 1=1{where}
        ORDER BY ss.name, ab.orbital_index
        LIMIT ?
        """,
        params,
    ).fetchall()
    truncated = len(rows) > SEARCH_RESULT_LIMIT
    return {"rows": [dict(r) for r in rows[:SEARCH_RESULT_LIMIT]], "truncated": truncated}


def search(conn, texts, tags, sizes=None):
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

    Returns:
        dict: `facets` (`{facet: [{"value","label","count","tooltip"}, ...]}`,
            one list per `SEARCH_TAG_FACETS` name), `autocomplete`
            (`sectors`/`systems`/`stars`/`planets`/`moons` -> list of
            distinct names), `facet_labels` (`{"facet:value": label}`,
            for rendering an active-filter chip without a second lookup),
            and `results` (`sectors`/`systems`/`stars`/`planets`/`moons`/
            `belts` -> `{"rows": [...], "truncated": bool}`, or `None`
            for a panel with no reason to run).
    """
    sizes = sizes or {}
    star_size, planet_size, moon_size = sizes.get("star"), sizes.get("planet"), sizes.get("moon")

    spectral_tags, luminosity_tags = tags.get("spectral", set()), tags.get("luminosity", set())
    class_tags, body_tags, life_tags = tags.get("class", set()), tags.get("body", set()), tags.get("life", set())
    moon_class_tags, moon_body_tags = tags.get("moon_class", set()), tags.get("moon_body", set())
    moon_life_tags = tags.get("moon_life", set())
    density_tags = tags.get("density", set())
    type_tags = tags.get("type", set())

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
    )
    facets = {name: options for name, options in facet_defs}
    facet_labels = {
        f"{name}:{opt['value']}": opt["label"]
        for name, options in facet_defs for opt in options
    }

    autocomplete = {
        "sectors": _search_name_list(conn, "sectors"),
        "systems": _search_name_list(conn, "star_systems"),
        "stars": _search_name_list(conn, "stars"),
        "planets": _search_name_list(conn, "planets"),
        "moons": _search_name_list(conn, "moons"),
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

    results = {"sectors": None, "systems": None, "stars": None, "planets": None, "moons": None, "belts": None}
    if texts.get("sector_q"):
        results["sectors"] = _search_result_sectors(conn, texts["sector_q"])
    if texts.get("system_q"):
        results["systems"] = _search_result_systems(conn, texts["system_q"])
    if stars_included:
        results["stars"] = _search_result_stars(
            conn, spectral_tags, luminosity_tags, texts.get("star_q", ""), size_range=star_size
        )
    if planets_included:
        results["planets"] = _search_result_planets(
            conn, class_tags, body_tags, life_tags, texts.get("planet_q", ""), size_range=planet_size
        )
    if moons_included:
        results["moons"] = _search_result_moons(
            conn, moon_class_tags, moon_body_tags, moon_life_tags, texts.get("moon_q", ""), size_range=moon_size
        )
    if belts_included:
        results["belts"] = _search_result_belts(conn, density_tags)

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
