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
`api/config.py`) rather than the read-write account `sectorGen.py`/
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
import math

import pymysql

from stellarObjects._db import add_mysql_connection_args, get_connection, mysql_config_from_args
from stellarObjects._version import VersionAction, version_banner
from stellarObjects.navGraph import build_knn_adjacency, shortest_path
from stellarObjects.navigation import course_between, warp_travel_times
from stellarObjects.program_constants import NAV_ADJACENCY_K
from stellarObjects.utils import milliparsecs_to_ly, pc_to_ly


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
        SELECT sec.id, sec.name, sec.edge_mpc, COUNT(ss.id) AS system_count
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
        {"id": r["id"], "name": r["name"], "edge_ly": milliparsecs_to_ly(r["edge_mpc"]),
         "system_count": r["system_count"]}
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

    if sector_id is not None:
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
        sector_id (int, optional): Restricts to one sector's systems.
        limit (int, optional): Caps the number of rows returned. `None`
            (the default -- what every existing caller of this function
            still gets) returns every matching system.
        offset (int, optional): Skips this many rows first. Ignored
            unless `limit` is also given; meaningless on its own.

    Returns:
        list[dict]: One row per matching system, with `id`, `name`,
                           `sector_id`, `is_binary`.
    """
    join_sql, where_sql, params = _systems_filter_clause(star_type_prefix, sector_id)
    query = f"SELECT DISTINCT ss.id, ss.name, ss.sector_id, ss.is_binary FROM star_systems ss{join_sql}{where_sql} ORDER BY ss.name"

    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params = params + [limit, offset or 0]

    return conn.execute(query, params).fetchall()


def count_systems(conn, star_type_prefix=None, sector_id=None):
    """
    Returns the total number of systems matching the same filters
    `list_systems` accepts, ignoring any pagination -- the denominator
    `list_systems(conn, limit=...)` callers (the API's `/api/systems`)
    need to report how many pages exist.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        star_type_prefix (str, optional): Same meaning as `list_systems`.
        sector_id (int, optional): Same meaning as `list_systems`.

    Returns:
        int: Total matching system count.
    """
    join_sql, where_sql, params = _systems_filter_clause(star_type_prefix, sector_id)
    query = f"SELECT COUNT(DISTINCT ss.id) AS n FROM star_systems ss{join_sql}{where_sql}"
    return conn.execute(query, params).fetchone()["n"]


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


def nav_between(conn, from_system_id, to_system_id, adjacency_k=NAV_ADJACENCY_K):
    """
    Resolves full NAV information between two systems: a direct course
    (distance/azimuth/altitude/warp travel times, from
    `stellarObjects.navigation`) plus an optimal route via adjacent
    systems (`stellarObjects.navGraph`), or raises if NAV doesn't apply to
    this pair.

    NAV availability rules (see docs/api.md's NAV section for the
    user-facing statement of these):
        - Either system not assigned to any sector -> unavailable.
        - Same sector -> available, scoped to that sector's own systems
          (sector-local positions).
        - Different sectors, both galaxy-placed -> available, scoped to
          every system in every galaxy-placed sector (absolute
          galaxy-frame positions).
        - Different sectors, either not galaxy-placed -> unavailable.

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        from_system_id (int): The `star_systems.id` to route from.
        to_system_id (int): The `star_systems.id` to route to.
        adjacency_k (int): Passed through to
            `navGraph.build_knn_adjacency` as `k`.

    Returns:
        dict: `scope` (`"sector"` or `"galaxy"`), `direct` (a
            `navigation.Course`), `warp_times` (a list of
            `navigation.WarpLeg`, for `direct.distance_ly`),
            `origin_position`/`destination_position` (the `(x, y, z)`
            light-year positions `direct` was computed from, in `scope`'s
            frame -- sector-local for `"sector"`, galaxy-frame for
            `"galaxy"`), and `route`: `None` if `from_system_id ==
            to_system_id` or no path exists through the adjacency graph,
            else `{"path": [...system ids...], "distance_ly": float,
            "positions": {system_id: (x, y, z), ...}}` (one entry per id in
            `path`, same frame as `origin_position`/`destination_position`
            -- for rendering the route, e.g. `html/lib/navmap.py`, without
            a second position lookup).

    Raises:
        ValueError: If either system id doesn't exist.
        NavUnavailable: If NAV doesn't apply to this pair, per the rules
            above. `str(exc)` explains why.
    """
    origin = _load_nav_endpoint(conn, from_system_id)
    destination = _load_nav_endpoint(conn, to_system_id)

    if origin["sector_id"] is None or destination["sector_id"] is None:
        raise NavUnavailable("NAV requires both systems to be assigned to a sector")

    if origin["sector_id"] == destination["sector_id"]:
        scope = "sector"
        positions = _sector_local_positions(conn, origin["sector_id"])
        origin_position, destination_position = origin["position_ly"], destination["position_ly"]
    else:
        if origin["galaxy_position_ly"] is None or destination["galaxy_position_ly"] is None:
            raise NavUnavailable(
                "NAV between different sectors requires both sectors to have a galaxy placement"
            )
        scope = "galaxy"
        positions = _galaxy_frame_positions(conn)
        origin_position, destination_position = origin["galaxy_position_ly"], destination["galaxy_position_ly"]

    direct = course_between(origin_position, destination_position)

    route = None
    if from_system_id != to_system_id:
        graph = build_knn_adjacency(positions, adjacency_k)
        found = shortest_path(graph, from_system_id, to_system_id)
        if found is not None:
            path, distance_ly = found
            route = {
                "path": path,
                "distance_ly": distance_ly,
                "positions": {system_id: positions[system_id] for system_id in path},
            }

    return {
        "scope": scope,
        "direct": direct,
        "warp_times": warp_travel_times(direct.distance_ly),
        "origin_position": origin_position,
        "destination_position": destination_position,
        "route": route,
    }


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

    # TODO: no subcommand here queries planets/moons directly -- `systems`
    # above only filters by star_type_prefix/sector_id, so "every Class D
    # planet smaller than Earth" or "sort these results by radius_km" can't
    # be asked of this CLI at all today. A `planets` subcommand (mirroring
    # `systems` above) would need --class (the existing planet_class
    # values, already exposed as a search facet in ../src/html/search.py)
    # plus a new --min-radius-km/--max-radius-km pair (or a --sort-by
    # radius_km flag) over the `planets`/`moons` tables' radius_km column.
    # See docs/TODO.md, "Investigate Further".

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
    finally:
        conn.close()


if __name__ == "__main__":
    main()
