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
import math

import pymysql

from stellarObjects._db import add_mysql_connection_args, get_connection, mysql_config_from_args
from stellarObjects._version import VersionAction, version_banner
from stellarObjects.navGraph import build_knn_adjacency, shortest_path
from stellarObjects.navigation import course_between, warp_travel_times
from stellarObjects.physical_constants import SPECTRAL_CLASS_COLORS
from stellarObjects.program_constants import NAV_ADJACENCY_K, PLANET_CLASSES
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
            `temperature_k`/`radius_km`/`luminosity_w`).

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
    }


def system_detail(conn, system_id):
    """
    Returns one system's full web-display detail: stars, planets (each
    with its own moons), asteroid belts, description text, and enough
    sector context to render `html/system.py` -- in one function, the
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
            `markdown_content`, `wikitext_content`, `stars` (id/role/name/
            star_type/mass_kg/radius_km/temperature_k/luminosity_w --
            `id` matches a `'wide'` binary's `planets`/`belts` rows' own
            `star_id`, disambiguating which star each orbits), `planets`
            (each a `planets` row, including its own `star_id`, plus its
            own `moons` list), `belts` (`asteroid_belts` rows, including
            `star_id`), and `sector_siblings` (`{id, name}` for every
            other system in the same sector, empty if standalone -- for
            linkifying `location`'s "nearest: ..." names without a second
            round trip).

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
        "markdown_content": system["markdown_content"], "wikitext_content": system["wikitext_content"],
        "stars": [dict(s) for s in stars],
        "planets": planets,
        "belts": [dict(b) for b in belts],
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
    use as a SQL LIKE pattern (paired with `ESCAPE '\\'` in the query)."""
    escaped = term.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")
    return f"%{escaped}%"


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
        "SELECT id, name, edge_mpc FROM sectors WHERE name LIKE ? ESCAPE '\\' ORDER BY name LIMIT ?",
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
        WHERE ss.name LIKE ? ESCAPE '\\'
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


def _search_result_stars(conn, spectral_tags, luminosity_tags, term):
    clauses, params = [], []
    if spectral_tags:
        clauses.append(f"SUBSTR(s.star_type, 1, 1) IN ({','.join('?' * len(spectral_tags))})")
        params.extend(sorted(spectral_tags))
    if luminosity_tags:
        clauses.append(f"s.yerkes_class IN ({','.join('?' * len(luminosity_tags))})")
        params.extend(sorted(luminosity_tags))
    if term:
        clauses.append("s.name LIKE ? ESCAPE '\\'")
        params.append(_search_like_pattern(term))
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    params.append(SEARCH_RESULT_LIMIT + 1)
    rows = conn.execute(
        f"""
        SELECT s.name, s.role, s.star_type, s.star_system_id, ss.name AS system_name, ss.sector_id
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


def _search_result_planets(conn, class_tags, body_tags, life_tags, term):
    # TODO: no way to filter or sort by planet size (planets.radius_km)
    # here -- see docs/TODO.md, "Open items" > "Search".
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
    if term:
        clauses.append("p.name LIKE ? ESCAPE '\\'")
        params.append(_search_like_pattern(term))
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    params.append(SEARCH_RESULT_LIMIT + 1)
    rows = conn.execute(
        f"""
        SELECT p.name, p.planet_class, p.body_type, p.life_chemical,
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


def _search_result_moons(conn, class_tags, body_tags, life_tags, term):
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
    if term:
        clauses.append("m.name LIKE ? ESCAPE '\\'")
        params.append(_search_like_pattern(term))
    where = (" AND " + " AND ".join(clauses)) if clauses else ""
    params.append(SEARCH_RESULT_LIMIT + 1)
    rows = conn.execute(
        f"""
        SELECT m.name, m.planet_class, m.body_type, m.life_chemical, p.name AS planet_name,
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


def search(conn, texts, tags):
    """
    Runs the faceted search behind `GET /api/search`/`html/search.py`:
    the same click-to-filter attribute tags (object type; star spectral/
    luminosity class; planet/moon class, body type, supported life
    chemistry; asteroid belt density) plus a per-entity name search this
    project's search page has always offered -- ported from
    `html/search.py`'s previous direct-SQL implementation unchanged (same
    queries, same "which result panels actually have a reason to run"
    logic: a panel only appears when one of its own tags/name field is
    active, or its object type is explicitly selected).

    Args:
        conn (stellarObjects._db.Connection): An open, read-only connection.
        texts (dict): `{"sector_q", "system_q", "star_q", "planet_q",
            "moon_q"}` -> search term (`""`/absent for "not searched").
        tags (dict): `{facet: set(value, ...)}`, one entry per name in
            `SEARCH_TAG_FACETS` -- an absent/empty set means no active
            filter for that facet.

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

    star_has_reason = bool(spectral_tags or luminosity_tags or texts.get("star_q"))
    planet_has_reason = bool(class_tags or body_tags or life_tags or texts.get("planet_q"))
    moon_has_reason = bool(moon_class_tags or moon_body_tags or moon_life_tags or texts.get("moon_q"))
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
        results["stars"] = _search_result_stars(conn, spectral_tags, luminosity_tags, texts.get("star_q", ""))
    if planets_included:
        results["planets"] = _search_result_planets(conn, class_tags, body_tags, life_tags, texts.get("planet_q", ""))
    if moons_included:
        results["moons"] = _search_result_moons(
            conn, moon_class_tags, moon_body_tags, moon_life_tags, texts.get("moon_q", "")
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

    # TODO: no subcommand here queries planets/moons directly -- `systems`
    # above only filters by star_type_prefix/sector_id, so "every Class D
    # planet smaller than Earth" or "sort these results by radius_km" can't
    # be asked of this CLI at all today. A `planets` subcommand (mirroring
    # `systems` above) would need --class (the existing planet_class
    # values, already exposed as a search facet in ../src/html/search.py)
    # plus a new --min-radius-km/--max-radius-km pair (or a --sort-by
    # radius_km flag) over the `planets`/`moons` tables' radius_km column.
    # See docs/TODO.md, "Open items" > "Search".

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
