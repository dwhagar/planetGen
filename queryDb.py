# queryDb.py

"""
List/query CLI for the planetGen database (stellarObjects/schema.sql).

A thin read-only front end over the tables `sectorGen.py`/`systemGen.py`
populate, for questions like "every G-type system," "everything within 50
light-years of a given system," and "what sectors exist" -- see TODO.md's
Phase 3 ("A way to list/query what's already stored"). Deliberately plain
SQL rather than routing through `stellarObjects._db`'s `load_star_system`/
`load_sector` (Phase 2's read path): these are simple, columnar listings,
not full object-graph reconstructions, so a raw query is the more direct
tool for the job -- the read path remains what a future richer tool (or a
re-upload/re-render workflow) would build on.

Opens the database strictly read-only (a `file:` URI with `mode=ro`) so
this can never accidentally write to a database another process is also
using.
"""

import argparse
import math
import sqlite3

from stellarObjects._db import DEFAULT_DB_PATH
from stellarObjects._version import VersionAction, version_banner
from stellarObjects.utils import milliparsecs_to_ly


def open_readonly(db_path):
    """
    Opens `db_path` strictly read-only -- fails outright if the file
    doesn't exist, rather than silently creating one (this tool never
    writes).

    Args:
        db_path (str): Path to the `.db` file.

    Returns:
        sqlite3.Connection: A read-only connection with `Row` row access.

    Raises:
        SystemExit: If the file doesn't exist or can't be opened.
    """
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        conn.row_factory = sqlite3.Row
        conn.execute("SELECT 1")  # forces the file to actually be opened now
    except sqlite3.OperationalError as exc:
        raise SystemExit(f"Error: could not open database at {db_path!r} ({exc}).")
    return conn


def list_sectors(conn):
    """
    Returns every sector, with its edge length (converted to light-years)
    and how many systems it contains.

    Args:
        conn (sqlite3.Connection): An open, read-only connection.

    Returns:
        list[sqlite3.Row]: One row per sector, with `id`, `name`,
                           `edge_ly`, `system_count`.
    """
    rows = conn.execute(
        """
        SELECT sec.id, sec.name, sec.edge_mpc, COUNT(ss.id) AS system_count
        FROM sectors sec
        LEFT JOIN star_systems ss ON ss.sector_id = sec.id
        GROUP BY sec.id
        ORDER BY sec.name
        """
    ).fetchall()
    return [
        {"id": r["id"], "name": r["name"], "edge_ly": milliparsecs_to_ly(r["edge_mpc"]),
         "system_count": r["system_count"]}
        for r in rows
    ]


def list_systems(conn, star_type_prefix=None, sector_id=None):
    """
    Returns systems, optionally filtered by star type and/or sector.

    Args:
        conn (sqlite3.Connection): An open, read-only connection.
        star_type_prefix (str, optional): Matches any system with at least
            one star (single, primary, or secondary) whose `star_type`
            starts with this text, e.g. `"G"` for every G-type system,
            `"G2V"` for an exact spectral/subclass/luminosity match.
            Case-sensitive, matching the stored spectral class letters.
        sector_id (int, optional): Restricts to one sector's systems.

    Returns:
        list[sqlite3.Row]: One row per matching system, with `id`, `name`,
                           `sector_id`, `is_binary`.
    """
    query = "SELECT DISTINCT ss.id, ss.name, ss.sector_id, ss.is_binary FROM star_systems ss"
    conditions = []
    params = []

    if star_type_prefix is not None:
        query += " JOIN stars s ON s.star_system_id = ss.id"
        conditions.append("s.star_type LIKE ?")
        params.append(f"{star_type_prefix}%")

    if sector_id is not None:
        conditions.append("ss.sector_id = ?")
        params.append(sector_id)

    if conditions:
        query += " WHERE " + " AND ".join(conditions)
    query += " ORDER BY ss.name"

    return conn.execute(query, params).fetchall()


def systems_within_radius(conn, system_id, radius_ly):
    """
    Finds every other system in the same sector as `system_id`, within
    `radius_ly` light-years, nearest first.

    Args:
        conn (sqlite3.Connection): An open, read-only connection.
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
    parser.add_argument('--db-path', type=str,
                        help=f"Path to the SQLite database file to query. Defaults to {DEFAULT_DB_PATH}.")

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

    return parser.parse_args()


def main():
    """
    The main entry point: dispatches to the requested subcommand and prints
    a plain-text listing of the results.
    """
    args = process_args()
    db_path = args.db_path or DEFAULT_DB_PATH
    conn = open_readonly(db_path)
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
