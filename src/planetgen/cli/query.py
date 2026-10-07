# planetgen.cli.query

"""
List/query what's already stored in the planetGen database, from the
command line; the queries themselves are `planetgen.db.query`'s. Run from
the checkout's `src/`:

    python3 -m planetgen.cli.query [--mysql-host HOST] ... sectors|systems|near|planets|moons
"""

import argparse

from planetgen._version import VersionAction, version_banner
from planetgen.db.query import list_moons, list_planets, list_sectors, list_systems, open_readonly, systems_within_radius
from planetgen.db.store import add_mysql_connection_args, mysql_config_from_args


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
    parser.add_argument('--version', action=VersionAction, banner=version_banner('planetgen.cli.query'))
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
