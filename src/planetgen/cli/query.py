# planetgen.cli.query

"""
List/query what's already stored in the planetGen database, from the
command line; the queries themselves are `planetgen.db.query`'s. Run from
the checkout's `src/`:

    python3 -m planetgen.cli.query [--mysql-host HOST] ... sectors|systems|near|planets|moons
"""

import argparse

from planetgen._version import VersionAction, version_banner
from planetgen.api import ids
from planetgen.db.query import list_moons, list_planets, list_sectors, list_systems, open_readonly
from planetgen.db import near
from planetgen.util import log
from planetgen.util.star_count import STAR_COUNT_NAMES, parse_star_count, star_count_label
from planetgen.util.format import format_distance_km, format_distance_ly, format_distance_pc
from planetgen.db.store import add_mysql_connection_args, mysql_config_from_args


STARS_HELP = ("Only systems with this many stars: a count (3 is exactly three), a count and a plus (3+ is three or more) or a name "
              f"({', '.join(name for name, _ in STAR_COUNT_NAMES)}; quintenary is five or more).")
"""str: The help of the `--stars` option of `systems` and `near` (UX.239): it names every choice."""

DEFAULT_PAGE = 100
"""int: How many sectors, systems, planets or moons a listing shows unless `--limit` says otherwise (PERF.207, PERF.211)."""


def _stars_argument(text):
    """The `--stars` value as `(least, most)` (`parse_star_count`); text it refuses is an argparse error."""
    try:
        return parse_star_count(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc))


def _add_paging(parser, kind):
    """The `--limit` and `--after` of a listing that can be very long (PERF.207)."""
    parser.add_argument('--limit', type=int, default=DEFAULT_PAGE,
                        help=f"How many {kind}s to show (default {DEFAULT_PAGE}; 0 shows every one, which on a big galaxy "
                             "reads them all).")
    parser.add_argument('--after', metavar='ID',
                        help=f"Carry on after this {kind} (the ID the previous listing ended on, which it prints).")


def _page(args):
    """`(rows to ask for, rows to show)` for a listing: one row more than shown, to know whether another page follows;
    both `None` for `--limit 0`."""
    if args.limit < 0:
        raise SystemExit("Error: --limit must be 0 or more.")
    return (None, None) if args.limit == 0 else (args.limit + 1, args.limit)


def _show(conn, kind, rows, shown):
    """The rows to print, and after them the `--after` to carry on with when more follow."""
    return rows if shown is None else rows[:shown], (
        None if shown is None or len(rows) <= shown else _printed(conn, kind, rows[shown - 1]["id"]))


def process_args():
    """
    Parses command-line arguments for the subcommands: `sectors`, `systems`, `near`, `planets` and `moons`.

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

    sectors_parser = subparsers.add_parser('sectors', help="List sectors, with their size and system count.")
    _add_paging(sectors_parser, "sector")

    systems_parser = subparsers.add_parser('systems', help="List systems, optionally filtered.")
    systems_parser.add_argument('--star-type', type=str,
                                help="Only systems with a star whose type starts with this "
                                     "(e.g. 'G' for every G-type system, 'G2V' for an exact match).")
    systems_parser.add_argument('--sector-id', help="Only systems in this sector.")
    systems_parser.add_argument('--stars', type=_stars_argument, help=STARS_HELP)
    _add_paging(systems_parser, "system")

    near_parser = subparsers.add_parser(
        'near', help="Find everything generated within a distance of a place.",
    )
    near_parser.add_argument('place', help="An object reference (system:<ID>, planet:<ID>, nebula:<ID> ...; a bare "
                                           "three-part ID is a system) or a galaxy-frame point 'x,y,z' in parsecs.")
    near_parser.add_argument('--distance', type=float,
                             help=f"Search distance in parsecs, up to {near.MAX_DISTANCE_PC:g}; with --nearest, how far "
                                  f"to look at most (default {near.MAX_DISTANCE_PC:g}).")
    near_parser.add_argument('--nearest', type=int,
                             help=f"The nearest N rows (up to {near.MAX_LIMIT}), found by trying {', '.join(f'{r:g}' for r in near.LADDER_PC)} "
                                  "pc in turn; quick even at the dense core.")
    near_parser.add_argument('--stars', type=_stars_argument,
                             help=STARS_HELP + " Phenomena and facilities outside a system are left out.")
    near_parser.add_argument('--kinds', help="Only these kinds, comma-separated: " + ", ".join(near.SEARCH_KINDS) + ".")
    near_parser.add_argument('--limit', type=int, default=near.DEFAULT_LIMIT, help="Rows to show (default 50).")
    near_parser.add_argument('--offset', type=int, default=0, help="Rows to skip.")

    planets_parser = subparsers.add_parser(
        'planets', help="List planets, optionally filtered by class, radius, sector, or system.",
    )
    planets_parser.add_argument('--class', dest='planet_class', type=str,
                                help="Only planets of this exact class (e.g. 'M').")
    planets_parser.add_argument('--min-radius-km', type=float, help="Only planets at least this large.")
    planets_parser.add_argument('--max-radius-km', type=float, help="Only planets at most this large.")
    planets_parser.add_argument('--sector-id', help="Only planets whose system is in this sector.")
    planets_parser.add_argument('--system-id', help="Only planets in this one system.")
    _add_paging(planets_parser, "planet")

    moons_parser = subparsers.add_parser(
        'moons', help="List moons, optionally filtered by class, radius, sector, or system.",
    )
    moons_parser.add_argument('--class', dest='planet_class', type=str,
                              help="Only moons of this exact class (e.g. 'M').")
    moons_parser.add_argument('--min-radius-km', type=float, help="Only moons at least this large.")
    moons_parser.add_argument('--max-radius-km', type=float, help="Only moons at most this large.")
    moons_parser.add_argument('--sector-id', help="Only moons whose system is in this sector.")
    moons_parser.add_argument('--system-id', help="Only moons in this one system.")
    _add_paging(moons_parser, "moon")

    log.add_debug_arguments(parser)
    args = parser.parse_args()
    if args.command == 'near' and args.distance is None and args.nearest is None:
        near_parser.error("give --distance (parsecs) or --nearest (a number of rows)")
    log.apply_debug_arguments(args, parser)
    return args


def _row(conn, kind, printed_id):
    """The row id of the object a printed ID names, or `None` for no filter; an ID that names nothing is an error."""
    if printed_id is None:
        return None
    found = ids.row_id(conn, kind, printed_id)
    if found is None:
        raise SystemExit(f"Error: no {kind} has the ID {printed_id}.")
    return found


def _printed(conn, kind, row_id):
    return ids.printed(conn, kind, row_id)


def main():
    """
    The main entry point: dispatches to the requested subcommand and prints
    a plain-text listing of the results.
    """
    args = process_args()
    conn = open_readonly(mysql_config_from_args(args))
    try:
        if args.command == 'sectors':
            ask, shown = _page(args)
            sectors = list_sectors(conn, limit=ask, after=_row(conn, "sector", args.after))
            if not sectors:
                print("No sectors stored." if args.after is None else "No more sectors.")
                return
            sectors, carry_on = _show(conn, "sector", sectors, shown)
            for sector in sectors:
                print(f"[{_printed(conn, 'sector', sector['id'])}] {sector['name']} "
                      f"(edge {format_distance_ly(sector['edge_ly'])}, {sector['system_count']} systems)")
            if carry_on:
                print(f"{len(sectors):,} sectors shown; more follow: carry on with --after {carry_on}")

        elif args.command == 'systems':
            ask, shown = _page(args)
            systems = list_systems(conn, star_type_prefix=args.star_type, sector_id=_row(conn, "sector", args.sector_id),
                                   limit=ask, after=_row(conn, "system", args.after), stars=args.stars)
            if not systems:
                print("No matching systems." if args.after is None else "No more systems.")
                return
            systems, carry_on = _show(conn, "system", systems, shown)
            for system in systems:
                stars = system.get("star_count") or (2 if system["is_binary"] else 1)
                kind = star_count_label((stars, stars)).lower()
                sector_note = (f"sector {_printed(conn, 'sector', system['sector_id'])}"
                               if system["sector_id"] is not None else "standalone")
                print(f"[{_printed(conn, 'system', system['id'])}] {system['name']} ({kind}, {sector_note})")
            if carry_on:
                print(f"{len(systems):,} systems shown; more follow: carry on with --after {carry_on}")

        elif args.command == 'near':
            try:
                if "," in args.place:
                    place = near.place_from_point(args.place.split(","))
                else:
                    resolved = ids.resolve_ref(conn, args.place)
                    if resolved is None:
                        raise near.NearError(f"{args.place!r} names nothing")
                    place = near.place_from_reference(conn, resolved)
                kinds = [kind for kind in args.kinds.split(",") if kind] if args.kinds else None
                if args.nearest is not None:
                    result = near.nearest(conn, place, args.nearest, kinds=kinds,
                                          max_distance_pc=args.distance or near.MAX_DISTANCE_PC, stars=args.stars)
                else:
                    result = near.objects_within(conn, place, args.distance, kinds=kinds,
                                                 limit=args.limit, offset=args.offset, stars=args.stars)
            except (near.NearError, ids.IdError) as exc:
                raise SystemExit(f"Error: {exc}")
            result = ids.translate(result, {"rows[].id": "by-kind"}, conn)
            if not result["rows"]:
                print(f"Nothing within {result['distance_pc']:g} pc of {place['name']}.")
            for row in result["rows"]:
                parent = f" -- {row['parent']['name']}" if row["parent"] else ""
                count = f" ({star_count_label((row['star_count'], row['star_count'])).lower()})" if "star_count" in row else ""
                print(f"[{row['ref']}] {row['name']}{count}{parent} -- {format_distance_pc(row['distance_pc'])}")
            capped = f" ({result['systems_in_range']:,}+ systems in range)" if result["systems_in_range_capped"] else ""
            print(f"{result['total']:,} objects found within {result['distance_pc']:g} pc{capped}; "
                  f"{result['not_generated_total']:,} more objects in range are not generated yet.")

        elif args.command == 'planets':
            ask, shown = _page(args)
            planets = list_planets(
                conn, planet_class=args.planet_class, min_radius_km=args.min_radius_km,
                max_radius_km=args.max_radius_km, sector_id=_row(conn, "sector", args.sector_id),
                system_id=_row(conn, "system", args.system_id), limit=ask, after=_row(conn, "planet", args.after),
            )
            if not planets:
                print("No matching planets." if args.after is None else "No more planets.")
                return
            planets, carry_on = _show(conn, "planet", planets, shown)
            for planet in planets:
                print(f"[{_printed(conn, 'planet', planet['id'])}] {planet['name']} (Class {planet['planet_class']}, "
                      f"{format_distance_km(planet['radius_km'])}) -- {planet['system_name']}")
            if carry_on:
                print(f"{len(planets):,} planets shown; more follow: carry on with --after {carry_on}")

        elif args.command == 'moons':
            ask, shown = _page(args)
            moons = list_moons(
                conn, planet_class=args.planet_class, min_radius_km=args.min_radius_km,
                max_radius_km=args.max_radius_km, sector_id=_row(conn, "sector", args.sector_id),
                system_id=_row(conn, "system", args.system_id), limit=ask, after=_row(conn, "moon", args.after),
            )
            if not moons:
                print("No matching moons." if args.after is None else "No more moons.")
                return
            moons, carry_on = _show(conn, "moon", moons, shown)
            for moon in moons:
                print(f"[{_printed(conn, 'moon', moon['id'])}] {moon['name']} (Class {moon['planet_class']}, "
                      f"{format_distance_km(moon['radius_km'])}) -- {moon['system_name']} / {moon['planet_name']}")
            if carry_on:
                print(f"{len(moons):,} moons shown; more follow: carry on with --after {carry_on}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
