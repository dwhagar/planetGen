#!/usr/bin/env python3
# planetgen.cli.species_check

"""
Consistency report (POP.14): lists every species placed outside its own environment limits.

A species may only be put where it could live (`planetgen.db.species_fit.check_placement`). This reads the stored
placements back and lists any that break that rule: a species whose homeworld no longer fits its limits (the world
was regenerated or its radiation refreshed, or a limit was edited). Nothing is changed. A species it lists is fixed
by regenerating its limits or the placement; there is no tolerance flag.

Run it as `python3 -m planetgen.cli.species_check`. It prints one line per species and exits 1 when it found any,
0 when every placement is consistent.

Usage:
    python3 -m planetgen.cli.species_check [--mysql-host HOST] [--mysql-port PORT] [--mysql-user USER]
                                           [--mysql-password PASSWORD] [--mysql-database DATABASE]
"""

import argparse
import sys

from planetgen.db import species_fit, store
from planetgen.db.query import open_readonly


def main(argv=None):
    parser = argparse.ArgumentParser(description="List every species placed outside its own environment limits.")
    store.add_mysql_connection_args(parser)
    args = parser.parse_args(argv)
    conn = open_readonly(store.mysql_config_from_args(args))
    try:
        found = species_fit.inconsistencies(conn)
        total = conn.execute("SELECT COUNT(*) AS n FROM species").fetchone()["n"]
    finally:
        conn.close()
    for item in found:
        print(f"{item['name']} (species {item['species_id']}): {item['reasons']}")
    print(f"{len(found)} of {total:,} species {'is' if len(found) == 1 else 'are'} placed outside "
          f"{'its' if len(found) == 1 else 'their'} limits." if found
          else f"All {total:,} species are placed inside their limits.")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
