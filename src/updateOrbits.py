#!/usr/bin/env python3
# src/updateOrbits.py

"""
Advances every planet's, moon's, star's, and binary system's orbital
position in the configured database based on real elapsed time since the
last run -- the "dedicated update script" `docs/TODO.md`'s orbital-motion
entry called for, meant to be run periodically (e.g. via cron, "once a
month or so") rather than on every generation run. This is the one
script that has to touch every table with a floating-point position/phase
column -- `planets`/`moons` (their own orbit) and `stars`/`star_systems`
(a star's galactic orbit, plus a binary pair's mutual orbit around each
other) -- so a single run brings the whole database's motion up to date
in one pass.

`stellarObjects._db.advance_orbital_phases` does the actual work: a
single set-based `UPDATE` per table (`planets`/`moons`), driven by each
body's own already-stored `period_years` and the elapsed real time since
`orbit_simulation_state.last_updated_at` (stored, not recomputed --
`stellarObjects._db.get_orbit_update_elapsed_years` measures it
server-side via `TIMESTAMPDIFF` rather than trusting this process' own
clock to agree with the database server's). A database this has never run
against yet has no `orbit_simulation_state` row -- the first run just
establishes that reference point (zero elapsed time, nothing to advance
yet) rather than guessing a start time.

`position_x/y/z_km` are recomputed in lockstep with `orbital_phase_deg`
(both are handled by the same `advance_orbital_phases` call -- position is
a pure function of distance/inclination/ascending-node/phase, so it has no
independent update of its own). A body is skipped entirely (no `UPDATE`
attempted at all, not just a no-op write) when this run's elapsed time is
below that body's own `min_update_interval_years` -- the point past which
the phase delta added would be smaller than `orbital_phase_deg`'s own
floating-point resolution and so is guaranteed to round back to the exact
value already stored (see `stellarObjects.utils.minimum_update_interval_years`).
`orbital_inclination_deg`/`orbital_ascending_node_deg`/`orbital_speed_kms`/
`min_update_interval_years` (fixed at generation time) and
`rotation_period_hours` (a static descriptive stat -- this generator
doesn't track rotational phase) are untouched; see
`stellarObjects.planetPhysics.generate_orbital_motion_properties`.

Every star's `galactic_orbital_phase_deg` (its position around the galactic
center) and, for a binary pair, `star_systems.binary_galactic_orbital_phase_deg`
and `binary_mutual_orbital_phase_deg` (the pair's own mutual orbit around
each other, entirely separate from their shared galactic orbit) are
advanced the same way, each guarded by its own `*_min_update_interval_years`
-- see `schema.sql`'s "v13" note and `stellarObjects._db.advance_orbital_phases`'s
docstring. `binary_mutual_position_x/y/z_km` (the secondary's position
relative to the primary) are recomputed in lockstep with
`binary_mutual_orbital_phase_deg`, the same "position has no independent
update of its own" treatment `position_x/y/z_km` gets above -- see
`schema.sql`'s "v14" note.

This file lives alongside `stellarObjects/` under `src/`, so Python's own
sys.path[0] (the running script's directory) already makes
`stellarObjects` importable -- no sys.path shim needed.

Usage:
    python3 src/updateOrbits.py [--mysql-host HOST] [--mysql-port PORT]
                                [--mysql-user USER] [--mysql-password PASSWORD]
                                [--mysql-database DATABASE]

    Every flag defaults to the same $PLANETGEN_MYSQL_* environment
    variable every other entry point in this project reads (see
    `stellarObjects._db.MySQLConfig`). Needs a read-write database
    account (like `sectorGen.py`/`systemGen.py`, not `queryDb.py`'s
    read-only one) -- this script mutates rows.
"""

import argparse
import sys

from stellarObjects._db import (
    add_mysql_connection_args,
    advance_orbital_phases,
    get_connection,
    get_orbit_update_elapsed_years,
    mysql_config_from_args,
)
from stellarObjects._version import VersionAction, version_banner


def main():
    parser = argparse.ArgumentParser(
        description="Advance every planet's, moon's, star's, and binary system's orbital position "
                    "based on real elapsed time.",
    )
    add_mysql_connection_args(parser)
    parser.add_argument('--version', action=VersionAction, banner=version_banner('updateOrbits.py'))
    args = parser.parse_args()

    config = mysql_config_from_args(args)

    try:
        conn = get_connection(config)
    except Exception as exc:
        print(f"error: could not open the database ({exc}).", file=sys.stderr)
        sys.exit(1)

    try:
        elapsed_years = get_orbit_update_elapsed_years(conn)
        if elapsed_years is None:
            print("No previous orbit update found for this database -- establishing a starting point now "
                  "(nothing to advance yet; run this again later to actually move anything).")
            elapsed_years = 0.0
        else:
            print(f"{elapsed_years:.6f} years elapsed since the last update -- advancing orbits.")

        planets_updated, moons_updated, stars_updated, binary_systems_updated = \
            advance_orbital_phases(conn, elapsed_years)
        print(
            f"Updated {planets_updated} planet(s), {moons_updated} moon(s), "
            f"{stars_updated} star(s), and {binary_systems_updated} binary system(s)."
        )
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
    finally:
        conn.close()


if __name__ == "__main__":
    main()
