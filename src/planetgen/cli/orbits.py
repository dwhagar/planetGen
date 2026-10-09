#!/usr/bin/env python3
# planetgen.cli.orbits

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

Each object keeps its own clock (GEN.106, schema v66): `epoch_unix`, when
its stored position holds, and an indexed `next_update_due`, when it will
have moved far enough to be worth storing again -- 0.01 mpc on a galactic
orbit (stars, black holes, neutron stars, every other star-like object
outside a system), 0.01 AU in a system (planets, comets, a pair's mutual
orbit, facilities round a star or in a belt), 100,000 km round a planet
(moons and the facilities orbiting a planet or moon) -- worked out from
its speed (`planetgen.physics.position.update_interval_s`). A run
(`planetgen.db.store.orbit_clock`) moves each object that is due from its
own epoch to now, by the database server's clock, sets its epoch to now
and works out its next due time; an object that isn't due is neither moved
nor counted. Rows without a due time yet (new, or edited since) get one
first, from the last run's time (`orbit_simulation_state.last_updated_at`,
which `finish_orbit_update` sets at the end of a run).

`planetgen.db.store.advance_orbital_phases` moves planets and moons (their
own orbit) and binary pairs' mutual orbits: a set-based `UPDATE` per table,
driven by each body's stored `period_years`. `position_x/y/z_km` are
recomputed in lockstep with `orbital_phase_deg` (position is a pure
function of distance/inclination/ascending-node/phase).
`orbital_inclination_deg`/`orbital_ascending_node_deg`/`orbital_speed_kms`
(fixed at generation time) and `rotation_period_hours` (a static
descriptive stat -- this generator doesn't track rotational phase) are
untouched; see `planetgen.physics.planets.generate_orbital_motion_properties`.

Every star's `galactic_orbital_phase_deg` (its position around the galactic
center) and a close pair's `star_systems.binary_galactic_orbital_phase_deg`
advance with their system's galactic move (`advance_galactic_positions`,
below), on the system's clock; a pair's `binary_mutual_orbital_phase_deg`
(its own mutual orbit, entirely separate from the shared galactic orbit)
has a clock of its own (`binary_epoch_unix`/`binary_next_update_due`).
`binary_mutual_position_x/y/z_km` (the secondary's position relative to
the primary) are recomputed in lockstep with
`binary_mutual_orbital_phase_deg`, the same "position has no independent
update of its own" treatment `position_x/y/z_km` gets above -- see
`schema.sql`'s "v14" note.

Every standalone exotic phenomenon (`phenomenonGen.py`, schema v16/v17) --
a black hole/neutron star with no owning `StarSystem`, a nebula, a
supernova remnant, a rogue planet, an interstellar comet, or a standalone
asteroid field -- gets its own `galactic_orbital_phase_deg` advanced with its galactic move,
on its own clock, the way a lone star's is: unbound from any specific STAR doesn't
mean unbound from the galaxy itself, so these still orbit the galactic
center on the same timescale. See `schema.sql`'s "v17" header note.

Star-bound comets (schema v18, `planetgen.db.store.advance_comet_orbits`)
are handled by a SEPARATE call, not folded into `advance_orbital_phases`
above: a comet's position isn't a linear function of elapsed time the way
a circular planet/moon orbit's is (Kepler's second law -- it moves far
faster near perihelion), so turning its advanced orbital anomaly into a
distance/position requires solving Kepler's or Barker's equation in
Python, not a set-based SQL `UPDATE`. See `advance_comet_orbits`'s own
docstring and `docs/design/comet-orbital-realism.md`.

Schema v20 layers a proper two-body (barycentric) "reflex offset"/
"wobble" on top of the above, for every relationship where the orbited
body isn't overwhelmingly more massive than what orbits it: a
planet-hosting star's own small displacement from its planets' combined
pull, a moon-hosting planet's from its moons', and (for a binary pair)
each star's own offset from the pair's shared barycenter rather than one
star sitting fixed while the other orbits it. None of the position
columns above change meaning -- these are new, additive columns
(`reflex_offset_x/y/z_km` on `stars`/`planets`,
`binary_primary/secondary_position_*_km` and
`binary_planetary_wobble_*_km` on `star_systems`), recomputed fresh on
every run from whatever the already-advanced children currently look
like, with no due time of their own. See
`schema.sql`'s "v20" header note and
`planetgen.physics.orbits.calculate_reflex_offset`'s docstring for the
formula.

Galactic motion (GEN.6): after the phases above,
`planetgen.db.store.advance_galactic_positions` turns every due star
system, standalone phenomenon and stand-alone facility about the galactic
axis by the angle its galactic phase advances (the galaxy's nucleus,
the quasar, stays put). Sectors are fixed cells, so an object that drifts
into another generated sector is refiled there: `sector_id`, its
sector-relative position and octant. `refresh_after_motion` then
recomputes containment (`refresh_containment`), octants and the stored
nearest systems (`refresh_nearest_systems`) for every placed sector, and
rewrites the `location` text of every system whose sector or nearest
systems changed. Last, the saved sector paths (GEN.123,
`planetgen.db.sector_paths`) are recomputed for every sector holding a star
system, rogue planet or comet. Orbital facilities advance their orbit phase like moons
(`advance_facility_orbits`). Page text is rendered from these rows
(schema v29), so nothing else names the old sector.

Run it as `python3 -m planetgen.cli.orbits` from anywhere (the editable
install makes `planetgen` importable).

Usage:
    python3 -m planetgen.cli.orbits [--mysql-host HOST] [--mysql-port PORT]
                                [--mysql-user USER] [--mysql-password PASSWORD]
                                [--mysql-database DATABASE]

    Every flag defaults to the same $PLANETGEN_MYSQL_* environment
    variable every other entry point in this project reads (see
    `planetgen.db.store.MySQLConfig`). Needs a read-write database
    account (like `sectorGen.py`/`systemGen.py`, not `planetgen.db.query`'s
    read-only one) -- this script mutates rows.
"""

import argparse
import sys

from planetgen.cli.stage_progress import StageProgress
from planetgen.db.sector_paths import compute_sector_paths
from planetgen.db.store import (
    add_mysql_connection_args,
    advance_comet_orbits,
    advance_facility_orbits,
    advance_galactic_positions,
    advance_orbital_phases,
    finish_orbit_update,
    get_connection,
    get_orbit_update_elapsed_years,
    mysql_config_from_args,
    orbit_clock,
    refresh_after_motion,
)
from planetgen._version import VersionAction, version_banner


TABLE_LABELS = {
    "planets": "planet(s)",
    "moons": "moon(s)",
    "binary_mutual_orbits": "binary mutual orbit(s)",
    "comets": "star-bound comet(s)",
    "facilities": "orbital facilit(ies)",
    "star_systems": "star system(s) on their galactic orbits",
    "binary_galactic_orbits": "close binary galactic orbit(s)",
    "black_holes": "standalone black hole(s)",
    "neutron_stars": "standalone neutron star(s)",
    "nebulae": "nebula(e)",
    "supernova_remnants": "supernova remnant(s)",
    "rogue_planets": "rogue planet(s)",
    "interstellar_comets": "interstellar comet(s)",
    "asteroid_fields": "asteroid field(s)",
    "space_facilities": "stand-alone facilit(ies)",
}
"""dict: The count keys of `advance_orbital_phases`, `advance_comet_orbits`,
`advance_facility_orbits` and `advance_galactic_positions` -> the label in
this script's summary line. Each counts only objects that were due (GEN.106);
the derived wobbles, recomputed every run, aren't listed."""


STAGES = (
    "Advancing orbital phases",
    "Advancing comets and facilities",
    "Moving along galactic orbits",
    "Refreshing containment, nearest systems and locations",
    "Recomputing sector paths",
)
"""tuple: The steps of a run, in order, as the progress bar names them."""


def _stage_label(index):
    return f"Step {index + 1} of {len(STAGES)}: {STAGES[index]}"


def main():
    parser = argparse.ArgumentParser(
        description="Advance every planet's, moon's, star's, binary system's, and standalone exotic "
                    "phenomenon's orbital position based on real elapsed time.",
    )
    add_mysql_connection_args(parser)
    parser.add_argument('--version', action=VersionAction, banner=version_banner('planetgen.cli.orbits'))
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
            print("No previous orbit update found for this database -- moving what is due since each "
                  "object was made.")
        elif elapsed_years < 0:
            # The server's clock went back since the last run: never move
            # anything backwards; this run restarts the clock from now.
            print(f"The last orbit update is {-elapsed_years:.6f} years in the future (the database server's "
                  "clock moved back?) -- moving nothing and restarting the clock from now.")
        else:
            print(f"{elapsed_years:.6f} years elapsed since the last update -- advancing what is due.")
        clock = orbit_clock(conn)

        # A bar over the steps, with the tables (or sectors) of the step
        # under way beneath it; the Generate page's job status draws the
        # same two bars.
        with StageProgress(len(STAGES)) as bar:
            bar.stage(_stage_label(0))
            counts = advance_orbital_phases(conn, clock, on_progress=bar.detail)
            # A comet's position needs Kepler's or Barker's equation, so it
            # moves in Python rather than in advance_orbital_phases' SQL.
            bar.stage(_stage_label(1))
            counts["comets"] = advance_comet_orbits(conn, clock)
            counts["facilities"] = advance_facility_orbits(conn, clock)

            # Galactic motion: move what is due along its galactic orbit,
            # refile what drifted into another sector, then recompute what
            # depends on position (containment, octants, nearest systems,
            # location text).
            bar.stage(_stage_label(2))
            motion = advance_galactic_positions(conn, clock, on_progress=bar.detail)
            counts.update(motion["counts"])
            finish_orbit_update(conn, clock)
            bar.stage(_stage_label(3))
            locations = refresh_after_motion(conn, motion["sectors"], on_progress=bar.detail)
            conn.commit()
            # Paths of every sector holding a star system, rogue planet or comet (they all begin
            # where their body is now, which just moved), and of the ones that emptied.
            bar.stage(_stage_label(4))
            holding = conn.execute(
                "SELECT sector_id FROM star_systems UNION SELECT sector_id FROM rogue_planets"
                " UNION SELECT sector_id FROM interstellar_comets").fetchall()
            sectors = sorted({row["sector_id"] for row in holding}.union(motion["sectors"]) - {None})
            for done, sector_id in enumerate(sectors, start=1):
                compute_sector_paths(conn, sector_id)
                bar.detail("sector paths", done, len(sectors))
            conn.commit()
        summary = ", ".join(f"{counts[table]} {label}" for table, label in TABLE_LABELS.items())
        print(f"Updated: {summary} (only what was due).")
        print(f"Moved {motion['moved']} object(s) along their galactic orbits; {motion['refiled']} changed "
              f"sector. Refreshed nearest systems and containment; rewrote {locations} location(s).")
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
    finally:
        conn.close()


if __name__ == "__main__":
    main()
