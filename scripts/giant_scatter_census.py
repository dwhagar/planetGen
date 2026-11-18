#!/usr/bin/env python3
"""
Giants, outer disks and comet reservoirs of a seeded sample of systems (GEN.1788): prints only, nothing is stored.

    python scripts/giant_scatter_census.py [--systems 300] [--seed 1000] [--jobs 4] [--star-type G2V] [--out rows.json]

Generates `--systems` systems from the seeds `--seed` upward (the star forced to `--star-type` when given, for a Sun-like sample) and reads, per system, its giants, its outer-disk belts
with their scattered share (`scatter.scattered_share`), whether the comet supply finds a reservoir (`comet_supply.family_weights`) and how many comets the system holds. Prints the rates:
a giant at any distance, at 2 to 8 AU and at 8 to 32 AU (observed 14.1 and 8.9 planets per 100 stars, Fulton et al. 2021), an outer disk, a giant with an outer disk, and the systems
with comets. Re-run it on a new main, and before and after any change to `physics/scatter.py` or `physics/belt_origin.py`.
"""

import argparse
import json
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import comet_supply, scatter
from planetgen.util import draw

EARTH_KG = 5.972e24
SOLAR_KG = 1.989e30


def measure(job):
    """One system's giants `(au, Earth masses, formation au)`, outer disks and comet reservoir."""
    seed, star_type = job
    config = SystemConfig()
    config.PLANETS = True
    if star_type:
        config.STAR_TYPE = star_type
    with draw.bound(seed):
        system = StarSystem(system_config=config)
    giants = sorted((p for p in system.all_planets if p.body_type == "g"), key=lambda p: p.distance)
    belts = [b for b in system.all_planets if b.body_type == "a"]
    disks = []
    for belt in belts:
        if getattr(belt, "origin", None) != "outer-disk":
            continue
        fates = getattr(belt, "scatter_fates", None)
        disks.append({"mass_earth": getattr(belt, "primordial_mass_earth", 0) or 0.0, "scattered": scatter.scattered_share(fates) if fates else None})
    return {"seed": seed, "star_msun": system.star.mass / SOLAR_KG,
            "giants": [(p.distance, p.mass / EARTH_KG, p.formation_distance_au or p.distance) for p in giants],
            "outer_disks": disks, "comets": len(system.comets) + len(getattr(system, "secondary_comets", None) or []),
            "reservoir": comet_supply.family_weights(belts) is not None}


def report(rows):
    """The rates of `rows`, as lines."""
    total = len(rows)
    with_giant = [row for row in rows if row["giants"]]
    both = [row for row in with_giant if row["outer_disks"]]

    def share(count):
        return f"{count} ({100.0 * count / total:.1f} percent)"

    def in_range(low, high):
        return sum(1 for row in rows if any(low <= giant[0] < high for giant in row["giants"]))
    return [f"systems {total}", f"any giant {share(len(with_giant))}", f"a giant at 2 to 8 AU {share(in_range(2.0, 8.0))} (observed 14.1 per 100 stars)",
            f"a giant at 8 to 32 AU {share(in_range(8.0, 32.0))} (observed 8.9 per 100 stars)",
            f"an outer disk {share(sum(1 for row in rows if row['outer_disks']))}", f"a giant and an outer disk {share(len(both))}",
            f"  of those, a scattered share above 0: {sum(1 for row in both if any(disk['scattered'] for disk in row['outer_disks']))}",
            f"systems with comets {share(sum(1 for row in rows if row['comets']))}", f"systems with a comet reservoir {share(sum(1 for row in rows if row['reservoir']))}"]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--systems", type=int, default=300)
    parser.add_argument("--seed", type=int, default=1000)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--star-type", default=None, help="force the star, e.g. G2V")
    parser.add_argument("--out", default=None, help="write the rows as JSON")
    args = parser.parse_args()
    with Pool(args.jobs) as pool:
        rows = pool.map(measure, [(seed, args.star_type) for seed in range(args.seed, args.seed + args.systems)], chunksize=1)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as handle:
            json.dump(rows, handle)
    print("\n".join(report(rows)))


if __name__ == "__main__":
    main()
