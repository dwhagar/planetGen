#!/usr/bin/env python3
"""
Class counts of a seeded sample of systems (GEN.1143): prints only, nothing is stored.

    python scripts/class_counts.py [--systems 300] [--seed 1000] [--jobs 4]

Generates `--systems` systems from the seeds `--seed` upward with the default settings and reads each planet and moon's class from the finished body
(`planet_classes.classify`). Prints:

1. the bodies of each class;
2. per class, the bodies between 250 and 330 K whose water is a liquid ocean (a surface ocean; a hycean ocean is counted beside it);
3. the hydrogen-rich N bodies under 15 Earth masses (hydrogen at least 30 percent of the air; helium is not a stored species), with their masses;
4. how many bodies reach E, H and L, and how many would match each row's own conditions if the rows before it did not take them.

The numbers go in the class reports (reports/class-index-rebuild-assessment-2026-10-30.md); re-run it on each new main.
"""

import argparse
import collections
import os
import statistics
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from planetgen.generation.config import SystemConfig  # noqa: E402
from planetgen.generation.system import StarSystem  # noqa: E402
from planetgen.physics import planet_classes  # noqa: E402
from planetgen.util import draw  # noqa: E402

TEMPERATE_K = (250.0, 330.0)
HYDROGEN_RICH_SHARE = 0.3
HYDROGEN_RICH_MAX_EARTH = 15.0
WATCHED = ("E", "H", "L")


def measure(seed):
    """The class letter, the water state, the temperature, the quantities that matter and the rows each body would match, per body of one system."""
    config = SystemConfig()
    config.PLANETS = True
    with draw.bound(seed):
        system = StarSystem(system_config=config)
    rows = []
    for planet in system.all_planets:
        if planet.body_type == "a":
            continue
        for body in (planet, *planet.moons):
            values = planet_classes.values_from_body(body)
            found = planet_classes.quantities(values)
            own = {c.letter for c in planet_classes.CLASSES if c.kind == found["body_type"] and c.when and c.matches(found)}
            rows.append({"letter": planet_classes.classify_values(values).letter, "hydrosphere": found["hydrosphere"],
                         "temperature_k": found["temperature_k"], "mass_earth": found["mass_earth"],
                         "hydrogen": found["share_h2"], "own": sorted(own & set(WATCHED))})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--systems", type=int, default=300)
    parser.add_argument("--seed", type=int, default=1000)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()

    with Pool(args.jobs) as pool:
        bodies = [b for rows in pool.map(measure, range(args.seed, args.seed + args.systems), chunksize=5) for b in rows]
    print(f"{args.systems} systems from seed {args.seed}: {len(bodies)} planets and moons\n")

    by_class = collections.Counter(b["letter"] for b in bodies)
    low, high = TEMPERATE_K
    print(f"{'class':<6}{'bodies':>8}{f'{low:g} to {high:g} K, surface ocean':>34}{'hycean':>8}")
    for letter in planet_classes.LETTERS:
        temperate = [b for b in bodies if b["letter"] == letter and b["temperature_k"] is not None and low <= b["temperature_k"] <= high]
        print(f"{letter:<6}{by_class[letter]:>8}{sum(b['hydrosphere'] == 'surface ocean' for b in temperate):>34}"
              f"{sum(b['hydrosphere'] == 'hycean' for b in temperate):>8}")
    total_ocean = sum(b["hydrosphere"] == "surface ocean" and b["temperature_k"] is not None and low <= b["temperature_k"] <= high for b in bodies)
    print(f"{'all':<6}{len(bodies):>8}{total_ocean:>34}\n")

    rich = [b for b in bodies if b["letter"] == "N" and b["hydrogen"] >= HYDROGEN_RICH_SHARE and b["mass_earth"] < HYDROGEN_RICH_MAX_EARTH]
    nines = by_class["N"]
    masses = [b["mass_earth"] for b in rich]
    print(f"N bodies: {nines}; hydrogen-rich (hydrogen at least {HYDROGEN_RICH_SHARE:.0%} of the air) and under {HYDROGEN_RICH_MAX_EARTH:g} Earth masses: {len(rich)}"
          + (f", median {statistics.median(masses):.2f}, range {min(masses):.2f} to {max(masses):.2f} Earth masses" if rich else ""))
    print(f"hydrogen-rich bodies of any class under {HYDROGEN_RICH_MAX_EARTH:g} Earth masses: "
          f"{sum(b['hydrogen'] >= HYDROGEN_RICH_SHARE and b['mass_earth'] < HYDROGEN_RICH_MAX_EARTH for b in bodies)}\n")

    print(f"{'class':<6}{'reached':>9}{'would match its own conditions':>34}")
    for letter in WATCHED:
        print(f"{letter:<6}{by_class[letter]:>9}{sum(letter in b['own'] for b in bodies):>34}")


if __name__ == "__main__":
    main()
