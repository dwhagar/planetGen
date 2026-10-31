#!/usr/bin/env python3
"""
src/tests/class_sweep_cli.py -- the class sweep of the condition-based classes (GEN.645, GEN.646). Not a pytest module (no `test_*` name).

Generates systems in memory and prints (and with `--report` writes) a Markdown report of the classes the generator derived: the share of
every letter among planets and among moons, the mass bins of the planets, and how many bodies change letter when `classify` is run
again on their stored conditions (none should). The distribution is what the thresholds in `physics/data/planet_classes.json` and the
mass law are tuned against.

Usage:
    python src/tests/class_sweep_cli.py
    python src/tests/class_sweep_cli.py --systems 2000 --seed 7 --report /mnt/project-files/reports/class-sweep.md
"""

import argparse
import collections
import sys

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, planet_classes
from planetgen.util import draw

MASS_BINS = ((0, 0.1), (0.1, 0.5), (0.5, 2), (2, 10), (10, 35), (35, 100), (100, 1e12))
"""tuple: The edges of the planet mass bins, in Earth masses."""


def sample_bodies(systems, seed, habitable_share=0.0):
    """`(body, is_moon)` for every planet and moon of `systems` generated systems from the run seed `seed`; the
    first `habitable_share` of the systems are asked for a habitable world (so the rare temperate classes are seen)."""
    draw.set_run_seed(seed)
    bodies = []
    for index in range(systems):
        config = SystemConfig()
        if index < habitable_share * systems:
            config.HABITABLE_WORLD = True
        system = StarSystem(system_config=config)
        for planet in system.planets:
            if not hasattr(planet, "planet_class"):
                continue  # an asteroid belt
            bodies.append((planet, False))
            bodies.extend((moon, True) for moon in planet.moons)
    return bodies


ASTEROID_LABEL = "Asteroid"


def changed_letters(bodies):
    """The bodies whose stored letter is not the one `classify` gives their stored conditions."""
    return [body for body, _ in bodies if planet_classes.classify(body) != body.planet_class]


def shares(bodies):
    """`Counter` of the letters of `bodies`; a body under a roundness ratio of 1 has no letter and is counted as `ASTEROID_LABEL` (GEN.1160)."""
    return collections.Counter(body.planet_class or ASTEROID_LABEL for body, _ in bodies)


def mass_bins(bodies):
    """`Counter` of the mass bin (index into `MASS_BINS`) of every planet in `bodies`."""
    counts = collections.Counter()
    for body, is_moon in bodies:
        if is_moon:
            continue
        earth = body.mass / constants.EARTH_MASS_TO_KG
        counts[next(i for i, (low, high) in enumerate(MASS_BINS) if low <= earth < high)] += 1
    return counts


def format_report(bodies, systems, seed):
    planets = [entry for entry in bodies if not entry[1]]
    moons = [entry for entry in bodies if entry[1]]
    planet_counts, moon_counts = shares(planets), shares(moons)
    letters = sorted(set(planet_counts) | set(moon_counts))
    lines = [
        "# Condition-based classes: sweep", "",
        f"Rule file version {planet_classes.VERSION}; {systems} generated systems from run seed {seed}; {len(planets)} planets and "
        f"{len(moons)} moons ({len(moons) / max(1, len(planets)):.2f} moons per planet, {len(planets) / max(1, systems):.2f} planets per system).", "",
        f"{len(changed_letters(bodies))} bodies change letter when classified again from their stored conditions.", "",
        "## Share of each letter", "", "| Letter | Planets | Planet share | Moons | Moon share |", "|---|---|---|---|---|",
    ]
    for letter in letters:
        lines.append(f"| {letter} | {planet_counts[letter]} | {100.0 * planet_counts[letter] / max(1, len(planets)):.2f}% | "
                     f"{moon_counts[letter]} | {100.0 * moon_counts[letter] / max(1, len(moons)):.2f}% |")
    bins = mass_bins(bodies)
    lines += ["", "## Planets by mass (Earth masses)", "", "| Mass | Planets per system |", "|---|---|"]
    for index, (low, high) in enumerate(MASS_BINS):
        lines.append(f"| {low} to {high if high < 1e9 else 'more'} | {bins[index] / max(1, systems):.3f} |")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--systems", type=int, default=800)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--habitable-share", type=float, default=0.0,
                        help="share of the systems generated with a forced habitable world (default 0)")
    parser.add_argument("--report", help="write the Markdown report to this file")
    args = parser.parse_args(argv)
    bodies = sample_bodies(args.systems, args.seed, args.habitable_share)
    report = format_report(bodies, args.systems, args.seed)
    print(report)
    if args.report:
        with open(args.report, "w", encoding="utf-8") as handle:
            handle.write(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
