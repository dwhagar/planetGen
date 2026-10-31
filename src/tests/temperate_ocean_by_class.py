#!/usr/bin/env python3
"""
src/tests/temperate_ocean_by_class.py -- for each class, how many bodies are temperate (250 to 330 K) with an ocean (GEN.1135, GEN.1136,
GEN.1142, asked for by Boss after the train 6 class epoch). Not a pytest module (no `test_*` name). Changes nothing; prints a table.

An ocean is counted two ways, in two columns: a surface ocean (the water is liquid at the surface: hydrosphere "surface ocean", "hycean" or "brine seep") and
an ice-covered ocean (liquid water under an ice shell, hydrosphere "ice-covered ocean"). The first version of this table counted both as one, so the glaciated
worlds of class P (ice over an ocean) showed as wet; `scripts/class_counts.py` counts the surface ocean only, and the two now agree on that column.

Uses the sampling of `class_sweep_cli` (a share of the systems is asked for a habitable world so the rare temperate classes are seen).

Usage:
    python src/tests/temperate_ocean_by_class.py
    python src/tests/temperate_ocean_by_class.py --systems 600 --seed 7 --habitable-share 0.25
"""

import argparse
import collections
import sys

from planetgen.physics import planet_classes
from tests import class_sweep_cli

SURFACE_OCEAN = ("surface ocean", "hycean", "brine seep")
"""tuple: The hydrospheres with liquid water at the surface (the same set `scripts/class_counts.py` counts as a surface ocean)."""

TEMPERATE_K = (250.0, 330.0)
"""tuple: The mean surface temperatures, K, that count as temperate here (Boss, 6:26 PM EDT 2026-10-30)."""


def temperate_ocean_counts(bodies):
    """`(all, temperate, surface_ocean, ice_covered)`: four `Counter`s of the class letters of `bodies` (a list of `(body, is_moon)`): every body, the temperate ones,
    the temperate ones with a surface ocean and the temperate ones with an ocean under ice."""
    every, temperate, surface, covered = (collections.Counter() for _ in range(4))
    for body, _is_moon in bodies:
        quantities = planet_classes.quantities(planet_classes.values_from_body(body))
        letter = body.planet_class
        every[letter] += 1
        temperature = quantities["temperature_k"]
        if temperature is not None and TEMPERATE_K[0] <= temperature <= TEMPERATE_K[1]:
            temperate[letter] += 1
            if quantities["hydrosphere"] in SURFACE_OCEAN:
                surface[letter] += 1
            elif quantities["hydrosphere"] == "ice-covered ocean":
                covered[letter] += 1
    return every, temperate, surface, covered


def format_table(bodies, systems, seed):
    every, temperate, surface, covered = temperate_ocean_counts(bodies)
    lines = [f"{systems} generated systems from run seed {seed}: {len(bodies)} bodies; temperate means {TEMPERATE_K[0]:.0f} to {TEMPERATE_K[1]:.0f} K.", "",
             "| Class | Name | Bodies | Temperate | Temperate with a surface ocean | Temperate with an ice-covered ocean |", "|---|---|---|---|---|---|"]
    for letter in sorted(every):
        name = next((c.name for c in planet_classes.CLASSES if c.letter == letter), "")
        lines.append(f"| {letter} | {name} | {every[letter]} | {temperate[letter]} | {surface[letter]} | {covered[letter]} |")
    lines.append(f"| all | | {sum(every.values())} | {sum(temperate.values())} | {sum(surface.values())} | {sum(covered.values())} |")
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--systems", type=int, default=300)
    parser.add_argument("--seed", type=int, default=5)
    parser.add_argument("--habitable-share", type=float, default=0.25)
    args = parser.parse_args(argv)
    bodies = class_sweep_cli.sample_bodies(args.systems, args.seed, args.habitable_share)
    print(format_table(bodies, args.systems, args.seed))
    return 0


if __name__ == "__main__":
    sys.exit(main())
