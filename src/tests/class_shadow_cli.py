#!/usr/bin/env python3
"""
src/tests/class_shadow_cli.py -- the shadow sweep of the condition-based classes (GEN.645, C2 of
docs/design/condition-based-planet-classes.md). Not a pytest module (no `test_*` name).

Generates systems in memory, classifies every planet and moon from its stored conditions
(`physics.planet_classes.classify`) and compares the letter with the one the generator drew. Prints (and with
`--report` writes) a Markdown report: how many bodies change letter, the confusion table drawn letter against
derived letter, the share of each letter before and after, and for every class the rule that is least like its
drawn bodies. The thresholds in `physics/data/planet_classes.json` are tuned until the changes are the fewest
the rules allow.

Usage:
    python src/tests/class_shadow_cli.py
    python src/tests/class_shadow_cli.py --systems 2000 --seed 7 --report /mnt/project-files/reports/class-shadow-sweep.md
"""

import argparse
import collections
import sys

from planetgen import tuning
from planetgen.generation import plausibility
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import planet_classes
from planetgen.util import draw


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


def sample_forced(per_class, seed):
    """`(body, is_moon)` for `per_class` planets (and their moons) of every class in every zone it allows, on host
    stars of every main-sequence type: the rare classes the natural mix hardly draws."""
    draw.set_run_seed(seed)
    bodies = []
    for letter, data in tuning.PLANET_CLASSES.items():
        for zone in plausibility.ZONE_CHARS:
            if not data[zone]:
                continue
            for _ in range(per_class):
                config = SystemConfig()
                config.STAR_TYPE = draw.choice(plausibility.HOST_STAR_TYPES)
                star = Star(config)
                planet = Planet(config, star, star.habitable_zone, plausibility.distance_for_zone(star, zone),
                                planet_class=letter)
                bodies.append((planet, False))
                bodies.extend((moon, True) for moon in planet.moons)
    return bodies


def compare(bodies):
    """`Counter` of `(drawn letter, derived letter)` over `bodies`, and the list of those two letters per body."""
    pairs = collections.Counter()
    for body, _ in bodies:
        pairs[(body.planet_class, planet_classes.classify(body))] += 1
    return pairs


def format_report(pairs, systems, seed, per_class=0):
    total = sum(pairs.values())
    same = sum(n for (drawn, derived), n in pairs.items() if drawn == derived)
    drawn_letters = sorted({d for d, _ in pairs})
    derived_letters = sorted({d for _, d in pairs})
    drawn_total = collections.Counter()
    derived_total = collections.Counter()
    for (drawn, derived), n in pairs.items():
        drawn_total[drawn] += n
        derived_total[derived] += n
    lines = [
        "# Condition-based classes: shadow sweep", "",
        f"Rule file version {planet_classes.VERSION}; "
        + (f"{per_class} planets of every class and zone it allows (plus their moons)" if per_class
           else f"{systems} generated systems") + f" from run seed {seed}; {total} planets and moons.", "",
        f"**{same} of {total} bodies ({100.0 * same / total:.1f} percent) keep their letter;** "
        f"{total - same} change.", "",
        "## Per drawn letter: how many keep it", "",
        "| Drawn | Bodies | Keep | Share kept | Becomes (most common) |", "|---|---|---|---|---|",
    ]
    for letter in drawn_letters:
        kept = pairs.get((letter, letter), 0)
        moved = sorted(((n, d) for (o, d), n in pairs.items() if o == letter and d != letter), reverse=True)
        becomes = ", ".join(f"{d} {n}" for n, d in moved[:4]) or "-"
        lines.append(f"| {letter} | {drawn_total[letter]} | {kept} | {100.0 * kept / drawn_total[letter]:.0f}% | {becomes} |")
    lines += ["", "## Share of each letter, drawn against derived", "", "| Letter | Drawn | Derived | Drawn share | Derived share |",
              "|---|---|---|---|---|"]
    for letter in sorted(set(drawn_letters) | set(derived_letters)):
        lines.append(f"| {letter} | {drawn_total[letter]} | {derived_total[letter]} | "
                     f"{100.0 * drawn_total[letter] / total:.2f}% | {100.0 * derived_total[letter] / total:.2f}% |")
    lines += ["", "## Confusion table (rows drawn, columns derived)", "",
              "| Drawn \\ Derived | " + " | ".join(derived_letters) + " |", "|---" * (len(derived_letters) + 1) + "|"]
    for letter in drawn_letters:
        lines.append(f"| {letter} | " + " | ".join(str(pairs.get((letter, d), 0) or "") for d in derived_letters) + " |")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--systems", type=int, default=800)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--habitable-share", type=float, default=0.0,
                        help="share of the systems generated with a forced habitable world (default 0)")
    parser.add_argument("--per-class", type=int, default=0,
                        help="instead of whole systems, this many planets of every class and zone it allows")
    parser.add_argument("--report", help="write the Markdown report to this file")
    args = parser.parse_args(argv)
    bodies = (sample_forced(args.per_class, args.seed) if args.per_class
              else sample_bodies(args.systems, args.seed, args.habitable_share))
    report = format_report(compare(bodies), args.systems, args.seed, args.per_class)
    print(report)
    if args.report:
        with open(args.report, "w", encoding="utf-8") as handle:
            handle.write(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
