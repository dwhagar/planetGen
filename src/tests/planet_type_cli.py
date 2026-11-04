#!/usr/bin/env python3
"""
src/tests/planet_type_cli.py -- the shadow report of the body types (GEN.1219). Not a pytest module (no `test_*` name). Read-only.

Types a sample of generated systems with `physics/planet_type.py` and prints: the type against the stored body_type (t or g) and the stored class letter,
the bodies within 20 percent of each cut, and the bodies whose radius disagrees with their type (a terrestrial body over 1.6 Earth radii with under 10
percent water is a rocky super-Earth; with 10 percent or more water it is a water-rich world).

Usage:
    python src/tests/planet_type_cli.py --systems 300 --seed 5 --habitable-share 0.25 --report /mnt/project-files/reports/planet-type-shadow.md
"""

import argparse
import collections
import sys

from planetgen.physics import planet_type as pt
from tests import class_sweep_cli

EARTH_RADIUS_KM = 6371.0


def rows_of(bodies):
    """One record per `(body, is_moon)`."""
    rows = []
    for body, is_moon in bodies:
        envelope = getattr(body, "envelope_mass_fraction", None)
        kind, known = pt.type_of(envelope, getattr(body, "body_type", None))
        rows.append({"type": kind, "known": known, "envelope": envelope, "body_type": getattr(body, "body_type", None), "class": body.planet_class,
                     "moon": is_moon, "radius": (getattr(body, "radius", None) or 0.0) / EARTH_RADIUS_KM,
                     "water": getattr(body, "water_mass_fraction", None) or 0.0,
                     "mass": (getattr(body, "mass", None) or 0.0) / 5.972e24})
    return rows


def format_report(rows, systems, seed, share):
    n = len(rows)
    lines = [f"# Body types shadow, {systems} systems from seed {seed} (habitable share {share})", "",
             f"{n} planets and moons; terrestrial under {pt.TERRESTRIAL_MAX_ENVELOPE:.1%} envelope, gas from {pt.GAS_MIN_ENVELOPE:.0%}.", "",
             "## Type against the stored body_type", "", "| type | body_type t | body_type g | envelope not stored | total |", "|---|---|---|---|---|"]
    for kind in pt.TYPES:
        sel = [r for r in rows if r["type"] == kind]
        lines.append(f"| {kind} | {sum(1 for r in sel if r['body_type'] == 't')} | {sum(1 for r in sel if r['body_type'] == 'g')} | {sum(1 for r in sel if not r['known'])} | {len(sel)} |")
    lines += ["", "## Type against the stored class letter", "", "| type | classes (count) |", "|---|---|"]
    for kind in pt.TYPES:
        c = collections.Counter(r["class"] for r in rows if r["type"] == kind)
        lines.append(f"| {kind} | " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())) + " |")
    lines += ["", "## Near a cut (within 20 percent of the line)", ""]
    for name, cut in (("terrestrial/transitional", pt.TERRESTRIAL_MAX_ENVELOPE), ("transitional/gas", pt.GAS_MIN_ENVELOPE)):
        near = [r for r in rows if r["envelope"] is not None and 0.8 * cut <= r["envelope"] <= 1.2 * cut]
        lines.append(f"- {name}: {len(near)} bodies" + (f" (classes {dict(collections.Counter(r['class'] for r in near))})" if near else ""))
    big = [r for r in rows if r["type"] == pt.TERRESTRIAL and not r["moon"] and r["radius"] > 1.6]
    wet = [r for r in big if r["water"] >= 0.1]
    lines += ["", "## Radius against type", "",
              f"- terrestrial planets over 1.6 Earth radii: {len(big)}; with 10 percent water or more: {len(wet)}; dry (rocky super-Earths up to 20 Earth masses): {len(big) - len(wet)}",
              f"- terrestrial-type bodies of class T (a hydrogen air but no stored envelope): {sum(1 for r in rows if r['type'] == pt.TERRESTRIAL and r['class'] == 'T')}",
              f"- gas-type bodies under 35 Earth masses: {sum(1 for r in rows if r['type'] == pt.GAS and r['mass'] < 35)}; transitional bodies over 35 Earth masses: {sum(1 for r in rows if r['type'] == pt.TRANSITIONAL and r['mass'] >= 35)}"]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--systems", type=int, default=300)
    parser.add_argument("--seed", type=int, default=5)
    parser.add_argument("--habitable-share", type=float, default=0.25)
    parser.add_argument("--report")
    args = parser.parse_args(argv)
    rows = rows_of(class_sweep_cli.sample_bodies(args.systems, args.seed, args.habitable_share))
    text = format_report(rows, args.systems, args.seed, args.habitable_share)
    print(text)
    if args.report:
        with open(args.report, "w", encoding="utf-8") as handle:
            handle.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
