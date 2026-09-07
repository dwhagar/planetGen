#!/usr/bin/env python
"""
physicalPlausibility.py -- physical-plausibility anomaly finder (CLI).

Batch-generates planets (and their moons) in memory -- no database
round-trip needed -- across every valid (planet class, zone) pair and a
broad grid of main-sequence host star spectral types, then reports:

  * hard physical-invariant violations (analytically-derived surface
    gravity bounds per class, plus basic finiteness/sign checks on
    temperature and pressure) -- these are unambiguous bugs;
  * statistical outliers (Tukey's fences on each class's own generated
    distribution) for surface_temperature, atmospheric_pressure, gravity,
    and scale_height -- flagged for human review, not asserted to be zero.

See `stellarObjects/plausibility.py` for the full design rationale (why two
tiers, why host stars are sampled across spectral types, why this isn't
hand-authored per-class numeric bounds) and TODO.md's "Physical-plausibility
test suite (anomaly finder)" future idea for the original ask.

Usage:
    python physicalPlausibility.py
    python physicalPlausibility.py --n 300 --no-moons
    python physicalPlausibility.py --classes M P --n 500
"""

import argparse
import logging
import sys

from stellarObjects import plausibility
from stellarObjects import program_constants as prog_c

# Suppress transformers warnings pulled in transitively via stellarObjects.
logging.getLogger("transformers").setLevel(logging.ERROR)


def process_args():
    parser = argparse.ArgumentParser(description="Physical-plausibility anomaly finder for planetGen.")
    parser.add_argument(
        "--n", type=int, default=150,
        help="Number of independent planets to generate per (class, zone) pair (default: 150).",
    )
    parser.add_argument(
        "--classes", nargs="+", default=None, metavar="CLASS",
        help="Restrict to these planet class codes (default: all classes with at least one valid zone).",
    )
    parser.add_argument(
        "--no-moons", action="store_true",
        help="Skip moon generation/analysis (faster; moons are included by default).",
    )
    parser.add_argument(
        "--k", type=float, default=3.0,
        help="Tukey's-fences multiplier for statistical outlier flagging (default: 3.0, the conventional "
             "'far outlier' threshold).",
    )
    parser.add_argument(
        "--no-examples", action="store_true",
        help="Omit individual outlier example lines from the report (summary stats only).",
    )
    return parser.parse_args()


def main():
    args = process_args()

    pairs = plausibility.VALID_CLASS_ZONE_PAIRS
    if args.classes:
        wanted = {c.upper() for c in args.classes}
        unknown = wanted - set(prog_c.PLANET_CLASSES)
        if unknown:
            print(f"Unknown planet class(es): {sorted(unknown)}", file=sys.stderr)
            return 2
        pairs = [(c, z) for c, z in pairs if c in wanted]
        if not pairs:
            print(f"None of {sorted(wanted)} have a valid zone to generate in.", file=sys.stderr)
            return 2

    total_bodies_estimate = args.n * len(pairs)
    print(
        f"Generating ~{total_bodies_estimate}+ bodies "
        f"({args.n} per (class, zone) pair x {len(pairs)} pairs, plus moons"
        f"{'' if not args.no_moons else ' [disabled]'}) ...",
        file=sys.stderr,
    )

    records = []
    for cls, zone in pairs:
        records.extend(plausibility.generate_sample(cls, zone, args.n, include_moons=not args.no_moons))

    report = plausibility.analyze(records, k=args.k)
    print(plausibility.format_report(report, show_outlier_examples=not args.no_examples))

    total_hard_violations = sum(len(v["hard_violations"]) for v in report.values())
    return 1 if total_hard_violations else 0


if __name__ == "__main__":
    sys.exit(main())
