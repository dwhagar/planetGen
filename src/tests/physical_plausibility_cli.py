#!/usr/bin/env python3
"""
src/tests/physical_plausibility_cli.py -- physical-plausibility anomaly
finder (CLI). Not a pytest test module itself (no `test_*` name, so pytest
won't collect it) -- it lives alongside `test_physical_plausibility.py`
because it shares that file's engine (`planetgen/generation/plausibility.py`)
and exists specifically to batch-run it for human-reviewed findings,
rather than as an automated pass/fail check.

Batch-generates planets (and their moons) in memory -- no database
round-trip needed -- across every zone and a
broad grid of main-sequence host star spectral types, then reports:

  * hard physical-invariant violations (analytically-derived surface
    gravity bounds per class, plus basic finiteness/sign checks on
    temperature and pressure) -- these are unambiguous bugs;
  * statistical outliers (Tukey's fences on each class's own generated
    distribution) for surface_temperature, atmospheric_pressure, gravity,
    and scale_height -- flagged for human review, not asserted to be zero.

See `planetgen/generation/plausibility.py` (same directory's parent) for the
full design rationale (why two tiers, why host stars are sampled across
spectral types, why this isn't hand-authored per-class numeric bounds)
and docs/TODO.md's "Physical-plausibility test suite (anomaly finder)"
future idea for the original ask.

Usage:
    python src/tests/physical_plausibility_cli.py
    python src/tests/physical_plausibility_cli.py --n 300 --no-moons
    python src/tests/physical_plausibility_cli.py --classes M P --n 500
"""

import argparse
import logging
import sys


from planetgen.generation import plausibility
from planetgen import tuning as prog_c

# Suppress transformers warnings pulled in transitively via planetgen.
logging.getLogger("transformers").setLevel(logging.ERROR)


def process_args():
    parser = argparse.ArgumentParser(description="Physical-plausibility anomaly finder for planetGen.")
    parser.add_argument(
        "--n", type=int, default=150,
        help="Number of independent planets to generate per zone (default: 150).",
    )
    parser.add_argument(
        "--classes", nargs="+", default=None, metavar="CLASS",
        help="Draw each sample again until it is one of these derived classes (default: every class as it comes).",
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

    wanted = [c.upper() for c in args.classes] if args.classes else [None]
    unknown = {c for c in wanted if c is not None} - set(prog_c.PLANET_CLASSES)
    if unknown:
        print(f"Unknown planet class(es): {sorted(unknown)}", file=sys.stderr)
        return 2

    zones = plausibility.ZONE_CHARS
    print(
        f"Generating ~{args.n * len(zones) * len(wanted)}+ bodies "
        f"({args.n} per zone x {len(zones)} zones"
        f"{'' if wanted == [None] else f' x classes {wanted}'}, plus moons"
        f"{'' if not args.no_moons else ' [disabled]'}) ...",
        file=sys.stderr,
    )

    records = []
    for wanted_class in wanted:
        for zone in zones:
            records.extend(plausibility.generate_sample(zone, args.n, include_moons=not args.no_moons, wanted_class=wanted_class))

    report = plausibility.analyze(records, k=args.k)
    print(plausibility.format_report(report, show_outlier_examples=not args.no_examples))

    total_hard_violations = sum(len(v["hard_violations"]) for v in report.values())
    return 1 if total_hard_violations else 0


if __name__ == "__main__":
    sys.exit(main())
