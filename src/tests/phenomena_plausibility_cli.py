#!/usr/bin/env python
"""
src/tests/phenomena_plausibility_cli.py -- exotic-phenomena plausibility
anomaly finder (CLI). Not a pytest test module itself (no `test_*` name,
so pytest won't collect it) -- it lives alongside
`test_phenomena_plausibility.py` because it shares that file's engine
(`stellarObjects/phenomenaPlausibility.py`) and exists specifically to
batch-run it for human-reviewed findings, rather than as an automated
pass/fail check. Mirrors `physical_plausibility_cli.py`'s own role for
the planet-focused engine.

Batch-generates standalone phenomena (black holes, neutron stars, nebulae,
supernova remnants, rogue planets, interstellar comets, asteroid fields)
in memory -- no database round-trip needed -- across all seven types,
then reports:

  * hard physical-invariant violations (analytically-derived formulas --
    e.g. a black hole's event horizon radius must match the Schwarzschild
    formula for its own mass) -- these are unambiguous bugs;
  * statistical outliers (Tukey's fences on each type's own generated
    distribution) for its continuous metrics, plus observed-vs-expected
    category-frequency comparisons (e.g. black hole accretion-disk
    proportion vs. its configured chance) -- flagged for human review,
    not asserted to be exactly the configured value.

See `stellarObjects/phenomenaPlausibility.py` for the full design
rationale.

This file lives at src/tests/, two levels under src/ where
`stellarObjects` actually lives -- one dirname() up from this file's own
directory reaches src/, added to sys.path below.

Usage:
    python src/tests/phenomena_plausibility_cli.py
    python src/tests/phenomena_plausibility_cli.py --n 500
    python src/tests/phenomena_plausibility_cli.py --types black-hole neutron-star
"""

import argparse
import logging
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from stellarObjects import phenomenaPlausibility as pp

# Suppress transformers warnings pulled in transitively via stellarObjects.
logging.getLogger("transformers").setLevel(logging.ERROR)


def process_args():
    parser = argparse.ArgumentParser(description="Exotic-phenomena plausibility anomaly finder for planetGen.")
    parser.add_argument(
        "--n", type=int, default=200,
        help="Number of independent bodies to generate per phenomenon type (default: 200).",
    )
    parser.add_argument(
        "--types", nargs="+", default=None, metavar="TYPE", choices=list(pp.PHENOMENON_TYPES),
        help="Restrict to these phenomenon types (default: all seven).",
    )
    parser.add_argument(
        "--k", type=float, default=3.0,
        help="Tukey's-fences multiplier for statistical outlier flagging (default: 3.0, the conventional "
             "'far outlier' threshold).",
    )
    parser.add_argument(
        "--category-margin", type=float, default=0.1,
        help="A category's observed proportion is flagged when it differs from its expected proportion by "
             "more than this (default: 0.1, i.e. 10 percentage points).",
    )
    parser.add_argument(
        "--no-examples", action="store_true",
        help="Omit individual outlier example lines from the report (summary stats only).",
    )
    return parser.parse_args()


def main():
    args = process_args()
    types = args.types or list(pp.PHENOMENON_TYPES)

    print(f"Generating {args.n} bodies per phenomenon type x {len(types)} types ...", file=sys.stderr)

    records = []
    for phenomenon_type in types:
        records.extend(pp.generate_sample(phenomenon_type, args.n))

    report = pp.analyze(records, k=args.k)
    print(pp.format_report(report, show_outlier_examples=not args.no_examples, category_flag_margin=args.category_margin))

    total_hard_violations = sum(len(v["hard_violations"]) for v in report.values())
    return 1 if total_hard_violations else 0


if __name__ == "__main__":
    sys.exit(main())
