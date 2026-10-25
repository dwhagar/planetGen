#!/usr/bin/env python3
"""
src/tests/climate_tuning_cli.py -- interactive per-class climate tuning tool.

Not a pytest test module itself (no `test_*` name, so pytest won't collect
it) -- it lives alongside `physical_plausibility_cli.py` and shares its
underlying engine (`planetgen/generation/plausibility.py`), but exists specifically
to let a human iterate on one class's `albedo_range`
(see `planetgen.tuning.PLANET_CLASSES`) without editing
source between runs: pass a candidate value as a CLI override, see the
resulting temperature/pressure distribution (with a delta against a
real-world reference where one exists -- Earth for Class M, Mars for K,
Venus for N), and only once satisfied, write the chosen value into
`PLANET_CLASSES` permanently. The pressure and the greenhouse of a body
are not class ranges any more (GEN.411, GEN.422): they come from its air
(`planetgen.physics.surface_air`), so only the albedo can be tuned here.

Overrides are applied as a temporary in-memory monkeypatch of
`PLANET_CLASSES[class]`, restored after the run -- nothing is written to
disk by this tool.

Usage:
    python src/tests/climate_tuning_cli.py --class M --n 300
    python src/tests/climate_tuning_cli.py --class M --n 300 --albedo 0.25 0.35
"""

import argparse
import contextlib
import logging
import sys


from planetgen.generation import plausibility
from planetgen import tuning as prog_c

logging.getLogger("transformers").setLevel(logging.ERROR)

# Real-world reference points for the classes with a direct physical analog.
# (mean_surface_temperature_k, mean_atmospheric_pressure_pa, source name).
# Classes without a direct analog (fictional/no-single-real-world-body
# classes) simply have no entry and the report omits the delta lines.
REAL_WORLD_TARGETS = {
    "M": (288.0, 101325.0, "Earth"),
    "K": (210.0, 610.0, "Mars"),
    "N": (737.0, 9.2e6, "Venus"),
}


def process_args():
    parser = argparse.ArgumentParser(description="Per-class climate tuning tool for planetGen.")
    parser.add_argument("--class", dest="planet_class", required=True, help="Planet class code to tune (e.g. M).")
    parser.add_argument("--n", type=int, default=300, help="Number of independent planets to generate (default: 300).")
    parser.add_argument("--zone", default="e", help="Zone character to generate in (default: 'e').")
    parser.add_argument("--no-moons", action="store_true", help="Skip moon generation (moons included by default).")
    parser.add_argument("--albedo", nargs=2, type=float, metavar=("LOW", "HIGH"), help="Override albedo_range for this run.")
    return parser.parse_args()


@contextlib.contextmanager
def temporary_class_overrides(cls, albedo=None):
    """
    Temporarily patches `program_constants.PLANET_CLASSES[cls]` with the
    given (low, high) albedo range override, restoring the original dict
    afterward regardless of how the block exits. Nothing is persisted.
    """
    original = prog_c.PLANET_CLASSES[cls]
    patched = dict(original)
    if albedo is not None:
        patched["albedo_range"] = tuple(albedo)
    prog_c.PLANET_CLASSES[cls] = patched
    try:
        yield patched
    finally:
        prog_c.PLANET_CLASSES[cls] = original


def _print_target_delta(cls, report):
    target = REAL_WORLD_TARGETS.get(cls)
    if not target or cls not in report:
        return
    target_temp, target_pressure, source = target
    metrics = report[cls]["metrics"]
    print(f"  --- delta from {source} reference (T={target_temp:.1f}K, P={target_pressure:.4g}Pa) ---")
    if "surface_temperature" in metrics:
        mean_t = metrics["surface_temperature"]["mean"]
        print(f"      surface_temperature: mean={mean_t:.1f}K, delta={mean_t - target_temp:+.1f}K ({(mean_t / target_temp - 1) * 100:+.1f}%)")
    if "atmospheric_pressure" in metrics:
        mean_p = metrics["atmospheric_pressure"]["mean"]
        print(f"      atmospheric_pressure: mean={mean_p:.4g}Pa, delta={mean_p - target_pressure:+.4g}Pa ({(mean_p / target_pressure - 1) * 100:+.1f}%)")


def main():
    args = process_args()
    cls = args.planet_class.upper()
    if cls not in prog_c.PLANET_CLASSES:
        print(f"Unknown planet class: {cls!r}", file=sys.stderr)
        return 2
    if not prog_c.PLANET_CLASSES[cls].get(args.zone):
        print(f"Class {cls!r} has no valid zone {args.zone!r}.", file=sys.stderr)
        return 2

    with temporary_class_overrides(cls, albedo=args.albedo) as patched:
        print(f"Generating {args.n} Class {cls} bodies in zone {args.zone!r} ...", file=sys.stderr)
        print(f"  albedo_range={patched.get('albedo_range', '(default)')}", file=sys.stderr)

        records = plausibility.generate_sample(cls, args.zone, args.n, include_moons=not args.no_moons)
        # Only this class's own records (moons may be a different class).
        records = [r for r in records if r["planet_class"] == cls]
        report = plausibility.analyze(records)

    print(plausibility.format_report(report, show_outlier_examples=False))
    _print_target_delta(cls, report)

    total_hard_violations = sum(len(v["hard_violations"]) for v in report.values())
    return 1 if total_hard_violations else 0


if __name__ == "__main__":
    sys.exit(main())
