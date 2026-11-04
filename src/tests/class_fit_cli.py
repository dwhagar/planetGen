#!/usr/bin/env python3
"""
src/tests/class_fit_cli.py -- the read-only shadow scorer of the best-fit classes (GEN.1206). Not a pytest module (no `test_*` name).

Generates systems in memory (the sampling of `class_sweep_cli`), scores every planet and moon both ways, with the first-match tree
(`planet_classes.classify`) and with the percent-match table (`physics/class_fit.py`, GEN.1205), and prints a Markdown report:

* the confusion table: the stored letter (rows) against the best-fit letter (columns);
* the share of bodies that score 100 percent for their best class, and the spread of the best match;
* the low-match list: bodies whose best match is under 50 percent (each is a hole in the class set);
* the close pairs: bodies whose runner-up is within 5 percentage points of the best;
* the temperate (250 to 330 K) surface-ocean count per class by both methods (M against O, V, F, G).

Nothing stored changes.

Usage:
    python src/tests/class_fit_cli.py --systems 300 --seed 1000
    python src/tests/class_fit_cli.py --systems 300 --seed 5 --habitable-share 0.25 --report /mnt/project-files/reports/class-fit-shadow.md
"""

import argparse
import collections
import sys
from collections.abc import Mapping

from planetgen.physics import class_fit, planet_classes
from tests import class_sweep_cli

LOW_MATCH = 0.5
CLOSE_MARGIN = 0.05
TEMPERATE_K = (250.0, 330.0)
SURFACE_WATER = ("surface ocean", "hycean", "brine seep")


def score_bodies(bodies):
    """One record per `(body, is_moon)`: the tree's letter, the fit and the quantities worth printing."""
    records = []
    for body, is_moon in bodies:
        values = body if isinstance(body, Mapping) else planet_classes.values_from_body(body)
        found = planet_classes.quantities(values)
        letter, fit, runner, runner_fit = class_fit.best_class(found)
        records.append({
            "tree": planet_classes.classify(values), "fit": letter or "-", "match": fit,
            "runner": runner or "-", "runner_match": runner_fit, "moon": is_moon, "found": found,
        })
    return records


def confusion(records):
    """`Counter` of `(stored letter, best-fit letter)`."""
    return collections.Counter((r["tree"], r["fit"]) for r in records)


def temperate_surface_ocean(records):
    """`(by_tree, by_fit)`: `Counter`s of the letters of the bodies at 250 to 330 K with a surface ocean."""
    by_tree, by_fit = collections.Counter(), collections.Counter()
    for r in records:
        t = r["found"]["temperature_k"]
        if t is not None and TEMPERATE_K[0] <= t <= TEMPERATE_K[1] and r["found"]["hydrosphere"] in SURFACE_WATER:
            by_tree[r["tree"]] += 1
            by_fit[r["fit"]] += 1
    return by_tree, by_fit


def describe_body(r):
    f = r["found"]
    parts = [f"{f['mass_earth']:.3g} Me", f"{f['gravity_g'] or 0:.2f} g", f"{f['temperature_k'] or 0:.0f} K", f"{f['pressure_pa']:.3g} Pa"]
    if f["pressure_pa"] > 0:
        parts.append(f"O2 {f['share_o2'] * 100:.1f}%, CO2 {f['share_co2'] * 100:.1f}%, reducing {f['reducing_share'] * 100:.1f}%")
    parts.append(f"water {f['ocean_fraction'] * 100:.0f}% ({f['hydrosphere']})")
    return ", ".join(parts)


def format_report(records, systems, seed, habitable_share):
    letters = sorted({r["tree"] for r in records} | {r["fit"] for r in records})
    table = confusion(records)
    total = len(records)
    full = sum(1 for r in records if r["match"] >= 1.0 - 1e-9)
    low = [r for r in records if r["match"] < LOW_MATCH]
    close = [r for r in records if r["runner"] != "-" and r["match"] - r["runner_match"] < CLOSE_MARGIN and r["match"] > 0.0]
    same = sum(n for (a, b), n in table.items() if a == b)
    lines = [f"# Best-fit class shadow, {systems} systems from seed {seed} (habitable share {habitable_share})", "",
             f"{total} planets and moons. First-match tree against the best-fit table (GEN.1205): the same letter for {same} ({same / total:.1%}); "
             f"{full} ({full / total:.1%}) score 100 percent for their best class; {len(low)} score under {LOW_MATCH:.0%}; "
             f"{len(close)} have a runner-up within {CLOSE_MARGIN * 100:.0f} points.", "",
             "## Confusion table (rows: tree letter, columns: best-fit letter)", "",
             "| tree \\ fit | " + " | ".join(letters) + " | total |", "|---|" + "---|" * (len(letters) + 1)]
    for a in letters:
        row = [table.get((a, b), 0) for b in letters]
        if any(row):
            lines.append(f"| {a} | " + " | ".join(str(n) if n else "" for n in row) + f" | {sum(row)} |")
    lines += ["", "## Where the two disagree (tree letter, best-fit letter, count, a sample body)", ""]
    for (a, b), n in sorted(table.items(), key=lambda kv: -kv[1]):
        if a != b:
            sample = next(r for r in records if r["tree"] == a and r["fit"] == b)
            lines.append(f"- {a} to {b}: {n}; e.g. {describe_body(sample)}; match {sample['match']:.2f}, runner-up {sample['runner']} {sample['runner_match']:.2f}")
    lines += ["", f"## Low-match list (best match under {LOW_MATCH:.0%}): holes in the class set", ""]
    holes = collections.Counter((r["tree"], r["fit"]) for r in low)
    if not low:
        lines.append("None.")
    for (a, b), n in holes.most_common():
        sample = next(r for r in low if r["tree"] == a and r["fit"] == b)
        lines.append(f"- tree {a}, best fit {b} ({sample['match']:.2f}): {n}; e.g. {describe_body(sample)}")
    lines += ["", f"## Close pairs (runner-up within {CLOSE_MARGIN * 100:.0f} points)", ""]
    pairs = collections.Counter(tuple(sorted((r["fit"], r["runner"]))) for r in close)
    for (a, b), n in pairs.most_common(12):
        lines.append(f"- {a} and {b}: {n}")
    if not pairs:
        lines.append("None.")
    by_tree, by_fit = temperate_surface_ocean(records)
    lines += ["", f"## Temperate ({TEMPERATE_K[0]:.0f} to {TEMPERATE_K[1]:.0f} K) bodies with a surface ocean, by class", "",
              "| class | tree | best fit |", "|---|---|---|"]
    for letter in sorted(set(by_tree) | set(by_fit)):
        lines.append(f"| {letter} | {by_tree[letter]} | {by_fit[letter]} |")
    lines.append(f"| all | {sum(by_tree.values())} | {sum(by_fit.values())} |")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--systems", type=int, default=300)
    parser.add_argument("--seed", type=int, default=1000)
    parser.add_argument("--habitable-share", type=float, default=0.0)
    parser.add_argument("--report", help="also write the report to this file")
    args = parser.parse_args(argv)
    records = score_bodies(class_sweep_cli.sample_bodies(args.systems, args.seed, args.habitable_share))
    text = format_report(records, args.systems, args.seed, args.habitable_share)
    print(text)
    if args.report:
        with open(args.report, "w", encoding="utf-8") as handle:
            handle.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
