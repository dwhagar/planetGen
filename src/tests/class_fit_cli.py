#!/usr/bin/env python3
"""
src/tests/class_fit_cli.py -- the best-fit class sweep (GEN.1206, GEN.1207). Not a pytest module (no `test_*` name).

Generates systems in memory (the sampling of `class_sweep_cli`), classifies every planet and moon with the best-fit classifier
(`planet_classes.classify_detail`) and prints a Markdown report:

* the share of each letter and of each class code (letter and digit);
* the spread of the match: how many bodies score 100 percent for their class, and the list of bodies under 50 percent (each is a hole in the class set:
  digit 9 is a body limited by its gravity or air amount, a match of 0 is a body that fits no class and shows its nearest);
* the close pairs: bodies whose runner-up is within 5 percentage points of the best;
* the temperate (250 to 330 K) surface-ocean count per class;
* with `--fudge`: the no-cap fudge sweep (GEN.1312, fudge-table-from-stand-in-confidence section 4c and 4d). For every body, the class that needs the fewest values moved
  (then the smallest move, in widths of the shoulder of the quantity) with no cap, how many bodies fit fully, and the share of bodies a cap would place, so the
  tiers of the cap are read off the knees of the curve;
* with `--compare REV`: how many bodies change class code against the fit table of git revision REV (GEN.1321: the oxygen and carbon dioxide edges became pressures),
  each change listed as from to with its count.

Nothing is stored.

Usage:
    python src/tests/class_fit_cli.py --systems 300 --seed 1000
    python src/tests/class_fit_cli.py --systems 300 --seed 5 --habitable-share 0.25 --report /mnt/project-files/reports/class-fit-sweep.md
    python src/tests/class_fit_cli.py --systems 300 --seed 1000 --habitable-share 0.25 --compare origin/main
    python src/tests/class_fit_cli.py --systems 150 --seed 1000 --fudge --report /mnt/project-files/reports/fudge-cap-sweep.md
"""

import argparse
import collections
import contextlib
import json
import subprocess
import sys
from collections.abc import Mapping

from planetgen.physics import class_code, class_fit, planet_classes
from tests import class_sweep_cli

LOW_MATCH = 0.5
CLOSE_MARGIN = 0.05
TEMPERATE_K = (250.0, 330.0)
SURFACE_WATER = ("surface ocean", "hycean", "brine seep")


def score_bodies(bodies):
    """One record per `(body, is_moon)`: the best-fit class and the quantities worth printing. An Asteroid has the letter `Ast`."""
    records = []
    for body, is_moon in bodies:
        values = body if isinstance(body, Mapping) else planet_classes.values_from_body(body)
        found = planet_classes.quantities(values)
        detail = planet_classes.classify_detail(values)
        records.append({
            "fit": "Ast" if detail is None else detail.letter, "digit": 0 if detail is None else detail.digit,
            "match": 1.0 if detail is None else detail.match, "runner": "-" if detail is None or detail.runner_up is None else detail.runner_up,
            "runner_match": 0.0 if detail is None else detail.runner_up_match, "moon": is_moon, "found": found,
            "full": class_fit._found(values),
        })
    return records


THERMAL = ("temperature_k", "hottest_region_k", "coldest_region_k", "day_temperature_k", "night_temperature_k")
"""tuple: The quantities that move together when the mean temperature moves; they count as one value in the fudge rule."""
CAPS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0)
"""tuple: The caps (in widths) the fudge sweep counts the placed bodies at."""


def overshoot(condition, found):
    """How far past the edge of one class condition the quantities `found` are, in widths of the quantity's shoulder: 0 inside, infinity for a quantity that is
    missing, an `equals` that differs, or an edge with no width."""
    quantity = condition["quantity"]
    value = found.get(quantity)
    if value is None:
        return float("inf")
    if "equals" in condition:
        return 0.0 if value == condition["equals"] else float("inf")
    if "min" in condition:
        edge = condition["min"]
        width = class_fit.shoulder_width(quantity, edge)
        return 0.0 if value >= edge else ((edge - value) / width if width > 0.0 else float("inf"))
    edge = condition["lt"]
    width = class_fit.shoulder_width(quantity, edge)
    return 0.0 if value < edge else ((value - edge) / width if width > 0.0 else float("inf"))


def class_options(found):
    """`[(letter, conditions, {value: widths})]` for every class whose gates the body passes: the values it would have to move (the thermal ones as one) and how far."""
    options = []
    for fit_class in class_fit.CLASSES:
        if not class_fit.passes_gates(fit_class, found):
            continue
        moves = {}
        for condition in fit_class.when:
            z = overshoot(condition, found)
            if z > 0.0:
                key = "thermal" if condition["quantity"] in THERMAL else condition["quantity"]
                moves[key] = max(moves.get(key, 0.0), z)
        options.append((fit_class.letter, len(fit_class.when), moves))
    return options


def minimal_fudge(options, cap=float("inf")):
    """`(values moved, total widths, letter, {value: widths})` of the class the rule takes, or `None` when no class is within `cap`: the fewest values, then the
    smallest move, then more conditions, then the letter (the tie-break of a full fit)."""
    eligible = [(len(moves), sum(moves.values()), -count, letter, moves) for letter, count, moves in options if all(z <= cap for z in moves.values())]
    if not eligible:
        return None
    n, total, _count, letter, moves = min(eligible, key=lambda t: (t[0], round(t[1], 9), t[2], t[3]))
    return n, total, letter, moves


def fudge_sweep(records):
    """The no-cap fudge sweep over the scored `records` (`score_bodies`): per body the options and the best choice with no cap."""
    rows = []
    for r in records:
        if r["fit"] == "Ast":
            continue
        options = class_options(r["full"])
        rows.append({"options": options, "best": minimal_fudge(options), "found": r["found"], "today": r["fit"]})
    return rows


def format_fudge(rows):
    total = len(rows)
    needing = [r for r in rows if r["best"] and r["best"][0] > 0]
    full = sum(1 for r in rows if r["best"] and r["best"][0] == 0)
    nothing = sum(1 for r in rows if r["best"] is None)
    lines = ["", "## The fudge sweep (no cap)", "",
             f"{total} bodies. {full} ({full / max(total, 1):.1%}) fit a class fully, so nothing is fudged; {len(needing)} need a fudge; {nothing} fit no class at any cap.", "",
             "Values moved by the bodies that need a fudge (the rule takes the fewest values first):", ""]
    lines += [f"- {n} value(s): {count}" for n, count in sorted(collections.Counter(r["best"][0] for r in needing).items())]
    lines += ["", "| cap (widths) | bodies placed | share | added |", "|---|---|---|---|"]
    previous = 0
    for cap in CAPS:
        placed = sum(1 for r in rows if minimal_fudge(r["options"], cap) is not None)
        lines.append(f"| {cap:g} | {placed} | {placed / max(total, 1):.1%} | {placed - previous:+d} |")
        previous = placed
    stretched = collections.Counter(key for r in needing for key in r["best"][3])
    lines += ["", "Quantities the rule stretches:", ""] + [f"- {key}: {count}" for key, count in stretched.most_common()]
    letters = collections.Counter(f"{r['today']} to {r['best'][2]}" for r in needing if r["best"][2] != r["today"])
    lines += ["", f"{sum(letters.values())} of the {len(needing)} fudged bodies take another letter than today's nearest-class choice:", ""]
    lines += [f"- {pair}: {count}" for pair, count in letters.most_common(12)]
    return "\n".join(lines) + "\n"


def table_at(revision):
    """The text of `planet_class_fit.json` at git revision `revision`."""
    path = "src/planetgen/physics/data/planet_class_fit.json"
    return subprocess.run(["git", "show", f"{revision}:{path}"], check=True, capture_output=True, text=True).stdout


@contextlib.contextmanager
def fit_table(text):
    """Swaps the classifier's table for the one in `text` (a `planet_class_fit.json`) while the block runs."""
    document = json.loads(text)
    saved = (class_fit.CLASSES, class_fit.SHOULDERS, class_fit.GATE_PRESSURE_PA)
    class_fit.CLASSES = tuple(class_fit.FitClass(row) for row in document["classes"])
    class_fit.SHOULDERS, class_fit.GATE_PRESSURE_PA = document["shoulders"], document["gate_pressure_pa"]
    try:
        yield
    finally:
        class_fit.CLASSES, class_fit.SHOULDERS, class_fit.GATE_PRESSURE_PA = saved


def code_of(found):
    """The class code (letter and digit) the best-fit classifier gives the stored quantities `found` with the table now in force."""
    detail = class_fit.classify(found)
    return class_code.format_code(detail.letter, detail.digit)


def changes_against(records, text):
    """`(changed, total, Counter)`: how many gaseous-or-solid bodies take another class code under the table in `text` than under the table now, and each
    `(before, after)` pair with its count. A record whose body the classifier skips (an asteroid) is left out."""
    pairs, total = collections.Counter(), 0
    for r in records:
        if r["fit"] == "Ast":
            continue
        total += 1
        after = class_code.format_code(r["fit"], r["digit"])
        with fit_table(text):
            before = code_of(r["full"])
        if before != after:
            pairs[(before, after)] += 1
    return sum(pairs.values()), total, pairs


def format_changes(revision, changed, total, pairs):
    lines = ["", f"## Class codes against the fit table at {revision}", "",
             f"{changed} of {total} bodies ({changed / max(total, 1):.1%}) take another class code than they did.", ""]
    if pairs:
        lines += ["| before | after | bodies |", "|---|---|---|"]
        lines += [f"| {a} | {b} | {n} |" for (a, b), n in pairs.most_common()]
    return "\n".join(lines) + "\n"


def temperate_surface_ocean(records):
    """`Counter` of the letters of the bodies at 250 to 330 K with a surface ocean."""
    counts = collections.Counter()
    for r in records:
        t = r["found"]["temperature_k"]
        if t is not None and TEMPERATE_K[0] <= t <= TEMPERATE_K[1] and r["found"]["hydrosphere"] in SURFACE_WATER:
            counts[r["fit"]] += 1
    return counts


def describe_body(r):
    f = r["found"]
    parts = [f"{f['mass_earth']:.3g} Me", f"{f['gravity_g'] or 0:.2f} g", f"{f['temperature_k'] or 0:.0f} K", f"{f['pressure_pa']:.3g} Pa"]
    if f["pressure_pa"] > 0:
        parts.append(f"O2 {f['o2_partial_pa'] / 1000.0:.2f} kPa ({f['share_o2'] * 100:.1f}%), CO2 {f['co2_partial_pa'] / 1000.0:.2f} kPa ({f['share_co2'] * 100:.1f}%), reducing {f['reducing_share'] * 100:.1f}%")
    parts.append(f"water {f['ocean_fraction'] * 100:.0f}% ({f['hydrosphere']})")
    return ", ".join(parts)


def format_report(records, systems, seed, habitable_share):
    total = len(records)
    full = sum(1 for r in records if r["match"] >= 1.0 - 1e-9)
    low = [r for r in records if r["match"] < LOW_MATCH]
    nothing = [r for r in records if r["match"] <= 0.0]
    limited = [r for r in records if r["digit"] == class_code.DIGITS[9]]
    close = [r for r in records if r["runner"] != "-" and r["match"] - r["runner_match"] < CLOSE_MARGIN and r["match"] > 0.0]
    lines = [f"# Best-fit class sweep, {systems} systems from seed {seed} (habitable share {habitable_share})", "",
             f"{total} planets and moons. {full} ({full / max(total, 1):.1%}) score 100 percent for their class; {len(low)} score under {LOW_MATCH:.0%}, of which "
             f"{len(nothing)} fit no class (the nearest is shown, at 0 percent); {len(limited)} are limited by another domain (digit 9); "
             f"{len(close)} have a runner-up within {CLOSE_MARGIN * 100:.0f} points.", "",
             "## Share of each class code", "", "| code | bodies | share |", "|---|---|---|"]
    codes = collections.Counter(class_code.format_code(r["fit"], r["digit"]) if r["fit"] != "Ast" else "Ast" for r in records)
    for code, n in sorted(codes.items()):
        lines.append(f"| {code} | {n} | {n / max(total, 1):.1%} |")
    lines += ["", f"## Low-match list (best match under {LOW_MATCH:.0%}): holes in the class set", ""]
    holes = collections.Counter((r["fit"], r["digit"]) for r in low)
    if not low:
        lines.append("None.")
    for (letter, digit), n in holes.most_common():
        sample = next(r for r in low if r["fit"] == letter and r["digit"] == digit)
        lines.append(f"- {class_code.format_code(letter, digit)} ({sample['match']:.2f}): {n}; e.g. {describe_body(sample)}")
    lines += ["", f"## Close pairs (runner-up within {CLOSE_MARGIN * 100:.0f} points)", ""]
    pairs = collections.Counter(tuple(sorted((r["fit"], r["runner"]))) for r in close)
    for (a, b), n in pairs.most_common(12):
        lines.append(f"- {a} and {b}: {n}")
    if not pairs:
        lines.append("None.")
    counts = temperate_surface_ocean(records)
    lines += ["", f"## Temperate ({TEMPERATE_K[0]:.0f} to {TEMPERATE_K[1]:.0f} K) bodies with a surface ocean, by class", "", "| class | bodies |", "|---|---|"]
    for letter in sorted(counts):
        lines.append(f"| {letter} | {counts[letter]} |")
    lines.append(f"| all | {sum(counts.values())} |")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--systems", type=int, default=300)
    parser.add_argument("--seed", type=int, default=1000)
    parser.add_argument("--habitable-share", type=float, default=0.0)
    parser.add_argument("--report", help="also write the report to this file")
    parser.add_argument("--fudge", action="store_true", help="also run the no-cap fudge sweep and print the share of bodies placed against the cap")
    parser.add_argument("--compare", metavar="REV", help="also count the bodies whose class code differs from the one the fit table at git revision REV gives")
    args = parser.parse_args(argv)
    records = score_bodies(class_sweep_cli.sample_bodies(args.systems, args.seed, args.habitable_share))
    text = format_report(records, args.systems, args.seed, args.habitable_share)
    if args.fudge:
        text += format_fudge(fudge_sweep(records))
    if args.compare:
        text += format_changes(args.compare, *changes_against(records, table_at(args.compare)))
    print(text)
    if args.report:
        with open(args.report, "w", encoding="utf-8") as handle:
            handle.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
