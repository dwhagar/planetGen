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
* with `--fudge`: what the fudge rule did (GEN.1314, fudge-table-from-stand-in-confidence section 4c): how many bodies fit a class fully, how many were placed at each cap tier by moving
  how many values, which quantities were stretched, and which bodies are holes;
* with `--dump FILE`: the class code of every body, written as JSON; with `--compare-dump FILE`: how many bodies change class code against such a dump (made by this tool on another
  revision from the same systems and seed), each change listed as from to with its count.

Nothing is stored.

Usage:
    python src/tests/class_fit_cli.py --systems 300 --seed 1000
    python src/tests/class_fit_cli.py --systems 300 --seed 5 --habitable-share 0.25 --report /mnt/project-files/reports/class-fit-sweep.md
    python src/tests/class_fit_cli.py --systems 300 --seed 1000 --habitable-share 0.25 --dump /tmp/before.json       (on the older revision), then --compare-dump /tmp/before.json
    python src/tests/class_fit_cli.py --systems 150 --seed 1000 --fudge --report /mnt/project-files/reports/fudge-sweep.md
"""

import argparse
import collections
import json
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
            "fudge_values": () if detail is None else getattr(detail, "fudge_values", ()), "fudge_amount": None if detail is None else getattr(detail, "fudge_amount", None),
            "full": class_fit._found(values),
        })
    return records


def fudge_summary(records):
    """The counts of the fudge rule over the scored `records`: bodies that fit fully, bodies placed with values moved (by how many values), digit 9, and holes."""
    rows = [r for r in records if r["fit"] != "Ast"]
    full = [r for r in rows if r["match"] > 0.0 and not r["fudge_values"]]
    fudged = [r for r in rows if r["match"] > 0.0 and r["fudge_values"]]
    holes = [r for r in rows if r["match"] <= 0.0]
    return {"total": len(rows), "full": full, "fudged": fudged, "holes": holes}


def format_fudge(records):
    summary = fudge_summary(records)
    total = summary["total"]
    fudged, holes = summary["fudged"], summary["holes"]
    lines = ["", "## The fudge rule", "",
             f"{total} bodies. {len(summary['full'])} ({len(summary['full']) / max(total, 1):.1%}) fit a class fully, so nothing is fudged; {len(fudged)} ({len(fudged) / max(total, 1):.1%}) were placed by moving "
             f"values (within the tiers {', '.join(f'{cap:g}' for cap in class_fit.CAP_TIERS)} widths); {len(holes)} ({len(holes) / max(total, 1):.1%}) are holes (no class within the largest cap, the nearest is shown at 0 percent).", "",
             "Values moved by the bodies that were fudged:", ""]
    lines += [f"- {n} value(s): {count}" for n, count in sorted(collections.Counter(len(r["fudge_values"]) for r in fudged).items())]
    lines += ["", "Total move of the fudged bodies, in widths:", ""]
    buckets = collections.Counter(min(int(r["fudge_amount"] / 0.25), 3) for r in fudged if r["fudge_amount"] is not None)
    lines += [f"- {('under ' + format((bucket + 1) * 0.25, 'g')) if bucket < 3 else '0.75 or more'}: {count}" for bucket, count in sorted(buckets.items())]
    stretched = collections.Counter(value for r in fudged for value in r["fudge_values"])
    lines += ["", "Quantities the rule stretches:", ""] + [f"- {key}: {count}" for key, count in stretched.most_common()]
    letters = collections.Counter(class_code.format_code(r["fit"], r["digit"]) for r in fudged)
    lines += ["", "Classes the fudged bodies take:", ""] + [f"- {code}: {count}" for code, count in letters.most_common(12)]
    if holes:
        near = collections.Counter(r["fit"] for r in holes)
        lines += ["", "Nearest class of the holes (each group of them is a missing class, GEN.1306):", ""] + [f"- {letter}: {count}" for letter, count in near.most_common(12)]
    return "\n".join(lines) + "\n"


def code_of_record(r):
    return "Ast" if r["fit"] == "Ast" else class_code.format_code(r["fit"], r["digit"])


def dump(records, path):
    """Writes the class code of every body, in order, with the mass and mean temperature that identify it, as JSON."""
    rows = [{"code": code_of_record(r), "mass": r["found"]["mass_earth"], "temperature": r["found"]["temperature_k"]} for r in records]
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(rows, handle)


def changes_against(records, before):
    """`(changed, total, Counter, details)`: how many bodies take another class code than the dump `before` (a list as `dump` writes it) records, each `(before, after)` pair with its count, and
    one line for each changed body (its quantities, the values the fudge rule moved, and the runner-up). The bodies must be the same ones in the same order; a mismatch of a mass or a
    temperature is an error."""
    if len(before) != len(records):
        raise SystemExit(f"the dump holds {len(before)} bodies and this run {len(records)}: use the same --systems, --seed and --habitable-share")
    pairs, total, details = collections.Counter(), 0, []
    for r, old in zip(records, before):
        if r["fit"] == "Ast":
            continue
        total += 1
        if abs((old["mass"] or 0.0) - (r["found"]["mass_earth"] or 0.0)) > 1e-9 * max(1.0, abs(old["mass"] or 0.0)):
            raise SystemExit("the bodies differ between the dump and this run (a generator change moved them): compare only a change that leaves the bodies alone")
        after = code_of_record(r)
        if old["code"] != after:
            pairs[(old["code"], after)] += 1
            details.append(f"{old['code']} to {after}: {describe_body(r)}; fudged {','.join(r['fudge_values']) or 'nothing'} by {r['fudge_amount'] or 0.0:.2f} widths, match {r['match']:.2f}, runner-up {r['runner']} {r['runner_match']:.2f}")
    return sum(pairs.values()), total, pairs, details


def format_changes(label, changed, total, pairs, details=()):
    lines = ["", f"## Class codes against {label}", "",
             f"{changed} of {total} bodies ({changed / max(total, 1):.1%}) take another class code than they did.", ""]
    if pairs:
        lines += ["| before | after | bodies |", "|---|---|---|"]
        lines += [f"| {a} | {b} | {n} |" for (a, b), n in pairs.most_common()]
        lines += ["", "Each changed body:", ""] + [f"- {line}" for line in details]
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
    parser.add_argument("--fudge", action="store_true", help="also print what the fudge rule did: the bodies fitted fully, fudged by how many values, and the holes")
    parser.add_argument("--dump", metavar="FILE", help="write the class code of every body as JSON")
    parser.add_argument("--compare-dump", metavar="FILE", help="also count the bodies whose class code differs from a dump made on another revision from the same systems")
    args = parser.parse_args(argv)
    records = score_bodies(class_sweep_cli.sample_bodies(args.systems, args.seed, args.habitable_share))
    text = format_report(records, args.systems, args.seed, args.habitable_share)
    if args.fudge:
        text += format_fudge(records)
    if args.dump:
        dump(records, args.dump)
    if args.compare_dump:
        with open(args.compare_dump, encoding="utf-8") as handle:
            text += format_changes(args.compare_dump, *changes_against(records, json.load(handle)))
    print(text)
    if args.report:
        with open(args.report, "w", encoding="utf-8") as handle:
            handle.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
