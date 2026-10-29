#!/usr/bin/env python3
"""
src/tests/life_class_shadow_cli.py -- the shadow sweep of the life classes (GEN.886, docs/design/habitability/life-classes.md section
10.3). Not a pytest module (no `test_*` name).

Generates systems in memory, scores every rocky planet and moon for every life class (`physics.life_classes.score_classes`) and
prints (and with `--report` writes) a Markdown report: for each class how many bodies it supports and how strongly, what it waited
for (a quantity the generator does not store yet, a placeholder supply), which class is the best of each body, how the best class
and the old life chemistry agree, and where Class 1a differs from today's `PHI_bio`. Nothing is stored and no draw is made.

Usage:
    python src/tests/life_class_shadow_cli.py
    python src/tests/life_class_shadow_cli.py --systems 800 --seed 1 --habitable-share 0.2 --report /mnt/project-files/reports/life-class-shadow-sweep.md
"""

import argparse
import collections
import sys

from planetgen import tuning
from planetgen.physics import life_classes
from tests import class_shadow_cli


def rocky(bodies):
    """The `(body, is_moon)` pairs a class can be scored for: a body with a surface temperature that is not a gas giant."""
    return [(b, m) for b, m in bodies if b.body_type != "g" and getattr(b, "surface_temperature", None) is not None]


def sweep(bodies):
    """`[(body, scores)]` for every rocky body of `bodies`."""
    return [(body, life_classes.score_classes(body)) for body, _ in rocky(bodies)]


def _letter(body):
    return (getattr(getattr(body, "star", None), "type", None) or "?")[:1]


def summarise(results):
    """Counters the report is made of."""
    supports = collections.Counter()
    strong = collections.Counter()
    waited = collections.defaultdict(collections.Counter)
    placeholder = collections.defaultdict(collections.Counter)
    best_of = collections.Counter()
    best_by_letter = collections.defaultdict(collections.Counter)
    for body, scores in results:
        top = life_classes.best(scores)
        best_of[top.id if top else "none"] += 1
        best_by_letter[_letter(body)][top.id if top else "none"] += 1
        for class_id, result in scores.items():
            if result.score > 0:
                supports[class_id] += 1
                for quantity in result.placeholders:
                    placeholder[class_id][quantity] += 1
            if result.score >= 5:
                strong[class_id] += 1
            for quantity in result.missing:
                waited[class_id][quantity] += 1
    return supports, strong, waited, placeholder, best_of, best_by_letter


def compare_with_index(results):
    """`(same, different, rows)`: bodies where Class 1a should be today's `PHI_bio` (a Sun-like star's band, sunlit water, a window
    that holds) and does or does not match it; `rows` lists up to ten that differ."""
    same = different = 0
    rows = []
    for body, scores in results:
        phi_bio = getattr(body, "phi_bio", None)
        one_a = scores["1a"]
        if phi_bio is None or one_a.failed or any(name.startswith("window:") and value < 1.0 for name, value in one_a.terms.items()):
            continue
        teff = getattr(getattr(body, "star", None), "temperature", None)
        if not teff or life_classes.band_fraction(teff, (400.0, 700.0)) < life_classes.band_fraction(life_classes.REFERENCE_STAR_K, (400.0, 700.0)):
            continue
        if life_classes.score_of(phi_bio) == one_a.score:
            same += 1
        else:
            different += 1
            if len(rows) < 10:
                rows.append((body, phi_bio, one_a))
    return same, different, rows


def format_report(results, systems, seed, habitable_share):
    supports, strong, waited, placeholder, best_of, best_by_letter = summarise(results)
    total = len(results)
    lines = [
        "# Life classes: shadow sweep", "",
        f"Rule file version {life_classes.VERSION}; {systems} generated systems from run seed {seed}"
        + (f", the first {habitable_share:.0%} asked for a habitable world" if habitable_share else "")
        + f"; {total} rocky planets and moons scored. Nothing is stored; no draw is made.", "",
        "## What each class supports", "",
        f"A body supports a class when its score is at least 1 (a likelihood of {tuning.HABITABILITY_LIFE_SCORE_MIN} or more); strong is 5 or more.", "",
        "| Class | Name | Supports | Share | Strong (5+) | Waits for | Placeholder supply |", "|---|---|---|---|---|---|---|",
    ]
    for row in life_classes.CLASSES:
        wait = ", ".join(f"{q} ({n})" for q, n in waited[row.id].most_common(2)) or "-"
        hold = ", ".join(f"{q} ({n})" for q, n in placeholder[row.id].most_common(2)) or "-"
        lines.append(f"| {row.id} | {row.name} | {supports[row.id]} | {100.0 * supports[row.id] / total:.1f}% | {strong[row.id]} | {wait} | {hold} |")
    lines += ["", "## The best class of each body", "", "| Best class | Bodies | Share |", "|---|---|---|"]
    for class_id, n in sorted(best_of.items(), key=lambda item: (-item[1], item[0])):
        lines.append(f"| {class_id} | {n} | {100.0 * n / total:.1f}% |")
    letters = sorted(best_by_letter)
    lines += ["", "## The best class by host star letter", "", "| Best \\ Star | " + " | ".join(letters) + " |", "|---" * (len(letters) + 1) + "|"]
    for class_id in [*life_classes.IDS, "none"]:
        counts = [best_by_letter[letter][class_id] for letter in letters]
        if any(counts):
            lines.append(f"| {class_id} | " + " | ".join(str(n or "") for n in counts) + " |")
    same, different, rows = compare_with_index(results)
    lines += ["", "## Class 1a against today's PHI_bio", "",
              f"Of the bodies where Class 1a should be exactly today's PHI_bio (sunlit water, all windows and gates holding, a star at least as rich "
              f"in 400 to 700 nm photons as the Sun): **{same} match and {different} differ**.", ""]
    if rows:
        lines += ["| Body | PHI_bio score | Class 1a score | Terms |", "|---|---|---|---|"]
        for body, phi_bio, one_a in rows:
            terms = ", ".join(f"{k} {v:.3f}" for k, v in one_a.terms.items())
            lines.append(f"| {getattr(body, 'name', '?')} | {life_classes.score_of(phi_bio)} | {one_a.score} | {terms} |")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--systems", type=int, default=400)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--habitable-share", type=float, default=0.2, help="share of the systems generated with a forced habitable world")
    parser.add_argument("--report", help="write the Markdown report to this file")
    args = parser.parse_args(argv)
    bodies = class_shadow_cli.sample_bodies(args.systems, args.seed, args.habitable_share)
    report = format_report(sweep(bodies), args.systems, args.seed, args.habitable_share)
    print(report)
    if args.report:
        with open(args.report, "w", encoding="utf-8") as handle:
            handle.write(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
