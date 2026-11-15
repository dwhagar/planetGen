#!/usr/bin/env python3
"""
The Class M meter: how many Class M and Blue worlds a Class M abundance tier gives per 1,000 Sun-like systems, and where the candidates are lost.

    python scripts/class_m_meter.py run DIR [--tiers science,sparse,...] [--systems 400] [--workers 4]
    python scripts/class_m_meter.py report DIR

`run` draws the `sunlike` sample of `scripts/class_census.py` (F7V to K2V in rotation, planets on, seeds from 30000) once per tier into `DIR/<tier>.jsonl`, with the
census's worker pool, memory and load gate and resume (a system already in a file is skipped), and prints the running counts (systems done of total, candidates,
Class M, Blue) as it goes. `report` reads every `<tier>.jsonl` of DIR and prints, for each tier, the funnel from the habitable-zone rocky candidates
(`abundance_tier.candidate`: a rocky planet, not a moon, in the ecosphere zone) to Class M, then the Class M and Blue worlds per 1,000 systems with a 95 percent
interval beside the tier's target (`abundance_tier.TARGET_SYSTEMS_PER_CLASS_M`). Funnel after `reports/class-m-blue-levers-2026-11-14/hzfunnel.py`.
After the funnel each tier has the table of the habitable family (GEN.1721; Boss, 8:24 PM EST 2026-11-14: classes H, O, P and R present beside M, M in front): for
each of M, H, O, P and R the worlds of the class per 1,000 systems, and of the rocky planets with the class's shape (its own limits of temperature, water cover
and gravity) how many pass its limits of pressure, oxygen and carbon dioxide and every limit it has (`family`; after the census's
`reports/train28-census-2026-11-14/five_classes.py`). The limits are read from the class rules (`class_fit.CLASSES`), never restated here.
The same files give the same tables.
"""

import argparse
import glob
import json
import math
import os
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import class_census  # noqa: E402 -- the census draw and its pool, reused

EARTH_MASS_KG = class_census.EARTH_MASS_KG
Z95 = 1.959964

FUNNEL = (
    ("candidates (rocky, not a moon, in the ecosphere)", lambda b: True),
    ("0.4 to 6 Earth masses", lambda b: 0.4 <= (b.get("mass") or 0.0) / EARTH_MASS_KG < 6),
    ("with ground (no hydrogen envelope hides it)", lambda b: b.get("hydrosphere") is not None),
    ("gravity 0.75 to 2 g", lambda b: 0.75 <= (b.get("gravity") or 0.0) < 2),
    ("water covers 20 to 95 percent", lambda b: 0.2 <= (b.get("ocean_fraction") or 0.0) < 0.95),
    ("liquid water (surface ocean)", lambda b: b.get("hydrosphere") == "surface ocean"),
    ("Class M", lambda b: b.get("planet_class") == "M"),
    ("habitability Blue", lambda b: b.get("overall_tier") == "Blue"),
    ("habitability Green", lambda b: b.get("overall_tier") == "Green"),
)
"""tuple: The rows of the funnel, (text, test of a stored body record). Each row's test is applied on top of all the rows above it down to Class M; the Blue and Green
rows are the candidates that pass the first six rows and carry that rating."""

LIQUID_ROW = 5
CLASS_M_ROW = 6


FAMILY = ("M", "H", "O", "P", "R")
"""tuple: The classes of the habitable family, the temperate one first (they share one air: oxygen, carbon dioxide, pressure; they differ by water and temperature)."""

SHAPE_QUANTITIES = ("temperature_k", "ocean_fraction", "gravity_g")
"""tuple: The quantities whose limits are the shape of a family class in the family table (mean temperature, water cover, gravity)."""

AIR_COLUMNS = (("pressure", "pressure_pa"), ("oxygen", "o2_partial_pa"), ("carbon dioxide", "co2_partial_pa"))
"""tuple: (column text, quantity) of the air limits the family table counts one by one."""


def is_candidate(body):
    """`abundance_tier.candidate` on a stored record: a rocky planet, not a moon, in the ecosphere."""
    return body.get("body_type") == "t" and not body.get("is_moon") and body.get("zone") == "e"


def wilson(hits, total):
    """The 95 percent Wilson interval of `hits` of `total` as a pair of shares (0, 0 for no systems)."""
    if not total:
        return 0.0, 0.0
    p = hits / total
    denominator = 1 + Z95 ** 2 / total
    center = (p + Z95 ** 2 / (2 * total)) / denominator
    half = Z95 * math.sqrt(p * (1 - p) / total + Z95 ** 2 / (4 * total ** 2)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def load(path):
    """The error-free system records of one `.jsonl` file."""
    records = []
    with open(path) as handle:
        for line in handle:
            try:
                record = json.loads(line)
            except ValueError:
                continue
            if "error" not in record:
                records.append(record)
    return records


def measure(records):
    """
    The numbers of one tier's sample: `systems`, the funnel counts (`funnel`, in `FUNNEL` order), the Class M and Blue worlds anywhere in the systems (every planet
    and moon: `class_m_worlds`, `blue_worlds`) and the systems holding at least one of each (`class_m_systems`, `blue_systems`).
    """
    funnel = [0] * len(FUNNEL)
    class_m_worlds = blue_worlds = class_m_systems = blue_systems = 0
    for record in records:
        has_m = has_blue = False
        for body in record.get("bodies", []):
            is_m, is_blue = body.get("planet_class") == "M", body.get("overall_tier") == "Blue"
            class_m_worlds += is_m
            blue_worlds += is_blue
            has_m, has_blue = has_m or is_m, has_blue or is_blue
            if not is_candidate(body):
                continue
            passed = 0
            for number, (_, test) in enumerate(FUNNEL[:CLASS_M_ROW + 1]):
                try:
                    ok = bool(test(body))
                except TypeError:
                    ok = False
                if not ok:
                    break
                funnel[number] += 1
                passed = number + 1
            if passed > LIQUID_ROW:
                funnel[CLASS_M_ROW + 1] += body.get("overall_tier") == "Blue"
                funnel[CLASS_M_ROW + 2] += body.get("overall_tier") == "Green"
        class_m_systems += has_m
        blue_systems += has_blue
    return {"systems": len(records), "funnel": funnel, "class_m_worlds": class_m_worlds, "blue_worlds": blue_worlds,
            "class_m_systems": class_m_systems, "blue_systems": blue_systems, "family": family(records)}


def family(records):
    """
    The habitable family of one tier's sample, `{letter: numbers}` in `FAMILY` order. `worlds` and `systems` count the worlds of the class (planets and moons) and
    the systems holding one. The rest count rocky planets (never moons) by the class's own rules (`class_fit.CLASSES`, each limit tested with `class_fit.overshoot`
    on the record's stored `found` quantities): `shape` passes every limit the class puts on `SHAPE_QUANTITIES`; of those, each column of `AIR_COLUMNS` passes the
    class's limits on that quantity, `all` passes every limit of the class, and `got` is the classes they were given (letter -> count, most common first). The
    shapes of two classes can overlap (a hot desert has the shape of H and of R), so the shape counts do not add up to a total.
    """
    from planetgen.physics import class_fit
    rules = {row.letter: row.when for row in class_fit.CLASSES if row.letter in FAMILY}
    out = {letter: {"worlds": 0, "systems": 0, "shape": 0, "all": 0, "got": {}, **{text: 0 for text, _ in AIR_COLUMNS}} for letter in FAMILY}

    def passes(conditions, found):
        return all(class_fit.overshoot(condition, found) == 0.0 for condition in conditions)

    for record in records:
        held = set()
        for body in record.get("bodies", []):
            letter = body.get("planet_class")
            if letter in out:
                out[letter]["worlds"] += 1
                held.add(letter)
            found = body.get("found")
            if body.get("body_type") != "t" or body.get("is_moon") or not isinstance(found, dict):
                continue
            for name, conditions in rules.items():
                if not passes([c for c in conditions if c["quantity"] in SHAPE_QUANTITIES], found):
                    continue
                numbers = out[name]
                numbers["shape"] += 1
                for text, quantity in AIR_COLUMNS:
                    numbers[text] += passes([c for c in conditions if c["quantity"] == quantity], found)
                numbers["all"] += passes(conditions, found)
                numbers["got"][letter or "none"] = numbers["got"].get(letter or "none", 0) + 1
        for letter in held:
            out[letter]["systems"] += 1
    for numbers in out.values():
        numbers["got"] = dict(sorted(numbers["got"].items(), key=lambda item: (-item[1], item[0])))
    return out


def family_table(numbers, systems):
    lines = ["| class | worlds | systems holding one | per 1,000 systems (95%) | planets with its shape | pass its pressure | pass its oxygen | pass its carbon dioxide "
             "| pass every limit | classes the shaped planets got |", "|---|---|---|---|---|---|---|---|---|---|"]
    for letter in FAMILY:
        row = numbers[letter]
        got = ", ".join(f"{name} {count}" for name, count in row["got"].items()) or "-"
        lines.append(f"| {letter} | {row['worlds']} | {row['systems']} | {per_thousand(row['systems'], systems)} | {row['shape']} | "
                     + " | ".join(str(row[text]) for text, _ in AIR_COLUMNS) + f" | {row['all']} | {got} |")
    return "\n".join(lines)


def per_thousand(hits, total):
    """Text `rate (low to high)` of `hits` systems of `total`, per 1,000 systems."""
    low, high = wilson(hits, total)
    return f"{1000.0 * hits / total if total else 0.0:.1f} ({1000.0 * low:.1f} to {1000.0 * high:.1f})"


def target_text(tier):
    from planetgen.generation import abundance_tier
    per = abundance_tier.TARGET_SYSTEMS_PER_CLASS_M.get(tier)
    return "none (what the laws give)" if per is None else f"{1000.0 / per:.3g} per 1,000 (1 in {per:,})"


def tier_table(tier, numbers):
    lines = [f"### {tier}: {numbers['systems']} Sun-like systems", "", "| step | planets | share of candidates |", "|---|---|---|"]
    top = numbers["funnel"][0] or 1
    for (text, _), count in zip(FUNNEL, numbers["funnel"]):
        lines.append(f"| {text} | {count} | {100.0 * count / top:.1f}% |")
    systems = numbers["systems"]
    lines += ["", "| per 1,000 systems (95% interval) | systems holding one | worlds | target |", "|---|---|---|---|",
              f"| Class M {per_thousand(numbers['class_m_systems'], systems)} | {numbers['class_m_systems']} | {numbers['class_m_worlds']} | {target_text(tier)} |",
              f"| Blue {per_thousand(numbers['blue_systems'], systems)} | {numbers['blue_systems']} | {numbers['blue_worlds']} | |"]
    lines += ["", "The habitable family (a planet has a class's shape when it passes the class's own limits of temperature, water cover and gravity; a class is also "
              "given to a planet a little past a limit, inside the caps of the class rules, so a class can be got without passing every limit):", "",
              family_table(numbers["family"], systems)]
    return "\n".join(lines)


def summary_table(rows):
    lines = ["| tier | systems | candidates | Class M per 1,000 (95%) | Blue per 1,000 (95%) | target |", "|---|---|---|---|---|---|"]
    for tier, numbers in rows:
        lines.append(f"| {tier} | {numbers['systems']} | {numbers['funnel'][0]} | {per_thousand(numbers['class_m_systems'], numbers['systems'])} "
                     f"| {per_thousand(numbers['blue_systems'], numbers['systems'])} | {target_text(tier)} |")
    return "\n".join(lines)


def tier_files(directory):
    from planetgen.generation import abundance_tier
    found = {pathlib.Path(p).stem: p for p in glob.glob(os.path.join(directory, "*.jsonl"))}
    return [(tier, found[tier]) for tier in abundance_tier.TIERS if tier in found]


def report(args):
    rows = [(tier, measure(load(path))) for tier, path in tier_files(args.directory)]
    if not rows:
        raise SystemExit(f"class_m_meter: no <tier>.jsonl in {args.directory}")
    print("## Class M meter\n\n" + summary_table(rows))
    for tier, numbers in rows:
        print("\n" + tier_table(tier, numbers))


def run(args):
    from planetgen.generation import abundance_tier
    tiers = [t for t in args.tiers.split(",") if t]
    for tier in tiers:
        if tier not in abundance_tier.TIERS:
            raise SystemExit(f"class_m_meter: {tier!r} is not one of {', '.join(abundance_tier.TIERS)}")
    os.makedirs(args.directory, exist_ok=True)
    for number, tier in enumerate(tiers, 1):
        print(f"tier {number} of {len(tiers)}: {tier}", flush=True)
        class_census.CLASS_M = None if tier == abundance_tier.DEFAULT else tier
        class_census.draw_sample(os.path.join(args.directory, tier + ".jsonl"), "sunlike", args.workers, args.systems, every=args.every)
    report(args)


def main(argv=None):
    from planetgen.generation import abundance_tier
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    runner = sub.add_parser("run", help="draw the Sun-like sample under each tier, then print the report")
    runner.add_argument("directory")
    runner.add_argument("--tiers", default=",".join(abundance_tier.TIERS))
    runner.add_argument("--systems", type=int, default=class_census.SUNLIKE_COUNT)
    runner.add_argument("--workers", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    runner.add_argument("--every", type=int, default=10, help="print the running counts every this many systems")
    runner.set_defaults(run=run)
    reporter = sub.add_parser("report", help="print the tables of a directory of <tier>.jsonl files")
    reporter.add_argument("directory")
    reporter.set_defaults(run=report)
    args = parser.parse_args(argv)
    args.run(args)


if __name__ == "__main__":
    main()
