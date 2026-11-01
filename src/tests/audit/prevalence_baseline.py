# tests/audit/prevalence_baseline.py

"""
Builds `tests/data/prevalence_expectations.json` from a measured census (TEST.430): the one command that re-baselines the
class census, so the commit, the draws per cell and the seeds are always recorded in the file.

    python -m tests.audit.prevalence_baseline measure --count 300 --out records.json   (run at the reference commit)
    python -m tests.audit.prevalence_baseline build records.json --commit <sha> --label "..." [--tag-current records2.json
        --anomaly GEN.855 --anomaly-note "..."]

`measure` draws every cell and stores the body records. `build` turns records into ranges: a share `p` measured over `n`
bodies is allowed `p +- (4 * sqrt(p (1 - p) (1 / n + 1 / n_run)) + PAD)`, where `n_run` is the body count of a later run, so a
repeat of the reference draw passes by a wide margin and a real shift trips a flag. `--tag-current` compares a second set of
records (the code as it is now) with the new ranges and marks every range it breaks with the anomaly item, so a known fault is
named and not hidden; the test fails on a flag with no such mark.
"""

import argparse
import json
import math
import sys

from tests.audit import prevalence as census

PAD = 0.04
"""float: The width added to every range beyond the sampling error: the classes are about to be reassigned from conditions
(GEN.689), so the starting ranges are wide and tighten as the classes settle (TEST.430)."""

MIN_EXPECTED_SHARE = 0.06
"""float: A class at least this common in the reference must appear in a run (its range has a `min` above zero)."""

SIGMAS = 4.0

MUST_APPEAR = 5
"""int: A colour or life share expected this many times or more in a run must appear at least once, whatever its band. A class
share does not get the rule: the classes are about to be reassigned from conditions (GEN.689), and a class's own band is wide."""


def _range(share, n_ref, n_run, must_appear=True):
    sigma = math.sqrt(max(share * (1.0 - share), 0.0) * (1.0 / max(n_ref, 1) + 1.0 / max(n_run, 1)))
    half = SIGMAS * sigma + PAD
    low = max(0.0, share - half) if share >= MIN_EXPECTED_SHARE else 0.0
    if must_appear and share * n_run >= MUST_APPEAR:
        low = max(low, 1.0 / max(n_run, 1))  # expected five or more times: none at all is a collapse, whatever the band says
    return round(low, 4), round(min(1.0, share + half), 4)


def zone_barred(zone, kind):
    """The classes barred from a zone: none since the class cut-over (GEN.646), which derives a body's class from its conditions, so
    any class may come out anywhere the conditions give it."""
    return []


GAS_DWARF_RANGE = {
    "min": 0.02, "max": 0.45, "kind": "stand-in",
    "source": "GEN.350 keeps class T (the sub-Neptune row); sub-Neptunes are the most common planets observed "
              "([R1459, R3499], recalled, unverified, tier T2/T3), so class T must appear among cold planets; the "
              "pre-cut-over spine drew it almost never"}
"""dict: The one range not taken from the reference: class T among cold planets (Boss: 'a class missing where expected, class T
today'). A stand-in at design tier, tweak level L3 (hand-set against the observed abundance of sub-Neptunes)."""


GAS_DWARF_NOTE = ("class T (gas dwarf) is never drawn: no 2 to 20 Earth-mass sub-Neptune population exists (GEN.829); the "
                  "pre-cut-over spine drew it almost never either")


def build(records, commit, label, count, origin="pre-cut-over"):
    cells = census.tally(records)
    out = {"reference": {"commit": commit, "label": label, "draws_per_cell": count,
                         "seeds": "crc32 of 'cell:<star>:<zone>' and 'orbit:<star>:<zone>' (tests.audit.prevalence.seed_for)",
                         "kind": f"stand-in: {origin} prevalence, no observed class rates exist",
                         "pad": PAD, "sigmas": SIGMAS},
           "anomalies": {}, "hangs": {"anomaly": ""}, "errors": {"anomaly": ""}, "cells": {}}
    for (star, zone, kind), cell in sorted(cells.items()):
        if not cell["n"] or kind not in ("planet", "moon"):
            continue
        spec = {"n_ref": cell["n"], "classes": {}, "tiers": {}, "barred": zone_barred(zone, kind),
                "source": f"{origin} measured at {commit[:9]}"}
        for letter, count_ in sorted(cell["cls"].items()):
            if letter in spec["barred"]:
                continue
            low, high = _range(count_ / cell["n"], cell["n"], cell["n"], must_appear=False)
            spec["classes"][letter] = {"min": low, "max": high, "kind": "stand-in",
                                       "source": f"{origin} share {count_ / cell['n']:.3f} at {commit[:9]}"}
        scored = sum(cell["tier"].values())
        for tier in census.TIERS:
            share = cell["tier"][tier] / scored if scored else 0.0
            low, high = _range(share, scored, scored)
            spec["tiers"][tier] = {"min": low, "max": high, "kind": "stand-in",
                                   "source": f"{origin} share {share:.3f} of {scored} scored bodies at {commit[:9]}"}
        spec["ceilings"] = {}
        for ceiling in census.CEILINGS:
            share = cell["ceiling"][ceiling] / scored if scored else 0.0
            low, high = _range(share, scored, scored)
            spec["ceilings"][ceiling] = {"min": low, "max": high, "kind": "stand-in",
                                         "source": f"{origin} share {share:.3f} of {scored} scored bodies at {commit[:9]}"}
        out["cells"][census.cell_key(star, zone, kind)] = spec
    for key, spec in out["cells"].items():
        if key.endswith("/cold/planet"):
            spec["classes"]["T"] = dict(GAS_DWARF_RANGE)
    return out


def tag_current(expectations, records, item, note, only=None):
    """Marks every untagged range the current code breaks with `item` (those `only` accepts, when given); the file keeps the
    note."""
    findings = census.compare(census.tally(records), expectations)
    tagged = 0
    for finding in findings:
        if finding.kind in ("hang", "error"):
            continue
        spec = expectations["cells"][finding.cell]
        prefix, name = finding.measure.split(" ", 1)
        group = next(g for g, p in census.GROUPS if p == prefix)
        if finding.kind == "barred" or name not in spec[group] or spec[group][name].get("anomaly"):
            continue
        if only and not only(finding):
            continue
        spec[group][name]["anomaly"] = item
        tagged += 1
    expectations["anomalies"][item] = note
    return tagged


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    measure = sub.add_parser("measure")
    measure.add_argument("--count", type=int, default=300)
    measure.add_argument("--out", required=True)
    build_p = sub.add_parser("build")
    build_p.add_argument("records")
    build_p.add_argument("--commit", required=True)
    build_p.add_argument("--label", default="")
    build_p.add_argument("--count", type=int, required=True)
    build_p.add_argument("--out", default=str(census.EXPECTATIONS_PATH))
    build_p.add_argument("--tag-current")
    build_p.add_argument("--anomaly")
    build_p.add_argument("--anomaly-note", default="")
    build_p.add_argument("--origin", default="pre-cut-over", help="what the records are, in the source text of every range")
    build_p.add_argument("--hang-anomaly", default="", help="the item that explains a draw over the time limit")
    build_p.add_argument("--error-anomaly", default="", help="the item that explains a draw that raised")
    args = parser.parse_args(argv)
    if args.command == "measure":
        records = census.draw_cells(args.count)
        with open(args.out, "w", encoding="utf-8") as handle:
            json.dump(records, handle)
        print(f"{len(records)} bodies")
        return 0
    with open(args.records, encoding="utf-8") as handle:
        records = json.load(handle)
    expectations = build(records, args.commit, args.label, args.count, args.origin)
    expectations["hangs"]["anomaly"] = args.hang_anomaly
    expectations["errors"]["anomaly"] = args.error_anomaly
    if args.tag_current:
        with open(args.tag_current, encoding="utf-8") as handle:
            current = json.load(handle)
        print("tagged T", tag_current(expectations, current, "GEN.829", GAS_DWARF_NOTE,
                                      only=lambda finding: finding.measure == "class T"))
        print("tagged", tag_current(expectations, current, args.anomaly, args.anomaly_note))
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(expectations, handle, indent=1, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
