#!/usr/bin/env python3
"""
GEN.1746 slice 4 census: how many worlds the mass extinctions change (measure only).

`collect OUT SET` draws the systems of a `class_census` set and, for each body, records what `life._reset_by_extinctions` did to the schedules of its life classes: the ages before and after, in the
record `extinction_census`. `report PATH` counts, for the dominant class of each body with life: the worlds that had animals by the star's age before the extinctions, the worlds the extinctions
changed, the worlds whose animals are gone now (still returning at the star's age), and the civilizations lost or delayed. Nothing in `src` is changed by it.
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "src"))

import class_census  # noqa: E402


def install_hook():
    """Wraps `life._reset_by_extinctions` so each body records the schedules before and after, and lets `class_census` record them."""
    from planetgen.generation import life
    original = life._reset_by_extinctions

    def wrapped(planet, cands):
        before = [(c.row.id, list(c.ages) if c.ages is not None else None) for c in cands]
        original(planet, cands)
        planet.extinction_census = {"star_age": planet.star.age, "classes": [
            {"class": class_id, "before": ages, "after": list(c.ages) if c.ages is not None else None} for (class_id, ages), c in zip(before, cands)]}

    life._reset_by_extinctions = wrapped
    class_census.BODY_ATTRS.append("extinction_census")


def collect(args):
    install_hook()
    class_census.CLASS_M = None
    class_census.draw_sample(args.out, args.set, args.workers, args.count)


def report(args):
    import json
    rows = []
    with open(args.path) as handle:
        for line in handle:
            try:
                record = json.loads(line)
            except ValueError:
                continue
            rows.extend(body["extinction_census"] for body in record.get("bodies", []) if isinstance(body.get("extinction_census"), dict))
    animals = changed = gone = civ_before = civ_lost = civ_later = 0
    for census in rows:
        star = census["star_age"]
        for one in census["classes"][:1]:       # the dominant class
            before, after = one["before"], one["after"]
            if not before or before[3] is None or before[3] > star:
                continue
            animals += 1
            if after != before:
                changed += 1
            gone += after[3] > star
            if before[4] is not None and before[4] <= star:
                civ_before += 1
                civ_lost += after[4] is None or after[4] > star
            elif before[4] is not None and after[4] != before[4]:
                civ_later += 1
    print(f"{len(rows)} bodies with a life schedule worked on; of them the dominant class had animals by the star's age on {animals}\n")
    print(f"- schedule changed by an extinction: {changed}\n- animals gone at the star's age (still returning): {gone}\n- civilizations that existed: {civ_before}, lost: {civ_lost}\n- civilizations that had not come, delayed: {civ_later}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    drawer = commands.add_parser("collect")
    drawer.add_argument("out")
    drawer.add_argument("set", choices=class_census.SETS)
    drawer.add_argument("--workers", type=int, default=3)
    drawer.add_argument("--count", type=int)
    reporter = commands.add_parser("report")
    reporter.add_argument("path")
    args = parser.parse_args()
    collect(args) if args.command == "collect" else report(args)


if __name__ == "__main__":
    main()
