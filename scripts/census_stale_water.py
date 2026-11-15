#!/usr/bin/env python3
"""
GEN.1590 check: how many terrestrial bodies of a class census file (`python scripts/class_census.py collect OUT.jsonl SET`) store a liquid ocean at a temperature
and pressure where water boils (a hydrosphere left as the life pass found it). Reports only; it changes nothing.

    python scripts/census_stale_water.py OUT.jsonl [OUT2.jsonl ...]
"""

import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

from planetgen.physics import hydrosphere  # noqa: E402

LIQUID = ("surface ocean", "brine seep", "hycean")


def stale_bodies(path):
    """`(bodies read, [(hydrosphere, K, kPa, class)])` for the terrestrial bodies of one census file whose stored liquid is above its boiling point."""
    seen, stale = 0, []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            for body in json.loads(line)["bodies"]:
                kelvin, pascal, state = body.get("surface_temperature"), body.get("atmospheric_pressure"), body.get("hydrosphere")
                if kelvin is None or pascal is None or body.get("body_type") != "t":
                    continue
                seen += 1
                if state in LIQUID and pascal > 0 and kelvin >= hydrosphere.water_boiling_k(pascal):
                    stale.append((state, round(kelvin), round(pascal / 1000), body.get("planet_class")))
    return seen, stale


def main(paths):
    for path in paths:
        seen, stale = stale_bodies(path)
        print(f"{path}: {seen} terrestrial bodies, {len(stale)} store a liquid above its boiling point {stale[:10]}")


if __name__ == "__main__":
    main(sys.argv[1:])
