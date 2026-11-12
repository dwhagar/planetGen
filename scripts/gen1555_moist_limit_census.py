#!/usr/bin/env python3
"""
GEN.1555 (prebuilt by Lane 07 for Lane 14): the score of a fix to the bistable surface air solve.

The air solve (`surface_air.solve`) has more than one fixed point near the moist limit and finds one by where it starts. Over a seeded sweep this script lists

* the bodies that sit EXACTLY on `surface_air.MOIST_LIMIT_K` (the pin of the `held_cool` rule, 375.0 K),
* the bodies near the boiling point of their own surface pressure (within `NEAR_BOILING_K`, which `test_surface_air` swings at a pinned 10 kPa), and
* for every body that kept at least 10 kPa of air or sits in one of those two groups, the start dependence: the same final `Setting` solved from `STARTS_K` with one draw
  stream, and the spread of the answers (a body is `bistable` (fresh probe stream) and `bistable_production` (the draw state of the production solve, so a pinned body is judged under its own design constants) when the temperatures differ by more than `BISTABLE_K` and `BISTABLE_RELATIVE` of their mean, or one start finds no air).

    python scripts/gen1555_moist_limit_census.py draw --out before.json [--systems 100] [--stars G2V K5V] [--workers 4]
    python scripts/gen1555_moist_limit_census.py report before.json [after.json]

The systems are the sweep of the GEN.1403 report (`reports/gen1403-sweep.py`): star type, planets on, comets and binaries off, seed 3000 + index. Run `draw` on the build
before the fix and on the build after it and `report` both: a fix is good when the bistable share falls, no body sits on 375.0 K only because the rule pinned it, and the
temperate count (273 to 330 K with 10 kPa) does not fall for the reason of the other branch. The final settings are the last solve of each body (a body is solved once at the
birth age and once at the final age; the last one is stored).
"""

import argparse
import json
import math
import multiprocessing
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

STARTS_K = (120.0, 300.0, 600.0)
"""tuple: The probe starts: cold, mid and hot, so a cold and a warm branch both show."""
BISTABLE_K = 20.0
BISTABLE_RELATIVE = 0.10
"""float: A body is `bistable` when its answers differ by more than `BISTABLE_K` AND this share of their mean (the steam and lava airs of 1,000 K and more wander by tens of
kelvin inside one branch and are not a second branch)."""
NEAR_BOILING_K = 10.0
PIN_TOLERANCE_K = 0.05
AIR_PA = 1.0e4
TEMPERATE_K = (273.0, 330.0)


def _boiling_k(pressure_pa, saturation):
    """The temperature at which water's saturation pressure is `pressure_pa`, by bisection on [250, 700] K; None outside the table."""
    low, high = 250.0, 700.0
    try:
        if saturation("H2O", low) > pressure_pa or saturation("H2O", high) < pressure_pa:
            return None
        for _ in range(60):
            mid = 0.5 * (low + high)
            if saturation("H2O", mid) < pressure_pa:
                low = mid
            else:
                high = mid
    except (ValueError, KeyError, ArithmeticError):
        return None
    return 0.5 * (low + high)


def _one(job):
    star_type, index = job
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import surface_air
    from planetgen.util import draw

    real = surface_air.solve
    last = {}

    def recording(setting, start_temperature_k):
        planet = sys._getframe(1).f_locals.get("planet")
        if planet is not None:
            last[id(planet)] = (planet, setting, draw.getstate())        # the draw state the production solve starts from: its design constants are drawn from it
        return real(setting, start_temperature_k)

    surface_air.solve = recording
    try:
        config = SystemConfig()
        config.STAR_TYPE = star_type
        config.PLANETS = True
        config.COMETS = False
        config.BINARY_SYSTEM = False
        with draw.bound(3000 + index):
            system = StarSystem(system_config=config)
    finally:
        surface_air.solve = real
    rows = []
    kept = {id(p): (p, s, state) for p, s, state in last.values()}
    for planet in [b for p in system.planets for b in [p] + list(getattr(p, "moons", None) or [])]:
        if getattr(planet, "body_type", None) != "t" or id(planet) not in kept:
            continue
        temperature = getattr(planet, "surface_temperature", None)
        pressure = getattr(planet, "atmospheric_pressure", None) or 0.0
        if temperature is None:
            continue
        boiling = _boiling_k(pressure, surface_air.saturation_limit_pa) if pressure > 0.0 else None
        pinned = abs(temperature - surface_air.MOIST_LIMIT_K) <= PIN_TOLERANCE_K
        near = boiling is not None and abs(temperature - boiling) <= NEAR_BOILING_K
        row = dict(star=star_type, index=index, au=round(planet.distance, 4), moon=bool(getattr(planet, "is_moon", False)), t=temperature, p=pressure,
                   cls=getattr(planet, "planet_class", None), pinned=pinned, near_boiling=near, boiling_k=boiling, answers=None, bistable=None,
                   answers_production=None, bistable_production=None)
        if pressure >= AIR_PA or pinned or near:
            setting, state = kept[id(planet)][1], kept[id(planet)][2]

            def probe(production):
                answers = []
                for start in STARTS_K:
                    with draw.bound("gen1555-probe"):
                        if production:
                            draw.setstate(state)         # the production draw state: the constants the body was designed with, not a fresh stream's
                        air = real(setting, start)
                    answers.append(None if air is None else [round(air.temperature_k, 2), round(air.pressure_pa, 1)])
                return answers

            def is_bistable(answers):
                temperatures = [a[0] for a in answers if a is not None]
                gap = max(temperatures) - min(temperatures) if temperatures else 0.0
                mean = sum(temperatures) / len(temperatures) if temperatures else 1.0
                return (len(temperatures) < len(answers) and bool(temperatures)) or (gap > BISTABLE_K and gap > BISTABLE_RELATIVE * mean)

            row["answers"] = probe(False)
            row["bistable"] = is_bistable(row["answers"])
            row["answers_production"] = probe(True)
            row["bistable_production"] = is_bistable(row["answers_production"])
        rows.append(row)
    return rows


def command_draw(args):
    jobs = [(star, i) for star in args.stars for i in range(args.systems)]
    with multiprocessing.Pool(args.workers) as pool:
        parts = pool.map(_one, jobs, chunksize=1)
    rows = [r for part in parts for r in part]
    pathlib.Path(args.out).write_text(json.dumps(rows))
    print(f"{len(jobs)} systems, {len(rows)} bodies -> {args.out}")


def _summary(rows):
    out = {}
    for star in sorted({r["star"] for r in rows}):
        group = [r for r in rows if r["star"] == star]
        air = [r for r in group if r["p"] >= AIR_PA]
        probed = [r for r in group if r["bistable"] is not None]
        temperate = [r for r in air if TEMPERATE_K[0] <= r["t"] <= TEMPERATE_K[1]]
        out[star] = dict(bodies=len(group), air10kpa=len(air), temperate=len(temperate), pinned_at_375=sum(r["pinned"] for r in group),
                         pinned_bistable_production=sum(1 for r in group if r["pinned"] and r.get("bistable_production")),
                         near_boiling=sum(r["near_boiling"] for r in group), probed=len(probed), bistable=sum(1 for r in probed if r["bistable"]))
    return out


def command_report(args):
    for path in args.sweeps:
        rows = json.loads(pathlib.Path(path).read_text())
        print(f"\n{path}")
        for star, s in _summary(rows).items():
            share = s["bistable"] / s["probed"] if s["probed"] else math.nan
            print(f"  {star}: {s['bodies']} bodies, {s['air10kpa']} with 10 kPa, {s['temperate']} temperate; pinned at 375.0 K {s['pinned_at_375']}; "
                  f"of those bistable under the production draw state {s['pinned_bistable_production']}; near boiling {s['near_boiling']}; probed {s['probed']}, bistable {s['bistable']} ({share:.0%})")
        for r in rows:
            if r["pinned"] or r["near_boiling"] or r["bistable"]:
                tag = ",".join(k for k, v in (("pinned", r["pinned"]), ("boiling", r["near_boiling"]), ("bistable", r["bistable"])) if v)
                print(f"    {r['star']} seed {3000 + r['index']} {r['au']} au T {r['t']:.1f} K P {r['p'] / 1e5:.2f} bar class {r['cls']} [{tag}] answers {r['answers']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    d = sub.add_parser("draw")
    d.add_argument("--out", required=True)
    d.add_argument("--stars", nargs="+", default=["G2V", "K5V"])
    d.add_argument("--systems", type=int, default=100)
    d.add_argument("--workers", type=int, default=4)
    d.set_defaults(func=command_draw)
    r = sub.add_parser("report")
    r.add_argument("sweeps", nargs="+")
    r.set_defaults(func=command_report)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
