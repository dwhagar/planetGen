#!/usr/bin/env python3
"""
GEN.1448: the before-and-after census of the system-interaction changes (GEN.1446 gravity from every body, GEN.1449, GEN.1454
obliquity variation, GEN.1450 secular chain, GEN.1447 Kozai moons).

Each change moves some bodies' ocean tide, tidal heating, locking, spin, tilt or orbit. This script draws the same seeded systems
on two builds of the code and counts how many bodies moved in each quantity, so a lane can say what its change did to the
world rather than that its tests passed.

    python scripts/census_interactions.py draw --out before.json [--per-type 40] [--workers 4]     # on the build to measure
    python scripts/census_interactions.py compare before.json after.json [--report reports/x.md]   # after the change lands

Every system is seeded from its star type and its index (`draw.set_run_seed`), so the same index draws the same system on both
builds unless the change itself moves the draws; a system whose body count differs between the two builds is counted apart
(`reshaped`) and its bodies are not matched. A body is matched by (star type, system index, planet index, moon index).

The quantities per body are in `QUANTITIES`; a quantity moved when its two values differ by more than `TOLERANCE` (a relative
tolerance, or the state name changing for the locking state). Any attribute whose name holds `kozai` is carried as a count of
what the body lost to it, so GEN.1447's moon losses are counted as soon as that lands (a build without it records none).
"""

import argparse
import json
import multiprocessing
import os
import pathlib
import signal
import sys
import zlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

STAR_TYPES = ("M2V", "K5V", "G2V", "F5V", "A0V", "K0III")
"""tuple: The census stars, the same set as `tests/audit/prevalence.py`: four dwarfs, an A dwarf and a giant."""

SYSTEM_TIME_LIMIT_S = 240
"""int: Seconds one system may take before it is abandoned and counted."""

TOLERANCE = 0.01
"""float: A numeric quantity has moved when it differs from its other value by more than this fraction of the larger of the two."""

QUANTITIES = (
    ("ocean_tide_spring_m", "ocean tide range (spring, m)"),
    ("tidal_flux_w_m2", "tidal heating (W/m2)"),
    ("lock_state", "locking state (free, locked, resonant)"),
    ("rotation_period_hours", "spin period (h)"),
    ("axial_tilt_deg", "obliquity (deg)"),
    ("eccentricity", "eccentricity"),
)
"""tuple: (key in a body record, what the report calls it)."""


def seed_for(star_type, index):
    """The seed of one system, the same on every run and different for each."""
    return zlib.crc32(f"interactions:{star_type}:{index}".encode()) or 1


class _Abandoned(BaseException):
    """Raised by the alarm when one system ran past `SYSTEM_TIME_LIMIT_S`."""


def _number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if number == number else None


def _body_record(body):
    from planetgen.physics import spin

    rotation = _number(getattr(body, "rotation_period_hours", None))
    period = _number(getattr(body, "period", None))
    state = None
    if rotation and period:
        state = spin.spin_state(rotation, period)[0]
    record = {
        "ocean_tide_spring_m": _number(getattr(body, "tide_range_spring_m", None)),
        "tidal_flux_w_m2": _number(getattr(body, "tidal_flux_w_m2", None)),
        "lock_state": state,
        "rotation_period_hours": rotation,
        "axial_tilt_deg": _number(getattr(body, "axial_tilt_deg", None)),
        "eccentricity": _number(getattr(body, "eccentricity", None)),
    }
    lost = {}
    for name, value in vars(body).items():
        if "kozai" in name.lower():
            count = len(value) if isinstance(value, (list, tuple, set, dict)) else _number(value)
            if count:
                lost[name] = count
    if lost:
        record["kozai"] = lost
    return record


def system_record(system):
    """`{"planets": [{"body": record, "moons": [record, ...]}, ...]}` for the primary-star planets of one system, in their order."""
    planets = []
    for planet in system.all_planets:
        if planet.body_type == "a":
            continue
        planets.append({"body": _body_record(planet), "moons": [_body_record(moon) for moon in planet.moons]})
    return {"planets": planets}


def draw_type(star_type, count, first=0):
    """`count` systems of `star_type`, as a list of `{"star", "index", ...}` records (an abandoned or crashed one says so)."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    def _alarm(*_):
        raise _Abandoned()

    previous = signal.signal(signal.SIGALRM, _alarm)
    records = []
    try:
        for index in range(first, first + count):
            draw.set_run_seed(seed_for(star_type, index))
            config = SystemConfig()
            config.STAR_TYPE = star_type
            signal.alarm(SYSTEM_TIME_LIMIT_S)
            entry = {"star": star_type, "index": index}
            try:
                entry.update(system_record(StarSystem(system_config=config)))
            except _Abandoned:
                entry["abandoned"] = True
            except Exception as error:  # a crash is a finding, not the end of the census
                entry["error"] = f"{type(error).__name__}: {error}"
            finally:
                signal.alarm(0)
            records.append(entry)
    finally:
        signal.signal(signal.SIGALRM, previous)
    return records


def _job(job):
    return draw_type(*job)


def draw_all(per_type, star_types=STAR_TYPES, workers=None, chunk=5):
    """Every star type's systems, split into chunks over worker processes so the cores stay busy."""
    jobs = [(star, min(chunk, per_type - first), first) for star in star_types for first in range(0, per_type, chunk)]
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, len(jobs))) as pool:
        results = pool.map(_job, jobs, chunksize=1)
    return [record for result in results for record in result]


def _moved(key, a, b):
    if a is None and b is None:
        return False
    if a is None or b is None:
        return True
    if key == "lock_state":
        return a != b
    return abs(a - b) > TOLERANCE * max(abs(a), abs(b))


def _bodies(record):
    out = {}
    for pi, planet in enumerate(record["planets"]):
        out[(pi, None)] = planet["body"]
        for mi, moon in enumerate(planet["moons"]):
            out[(pi, mi)] = moon
    return out


def compare(before, after):
    """Per star type and overall: systems matched, reshaped, abandoned; bodies matched; bodies moved per quantity; Kozai losses."""
    index = {(r["star"], r["index"]): r for r in after}
    rows = {}
    for old in before:
        row = rows.setdefault(old["star"], {"systems": 0, "matched": 0, "reshaped": 0, "failed": 0, "bodies": 0, "moved_any": 0,
                                            "moved": {key: 0 for key, _ in QUANTITIES}, "of": {key: 0 for key, _ in QUANTITIES},
                                            "kozai_before": 0.0, "kozai_after": 0.0, "examples": {key: [] for key, _ in QUANTITIES}})
        row["systems"] += 1
        new = index.get((old["star"], old["index"]))
        if new is None or any(r.get(flag) for r in (old, new) for flag in ("abandoned", "error")):
            row["failed"] += 1
            continue
        old_bodies, new_bodies = _bodies(old), _bodies(new)
        for source, key in ((old_bodies, "kozai_before"), (new_bodies, "kozai_after")):
            row[key] += sum(sum(b.get("kozai", {}).values()) for b in source.values())
        if old_bodies.keys() != new_bodies.keys():
            row["reshaped"] += 1
            continue
        row["matched"] += 1
        for where, a in old_bodies.items():
            b = new_bodies[where]
            row["bodies"] += 1
            any_moved = False
            for key, _label in QUANTITIES:
                if a.get(key) is None and b.get(key) is None:
                    continue
                row["of"][key] += 1
                if _moved(key, a.get(key), b.get(key)):
                    row["moved"][key] += 1
                    any_moved = True
                    if len(row["examples"][key]) < 3:
                        row["examples"][key].append([old["star"], old["index"], list(where), a.get(key), b.get(key)])
            row["moved_any"] += any_moved
    return rows


def _share(n, of):
    return f"{n} of {of} ({n / of:.1%})" if of else "none"


def report(rows, before_name, after_name):
    lines = ["# System interaction census (GEN.1448)", "",
             f"Before: `{before_name}`. After: `{after_name}`. A quantity counts as moved when it differs by more than "
             f"{TOLERANCE:.0%} of the larger value (the locking state when the name changes). Bodies are matched by star type, "
             "system index, planet index and moon index; a system whose body count changed is `reshaped` and not matched.", ""]
    total = {"systems": 0, "matched": 0, "reshaped": 0, "failed": 0, "bodies": 0, "moved_any": 0, "moved": {k: 0 for k, _ in QUANTITIES},
             "of": {k: 0 for k, _ in QUANTITIES}, "kozai_before": 0.0, "kozai_after": 0.0}
    for row in rows.values():
        for key in ("systems", "matched", "reshaped", "failed", "bodies", "moved_any", "kozai_before", "kozai_after"):
            total[key] += row[key]
        for key, _ in QUANTITIES:
            total["moved"][key] += row["moved"][key]
            total["of"][key] += row["of"][key]
    lines += ["## Systems", "", "| Star | Systems | Matched | Reshaped | Failed | Bodies matched | Bodies with any change |", "|---|---|---|---|---|---|---|"]
    for star, row in {**rows, "all": total}.items():
        lines.append(f"| {star} | {row['systems']} | {row['matched']} | {row['reshaped']} | {row['failed']} | {row['bodies']} | {_share(row['moved_any'], row['bodies'])} |")
    lines += ["", "## Bodies moved, by quantity", "", "| Quantity | " + " | ".join(list(rows) + ["all"]) + " |", "|---|" + "---|" * (len(rows) + 1)]
    for key, label in QUANTITIES:
        cells = [_share(row["moved"][key], row["of"][key]) for row in rows.values()] + [_share(total["moved"][key], total["of"][key])]
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    lines += ["", "## Kozai losses (count of moons or bodies a build records as lost to the Kozai cycle)", "",
              f"Before: {total['kozai_before']:g}. After: {total['kozai_after']:g}.", "",
              "## Examples (star, system index, (planet, moon or none), before, after)", ""]
    for key, label in QUANTITIES:
        shown = [example for row in rows.values() for example in row["examples"][key]][:6]
        if shown:
            lines.append(f"- {label}: " + "; ".join(f"{e[0]} #{e[1]} {tuple(e[2])}: {e[3]} to {e[4]}" for e in shown))
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    d = sub.add_parser("draw", help="draw the seeded systems and save their records")
    d.add_argument("--out", required=True)
    d.add_argument("--per-type", type=int, default=40)
    d.add_argument("--workers", type=int, default=None)
    d.add_argument("--stars", nargs="*", default=list(STAR_TYPES))
    c = sub.add_parser("compare", help="compare two saved draws")
    c.add_argument("before")
    c.add_argument("after")
    c.add_argument("--report", default=None)
    args = parser.parse_args(argv)
    if args.command == "draw":
        records = draw_all(args.per_type, tuple(args.stars), args.workers)
        pathlib.Path(args.out).write_text(json.dumps(records), encoding="utf-8")
        failed = sum(1 for r in records if r.get("abandoned") or r.get("error"))
        print(f"{len(records)} systems drawn, {failed} abandoned or crashed -> {args.out}")
        return 0
    rows = compare(json.loads(pathlib.Path(args.before).read_text()), json.loads(pathlib.Path(args.after).read_text()))
    text = report(rows, args.before, args.after)
    if args.report:
        pathlib.Path(args.report).write_text(text, encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
