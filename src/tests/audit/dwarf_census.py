# tests/audit/dwarf_census.py

"""
The dwarf census (GEN.1087, docs/design/systems/belt-origins.md): how many Class D bodies, dwarf-size bodies, small moons and belts the
generator makes per system, per census star type. It measures; it asserts nothing. The numbers seed the test targets of the dwarf population
and scatter modules and say whether a system needs more Class D bodies than it makes.

A *dwarf-size* body is one of 400 to 2,500 km radius that is not a giant (the dwarf planets of the Solar System run from Ceres, 470 km, to
Pluto, 1,190 km, and Eris, 1,163 km); a *small moon* is a moon of 50 to 500 km radius (the lower figure of `planet_classes.json`'s older wording,
`mass-object-class.md`). A body counts as Class D when `planet_class` is `D`, whether it is a planet, a moon or a belt member.

    PYTHONPATH=. python3 -m tests.audit.dwarf_census --count 60 --out /tmp/dwarf-census.md
"""

import argparse
import json
import math
import multiprocessing
import os
import signal

from tests.audit import prevalence as census

SYSTEM_TIME_LIMIT_S = 240
"""int: Seconds one system may take before it is abandoned and counted."""

STAR_TYPES = census.STAR_TYPES + ("M5V",)
"""tuple: The census star types, and a late M dwarf."""

DWARF_RADIUS_KM = (400.0, 2500.0)
SMALL_MOON_RADIUS_KM = (50.0, 500.0)
"""tuple: The radius ranges, km, of a dwarf-size body and of a small moon."""

RADIUS_BINS_KM = (0.0, 50.0, 200.0, 500.0, 1000.0, 2500.0, math.inf)
"""tuple: The bin edges of the radius histogram of every non-giant body, km."""


class _Abandoned(BaseException):
    """Raised by the alarm when one system ran past `SYSTEM_TIME_LIMIT_S`."""


def _record(system):
    """One system's counts: planets, moons, belts, Class D bodies by kind, small moons, dwarf-size bodies, and the radius of each moon and planet."""
    planets = [planet for planet in system.all_planets if planet.body_type != "a"]
    belts = [planet for planet in system.all_planets if planet.body_type == "a"]
    moons = [moon for planet in planets for moon in planet.moons]
    solid = [body for body in (*planets, *moons) if body.body_type != "g"]
    radii = [round(body.radius, 1) for body in solid if getattr(body, "radius", None)]
    return {
        "planets": len(planets), "moons": len(moons), "belts": len(belts),
        "belt_widths_au": [round(belt.upper_limit - belt.lower_limit, 3) for belt in belts],
        "class_d_planets": sum(1 for p in planets if getattr(p, "planet_class", None) == "D"),
        "class_d_moons": sum(1 for m in moons if getattr(m, "planet_class", None) == "D"),
        "small_moons": sum(1 for m in moons if SMALL_MOON_RADIUS_KM[0] <= (m.radius or 0.0) <= SMALL_MOON_RADIUS_KM[1]),
        "dwarf_size": sum(1 for r in radii if DWARF_RADIUS_KM[0] <= r <= DWARF_RADIUS_KM[1]),
        "radii_km": radii,
        "planet_radii_km": [round(p.radius, 1) for p in planets if p.body_type != "g" and getattr(p, "radius", None)],
        "moon_radii_km": [round(m.radius, 1) for m in moons if getattr(m, "radius", None)],
        "class_d_radii_km": [round(b.radius, 1) for b in (*planets, *moons) if getattr(b, "planet_class", None) == "D" and getattr(b, "radius", None)],
        "giants": sum(1 for p in planets if p.body_type == "g"),
        "classes": sorted({getattr(body, "planet_class", "?") for body in solid}),
        "binary": system.binary_type is not None,
    }


def draw_type(star_type, count):
    """`count` systems round `star_type` as records; a system over the time limit is `{"abandoned": True}`, a crash `{"error": ...}`."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    draw.set_run_seed(census.seed_for("dwarf-census", star_type))
    records = []

    def _alarm(*_):
        raise _Abandoned()

    previous = signal.signal(signal.SIGALRM, _alarm)
    try:
        for _ in range(count):
            config = SystemConfig()
            config.STAR_TYPE = star_type
            signal.alarm(SYSTEM_TIME_LIMIT_S)
            try:
                records.append(_record(StarSystem(system_config=config)))
            except _Abandoned:
                records.append({"abandoned": True})
            except Exception as error:
                records.append({"error": f"{type(error).__name__}: {error}"})
            finally:
                signal.alarm(0)
    finally:
        signal.signal(signal.SIGALRM, previous)
    return [{"star": star_type, **record} for record in records]


def _job(job):
    return draw_type(*job)


def draw_types(count, star_types=STAR_TYPES, workers=None):
    jobs = [(star, count) for star in star_types]
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, len(jobs))) as pool:
        results = pool.map(_job, jobs, chunksize=1)
    return [record for result in results for record in result]


def summarize(records):
    """`{star: {n, abandoned, errors, per-system means and shares}}`."""
    out = {}
    for star in sorted({record["star"] for record in records}):
        mine = [r for r in records if r["star"] == star]
        good = [r for r in mine if "planets" in r]
        n = len(good) or 1
        mean = lambda key: sum(r[key] for r in good) / n   # noqa: E731
        share = lambda key: sum(1 for r in good if r[key] > 0) / n   # noqa: E731
        bins = [0] * (len(RADIUS_BINS_KM) - 1)
        for r in good:
            for radius in r["radii_km"]:
                for i in range(len(bins)):
                    if RADIUS_BINS_KM[i] <= radius < RADIUS_BINS_KM[i + 1]:
                        bins[i] += 1
        def histogram(key):
            counts = [0] * (len(RADIUS_BINS_KM) - 1)
            for r in good:
                for radius in r.get(key, ()):
                    for i in range(len(counts)):
                        if RADIUS_BINS_KM[i] <= radius < RADIUS_BINS_KM[i + 1]:
                            counts[i] += 1
            return [c / n for c in counts]

        out[star] = {
            "planet_bins": histogram("planet_radii_km"), "moon_bins": histogram("moon_radii_km"), "class_d_bins": histogram("class_d_radii_km"),
            "giants": mean("giants"),
            "n": len(good), "abandoned": sum(1 for r in mine if r.get("abandoned")), "errors": [r["error"] for r in mine if "error" in r][:3],
            "planets": mean("planets"), "moons": mean("moons"), "belts": mean("belts"),
            "belts_share": share("belts"), "class_d_planets": mean("class_d_planets"), "class_d_moons": mean("class_d_moons"),
            "class_d_share": sum(1 for r in good if r["class_d_planets"] + r["class_d_moons"] > 0) / n,
            "small_moons": mean("small_moons"), "small_moons_share": share("small_moons"),
            "dwarf_size": mean("dwarf_size"), "dwarf_size_share": share("dwarf_size"),
            "belt_width_au": (sum(w for r in good for w in r["belt_widths_au"]) / max(1, sum(len(r["belt_widths_au"]) for r in good))),
            "radius_bins": [b / n for b in bins],
            "classes": sorted({c for r in good for c in r["classes"]}),
        }
    return out


def format_report(summary, count):
    lines = [f"# Dwarf census ({count} systems per star type requested)", "",
             "| star | n | planets | moons | belts | share with a belt | Class D planets | Class D moons | share with a Class D | small moons (50 to 500 km) | share with one | dwarf-size (400 to 2500 km) | share with one | mean belt width AU |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for star, v in summary.items():
        lines.append(f"| {star} | {v['n']} ({v['abandoned']} abandoned) | {v['planets']:.2f} | {v['moons']:.2f} | {v['belts']:.2f} | {v['belts_share']:.2f} | "
                     f"{v['class_d_planets']:.3f} | {v['class_d_moons']:.3f} | {v['class_d_share']:.2f} | {v['small_moons']:.2f} | {v['small_moons_share']:.2f} | "
                     f"{v['dwarf_size']:.3f} | {v['dwarf_size_share']:.2f} | {v['belt_width_au']:.2f} |")
    lines += ["", "Radius histogram, bodies per system (non-giant planets and moons), bins in km " + ", ".join(str(b) for b in RADIUS_BINS_KM), ""]
    for star, v in summary.items():
        lines.append(f"- {star}: all " + ", ".join(f"{b:.2f}" for b in v["radius_bins"]) + f"; classes seen {''.join(v['classes'])}")
        lines.append("  - planets only: " + ", ".join(f"{b:.2f}" for b in v["planet_bins"]) + f" (giants {v['giants']:.2f}); moons only: "
                     + ", ".join(f"{b:.2f}" for b in v["moon_bins"]) + "; Class D bodies: " + ", ".join(f"{b:.2f}" for b in v["class_d_bins"]))
        if v["errors"]:
            lines.append(f"  - errors: {v['errors']}")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--count", type=int, default=60, help="systems per star type")
    parser.add_argument("--out", default="/tmp/dwarf-census.md")
    parser.add_argument("--json", default=None, help="also write the raw records here")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args(argv)
    records = draw_types(args.count, workers=args.workers)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump(records, handle)
    report = format_report(summarize(records), args.count)
    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write(report)
    print(report)


if __name__ == "__main__":
    main()
