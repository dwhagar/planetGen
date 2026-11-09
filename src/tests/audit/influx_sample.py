# tests/audit/influx_sample.py

"""
A sample of the belt dust and comet influx of rocky planets (GEN.1440, GEN.1453 slice A, `physics/influx.py`): whole seeded systems are generated, and each planet that is not a giant is given the influx of the
belts and comets of its system, with and without the shielding of the giants in the way. The sample shows that the influx varies with the system's belts (none, an inner belt, a debris disk's outer belt), and how
many planets a giant shields.
"""

import collections
import multiprocessing
import os
import pathlib
import statistics

SEEDS = 120
REPORT_NAME = "influx-sample-latest.md"


def system_rows(seed):
    """`[(kind, influx, unshielded, belts, giants)]` for every non-giant planet of the system of `seed`; `kind` is `outer` if the system has an outer-disk belt, `inner` if only other belts, else `none`."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import constants, influx
    from planetgen.util import draw

    draw.set_run_seed(seed)
    system = StarSystem(system_config=SystemConfig())
    star = system.star
    bodies = list(system.planets)
    belts = [b for b in bodies if b.body_type == 'a' and b.mass_kg]
    reservoirs = [influx.Reservoir("belt", b.mass_kg, b.lower_limit, b.upper_limit, star.age * 1.0e9) for b in belts]
    reservoirs += [influx.Reservoir("comet", c.mass_kg, c.perihelion_distance_au, c.perihelion_distance_au) for c in system.comets if c.mass_kg]
    giants = [influx.Giant(b.mass / constants.EARTH_MASS_TO_KG, b.radius / 6371.0, b.distance) for b in bodies if b.body_type == 'g']
    kind = "outer" if any(getattr(b, "origin", None) == "outer-disk" for b in belts) else ("inner" if belts else "none")
    rows = []
    for planet in bodies:
        if planet.body_type in ('a', 'g'):
            continue
        mass_solar = star.mass / constants.SOLAR_MASS_TO_KG
        shielded = influx.planet_influx_kg_per_yr(planet.radius, planet.distance, mass_solar, reservoirs, giants)[0]
        bare = influx.planet_influx_kg_per_yr(planet.radius, planet.distance, mass_solar, reservoirs, [])[0]
        rows.append((kind, shielded, bare, len(belts), len(giants)))
    return rows


def run(seeds=SEEDS, workers=None):
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(system_rows, range(1, seeds + 1), chunksize=1)
    return [row for part in parts for row in part]


def format_report(rows, header=""):
    lines = ["# Belt dust and comet influx sample (GEN.1440, GEN.1453)", "", header, "",
             "| system kind | planets | with any influx | median kg a year | 90th percentile kg a year | shielded by a giant (share left, median) |", "|---|---|---|---|---|---|"]
    by_kind = collections.defaultdict(list)
    for row in rows:
        by_kind[row[0]].append(row)
    for kind in ("none", "inner", "outer"):
        group = by_kind.get(kind, [])
        positive = sorted(row[1] for row in group if row[1] > 0.0)
        median = statistics.median(positive) if positive else 0.0
        top = positive[int(0.9 * (len(positive) - 1))] if positive else 0.0
        left = sorted(row[1] / row[2] for row in group if row[2] > 0.0 and row[1] < 0.999 * row[2])
        shielded = f"{len(left)} ({statistics.median(left):.2f})" if left else "0"
        lines.append(f"| {kind} belts | {len(group)} | {len(positive)} | {median:.3g} | {top:.3g} | {shielded} |")
    return "\n".join(lines) + "\n"


def write_report(text, directory=None, name=REPORT_NAME):
    directory = pathlib.Path(directory or os.environ.get("PLANETGEN_REPORT_DIR") or ".test-reports")
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(text, encoding="utf-8")
    return path


if __name__ == "__main__":
    print(write_report(format_report(run(), f"Seeds 1 to {SEEDS}, rocky planets of whole systems.")))
