# tests/audit/late_supply_census.py

"""
The late supply census (GEN.1493 slice 0): every seed is generated as a whole system and each rocky planet's veneer is worked out from the system's own belts and planets
(`physics/late_supply.py`) and set against the mean the present law gives it (`impact_history.late_accreted_percent_mean`, GEN.739).

It reports, for rocky planets of 0.3 to 3 Earth masses of single stars, the share that receive nothing (no belt reaches them), the median and the spread of the ratio
(supply at INWARD_SHARE 1) / (the law's mean), and from it the INWARD_SHARE at which the median planet takes the law's mean. That number is the one fitted constant (row LS1, tweak L3).
The census reads, it does not wire: nothing stored or drawn changes.
"""

import collections
import multiprocessing
import os
import statistics

SEEDS = 480
REPORT_NAME = "late-supply-census-latest.md"
EARTH_KG = 5.9722e24


def system_rows(seed):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import impact_history, influx, late_supply
    from planetgen.util import draw

    draw.set_run_seed(seed)
    system = StarSystem(system_config=SystemConfig())
    if len(system.stars) > 1:
        return []
    star_mass = system.star.mass / 1.98841e30
    bodies = [(body.radius, body.distance, body.mass) for body in system.planets if body.body_type in ("t", "g") and body.mass and body.radius]
    reservoirs = [influx.Reservoir("belt", belt.primordial_mass_earth * EARTH_KG, belt.lower_limit, belt.upper_limit)
                  for belt in system.planets if belt.body_type == "a" and belt.primordial_mass_earth]
    rows = []
    for body in system.planets:
        if body.body_type != "t" or not body.radius or not body.mass:
            continue
        mass_earth = body.mass / EARTH_KG
        if not (0.3 <= mass_earth <= 3.0):
            continue
        mean = impact_history.late_accreted_percent_mean(mass_earth) / 100.0
        supply = late_supply.veneer_fraction((body.radius, body.distance, body.mass), bodies, star_mass, reservoirs) / late_supply.INWARD_SHARE
        giants = any(other.body_type == "g" for other in system.planets)
        rows.append((seed, mass_earth, body.distance, bool(reservoirs), giants, supply / mean))
    return rows


def run(seeds=SEEDS, workers=None):
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(system_rows, range(1, seeds + 1), chunksize=1)
    return [row for part in parts for row in part]


def summarize(rows):
    """`{group: (count, none, median, p10, p90)}` of the ratio (supply at share 1 over the law's mean) for all planets, and for those with and without a giant."""
    groups = {"all": rows, "with a giant": [r for r in rows if r[4]], "without a giant": [r for r in rows if not r[4]]}
    out = {}
    for name, part in groups.items():
        ratios = sorted(r[5] for r in part)
        taking = [x for x in ratios if x > 0.0]
        if taking:
            q = statistics.quantiles(taking, n=10) if len(taking) >= 10 else [min(taking)] * 9
            out[name] = (len(part), len(part) - len(taking), statistics.median(taking), q[0], q[-1])
        else:
            out[name] = (len(part), len(part), 0.0, 0.0, 0.0)
    return out


def format_report(rows, head):
    table = summarize(rows)
    lines = [f"# Late supply census (GEN.1493 slice 0), head {head}", "",
             "Rocky planets of 0.3 to 3 Earth masses of single stars, the supply of the system's own belts at INWARD_SHARE 1 over the mean of the present law.", "",
             "| group | planets | taking nothing | median ratio | 10th percentile | 90th percentile |", "|---|---|---|---|---|---|"]
    for name, (count, none, median, low, high) in table.items():
        lines.append(f"| {name} | {count} | {none} | {median:.4g} | {low:.4g} | {high:.4g} |")
    median = table["all"][2]
    lines += ["", f"INWARD_SHARE at which the median planet that takes anything takes the law's mean: {1.0 / median:.3g}" if median > 0 else "No planet took anything."]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    import pathlib
    import subprocess

    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    report = format_report(run(), head)
    out = pathlib.Path(os.environ.get("PLANETGEN_REPORT_DIR", ".")) / REPORT_NAME
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(report)
    print(report)
