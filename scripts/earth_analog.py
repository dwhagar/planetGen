#!/usr/bin/env python3
"""
The Earth analog as a standing script (TEST.726, moved from reports/class-m-audit-2026-11-12/earth_analog.py): an exact Earth put through the real generator many times, and
what comes out held against a stored baseline, so a change to the air, water, life or class code shows its effect on one world we know the answer for.

    python scripts/earth_analog.py run BLOCK N [--age GY] [--water W] [--out FILE]   one JSON line per trial (the analysis form of the audit)
    python scripts/earth_analog.py baseline [--n 800] [--write]                      measure the baseline on main (the stored file is written only with --write)
    python scripts/earth_analog.py compare [--n 400] [--workers K]                   the numbers now against the baseline, for the train PR's body

The analog: a generated G2V system (seed 4001) supplies a rocky planet object; the star is set to the Sun (1 Msun, 3.828e26 W, 5772 K, age `--age`, default 4.6 Gy, snow line 2.7 AU) and
the planet to Earth (late accretion 0.7 percent, `impact_history.LATE_ACCRETED_PERCENT_EARTH`; 5.972e24 kg, 6371 km, 1 g, 5514 kg/m3, 1 AU, e 0.0167, metal 0.32, tilt 23.4, day 24 h). Each trial draws the bulk water, then calls
`planets.calculate_atmospheric_conditions`, `planets.apply_class`, `life.apply_life_data` and `planets.apply_class` as the system build does, on its own draw seed `BLOCK * 100000 + i`.
`--water` is `sample` (default: `delivered_water.sample_fraction` at 1 AU and a 2.7 AU snow line), a Bern class (`low`, `high`, `dry`, `rich`), or a fixed bulk water mass fraction.

The comparison is paired: trial i draws the same seed in the baseline and now, so with no change to the generator every trial comes out the same and the table shows no movement. It
is not a gate and never blocks a train; its numbers are what the items that move the Earth analog (GEN.1581, 1584, 1585, 1586) quote before and after. The baseline file is
`scripts/data/earth_analog_baseline.json`; a train that changes the numbers on purpose renews it with `baseline --write`.
"""

import argparse
import collections
import json
import multiprocessing
import os
import statistics
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

BASELINE_PATH = os.path.join(ROOT, "scripts", "data", "earth_analog_baseline.json")
BLOCK = 1
AGE_GY = 4.6
DISTANCE_STEP_AU = 1e-9
"""float: How far each trial moves the planet's orbit (AU), so its oddity streams differ from trial to trial (GEN.1646); about 1e-9 of the orbit, far under anything the air or the climate reads."""
CHUNK = 25
"""int: Trials per work unit when the run is spread over workers."""

GROUPS = {"gravity": ["gravity_g"], "pressure": ["pressure_pa"], "O2": ["o2_partial_pa"], "CO2": ["co2_partial_pa"], "gases": ["reducing_share", "share_h2", "share_so2"],
          "T_mean": ["temperature_k"], "T_extremes": ["hottest_region_k", "coldest_region_k", "night_temperature_k", "wet_bulb_c"], "cover": ["ocean_fraction"],
          "water_phase": ["phase_water"]}
"""dict: The conditions of Class M, by the name of the group a failing one is counted under."""


def _setup(age_gy):
    """The Earth analog planet and Class M, built once per process."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import class_fit, impact_history
    from planetgen.util import draw

    config = SystemConfig()
    config.PLANETS = True
    config.STAR_TYPE = "G2V"
    with draw.bound(4001):
        system = StarSystem(system_config=config)
    star = system.stars[0]
    planet = [body for body in system.all_planets if body.body_type == "t" and not body.is_moon][0]
    star.age, star.mass, star.luminosity, star.temperature = age_gy, 1.989e30, 3.828e26, 5772
    for name, value in (("radius", 6.957e5), ("snow_line_au", 2.7)):
        if hasattr(star, name):
            setattr(star, name, value)
    planet.mass, planet.radius, planet.gravity, planet.density, planet.distance, planet.eccentricity = 5.972e24, 6371.0, 1.0, 5514, 1.0, 0.0167
    planet.ice_mass_fraction = planet.envelope_mass_fraction = 0.0
    planet.primary_mass_kg, planet._star_distance_au, planet.formation_distance_au = star.mass, 1.0, 1.0
    planet.rotation_period_hours, planet.axial_tilt_deg = 24.0, 23.4
    planet._late_accreted_percent = impact_history.LATE_ACCRETED_PERCENT_EARTH     # the generated system's own belts set 3.0, six times the Earth's veneer (0.7, see LATE_ACCRETED_PERCENT_EARTH; [R4598] gives 0.5 to 0.8)
    return planet, [c for c in class_fit.CLASSES if c.letter == "M"][0]


def _water(water):
    from planetgen.physics import delivered_water
    from planetgen.util import draw
    from planetgen.util.random import log_uniform

    if water == "sample":
        return delivered_water.sample_fraction(1.0, 1.0 / 2.7)
    if water in delivered_water.CLASSES:
        low, high = delivered_water.bracket(water)
        return log_uniform(low, high, rng=draw)
    return float(water)


def trials(block, start, count, age_gy=AGE_GY, water="sample", full=False):
    """Trials `start` to `start + count - 1` of `block` as rows. A row holds the water, the class, the hydrosphere, the temperature (K), the pressure (Pa), the ocean cover, the life
    class and the groups of Class M's conditions that fail; `full` adds the classifier's quantities (the form of `run`). A trial that raises holds `error`."""
    from planetgen.generation import life
    from planetgen.physics import atmosphere, class_fit, planets
    from planetgen.util import draw

    planet, class_m = _setup(age_gy)
    to_group = {quantity: group for group, quantities in GROUPS.items() for quantity in quantities}
    rows = []
    for i in range(start, start + count):
        with draw.bound(block * 100000 + i):
            planet.distance = 1.0 + i * DISTANCE_STEP_AU      # GEN.1646: the oddities' streams are keyed on the body (`oddities.body_key`: star, distance, mass), so one reused planet drew every trial's second climate branch from one stream
            planet.water_mass_fraction = _water(water)
            planet.__dict__.pop("_formation_fractions", None)
            planet.metal_mass_fraction, planet.ice_mass_fraction, planet.envelope_mass_fraction = 0.32, 0.0, 0.0
            planet.rock_mass_fraction = 1 - 0.32 - planet.water_mass_fraction
            row = {"i": i, "w": planet.water_mass_fraction}
            try:
                planets.calculate_atmospheric_conditions(planet)
                planets.apply_class(planet)
                life.apply_life_data(planet)
                planets.apply_class(planet)
            except Exception as error:      # a trial that raises is a finding, kept in the row
                row["error"] = f"{type(error).__name__}: {error}"[:120]
                rows.append(row)
                continue
        found = class_fit._found(planet)
        chosen = class_fit.classify(planet)
        failing = sorted({to_group[c["quantity"]] for c in class_m.when if class_fit.membership(c, found) < 1.0})
        row.update(cls=chosen.letter, hyd=planet.hydrosphere, T=round(planet.surface_temperature, 1), P=round(planet.atmospheric_pressure, 1), oc=round(planet.ocean_fraction, 3),
                   life_class=planet.life_class, fails=failing)
        if full:
            row.update(match=chosen.match, runner=chosen.runner_up, stage=planet.life_stage, o2=found["o2_partial_pa"], co2=found["co2_partial_pa"], g=planet.gravity,
                       hot=found["hottest_region_k"], cold=found["coldest_region_k"], night=found["night_temperature_k"], wetbulb=found["wet_bulb_c"], albedo=planet.bond_albedo,
                       G=planet.greenhouse_factor, gas={k: round(v, 1) for k, v in atmosphere.gas_partials_pa(planet).items() if v >= 1.0})
        rows.append(row)
    return rows


def _job(job):
    return trials(*job)


def measure(n, workers=None, block=BLOCK, age_gy=AGE_GY, water="sample"):
    """The first `n` trials of `block`, over `workers` processes (default: the cores), in trial order."""
    jobs = [(block, start, min(CHUNK, n - start), age_gy, water) for start in range(0, n, CHUNK)]
    if len(jobs) == 1 or workers == 1:
        parts = [_job(job) for job in jobs]
    else:
        with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, len(jobs))) as pool:
            parts = pool.map(_job, jobs, chunksize=1)
    return [row for part in parts for row in part]


def summarize(rows):
    """The numbers of a set of rows: the trials, the classes, the hydrospheres, the temperature bins, the medians and how often each group of Class M's conditions fails."""
    good = [r for r in rows if "error" not in r]
    temperature_bin = lambda t: "<258 K" if t < 258 else "258-343 K" if t <= 343 else ">343 K"   # noqa: E731
    return {"trials": len(rows), "errors": len(rows) - len(good),
            "classes": dict(collections.Counter(r["cls"] for r in good)), "hydrosphere": dict(collections.Counter(r["hyd"] for r in good)),
            "temperature": dict(collections.Counter(temperature_bin(r["T"]) for r in good)),
            "median_pressure_kpa": round(statistics.median(r["P"] for r in good) / 1000.0) if good else None,
            "median_temperature_k": round(statistics.median(r["T"] for r in good)) if good else None,
            "m_fails": {group: sum(1 for r in good if group in r["fails"]) for group in GROUPS}}


def git_head():
    done = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=False)
    return done.stdout.strip() or "unknown"


def load_baseline(path=None):
    with open(path or BASELINE_PATH, encoding="utf-8") as handle:
        return json.load(handle)


def compare_text(baseline, rows):
    """The Markdown for the train PR's body: the baseline's first `len(rows)` trials against `rows`, then the trials whose class changed."""
    n = len(rows)
    then_rows = [r for r in baseline["rows"] if r["i"] < n]
    then, now = summarize(then_rows), summarize(rows)

    def counts(key, labels):
        return [(label, then[key].get(label, 0), now[key].get(label, 0)) for label in labels]

    lines = [f"Earth analog: {n} trials (seeds {baseline['block']}xxxxx, star age {baseline['age_gy']} Gy), baseline made on {baseline['made_on']}; paired draws, so no change to the generator moves nothing.", "",
             "| measure | baseline | now |", "|---|---|---|"]
    for key, title in (("classes", "class"), ("hydrosphere", "hydrosphere"), ("temperature", "surface temperature")):
        labels = sorted(set(then[key]) | set(now[key]), key=lambda label: -(then[key].get(label, 0) + now[key].get(label, 0)))
        for label, before, after in counts(key, labels):
            lines.append(f"| {title} {label} | {before} | {after}{'' if before == after else f' ({after - before:+d})'} |")
    lines.append(f"| median pressure kPa | {then['median_pressure_kpa']} | {now['median_pressure_kpa']} |")
    lines.append(f"| median temperature K | {then['median_temperature_k']} | {now['median_temperature_k']} |")
    for group in GROUPS:
        before, after = then["m_fails"][group], now["m_fails"][group]
        if before or after:
            lines.append(f"| trials failing Class M's {group} | {before} | {after}{'' if before == after else f' ({after - before:+d})'} |")
    if now["errors"] != then["errors"]:
        lines.append(f"| trials that raised | {then['errors']} | {now['errors']} |")
    before_by_i = {r["i"]: r for r in then_rows}
    moves = collections.Counter((before_by_i[r["i"]]["cls"], r["cls"]) for r in rows
                                if "error" not in r and r["i"] in before_by_i and "error" not in before_by_i[r["i"]] and before_by_i[r["i"]]["cls"] != r["cls"])
    lines.append("")
    lines.append("Trials that changed class: " + (", ".join(f"{a} to {b} {count}" for (a, b), count in moves.most_common()) if moves else "none") + ".")
    return "\n".join(lines)


def report(n=400, workers=None):
    """`compare_text` for the stored baseline, measured now; a missing baseline or any failure is a sentence, never an exception (this never blocks a train)."""
    try:
        baseline = load_baseline()
        return compare_text(baseline, measure(min(n, len(baseline["rows"])), workers, baseline["block"], baseline["age_gy"], baseline["water"]))
    except Exception as error:
        return f"Earth analog: not run ({type(error).__name__}: {error})"


def command_run(args):
    for row in trials(args.block, 0, args.n, args.age, args.water, full=True):
        print(json.dumps(row), flush=True)


def command_baseline(args):
    rows = measure(args.n, args.workers)
    document = {"made_on": git_head(), "block": BLOCK, "age_gy": AGE_GY, "water": "sample", "rows": rows}
    print(json.dumps(summarize(rows), indent=1))
    if args.write:
        os.makedirs(os.path.dirname(BASELINE_PATH), exist_ok=True)
        with open(BASELINE_PATH, "w", encoding="utf-8") as handle:
            json.dump(document, handle, separators=(",", ":"))
            handle.write("\n")
        print(f"wrote {BASELINE_PATH}")


def command_compare(args):
    print(report(args.n, args.workers))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run")
    run.add_argument("block", type=int)
    run.add_argument("n", type=int)
    run.add_argument("--age", type=float, default=AGE_GY)
    run.add_argument("--water", default="sample")
    run.set_defaults(func=command_run)
    base = sub.add_parser("baseline")
    base.add_argument("--n", type=int, default=800)
    base.add_argument("--workers", type=int, default=None)
    base.add_argument("--write", action="store_true")
    base.set_defaults(func=command_baseline)
    compare = sub.add_parser("compare")
    compare.add_argument("--n", type=int, default=400)
    compare.add_argument("--workers", type=int, default=None)
    compare.set_defaults(func=command_compare)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
