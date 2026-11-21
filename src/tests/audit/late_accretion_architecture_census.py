# tests/audit/late_accretion_architecture_census.py

"""
The late accretion census by architecture (GEN.1493 slice 3, read-only): every seed is generated as a whole system and each rocky planet's late accretion percent, the one the pass set
(`generation/late_accretion.py`), is set against the mean the law gives its mass (`impact_history.late_accreted_percent_mean`), grouped by the architecture of the planet's own star:

* `no giant`: no giant planet of the star;
* `giant inside the belt zone`: a giant, none beyond 10 AU (hot and warm giants);
* `giant beyond 10 AU`: the outermost giant at 10 AU or more (a Jupiter or a Neptune among them);
* `with a companion star`: the system has more than one star (the planets of the primary are counted), whatever its giants;
* `comets and no belt`: the star has comets and no belt (a split of the groups above, shown beside them).

The baseline is main before slice 2 (the inward share read from scattering, stirring and comets). Run it again after slice 2 to see what moved. Nothing here changes any generated value.
"""

import multiprocessing
import os
import statistics

SEEDS = 240
REPORT_NAME = "late-accretion-architecture-census-latest.md"  # PLANETGEN_REPORT_NAME overrides it
EARTH_KG = 5.9722e24
GIANT_FAR_AU = 10.0
GROUPS = ("no giant", "giant inside the belt zone", "giant beyond 10 AU", "with a companion star", "comets and no belt")


def groups_of(star_count, giant_au, has_belt, has_comets):
    """The groups one planet falls in: its architecture, and `comets and no belt` beside it when that holds."""
    if star_count > 1:
        found = ["with a companion star"]
    elif giant_au is None:
        found = ["no giant"]
    elif giant_au >= GIANT_FAR_AU:
        found = ["giant beyond 10 AU"]
    else:
        found = ["giant inside the belt zone"]
    if has_comets and not has_belt:
        found.append("comets and no belt")
    return found


def system_rows(seed, star_type=None):
    """`[(group, mass_earth, percent, mean_percent, floored, capped, no_belt_reaches)]` of the rocky planets of 0.3 to 3 Earth masses of the primary star of the system made from `seed` (one row for each group the planet is in)."""
    from planetgen.generation import late_accretion
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import impact_history, influx
    from planetgen.util import draw

    draw.set_run_seed(seed)
    config = SystemConfig()
    if star_type:
        config.STAR_TYPE = star_type
        config.PLANETS = True
    system = StarSystem(system_config=config)
    star = system.star
    bodies = list(system.all_planets)
    giants = [body.distance for body in bodies if body.body_type == "g" and not body.is_moon and body.star is star]
    has_belt = any(body.body_type == "a" for body in bodies)
    has_comets = bool(system.comets)
    far = max(giants) if giants else None
    rows = []
    for body in bodies:
        if body.body_type != "t" or body.is_moon or body.star is not star:
            continue
        percent = body.build.late_accreted_percent
        if percent is None or not body.mass:
            continue
        mass_earth = body.mass / EARTH_KG
        if not 0.3 <= mass_earth <= 3.0:
            continue
        mean = impact_history.late_accreted_percent_mean(mass_earth)
        floored = percent <= late_accretion.FLOOR_SHARE * mean * (1.0 + 1e-9)
        capped = percent >= impact_history.LATE_ACCRETED_PERCENT_MAX * (1.0 - 1e-9)
        unreached = not any(influx.reach(influx.Reservoir("belt", 1.0, belt.lower_limit, belt.upper_limit), body.distance) > 0.0
                            for belt in bodies if belt.body_type == "a" and late_accretion._same_star(belt, body) and late_accretion._primordial_kg(belt) > 0.0)
        for group in groups_of(len(system.stars), far, has_belt, has_comets):
            rows.append((group, mass_earth, percent, mean, floored, capped, unreached))
    return rows


def _one(args):
    return system_rows(*args)


def run(seeds=SEEDS, workers=None, first_seed=1, star_type=None):
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(_one, [(seed, star_type) for seed in range(first_seed, first_seed + seeds)], chunksize=1)
    return [row for part in parts for row in part]


def summarize(rows):
    """`{group: {planets, floor, cap, median, p10, p90, median_percent}}` from the rows `system_rows` gives."""
    found = {}
    for group in GROUPS:
        mine = [row for row in rows if row[0] == group]
        ratios = sorted(row[2] / row[3] for row in mine)
        n = len(ratios)

        def pick(share, ratios=ratios, n=n):
            return ratios[min(n - 1, int(share * n))] if n else float("nan")

        found[group] = {"planets": n, "floor": sum(1 for row in mine if row[4]), "unreached": sum(1 for row in mine if row[6]), "cap": sum(1 for row in mine if row[5]),
                        "median": statistics.median(ratios) if n else float("nan"), "p10": pick(0.10), "p90": pick(0.90),
                        "median_percent": statistics.median(row[2] for row in mine) if n else float("nan")}
    return found


def format_report(rows, heading):
    found = summarize(rows)
    lines = [f"# {heading}", "", f"Rocky planets of 0.3 to 3 Earth masses of the primary star, {SEEDS} seeds. Ratio is the late accretion over the law's mean for the planet's mass.", "",
             "| Architecture | Planets | Floor | of them no belt outside | Ceiling | Median percent | Median ratio | 10th and 90th ratio |", "|---|---|---|---|---|---|---|---|"]
    for group in GROUPS:
        s = found[group]
        n = max(s["planets"], 1)
        lines.append(f"| {group} | {s['planets']} | {100.0 * s['floor'] / n:.0f} percent | {100.0 * s['unreached'] / n:.0f} percent of planets | {100.0 * s['cap'] / n:.0f} percent | {s['median_percent']:.3g} | {s['median']:.3g} | {s['p10']:.3g} and {s['p90']:.3g} |")
    lines.append("")
    return "\n".join(lines)


def main():
    """Environment: `ARCH_STAR_TYPE` (for example G2V, a Sun-like star: giants and comets are rare around the default stars), `ARCH_FIRST_SEED`, `ARCH_WORKERS` and `ARCH_LABEL` (the heading's label)."""
    star_type = os.environ.get("ARCH_STAR_TYPE") or None
    text = format_report(run(first_seed=int(os.environ.get("ARCH_FIRST_SEED", "1")), star_type=star_type, workers=int(os.environ.get("ARCH_WORKERS", "0")) or None),
                         f"The late accretion census by architecture (GEN.1493 slice 3, {os.environ.get('ARCH_LABEL', 'baseline before slice 2')}, {star_type or 'default stars'})")
    directory = os.environ.get("PLANETGEN_REPORT_DIR", ".")
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, os.environ.get("PLANETGEN_REPORT_NAME", REPORT_NAME)), "w", encoding="utf-8") as handle:
        handle.write(text)
    print(text)


if __name__ == "__main__":
    main()
