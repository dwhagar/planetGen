# tests/audit/late_accretion_census.py

"""
The late accretion census (GEN.1493 slice 1): every seed is generated as a whole system with the late accretion pass on (`generation/late_accretion.py`) and each rocky planet's
percent, the one the pass set, is set against the mean the law gives its mass (`impact_history.late_accreted_percent_mean`).

It reports, for rocky planets of 0.3 to 3 Earth masses of single stars, the share that took the floor (no belt reached them), the share held at the ceiling, the median and the spread of
(percent) / (the law's mean), what the pass cost with the gate (`planets.LATE_ACCRETION_GATE`, row LA3) against the pass for every planet, and how far the gated air differs from the full
replay (every seed is made twice). Rows LA1 (floor), LA2 (no scatter of its own) and LA3 (the gate) are read from it.
"""

import multiprocessing
import os
import statistics
import time

SEEDS = 240
REPORT_NAME = "late-accretion-census-latest.md"
EARTH_KG = 5.9722e24


def _generate(seed, gate):
    """`(system, seconds of the pass, planets made again, seconds of the whole system)` of the seed with the pass at `gate`."""
    from planetgen.generation import late_accretion
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import planets as planet_physics
    from planetgen.util import draw

    spent, made = [0.0], [0]
    real_apply, real_remake = late_accretion.apply, planet_physics.remake_air

    def timed(planets, gate_=gate):
        began = time.perf_counter()
        try:
            return real_apply(planets, gate)
        finally:
            spent[0] += time.perf_counter() - began

    def counted(planet, percent, gate_=gate):
        remade = real_remake(planet, percent, gate_)
        made[0] += bool(remade)
        return remade

    late_accretion.apply, planet_physics.remake_air = timed, counted
    try:
        draw.set_run_seed(seed)
        began = time.perf_counter()
        system = StarSystem(system_config=SystemConfig())
        total = time.perf_counter() - began
    finally:
        late_accretion.apply, planet_physics.remake_air = real_apply, real_remake
    return system, spent[0], made[0], total


def system_rows(seed):
    """
    `[(seed, mass_earth, percent, mean_percent, floored, capped, pressure_change, temperature_change)]` of the rocky planets in range of the system made with the gate (the changes against the
    same system made again for every planet), and the last row `("time", gated pass, full pass, whole gated system, planets made again, planets that apply)`.
    """
    from planetgen.generation import late_accretion
    from planetgen.physics import impact_history

    system, gated_s, made, total = _generate(seed, late_accretion.planet_physics.LATE_ACCRETION_GATE)
    full, full_s, applies, _total = _generate(seed, 0.0)
    rows = []
    if len(system.stars) == 1:
        for body, twin in zip(system.all_planets, full.all_planets):
            percent = getattr(body, "_late_accreted_percent", None)
            if body.body_type != "t" or body.is_moon or percent is None or not body.mass:
                continue
            mass_earth = body.mass / EARTH_KG
            if not (0.3 <= mass_earth <= 3.0):
                continue
            mean = impact_history.late_accreted_percent_mean(mass_earth)
            pressure = abs(body.atmospheric_pressure - twin.atmospheric_pressure) / twin.atmospheric_pressure if twin.atmospheric_pressure else 0.0
            temperature = abs(body.surface_temperature - twin.surface_temperature) / twin.surface_temperature
            rows.append((seed, mass_earth, percent, mean, percent <= late_accretion.FLOOR_SHARE * mean * (1.0 + 1e-9),
                         percent >= impact_history.LATE_ACCRETED_PERCENT_MAX * (1.0 - 1e-9), pressure, temperature))
    rows.append(("time", gated_s, full_s, total, made, applies))
    return rows


def run(seeds=SEEDS, workers=None):
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(system_rows, range(1, seeds + 1), chunksize=1)
    return [row for part in parts for row in part]


def summarize(rows):
    """The numbers of the report from the rows `system_rows` gives."""
    planets = [row for row in rows if row[0] != "time"]
    times = [row for row in rows if row[0] == "time"]
    ratios = sorted(row[2] / row[3] for row in planets)
    pressure = sorted(row[6] for row in planets)
    temperature = sorted(row[7] for row in planets)
    n = len(ratios)

    def pick(values, share):
        return values[min(n - 1, int(share * n))] if n else float("nan")

    return {"planets": n, "floor": sum(1 for row in planets if row[4]), "cap": sum(1 for row in planets if row[5]),
            "median": statistics.median(ratios) if n else float("nan"), "p10": pick(ratios, 0.10), "p90": pick(ratios, 0.90),
            "same": sum(1 for row in planets if row[6] == 0.0 and row[7] == 0.0),
            "over_one_percent": sum(1 for row in planets if row[6] > 0.01 or row[7] > 0.01),
            "dp_p90": pick(pressure, 0.90), "dp_max": pressure[-1] if n else float("nan"), "dt_p90": pick(temperature, 0.90), "dt_max": temperature[-1] if n else float("nan"),
            "pass_s": sum(row[1] for row in times), "full_s": sum(row[2] for row in times), "total_s": sum(row[3] for row in times),
            "made": sum(row[4] for row in times), "applies": sum(row[5] for row in times)}


def format_report(rows, heading):
    s = summarize(rows)
    n = max(s["planets"], 1)
    return "\n".join([
        f"# {heading}", "",
        f"Rocky planets of 0.3 to 3 Earth masses of single stars, {SEEDS} seeds: {s['planets']}.", "",
        "| What | Value |", "|---|---|",
        f"| took the floor (no belt reached them) | {s['floor']} ({100.0 * s['floor'] / n:.0f} percent) |",
        f"| held at the ceiling | {s['cap']} ({100.0 * s['cap'] / n:.0f} percent) |",
        f"| percent over the law's mean: median | {s['median']:.3g} |",
        f"| percent over the law's mean: 10th and 90th | {s['p10']:.3g} and {s['p90']:.3g} |",
        f"| airs made again, of the rocky planets with one to make | {s['made']} of {s['applies']} |",
        f"| identical to the full replay | {s['same']} of {s['planets']} |",
        f"| changed by more than 1 percent in pressure or temperature | {s['over_one_percent']} of {s['planets']} |",
        f"| pressure against the full replay: 90th percentile and largest change | {s['dp_p90']:.2g} and {s['dp_max']:.2g} |",
        f"| temperature against the full replay: 90th percentile and largest change | {s['dt_p90']:.2g} and {s['dt_max']:.2g} |",
        f"| the pass with the gate, seconds over all seeds | {s['pass_s']:.1f} ({100.0 * s['pass_s'] / max(s['total_s'], 1e-9):.1f} percent of the whole systems) |",
        f"| the pass for every planet (gate 0), seconds over all seeds | {s['full_s']:.1f} |",
        f"| whole systems with the gate, seconds over all seeds | {s['total_s']:.1f} |", ""])


def main():
    rows = run()
    text = format_report(rows, "The late accretion census (GEN.1493 slice 1)")
    directory = os.environ.get("PLANETGEN_REPORT_DIR", ".")
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, REPORT_NAME), "w", encoding="utf-8") as handle:
        handle.write(text)
    print(text)


if __name__ == "__main__":
    main()
