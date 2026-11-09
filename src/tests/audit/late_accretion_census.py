# tests/audit/late_accretion_census.py

"""
The late accretion census (GEN.1493 slice 1): every seed is generated as a whole system with the late accretion pass on (`generation/late_accretion.py`) and each rocky planet's
percent, the one the pass set, is set against the mean the law gives its mass (`impact_history.late_accreted_percent_mean`).

It reports, for rocky planets of 0.3 to 3 Earth masses of single stars, the share that took the floor (no belt reached them), the share held at the ceiling, the median and the spread of
(percent) / (the law's mean), and what the pass cost: the seconds it took against the seconds of the whole system. Rows LA1 (floor) and LA2 (no scatter of its own) are read from it.
"""

import multiprocessing
import os
import statistics
import time

SEEDS = 240
REPORT_NAME = "late-accretion-census-latest.md"
EARTH_KG = 5.9722e24


def system_rows(seed):
    """`[(seed, mass_earth, percent, mean_percent, floored, capped)]` of the rocky planets in range, and the seconds `(pass, whole system)` as the last row `("time", pass, total)`."""
    from planetgen.generation import late_accretion
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import impact_history
    from planetgen.util import draw

    spent = [0.0]
    real = late_accretion.apply

    def timed(planets):
        began = time.perf_counter()
        try:
            return real(planets)
        finally:
            spent[0] += time.perf_counter() - began

    late_accretion.apply = timed
    try:
        draw.set_run_seed(seed)
        began = time.perf_counter()
        system = StarSystem(system_config=SystemConfig())
        total = time.perf_counter() - began
    finally:
        late_accretion.apply = real
    rows = []
    if len(system.stars) == 1:
        for body in system.all_planets:
            percent = getattr(body, "_late_accreted_percent", None)
            if body.body_type != "t" or body.is_moon or percent is None or not body.mass:
                continue
            mass_earth = body.mass / EARTH_KG
            if not (0.3 <= mass_earth <= 3.0):
                continue
            mean = impact_history.late_accreted_percent_mean(mass_earth)
            rows.append((seed, mass_earth, percent, mean, percent <= late_accretion.FLOOR_SHARE * mean * (1.0 + 1e-9),
                         percent >= impact_history.LATE_ACCRETED_PERCENT_MAX * (1.0 - 1e-9)))
    rows.append(("time", spent[0], total))
    return rows


def run(seeds=SEEDS, workers=None):
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(system_rows, range(1, seeds + 1), chunksize=1)
    return [row for part in parts for row in part]


def summarize(rows):
    """`{"planets", "floor", "cap", "median", "p10", "p90", "pass_s", "total_s"}` from the rows `system_rows` gives."""
    planets = [row for row in rows if row[0] != "time"]
    times = [row for row in rows if row[0] == "time"]
    ratios = sorted(row[2] / row[3] for row in planets)
    n = len(ratios)

    def pick(share):
        return ratios[min(n - 1, int(share * n))] if n else float("nan")

    return {"planets": n, "floor": sum(1 for row in planets if row[4]), "cap": sum(1 for row in planets if row[5]),
            "median": statistics.median(ratios) if n else float("nan"), "p10": pick(0.10), "p90": pick(0.90),
            "pass_s": sum(row[1] for row in times), "total_s": sum(row[2] for row in times)}


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
        f"| the pass, seconds over all seeds | {s['pass_s']:.1f} |",
        f"| whole systems, seconds over all seeds | {s['total_s']:.1f} ({100.0 * s['pass_s'] / max(s['total_s'], 1e-9):.0f} percent is the pass) |", ""])


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
