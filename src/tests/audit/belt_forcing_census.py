# tests/audit/belt_forcing_census.py

"""
The belt forcing census (GEN.1492 slice 0): every seed is generated as a pair of stars with a belt, and each belt of a star with a companion has the forced eccentricity its companions give it
(`physics/belt_forcing.py`) set against the eccentricity its members are drawn with (`belt_roster.rms_eccentricity`), with and without the star's planets.

It reports the belts a companion reaches (the middle of the belt inside `ALPHA_LIMIT` of the companion's axis), the median and spread of (forced eccentricity) / (the drawn rms), the share
whose forcing is larger than the drawn rms, and how often the star's planets dilute the forcing by half or more. The census reads, it does not wire: nothing stored or drawn changes.
"""

import math
import multiprocessing
import os
import statistics

SEEDS = 300
REPORT_NAME = "belt-forcing-census-latest.md"


def system_rows(seed):
    """`[(seed, alpha, rms, bare, with_planets, held)]` of the belts of the stars with companions of one pair of stars."""
    from planetgen.generation import belt as belt_module
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import belt_forcing, belt_roster, radiation_field
    from planetgen.util import draw

    cfg = SystemConfig()
    cfg.STARS, cfg.BINARY_SYSTEM, cfg.ASTEROID_BELT = 2, True, True
    draw.set_run_seed(seed)
    system = StarSystem(system_config=cfg)
    rows = []
    for body in system.all_planets:
        if body.body_type != "a":
            continue
        star = belt_module._STARS.get(body)
        companions = radiation_field.companions_of(star) if star is not None else ()
        if not companions or not getattr(star, "mass", None):
            continue
        middle = math.sqrt(body.lower_limit * body.upper_limit)
        sources = [(c.star.mass, c.axis_au, c.eccentricity, belt_forcing.periapsis_rad(star.mass, c.member_no, c.axis_au)) for c in companions if getattr(c.star, "mass", None)]
        siblings = [p for p in system.all_planets if p.body_type in ("t", "g") and not p.is_moon and p.mass and getattr(p, "star", None) is star]
        rate = belt_forcing.planet_rate(middle, [p.distance for p in siblings], [p.mass for p in siblings], star.mass)
        bare = abs(belt_forcing.forced_vector(middle, star.mass, sources, belt_forcing.planet_rate(middle, [], [], star.mass)))
        planets = abs(belt_forcing.forced_vector(middle, star.mass, sources, rate))
        alpha = middle / min(source[1] for source in sources)
        rows.append((seed, alpha, belt_roster.rms_eccentricity(body.lower_limit, body.upper_limit), bare, planets, alpha > belt_forcing.ALPHA_LIMIT))
    return rows


def run(seeds=SEEDS, workers=None):
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(system_rows, range(1, seeds + 1), chunksize=1)
    return [row for part in parts for row in part]


def summarize(rows):
    """The numbers of the report from the rows `system_rows` gives."""
    reached = [row for row in rows if not row[5] and row[3] > 0.0]
    ratios = sorted(row[3] / row[2] for row in reached)
    n = len(ratios)
    diluted = [row for row in reached if row[4] <= 0.5 * row[3]]
    return {"belts": len(rows), "held": sum(1 for row in rows if row[5]), "reached": n,
            "median": statistics.median(ratios) if n else float("nan"), "p10": ratios[int(0.1 * n)] if n else float("nan"), "p90": ratios[min(n - 1, int(0.9 * n))] if n else float("nan"),
            "over_rms": sum(1 for ratio in ratios if ratio > 1.0), "diluted": len(diluted),
            "median_alpha": statistics.median(row[1] for row in reached) if n else float("nan")}


def format_report(rows, heading):
    s = summarize(rows)
    n = max(s["reached"], 1)
    return "\n".join([
        f"# {heading}", "",
        f"Belts of stars with a companion in {SEEDS} generated pairs: {s['belts']}.", "",
        "| What | Value |", "|---|---|",
        f"| belts outside the theory (middle beyond `ALPHA_LIMIT` of the companion's axis) | {s['held']} |",
        f"| belts the companion forces | {s['reached']} (median alpha {s['median_alpha']:.3g}) |",
        f"| forced eccentricity over the drawn rms: median | {s['median']:.3g} |",
        f"| forced eccentricity over the drawn rms: 10th and 90th | {s['p10']:.3g} and {s['p90']:.3g} |",
        f"| forcing larger than the drawn rms | {s['over_rms']} ({100.0 * s['over_rms'] / n:.0f} percent) |",
        f"| forcing cut by half or more by the star's planets | {s['diluted']} ({100.0 * s['diluted'] / n:.0f} percent) |", ""])


def main():
    rows = run()
    text = format_report(rows, "The belt forcing census (GEN.1492 slice 0)")
    directory = os.environ.get("PLANETGEN_REPORT_DIR", ".")
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, REPORT_NAME), "w", encoding="utf-8") as handle:
        handle.write(text)
    print(text)


if __name__ == "__main__":
    main()
