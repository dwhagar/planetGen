# tests/audit/planet_forced_census.py

"""
The forced planets census (GEN.1499 slice 1): every seed is generated as a pair of stars with planets, twice, with the companion's secular terms (`belt_forcing.stellar_terms`) wired into the
planets' exchange and with them switched off. It reports, for the stars with a companion, how many planets the terms move, the mean eccentricity and inclination before and after, and how
many exchanges the terms turn into "outside the theory" (a forced eccentricity or tilt past the limits keeps the planets as they were). A lone star's planets are not read: they give today's
numbers by construction (`tests/test_secular_chain.py`, the first test of GEN.1499).
"""

import multiprocessing
import os
import statistics

SEEDS = 200
REPORT_NAME = "planet-forced-census-latest.md"


def _planets(system):
    from planetgen.physics import radiation_field
    out = {}
    for star, planets in ((host.star, host.orbiters) for host in system.hosts):
        if not radiation_field.companions_of(star):
            continue
        out[star.name] = (system.secular_outcomes.get(star.name),
                          [(p.eccentricity, p.orbital_inclination_deg or 0.0) for p in planets if p.body_type != 'a'])
    return out


def system_rows(seed):
    """`[(seed, planets, moved, mean_e_before, mean_e_after, mean_i_before, mean_i_after, outcome_before, outcome_after)]`, one for each star with a companion and planets."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import belt_forcing
    from planetgen.util import draw

    def build(wired):
        cfg = SystemConfig()
        cfg.STARS, cfg.BINARY_SYSTEM, cfg.PLANETS, cfg.MAX_PLANETS = 2, True, True, True
        draw.set_run_seed(seed)
        original = belt_forcing.stellar_terms
        if not wired:
            belt_forcing.stellar_terms = lambda *args, **kwargs: None
        try:
            return _planets(StarSystem(system_config=cfg))
        finally:
            belt_forcing.stellar_terms = original

    with_terms, without = build(True), build(False)
    rows = []
    for name, (outcome, planets) in with_terms.items():
        before_outcome, before = without.get(name, (None, []))
        if not planets or len(before) != len(planets):
            continue
        moved = sum(1 for a, b in zip(before, planets) if abs(a[0] - b[0]) > 1e-9 or abs(a[1] - b[1]) > 1e-6)
        rows.append((seed, len(planets), moved, statistics.fmean(a[0] for a in before), statistics.fmean(b[0] for b in planets),
                     statistics.fmean(a[1] for a in before), statistics.fmean(b[1] for b in planets), before_outcome, outcome))
    return rows


def run(seeds=SEEDS, workers=None):
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(system_rows, range(1, seeds + 1), chunksize=1)
    return [row for part in parts for row in part]


def summarize(rows):
    """The numbers of the report from the rows `system_rows` gives."""
    moved = [row for row in rows if row[2] > 0]
    return {"stars": len(rows), "planets": sum(row[1] for row in rows), "moved_stars": len(moved), "moved_planets": sum(row[2] for row in rows),
            "e_before": statistics.fmean(row[3] for row in moved) if moved else float("nan"),
            "e_after": statistics.fmean(row[4] for row in moved) if moved else float("nan"),
            "i_before": statistics.fmean(row[5] for row in moved) if moved else float("nan"),
            "i_after": statistics.fmean(row[6] for row in moved) if moved else float("nan"),
            "turned_outside": sum(1 for row in rows if row[8] == "outside the theory" and row[7] != "outside the theory")}


def format_report(rows, heading):
    s = summarize(rows)
    return "\n".join([
        f"# {heading}", "",
        f"Stars with a companion and planets in {SEEDS} generated pairs: {s['stars']} ({s['planets']} planets).", "",
        "| What | Value |", "|---|---|",
        f"| stars whose planets the companion terms move | {s['moved_stars']} |",
        f"| planets moved | {s['moved_planets']} of {s['planets']} |",
        f"| mean eccentricity of a moved star's planets: before and after | {s['e_before']:.3g} and {s['e_after']:.3g} |",
        f"| mean inclination of a moved star's planets, degrees: before and after | {s['i_before']:.3g} and {s['i_after']:.3g} |",
        f"| stars the terms turn to \"outside the theory\" | {s['turned_outside']} |", ""])


def main():
    text = format_report(run(), "The forced planets census (GEN.1499 slice 1)")
    directory = os.environ.get("PLANETGEN_REPORT_DIR", ".")
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, REPORT_NAME), "w", encoding="utf-8") as handle:
        handle.write(text)
    print(text)


if __name__ == "__main__":
    main()
