# tests/audit/belt_forced_members_census.py

"""
The forced members census (GEN.1492 slice 1): every seed is generated as a pair of stars with a belt, and each belt of a star with a companion has its stored members drawn twice, with the
companion's secular forcing wired in (`AsteroidBelt._force_by_companions`) and with it switched off. It reports how many belts and members the forcing moves, the mean eccentricity before and
after, and the members it leaves alone because they stand outside the theory. A belt of a lone star is not read: it gives today's numbers by construction (the first test of
`tests/test_belt_forced_members.py`).
"""

import multiprocessing
import os
import statistics

SEEDS = 300
REPORT_NAME = "belt-forced-members-census-latest.md"


def system_rows(seed):
    """`[(seed, members, moved, mean_e_before, mean_e_after, over_half)]` of the belts of the stars with companions of one pair of stars."""
    from planetgen.generation import belt as belt_module
    from planetgen.generation.belt import AsteroidBelt
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import radiation_field
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
        if star is None or not radiation_field.companions_of(star):
            continue
        forced = [dict(r) for r in body.members]
        body._members = None
        original = AsteroidBelt._force_by_companions
        AsteroidBelt._force_by_companions = lambda self, star, rows: rows
        try:
            bare = [dict(r) for r in body.members]
        finally:
            AsteroidBelt._force_by_companions = original
            body._members = None
        if not bare:
            continue
        moved = sum(1 for a, b in zip(bare, forced) if a["eccentricity"] != b["eccentricity"])
        before = statistics.fmean(r["eccentricity"] for r in bare)
        after = statistics.fmean(r["eccentricity"] for r in forced)
        rows.append((seed, len(bare), moved, before, after, after > 1.5 * before))
    return rows


def run(seeds=SEEDS, workers=None):
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(system_rows, range(1, seeds + 1), chunksize=1)
    return [row for part in parts for row in part]


def summarize(rows):
    """The numbers of the report from the rows `system_rows` gives."""
    reached = [row for row in rows if row[2] > 0]
    return {"belts": len(rows), "moved_belts": len(reached), "members": sum(row[1] for row in rows), "moved_members": sum(row[2] for row in rows),
            "before": statistics.fmean(row[3] for row in reached) if reached else float("nan"),
            "after": statistics.fmean(row[4] for row in reached) if reached else float("nan"),
            "over_half": sum(1 for row in reached if row[5])}


def format_report(rows, heading):
    s = summarize(rows)
    n = max(s["moved_belts"], 1)
    return "\n".join([
        f"# {heading}", "",
        f"Belts of stars with a companion in {SEEDS} generated pairs: {s['belts']}.", "",
        "| What | Value |", "|---|---|",
        f"| belts whose members the forcing moves | {s['moved_belts']} |",
        f"| members of those belts moved | {s['moved_members']} of {s['members']} (in all {s['belts']} belts) |",
        f"| mean eccentricity of a moved belt: before and after | {s['before']:.3g} and {s['after']:.3g} |",
        f"| belts whose mean eccentricity rises by half or more | {s['over_half']} ({100.0 * s['over_half'] / n:.0f} percent of the moved) |", ""])


def main():
    rows = run()
    text = format_report(rows, "The forced members census (GEN.1492 slice 1)")
    directory = os.environ.get("PLANETGEN_REPORT_DIR", ".")
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, REPORT_NAME), "w", encoding="utf-8") as handle:
        handle.write(text)
    print(text)


if __name__ == "__main__":
    main()
