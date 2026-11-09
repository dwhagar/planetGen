# tests/audit/trojan_census.py

"""
The trojan swarm census (GEN.1455, system-interactions.md rule 5): every seed is generated as a whole system and every outer-disk swarm the generator draws is recorded with the share of it that the neighbours leave
(`trojans.surviving_share`). The survival does not change a draw, so the swarm of before is the swarm now over the survival and one run per seed serves.

The rule: a belt of a system with one giant and no companion star never changes.
"""

import collections
import multiprocessing
import os
import pathlib

SEEDS = 960
REPORT_NAME = "trojan-census-latest.md"


def _swarms(seed):
    """`{"giants", "bound", "calls"}` of one seed: every outer-disk swarm the generator drew, as `(survival, neighbours, trojans share, share of before)`."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import trojans
    from planetgen.util import draw

    calls = []
    real = trojans.surviving_share

    def recorded(*args, **kwargs):
        survival = real(*args, **kwargs)
        calls.append((survival, len(args[3])))
        return survival

    from planetgen.physics import scatter

    real_shares = scatter.scatter_shares
    before = []

    def recorded_shares(*args, **kwargs):
        before.append(real_shares(*args, **{**kwargs, "trojan_survival": 1.0})["trojans"])     # the swarm of before: the same call with every neighbour left out
        return real_shares(*args, **kwargs)

    trojans.surviving_share = recorded
    scatter.scatter_shares = recorded_shares
    try:
        draw.set_run_seed(seed)
        system = StarSystem(system_config=SystemConfig())
    finally:
        trojans.surviving_share = real
        scatter.scatter_shares = real_shares
    giants = [b for b in system.all_planets if b.body_type == 'g']
    shares = [b.scatter_fates["trojans"] for b in system.all_planets if b.body_type == 'a' and getattr(b, "scatter_fates", None)]
    return {"giants": len(giants), "bound": len(system.stars) > 1 or bool(system.dark_members) or system.star.a_crit_au is not None,
            "calls": [(survival, neighbours, share, old) for (survival, neighbours), share, old in zip(calls, shares, before)]}


def system_counts(seed):
    """Counts for one seed: each swarm's survival is the share it keeps against the swarm of before (survival 1.0)."""
    seen = _swarms(seed)
    group = ("companion" if seen["bound"] else "single") + ("-multi-giant" if seen["giants"] > 1 else "")
    counts = collections.Counter()
    counts[(group, "systems")] += 1
    for survival, _neighbours, share, old in seen["calls"]:
        counts[(group, "swarms")] += 1
        counts[(group, "swarm_share_after")] += share
        counts[(group, "swarm_share_before")] += old
        if survival < 1.0:
            counts[(group, "swarms_changed")] += 1
            if survival < 0.5:
                counts[(group, "swarms_over_half_lost")] += 1
    return counts


def run(seeds=SEEDS, workers=None):
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(system_counts, range(1, seeds + 1), chunksize=1)
    total = collections.Counter()
    for part in parts:
        total.update(part)
    return total


GROUPS = ("single", "single-multi-giant", "companion", "companion-multi-giant")
ROWS = (("systems", "Systems"), ("swarms", "Outer-disk swarms drawn"), ("swarms_changed", "Swarms that changed"), ("swarms_over_half_lost", "... that lost over half"),
        ("swarm_share_before", "Sum of the swarm shares of before"), ("swarm_share_after", "Sum of the swarm shares now"))


def violations(counts):
    """A system with one giant and no companion never changes."""
    return [f"{counts[('single', 'swarms_changed')]} swarms changed among single-giant systems of lone stars"] if counts[("single", "swarms_changed")] else []


def format_report(counts, header=""):
    lines = ["# Trojan swarm census (GEN.1455)", "", header, "", "| measure | " + " | ".join(GROUPS) + " |", "|---|" + "---|" * len(GROUPS)]
    lines += [f"| {label} | " + " | ".join(str(counts[(group, key)]) for group in GROUPS) + " |" for key, label in ROWS]
    found = violations(counts)
    lines += ["", "Rule violations: " + ("none" if not found else ""), *[f"- {item}" for item in found], ""]
    return "\n".join(lines)


def write_report(text, directory=None, name=REPORT_NAME):
    directory = pathlib.Path(directory or os.environ.get("PLANETGEN_REPORT_DIR") or ".test-reports")
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(text, encoding="utf-8")
    return path


if __name__ == "__main__":
    print(write_report(format_report(run(), f"Seeds 1 to {SEEDS}, every outer-disk swarm with the share the neighbours leave.")))
