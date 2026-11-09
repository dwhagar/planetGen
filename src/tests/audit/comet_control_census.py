# tests/audit/comet_control_census.py

"""
The comet control census (GEN.1457, system-interactions.md rule 5, Comet orbits): every seed is generated as a whole system twice, *before* (the controlling giant is held at None, so the comets are drawn as they were) and *after*
(the comets of a star with a giant are drawn about its most massive one). The comets of a star are the first draws after its planets, so the two runs agree on everything up to them.

The rule: a star with no giant never changes; the families other than the Jupiter family never change in their elements; every Jupiter-family comet of a star with a giant has a Tisserand parameter about it in the family's range.
"""

import collections
import multiprocessing
import os

SEEDS = 480
REPORT_NAME = "comet-control-census-latest.md"


def _comets(system):
    return [(star_name, comet) for star_name, comets in (("primary", system.comets), ("secondary", system.secondary_comets)) for comet in comets]


def _generate(seed, held):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import comet_control
    from planetgen.util import draw

    real = comet_control.controlling_giant_au
    if held:
        comet_control.controlling_giant_au = lambda planets: None
    try:
        draw.set_run_seed(seed)
        return StarSystem(system_config=SystemConfig())
    finally:
        comet_control.controlling_giant_au = real


def system_counts(seed):
    from planetgen.physics import comet_control

    counts = collections.Counter()
    before = _generate(seed, True)
    after = _generate(seed, False)
    giant = comet_control.controlling_giant_au(after.planets)
    counts["systems"] += 1
    counts["systems with a giant"] += giant is not None
    if len(before.comets) != len(after.comets):
        counts["comet count differs"] += 1
        return counts
    moved = False    # a changed comet shifts the draws of the comets after it on the same star, so those are not compared
    for old, new in zip(before.comets, after.comets):
        counts["primary comets"] += 1
        key = (old.period_class, old.perihelion_distance_au, old.eccentricity, old.inclination_deg)
        same = key == (new.period_class, new.perihelion_distance_au, new.eccentricity, new.inclination_deg)
        if giant is None:
            counts["comets of stars with no giant"] += 1
            counts["changed with no giant"] += not same
            continue
        counts["comets of stars with a giant"] += 1
        if old.period_class == "jupiter_family":
            counts["jupiter-family before"] += 1
            counts["jupiter-family before: perihelion sum"] += old.perihelion_distance_au
            counts["jupiter-family before: in range"] += comet_control.controls(giant, old.perihelion_distance_au, old.eccentricity, old.inclination_deg)
        elif not same and not moved:
            counts["other family changed"] += 1
        if new.period_class == "jupiter_family":
            counts["jupiter-family after"] += 1
            counts["jupiter-family after: perihelion sum"] += new.perihelion_distance_au
            counts["jupiter-family after: out of range"] += not comet_control.controls(giant, new.perihelion_distance_au, new.eccentricity, new.inclination_deg)
        if old.period_class == "jupiter_family" and new.period_class != "jupiter_family":
            counts["jupiter-family reclassified"] += 1
        moved = moved or not same
    return counts


def run(seeds=SEEDS, workers=None):
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(system_counts, range(1, seeds + 1), chunksize=1)
    total = collections.Counter()
    for part in parts:
        total.update(part)
    return total


ROWS = ("systems", "systems with a giant", "comet count differs", "primary comets", "comets of stars with no giant", "changed with no giant", "comets of stars with a giant", "jupiter-family before", "jupiter-family before: in range",
        "jupiter-family after", "jupiter-family after: out of range", "jupiter-family reclassified", "other family changed")


def violations(counts):
    """The census findings that break the rule: any change with no giant, any other family changed, any Jupiter-family comet out of range, any comet count that moved."""
    return [name for name in ("changed with no giant", "other family changed", "jupiter-family after: out of range", "comet count differs") if counts[name]]


def format_report(counts, head):
    lines = [f"# Comet control census (GEN.1457), head {head}", "", f"{counts['systems']} systems, comets drawn before and after the controlling giant.", "", "| row | count |", "|---|---|"]
    lines += [f"| {row} | {counts[row]} |" for row in ROWS]
    for label in ("before", "after"):
        n = counts[f"jupiter-family {label}"]
        if n:
            lines.append(f"| jupiter-family {label}: mean perihelion (AU) | {counts[f'jupiter-family {label}: perihelion sum'] / n:.3f} |")
    broken = violations(counts)
    lines += ["", "Rule violations: " + (", ".join(broken) if broken else "none")]
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
