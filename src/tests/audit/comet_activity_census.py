# tests/audit/comet_activity_census.py

"""
The comet activity census (GEN.1491, system-interactions.md rule 5, Comet activity): every seed is generated as a whole system and each comet's chance of activity is read twice, *before* (the Sun's 3 AU read against the
perihelion as lit by every star, GEN.1444) and *after* (that distance over the square root of the host's luminosity in suns).

The rule: a comet of a host of exactly the Sun's luminosity, or of no star, never changes; a brighter host gives a higher chance and a fainter one a lower chance at the same perihelion.
"""

import collections
import multiprocessing
import os

SEEDS = 480
REPORT_NAME = "comet-activity-census-latest.md"
CLASSES = (("under 0.1 Lsun", 0.0, 0.1), ("0.1 to 3 Lsun", 0.1, 3.0), ("over 3 Lsun", 3.0, float("inf")))


def _class_of(solar_luminosities):
    for name, low, high in CLASSES:
        if low <= solar_luminosities < high:
            return name
    return CLASSES[-1][0]


def system_counts(seed):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.generation.comet import _activity_chance
    from planetgen.physics import comet_activity, constants, radiation_field
    from planetgen.util import draw

    draw.set_run_seed(seed)
    system = StarSystem(system_config=SystemConfig())
    counts = collections.Counter()
    for star, comets in ((system.star, system.comets), (getattr(system, "secondary_star", None), system.secondary_comets)):
        if star is None:
            continue
        luminosity = getattr(star, "luminosity", None)
        for comet in comets:
            lit = radiation_field.lit_distance_au(star, comet.perihelion_distance_au)
            before = _activity_chance(lit)
            after = _activity_chance(comet_activity.solar_equivalent_au(lit, luminosity))
            group = _class_of((luminosity or constants.SOLAR_LUMINOSITY) / constants.SOLAR_LUMINOSITY)
            counts[(group, "comets")] += 1
            counts[(group, "chance before")] += before
            counts[(group, "chance after")] += after
            counts[(group, "active")] += bool(comet.is_active)
            if abs(before - after) > 1e-12:
                counts[(group, "chance changed")] += 1
            if luminosity == constants.SOLAR_LUMINOSITY and before != after:
                counts[(group, "changed at exactly the Sun's luminosity")] += 1
    return counts


def run(seeds=SEEDS, workers=None):
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(system_counts, range(1, seeds + 1), chunksize=1)
    total = collections.Counter()
    for part in parts:
        total.update(part)
    return total


ROWS = ("comets", "chance changed", "active", "chance before", "chance after", "changed at exactly the Sun's luminosity")


def violations(counts):
    """The findings that break the rule: any change at exactly the Sun's luminosity, or a brighter-host class whose mean chance fell or a fainter one whose rose."""
    broken = [group for group, _low, _high in CLASSES if counts[(group, "changed at exactly the Sun's luminosity")]]
    bright = CLASSES[2][0]
    faint = CLASSES[0][0]
    if counts[(bright, "chance after")] < counts[(bright, "chance before")]:
        broken.append(bright + " fell")
    if counts[(faint, "chance after")] > counts[(faint, "chance before")]:
        broken.append(faint + " rose")
    return broken


def format_report(counts, head):
    lines = [f"# Comet activity census (GEN.1491), head {head}", "", "Mean chance of activity per comet and the comets drawn active, by host luminosity.", "",
             "| row | " + " | ".join(name for name, _low, _high in CLASSES) + " |", "|---|" + "---|" * len(CLASSES)]
    for row in ROWS:
        cells = []
        for name, _low, _high in CLASSES:
            value = counts[(name, row)]
            if row.startswith("chance ") and row != "chance changed" and counts[(name, "comets")]:
                value = f"{value / counts[(name, 'comets')]:.3f}"
            cells.append(str(value))
        lines.append(f"| {row} | " + " | ".join(cells) + " |")
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
