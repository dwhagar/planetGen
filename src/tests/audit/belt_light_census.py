# tests/audit/belt_light_census.py

"""
The belt and comet light census (GEN.1444, and the rule of docs/design/generation/system-interactions.md that a change of the radiation field ships with a before and after).

Each seed is generated twice as a whole system: *before*, with the belts and the comets reading their host star's light alone (the ring's light ratio held at 1.0 and a comet's lit
perihelion held at its perihelion, which is what the generator did before GEN.1444), and *after*, as it is now. Both runs draw from the same seed, so a body that is the same in both has
the same draws; a system whose layout differs (a belt that was left out or put elsewhere because its size law moved the mass the belt keeps) is counted apart and not compared body by body.

For the systems whose belts and comets are the same in both runs, every belt is compared on its mean temperature (the temperature its bodies' makeup was drawn at), its size law, its mass and
its taxon makeup, and every comet on its chance of showing activity and on whether it does.

The rule the census checks (`violations`): a belt or a comet of a lone star never changes.

`python -m tests.audit.belt_light_census truncation` runs the same kind of before and after for GEN.1451 (`truncation_counts`): each seed with and without the companion's cut of the disk
(`disk_truncation.star_disk_limit_au`), the belts matched by their inner limit, since the cut moves the outer one and can end a belt.
"""

import collections
import multiprocessing
import os
import pathlib
from types import SimpleNamespace

SEEDS = 240
"""int: Seeds 1 to this are generated."""
REPORT_NAME = "belt-light-census-latest.md"
PERCENTS = (0.1, 0.5, 1, 5)
"""tuple: The temperature changes (percent) the report counts belts over."""


def _before(module_belt, module_comet):
    """Puts the generator back as it read the host alone: returns what to put back."""
    saved = (module_belt.ring_light_ratio, module_comet.radiation_field)
    module_belt.ring_light_ratio = lambda star, lower, upper: 1.0
    module_comet.radiation_field = SimpleNamespace(lit_distance_au=lambda star, distance: distance)
    return saved


def _undo(module_belt, module_comet, saved):
    module_belt.ring_light_ratio, module_comet.radiation_field = saved


def _bodies(system):
    belts = [body for body in system.all_planets if getattr(body, "body_type", None) == "a"]
    return belts, list(system.comets)


def _snapshot(seed):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    draw.set_run_seed(seed)
    system = StarSystem(system_config=SystemConfig())
    belts, comets = _bodies(system)
    return {
        "stars": len(system.stars),
        "bound": len(system.stars) > 1 or bool(system.dark_members) or system.star.a_crit_au is not None,     # a companion of any kind, a dark one too, limits the star's orbits
        "belts": [(belt.distance, belt.lower_limit, belt.upper_limit, belt.mass_kg, belt.largest_body_km, belt.makeup.detail["temperature_k"], belt.density,
                   tuple(belt.makeup.detail.get("taxon_shares", {}).items()) if isinstance(belt.makeup.detail.get("taxon_shares"), dict) else ()) for belt in belts],
        "comets": [(comet.perihelion_distance_au, comet.is_active) for comet in comets],
        "law": [_law(system, belt) for belt in belts],
    }


def _law(system, belt):
    from planetgen.generation import belt as belt_generation

    star = belt_generation._STARS.get(belt)
    return None if star is None else belt_generation.size_law(belt.distance, star.snow_line_au, belt_generation.ring_light_ratio(star, belt.lower_limit, belt.upper_limit))


def system_counts(seed):
    """Counts for one seed (see the module): `{key: number}` with the keys named in `format_report`."""
    from planetgen.generation import belt as module_belt
    from planetgen.generation import comet as module_comet

    after = _snapshot(seed)
    saved = _before(module_belt, module_comet)
    try:
        before = _snapshot(seed)
    finally:
        _undo(module_belt, module_comet, saved)
    group = "lone" if after["stars"] < 2 else "multi"
    counts = collections.Counter()
    counts[(group, "systems")] += 1
    if len(before["belts"]) != len(after["belts"]) or len(before["comets"]) != len(after["comets"]) or \
            [belt[:3] for belt in before["belts"]] != [belt[:3] for belt in after["belts"]]:
        counts[(group, "layout_changed")] += 1
        return counts
    for old, new, law in zip(before["belts"], after["belts"], after["law"]):
        counts[(group, "belts")] += 1
        warmer = new[5] / old[5] - 1.0
        if new[:3] != old[:3] or new[3:] != old[3:]:
            counts[(group, "belts_changed")] += 1
        if new[5] != old[5]:
            counts[(group, "belts_temperature_changed")] += 1
            for percent in PERCENTS:
                if abs(warmer) > percent / 100.0:
                    counts[(group, f"belts_temperature_over_{percent}")] += 1
        if new[4] != old[4] or new[3] != old[3]:
            counts[(group, "belts_size_or_mass_changed")] += 1
        if new[7] != old[7]:
            counts[(group, "belts_taxa_changed")] += 1
    for old, new in zip(before["comets"], after["comets"]):
        counts[(group, "comets")] += 1
        if new != old:
            counts[(group, "comets_changed")] += 1
    return counts


def _no_truncation(module_birth):
    saved = module_birth.disk_truncation
    module_birth.disk_truncation = SimpleNamespace(star_disk_limit_au=lambda star: float("inf"))
    return saved


def truncation_counts(seed):
    """GEN.1451, the companion's cut of a disk: counts for one seed, each seed generated with and without the cut. Belts are matched by their inner limit (the cut moves the outer one)."""
    import collections as _collections

    from planetgen.generation import belt_birth

    after = _snapshot(seed)
    saved = _no_truncation(belt_birth)
    try:
        before = _snapshot(seed)
    finally:
        belt_birth.disk_truncation = saved
    group = "multi" if after["bound"] else "lone"     # a star with a dark companion is not alone: its stability limit cuts its belts
    counts = _collections.Counter()
    counts[(group, "systems")] += 1
    old = {round(belt[1], 9): belt for belt in before["belts"]}
    new = {round(belt[1], 9): belt for belt in after["belts"]}
    counts[(group, "belts_before")] += len(old)
    counts[(group, "belts_after")] += len(new)
    counts[(group, "belts_lost")] += len(set(old) - set(new))
    counts[(group, "belts_new")] += len(set(new) - set(old))
    for key in set(old) & set(new):
        if new[key][2] < old[key][2]:
            counts[(group, "belts_cut")] += 1
        if new[key][3] != old[key][3]:
            counts[(group, "belts_mass_changed")] += 1
    return counts


TRUNCATION_ROWS = (("systems", "Systems"), ("belts_before", "Belts without the cut"), ("belts_after", "Belts with the cut"), ("belts_lost", "Belts that no longer exist"),
                   ("belts_new", "Belts that did not exist"), ("belts_cut", "Belts whose outer limit is cut"), ("belts_mass_changed", "Belts whose mass changed"))


def truncation_violations(counts):
    """A lone star's belts never change with the cut."""
    return [f"{counts[('lone', key)]} {key.replace('_', ' ')} among the lone stars" for key in ("belts_lost", "belts_new", "belts_cut", "belts_mass_changed") if counts[("lone", key)]]


def format_truncation_report(counts, header=""):
    lines = ["# Companion disk-cut census (GEN.1451)", "", header, "", "| measure | lone stars | multi-star systems |", "|---|---|---|"]
    lines += [f"| {label} | {counts[('lone', key)]} | {counts[('multi', key)]} |" for key, label in TRUNCATION_ROWS]
    found = truncation_violations(counts)
    lines += ["", "Rule violations: " + ("none" if not found else ""), *[f"- {item}" for item in found], ""]
    return "\n".join(lines)


def run(seeds=SEEDS, workers=None, job=None):
    """The summed counts of seeds 1 to `seeds`."""
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(job or system_counts, range(1, seeds + 1), chunksize=1)
    total = collections.Counter()
    for part in parts:
        total.update(part)
    return total


def violations(counts):
    """What breaks the rule over `counts` (sentences; empty when it holds)."""
    found = []
    for key in ("layout_changed", "belts_changed", "comets_changed"):
        if counts[("lone", key)]:
            found.append(f"{counts[('lone', key)]} {key.replace('_', ' ')} among the lone stars")
    return found


ROWS = (("systems", "Systems"), ("layout_changed", "Systems whose belts or comets differ in layout (left out below)"), ("belts", "Belts compared"),
        ("belts_changed", "Belts with any stored value changed"), ("belts_temperature_changed", "Belts whose temperature changed"),
        *((f"belts_temperature_over_{percent}", f"... by over {percent} percent") for percent in PERCENTS),
        ("belts_size_or_mass_changed", "Belts whose mass or largest body changed (the size law moved)"), ("belts_taxa_changed", "Belts whose taxon makeup changed"),
        ("comets", "Comets compared"), ("comets_changed", "Comets whose activity changed"))


def format_report(counts, header=""):
    """The census as Markdown."""
    lines = ["# Belt and comet light census (GEN.1444)", "", header, "", "| measure | lone stars | multi-star systems |", "|---|---|---|"]
    lines += [f"| {label} | {counts[('lone', key)]} | {counts[('multi', key)]} |" for key, label in ROWS]
    found = violations(counts)
    lines += ["", "Rule violations: " + ("none" if not found else ""), *[f"- {item}" for item in found], ""]
    return "\n".join(lines)


def write_report(text, directory=None, name=REPORT_NAME):
    """Writes the report to `directory`, `$PLANETGEN_REPORT_DIR` or `.test-reports/`."""
    directory = pathlib.Path(directory or os.environ.get("PLANETGEN_REPORT_DIR") or ".test-reports")
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(text, encoding="utf-8")
    return path


if __name__ == "__main__":
    import sys

    if "truncation" in sys.argv[1:]:
        print(write_report(format_truncation_report(run(job=truncation_counts), f"Seeds 1 to {SEEDS}, each generated with and without the companion's cut of the disk."), name="belt-truncation-census-latest.md"))
    else:
        print(write_report(format_report(run(), f"Seeds 1 to {SEEDS}, each generated before and after GEN.1444.")))
