# tests/audit/scatterer_census.py

"""
The scatterer census (GEN.1456, system-interactions.md rule 5): every seed is generated as a whole system and the losses of each star's planet architecture are read twice, *before* (the scatterer of a thrown or
swept planet is the neighbour the clearing recorded, `planet_architecture._reaches` held False) and *after* (the planet of the highest Safronov number among all whose reach holds the body).

The rule: a star with no giant and a system of one planet never change; a body no other planet reaches never changes.
"""

import collections
import multiprocessing
import os
import pathlib

SEEDS = 480
REPORT_NAME = "scatterer-census-latest.md"


def _read(architecture):
    from planetgen.physics import planet_architecture

    return [(loss.kind, loss.perturber_mass_kg, round(loss.ejected, 12)) for loss in planet_architecture.losses(architecture) if loss.kind in ("thrown", "swept")]


def system_counts(seed):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import planet_architecture
    from planetgen.util import draw

    draw.set_run_seed(seed)
    system = StarSystem(system_config=SystemConfig())
    counts = collections.Counter()
    real = planet_architecture._reaches
    for _star, architecture in getattr(system, "architectures", ()):
        group = "companion" if len(system.stars) > 1 else "single"
        giants = any(planet.runaway for planet in architecture.planets)
        after = _read(architecture)
        planet_architecture._reaches = lambda *args: False
        try:
            before = _read(architecture)
        finally:
            planet_architecture._reaches = real
        counts[(group, "architectures")] += 1
        counts[(group, "with giants" if giants else "without giants")] += 1
        counts[(group, "bodies scattered")] += len(after)
        for old, new in zip(before, after):
            if old != new:
                counts[(group, "scatterer changed")] += 1
                counts[(group, "ejected share before")] += old[2]
                counts[(group, "ejected share after")] += new[2]
                if not giants:
                    counts[(group, "changed with no giant")] += 1
        if len(architecture.planets) <= 1 and before != after:
            counts[(group, "changed with one planet")] += 1
    return counts


def run(seeds=SEEDS, workers=None):
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(system_counts, range(1, seeds + 1), chunksize=1)
    total = collections.Counter()
    for part in parts:
        total.update(part)
    return total


GROUPS = ("single", "companion")
ROWS = ("architectures", "with giants", "without giants", "bodies scattered", "scatterer changed", "ejected share before", "ejected share after", "changed with no giant", "changed with one planet")


def violations(counts):
    """A star with one planet, or with none that is a giant, never changes."""
    return [f"{counts[(group, key)]} {key} among the {group} stars" for group in GROUPS for key in ("changed with one planet",) if counts[(group, key)]]


def format_report(counts, header=""):
    lines = ["# Scatterer census (GEN.1456)", "", header, "", "| measure | " + " | ".join(GROUPS) + " |", "|---|" + "---|" * len(GROUPS)]
    lines += [f"| {key} | " + " | ".join(f"{counts[(group, key)]:.3g}" if isinstance(counts[(group, key)], float) else str(counts[(group, key)]) for group in GROUPS) + " |" for key in ROWS]
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
    print(write_report(format_report(run(), f"Seeds 1 to {SEEDS}, the losses of every star's planet architecture read with and without the planets that reach the body.")))
