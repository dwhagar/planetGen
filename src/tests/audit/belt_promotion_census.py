# tests/audit/belt_promotion_census.py

"""
The census of belt bodies made planets (GEN.1511, the report asked for by GEN.1509; the rule is GEN.1510, `generation/belt_planets.py`).

Seeds 1 to `SEEDS` are generated as whole systems. For every belt that has at least one promoted planet the census counts the planets by class and by the spectral class of the star
the belt goes round, the planets per belt, the generation time, and checks the rules of the promotion:

* no belt has more than `tuning.BELT_PLANETS_PER_BELT` promoted planets;
* every promoted planet has a class, no envelope, and a roundness ratio of 1 or more;
* a promoted planet's rank is no longer a member of its belt;
* the planets of a star stand in order of distance.

`python -m tests.audit.belt_promotion_census` writes `belt-promotion-census-latest.md` (to `$PLANETGEN_REPORT_DIR` or `.test-reports/`).
"""

import collections
import multiprocessing
import os
import pathlib
import time

SEEDS = 120
"""int: Seeds 1 to this are generated."""
REPORT_NAME = "belt-promotion-census-latest.md"


def _star_class(star):
    """The spectral class letter of a star ("M", "G", "WD" for a white dwarf, "other")."""
    text = str(getattr(star, "type", "") or "")
    if "white dwarf" in text.lower() or text.upper().startswith("D"):
        return "WD"
    return text[:1].upper() if text[:1].upper() in "OBAFGKM" else "other"


def system_counts(seed):
    """Counts for one seed (see the module): `{key: number}`; the keys are tuples `(kind, what)`."""
    from planetgen import tuning
    from planetgen.generation import irregular_moons, validation
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import body_strength
    from planetgen.util import draw

    draw.set_run_seed(seed)
    started = time.monotonic()
    system = StarSystem(system_config=SystemConfig())
    counts = collections.Counter()
    counts[("total", "systems")] += 1
    counts[("total", "seconds_x1000")] += int((time.monotonic() - started) * 1000)
    for star, bodies in validation.star_lists(system):
        kind = _star_class(star)
        belts = irregular_moons.belts_of(bodies)
        promoted = [body for body in bodies if getattr(body, "belt_rank", None) is not None]
        counts[("belts", kind)] += len(belts)
        distances = [body.distance for body in bodies]
        if distances != sorted(distances):
            counts[("violation", "planets out of order of distance")] += 1
        for order, belt in belts.items():
            mine = [body for body in promoted if body.belt_order == order]
            if mine:
                counts[("belts_with_planets", kind)] += 1
            if len(mine) > tuning.BELT_PLANETS_PER_BELT:
                counts[("violation", "a belt with more promoted planets than the limit")] += 1
            held = {member["rank"] for member in belt.members}
            if any(body.belt_rank in held for body in mine):
                counts[("violation", "a promoted rank is still a member")] += 1
        for body in promoted:
            counts[("planets", kind)] += 1
            counts[("class", body.planet_class or "none")] += 1
            counts[("class_by_star", f"{kind}/{body.planet_class or 'none'}")] += 1
            if body.planet_class is None:
                counts[("violation", "a promoted planet with no class")] += 1
            if body.envelope_mass_fraction:
                counts[("violation", "a promoted planet with an envelope")] += 1
            if body_strength.ratio_for(body.mass, min(body.ice_mass_fraction + (body.water_mass_fraction or 0.0), 1.0)) < 1.0:       # the member's status counts its water ice with the other ice (`belt_page.body_entry`)
                counts[("violation", "a promoted planet under the round mass")] += 1
    return counts


def run(seeds=SEEDS, workers=None, job=None):
    """The summed counts of seeds 1 to `seeds`."""
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, seeds)) as pool:
        parts = pool.map(job or system_counts, range(1, seeds + 1), chunksize=1)
    total = collections.Counter()
    for part in parts:
        total.update(part)
    return total


def violations(counts):
    """What breaks the rules over `counts` (sentences; empty when they hold)."""
    return [f"{number} x {what}" for (kind, what), number in sorted(counts.items()) if kind == "violation" and number]


def format_report(counts, header=""):
    """The census as Markdown: planets by class, by star class, per belt, and the rule violations."""
    systems = counts[("total", "systems")]
    star_classes = sorted({what for (kind, what) in counts if kind == "belts"})
    classes = sorted({what for (kind, what) in counts if kind == "class"})
    lines = ["# Belt bodies made planets (GEN.1510 census, GEN.1511)", "", header, "",
             f"Systems: {systems}. Generation time per system: {counts[('total', 'seconds_x1000')] / 1000.0 / max(systems, 1):.2f} s.", "",
             "| star class | belts | belts with planets | planets made | planets per belt with any |", "|---|---|---|---|---|"]
    for star in star_classes:
        have, made = counts[("belts_with_planets", star)], counts[("planets", star)]
        lines.append(f"| {star} | {counts[('belts', star)]} | {have} | {made} | {made / have:.1f} |" if have else f"| {star} | {counts[('belts', star)]} | 0 | 0 | |")
    lines += ["", "| class | planets | " + " | ".join(star_classes) + " |", "|---|---|" + "---|" * len(star_classes)]
    for letter in classes:
        lines.append(f"| {letter} | {counts[('class', letter)]} | " + " | ".join(str(counts[('class_by_star', f'{star}/{letter}')]) for star in star_classes) + " |")
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
    print(write_report(format_report(run(), f"Seeds 1 to {SEEDS}, default configuration.")))
