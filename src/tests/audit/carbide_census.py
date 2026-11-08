# tests/audit/carbide_census.py

"""
The carbide census (GEN.1423, the Done line: "both minerals appear on carbon-rich worlds in a census sample and never on oxygen-rich ones").

Whole seeded systems are generated, `PER_GROUP` around stars whose C/O is at or above `bulk_draw.CARBON_OVER_OXYGEN` (the carbon-rich group) and `PER_GROUP` around the rest
(the oxygen-rich group). Seeds are screened on the star alone (a star is drawn before any planet, so the cheap screen names the group of the whole system), but each system
is put in the group of its own finished star, so a screen that missed would only move a system between the two columns. Every planet and moon the ore pass assessed is counted by
whether its rock is carbon rock and what status `physics/carbide_deposits.py` gave each of its two rows (`moissanite_crust`, `cohenite_metal`).

The rule the census checks (`violations`): a body whose rock is silicate has both rows `absent_by_state`; a body around a star under C/O 0.8 never has carbon rock; every carbon-rock
body the ore pass assessed has both rows; and, over the sample, each mineral is present on at least one carbon-rock body.
"""

import collections
import multiprocessing
import os
import pathlib

PER_GROUP = 40
"""int: Systems drawn in each group."""
SCREEN_LIMIT = 6000
"""int: The most seeds looked at to fill the groups."""
SYSTEM_TIME_LIMIT_S = 240
REPORT_NAME = "carbide-census-latest.md"
MINERALS = ("moissanite_crust", "cohenite_metal")
STATUSES = ("present", "absent_by_chance", "absent_by_state")


def _star_c_over_o(seed):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import bulk_draw
    from planetgen.util import draw

    draw.set_run_seed(seed)
    config = SystemConfig()
    config.PLANETS = False
    return bulk_draw.host_ratios(StarSystem(system_config=config).primary_star).c_over_o


def screen(per_group=PER_GROUP, limit=SCREEN_LIMIT):
    """`(carbon seeds, oxygen seeds)`: the first seeds, counting up from 1, whose star is at or above C/O 0.8 and the first that are not, `per_group` of each."""
    from planetgen.physics import bulk_draw

    carbon, oxygen = [], []
    for seed in range(1, limit + 1):
        if len(carbon) >= per_group and len(oxygen) >= per_group:
            break
        rich = _star_c_over_o(seed) >= bulk_draw.CARBON_OVER_OXYGEN
        group = carbon if rich else oxygen
        if len(group) < per_group:
            group.append(seed)
    return carbon, oxygen


def system_counts(seed):
    """Counts for one finished system: `{(group, carbon_rock, mineral, status): bodies, (group, carbon_rock, "bodies"): bodies}` with `group` `"carbon"` or `"oxygen"` by its own star."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import bulk_draw, ore_deposits
    from planetgen.util import draw

    draw.set_run_seed(seed)
    system = StarSystem(system_config=SystemConfig())
    rich = bulk_draw.host_ratios(system.primary_star).c_over_o >= bulk_draw.CARBON_OVER_OXYGEN
    group = "carbon" if rich else "oxygen"
    counts = collections.Counter()
    for planet in system.planets:
        for body in [planet, *(getattr(planet, "moons", None) or [])]:
            rows = getattr(body, ore_deposits.ATTRIBUTE, None)
            if rows is None:
                continue
            rock = bool(getattr(body, "carbon_rock", False))
            counts[(group, rock, "bodies")] += 1
            for row in rows:
                if row["deposit_type"] in MINERALS:
                    counts[(group, rock, row["deposit_type"], row["status"])] += 1
    counts[(group, None, "systems")] += 1
    return counts


def _job(seed):
    return system_counts(seed)


def run(per_group=PER_GROUP, workers=None, limit=SCREEN_LIMIT):
    """The summed counts of `per_group` systems in each group (see `system_counts`)."""
    carbon, oxygen = screen(per_group, limit)
    seeds = carbon + oxygen
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, max(len(seeds), 1))) as pool:
        parts = pool.map(_job, seeds, chunksize=1)
    total = collections.Counter()
    for part in parts:
        total.update(part)
    return total


def violations(counts):
    """What breaks the rule over `counts` (a list of sentences; empty when it holds)."""
    found = []
    for key, number in sorted(counts.items(), key=str):
        if len(key) == 4:
            group, rock, mineral, status = key
            if not rock and status != "absent_by_state":
                found.append(f"{number} silicate bod{'y' if number == 1 else 'ies'} in the {group} group had {mineral} {status}")
    for group, rock, kind in [key for key in counts if len(key) == 3 and key[2] == "bodies"]:
        if group == "oxygen" and rock:
            found.append(f"{counts[(group, rock, kind)]} carbon-rock bodies around stars under C/O 0.8")
        if rock:
            for mineral in MINERALS:
                rows = sum(counts[(group, rock, mineral, status)] for status in STATUSES)
                if rows != counts[(group, rock, kind)]:
                    found.append(f"{counts[(group, rock, kind)] - rows} carbon-rock bodies without a {mineral} row")
    for mineral in MINERALS:
        if not any(number for key, number in counts.items() if len(key) == 4 and key[1] and key[2] == mineral and key[3] == "present"):
            found.append(f"{mineral} is present on no carbon-rock body of the sample")
    return found


def format_report(counts, header=""):
    """The census as Markdown."""
    lines = ["# Carbide census (GEN.1423)", "", header, "", "| group | systems | assessed bodies | carbon-rock bodies | moissanite present | cohenite present | silicate bodies with a carbide present |",
             "|---|---|---|---|---|---|---|"]
    for group in ("carbon", "oxygen"):
        rock = lambda flag, *rest: counts[(group, flag, *rest)]
        bodies = rock(True, "bodies") + rock(False, "bodies")
        lines.append(f"| {group}-rich stars | {counts[(group, None, 'systems')]} | {bodies} | {rock(True, 'bodies')} | {rock(True, 'moissanite_crust', 'present')} of {rock(True, 'bodies')} | "
                     f"{rock(True, 'cohenite_metal', 'present')} of {rock(True, 'bodies')} | {sum(rock(False, m, 'present') for m in MINERALS)} |")
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
