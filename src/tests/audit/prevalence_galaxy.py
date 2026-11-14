# tests/audit/prevalence_galaxy.py

"""
Planet occurrence on the galaxy's own stars (TEST.395, the galaxy-placed repeat of part E of the prevalence suite).

`prevalence_occurrence` draws a fixed number of systems round each census star type, so its host mix is the census's, not the sky's.
This module draws seeded sectors at the reference solar position the way a galaxy sector is filled (`census.fill_region_sector`,
`run_sector.generate_sector` with the position's profile), so the host mix is whatever the galaxy puts there, and holds the planets of
the stars it placed to the same observed bands (`planet_occurrence`, each with its paper and tier) on their final classes and stored
life.

How a band is read for a mix that is not the census's:

* A host group (M, K, G, F by the star's mass, `planet_occurrence.host_group`) with `MIN_GROUP` stars or more is a cell of its own, held
  to the group's bands widened for the sampling error of its own `n` (`prevalence_occurrence.build_expectations`).
* The `galaxy` cell pools every group that has a band for the measure: its value is the pooled figure of those groups' stars, and its range
  is the groups' ranges averaged with the weight of their stars. So a mix that is mostly M dwarfs is held mostly to the M dwarf rates, and a
  group too small to bind alone still counts in the pool in proportion to its stars.
* A flag whose group carries a known-gap item in the stored occurrence file carries the same item here (the cause is the generator's, not
  the sample's); a flag in the pool is known when any group it pools carries one for that measure.

A flag names its cell and measure (the bin that left its band). The report also lists the class mix of the sample (planets and moons by final
class, the habitable-class bodies per star counted from the stored life class), report only: no survey measures a class.
Classes O, P and Q have no target rate (Boss, 2026-10-30: physics drives, no game rates); they are tested for being possible when conditions
allow them, not for a share.
"""

import collections
import math
import multiprocessing
import os
import signal

from planetgen.physics import planet_occurrence as po
from tests.audit import census as audit
from tests.audit import prevalence as census
from tests.audit import prevalence_occurrence as occurrence
from tests.audit import prevalence_stars as stars
from tests.audit import prevalence_systems as systems

REPORT_NAME = "occurrence-galaxy-latest.md"
POOL = "galaxy"
"""str: The cell name of the pooled figure."""
SECTORS = 100
"""int: Seeded sectors drawn at the solar position (about three systems each at the default galaxy)."""
MIN_GROUP = 25
"""int: Stars a host group needs to be a cell of its own; fewer only count in the pool."""
SECTOR_TIME_LIMIT_S = 900
"""int: Seconds one sector may take before the census abandons it and counts it (a generation hang)."""

REPRESENTATIVE = {group: members[0] for group, members in occurrence.GROUPS.items()}
"""dict: The census star type a group's records are filed under, so `prevalence_occurrence.summarize` reads them."""


class _Abandoned(BaseException):
    """Raised by the alarm when one sector ran past `SECTOR_TIME_LIMIT_S`."""


def sector_case():
    """The one region the census fills: a sector at the reference solar position, its system count drawn from the density."""
    from tests.audit import representation as rep

    return audit.RegionCase("solar", stars.solar_position_pc(), rep.WITH_COMPONENTS, rep.EDGE, rep.E1, sectors=SECTORS)


def class_counts(system):
    """The final class of every planet and moon of `system` as a `{letter: count}` (a belt counts as `asteroid`)."""
    found = collections.Counter()
    for planet in system.all_planets:
        if planet.body_type == "a":
            found["asteroid"] += 1
            continue
        for body in (planet, *planet.moons):
            found[body.planet_class or "none"] += 1
    return dict(found)


def record_of(entry):
    """The record of one placed system: the system census's record, its host group, its final classes and the bodies with life."""
    system = entry.star_system
    record = systems._record(system)
    group = po.host_group(record["star_mass"])
    record.update(group=group, star=REPRESENTATIVE.get(group, "other"), classes=class_counts(system),
                  life_bodies=system.count_habitable_class_worlds())
    return record


def draw_sector(index):
    """The records of the systems of sector `index`; a sector that ran over the limit is one `{"abandoned": True}`."""
    def _alarm(*_):
        raise _Abandoned()

    previous = signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(SECTOR_TIME_LIMIT_S)
    try:
        sector = audit.fill_region_sector(sector_case(), audit.seed_for("solar", index))
        return [record_of(entry) for entry in sector.entries if getattr(entry, "star_system", None) is not None]
    except _Abandoned:
        return [{"abandoned": True, "star": "other", "group": None}]
    except Exception as error:  # a crash is a finding, not the end of the census
        import pathlib
        import traceback

        frame = traceback.extract_tb(error.__traceback__)[-1]
        return [{"error": f"{type(error).__name__}: {error} ({pathlib.Path(frame.filename).name}:{frame.lineno})", "star": "other", "group": None}]
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def draw_sectors(count=SECTORS, workers=None):
    """The records of sectors 0 to `count - 1`, on worker processes (the sector alarm must not be shared with a test run's own timeout)."""
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, count)) as pool:
        results = pool.map(draw_sector, range(count), chunksize=1)
    return [record for result in results for record in result]


def summarize(records):
    """Per host group and for the pool: the figures of `prevalence_occurrence.summarize`, plus how many stars fell outside every group."""
    grouped = occurrence.summarize([r for r in records if r.get("star") in REPRESENTATIVE.values()])
    out = dict(grouped)
    out[POOL] = _pooled(grouped)
    out["other"] = {"n": sum(1 for r in records if r.get("star") == "other" and not r.get("abandoned") and not r.get("error"))}
    return out


def _pooled(grouped):
    """The pooled cell: each measure over the groups that have a band for it, the figure of those groups' stars together."""
    values, weights = {}, {}
    for measure in occurrence.MEASURES:
        members = [(name, spec) for name, spec in grouped.items() if measure in spec["values"] and _has_band(name, measure)]
        if not members:
            continue
        weight = sum(spec["n"] for _name, spec in members)
        values[measure] = sum(spec["values"][measure] * spec["n"] for _name, spec in members) / weight
        weights[measure] = {name: spec["n"] for name, spec in members}
    return {"n": sum(spec["n"] for spec in grouped.values()), "abandoned": sum(spec["abandoned"] for spec in grouped.values()),
            "errors": sorted({e for spec in grouped.values() for e in spec["errors"]}),
            "error_count": sum(spec["error_count"] for spec in grouped.values()), "values": values, "sd": {}, "weights": weights}


def _has_band(group, measure):
    """Whether the observed bands give `group` a range for `measure` (the same test `prevalence_occurrence.build_expectations` makes)."""
    return measure in occurrence.build_expectations(1)["types"].get(group, {"measures": {}})["measures"]


def expectations_for(summary, stored):
    """The ranges the sample is held to: a cell per group with `MIN_GROUP` stars or more at its own `n`, and the pooled cell with the groups'
    ranges averaged by their stars; known-gap items copied from the `stored` occurrence file."""
    out = {"reference": {"kind": "observed rates of planet_occurrence.json, galaxy-placed stars", "sectors": SECTORS, "min_group": MIN_GROUP},
           "anomalies": dict(stored["anomalies"]), "hangs": {"anomaly": stored["hangs"]["anomaly"]},
           "errors": {"anomaly": stored["errors"]["anomaly"]}, "types": {}}
    ranges = {}
    for group in occurrence.GROUPS:
        n = summary[group]["n"]
        if not n:
            continue
        built = occurrence.build_expectations(n)["types"].get(group, {"measures": {}})["measures"]
        ranges[group] = built
        if n >= MIN_GROUP:
            measures = {}
            for name, bounds in built.items():
                kept = dict(bounds)
                gap = stored["types"].get(group, {"measures": {}})["measures"].get(name, {}).get("anomaly")
                if gap:
                    kept["anomaly"] = gap
                measures[name] = kept
            out["types"][group] = {"measures": measures}
    pooled = {}
    for measure, weights in summary[POOL].get("weights", {}).items():
        have = [(group, n) for group, n in weights.items() if measure in ranges.get(group, {})]
        if not have:
            continue
        whole = sum(n for _g, n in have)
        bounds = {"min": round(sum(ranges[g][measure]["min"] * n for g, n in have) / whole, 4),
                  "max": round(sum(ranges[g][measure]["max"] * n for g, n in have) / whole, 4), "kind": "observed",
                  "source": "the host groups' bands averaged by their stars: " + "; ".join(f"{g} {n}" for g, n in have)}
        gaps = sorted({stored["types"].get(g, {"measures": {}})["measures"].get(measure, {}).get("anomaly") for g, _n in have} - {None})
        if gaps:
            bounds["anomaly"] = gaps[0]
        pooled[measure] = bounds
    if pooled:
        out["types"][POOL] = {"measures": pooled}
    return out


def compare(summary, expectations):
    """Flags over the sample: the shared comparison, plus a group-less sample (every star outside the groups) as a missing cell."""
    return systems.compare({k: v for k, v in summary.items() if k != "other"}, expectations)


def class_table(records):
    """`(class letter, bodies, share)` rows of the sample's final classes, most common first, and the stars they are of."""
    found = collections.Counter()
    stars_counted = 0
    for record in records:
        if record.get("abandoned") or record.get("error"):
            continue
        stars_counted += 1
        found.update(record["classes"])
    whole = sum(found.values()) or 1
    return [(letter, count, count / whole) for letter, count in found.most_common()], stars_counted


def format_report(summary, records, findings, header=""):
    lines = ["# Planet occurrence on galaxy-placed stars (TEST.395), latest run", ""]
    if header:
        lines += [header, ""]
    unexplained = [f for f in findings if not f.known]
    explained = [f for f in findings if f.known]
    lines += [f"Flags: {len(unexplained)} new, {len(explained)} known gaps.", ""]
    lines += [f"- {finding.line()}" for finding in unexplained]
    if explained:
        lines += ["", "Known gaps against observation:", ""]
        lines += [f"- {finding.line()}" for finding in explained]
    lines += ["", "## Figures per host group and pooled", ""]
    for name in (*occurrence.GROUPS, POOL):
        mine = summary.get(name)
        if not mine or "values" not in mine:
            continue
        figures = ", ".join(f"{m} {mine['values'][m]:.3f}" for m in occurrence.MEASURES if m in mine["values"])
        note = "" if name == POOL or mine["n"] >= MIN_GROUP else f" (under {MIN_GROUP}: counted in the pool only)"
        lines.append(f"- {name} (n={mine['n']}{note}): {figures}")
    lines.append(f"- stars outside the host groups: {summary.get('other', {}).get('n', 0)}")
    table, counted = class_table(records)
    life = sum(r["life_bodies"] for r in records if "life_bodies" in r)
    lines += ["", f"## Final classes of the {counted} systems (report only)", "",
              f"Bodies with life, from the stored life class: {life} ({life / max(counted, 1):.3f} per star).", "",
              "| class | bodies | share |", "|---|---|---|"]
    lines += [f"| {letter} | {count} | {share * 100:.1f} % |" for letter, count, share in table]
    return "\n".join(lines) + "\n"
