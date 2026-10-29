# tests/audit/prevalence_occurrence.py

"""
Planet occurrence against observation (TEST.395, part E of the prevalence suite): the planets the generator makes per star, by
mass and orbit, against the rates surveys measured, read from the same seeded systems the system census draws
(`tests.audit.prevalence_systems`), grouped by host (M dwarfs, K/G/F dwarfs).

Unlike the other parts the ranges are not taken from a reference run: each is an observed rate with its source and tier, widened
for the sample (a host group is 100 to 300 systems). A range the generator misses is a real gap against observation and carries
the item that will close it (GEN.828 for giants and hot Jupiters, GEN.829 for the sub-Neptune population). The observed rates
here were recalled by Aux 16 (reports/planet-class-distribution-check-2026-10-25.md) and are marked T2/T3, unverified, until a
lane has read the papers; Fressin 2013, Petigura 2018, Zhu 2018, Fulton 2021, Bryson 2021, Wright 2012, Johnson 2010,
Bonfils 2013, Dressing and Charbonneau 2015.

Measures are per star of the group unless named a ratio: `giants` (100 Earth masses or more, any orbit), `hot_jupiters` (the same
inside 0.1 AU), `cold_giants` (1 to 10 AU), `inner_planets` (1 to 20 Earth masses inside 0.5 AU), `earth_mass_in_hz` (0.5 to 2 Earth
masses inside the star's ecosphere) and `super_earths_per_neptune` (planets of 2 to 10 over planets of 10 to 50 Earth masses, which
observation has above one).
"""

import json
import pathlib

from tests.audit import prevalence as census
from tests.audit import prevalence_systems as systems

EXPECTATIONS_PATH = pathlib.Path(__file__).resolve().parent.parent / "data" / "prevalence_occurrence_expectations.json"
REPORT_NAME = "occurrence-census-latest.md"

GROUPS = {"M": ("M2V",), "KGF": ("K5V", "G2V", "F5V")}
"""dict: The host groups: M dwarfs and K, G and F dwarfs, each a pool of the census star types."""

MEASURES = ("giants", "hot_jupiters", "cold_giants", "inner_planets", "earth_mass_in_hz", "super_earths_per_neptune")

GIANT_ME = 100.0

OBSERVED = {
    "M": {
        "giants": (0.0, 0.12, "giants about 0.03 to 0.06 per M dwarf (Johnson et al. 2010, Bonfils et al. 2013), widened for 100 stars"),
        "hot_jupiters": (0.0, 0.01, "hot Jupiters under 0.002 per M dwarf (Johnson et al. 2010; transit surveys), widened for 100 stars"),
    },
    "KGF": {
        "giants": (0.08, 0.35, "giants of 100 Earth masses or more about 0.15 to 0.25 per FGK star (Wright et al. 2012, Fulton et al. 2021)"),
        "hot_jupiters": (0.0, 0.02, "hot Jupiters about 0.005 to 0.01 per FGK star (Wright et al. 2012; Fressin et al. 2013: 0.0043), "
                                    "widened for 300 stars"),
        "cold_giants": (0.05, 0.25, "giants at 1 to 10 AU about 0.10 to 0.15 per FGK star (Fulton et al. 2021, 0.3 to 10 AU)"),
        "inner_planets": (0.4, 1.6, "1 to 20 Earth-mass planets inside 0.5 AU about 0.5 to 1 per FGK star (Fressin et al. 2013, "
                                    "Zhu et al. 2018, Petigura et al. 2018)"),
        "earth_mass_in_hz": (0.12, 0.8, "Earth-mass planets in the ecosphere (eta-Earth) 0.2 to 0.6 per Sun-like star (Bryson et al. "
                                        "2021, Petigura et al. 2018), widened for 300 stars"),
        "super_earths_per_neptune": (1.0, 12.0, "planets of 2 to 10 Earth masses outnumber those of 10 to 50: the mass function falls "
                                                "steadily above about 10 Earth masses (Fressin et al. 2013, Petigura et al. 2018)"),
    },
}
"""dict: `group -> measure -> (min, max, source)`: observed rates, kind `observed`, tier T2/T3 (recalled, to be checked against the
papers). A host group with no entry for a measure is not checked for it."""


def _group_of(star):
    return next((name for name, members in GROUPS.items() if star in members), None)


def summarize(records):
    """Per host group, the per-star figures of the pooled systems, in the shape `prevalence_systems.compare` reads."""
    out = {}
    for group, members in GROUPS.items():
        mine = [r for r in records if r["star"] in members]
        kept = [r for r in mine if not r.get("abandoned") and not r.get("error")]
        n = len(kept)
        values = {}
        if n:
            planets = [(mass, a, cls, r["hz"]) for r in kept for mass, a, cls in r["orbits"]]
            giants = [p for p in planets if p[0] >= GIANT_ME]
            values["giants"] = len(giants) / n
            values["hot_jupiters"] = sum(1 for p in giants if p[1] < 0.1) / n
            values["cold_giants"] = sum(1 for p in giants if 1.0 <= p[1] <= 10.0) / n
            values["inner_planets"] = sum(1 for p in planets if 1.0 <= p[0] <= 20.0 and p[1] < 0.5) / n
            values["earth_mass_in_hz"] = sum(1 for m, a, c, hz in planets if 0.5 <= m <= 2.0 and hz[0] <= a <= hz[1]) / n
            nep = sum(1 for p in planets if 10.0 < p[0] <= 50.0)
            values["super_earths_per_neptune"] = (sum(1 for p in planets if 2.0 <= p[0] <= 10.0) / nep) if nep else 99.0
        out[group] = {"n": n, "abandoned": sum(1 for r in mine if r.get("abandoned")),
                      "errors": sorted({r["error"] for r in mine if r.get("error")}),
                      "error_count": sum(1 for r in mine if r.get("error")), "values": values, "sd": {}}
    return out


def build_expectations(systems_per_type):
    """The expectations file: the observed ranges with their kind and source, the anomalies empty until `tag_current`."""
    out = {"reference": {"kind": "observed rates (recalled by Aux 16, tier T2/T3, unverified)", "systems_per_type": systems_per_type,
                         "seeds": "the system census's: crc32 of 'system:<star>' (tests.audit.prevalence.seed_for)",
                         "groups": {name: list(members) for name, members in GROUPS.items()}},
           "anomalies": {}, "hangs": {"anomaly": ""}, "errors": {"anomaly": ""}, "types": {}}
    for group, measures in OBSERVED.items():
        out["types"][group] = {"measures": {name: {"min": low, "max": high, "kind": "observed", "source": source}
                                            for name, (low, high, source) in measures.items()}}
    return out


def compare(summary, expectations):
    return systems.compare(summary, expectations)


def tag_current(expectations, records, item, note, only=None):
    tagged = 0
    for finding in compare(summarize(records), expectations):
        if finding.kind in ("hang", "error"):
            continue
        bounds = expectations["types"][finding.cell]["measures"][finding.measure]
        if bounds.get("anomaly") or (only and not only(finding)):
            continue
        bounds["anomaly"] = item
        tagged += 1
    expectations["anomalies"][item] = note
    return tagged


def load_expectations(path=EXPECTATIONS_PATH):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def format_report(summary, findings, header=""):
    lines = ["# Planet occurrence against observation (TEST.395), latest run", ""]
    if header:
        lines += [header, ""]
    unexplained = [f for f in findings if not f.known]
    explained = [f for f in findings if f.known]
    lines += [f"Flags: {len(unexplained)} new, {len(explained)} known gaps.", ""]
    lines += [f"- {finding.line()}" for finding in unexplained]
    if explained:
        lines += ["", "Known gaps against observation:", ""]
        lines += [f"- {finding.line()}" for finding in explained]
    lines += ["", "## Figures per host group", ""]
    for group, mine in summary.items():
        figures = ", ".join(f"{name} {mine['values'][name]:.3f}" for name in MEASURES if name in mine["values"])
        lines.append(f"- {group} ({', '.join(GROUPS[group])}; n={mine['n']}): {figures}")
    return "\n".join(lines) + "\n"


GAP_ITEMS = {
    ("M", "giants"): ("GEN.828", "giants per M dwarf are about 0.3 to 0.5 against 0.03 to 0.06 observed: giant occurrence does not fall "
                                 "with the host mass"),
    ("M", "hot_jupiters"): ("GEN.828", "hot Jupiters per M dwarf are about 0.07 to 0.12 against under 0.002 observed"),
    ("KGF", "inner_planets"): ("GEN.829", "1 to 20 Earth-mass planets inside 0.5 AU per FGK star are about 0.1 against 0.5 to 1 "
                                          "observed: the sub-Neptune population is missing"),
    ("KGF", "super_earths_per_neptune"): ("GEN.829", "planets of 10 to 50 Earth masses outnumber those of 2 to 10: the mass function "
                                                      "is inverted against observation"),
}
"""dict: The gaps against observation known today, each with the item that will close it (the first run found them)."""


def main(argv=None):
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Writes tests/data/prevalence_occurrence_expectations.json (TEST.395).")
    parser.add_argument("records", help="system records of the code now (prevalence_systems measure)")
    parser.add_argument("--systems-per-type", type=int, required=True)
    parser.add_argument("--out", default=str(EXPECTATIONS_PATH))
    parser.add_argument("--hang-anomaly", default="")
    parser.add_argument("--error-anomaly", default="")
    args = parser.parse_args(argv)
    with open(args.records, encoding="utf-8") as handle:
        records = json.load(handle)
    expectations = build_expectations(args.systems_per_type)
    expectations["hangs"]["anomaly"] = args.hang_anomaly
    expectations["errors"]["anomaly"] = args.error_anomaly
    for (group, measure), (item, note) in GAP_ITEMS.items():
        tagged = tag_current(expectations, records, item, note, only=lambda f, g=group, m=measure: (f.cell, f.measure) == (g, m))
        print(group, measure, "tagged" if tagged else "meets the observed range now")
        if not tagged:
            expectations["anomalies"].pop(item, None)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(expectations, handle, indent=1, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
