# tests/audit/prevalence_occurrence.py

"""
Planet occurrence against observation (TEST.395, part E of the prevalence suite): the planets the generator makes per star against the
rates surveys measured, read from the same seeded systems the system census draws (`tests.audit.prevalence_systems`), one group per
host type (M, K, G and F dwarfs of the census).

The ranges are not taken from a reference run. Each is an observed figure of `physics/data/planet_occurrence.json` (the paper it comes
from, READ in the source or 2ND-HAND in a paper that quotes it, with its tier, GEN.728) widened for the sample: a figure from `n` systems is
held to the published 1-sigma interval and the sampling error of `n` systems added in quadrature, two sigma either side. A range the
generator misses is a real gap against observation and carries the item that will close it. The checks per host group:

* `inner_fraction` the share of stars with a planet of 1 to 20 Earth radii and a period under 400 days [R4812, R4898];
* `inner_multiplicity` the mean number of such planets among those stars (the same papers);
* `planets_per_star` the mean number of such planets per star [R4897];
* `small_planets_per_star` planets of 1 to 4 Earth radii with a period under 200 days, per star ([R1156], M dwarfs);
* `giants_2_8_au` the share of stars with a giant (0.1 to 20 Jupiter masses) at 2 to 8 AU [R1479];
* `hot_jupiter_rate` the share with a planet of 8 to 20 Earth radii and a period under 10 days [R1504, R4779, R4897];
* `earth_in_hz` planets of 0.5 to 1.5 Earth radii inside the star's habitable zone, per star [R0572];
* `pairs_over_10_hill_radii` the share of adjacent pairs at least 10 mutual Hill radii apart [R4660].

A planet's period comes from its distance and the star's mass, and its radius from the planet's own, so the figures follow the same
cuts the surveys used. The count of planets per star by mass bin of the class census stays in `tests/test_prevalence_census.py`.
"""

import json
import math
import pathlib

from planetgen.physics import constants, planet_occurrence as po
from tests.audit import prevalence as census
from tests.audit import prevalence_systems as systems

EXPECTATIONS_PATH = pathlib.Path(__file__).resolve().parent.parent / "data" / "prevalence_occurrence_expectations.json"
REPORT_NAME = "occurrence-census-latest.md"

GROUPS = {"M": ("M2V",), "K": ("K5V",), "G": ("G2V",), "F": ("F5V",)}
"""dict: The host groups, each a pool of the census star types (the keys are the groups of `planet_occurrence.HOSTS`)."""

SHARE_MEASURES = ("inner_fraction", "giants_2_8_au", "hot_jupiter_rate", "pairs_over_10_hill_radii")
"""tuple: The measures that are a share of something (a binomial error); the rest are a mean (a Poisson error)."""

MEASURES = ("inner_fraction", "inner_multiplicity", "planets_per_star", "small_planets_per_star", "giants_2_8_au", "hot_jupiter_rate",
            "earth_in_hz", "pairs_over_10_hill_radii")

BAND_MEASURES = {
    "inner_fraction": "inner_fraction",
    "inner_multiplicity": "inner_multiplicity",
    "planets_per_star": "planets_per_star_1_20re_p400",
    "small_planets_per_star": "planets_per_star_1_4re_p200",
    "giants_2_8_au": "giants_2_8_au",
    "hot_jupiter_rate": "hot_jupiter_rate",
    "earth_in_hz": "earth_mass_in_conservative_hz",
    "pairs_over_10_hill_radii": "pairs_at_least_10_hill_radii_apart",
}
"""dict: Each measure of this census and the measure of `planet_occurrence` whose bands it is checked against."""

GIANT_MIN_EARTH_MASSES = 31.8
"""float: 0.1 Jupiter masses, the lightest planet of [R1479]'s giants."""
SIGMAS = 2.0
SAMPLING_PAD = 0.05
"""float: The share of the observed value added to the range on both sides (the generator is not the survey: its planets are the
intrinsic ones, the surveys' were corrected for detection)."""

GAP_ITEMS = {
    "inner_fraction": ("GEN.728", "stars with a Kepler-like planet are 60 to 70 percent against 30 observed (the count is a uniform integer, "
                                   "not an occurrence law)"),
    "inner_multiplicity": ("GEN.728", "the number of Kepler-like planets of a host is not drawn from a multiplicity law"),
    "planets_per_star": ("GEN.728", "planets of 1 to 20 Earth radii per star are not the observed number"),
    "small_planets_per_star": ("GEN.728", "planets of 1 to 4 Earth radii per M dwarf are not the observed number"),
    "giants_2_8_au": ("GEN.920", "giants at 2 to 8 AU round FGK stars are rarer than the California Legacy Survey finds"),
    "hot_jupiter_rate": ("GEN.828", "no giant moves in to a hot Jupiter orbit; the observed rate is 0.3 to 1.2 percent of stars"),
    "earth_in_hz": ("GEN.728", "Earth-size planets in the habitable zone are not the observed number"),
    "pairs_over_10_hill_radii": ("GEN.729", "neighbours are spaced by a uniform ratio of 1.3 to 1.8, not by the observed spacing"),
}
"""dict: The gaps against observation known today, each with the item that will close it (tagged when the stored file is built)."""


def _group_of(star):
    return next((name for name, members in GROUPS.items() if star in members), None)


def _period_days(distance_au, star_mass_solar):
    return 365.25 * math.sqrt(distance_au ** 3 / star_mass_solar)


def _planets(record):
    """The planets of a record as `(mass in Earth masses, distance in AU, radius in Earth radii, period in days)`."""
    mass = record.get("star_mass") or 1.0
    return [(m, a, r, _period_days(a, mass)) for m, a, _cls, r in record["orbits"]]


def _hill_spacings(record):
    mass_kg = (record.get("star_mass") or 1.0) * constants.SOLAR_MASS_TO_KG
    found = sorted((a, m * constants.EARTH_MASS_TO_KG) for m, a, _cls, _r in record["orbits"])
    out = []
    for (a1, m1), (a2, m2) in zip(found, found[1:]):
        r_h = 0.5 * (a1 + a2) * ((m1 + m2) / (3.0 * mass_kg)) ** (1.0 / 3.0)
        out.append((a2 - a1) / r_h if r_h > 0 else math.inf)
    return out


def summarize(records):
    """Per host group, the per-star figures of the pooled systems, in the shape `prevalence_systems.compare` reads."""
    out = {}
    for group, members in GROUPS.items():
        mine = [r for r in records if r["star"] in members]
        kept = [r for r in mine if not r.get("abandoned") and not r.get("error")]
        n = len(kept)
        values = {}
        if n:
            inner = [[p for p in _planets(r) if 1.0 <= p[2] <= 20.0 and p[3] < 400.0] for r in kept]
            hosts = [x for x in inner if x]
            values["inner_fraction"] = len(hosts) / n
            values["inner_multiplicity"] = sum(len(x) for x in hosts) / len(hosts) if hosts else 0.0
            values["planets_per_star"] = sum(len(x) for x in inner) / n
            values["small_planets_per_star"] = sum(1 for r in kept for p in _planets(r) if 1.0 <= p[2] <= 4.0 and p[3] < 200.0) / n
            values["giants_2_8_au"] = sum(1 for r in kept if any(p[0] >= GIANT_MIN_EARTH_MASSES and 2.0 <= p[1] <= 8.0 for p in _planets(r))) / n
            values["hot_jupiter_rate"] = sum(1 for r in kept if any(8.0 <= p[2] <= 20.0 and p[3] < 10.0 for p in _planets(r))) / n
            values["earth_in_hz"] = sum(1 for r in kept for p in _planets(r)
                                        if 0.5 <= p[2] <= 1.5 and r["hz"][0] <= p[1] <= r["hz"][1]) / n
            pairs = [k for r in kept for k in _hill_spacings(r)]
            if pairs:
                values["pairs_over_10_hill_radii"] = sum(1 for k in pairs if k >= 10.0) / len(pairs)
        out[group] = {"n": n, "abandoned": sum(1 for r in mine if r.get("abandoned")),
                      "errors": sorted({r["error"] for r in mine if r.get("error")}),
                      "error_count": sum(1 for r in mine if r.get("error")), "values": values, "sd": {}}
    return out


def _sampling_sigma(measure, value, n):
    if measure in SHARE_MEASURES:
        return po.share_sigma(min(max(value, 0.0), 1.0), n)
    return math.sqrt(max(value, 1.0 / n) / n)


def _range(measure, bands, n):
    """The `(low, high, source)` a run of `n` systems must fall in: the widest of the bands' two-sigma intervals, each with the sampling
    error of `n` systems added in quadrature, and `SAMPLING_PAD` of the value either side."""
    low, high, sources = math.inf, -math.inf, []
    for band in bands:
        sampling = _sampling_sigma(measure, band.value, n)
        down, up = math.hypot(band.minus, sampling), math.hypot(band.plus, sampling)
        pad = SAMPLING_PAD * band.value
        low = min(low, band.value - SIGMAS * down - pad)
        high = max(high, band.value + SIGMAS * up + pad)
        sources.append(f"{band.definition} ({band.source}; {band.label}, {band.tier})")
    low = max(0.0, low)
    if measure in SHARE_MEASURES:
        high = min(1.0, high)
    return round(low, 4), round(high, 4), "; or ".join(sources)


def build_expectations(systems_per_type):
    """The expectations file: the observed ranges with their kind and source for a run of `systems_per_type` systems per host group,
    the anomalies empty until `tag_current`."""
    out = {"reference": {"kind": "observed rates of planet_occurrence.json (READ or 2ND-HAND, tier T1 or T2)",
                         "systems_per_type": systems_per_type,
                         "seeds": "the system census's: crc32 of 'system:<star>' (tests.audit.prevalence.seed_for)",
                         "groups": {name: list(members) for name, members in GROUPS.items()}},
           "anomalies": {}, "hangs": {"anomaly": ""}, "errors": {"anomaly": ""}, "types": {}}
    for group in GROUPS:
        measures = {}
        for measure in MEASURES:
            bands = [b for b in po.bands_of(BAND_MEASURES[measure], group) if b.host == group or (b.host == "all" and group in ("G", "K"))]
            if measure == "giants_2_8_au" and group not in ("G", "K"):
                continue
            if bands:
                low, high, source = _range(measure, bands, systems_per_type)
                measures[measure] = {"min": low, "max": high, "kind": "observed", "source": source}
        if measures:
            out["types"][group] = {"measures": measures}
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


def main(argv=None):
    import argparse

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
    notes, tagged_items = {}, set()
    for measure, (item, note) in GAP_ITEMS.items():
        tagged = tag_current(expectations, records, item, note, only=lambda f, m=measure: f.measure == m)
        print(measure, f"{tagged} host groups tagged {item}" if tagged else "meets the observed range now")
        if tagged:
            tagged_items.add(item)
            notes.setdefault(item, []).append(note)
    expectations["anomalies"] = {item: "; ".join(parts) for item, parts in notes.items()}
    assert set(expectations["anomalies"]) == tagged_items
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(expectations, handle, indent=1, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
