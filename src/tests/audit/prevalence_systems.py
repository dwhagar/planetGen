# tests/audit/prevalence_systems.py

"""
The system census (TEST.431, part B of the prevalence suite): seeded full systems of each census star type, summed into the
figures a whole system decides and a planet draw cannot: planets, moons, belts and comets per system, the share with no planet,
the share that are pairs, the share with a habitable-class world or a Blue or Green one, and the share where life reaches each
milestone. Checked against `tests/data/prevalence_systems_expectations.json`, whose ranges are marked `observed`, `norm` or
`stand-in` exactly as the class census's are (`tests/audit/prevalence.py`).

A measure is a `share` (a fraction of the systems) or a `mean` (per system, or per planet for moons); a `share` range is a
binomial band and a `mean` range is a band on the mean from the reference's own spread, so a repeat of the reference passes by a
wide margin and a real shift trips a flag. The run is seeded from the star type, a system that runs past the time limit is
counted (`hang`), and the draws are made in worker processes like the class census's.
"""

import json
import math
import multiprocessing
import os
import pathlib
import signal

from tests.audit import prevalence as census

EXPECTATIONS_PATH = pathlib.Path(__file__).resolve().parent.parent / "data" / "prevalence_systems_expectations.json"
REPORT_NAME = "system-census-latest.md"
SYSTEM_TIME_LIMIT_S = 180
"""int: Seconds one system may take before the census abandons it and counts it."""

EARTH_MASS_KG = 5.9722e24
EARTH_RADIUS_KM = 6371.0
SOLAR_MASS_KG = 1.98847e30

MILESTONES = ("abiogenesis", "photosynthesis", "complex_cells", "multicellularity", "technological_civilization")

SHARE_MEASURES = ("zero_planets", "binary", "comets", "belts", "habitable_class", "blue_green", "life_abiogenesis",
                  "life_photosynthesis", "life_complex_cells", "life_multicellularity", "life_technological_civilization")
MEAN_MEASURES = ("planets", "belts_count", "moons_per_planet", "comets_count")

PAD_SHARE = 0.05
PAD_MEAN = 0.15
"""The widths added beyond the sampling error: an absolute 0.05 for a share, 15 percent of the mean for a mean. The class
cut-over reassigns the classes, so the ranges start wide and tighten as it settles."""
SIGMAS = 4.0
MUST_APPEAR = 5
"""int: A share expected this many times or more in a run must appear at least once (none in a run is under 1 chance in 100 when the
reference share holds), whatever its band: the band is wide, and a collapse to zero is the fault it must not miss."""


class _Abandoned(BaseException):
    """Raised by the alarm when one system ran past `SYSTEM_TIME_LIMIT_S`."""


def _record(system):
    from planetgen.generation import life  # noqa: F401  (the stage names live on the bodies)

    bodies = [body for planet in system.all_planets if planet.body_type != "a" for body in (planet, *planet.moons)]
    planets = [planet for planet in system.all_planets if planet.body_type != "a"]
    stages = {getattr(body, "life_stage", None) for body in bodies}
    reached = {}
    for milestone in MILESTONES:
        reached[milestone] = any(MILESTONES.index(stage) >= MILESTONES.index(milestone)
                                 for stage in stages if stage in MILESTONES)
    belts = sum(1 for planet in system.all_planets if planet.body_type == "a")
    comets = len(system.comets) + len(getattr(system, "secondary_comets", ()))
    hab_count, _m = system.count_habitable()
    primary = [planet for planet in system.planets if planet.body_type != "a"]
    orbits = [[round((planet.mass or 0.0) / EARTH_MASS_KG, 4), round(planet.distance, 4), planet.planet_class,
               round((getattr(planet, "radius", 0.0) or 0.0) / EARTH_RADIUS_KM, 4)] for planet in primary]
    return {"orbits": orbits, "star_mass": round(system.primary_star.mass / SOLAR_MASS_KG, 4), "hz": [round(edge, 4) for edge in system.primary_star.habitable_zone], "planets": len(planets), "belts": belts, "moons": sum(len(planet.moons) for planet in planets), "comets": comets,
            "binary": system.binary_type is not None, "habitable_class": system.count_habitable_class_worlds() > 0,
            "blue_green": hab_count > 0, "reached": reached}


def draw_type(star_type, count):
    """`count` systems around `star_type`, as a list of records; a system over the time limit is `{"abandoned": True}`."""
    from planetgen.util import draw
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem

    draw.set_run_seed(census.seed_for("system", star_type))
    records = []

    def _alarm(*_):
        raise _Abandoned()

    previous = signal.signal(signal.SIGALRM, _alarm)
    try:
        for _ in range(count):
            config = SystemConfig()
            config.STAR_TYPE = star_type
            signal.alarm(SYSTEM_TIME_LIMIT_S)
            try:
                records.append(_record(StarSystem(system_config=config)))
            except _Abandoned:
                records.append({"abandoned": True})
            except Exception as error:  # a crash is a finding, not the end of the census
                import traceback

                frame = traceback.extract_tb(error.__traceback__)[-1]
                records.append({"error": f"{type(error).__name__}: {error} ({pathlib.Path(frame.filename).name}:{frame.lineno})"})
            finally:
                signal.alarm(0)
    finally:
        signal.signal(signal.SIGALRM, previous)
    return [{"star": star_type, **record} for record in records]


def _job(job):
    return draw_type(*job)


def draw_types(count, star_types=census.STAR_TYPES, workers=None):
    """Every star type's records, on worker processes (the time guard is an alarm, which a test run's own timeout must not share)."""
    jobs = [(star, count) for star in star_types]
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, len(jobs))) as pool:
        results = pool.map(_job, jobs, chunksize=1)
    return [record for result in results for record in result]


def summarize(records):
    """Per star type: `{"n": systems, "abandoned": n, "values": {measure: value}, "sd": {mean measure: sd}}`."""
    out = {}
    for star in sorted({record["star"] for record in records}):
        mine = [r for r in records if r["star"] == star]
        kept = [r for r in mine if not r.get("abandoned") and not r.get("error")]
        errors = sorted({r["error"] for r in mine if r.get("error")})
        n = len(kept)
        values, spread = {}, {}
        if n:
            planets = sum(r["planets"] for r in kept)
            values["zero_planets"] = sum(1 for r in kept if r["planets"] + r["belts"] == 0) / n
            values["binary"] = sum(r["binary"] for r in kept) / n
            values["comets"] = sum(1 for r in kept if r["comets"] > 0) / n
            values["belts"] = sum(1 for r in kept if r["belts"] > 0) / n
            values["habitable_class"] = sum(r["habitable_class"] for r in kept) / n
            values["blue_green"] = sum(r["blue_green"] for r in kept) / n
            for milestone in MILESTONES:
                values[f"life_{milestone}"] = sum(r["reached"][milestone] for r in kept) / n
            values["planets"] = planets / n
            values["belts_count"] = sum(r["belts"] for r in kept) / n
            values["comets_count"] = sum(r["comets"] for r in kept) / n
            values["moons_per_planet"] = sum(r["moons"] for r in kept) / planets if planets else 0.0
            for measure, key in (("planets", "planets"), ("belts_count", "belts"), ("comets_count", "comets")):
                mean = values[measure]
                spread[measure] = math.sqrt(sum((r[key] - mean) ** 2 for r in kept) / max(n - 1, 1))
            moon_mean = values["moons_per_planet"]
            per_system = [(r["moons"] / r["planets"]) if r["planets"] else 0.0 for r in kept]
            spread["moons_per_planet"] = math.sqrt(sum((x - moon_mean) ** 2 for x in per_system) / max(n - 1, 1))
        out[star] = {"n": n, "abandoned": sum(1 for r in mine if r.get("abandoned")), "errors": errors,
                     "error_count": sum(1 for r in mine if r.get("error")), "values": values, "sd": spread}
    return out


def measure_range(measure, value, n_ref, n_run, sd=0.0):
    """The band a reference `value` over `n_ref` systems allows a run of `n_run`: a binomial band for a share, a band on the mean
    from the reference's spread for a mean."""
    if measure in SHARE_MEASURES:
        sigma = math.sqrt(max(value * (1.0 - value), 0.0) * (1.0 / max(n_ref, 1) + 1.0 / max(n_run, 1)))
        half = SIGMAS * sigma + PAD_SHARE
        low = max(0.0, value - half)
        if value * n_run >= MUST_APPEAR:
            low = max(low, 1.0 / n_run)  # expected five or more times: none at all is a collapse, whatever the band says
        return round(low, 4), round(min(1.0, value + half), 4)
    sigma = sd * math.sqrt(1.0 / max(n_ref, 1) + 1.0 / max(n_run, 1))
    half = SIGMAS * sigma + PAD_MEAN * abs(value)
    return round(max(0.0, value - half), 4), round(value + half, 4)



def load_expectations(path=EXPECTATIONS_PATH):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def compare(summary, expectations):
    """Every flag over the census: each expected measure of each star type against its range, and any abandoned system."""
    findings = []
    for star, spec in expectations.get("types", {}).items():
        mine = summary.get(star)
        if not mine or not mine["n"]:
            findings.append(census.Finding("missing", star, "type", 0.0, 1.0, 1.0, "no systems drawn"))
            continue
        for measure, bounds in spec.get("measures", {}).items():
            value = mine["values"].get(measure, 0.0)
            low, high = bounds["min"], bounds["max"]
            if value < low or value > high:
                kind = "missing" if value == 0.0 and low > 0 else "out-of-range"
                findings.append(census.Finding(kind, star, measure, value, low, high,
                                               f"{bounds.get('kind', '?')}: {bounds.get('source', '')}",
                                               bounds.get("anomaly", "")))
        if mine["error_count"]:
            findings.append(census.Finding("error", star, "systems that raised", float(mine["error_count"]), 0.0, 0.0,
                                           "; ".join(mine["errors"])[:300], expectations.get("errors", {}).get("anomaly", "")))
        if mine["abandoned"]:
            findings.append(census.Finding("hang", star, "systems over the time limit", float(mine["abandoned"]), 0.0, 0.0,
                                           f"a system ran past {SYSTEM_TIME_LIMIT_S} s",
                                           expectations.get("hangs", {}).get("anomaly", "")))
    return findings


def stale_anomalies(expectations, findings):
    broken = {(f.cell, f.measure) for f in findings}
    return [f"{star} {measure} ({bounds['anomaly']})" for star, spec in expectations.get("types", {}).items()
            for measure, bounds in spec.get("measures", {}).items()
            if bounds.get("anomaly") and (star, measure) not in broken]


def format_report(summary, findings, header=""):
    lines = ["# System census (TEST.431), latest run", ""]
    if header:
        lines += [header, ""]
    unexplained = [f for f in findings if not f.known]
    explained = [f for f in findings if f.known]
    lines += [f"Flags: {len(unexplained)} new, {len(explained)} known anomalies.", ""]
    lines += [f"- {finding.line()}" for finding in unexplained]
    if explained:
        lines += ["", "Known anomalies:", ""]
        lines += [f"- {finding.line()}" for finding in explained]
    lines += ["", "## Figures per star type", ""]
    for star in census.STAR_TYPES:
        mine = summary.get(star)
        if not mine:
            continue
        figures = ", ".join(f"{name} {mine['values'][name]:.3f}" for name in (*MEAN_MEASURES, *SHARE_MEASURES)
                            if name in mine["values"])
        lines.append(f"- {star} (n={mine['n']}, abandoned {mine['abandoned']}, raised {mine['error_count']}): {figures}")
    return "\n".join(lines) + "\n"


# --- The reference and the ranges ------------------------------------------------------------------------------------

def build_expectations(reference_records, commit, label, n_run):
    """Ranges from the reference records (the code before the change being tested): shares and means as `measure_range`, all
    stand-ins. A researched norm needs a sample large enough to bind (a binary share at 100 systems is good to about 0.2 either
    way), so the norms of multiplicity, star type and the like are in the star census (TEST.432), which draws thousands of stars."""
    summary = summarize(reference_records)
    out = {"reference": {"commit": commit, "label": label, "systems_per_type": n_run,
                         "seeds": "crc32 of 'system:<star>' (tests.audit.prevalence.seed_for)",
                         "kind": "stand-in: the prevalence before the change measured",
                         "pad_share": PAD_SHARE, "pad_mean": PAD_MEAN, "sigmas": SIGMAS},
           "anomalies": {}, "hangs": {"anomaly": ""}, "errors": {"anomaly": ""}, "types": {}}
    for star, mine in summary.items():
        spec = {"n_ref": mine["n"], "measures": {}}
        for measure in (*MEAN_MEASURES, *SHARE_MEASURES):
            value = mine["values"][measure]
            low, high = measure_range(measure, value, mine["n"], n_run, mine["sd"].get(measure, 0.0))
            bounds = {"min": low, "max": high, "kind": "stand-in",
                      "source": f"pre-cut-over {value:.3f} over {mine['n']} systems at {commit[:9]}"}
            spec["measures"][measure] = bounds
        out["types"][star] = spec
    return out


def tag_current(expectations, records, item, note, only=None):
    """Marks every untagged range the current code breaks with `item` (those `only` accepts)."""
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


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description="Re-baselines tests/data/prevalence_systems_expectations.json (TEST.431).")
    sub = parser.add_subparsers(dest="command", required=True)
    measure = sub.add_parser("measure", help="draw the systems at the code being measured and store the records")
    measure.add_argument("--count", type=int, default=100)
    measure.add_argument("--out", required=True)
    build = sub.add_parser("build", help="turn a reference's records into ranges")
    build.add_argument("records")
    build.add_argument("--commit", required=True)
    build.add_argument("--label", default="")
    build.add_argument("--run-count", type=int, required=True, help="systems per star type a later run draws")
    build.add_argument("--out", default=str(EXPECTATIONS_PATH))
    build.add_argument("--tag-current", help="records of the code as it is now: ranges it breaks are tagged with --anomaly")
    build.add_argument("--anomaly")
    build.add_argument("--anomaly-note", default="")
    build.add_argument("--hang-anomaly", default="")
    build.add_argument("--error-anomaly", default="")
    args = parser.parse_args(argv)
    if args.command == "measure":
        with open(args.out, "w", encoding="utf-8") as handle:
            json.dump(draw_types(args.count), handle)
        return 0
    with open(args.records, encoding="utf-8") as handle:
        expectations = build_expectations(json.load(handle), args.commit, args.label, args.run_count)
    expectations["hangs"]["anomaly"] = args.hang_anomaly
    expectations["errors"]["anomaly"] = args.error_anomaly
    if args.tag_current:
        with open(args.tag_current, encoding="utf-8") as handle:
            print("tagged", tag_current(expectations, json.load(handle), args.anomaly, args.anomaly_note))
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(expectations, handle, indent=1, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
