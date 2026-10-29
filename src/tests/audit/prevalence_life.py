# tests/audit/prevalence_life.py

"""
The life census (TEST.435, part G of the prevalence suite): how far the life of a generated world gets, against the targets Boss set for
the life class system (docs/design/habitability/life-classes.md 5.2).

`SYSTEMS` seeded systems are generated whole. Every terrestrial planet and moon with a life chemistry is a life world, and it is put on
the highest rung it reached: `origin` (abiogenesis), `producers` (photosynthesis), `complex` (complex cells), `multicellular` (macroscopic
life with no animals), `animals` (fauna above 0) or `sentient` (a technological civilization). The shares of life worlds are compared with
the bands of `tests/data/prevalence_life_expectations.json` (stand-ins, tier T9 and tweak level L3, set from Boss's rule and Core 1's timing
review), widened by four binomial sigma of the run's own count of life worlds, so a small run never tightens a band.

The animal and sentient shares are taken per oxygenic life world (GEN.907, Boss 6:04 AM EDT 10-29): Core 1's timing figures describe worlds
with oxygenic life, and the default galaxy holds almost none (chlorophyll chemistries are listed for G and K stars and need a sunlit sea),
so they are measured on `HOST_SYSTEMS` systems made around each of `OXYGENIC_HOSTS`. The share that stops before animals is still taken over
every life world of the default systems. The class of each life world by the shadow classifier (`physics.life_classes.best`, GEN.886; the class is stored from GEN.889) is counted
too. Per-class bands are an outcome to measure, not to force (the design's words), so the class shares are in the report and not checked.
"""

import json
import math
import multiprocessing
import os
import pathlib
import signal
import traceback
from collections import Counter

from tests.audit import prevalence as census

EXPECTATIONS_PATH = pathlib.Path(__file__).resolve().parent.parent / "data" / "prevalence_life_expectations.json"
REPORT_NAME = "life-census-latest.md"

SYSTEMS = 240
OXYGENIC_HOSTS = ("G2V", "K5V")
"""tuple: The host stars of the oxygenic population. The default galaxy gives 4 oxygenic worlds in 240 systems (chlorophyll chemistries are
listed for G and K stars only and need a sunlit sea), so a share per oxygenic world is measured on systems made around these stars."""
HOST_SYSTEMS = 300
"""int: Systems drawn around each of `OXYGENIC_HOSTS`."""
SYSTEM_TIME_LIMIT_S = 60
SIGMAS = 4.0

RUNGS = ("origin", "producers", "complex", "multicellular", "animals", "sentient")
"""tuple: The highest rung a life world reached, lowest first."""

MEASURES = ("stops_before_animals", "reaches_animals", "reaches_sentient", "sentient_of_animal_worlds")
"""tuple: What the bands bind: the share of life worlds that stop at origin, producers or complex cells, that reach animals (sentient
worlds included), that reach a sentient civilization, and the share of the worlds that have animals that are also sentient."""

OXYGENIC_MEASURES = ("reaches_animals", "reaches_sentient", "sentient_of_animal_worlds")
"""tuple: The measures taken per oxygenic life world (GEN.907: Core 1's timing figures describe worlds with oxygenic life). The others are
taken over every life world, because most life in a generated galaxy stopping before animals is the point of that band."""


class _Abandoned(BaseException):
    """Raised by the alarm when one system ran past `SYSTEM_TIME_LIMIT_S`."""


def rung(body):
    """The highest rung of a life world (a body with a life chemistry)."""
    if getattr(body, "life_sentient", None):
        return "sentient"
    if (getattr(body, "life_fauna", None) or 0) > 0:
        return "animals"
    stage = getattr(body, "life_stage", None)
    if (getattr(body, "life_flora", None) or 0) > 0 and stage not in ("abiogenesis", "photosynthesis", "complex_cells"):
        return "multicellular"
    return {"photosynthesis": "producers", "complex_cells": "complex", "multicellularity": "multicellular",
            "technological_civilization": "multicellular"}.get(stage, "origin")


def class_of(body):
    """The shadow life class of a body (`None` for no class at or above the line)."""
    from planetgen.physics import life_classes

    best = life_classes.best(life_classes.score_classes(body))
    return best.id if best else None


def draw_life(indices, star_type=None):
    """`{"systems", "worlds", "life", "rungs": {rung: n}, "classes": {class: n}, "oxygenic": {rung: n}, "abandoned", "errors",
    "error_count"}` for the systems of `indices` (seeded by their index; made around `star_type` when given, with seeds of their own)."""
    from planetgen.generation import evolution
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    rungs, classes, oxygenic, errors = Counter(), Counter(), Counter(), []
    worlds = life = abandoned = error_count = 0

    def _alarm(*_):
        raise _Abandoned()

    previous = signal.signal(signal.SIGALRM, _alarm)
    try:
        for index in indices:
            draw.set_run_seed(census.seed_for("life-census", index) if star_type is None else census.seed_for("life-census", star_type, index))
            signal.alarm(SYSTEM_TIME_LIMIT_S)
            try:
                config = SystemConfig()
                config.STAR_TYPE = star_type
                system = StarSystem(system_config=config)
                for planet in system.all_planets:
                    if planet.body_type == "a":
                        continue
                    for body in (planet, *planet.moons):
                        if body.body_type != "t":
                            continue
                        worlds += 1
                        if body.life_chemical:
                            life += 1
                            rungs[rung(body)] += 1
                            if evolution.releases_oxygen(body.life_chemical):
                                oxygenic[rung(body)] += 1
                            classes[str(class_of(body))] += 1
            except _Abandoned:
                abandoned += 1
            except Exception as error:
                frame = traceback.extract_tb(error.__traceback__)[-1]
                errors.append(f"{type(error).__name__}: {error} ({pathlib.Path(frame.filename).name}:{frame.lineno})")
                error_count += 1
            finally:
                signal.alarm(0)
    finally:
        signal.signal(signal.SIGALRM, previous)
    return {"systems": len(list(indices)), "worlds": worlds, "life": life, "rungs": dict(rungs), "classes": dict(classes),
            "oxygenic": dict(oxygenic), "abandoned": abandoned, "errors": sorted(set(errors)), "error_count": error_count}


def _job(job):
    return draw_life(*job)


def run(systems=SYSTEMS, workers=None, host_systems=HOST_SYSTEMS):
    """The census on worker processes (the time guard is an alarm), merged: `systems` default systems for the shares of every life
    world, and `host_systems` systems around each of `OXYGENIC_HOSTS` for the shares of the oxygenic ones (`host_oxygenic`)."""
    workers = min(workers or os.cpu_count() or 1, max(systems, host_systems, 1))
    jobs = [(list(range(start, systems, workers)), None) for start in range(workers)]
    jobs += [(list(range(start, host_systems, workers)), host) for host in OXYGENIC_HOSTS for start in range(workers)]
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        parts = pool.map(_job, jobs, chunksize=1)
    merged = {"systems": 0, "worlds": 0, "life": 0, "rungs": Counter(), "classes": Counter(), "oxygenic": Counter(),
              "host_systems": 0, "host_oxygenic": Counter(), "abandoned": 0, "errors": [], "error_count": 0}
    for (_, host), part in zip(jobs, parts):
        for key in ("abandoned", "error_count"):
            merged[key] += part[key]
        merged["errors"] += part["errors"]
        if host is None:
            for key in ("systems", "worlds", "life"):
                merged[key] += part[key]
            merged["rungs"].update(part["rungs"])
            merged["classes"].update(part["classes"])
            merged["oxygenic"].update(part["oxygenic"])
        else:
            merged["host_systems"] += part["systems"]
            merged["host_oxygenic"].update(part["oxygenic"])
    for key in ("rungs", "classes", "oxygenic", "host_oxygenic"):
        merged[key] = dict(merged[key])
    merged["errors"] = sorted(set(merged["errors"]))
    return merged


def summarize(result):
    """The four measured shares (see `MEASURES`), and the worlds they rest on: the first over every life world of the default
    systems, the others over the oxygenic life worlds of the systems made around `OXYGENIC_HOSTS` (GEN.907)."""
    life = result["life"]
    rungs = result["rungs"]
    oxygenic = result.get("host_oxygenic", {})
    oxygenic_worlds = sum(oxygenic.values())
    animal_worlds = oxygenic.get("animals", 0) + oxygenic.get("sentient", 0)
    values = {
        "stops_before_animals": sum(rungs.get(name, 0) for name in ("origin", "producers", "complex")) / life if life else 0.0,
        "reaches_animals": animal_worlds / oxygenic_worlds if oxygenic_worlds else 0.0,
        "reaches_sentient": oxygenic.get("sentient", 0) / oxygenic_worlds if oxygenic_worlds else 0.0,
        "sentient_of_animal_worlds": oxygenic.get("sentient", 0) / animal_worlds if animal_worlds else 0.0,
    }
    counts = {"stops_before_animals": life, "reaches_animals": oxygenic_worlds, "reaches_sentient": oxygenic_worlds,
              "sentient_of_animal_worlds": animal_worlds}
    return values, counts


def load_expectations(path=EXPECTATIONS_PATH):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def allowed(bounds, n):
    """The band a run resting on `n` worlds allows: the stored one widened by `SIGMAS` binomial sigma at each edge."""
    def slack(share):
        return SIGMAS * math.sqrt(max(share * (1.0 - share), 0.0) / max(n, 1))
    return max(0.0, bounds["min"] - slack(bounds["min"])), min(1.0, bounds["max"] + slack(bounds["max"]))


def compare(result, expectations):
    """Every flag: each share against its band, and any system that raised or ran past the time limit. A measure that rests on no
    world at all (no animal world for the sentient-per-animal share) is not judged."""
    findings = []
    if not result["life"]:
        return [census.Finding("missing", "life", "life worlds", 0.0, 1.0, 1.0, "no world had life")]
    values, counts = summarize(result)
    if result.get("host_systems") and not sum(result.get("host_oxygenic", {}).values()):
        findings.append(census.Finding("missing", "life", "oxygenic life worlds", 0.0, 1.0, 1.0,
                                       f"none in {result['host_systems']} systems around {', '.join(OXYGENIC_HOSTS)}"))
    for measure, bounds in expectations["measures"].items():
        if not counts[measure]:
            continue
        low, high = allowed(bounds, counts[measure])
        if values[measure] < low or values[measure] > high:
            kind = "missing" if values[measure] == 0.0 and low > 0 else "out-of-range"
            findings.append(census.Finding(kind, "life", measure, values[measure], low, high,
                                           f"{bounds['kind']}: {bounds['source']}", bounds.get("anomaly", "")))
    if result["error_count"]:
        findings.append(census.Finding("error", "life", "systems that raised", float(result["error_count"]), 0.0, 0.0,
                                       "; ".join(result["errors"])[:300], expectations.get("errors", {}).get("anomaly", "")))
    if result["abandoned"]:
        findings.append(census.Finding("hang", "life", "systems over the time limit", float(result["abandoned"]), 0.0, 0.0,
                                       f"a system ran past {SYSTEM_TIME_LIMIT_S} s", expectations.get("hangs", {}).get("anomaly", "")))
    return findings


def stale_anomalies(expectations, findings):
    broken = {finding.measure for finding in findings}
    return [f"{measure} ({bounds['anomaly']})" for measure, bounds in expectations["measures"].items()
            if bounds.get("anomaly") and measure not in broken]


def format_report(result, findings, header=""):
    values, counts = summarize(result)
    lines = ["# Life census (TEST.435), latest run", ""]
    if header:
        lines += [header, ""]
    unexplained = [finding for finding in findings if not finding.known]
    explained = [finding for finding in findings if finding.known]
    lines += [f"Flags: {len(unexplained)} new, {len(explained)} known anomalies.", ""]
    lines += [f"- {finding.line()}" for finding in unexplained]
    if explained:
        lines += ["", "Known anomalies:", ""]
        lines += [f"- {finding.line()}" for finding in explained]
    lines += ["", f"## {result['systems']} systems, {result['worlds']} terrestrial worlds, {result['life']} with life", ""]
    lines += [f"- {measure}: {values[measure]:.4f} (of {counts[measure]})" for measure in MEASURES]
    lines += ["", "Highest rung: " + ", ".join(f"{name} {result['rungs'].get(name, 0)}" for name in RUNGS)]
    oxygenic = result.get("oxygenic", {})
    lines += ["", f"Oxygenic life worlds of the default systems: {sum(oxygenic.values())} of {result['life']}"
              + (", highest rung: " + ", ".join(f"{name} {oxygenic.get(name, 0)}" for name in RUNGS) if oxygenic else "")]
    host = result.get("host_oxygenic", {})
    lines += ["", f"Oxygenic life worlds in {result.get('host_systems', 0)} systems around {', '.join(OXYGENIC_HOSTS)} (the population of "
              f"reaches_animals, reaches_sentient and sentient_of_animal_worlds): {sum(host.values())}"
              + (", highest rung: " + ", ".join(f"{name} {host.get(name, 0)}" for name in RUNGS) if host else "")]
    lines += ["", "Shadow life class of the life worlds (an outcome to measure, not checked): "
              + ", ".join(f"{name} {count}" for name, count in sorted(result["classes"].items(), key=lambda item: -item[1]))]
    return "\n".join(lines) + "\n"


def main(argv=None):
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="How far the life of generated worlds gets (TEST.435).")
    parser.add_argument("--systems", type=int, default=SYSTEMS)
    args = parser.parse_args(argv)
    result = run(args.systems)
    sys.stdout.write(format_report(result, compare(result, load_expectations())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
