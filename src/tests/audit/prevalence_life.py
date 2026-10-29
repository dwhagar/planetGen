# tests/audit/prevalence_life.py

"""
The life census (TEST.435, part G of the prevalence suite): how far the life of a generated world gets, against the targets Boss set for
the life class system (docs/design/habitability/life-classes.md 5.2).

`SYSTEMS` seeded systems are generated whole. Every terrestrial planet and moon with a life chemistry is a life world, and it is put on
the highest rung it reached: `origin` (abiogenesis), `producers` (photosynthesis), `complex` (complex cells), `multicellular` (macroscopic
life with no animals), `animals` (fauna above 0) or `sentient` (a technological civilization). The shares of life worlds are compared with
the bands of `tests/data/prevalence_life_expectations.json` (stand-ins, tier T9 and tweak level L3, set from Boss's rule and Core 1's timing
review), widened by four binomial sigma of the run's own count of life worlds, so a small run never tightens a band.

The class of each life world by the shadow classifier (`physics.life_classes.best`, GEN.886; the class is stored from GEN.889) is counted
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
SYSTEM_TIME_LIMIT_S = 60
SIGMAS = 4.0

RUNGS = ("origin", "producers", "complex", "multicellular", "animals", "sentient")
"""tuple: The highest rung a life world reached, lowest first."""

MEASURES = ("stops_before_animals", "reaches_animals", "reaches_sentient", "sentient_of_animal_worlds")
"""tuple: What the bands bind: the share of life worlds that stop at origin, producers or complex cells, that reach animals (sentient
worlds included), that reach a sentient civilization, and the share of the worlds that have animals that are also sentient."""


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


def draw_life(indices):
    """`{"systems", "worlds", "life", "rungs": {rung: n}, "classes": {class: n}, "abandoned", "errors", "error_count"}` for the systems
    of `indices` (seeded by their index)."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    rungs, classes, errors = Counter(), Counter(), []
    worlds = life = abandoned = error_count = 0

    def _alarm(*_):
        raise _Abandoned()

    previous = signal.signal(signal.SIGALRM, _alarm)
    try:
        for index in indices:
            draw.set_run_seed(census.seed_for("life-census", index))
            signal.alarm(SYSTEM_TIME_LIMIT_S)
            try:
                system = StarSystem(system_config=SystemConfig())
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
            "abandoned": abandoned, "errors": sorted(set(errors)), "error_count": error_count}


def run(systems=SYSTEMS, workers=None):
    """The census on worker processes (the time guard is an alarm), merged."""
    workers = min(workers or os.cpu_count() or 1, max(systems, 1))
    chunks = [list(range(start, systems, workers)) for start in range(workers)]
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        parts = pool.map(draw_life, chunks, chunksize=1)
    merged = {"systems": 0, "worlds": 0, "life": 0, "rungs": Counter(), "classes": Counter(), "abandoned": 0, "errors": [], "error_count": 0}
    for part in parts:
        for key in ("systems", "worlds", "life", "abandoned", "error_count"):
            merged[key] += part[key]
        merged["rungs"].update(part["rungs"])
        merged["classes"].update(part["classes"])
        merged["errors"] += part["errors"]
    merged["rungs"], merged["classes"], merged["errors"] = dict(merged["rungs"]), dict(merged["classes"]), sorted(set(merged["errors"]))
    return merged


def summarize(result):
    """The four measured shares of life worlds (see `MEASURES`), and the worlds they rest on."""
    life = result["life"]
    rungs = result["rungs"]
    animal_worlds = rungs.get("animals", 0) + rungs.get("sentient", 0)
    values = {
        "stops_before_animals": sum(rungs.get(name, 0) for name in ("origin", "producers", "complex")) / life if life else 0.0,
        "reaches_animals": animal_worlds / life if life else 0.0,
        "reaches_sentient": rungs.get("sentient", 0) / life if life else 0.0,
        "sentient_of_animal_worlds": rungs.get("sentient", 0) / animal_worlds if animal_worlds else 0.0,
    }
    counts = {"stops_before_animals": life, "reaches_animals": life, "reaches_sentient": life, "sentient_of_animal_worlds": animal_worlds}
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
