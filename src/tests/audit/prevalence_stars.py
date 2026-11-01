# tests/audit/prevalence_stars.py

"""
The star census (TEST.432, part C of the prevalence suite): what the generator makes of stars against what surveys found.

Two draws, both large because a researched norm binds only when the sample is:

- The star type mix. Stars drawn from the stellar population model (`stellar_evolution.sample_living_star`) at the default
  galaxy's reference solar position, each population picked by its density there, counted into the spectral letters of
  main-sequence stars, white dwarfs and evolved stars (subgiants and giants). Compared with the local census: the RECONS
  10 parsec census (Henry et al., data of 2018-04-12: 378 stars, 21 of them white dwarfs) and Mamajek's number densities
  of stars of each type (read 2026-10-28), as shares of all stars.
- Multiplicity by primary type. Whole systems with no planets (`SystemConfig.PLANETS = False`; planets are the system
  census's, TEST.431) around M2V, K5V, G2V, F5V and A0V, counted into the stars each holds: the share with a companion
  star, and for a solar-type primary the share with three or more. Compared with [R1175], Raghavan et al.
  2010, and [R3108], whose figures `tuning.py` already quotes.

`tests/data/prevalence_stars_expectations.json` holds the norm ranges (a range, not a point, wherever two sources differ:
the widest range stands) and the item that explains each range the generator is known to break. The sampling slack of a
run (four binomial sigma at the run's own draw count) is added when the run is compared, so a norm is never tightened by
a small sample. A failure is a flag no known item explains; the known ones print in the latest-run report, and a tag
whose range holds again is listed as stale there.
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

EXPECTATIONS_PATH = pathlib.Path(__file__).resolve().parent.parent / "data" / "prevalence_stars_expectations.json"
REPORT_NAME = "star-census-latest.md"

MIX_DRAWS = 40000
"""int: Stars drawn for the type mix (the smallest norm worth binding, 0.04 percent of B stars, is then 16 draws)."""

SYSTEMS_PER_TYPE = 2500
"""int: Systems drawn around each primary type (a share near one half is then good to about 1 percent)."""

MULTIPLICITY_TYPES = ("M2V", "K5V", "G2V", "F5V", "A0V")
SYSTEM_TIME_LIMIT_S = 60
SIGMAS = 4.0

MIX_MEASURES = ("M", "K", "G", "F", "A", "B", "O", "wd", "evolved")
"""tuple: The kinds of star the mix counts: a spectral letter for a main-sequence star, `wd` for a white dwarf, `evolved` for a
subgiant or a giant (Yerkes IV, III, II, I, 0)."""

BROWN_DWARF_LETTERS = "LTY"
"""str: A member whose type begins with one of these is a brown dwarf, not a star (the brown dwarf census is part D)."""


class _Abandoned(BaseException):
    """Raised by the alarm when one system ran past `SYSTEM_TIME_LIMIT_S`."""


def solar_position_pc():
    """The default galaxy's reference solar position (the radius and angle the reference place of the profiles uses), in the
    plane."""
    from tests import test_satellites as satellites

    return (satellites.SOLAR_RADIUS * math.cos(satellites.SOLAR_ANGLE),
            satellites.SOLAR_RADIUS * math.sin(satellites.SOLAR_ANGLE), 0.0)


def solar_populations():
    """`{population: density}` at the reference solar position."""
    from planetgen.galaxy import density
    from tests import test_satellites as satellites

    return density.population_densities(solar_position_pc(), satellites.SHAPE)


def _kind(mass, age, state):
    from planetgen.physics.stellar_evolution import star_params

    yerkes = state["yerkes_class"]
    if yerkes == "VII":
        return "wd"
    if yerkes in ("IV", "III", "II", "I", "Ia", "Ib", "0"):
        return "evolved"
    return star_params(mass, age, state)["type"][0]


def draw_mix(count):
    """`{"n": draws, "counts": {kind: n}, "errors": [...]}` for `count` stars of the reference solar position."""
    from planetgen.generation import star_population
    from planetgen.physics.stellar_evolution import sample_living_star
    from planetgen.util import draw

    draw.set_run_seed(census.seed_for("star-mix"))
    populations = solar_populations()
    counts, errors = Counter(), []
    for _ in range(count):
        try:
            population = star_population.pick_population(populations, rng=draw)
            counts[_kind(*sample_living_star(population=population))] += 1
        except Exception as error:  # a crash is a finding, not the end of the census
            frame = traceback.extract_tb(error.__traceback__)[-1]
            errors.append(f"{type(error).__name__}: {error} ({pathlib.Path(frame.filename).name}:{frame.lineno})")
    return {"n": count, "counts": dict(counts), "errors": sorted(set(errors)), "error_count": len(errors)}


def _stars_in(system):
    return sum(1 for star in system.stars if not str(star.type).startswith(tuple(BROWN_DWARF_LETTERS)))


def draw_multiplicity(star_type, count):
    """`{"n", "stars": {count: systems}, "abandoned", "errors", "error_count"}` for `count` planetless systems around
    `star_type`."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    draw.set_run_seed(census.seed_for("multiplicity", star_type))
    stars, errors, abandoned, error_count = Counter(), [], 0, 0

    def _alarm(*_):
        raise _Abandoned()

    previous = signal.signal(signal.SIGALRM, _alarm)
    try:
        for _ in range(count):
            config = SystemConfig()
            config.STAR_TYPE = star_type
            config.PLANETS = False
            signal.alarm(SYSTEM_TIME_LIMIT_S)
            try:
                stars[_stars_in(StarSystem(system_config=config))] += 1
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
    return {"n": count, "stars": {str(key): value for key, value in sorted(stars.items())}, "abandoned": abandoned,
            "errors": sorted(set(errors)), "error_count": error_count}


def _job(job):
    kind, *arguments = job
    return (kind, arguments[0] if kind == "multiplicity" else "mix",
            draw_mix(*arguments) if kind == "mix" else draw_multiplicity(*arguments))


def run(mix_draws=MIX_DRAWS, systems=SYSTEMS_PER_TYPE, types=MULTIPLICITY_TYPES, workers=None):
    """Both draws on worker processes (the time guard is an alarm, which a test run's own timeout must not share):
    `{"mix": result, "multiplicity": {star type: result}}`."""
    jobs = [("mix", mix_draws)] + [("multiplicity", star, systems) for star in types]
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, len(jobs))) as pool:
        results = pool.map(_job, jobs, chunksize=1)
    return {"mix": next(result for kind, _name, result in results if kind == "mix"),
            "multiplicity": {name: result for kind, name, result in results if kind == "multiplicity"}}


def summarize(draws):
    """The measured shares: `{"mix": {"n", "values", ...}, "multiplicity": {star: {"n", "values", ...}}}`."""
    mix = draws["mix"]
    kept = mix["n"] - mix["error_count"]
    mix_values = {kind: mix["counts"].get(kind, 0) / kept if kept else 0.0 for kind in MIX_MEASURES}
    out = {"mix": {"n": kept, "values": mix_values, "errors": mix["errors"], "error_count": mix["error_count"],
                   "abandoned": 0},
           "multiplicity": {}}
    for star, result in draws["multiplicity"].items():
        kept = result["n"] - result["abandoned"] - result["error_count"]
        stars = {int(key): value for key, value in result["stars"].items()}
        values = {}
        if kept:
            values["companion"] = sum(v for k, v in stars.items() if k >= 2) / kept
            values["triple_or_more"] = sum(v for k, v in stars.items() if k >= 3) / kept
        out["multiplicity"][star] = {"n": kept, "values": values, "errors": result["errors"],
                                     "error_count": result["error_count"], "abandoned": result["abandoned"]}
    return out


def load_expectations(path=EXPECTATIONS_PATH):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def allowed(bounds, n):
    """The range a run of `n` draws allows: the norm range widened by `SIGMAS` binomial sigma at each edge."""
    def slack(share):
        return SIGMAS * math.sqrt(max(share * (1.0 - share), 0.0) / max(n, 1))
    low = max(0.0, bounds["norm_min"] - slack(bounds["norm_min"])) if bounds["norm_min"] > 0 else 0.0
    return low, min(1.0, bounds["norm_max"] + slack(bounds["norm_max"]))


def _groups(summary, expectations):
    yield "mix", summary["mix"], expectations["mix"]
    for star, measures in expectations["multiplicity"].items():
        yield star, summary["multiplicity"].get(star), measures


def compare(summary, expectations):
    """Every flag over the census: each measure of the mix and of each primary type against its range, any draw that raised
    and any system over the time limit."""
    findings = []
    for name, mine, measures in _groups(summary, expectations):
        if not mine or not mine["n"]:
            findings.append(census.Finding("missing", name, "draws", 0.0, 1.0, 1.0, "nothing drawn"))
            continue
        for measure, bounds in measures.items():
            low, high = allowed(bounds, mine["n"])
            value = mine["values"].get(measure, 0.0)
            if value < low or value > high:
                kind = "missing" if value == 0.0 and low > 0 else "out-of-range"
                findings.append(census.Finding(kind, name, measure, value, low, high,
                                               f"{bounds['kind']}: {bounds['source']}", bounds.get("anomaly", "")))
        if mine["error_count"]:
            findings.append(census.Finding("error", name, "draws that raised", float(mine["error_count"]), 0.0, 0.0,
                                           "; ".join(mine["errors"])[:300], expectations.get("errors", {}).get("anomaly", "")))
        if mine["abandoned"]:
            findings.append(census.Finding("hang", name, "systems over the time limit", float(mine["abandoned"]), 0.0, 0.0,
                                           f"a system ran past {SYSTEM_TIME_LIMIT_S} s",
                                           expectations.get("hangs", {}).get("anomaly", "")))
    return findings


def stale_anomalies(expectations, findings):
    broken = {(f.cell, f.measure) for f in findings}
    groups = {"mix": expectations["mix"], **expectations["multiplicity"]}
    return [f"{name} {measure} ({bounds['anomaly']})" for name, measures in groups.items()
            for measure, bounds in measures.items() if bounds.get("anomaly") and (name, measure) not in broken]


def format_report(summary, findings, header=""):
    lines = ["# Star census (TEST.432), latest run", ""]
    if header:
        lines += [header, ""]
    unexplained = [f for f in findings if not f.known]
    explained = [f for f in findings if f.known]
    lines += [f"Flags: {len(unexplained)} new, {len(explained)} known anomalies.", ""]
    lines += [f"- {finding.line()}" for finding in unexplained]
    if explained:
        lines += ["", "Known anomalies:", ""]
        lines += [f"- {finding.line()}" for finding in explained]
    mix = summary["mix"]
    lines += ["", "## Star type mix at the reference solar position", "",
              f"{mix['n']} stars: " + ", ".join(f"{name} {mix['values'].get(name, 0.0) * 100:.2f}%" for name in MIX_MEASURES), "",
              "## Multiplicity by primary type", ""]
    for star in MULTIPLICITY_TYPES:
        mine = summary["multiplicity"].get(star)
        if mine:
            figures = ", ".join(f"{name} {value:.3f}" for name, value in mine["values"].items())
            lines.append(f"- {star} (n={mine['n']}, abandoned {mine['abandoned']}, raised {mine['error_count']}): {figures}")
    return "\n".join(lines) + "\n"


def main(argv=None):
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Draws the star census and prints its report (TEST.432).")
    parser.add_argument("--stars", type=int, default=MIX_DRAWS)
    parser.add_argument("--systems", type=int, default=SYSTEMS_PER_TYPE)
    args = parser.parse_args(argv)
    summary = summarize(run(args.stars, args.systems))
    findings = compare(summary, load_expectations())
    sys.stdout.write(format_report(summary, findings))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
