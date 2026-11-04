# tests/audit/redraw_rates.py

"""
The redraw-rate census (TEST.476, GEN.1006): how often the generator throws a drawn item away, or moves or cuts it, to make it
fit. Seeded full systems of each star-count group (ordinary systems and systems forced to 1 to 6 stars with `STARS`, so the rare
triples and quadruples appear), each built with `draw.counting()` on, summed into a rate per kind per generated system
(`planetgen.generation.redraw_kinds.KINDS`), by group, by the number of stars the system came out with and by the class of its
placed star. The counters read no stream, so a counted system is the same system.

Two gates read the figures (`tests/test_redraw_rates.py`):

* a **ceiling** per kind and group, from `tests/data/redraw_rate_expectations.json`: twice the baseline rate measured on the code
  before the change (plus one event of the reference's size), failing only when the run's rate is above it by more than its own
  sampling error, so a change that suddenly redraws more is caught and a repeat of the reference is not;
* an **invariant** for the overlap kinds (`redraw_kinds.OVERLAP_KINDS`): their count is exactly zero once the overlap-free build
  (GEN.1009 to GEN.1015) has deleted the GEN.1004 quick fix. It is held as an expected failure until then.

`python -m tests.audit.redraw_rates measure --count 2500 --out records.json` draws the big sample (the research thread's baseline
was 2,500 forced systems per star count); `build` turns records into the expectations file; `report` prints the table.
"""

import argparse
import json
import math
import multiprocessing
import os
import pathlib
import signal
import sys

from tests.audit import prevalence as census

EXPECTATIONS_PATH = pathlib.Path(__file__).resolve().parent.parent / "data" / "redraw_rate_expectations.json"
REPORT_NAME = "redraw-rate-latest.md"
SYSTEM_TIME_LIMIT_S = 180
"""int: Seconds one system may take before the census abandons it and counts it."""

GROUPS = {"ordinary": {}, "stars-1": {"STARS": 1}, "stars-2": {"STARS": 2}, "stars-3": {"STARS": 3}, "stars-4": {"STARS": 4},
          "stars-5": {"STARS": 5}, "stars-6": {"STARS": 6},
          "world-requested": {"HABITABLE_WORLD": True, "ASTEROID_BELT": True}, "world-barred": {"HABITABLE_WORLD": False}}
"""dict: The groups of the census, by name, with the `SystemConfig` settings each forces: the star count (`STARS`, so the rare triples
and quadruples appear), a habitable world and a belt asked for (the request retries) and habitable worlds barred (the redraw of a body that
came out habitable). An empty setting leaves it to the tables."""

BUILT_BUCKETS = ("1", "2", "3", "4", "5+")
"""tuple: The star counts a system is reported under, by the number of stars it came out with."""

CEILING_FACTOR = 2.0
"""float: A ceiling is this many times the baseline rate (a rate that doubles trips it)."""

SIGMAS = 2.0
"""float: How many standard errors of the run's own count the rate must exceed the ceiling by to fail (a Poisson count's)."""


class _Abandoned(BaseException):
    """Raised by the alarm when one system ran past `SYSTEM_TIME_LIMIT_S`."""


def bucket_of(star_count):
    """The report bucket (`BUILT_BUCKETS`) of a system with this many stars."""
    return "5+" if star_count >= 5 else str(max(int(star_count), 1))


def placed_class(system):
    """The class of the placed star as its first spectral token (`G2V` gives `GV`, a giant's `K0III` gives `KIII`)."""
    token = (system.primary_star.type or "?").split()[0]
    letter = token[0] if token else "?"
    return letter + "".join(ch for ch in token if ch in "IVX")


def draw_group(group, count, first=0):
    """`count` systems of `group`, as records `{"group", "number", "stars", "class", "counts"}`; a system over the time limit is
    `{"group", "number", "abandoned": True}` and one that raised is `{"group", "number", "error": text}`. Seeded from the group and the system's number,
    so a second run draws the same systems."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    from tests.fuzz_support import deterministic_entropy

    def _alarm(*_):
        raise _Abandoned()

    records = []
    previous = signal.signal(signal.SIGALRM, _alarm)
    try:
        for number in range(first, first + count):
            seed = census.seed_for("redraw", group, number)
            config = SystemConfig()
            for setting, value in GROUPS[group].items():
                setattr(config, setting, value)
            signal.alarm(SYSTEM_TIME_LIMIT_S)
            try:
                with deterministic_entropy(seed):
                    draw.set_run_seed(seed)
                    with draw.counting() as tally:
                        system = StarSystem(system_config=config)
                records.append({"group": group, "number": number, "stars": len(system.stars), "class": placed_class(system),
                                "counts": dict(tally)})
            except _Abandoned:
                records.append({"group": group, "number": number, "abandoned": True})
            except Exception as error:  # a crash is a finding, not the end of the census
                frame = __import__("traceback").extract_tb(error.__traceback__)[-1]
                records.append({"group": group, "number": number, "error": f"{type(error).__name__}: {error} ({pathlib.Path(frame.filename).name}:{frame.lineno})"})
            finally:
                signal.alarm(0)
    finally:
        signal.signal(signal.SIGALRM, previous)
    return records


def _job(job):
    return draw_group(*job)


def draw_groups(count, groups=tuple(GROUPS), workers=None, chunk=10):
    """The records of every group, on worker processes, in chunks of `chunk` systems (the time guard is an alarm, which a test run's own
    timeout must not share)."""
    jobs = [(group, min(chunk, count - first), first) for group in groups for first in range(0, count, chunk)]
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, len(jobs))) as pool:
        results = pool.map(_job, jobs, chunksize=1)
    return [record for result in results for record in result]


def summarize(records):
    """Per group: `{"n": systems, "abandoned": n, "errors": [...], "events": {kind: count}, "systems": {kind: systems with at least one},
    "squares": {kind: sum of the squared per-system counts}}`;
    and the same per built star bucket (`by_stars`) and per placed-star class (`by_class`)."""
    def tally(rows):
        events, systems, squares = {}, {}, {}
        for row in rows:
            for kind, count in row["counts"].items():
                events[kind] = events.get(kind, 0) + count
                systems[kind] = systems.get(kind, 0) + 1
                squares[kind] = squares.get(kind, 0) + count * count
        return {"n": len(rows), "events": events, "systems": systems, "squares": squares}

    out = {"groups": {}, "by_stars": {}, "by_class": {}}
    for group in GROUPS:
        mine = [r for r in records if r["group"] == group]
        if not mine:
            continue
        kept = [r for r in mine if "counts" in r]
        out["groups"][group] = {**tally(kept), "abandoned": sum(1 for r in mine if r.get("abandoned")),
                                "errors": sorted({r["error"] for r in mine if r.get("error")}),
                                "error_count": sum(1 for r in mine if r.get("error"))}
    kept = [r for r in records if "counts" in r]
    for bucket in BUILT_BUCKETS:
        rows = [r for r in kept if bucket_of(r["stars"]) == bucket]
        if rows:
            out["by_stars"][bucket] = tally(rows)
    for letter in sorted({r["class"] for r in kept}):
        out["by_class"][letter] = tally([r for r in kept if r["class"] == letter])
    return out


def wilson(successes, n, z=1.96):
    """The Wilson score interval `(low, high)` of a share of `n` systems."""
    if n <= 0:
        return 0.0, 1.0
    p = successes / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, centre - half), min(1.0, centre + half)


def rate_lower_bound(events, n, squares=None):
    """The lower bound of an events-per-system rate: the mean less `SIGMAS` standard errors of the mean. The standard error is the spread
    of the per-system counts when their sum of squares is given (the events of a system come in clusters, so a Poisson count would be far
    too tight); without it, that of a Poisson count."""
    if n <= 0:
        return 0.0
    if squares is None or n < 2:
        return max(0.0, events - SIGMAS * math.sqrt(events)) / n
    mean = events / n
    variance = max(0.0, (squares - n * mean * mean) / (n - 1))
    return max(0.0, mean - SIGMAS * math.sqrt(variance / n))


def load_expectations(path=EXPECTATIONS_PATH):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def ceilings_broken(summary, expectations):
    """Every (group, kind) whose run rate is above its ceiling by more than the run's own sampling error:
    `[(group, kind, events, n, rate, ceiling)]`. The overlap kinds are the invariant's, not the ceilings'."""
    from planetgen.generation import redraw_kinds

    broken = []
    for group, spec in expectations.get("groups", {}).items():
        mine = summary["groups"].get(group)
        if not mine or not mine["n"]:
            continue
        for kind, bounds in spec.get("kinds", {}).items():
            if kind in redraw_kinds.OVERLAP_KINDS:
                continue
            events = mine["events"].get(kind, 0)
            if rate_lower_bound(events, mine["n"], mine["squares"].get(kind, 0) if "squares" in mine else None) > bounds["ceiling"]:
                broken.append((group, kind, events, mine["n"], events / mine["n"], bounds["ceiling"]))
    return broken


def unregistered_kinds(summary):
    """The kinds that fired in the run and are not in `redraw_kinds.KINDS`, or have no ceiling in the expectations for their group."""
    from planetgen.generation import redraw_kinds

    return sorted({kind for mine in summary["groups"].values() for kind in mine["events"]} - set(redraw_kinds.KINDS))


def overlap_events(summary):
    """The overlap-kind events of the run: `{kind: count}` over every group."""
    from planetgen.generation import redraw_kinds

    found = {}
    for mine in summary["groups"].values():
        for kind in redraw_kinds.OVERLAP_KINDS:
            if mine["events"].get(kind):
                found[kind] = found.get(kind, 0) + mine["events"][kind]
    return found


def format_report(summary, header=""):
    from planetgen.generation import redraw_kinds

    lines = ["# Redraw-rate census (TEST.476), latest run", ""]
    if header:
        lines += [header, ""]
    groups = list(summary["groups"])
    lines += ["## Events per generated system, by kind and group", "",
              "| kind | scope | overlap | " + " | ".join(f"{g} (n={summary['groups'][g]['n']})" for g in groups) + " |",
              "|---|---|---|" + "---|" * len(groups)]
    for kind, spec in redraw_kinds.KINDS.items():
        cells = []
        for group in groups:
            mine = summary["groups"][group]
            events, systems = mine["events"].get(kind, 0), mine["systems"].get(kind, 0)
            low, high = wilson(systems, mine["n"])
            cells.append(f"{events / mine['n']:.4f} ({events}; {low:.3f} to {high:.3f})" if mine["n"] else "-")
        lines.append(f"| {kind} | {spec.scope} | {'yes' if spec.overlap else ''} | " + " | ".join(cells) + " |")
    lines += ["", "Each cell: events per system (events; the 95% interval of the share of systems with at least one).", ""]
    for title, key in (("by number of stars the system came out with", "by_stars"), ("by class of the placed star", "by_class")):
        parts = summary[key]
        lines += [f"## Events per generated system, {title}", "",
                  "| kind | " + " | ".join(f"{name} (n={parts[name]['n']})" for name in parts) + " |", "|---|" + "---|" * len(parts)]
        for kind in redraw_kinds.KINDS:
            if not any(parts[name]["events"].get(kind) for name in parts):
                continue
            lines.append(f"| {kind} | " + " | ".join(f"{parts[name]['events'].get(kind, 0) / parts[name]['n']:.4f}" for name in parts) + " |")
        lines.append("")
    problems = [(g, m["abandoned"], m["error_count"]) for g, m in summary["groups"].items() if m["abandoned"] or m["error_count"]]
    if problems:
        lines += ["## Systems that hung or raised", ""] + [f"- {g}: {a} over the time limit, {e} raised" for g, a, e in problems] + [""]
    return "\n".join(lines)


def build_expectations(records, commit, label):
    """The expectations file from a reference run: per group and kind the baseline rate and its ceiling (twice the baseline plus one
    event of the reference's size), all registered as stand-ins of the tier and tweak level below."""
    summary = summarize(records)
    out = {"reference": {"commit": commit, "label": label, "seeds": "crc32 of 'redraw:<group>:<number>' (tests.audit.prevalence.seed_for)",
                         "kind": "stand-in: the redraw rate before the overlap-free build measured", "ceiling_factor": CEILING_FACTOR,
                         "tier": "T9", "tweak": "L3",
                         "note": "T9 (from needed output), L3 (hand-set against the baseline run): the ceiling is the rate that would double, "
                                 "not a physical limit; GEN.1006's parts lower the baselines, and each part re-measures and tightens them"},
           "groups": {}}
    for group, mine in summary["groups"].items():
        kinds = {}
        for kind, events in sorted(mine["events"].items()):
            rate = events / mine["n"]
            kinds[kind] = {"baseline": round(rate, 6), "events": events,
                           "ceiling": round(CEILING_FACTOR * rate + 1.0 / mine["n"], 6)}
        out["groups"][group] = {"n_ref": mine["n"], "kinds": kinds}
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    measure = sub.add_parser("measure", help="draw the systems and write their records")
    measure.add_argument("--count", type=int, default=100, help="systems per group")
    measure.add_argument("--out", required=True)
    measure.add_argument("--groups", default=",".join(GROUPS))
    build = sub.add_parser("build", help="write the expectations file from records")
    build.add_argument("records")
    build.add_argument("--commit", required=True)
    build.add_argument("--label", default="")
    build.add_argument("--out", default=str(EXPECTATIONS_PATH))
    report = sub.add_parser("report", help="print the table of a records file")
    report.add_argument("records")
    args = parser.parse_args(argv)
    if args.command == "measure":
        records = draw_groups(args.count, tuple(g for g in args.groups.split(",") if g))
        with open(args.out, "w", encoding="utf-8") as handle:
            json.dump(records, handle)
        print(f"{len(records)} records written to {args.out}")
    elif args.command == "build":
        with open(args.records, encoding="utf-8") as handle:
            records = json.load(handle)
        expectations = build_expectations(records, args.commit, args.label)
        with open(args.out, "w", encoding="utf-8") as handle:
            json.dump(expectations, handle, indent=1, sort_keys=True)
            handle.write("\n")
        print(f"expectations for {len(expectations['groups'])} groups written to {args.out}")
    else:
        with open(args.records, encoding="utf-8") as handle:
            print(format_report(summarize(json.load(handle))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
