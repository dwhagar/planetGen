#!/usr/bin/env python3
"""
Sensitivity of the class quantities to the stand-in values (GEN.1311, part 1 of the fudge table).

    python scripts/stand_in_sensitivity.py rows [--json OUT.json]                    list every registry row the harness can and cannot override
    python scripts/stand_in_sensitivity.py screen [--systems 5] [--seed 20261108] [--jobs 4] [--only IDS] [--limit N] [--out OUT.json] [--report OUT.md]

The fudge table (design-drafts/fudge-table-from-stand-in-confidence-2026-11-06.md) sets the width of a class edge from how far the quantity could move if the
stand-ins that make it were wrong by what the registry says. This is the instrument for that: it moves **one stand-in at a time** across its registered range and records
how far each quantity the class fit reads moves.

1. **Rows.** The registry (`docs/design/research/gap-filling-plan.md`, read the way `build_stand_in_values.py` reads it) has a `constant` cell that names a value as
   `module.NAME` and a `knobs` cell that gives its range. A row is *overridable* when its first `module.NAME` is a module constant that exists and is a plain number, and its
   knobs hold a range the harness can read (`A to B`, `a factor of N`). Every other row is *not overridable* and `rows` says why (no constant named, no such constant, not a
   number, no range). A row that names several constants moves the first only (the others are listed in `notes`).
2. **Override.** `Override` sets the module constant for the length of a `with` block and puts it back, and clears every `functools` cache of the package so a value
   computed from the old constant is not reused. A constant another module copied with `from x import NAME` does not move; the screening pass finds that as a row that moves
   nothing.
3. **Sample.** A fixed sample of systems, each from its own run seed (`seed + index`, so a change in the number of draws of one system shifts no other), read body by body
   through `class_fit._found`: the very quantities the class fit reads (gravity, mean temperature, pressure, gas shares, ocean cover, the regional, day and night
   temperatures, the wet bulb). A run is paired with the baseline body by body (system, planet, moon); a system whose number of bodies changed is counted and left out of
   the pairs.
4. **Screening pass.** `screen` runs the baseline, then for each overridable row the low end and the high end of its range, and records per quantity the root mean square of
   the movement against the baseline (kelvin for a temperature, the natural log of the ratio for a pressure, a gravity or a mass, the plain difference for a share or a
   fraction) and the share of bodies that moved. A row whose both ends move no quantity is listed as one that moves nothing (a row to drop from the full pass).

Everything it draws goes through `draw` (the one random wrapper); nothing here changes a generated world.
"""

import argparse
import collections
import contextlib
import functools
import importlib
import json
import math
import os
import re
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)

PLAN = os.path.join(HERE, "..", "docs", "design", "research", "gap-filling-plan.md")
DEFAULT_SEED = 20261108
DEFAULT_SYSTEMS = 5

QUANTITY_GROUPS = (
    ("gravity", ("gravity_g",)),
    ("mean temperature", ("temperature_k",)),
    ("pressure", ("pressure_pa",)),
    ("gas shares", ("share_", "o2_partial_pa", "co2_partial_pa", "reducing_share")),
    ("ocean and water", ("ocean_fraction", "water_mass_fraction", "ocean_ph", "ammonia_water_mass_fraction", "hydrosphere")),
    ("band temperatures", ("coldest_region_k", "hottest_region_k", "day_temperature_k", "night_temperature_k")),
    ("wet bulb and humidity", ("wet_bulb_c", "relative_humidity")),
    ("mass", ("mass_earth",)),
)
"""tuple: `(group, prefixes or names)` that sort the quantities of `class_fit._found` for the report; a quantity in none of them is `other`."""
LOG_QUANTITIES = ("gravity_g", "pressure_pa", "mass_earth", "o2_partial_pa", "co2_partial_pa", "flux_variation")
"""tuple: The quantities that move by orders of magnitude: their movement is the natural log of the ratio. Every other number moves by its plain difference."""
LOG_FLOOR = 1.0e-6
"""float: The value a log quantity is held above, so a zero is a number."""
SHARE_EPSILON = 1.0e-9
"""float: A movement under this is no movement (rounding)."""


# --- The rows of the registry -----------------------------------------------------------------------------------------------------------------------------------------

NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
CONSTANT = re.compile(r"\b([a-z][a-z0-9_]*)\.([A-Z][A-Z0-9_]{2,})\b")
CITATION = re.compile(r"\[[^\]]*\]")
RANGE = re.compile(rf"({NUMBER})\s*(?:[A-Za-z/%°]+\s*)?(?:to|-|–)\s*({NUMBER})")
FACTOR = re.compile(rf"factor of ({NUMBER})")


def read_rows(plan_path=PLAN):
    """
    Every registry row that has a `constant` cell and a `knobs` cell, as dicts (`fill_id`, `section`, `tier`, `tweak`, `constant`, `knobs`). The rows are the plan's, read by
    `build_stand_in_values.registry_tables`, so the harness and the stand-in list cannot disagree about what a row is.
    """
    import build_stand_in_values

    with open(plan_path, encoding="utf-8") as handle:
        tables = build_stand_in_values.registry_tables(handle.read())
    rows = []
    for section, (header, cells) in tables.items():
        if "constant" not in header or "knobs" not in header:
            continue
        for cell in cells:
            if len(cell) != len(header):
                continue
            row = dict(zip(header, cell))
            rows.append({"fill_id": row["fill_id"], "section": section, "tier": row.get("tier"), "tweak": row.get("tweak"), "constant": row["constant"], "knobs": row["knobs"]})
    return rows


def parse_range(knobs, value):
    """
    `(low, high)` of the registered range of a constant whose current value is `value`, or `None`. Reads the first of `A to B` (units after the number are allowed, citations
    are dropped), `a factor of N` (`value / N` to `value * N`); the text after the first `;` is another knob and is left out.
    """
    text = CITATION.sub("", knobs).split(";")[0]
    match = RANGE.search(text)
    if match:
        low, high = float(match.group(1)), float(match.group(2))
        if low == high:
            return None
        return (low, high) if low < high else (high, low)
    match = FACTOR.search(text)
    if match and value:
        factor = float(match.group(1))
        if factor > 1.0:
            return (value / factor, value * factor) if value > 0 else (value * factor, value / factor)
    return None


def module_index():
    """`{file name without .py: [dotted module names]}` of every module under `src/planetgen`, so `module.NAME` finds the module wherever in the package it lives."""
    root = os.path.join(HERE, "..", "src")
    index = collections.defaultdict(list)
    for folder, _, files in os.walk(os.path.join(root, "planetgen")):
        for file in files:
            if file.endswith(".py") and file != "__init__.py":
                dotted = os.path.relpath(os.path.join(folder, file[:-3]), root).replace(os.sep, ".")
                index[file[:-3]].append(dotted)
    return index


_INDEX = None


def resolve(module_name, name):
    """`(module, value)` of the module constant `name` of the module called `module_name` (anywhere under `planetgen`), or `(None, reason)`."""
    global _INDEX
    if _INDEX is None:
        _INDEX = module_index()
    candidates = _INDEX.get(module_name, [])
    if not candidates:
        return None, f"no module {module_name}"
    importlib.import_module("planetgen.generation.system")     # the normal import order first; some modules import each other in a ring
    problem = None
    for qualified in candidates:
        try:
            module = importlib.import_module(qualified)
        except Exception as error:
            problem = f"module {module_name} does not import ({type(error).__name__})"
            continue
        if hasattr(module, name):
            return module, getattr(module, name)
    return None, problem or f"no module constant {module_name}.{name}"


def classify_rows(rows):
    """
    Each row with its status: `target` (`module`, `name`, `value`, `low`, `high`) where it can be overridden, else `reason` (`no constant named`, `no such constant`, `not a
    plain number`, `no range to read`). `notes` holds the other constants a row names, which are not moved.
    """
    found = []
    for row in rows:
        names = CONSTANT.findall(CITATION.sub("", row["constant"]))
        entry = dict(row, notes=[f"{m}.{n}" for m, n in names[1:]])
        if not names:
            found.append(dict(entry, reason="no constant named"))
            continue
        module_name, name = names[0]
        module, value = resolve(module_name, name)
        if module is None:
            found.append(dict(entry, reason=value))
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            found.append(dict(entry, reason="not a plain number"))
            continue
        span = parse_range(row["knobs"], float(value))
        if span is None:
            found.append(dict(entry, reason="no range to read"))
            continue
        found.append(dict(entry, target={"module": module.__name__, "name": name, "value": value, "low": span[0], "high": span[1]}))
    return found


def clear_caches():
    """Clears every `functools` cache of the package, so a value worked out from the old constant is not reused."""
    for name, module in list(sys.modules.items()):
        if not name.startswith("planetgen") or module is None:
            continue
        for attribute in list(vars(module).values()):
            if type(attribute).__name__ == "_lru_cache_wrapper":        # a type check: a proxy (Flask's `current_app`) raises on any attribute read
                attribute.cache_clear()


def holders(module_name, name):
    """The modules of the package that hold the same-named attribute with the very value of `module_name.name` (a `from x import NAME` copy), the home module first."""
    home = importlib.import_module(module_name)
    original = getattr(home, name)
    found = [home]
    for key, module in list(sys.modules.items()):
        if module is None or module is home or not key.startswith("planetgen"):
            continue
        if name in vars(module) and vars(module)[name] is original:
            found.append(module)
    return found


@contextlib.contextmanager
def override(module_name, name, value):
    """
    Sets the module constant `module_name.name` to `value` (an int stays an int when the value is whole) for the length of the block, then puts it back. A module that copied
    the constant with `from x import NAME` (it holds the very same object) is set too, so the override reaches every reader.
    """
    modules = holders(module_name, name)
    original = getattr(modules[0], name)
    new = int(round(value)) if isinstance(original, int) and float(value).is_integer() else value
    for module in modules:
        setattr(module, name, new)
    clear_caches()
    try:
        yield len(modules)
    finally:
        for module in modules:
            setattr(module, name, original)
        clear_caches()


# --- The sample ------------------------------------------------------------------------------------------------------------------------------------------------------------

def group_of(quantity):
    for group, names in QUANTITY_GROUPS:
        if any(quantity == item or (item.endswith("_") and quantity.startswith(item)) for item in names):
            return group
    return "other"


def read_system(seed):
    """`{(planet, moon): {quantity: number}}` of the bodies of the system made from run seed `seed` (`moon` is -1 for the planet itself); the quantities of `class_fit._found`."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import class_fit
    from planetgen.util import draw

    config = SystemConfig()
    config.PLANETS = True
    with draw.bound(seed):
        system = StarSystem(system_config=config)
    bodies = {}
    for planet_index, planet in enumerate([p for p in system.all_planets if p.body_type != "a"]):
        for moon_index, body in [(-1, planet), *enumerate(planet.moons)]:
            try:
                found = class_fit._found(body)
            except Exception as error:      # a body the class fit cannot read is recorded as such and left out of the pairs
                bodies[(planet_index, moon_index)] = {"error": type(error).__name__}
                continue
            bodies[(planet_index, moon_index)] = {key: float(value) for key, value in found.items() if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)}
    return bodies


def read_sample(systems, seed):
    """The bodies of `systems` systems from the run seeds `seed`, `seed + 1`, ..."""
    return [read_system(seed + index) for index in range(systems)]


def movement(quantity, base, moved):
    """How far `quantity` moved from `base` to `moved`: the natural log of the ratio for a `LOG_QUANTITIES`, else the plain difference."""
    if quantity in LOG_QUANTITIES:
        return math.log(max(moved, LOG_FLOOR) / max(base, LOG_FLOOR))
    return moved - base


def compare(base_sample, sample):
    """
    Per quantity, over the bodies present in both samples: `{quantity: [sum of squared movements, bodies, bodies that moved]}`, and the number of systems whose set of bodies
    changed (`structure`).
    """
    totals = collections.defaultdict(lambda: [0.0, 0, 0])
    structure = 0
    for base_bodies, bodies in zip(base_sample, sample):
        if set(base_bodies) != set(bodies):
            structure += 1
            continue
        for key, base in base_bodies.items():
            moved = bodies[key]
            if "error" in base or "error" in moved:
                continue
            for quantity, value in base.items():
                if quantity not in moved:
                    continue
                step = movement(quantity, value, moved[quantity])
                entry = totals[quantity]
                entry[0] += step * step
                entry[1] += 1
                entry[2] += abs(step) > SHARE_EPSILON
    return {quantity: tuple(entry) for quantity, entry in totals.items()}, structure


def summarize(totals):
    """`{quantity: {"rms": ..., "moved": share of bodies}}` of a `compare` result."""
    return {quantity: {"rms": math.sqrt(entry[0] / entry[1]) if entry[1] else 0.0, "moved": entry[2] / entry[1] if entry[1] else 0.0} for quantity, entry in totals.items()}


# --- The pass ---------------------------------------------------------------------------------------------------------------------------------------------------------------

_BASELINE = None
_SETTINGS = None


def _evaluate(task):
    """Runs the sample with one constant at its low end and at its high end against the baseline (a worker of `screen`)."""
    fill_id, target = task
    systems, seed = _SETTINGS
    started = time.time()
    result = {"fill_id": fill_id}
    errors = []
    for end in ("low", "high"):
        try:
            with override(target["module"], target["name"], target[end]):
                totals, structure = compare(_BASELINE, read_sample(systems, seed))
            result[end] = {"quantities": summarize(totals), "structure": structure}
        except Exception as error:      # an end the generator cannot run (a zero divisor, a value out of its domain) is recorded and the other end still counts
            errors.append(f"{end} end {target[end]:.4g}: {type(error).__name__}: {error}")
            result[end] = {"quantities": {}, "structure": 0}
    if errors:
        result["errors"] = errors
        if len(errors) == 2:
            result["error"] = "; ".join(errors)
    result["seconds"] = round(time.time() - started, 1)
    return result


def screen(rows, systems, seed, jobs, progress=print, checkpoint=None):
    """The screening pass over the overridable `rows` (`classify_rows`): the baseline, then every row at both ends. Returns `{"baseline": ..., "results": [...]}`."""
    global _BASELINE, _SETTINGS
    _BASELINE, _SETTINGS = read_sample(systems, seed), (systems, seed)
    tasks = [(row["fill_id"], row["target"]) for row in rows if "target" in row]
    order = {fill_id: n for n, (fill_id, _) in enumerate(tasks)}
    results = []
    if checkpoint and os.path.exists(checkpoint):
        with open(checkpoint, encoding="utf-8") as handle:
            results = [json.loads(line) for line in handle if line.strip()]
        done = {result["fill_id"] for result in results}
        tasks = [task for task in tasks if task[0] not in done]
        progress(f"  resuming: {len(done)} rows already in {checkpoint}")
    progress(f"baseline: {sum(len(s) for s in _BASELINE)} bodies in {systems} systems; {len(tasks)} constants, two ends each, {jobs} jobs")
    with Pool(jobs) as pool:
        for count, result in enumerate(pool.imap_unordered(_evaluate, tasks, chunksize=1), 1):
            results.append(result)
            if checkpoint:
                with open(checkpoint, "a", encoding="utf-8") as handle:
                    handle.write(json.dumps(result) + "\n")
            if count % 20 == 0 or count == len(tasks):
                progress(f"  {count} of {len(tasks)} done")
    results.sort(key=lambda result: order[result["fill_id"]])
    return {"bodies": sum(len(s) for s in _BASELINE), "systems": systems, "seed": seed, "results": results}


def strongest(result):
    """`(movement, quantity)` of the largest root mean square movement of a screening result, over both ends and every quantity, with the quantity's unit scaled so a kelvin and a
    log ratio compare only within one quantity; used to rank rows within a group, not across."""
    best = (0.0, None)
    for end in ("low", "high"):
        for quantity, entry in result.get(end, {}).get("quantities", {}).items():
            if entry["rms"] > best[0]:
                best = (entry["rms"], quantity)
    return best


def moves_nothing(result):
    return "error" not in result and "errors" not in result and all(entry["rms"] <= SHARE_EPSILON for end in ("low", "high") for entry in result[end]["quantities"].values()) \
        and result["low"]["structure"] == 0 and result["high"]["structure"] == 0


def report(rows, screened):
    """The markdown of a screening pass: how many rows can and cannot be overridden and why, the ones that move nothing, the ones that fail, and the strongest movers by group."""
    by_id = {row["fill_id"]: row for row in rows}
    movers, nothing, failed = [], [], []
    for result in screened["results"]:
        if "error" in result:
            failed.append(result)
        elif moves_nothing(result):
            nothing.append(result)
        else:
            movers.append(result)
    reasons = collections.Counter(row["reason"] for row in rows if "reason" in row)
    lines = [f"# Stand-in sensitivity, screening pass ({screened['systems']} systems from seed {screened['seed']}, {screened['bodies']} bodies)", "",
             f"Registry rows with a constant and a knobs cell: {len(rows)}. Overridable: {sum(1 for row in rows if 'target' in row)}. Not overridable: {sum(reasons.values())}.", "",
             "## Not overridable", "", "| Reason | Rows |", "|---|---|"]
    lines += [f"| {reason} | {count} |" for reason, count in reasons.most_common()]
    lines += ["", f"## Overridable, run: {len(screened['results'])}", "",
              f"* Move at least one class quantity: {len(movers)}", f"* Move nothing at either end (drop from the full pass): {len(nothing)}",
              f"* Failed to run: {len(failed)}", ""]
    if nothing:
        lines += ["### Move nothing", "", ", ".join(result["fill_id"] for result in nothing), ""]
    if failed:
        lines += ["### Failed", ""] + [f"* {result['fill_id']}: {result['error']}" for result in failed] + [""]
    partial = [result for result in movers if "errors" in result]
    if partial:
        lines += ["### One end failed (the other end counts)", ""] + [f"* {result['fill_id']}: {'; '.join(result['errors'])}" for result in partial] + [""]
    changed = [result["fill_id"] for result in movers if result["low"]["structure"] or result["high"]["structure"]]
    if changed:
        lines += ["### Changed the number of bodies of a system (those systems are left out of the pairs)", "", ", ".join(changed), ""]
    lines += ["## Strongest movers by group", "", "(root mean square movement over the paired bodies, the larger of the two ends; kelvin for a temperature, natural log of the ratio for a log quantity, plain difference for a share)", ""]
    for group, _ in QUANTITY_GROUPS:
        ranked = []
        for result in movers:
            for end in ("low", "high"):
                for quantity, entry in result[end]["quantities"].items():
                    if group_of(quantity) == group and entry["rms"] > SHARE_EPSILON:
                        ranked.append((entry["rms"], result["fill_id"], quantity))
        best = {}
        for rms, fill_id, quantity in ranked:
            key = (fill_id, quantity)
            best[key] = max(best.get(key, 0.0), rms)
        top = sorted(((rms, key) for key, rms in best.items()), reverse=True)[:8]
        lines += [f"### {group}", ""]
        lines += [f"* {key[0]} on {key[1]}: {rms:.3g} ({by_id[key[0]]['constant'][:70].strip()})" for rms, key in top] or ["* none moved"]
        lines.append("")
    return "\n".join(lines)


def combine(screens):
    """
    The stand-in by end-result matrix from several screening passes (`--out` files, one per seed): `{fill_id: {"quantities": {quantity: {"rms", "moved", "bodies"}}, "structure": n, "errors": [...]}}`.
    A row's movement in a quantity is the root mean square over every pair of every pass, taken at the end that moves it more; `moved` is the share of paired bodies that moved; a row that
    appears twice in one pass (a resumed run) counts once, the later entry. `bodies` is how many pairs the row's passes held.
    """
    matrix = {}
    for screened in screens:
        by_row = {}
        for result in screened["results"]:
            by_row[result["fill_id"]] = result
        for fill_id, result in by_row.items():
            entry = matrix.setdefault(fill_id, {"quantities": {}, "structure": 0, "errors": [], "passes": 0})
            entry["passes"] += 1
            entry["structure"] += result["low"]["structure"] + result["high"]["structure"]
            entry["errors"] += result.get("errors", [])
            for quantity in set(result["low"]["quantities"]) | set(result["high"]["quantities"]):
                ends = [result[end]["quantities"].get(quantity) for end in ("low", "high")]
                ends = [end for end in ends if end]
                best = max(ends, key=lambda end: end["rms"])
                held = entry["quantities"].setdefault(quantity, {"sum_sq": 0.0, "moved": 0.0, "n": 0})
                held["sum_sq"] += best["rms"] ** 2
                held["moved"] += best["moved"]
                held["n"] += 1
    for entry in matrix.values():
        entry["quantities"] = {quantity: {"rms": math.sqrt(held["sum_sq"] / held["n"]), "moved": held["moved"] / held["n"], "passes": held["n"]}
                               for quantity, held in entry["quantities"].items()}
    return matrix


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    rows_parser = sub.add_parser("rows", help="list the rows and whether each can be overridden")
    rows_parser.add_argument("--json")
    combine_parser = sub.add_parser("combine", help="pool several screening passes (one per seed) into the stand-in by end-result matrix")
    combine_parser.add_argument("passes", nargs="+", help="the --out files of the passes")
    combine_parser.add_argument("--out", required=True, help="write the matrix as JSON")
    screen_parser = sub.add_parser("screen", help="the screening pass")
    screen_parser.add_argument("--systems", type=int, default=DEFAULT_SYSTEMS)
    screen_parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    screen_parser.add_argument("--jobs", type=int, default=4)
    screen_parser.add_argument("--only", help="comma-separated fill ids")
    screen_parser.add_argument("--movers-of", help="a screening JSON (--out of an earlier pass): run only the rows that moved a quantity or changed a system's bodies there")
    screen_parser.add_argument("--limit", type=int, help="run only the first N overridable rows (a trial)")
    screen_parser.add_argument("--out", help="write the results as JSON")
    screen_parser.add_argument("--checkpoint", help="a JSON-lines file each finished row is appended to; a rerun with the same file skips the rows already in it")
    screen_parser.add_argument("--report", help="write the markdown report")
    args = parser.parse_args(argv)

    if args.command == "combine":
        screens = []
        for path in args.passes:
            with open(path, encoding="utf-8") as handle:
                screens.append(json.load(handle)["screened"])
        matrix = combine(screens)
        with open(args.out, "w", encoding="utf-8") as handle:
            json.dump(matrix, handle, indent=1, sort_keys=True)
        print(f"{len(matrix)} rows from {len(screens)} passes")
        return 0
    rows = classify_rows(read_rows())
    if args.command == "rows":
        counts = collections.Counter("overridable" if "target" in row else row["reason"] for row in rows)
        for reason, count in counts.most_common():
            print(f"{count:5d}  {reason}")
        if args.json:
            with open(args.json, "w", encoding="utf-8") as handle:
                json.dump(rows, handle, indent=1)
        return 0
    chosen = rows
    if args.only:
        wanted = set(args.only.split(","))
        chosen = [row for row in rows if row["fill_id"] in wanted]
    if args.movers_of:
        with open(args.movers_of, encoding="utf-8") as handle:
            earlier = json.load(handle)["screened"]["results"]
        movers = {result["fill_id"] for result in earlier if not moves_nothing(result)}
        chosen = [row for row in chosen if row["fill_id"] in movers]
    overridable = [row for row in chosen if "target" in row]
    if args.limit:
        overridable = overridable[:args.limit]
    screened = screen(overridable, args.systems, args.seed, args.jobs, checkpoint=args.checkpoint)
    text = report(rows if not (args.only or args.limit or args.movers_of) else overridable, screened)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as handle:
            json.dump({"rows": rows, "screened": screened}, handle, indent=1)
    if args.report:
        with open(args.report, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
