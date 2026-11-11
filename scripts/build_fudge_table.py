#!/usr/bin/env python3
"""
Builds the fudge table (GEN.1313, part 3 of the fudge table) from the stand-in by end-result matrix and the registry.

    python scripts/build_fudge_table.py            # write src/planetgen/physics/data/fudge_table.json
    python scripts/build_fudge_table.py --check    # build in memory and fail if the committed file differs
    python scripts/build_fudge_table.py --update   # what update.sh runs (GEN.1317): build the table with its header into PLANETGEN_DATA_DIR unless it is current
    python scripts/build_fudge_table.py --status   # say whether the table in PLANETGEN_DATA_DIR is current (exit 0), stale or missing (exit 1)

The width of a class quantity is how far it could move if the stand-ins that make it were wrong by what the registry says they could be (design-drafts/fudge-table-from-stand-in-confidence-2026-11-06.md,
section 3):

1. **The matrix** (`fudge_matrix.json`, made by `scripts/stand_in_sensitivity.py screen --frozen` and `combine`, kept in the repo so this script is quick and exact): per registry row and per quantity the root mean
   square movement of the quantity over the paired bodies when the row's constant goes to either end of its range (the end that moves it more), pooled over the passes.
2. **Combine.** Each cell is first clipped: a body that moves farther than half the gap between two class edges of the quantity has changed class whatever the distance (a runaway greenhouse, a
   tectonic switch), so no movement counts for more than that half-gap (`clip`; a quantity with one edge clips at `UNCAPPED_CLIP_SHOULDERS` hand shoulders). `sigma` of a quantity is the root sum of squares of its column, each cell multiplied by the factor of the row's tier (`TIER_FACTORS`: a range of tier T1 to T3 is used as printed, T4 to T8 widened by
   1.5, T9 by 2). The shoulder width is `SIGMA_MULTIPLE` (2) sigma: a value one sigma past a class edge keeps half its membership.
3. **Cap.** No width may exceed half the gap between a class edge and the nearest other edge of the same quantity in the class table, or the limit the quantity cannot pass (0, and 1 for a share); where the cap binds the table says so,
   because that is a place a class is missing, not a place for a fudge.
4. **Quantities with no measurement.** The frozen sample makes the end stages of rocky planets and moons again (air, water, surface, class, spin), so a quantity fixed before those stages (the mass and
   the gravity) has no column, and a quantity no registered constant moves has none either. Such a quantity takes the shoulder it has today (the hand-chosen CF/DM row of `planet_class_fit.json`) scaled by
   `1 + score / 12`, where `score` is the mean registry unreliability score (tier points plus tweak points, 0 best, 12 worst) of the area that makes it (`FALLBACK_SCORES`, from reports/stand-in-reliability-ranking-2026-11-06.md).
   This is the stand-in of registry row FT1; the table marks the quantity `fallback` so the next measured run replaces it.
5. **Rows not in the table.** A row whose perturbation changed which bodies exist in the old regenerate-the-system run but moved nothing when the bodies were held (`population_rows`) acts on the draw of the
   population, not on the quantity of a given body; it is listed, not summed (GEN.1318).

The table is derived: it is built from the matrix, the class table and the registry only, so two hosts get the same file, and `src/tests/test_fudge_table.py` fails when the committed file differs.
"""

import argparse
import collections
import json
import math
import os
import re
import statistics
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)

DATA = os.path.join(HERE, "..", "src", "planetgen", "physics", "data")
MATRIX = os.path.join(DATA, "fudge_matrix.json")
OUT = os.path.join(DATA, "fudge_table.json")

TIER_FACTORS = ((3, 1.0), (8, 1.5), (9, 2.0))
"""tuple: `(highest tier number, factor)`: the factor a row's range is widened by (registry row FT2)."""
SIGMA_MULTIPLE = 2.0
"""float: The shoulder is this many sigma wide: a value one sigma past the edge keeps half its membership (registry row FT3)."""
UNCAPPED_CLIP_LN = math.log(10.0)
"""float: Where a log quantity has a single class edge, a movement is clipped at a factor of ten (registry row FT4)."""
THERMAL = ("temperature_k", "hottest_region_k", "coldest_region_k", "day_temperature_k", "night_temperature_k")
"""tuple: The quantities that move together with the mean temperature; one of them with a single edge takes the median half-gap of the others as its clip (registry row FT4)."""
TOP_CONTRIBUTORS = 5
"""int: How many rows are named as the makers of a width."""
FALLBACK_SCORES = {
    "mass_earth": 5.2, "gravity_g": 3.3, "temperature_k": 3.9, "pressure_pa": 3.3, "ocean_fraction": 3.6, "water_mass_fraction": 3.6, "ammonia_water_mass_fraction": 3.6,
    "o2_partial_pa": 3.5, "co2_partial_pa": 3.5, "share_co2": 3.5, "share_h2": 3.5, "share_so2": 3.5, "share_ch4": 3.5, "reducing_share": 3.5,
    "hottest_region_k": 5.1, "coldest_region_k": 5.1, "day_temperature_k": 5.1, "night_temperature_k": 5.1, "wet_bulb_c": 5.1, "relative_humidity": 5.1,
}
"""dict: The mean unreliability score of the area that makes each quantity (the ranking report of 11-06, section 2b): the stand-in of registry row FT1."""
FALLBACK_BASE = {"mass_earth": 0.510826, "gravity_g": 0.105361}
"""dict: The hand shoulders the class fit had before the table (a mass, 40 percent: the log of 1 / 0.6; a gravity, 10 percent: the log of 1 / 0.9), in the width's unit. The frozen
sample cannot move the mass or the gravity (both are fixed before the end stages it makes again), so they keep these, scaled by `1 + score / 12` (registry row FT1); a measured run that moves the
mass law and the radius law with the population held (GEN.1318) replaces them."""
TIER = "T6"
TWEAK = "L0"
"""str: The tier and tweak of the widths themselves: propagated registry ranges with a stated factor, the physics argued (registry row FT3)."""
LOG_QUANTITIES = ("gravity_g", "pressure_pa", "mass_earth", "o2_partial_pa", "co2_partial_pa", "flux_variation")
"""tuple: The quantities that move by orders of magnitude: their movement and width are natural logs of a ratio (the same list as `stand_in_sensitivity.LOG_QUANTITIES`)."""
EPSILON = 1.0e-9
SHARES = ("ocean_fraction", "water_mass_fraction", "ammonia_water_mass_fraction", "share_co2", "share_h2", "share_so2", "share_ch4", "reducing_share", "relative_humidity")
"""tuple: The quantities that are shares or fractions: they stop at 0 and at 1, and the nearer stop acts as an edge when a cap is worked out (registry row FT4)."""


def tier_factor(tier):
    number = int(re.search(r"\d+", tier or "T9").group())
    for highest, factor in TIER_FACTORS:
        if number <= highest:
            return factor
    return TIER_FACTORS[-1][1]


def tier_number(tier):
    return int(re.search(r"\d+", tier or "T9").group())


def tuned(tweak):
    return (tweak or "") in ("L3", "L4")


def edges_of(quantity, classes):
    """The distinct class edges of `quantity` (`min` and `lt` conditions of the class table), sorted."""
    return sorted({condition.get("min", condition.get("lt")) for cls in classes for condition in cls.when if condition["quantity"] == quantity and ("min" in condition or "lt" in condition)})


def cap_of(quantity, edge, edges):
    """
    Half the gap from `edge` to the nearest other edge of the same quantity in the unit of the width (the log for a log quantity), or `None` when there is none. The limit a quantity cannot
    pass counts as an edge: a share stops at 0 and 1, so a shoulder never reaches past half the way to the limit. A temperature, a pressure and a mass have no limit that matters to a class.
    """
    others = [other for other in edges if other != edge]
    if quantity in SHARES:
        others = [other for other in others + [0.0, 1.0] if other != edge]
    if not others:
        return None
    if quantity in LOG_QUANTITIES:
        return min(abs(math.log(max(edge, 1.0e-6) / max(other, 1.0e-6))) for other in others) / 2.0
    return min(abs(edge - other) for other in others) / 2.0


def clipped_rms(cell, clip):
    """The root mean square movement of a matrix `cell` over its bodies with every movement held at `clip` (none when `clip` is `None` or the cell has no per-body movements): a body that moves
    farther than half the gap between two class edges has changed class whatever the distance, so the size past that carries no information about a shoulder."""
    steps = cell.get("steps")
    if clip is None or not steps or not cell.get("bodies"):
        return cell["rms"]
    return math.sqrt(sum(min(step, clip) ** 2 for step in steps) / cell["bodies"])


def fallback_base(quantity):
    """The width a quantity takes when no row of the matrix moves it, before the score scales it: the hand shoulder the class fit had before the table (`FALLBACK_BASE`), or `None`."""
    return FALLBACK_BASE.get(quantity)


def column(quantity, matrix, by_id, clip):
    """The rows of the matrix that move `quantity`, each with its tier factor and its movement clipped at `clip` (`clipped_rms`), and the sigma they sum to (root sum of squares of the widened movements)."""
    contributions = []
    for fill_id, entry in sorted(matrix["rows"].items()):
        cell = entry["quantities"].get(quantity)
        if not cell or cell["rms"] <= EPSILON:
            continue
        row = by_id.get(fill_id, {})
        factor = tier_factor(row.get("tier"))
        movement = clipped_rms(cell, clip)
        if movement > EPSILON:
            contributions.append({"fill_id": fill_id, "tier": row.get("tier"), "tweak": row.get("tweak"), "factor": factor, "movement": movement, "variance": (factor * movement) ** 2})
    return contributions, math.sqrt(sum(item["variance"] for item in contributions))


def class_conditions(path=None):
    """The class table (`planet_class_fit.json`) as objects with a `when` list, read straight from the file so the build does not need a table to exist."""
    import types

    with open(path or os.path.join(DATA, "planet_class_fit.json"), encoding="utf-8") as handle:
        document = json.load(handle)
    return [types.SimpleNamespace(letter=row["letter"], when=row["when"]) for row in document["classes"]]


def build(matrix, rows, classes):
    """The table as a dict, from the matrix document, the registry `rows` (`stand_in_sensitivity.read_rows`) and the class table (`class_conditions`)."""
    by_id = {row["fill_id"]: row for row in rows}
    quantities = {}
    names = sorted({condition["quantity"] for cls in classes for condition in cls.when if "equals" not in condition})
    thermal_caps = [cap for name in names if name in THERMAL for cap in (cap_of(name, edge, edges_of(name, classes)) for edge in edges_of(name, classes)) if cap is not None]
    thermal_clip = statistics.median(thermal_caps) if thermal_caps else None
    for quantity in names:
        edges = edges_of(quantity, classes)
        hand = fallback_base(quantity)
        uncapped = UNCAPPED_CLIP_LN if quantity in LOG_QUANTITIES else (thermal_clip if quantity in THERMAL else None)
        caps = [cap_of(quantity, edge, edges) for edge in edges]
        caps = [uncapped if cap is None else cap for cap in caps]        # an edge with no neighbour is capped at the clip of a lone edge (a factor of ten, or the median half-gap of the thermal values)
        clip = statistics.median([cap for cap in caps if cap is not None]) if any(cap is not None for cap in caps) else uncapped        # the quantity's own figure clips at the middle half-gap; each edge's width below at its own
        contributions, sigma = column(quantity, matrix, by_id, clip)
        total = sigma * sigma
        unit = "ln" if quantity in LOG_QUANTITIES else "plain"
        measured = bool(contributions)
        if not measured and hand is None:
            raise SystemExit(f"{quantity}: no row of the matrix moves it and it has no fallback base (FALLBACK_BASE)")
        fallback = hand * (1.0 + FALLBACK_SCORES[quantity] / 12.0) if hand else None
        width, basis = (SIGMA_MULTIPLE * sigma, "measured") if measured else (fallback, "fallback")
        top = sorted(contributions, key=lambda item: (-item["variance"], item["fill_id"]))[:TOP_CONTRIBUTORS]
        edge_rows = []
        for edge, cap in zip(edges, caps):
            # the width at one edge: the sigma of the movements clipped at that edge's own half-gap, twice, and never past the half-gap itself
            edge_sigma = column(quantity, matrix, by_id, cap if cap is not None else uncapped)[1] if measured else None
            edge_width = SIGMA_MULTIPLE * edge_sigma if measured else fallback
            edge_rows.append({"edge": edge, "cap": round(cap, 6) if cap is not None else None, "binds": cap is not None and edge_width > cap,
                              "width": round(edge_width if cap is None else min(edge_width, cap), 6)})
        quantities[quantity] = {
            "unit": unit, "basis": basis, "rows": len(contributions), "clip": round(clip, 6) if clip is not None else None, "sigma": round(sigma, 6), "width": round(width, 6),
            "fallback_base": round(hand, 6) if hand else None, "tier": TIER, "tweak": TWEAK,
            "tuned_share": round(sum(item["variance"] for item in contributions if tuned(item["tweak"])) / total, 4) if total else 0.0,
            "contributors": [{"fill_id": item["fill_id"], "tier": item["tier"], "tweak": item["tweak"], "movement": round(item["movement"], 6), "share": round(item["variance"] / total, 4)} for item in top],
            "edges": edge_rows,
        }
    return {
        "title": "The fudge table: the width of a class quantity from the confidence of the stand-ins that make it (GEN.1313)",
        "version": 1,
        "status": "Derived by scripts/build_fudge_table.py from fudge_matrix.json, the class table and the stand-in registry. Do not edit; rebuild.",
        "basis": ("Width = 2 sigma, sigma the root sum of squares of the movements of the quantity when each registered stand-in is moved across its range (tier T1 to T3 range as printed, T4 to T8 widened by 1.5, T9 by 2), "
                  "capped at half the gap to the next class edge. Quantities with no measured row take the hand shoulder scaled by 1 + score / 12 (registry row FT1). Correlated stand-ins are summed as independent, "
                  "which slightly understates the width."),
        "sample": matrix["sample"],
        "tier": TIER,
        "tweak": TWEAK,
        "factors": {"tier_factors": [list(pair) for pair in TIER_FACTORS], "sigma_multiple": SIGMA_MULTIPLE},
        "population_rows": sorted(matrix.get("population_rows", [])),
        "quantities": quantities,
    }


def assemble(passes, regenerated=None):
    """
    The matrix document from the `--out` files of the frozen screening passes (`stand_in_sensitivity.py screen --frozen`, one per seed; the first the full pass, the others of its movers) and,
    for the population rows, the matrix of the old regenerate-the-system runs: `rows` the paired rows, with the shoulder quantities only; `unstable_rows` those that never came back to the baseline after the override (left out), `partly_unstable_rows` those that did not in some pass (that pass left out);
    `population_rows` those the regenerate runs moved and the frozen runs did not.
    """
    import stand_in_sensitivity
    from planetgen.physics import class_fit

    screens = []
    for path in passes:
        with open(path, encoding="utf-8") as handle:
            screens.append(json.load(handle)["screened"])
    flagged = collections.defaultdict(int)
    for screened in screens:
        kept = []
        for result in screened["results"]:
            if any("did not return" in error for error in result.get("errors", [])):
                flagged[result["fill_id"]] += 1        # the bodies did not come back to the baseline after this row: its numbers in this pass are not from one state
            else:
                kept.append(result)
        screened["results"] = kept
    combined = stand_in_sensitivity.combine(screens)
    unstable = sorted(fill_id for fill_id in flagged if fill_id not in combined)
    partly = sorted(fill_id for fill_id in flagged if fill_id in combined)
    rows = {}
    for fill_id, entry in sorted(combined.items()):
        cells = {quantity: {"rms": round(cell["rms"], 6), "moved": round(cell["moved"], 4), "bodies": cell["bodies"], "steps": cell["steps"]} for quantity, cell in sorted(entry["quantities"].items())
                 if quantity in class_fit.SHOULDERS and cell["rms"] > EPSILON}
        if cells:
            rows[fill_id] = {"passes": entry["passes"], "quantities": cells}
    population = []
    if regenerated:
        with open(regenerated, encoding="utf-8") as handle:
            population = sorted(set(json.load(handle)) - set(rows) - set(unstable))
    from planetgen.util import generator_epoch

    return {"epoch": generator_epoch.generator_epoch(),      # the epoch of the checkout the passes were assembled on (GEN.1317: the table's header says how old the matrix is)
            "sample": {"mode": "frozen end stages (planets.remake_end_stages) of the rocky planets and moons of the systems below",
                       "passes": [{"seed": screened["seed"], "systems": screened["systems"], "bodies": screened["bodies"], "rows_run": len(screened["results"])} for screened in screens]},
            "rows": rows, "unstable_rows": unstable, "partly_unstable_rows": partly, "population_rows": population}


def render(table):
    return json.dumps(table, indent=1, sort_keys=False) + "\n"


def make():
    import stand_in_sensitivity

    with open(MATRIX, encoding="utf-8") as handle:
        matrix = json.load(handle)
    return render(build(matrix, stand_in_sensitivity.read_rows(), class_conditions()))


def update(force=False, say=print):
    """
    Update's step (GEN.1317): the table in `PLANETGEN_DATA_DIR` is left alone when its header is current (a hash of three small files, under a second, update stays a no-op on a
    current server) and built again from the stored matrix when it is missing or stale: the fast mode, a fraction of a second. Returns 0, or 1 when the build failed.
    """
    from planetgen.physics import fudge_table

    found, reason = fudge_table.state()
    if found == "current" and not force:
        say(f"fudge table: current (epoch {fudge_table.header_of(fudge_table.data_file())['epoch']}), nothing to do")
        return 0
    say(f"fudge table: {found} ({reason}); building it from the stored matrix")
    try:
        started = time.time()
        table = json.loads(make())
        with open(MATRIX, encoding="utf-8") as handle:
            measured = json.load(handle).get("epoch")
        table["header"] = {**fudge_table.current_inputs(), "matrix_epoch": measured}
        path = fudge_table.data_file()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        temporary = path + ".new"
        with open(temporary, "w", encoding="utf-8") as handle:
            handle.write(render(table))
        os.replace(temporary, path)
        os.chmod(path, 0o644)
    except (OSError, ValueError, KeyError) as problem:
        say(f"fudge table: the build failed ({problem}); bulk generation refuses to start until it is built")
        return 1
    say(f"fudge table: built into {path} in {time.time() - started:.2f} s")
    age = fudge_table.matrix_age_epochs(table["header"])
    if age is None:
        say("fudge table: the stored sensitivity matrix does not say when it was measured; refresh it with scripts/stand_in_sensitivity.py when the stand-ins have moved (about 1.5 hours on four cores)")
    elif age > fudge_table.MATRIX_MAX_AGE_EPOCHS:
        say(f"fudge table: the stored sensitivity matrix was measured {age} epochs ago; refresh it (scripts/stand_in_sensitivity.py screen --frozen, about 1.5 hours on four cores)")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--check", action="store_true", help="fail if the committed table differs from the matrix")
    parser.add_argument("--update", action="store_true", help="build the table with its header into PLANETGEN_DATA_DIR unless it is current (what update.sh runs)")
    parser.add_argument("--force", action="store_true", help="with --update: build even when the table is current")
    parser.add_argument("--status", action="store_true", help="exit 0 when the table in PLANETGEN_DATA_DIR is current, 1 when stale or missing")
    parser.add_argument("--assemble", nargs="+", metavar="PASS.json", help="write fudge_matrix.json from the frozen screening passes (the first the full pass), then build the table")
    parser.add_argument("--regenerated", help="with --assemble: the matrix JSON of the old regenerate-the-system runs, for the population rows")
    args = parser.parse_args(argv)
    if args.status:
        from planetgen.physics import fudge_table

        found, reason = fudge_table.state()
        print(f"fudge table: {found}" + (f" ({reason})" if reason else ""))
        return 0 if found == "current" else 1
    if args.update:
        return update(force=args.force)
    if args.assemble:
        with open(MATRIX, "w", encoding="utf-8") as handle:
            handle.write(json.dumps(assemble(args.assemble, args.regenerated), indent=1) + "\n")
        print(f"wrote {os.path.relpath(MATRIX)}")
    text = make()
    if args.check:
        current = open(OUT, encoding="utf-8").read() if os.path.exists(OUT) else ""
        if current != text:
            print("src/planetgen/physics/data/fudge_table.json is out of date; run python scripts/build_fudge_table.py", file=sys.stderr)
            return 1
        return 0
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(text)
    print(f"wrote {os.path.relpath(OUT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
