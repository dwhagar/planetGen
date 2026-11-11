# tests/audit/ledger_light_census.py

"""
GEN.1465: what the stored light and ultraviolet rows of the effect ledger say, measured on generated systems.

Light and UV are the two effects the generator fills today (`StarSystem._fill_light_ledger`); a body with one star sending it light writes no row (the source is
implied), so every row there is a body with two or more stars. This audit reads the rows `EffectLedger.rows` would store (the same call `store.insert_star_system`
makes) for each generated system and tabulates, for each effect: how many bodies have a second source, which star is the largest source (the body's host or a
companion), how large the second source is as a share of the largest, and how many bodies have a source over 1 percent, 10 percent and 50 percent. The database
holds the same numbers: a (body, effect) row with `share` under 1.0 is a second source, and the summary row carries `n_below` and `largest_value`. Run:

    python -m tests.audit.ledger_light_census 300 [first_seed] [binary]

where `binary` forces a multiple-star system (the natural mix is mostly lone stars, which store nothing) and the default config is used otherwise.
"""

import collections
import statistics
import sys

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import effect_ledger
from planetgen.util import draw

THRESHOLDS = (0.001, 0.01, 0.1, 0.5)


def ids_of(system):
    """`id(object)` to `(table, id)` for a system's stars, planets and moons, numbered in order (what the save would give)."""
    ids, counter = {}, 0
    for star in system.stars:
        counter += 1
        ids[id(star)] = ("stars", counter)
    for obj in system.all_planets:
        counter += 1
        ids[id(obj)] = ("planets", counter)
        for moon in getattr(obj, "moons", ()):
            counter += 1
            ids[id(moon)] = ("moons", counter)
    return ids


def measure(systems):
    """The tallies of `systems` (a list of `StarSystem`): per effect, bodies with rows, second-source shares, largest-source role."""
    codes = {code: name for name, code, *_rest in effect_ledger.EFFECTS}
    tally = {name: {"bodies": 0, "shares": [], "largest_host": 0, "largest_companion": 0, "by_kind": collections.Counter()} for name in codes.values()}
    totals = {"systems": len(systems), "multiple": 0, "bodies": 0, "rows": 0}
    for system in systems:
        ids = ids_of(system)
        stars = [ids[id(star)][1] for star in system.stars]
        totals["multiple"] += len(stars) > 1
        totals["bodies"] += sum(1 for table, _ in ids.values() if table != "stars")
        rows = system.effect_ledger.rows(1, ids)
        totals["rows"] += len(rows)
        grouped = collections.defaultdict(list)
        for row in rows:
            grouped[(row["effect"], row["body_type"], row["body_id"])].append(row)
        for (effect, body_type, _body), group in grouped.items():
            row = tally[codes[effect]]
            row["bodies"] += 1
            row["by_kind"][effect_ledger.BODY_TABLES[body_type - 1]] += 1
            sources = [r for r in group if r["source_type"] != 0]
            top = next(r for r in sources if r["share"] == 1.0)
            if top["source_id"] == stars[0]:
                row["largest_host"] += 1
            else:
                row["largest_companion"] += 1
            for source in sources:
                if source is not top:
                    row["shares"].append(source["share"])
            summary = next(r for r in group if r["source_type"] == 0)
            row["shares"].extend([0.0] * (summary["n_below"] or 0))      # a source under the floor: counted, its size is under 0.1 percent
    return totals, tally


def report(totals, tally):
    lines = [f"{totals['systems']} systems, {totals['multiple']} with more than one star, {totals['bodies']} planets and moons, {totals['rows']} rows "
             f"({totals['rows'] / max(totals['systems'], 1):.2f} a system)."]
    for effect, row in tally.items():
        if not row["bodies"]:
            continue
        shares = row["shares"]
        cross = ", ".join(f"over {t * 100:g} percent: {sum(1 for s in shares if s >= t)}" for t in THRESHOLDS)
        lines.append(f"{effect}: {row['bodies']} bodies with a second source; largest source is the host {row['largest_host']} times and a companion {row['largest_companion']} times; "
                     f"second sources by share ({len(shares)} in all; {cross}); median share {statistics.median(shares) if shares else 0:.4f}; by receiver {dict(row['by_kind'])}.")
    return "\n".join(lines)


def main(argv):
    count = int(argv[1]) if len(argv) > 1 else 100
    first = int(argv[2]) if len(argv) > 2 and argv[2].isdigit() else 0
    config = SystemConfig()
    config.PLANETS = True
    if "binary" in argv:
        config.BINARY_SYSTEM = True
    systems = []
    for seed in range(first, first + count):
        with draw.bound(seed):
            systems.append(StarSystem(system_config=config))
    print(report(*measure(systems)))


if __name__ == "__main__":
    main(sys.argv)
