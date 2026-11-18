# tests/golden_stages.py

"""
One digest per generation stage of the golden galaxy (OPS.147), so a train's declared stages can be checked against what it moved.

`golden_seed.galaxy_record` digests the plan, each sector and the region, but the plan digest is taken after the whole run, so it moves
with what the later stages write (the audit of OPS.133, `reports/lane05/affects-audit-2026-11-11.md`: it moved at seven epochs while the
plan stage alone moved at one). This module digests each stage's own output instead, from the same stored galaxy:

- `plan`: the plan's own rows, `galaxy_shape` without the columns the later stages write (`PLAN_OWNERS`) and `galaxy_shell`;
- `scatter`: the scatter's tables (`bright_stars`, `phenomenon_scatter`, `scatter_members`) and the galaxy-shape columns the scatter writes;
- `fill`: every sector's own sections but the scatter's and the facilities' (the systems, stars, planets, moons, belts, comets);
- `population`: the facilities. The golden galaxy has none today, so this one is `None` until it does (a stage that cannot be measured is
  never reported as moved);
- `derived` (OPS.133): values that are a pure function of stored ones. They are stored inside the fill's rows, so the digest of its own cannot be told
  from the fill's; it is `None` (not measured) and a note that declares it is never judged against the golden.

Each is `{"full": ..., "structural": ..., "rows": n}`: a stage moved when its `full` digest differs; `structural` says it changed shape and not
only numbers. The golden test writes them to the file named by `PLANETGEN_STAGE_DIGESTS_OUT` when it builds the galaxy (no second build), and
`python scripts/train.py stages --before A --after B` compares two such files with the stages a train declared.
"""

import json
from pathlib import Path

from planetgen.db import fingerprint
from planetgen.util import fingerprint_encoding as encoding

STAGES = ("plan", "scatter", "fill", "population", "derived")
"""tuple: The stages, in the order each is redone from (`bump_version.GENERATION_STAGES`)."""

OUT_ENVIRONMENT = "PLANETGEN_STAGE_DIGESTS_OUT"
"""str: The environment variable that names the file the golden test writes the stage digests to."""

SCATTER_TABLES = ("bright_stars", "phenomenon_scatter", "scatter_members")
"""tuple: The scatter stage's tables."""

SCATTER_OWNED_SECTIONS = frozenset({"phenomenon", "population"})
"""frozenset: The sector sections the fill stage does not own: the scatter's bright stars and phenomena, and the facilities."""

SHAPE_COLUMNS_OF_LATER_STAGES = {
    "scatter": frozenset({"bright_star_min_luminosity_sol", "bright_star_seed", "bright_star_mass_limit_sol",
                          "phenomenon_scatter_seed", "phenomenon_min_mass_solar"}),
    "fill": frozenset({"density_ratio_avg"}),
}
"""dict: The `galaxy_shape` columns a stage after the plan writes, by the stage (`db/store.py`): the plan's own digest leaves them out, and
the scatter's holds the ones it writes."""


def _digest(rows):
    found = encoding.row_set_digests(rows)
    return {"full": found.full, "structural": found.structural, "rows": len(rows)}


def _without(value, columns):
    row = {name: item for name, item in value["row"].items() if name not in columns}
    return {"table": value["table"], "row": row}


def _only(value, columns):
    return {"table": value["table"], "row": {name: item for name, item in value["row"].items() if name in columns}}


def stage_digests(conn):
    """
    The digest of each stage's own output of the stored galaxy.

    Returns:
        dict: `{stage: {"full", "structural", "rows"} or None}` for each stage in `STAGES`.
    """
    reader = fingerprint._Reader(conn)
    later = frozenset().union(*SHAPE_COLUMNS_OF_LATER_STAGES.values())
    shape = [reader.row_value("galaxy_shape", row) for row in reader.select("galaxy_shape", "1 = 1")]
    plan = [_without(value, later) for value in shape]
    plan += [reader.row_value("galaxy_shell", row) for row in reader.select("galaxy_shell", "1 = 1")]
    scatter = [_only(value, SHAPE_COLUMNS_OF_LATER_STAGES["scatter"]) for value in shape]
    for table in SCATTER_TABLES:
        scatter += [reader.row_value(table, row) for row in reader.select(table, "1 = 1")]
    fill = []
    for (label,) in sorted((f"{row['shell_index']} {row['band_index']} {row['wedge_index']}",) for row in conn.execute(
            "SELECT shell_index, band_index, wedge_index FROM sectors").fetchall()):
        address = tuple(int(part) for part in label.split())
        for name, value in sorted(fingerprint.sector_sections(conn, address=address).items()):
            if name not in SCATTER_OWNED_SECTIONS:
                fill.append({"sector": label, "section": name, "full": value.full, "structural": value.structural})
    facilities = [reader.row_value("facilities", row) for row in reader.select("facilities", "1 = 1")]
    return {"plan": _digest(plan), "scatter": _digest(scatter), "fill": _digest(fill),
            "population": _digest(facilities) if facilities else None, "derived": None}


def write(conn, path):
    """Writes `stage_digests(conn)` to `path` as JSON."""
    Path(path).write_text(json.dumps(stage_digests(conn), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def read(path):
    """The stage digests in the file at `path`."""
    return json.loads(Path(path).read_text(encoding="utf-8"))


def moved(before, after):
    """The stages whose digest differs between two records, in stage order. A stage either record could not measure (`None`) is not moved."""
    found = []
    for stage in STAGES:
        old, new = before.get(stage), after.get(stage)
        if old is not None and new is not None and old["full"] != new["full"]:
            found.append(stage)
    return found
