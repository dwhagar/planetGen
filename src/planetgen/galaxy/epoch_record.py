# planetgen/galaxy/epoch_record.py

"""
The generator epoch each stage of a stored galaxy was made under (GEN.764).

A galaxy is made in stages (`util/generator_epoch.stages()`), and each keeps its own record of the epoch: the plan in
`galaxy_shape.generator_epoch` (v131), each scatter in `scatter_epochs` (v164), each sector in `sectors.generator_epoch`
(v121) and each species in `species.generator_epoch` (v164). `recorded(conn)` reads them all; `behind(conn)` says, per
stage, what is older than the running epoch and which stages that makes a galaxy redo
(`util/generator_epoch.stages_between`). `None` is an epoch the record does not know (a galaxy made before it was kept).
"""

from planetgen.db import store
from planetgen.util import generator_epoch


def recorded(conn):
    """
    The epochs the galaxy's stages were made under.

    Returns:
        dict: `plan` (the plan's epoch, or None), `scatter` (`{kind: epoch}` of the scatters drawn), `fill`
        (`{epoch: sectors}`), `population` (`{epoch: species}`). The epoch keys of `fill` and `population` can be None.
    """
    plan = conn.execute("SELECT generator_epoch FROM galaxy_shape WHERE id = 1").fetchone()
    return {
        "plan": None if plan is None else plan["generator_epoch"],
        "scatter": store.scatter_epochs(conn),
        "fill": _counts(conn, "sectors"),
        "population": _counts(conn, "species"),
    }


def _counts(conn, table):
    rows = conn.execute(f"SELECT generator_epoch, COUNT(*) AS n FROM {table} GROUP BY generator_epoch").fetchall()
    return {row["generator_epoch"]: row["n"] for row in sorted(
        rows, key=lambda row: (row["generator_epoch"] is None, row["generator_epoch"] or 0))}


def needs(stage, epoch, current):
    """Whether a result of `stage` made under `epoch` must be redone for `current`: some epoch since changed that
    stage. An unknown epoch (None) counts as the first one, so whatever any epoch since changed."""
    return stage in generator_epoch.stages_between(1 if epoch is None else min(epoch, current), current)


def behind(conn, current=None):
    """
    What must be redone to bring the galaxy to epoch `current` (default the running one).

    Returns:
        dict: `stale`, the records that are out of date for their own stage: `plan` (its epoch, only when the plan is),
        `scatter` (`{kind: epoch}`), `fill` and `population` (`{epoch: count}`); and `stages`, the stages to redo in the
        order of `generator_epoch.stages()`: those with a stale record, and everything that is made from a stage
        redone (a new plan redoes all, a scatter redoes the fill, a fill redoes the population).
    """
    current = generator_epoch.generator_epoch() if current is None else current
    found = recorded(conn)
    stale = {"scatter": {kind: epoch for kind, epoch in found["scatter"].items() if needs("scatter", epoch, current)},
             "fill": {epoch: n for epoch, n in found["fill"].items() if needs("fill", epoch, current)},
             "population": {epoch: n for epoch, n in found["population"].items() if needs("population", epoch, current)}}
    if found["plan"] is not None and needs("plan", found["plan"], current):
        stale["plan"] = found["plan"]
    stages = generator_epoch.stages()
    redo = {stage for stage in stages if stale.get(stage)}
    for stage, built_from in (("scatter", "plan"), ("fill", "scatter"), ("population", "fill")):
        if built_from in redo:
            redo.add(stage)
    if "plan" in redo:
        redo.update(stages)
    return {"stale": stale, "stages": tuple(stage for stage in stages if stage in redo)}
