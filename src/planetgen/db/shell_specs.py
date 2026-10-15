# planetgen/db/shell_specs.py

"""
The galaxy's outline as the Generate page shows it (ADM.28): how many shells
it has, how far out they reach, and how many sectors are charted in each.
Read from `galaxy_shell` (one row per shell with content; see `schema.sql`).
"""

from planetgen.galaxy import spherical

SHELL_SPEC_ROWS = 41
"""int: The most shell rows `galaxy_shell_specs` returns; more shells are an even sample."""


def galaxy_shell_specs(conn, edge_pc):
    """
    The galaxy's shells as the Generate page shows them.

    Returns:
        dict | None: `None` before any shell is stored; else `count` (shells
            with content, the core included), `outermost` (the highest shell
            index), `depth_pc` (one shell's depth, always the sector edge),
            `radius_pc` (the outermost shell's outer edge), `cells` (sectors
            inside the outline, all shells), `charted` (generated
            galaxy-placed sectors, all shells), `sampled` (whether `rows` is
            a sample) and `rows`, outermost shell first: `shell`,
            `inner_pc`, `outer_pc`, `bands` (how many bands the outline
            reaches in it, the caps included), `cells` and `charted`.
    """
    shells = conn.execute(
        "SELECT shell_index, high_band_index, cell_count FROM galaxy_shell ORDER BY shell_index DESC").fetchall()
    if not shells:
        return None
    charted = {row["shell_index"]: row["n"] for row in conn.execute(
        "SELECT shell_index, COUNT(*) AS n FROM sectors WHERE center_x_pc IS NOT NULL GROUP BY shell_index"
    ).fetchall()}
    keep = shells
    if len(shells) > SHELL_SPEC_ROWS:
        picks = {round(i * (len(shells) - 1) / (SHELL_SPEC_ROWS - 1)) for i in range(SHELL_SPEC_ROWS)}
        picks.add(len(shells) - 1)
        keep = [shells[i] for i in sorted(picks)]
    rows = []
    for row in keep:
        inner, outer = spherical.shell_bounds_pc(row["shell_index"], edge_pc)
        rows.append({"shell": row["shell_index"], "inner_pc": inner, "outer_pc": outer,
                     "bands": 2 * row["high_band_index"] + 1, "cells": row["cell_count"],
                     "charted": charted.get(row["shell_index"], 0)})
    outermost = shells[0]["shell_index"]
    return {"count": len(shells), "outermost": outermost, "depth_pc": edge_pc,
            "radius_pc": spherical.shell_bounds_pc(outermost, edge_pc)[1],
            "cells": sum(row["cell_count"] for row in shells), "charted": sum(charted.values()),
            "sampled": len(rows) < len(shells), "rows": rows}
