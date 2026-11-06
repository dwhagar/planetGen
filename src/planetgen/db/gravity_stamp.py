# planetgen/db/gravity_stamp.py

"""
Which stored sectors a gravity grid reads, and the stamp that names its cache file (MAP.387, API.24; docs/design/orbits/gravity-map.md)
=======================================================================================================================================

The part of the sector grid job (`db.gravity_grid`) that needs no compiled code, so the web app can tell whether a sector's grid is cached and current
(`api/gravity`) without importing numba: the settings of the grid, the sectors in a target's reach and the stamp of the rule, the settings, the potential and
those sectors' rows (`galaxy.gravity_cache` stores and reads the files).
"""

import numpy as np

from planetgen.galaxy import gravity_cache, gravity_settings

ZONES = 16
"""int: Zones along each edge of a sector (0.25 pc at 4 pc; a setting of the map)."""
REACH_EDGES = 4
"""int: How many sector edges away, along each axis, a sector still adds its stars to the grid."""
SUBDIVISION = 8
"""int: Point masses along each edge of the cube the smooth galaxy's own stars of one sector are spread through."""
OWN_SOFTENING_PC = 0.25
"""float: The Plummer softening of those point masses, parsecs (a zone edge at the default grid): they stand for a continuous density."""
_TABLES = ("star_systems", "black_holes", "neutron_stars")


def sector_row(conn, sector_id):
    row = conn.execute("SELECT id, edge_mpc, shell_index, band_index, wedge_index, center_x_pc, center_y_pc, center_z_pc FROM sectors WHERE id = ?",
                       (sector_id,)).fetchone()
    if row is None:
        raise ValueError(f"no sector {sector_id}")
    return row


def sectors_in_reach(conn, target, reach_edges=REACH_EDGES):
    """The rows (id, edge, address, center) of the stored sectors in the cube of `reach_edges` edges round `target` (a sector row), in id order."""
    reach = reach_edges * target["edge_mpc"] / 1000.0
    cx, cy, cz = target["center_x_pc"], target["center_y_pc"], target["center_z_pc"]
    return conn.execute(
        "SELECT id, edge_mpc, shell_index, band_index, wedge_index, center_x_pc, center_y_pc, center_z_pc FROM sectors"
        " WHERE center_x_pc BETWEEN ? AND ? AND center_y_pc BETWEEN ? AND ? AND center_z_pc BETWEEN ? AND ? ORDER BY id",
        (cx - reach, cx + reach, cy - reach, cy + reach, cz - reach, cz + reach)).fetchall()


def source_stamp(conn, sectors, zones=ZONES, reach_edges=REACH_EDGES, pot=None):
    """The stamp of the grid of a target whose reach holds `sectors` (rows from `sectors_in_reach`) under `pot` (see the module doc)."""
    ids = [row["id"] for row in sectors]
    marks = ", ".join("?" for _ in ids)
    parts = [gravity_cache.RULE_VERSION, zones, reach_edges, SUBDIVISION, OWN_SOFTENING_PC, gravity_settings.LEAF_PC, gravity_settings.PAIR_BUDGET,
             gravity_settings.MIN_THETA, gravity_settings.MAX_THETA, gravity_settings.THETA_STEP, repr(pot)]
    modified = {row["id"]: str(row["modified_at"]) for row in conn.execute(f"SELECT id, modified_at FROM sectors WHERE id IN ({marks})", ids).fetchall()} if ids else {}
    rows = {}
    for table in _TABLES:
        sector_column = "sector_id"
        if ids:
            for row in conn.execute(
                    f"SELECT {sector_column} AS sector_id, COUNT(*) AS n, MAX(modified_at) AS latest FROM {table} WHERE {sector_column} IN ({marks})"
                    f" GROUP BY {sector_column}", ids).fetchall():
                rows[(row["sector_id"], table)] = (row["n"], str(row["latest"]))
    for sector_id in sorted(ids):
        parts.append((sector_id, modified.get(sector_id), tuple(rows.get((sector_id, table)) for table in _TABLES)))
    return gravity_cache.stamp(tuple(parts))


def zone_points(center, edge_pc, zones=ZONES):
    """The centers of the `zones^3` zones of the cube of `edge_pc` about `center`, `(zones^3, 3)`, zone `(i, j, k)` at row `(i * zones + j) * zones + k`."""
    c = (np.arange(zones) + 0.5) / zones * edge_pc - edge_pc / 2.0
    x, y, z = np.meshgrid(c, c, c, indexing="ij")
    return np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1) + np.asarray(center, dtype=np.float64)
