# planetgen/db/gravity_grid.py

"""
The gravity map's sector grid: the local field of one sector on 16 x 16 x 16 zones, cached (MAP.387; docs/design/orbits/gravity-map.md)
===========================================================================================================================================

`grid(conn, sector_id, directory)` returns the *local* gravity field of a sector at the centre of each of its zones, from the rule of MAP.383: **the generated
stars of the sectors round it, less the smooth galaxy's own stars in those sectors**. The total field is that plus the smooth galaxy at the zone
(`galaxy.potential`), which the reader adds because it is a closed form and keeps the file small. The rule makes the double count of the smooth potential and
the generated stars vanish: in the tidal strength it was 36 percent at the median (`reports/map383-background-rule-2026-11-05.md`).

- **Sectors.** Every stored sector whose centre lies within `REACH_EDGES` sector edges of the target's along each axis (a cube of `2 * REACH_EDGES + 1` edges, the
  target in the middle). A sector that is not stored adds and subtracts nothing. Farther generated sectors are left out; their mass differs from the smooth
  model's by the shot noise of their stars, which falls with the square of the distance.
- **Generated stars.** `gravity_sources.system_point_masses` (one point per system, MAP.386) summed by `gravity_tree` (the 0.5 pc hierarchy of MAP.385) and
  `gravity_kernel` (MAP.384), the opening angle chosen to the pair budget.
- **The smooth galaxy's own stars.** For each sector the mass of the potential's own stellar density at its centre times the cell's volume
  (`potential.stellar_density`), spread through a cube of that volume as `SUBDIVISION^3` equal point masses softened by `OWN_SOFTENING_PC` (the zone
  edge) and summed the same way. It is a continuous density where the stars are points, so it is softened and the stars are not.
- **The zones.** A cube of the sector's edge about its centre. The cells of the galaxy are shells, bands and wedges, not cubes; the map draws the cube.
- **Cache.** One file per sector (`gravity_cache`), keyed by a stamp of the rule, the settings, the potential and, for every sector in reach, its `modified_at` and
  the count and latest `modified_at` of its rows; a moved star changes the stamp of the grids that reach it and of no other.

The module imports numba through `gravity_tree`, so the web app must not import it; the web reads the files through `galaxy.gravity_cache` (OPS.42). Units:
positions in parsecs, masses in solar masses, the field in `gravity_kernel`'s.
"""

import dataclasses
import math

import numpy as np

from planetgen.db import gravity_sources
from planetgen.galaxy import gravity_cache, gravity_kernel, gravity_tree, potential, spherical

ZONES = 16
"""int: Zones along each edge of a sector (0.25 pc at 4 pc; a setting of the map)."""
REACH_EDGES = 4
"""int: How many sector edges away, along each axis, a sector still adds its stars to the grid."""
SUBDIVISION = 8
"""int: Point masses along each edge of the cube the smooth galaxy's own stars of one sector are spread through."""
OWN_SOFTENING_PC = 0.25
"""float: The Plummer softening of those point masses, parsecs (a zone edge at the default grid): they stand for a continuous density."""
_TABLES = ("star_systems", "black_holes", "neutron_stars")


@dataclasses.dataclass(frozen=True)
class GridResult:
    """A sector's grid: `values` `(zones^3, 10)` float32, `hit` whether the file was there, `stamp` the key of its file."""

    values: np.ndarray
    hit: bool
    stamp: str


def _sector(conn, sector_id):
    row = conn.execute("SELECT id, edge_mpc, shell_index, band_index, wedge_index, center_x_pc, center_y_pc, center_z_pc FROM sectors WHERE id = ?",
                       (sector_id,)).fetchone()
    if row is None:
        raise ValueError(f"no sector {sector_id}")
    return row


def sectors_in_reach(conn, target, reach_edges=REACH_EDGES):
    """The rows (id, edge, address, centre) of the stored sectors in the cube of `reach_edges` edges round `target` (a sector row), in id order."""
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
    parts = [gravity_cache.RULE_VERSION, zones, reach_edges, SUBDIVISION, OWN_SOFTENING_PC, gravity_tree.LEAF_PC, gravity_tree.PAIR_BUDGET,
             gravity_tree.MIN_THETA, gravity_tree.MAX_THETA, gravity_tree.THETA_STEP, repr(pot)]
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


def zone_points(centre, edge_pc, zones=ZONES):
    """The centres of the `zones^3` zones of the cube of `edge_pc` about `centre`, `(zones^3, 3)`, zone `(i, j, k)` at row `(i * zones + j) * zones + k`."""
    c = (np.arange(zones) + 0.5) / zones * edge_pc - edge_pc / 2.0
    x, y, z = np.meshgrid(c, c, c, indexing="ij")
    return np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1) + np.asarray(centre, dtype=np.float64)


def own_stars(sectors, pot, subdivision=SUBDIVISION):
    """
    The smooth galaxy's own stars in `sectors`: `(positions, masses)` of `subdivision^3` equal point masses per sector, in a cube of the sector's volume about its
    centre, the mass being the potential's stellar density at the centre times the cell's volume.
    """
    offsets = (np.arange(subdivision) + 0.5) / subdivision - 0.5
    x, y, z = np.meshgrid(offsets, offsets, offsets, indexing="ij")
    unit = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
    positions, masses = [], []
    for row in sectors:
        edge_pc = row["edge_mpc"] / 1000.0
        volume = spherical.cell_volume(row["shell_index"], row["band_index"], row["wedge_index"], edge_pc)
        centre = np.array([row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]])
        density = potential.stellar_density(pot, *centre)
        mass = float(density) * volume if math.isfinite(float(density)) else 0.0
        positions.append(unit * volume ** (1.0 / 3.0) + centre)
        masses.append(np.full(len(unit), mass / len(unit)))
    if not positions:
        return np.zeros((0, 3)), np.zeros(0)
    return np.concatenate(positions), np.concatenate(masses)


def compute(conn, sector_id, pot, zones=ZONES, reach_edges=REACH_EDGES, threads=None):
    """The local field of the sector on `zones^3` zones, float64 `(zones^3, 10)`: generated stars less the smooth galaxy's own (see the module doc)."""
    target = _sector(conn, sector_id)
    sectors = sectors_in_reach(conn, target, reach_edges)
    edge_pc = target["edge_mpc"] / 1000.0
    points = zone_points((target["center_x_pc"], target["center_y_pc"], target["center_z_pc"]), edge_pc, zones)
    stars = gravity_sources.system_point_masses(conn, [row["id"] for row in sectors])
    generated = gravity_tree.fields_at(gravity_tree.build(stars.positions_pc, stars.masses_solar), points, threads=threads).values
    own_positions, own_masses = own_stars(sectors, pot)
    own = gravity_tree.fields_at(gravity_tree.build(own_positions, own_masses), points, softening_pc=OWN_SOFTENING_PC, threads=threads).values
    return generated - own


def grid(conn, sector_id, directory, pot=None, zones=ZONES, reach_edges=REACH_EDGES, threads=None):
    """
    The `GridResult` of the sector: read from `directory` when a file of the current stamp is there, else worked out and stored. `pot` is the galaxy's
    `potential.Potential` (`potential.current()` when not given).
    """
    pot = potential.current() if pot is None else pot
    target = _sector(conn, sector_id)
    stamp = source_stamp(conn, sectors_in_reach(conn, target, reach_edges), zones, reach_edges, pot)
    cached = gravity_cache.read(directory, sector_id, zones, stamp)
    if cached is not None:
        return GridResult(cached, True, stamp)
    values = compute(conn, sector_id, pot, zones, reach_edges, threads).astype(np.float32)
    gravity_cache.write(directory, sector_id, zones, stamp, values)
    return GridResult(values, False, stamp)
