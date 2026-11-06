# planetgen/db/gravity_grid.py

"""
The gravity map's sector grid: the local field of one sector on 16 x 16 x 16 zones, cached (MAP.387; docs/design/orbits/gravity-map.md)
===========================================================================================================================================

`grid(conn, sector_id, directory)` returns the *local* gravity field of a sector at the center of each of its zones, from the rule of MAP.383: **the generated
stars of the sectors round it, less the smooth galaxy's own stars in those sectors**. The total field is that plus the smooth galaxy at the zone
(`galaxy.potential`), which the reader adds because it is a closed form and keeps the file small. The rule makes the double count of the smooth potential and
the generated stars vanish: in the tidal strength it was 36 percent at the median (`reports/map383-background-rule-2026-11-05.md`).

- **Sectors.** Every stored sector whose center lies within `REACH_EDGES` sector edges of the target's along each axis (a cube of `2 * REACH_EDGES + 1` edges, the
  target in the middle). A sector that is not stored adds and subtracts nothing. Farther generated sectors are left out; their mass differs from the smooth
  model's by the shot noise of their stars, which falls with the square of the distance.
- **Generated stars.** `gravity_sources.system_point_masses` (one point per system, MAP.386) summed by `gravity_tree` (the 0.5 pc hierarchy of MAP.385) and
  `gravity_kernel` (MAP.384), the opening angle chosen to the pair budget.
- **The smooth galaxy's own stars.** For each sector the mass of the potential's own stellar density at its center times the cell's volume
  (`potential.stellar_density`), spread through a cube of that volume as `SUBDIVISION^3` equal point masses softened by `OWN_SOFTENING_PC` (the zone
  edge) and summed the same way. It is a continuous density where the stars are points, so it is softened and the stars are not.
- **The zones.** A cube of the sector's edge about its center. The cells of the galaxy are shells, bands and wedges, not cubes; the map draws the cube.
- **Cache.** One file per sector (`gravity_cache`), keyed by a stamp of the rule, the settings, the potential and, for every sector in reach, its `modified_at` and
  the count and latest `modified_at` of its rows; a moved star changes the stamp of the grids that reach it and of no other.

The module imports numba through `gravity_tree`, so the web app must not import it; the web reads the files through `galaxy.gravity_cache` (OPS.42). Units:
positions in parsecs, masses in solar masses, the field in `gravity_kernel`'s.
"""

import dataclasses
import math

import numpy as np

from planetgen.db import gravity_sources
from planetgen.db.gravity_stamp import OWN_SOFTENING_PC, REACH_EDGES, SUBDIVISION, ZONES, sector_row, sectors_in_reach, source_stamp, zone_points
from planetgen.galaxy import gravity_cache, gravity_kernel, gravity_tree, potential, spherical

@dataclasses.dataclass(frozen=True)
class GridResult:
    """A sector's grid: `values` `(zones^3, 10)` float32, `hit` whether the file was there, `stamp` the key of its file."""

    values: np.ndarray
    hit: bool
    stamp: str


def own_stars(sectors, pot, subdivision=SUBDIVISION):
    """
    The smooth galaxy's own stars in `sectors`: `(positions, masses)` of `subdivision^3` equal point masses per sector, in a cube of the sector's volume about its
    center, the mass being the potential's stellar density at the center times the cell's volume.
    """
    offsets = (np.arange(subdivision) + 0.5) / subdivision - 0.5
    x, y, z = np.meshgrid(offsets, offsets, offsets, indexing="ij")
    unit = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
    positions, masses = [], []
    for row in sectors:
        edge_pc = row["edge_mpc"] / 1000.0
        volume = spherical.cell_volume(row["shell_index"], row["band_index"], row["wedge_index"], edge_pc)
        center = np.array([row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]])
        density = potential.stellar_density(pot, *center)
        mass = float(density) * volume if math.isfinite(float(density)) else 0.0
        positions.append(unit * volume ** (1.0 / 3.0) + center)
        masses.append(np.full(len(unit), mass / len(unit)))
    if not positions:
        return np.zeros((0, 3)), np.zeros(0)
    return np.concatenate(positions), np.concatenate(masses)


def compute(conn, sector_id, pot, zones=ZONES, reach_edges=REACH_EDGES, threads=None):
    """The local field of the sector on `zones^3` zones, float64 `(zones^3, 10)`: generated stars less the smooth galaxy's own (see the module doc)."""
    target = sector_row(conn, sector_id)
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
    target = sector_row(conn, sector_id)
    stamp = source_stamp(conn, sectors_in_reach(conn, target, reach_edges), zones, reach_edges, pot)
    cached = gravity_cache.read(directory, sector_id, zones, stamp)
    if cached is not None:
        return GridResult(cached, True, stamp)
    values = compute(conn, sector_id, pot, zones, reach_edges, threads).astype(np.float32)
    gravity_cache.write(directory, sector_id, zones, stamp, values)
    return GridResult(values, False, stamp)
