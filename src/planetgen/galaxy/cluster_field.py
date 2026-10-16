# planetgen/galaxy/cluster_field.py

"""
The globular clusters' density in the sector gate (GEN.160; docs/design/globular-clusters.md sections 3 and 4).

A `ClusterField` holds a galaxy's King-profile cluster components (`galaxy.components.Component`, the
`galaxy_components` table) and answers, in the units of `galaxy.density.relative_density`, how many more systems a cell
holds for the clusters in it: the systems the clusters put inside the cell's cube over the systems a cell holds at
density 1 (`expected_system_count_at_density_1`). It is carried by `GalaxyShape.clusters` (`None` for a galaxy without),
so every path that asks `relative_density`, `population_densities` or `predicted_star_count` gets it with no new
argument. A grid of 1 kpc cells names, for any point, the few clusters whose tidal sphere reaches it, so a cell far
from every cluster costs one dictionary lookup.

The field is immutable and compares and hashes by identity (the density model caches its terms by shape).
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.physics import king


class ClusterField:
    """
    The King-profile clusters of a galaxy, for the sector gate.

    Args:
        found (iterable[galaxy.components.Component]): The galaxy's components; only the `cluster` kind with a `king`
            profile count here.
        edge_pc (float): The sector edge, parsecs: the cube a cell's density is integrated over.
        expected_at_density_1 (float): The systems a cell holds at `relative_density` 1.
    """

    __slots__ = ("edge_pc", "expected_at_density_1", "names", "_centers", "_reach", "_scale_pc", "_central",
                 "_tables", "_index")

    def __init__(self, found, edge_pc, expected_at_density_1):
        kings = [c for c in found if c.kind == "cluster" and c.profile == "king"]
        self.edge_pc = float(edge_pc)
        self.expected_at_density_1 = float(expected_at_density_1)
        self.names = tuple(c.name for c in kings)
        self._centers = np.array([c.center_pc for c in kings], dtype=float).reshape(-1, 3)
        self._reach = np.array([c.tidal_radius_pc for c in kings], dtype=float)
        self._scale_pc = np.array([king.scale_radius_pc(c.king_w0, c.half_mass_radius_pc) for c in kings], dtype=float)
        # Systems per pc3 at the centre: the central density of the mass over the mean system mass.
        self._central = np.array(
            [king.central_density_per_pc3(c.king_w0, c.mass_solar, c.half_mass_radius_pc) for c in kings], dtype=float
        ) / tuning.CLUSTER_MEAN_SYSTEM_MASS_SOL
        self._tables = tuple(king.king_table(int(c.king_w0)) for c in kings)
        self._index = self._build_index()

    __hash__ = object.__hash__

    def __eq__(self, other):
        return self is other

    def __len__(self):
        return len(self.names)

    def _grid_cell(self, point):
        grid = tuning.CLUSTER_GRID_PC
        return (math.floor(point[0] / grid), math.floor(point[1] / grid), math.floor(point[2] / grid))

    def _build_index(self):
        """Grid cell -> the clusters whose tidal sphere, grown by the cell cube's reach, touches it."""
        grid = tuning.CLUSTER_GRID_PC
        margin = self.edge_pc * (tuning.CLUSTER_NEAR_EDGES + 1.0)
        index = {}
        for i, (center, reach) in enumerate(zip(self._centers, self._reach)):
            low = self._grid_cell(center - (reach + margin))
            high = self._grid_cell(center + (reach + margin))
            for a in range(low[0], high[0] + 1):
                for b in range(low[1], high[1] + 1):
                    for c in range(low[2], high[2] + 1):
                        index.setdefault((a, b, c), []).append(i)
        return {key: tuple(value) for key, value in index.items()}

    def _systems_per_pc3(self, i, radius_pc):
        """Cluster `i`'s systems per cubic parsec at `radius_pc` from its centre (0 beyond the tidal radius)."""
        table = self._tables[i]
        x = radius_pc / self._scale_pc[i]
        if x >= table.tidal_x:
            return 0.0
        return float(self._central[i] * np.interp(x, table.x, table.density))

    def _cube_systems(self, i, position):
        """The systems of cluster `i` inside the cube of the sector edge about `position`."""
        edge = self.edge_pc
        center = self._centers[i]
        distance = float(np.linalg.norm(position - center))
        if distance - 0.87 * edge >= self._reach[i]:
            return 0.0
        if distance > tuning.CLUSTER_NEAR_EDGES * edge:
            return self._systems_per_pc3(i, distance) * edge ** 3
        n = tuning.CLUSTER_CELL_SUBSAMPLES
        if self._scale_pc[i] < edge / 4.0:
            n *= 2          # a core much smaller than the cell: a finer grid for the few cells it sits in
        offsets = (np.arange(n) + 0.5) / n * edge - edge / 2.0
        gx, gy, gz = np.meshgrid(offsets, offsets, offsets, indexing="ij")
        radius = np.sqrt((position[0] + gx - center[0]) ** 2 + (position[1] + gy - center[1]) ** 2
                         + (position[2] + gz - center[2]) ** 2)
        table = self._tables[i]
        x = radius / self._scale_pc[i]
        density = np.where(x < table.tidal_x, np.interp(x, table.x, table.density), 0.0)
        return float(self._central[i] * density.sum() * (edge / n) ** 3)

    def relative_density(self, position_pc):
        """
        The clusters' part of `relative_density` at `position_pc`: the systems they put in the cell's cube over
        the systems a cell holds at density 1. `0.0` away from every cluster.
        """
        near = self._index.get(self._grid_cell(position_pc))
        if not near:
            return 0.0
        position = np.asarray(position_pc, dtype=float)
        return sum(self._cube_systems(i, position) for i in near) / self.expected_at_density_1

    def bound(self, r_cyl, z):
        """
        An upper bound of `relative_density` over every point at cylindrical radius `r_cyl` and height `z`, whatever
        its angle: each cluster at the nearest the ring comes to it, less the cell's reach (`CLUSTER_BOUND_SLACK_EDGES`
        edges), the clusters added. `0.0` for a field without clusters.
        """
        if not len(self):
            return 0.0
        cylindrical = np.hypot(self._centers[:, 0], self._centers[:, 1])
        nearest = np.hypot(r_cyl - cylindrical, z - self._centers[:, 2]) - tuning.CLUSTER_BOUND_SLACK_EDGES * self.edge_pc
        total = 0.0
        for i in np.nonzero(nearest < self._reach)[0]:
            total += self._systems_per_pc3(int(i), max(float(nearest[i]), 0.0)) * self.edge_pc ** 3
        return total / self.expected_at_density_1

    def high_band(self, shell, band_at, edge_pc):
        """
        The highest `|band|` of shell `shell` any cluster reaches (grown by the bound's slack), or `-1` when none
        does: a cluster far from the plane is no part of the run of bands the disk and bulge qualify, so the outline
        takes it in by hand. `band_at` is `galaxy.spherical.band_at`.
        """
        high = -1
        radius = shell * edge_pc
        slack = tuning.CLUSTER_BOUND_SLACK_EDGES * edge_pc
        for center, reach in zip(self._centers, self._reach):
            distance = float(np.linalg.norm(center))
            grown = reach + slack
            if distance <= grown:
                return band_at(shell, 0.0) if shell else 0   # the sphere holds the pole: the shell's top band
            if abs(distance - radius) > grown:
                continue
            latitude = math.asin(max(-1.0, min(1.0, center[2] / distance)))
            half_angle = math.asin(min(1.0, grown / distance))
            highest = min(math.pi / 2.0, abs(latitude) + half_angle)
            high = max(high, abs(band_at(shell, math.pi / 2.0 - highest)))
        return high
