# planetgen/galaxy/component_field.py

"""
The components' density in the sector gate (GEN.160, GEN.237; docs/design/globular-clusters.md sections 3 and 4).

A `ComponentField` holds a galaxy's components (`galaxy.components.Component`, the `galaxy_components` table): its
globular clusters (a King profile each) and its satellite galaxies (an exponential or Plummer profile each; the
profiles are `physics/component_profiles.py`). It answers, in the units of `galaxy.density.relative_density`, how many
more systems a cell holds for the components in it: the systems they put inside the cell's cube over the systems a
cell holds at density 1 (`expected_system_count_at_density_1`), by kind (`densities`: the `cluster` and `satellite`
populations of the profile's mix). It is carried by `GalaxyShape.components` (`None` for a galaxy without), so every
path that asks `relative_density`, `population_densities` or `predicted_star_count` gets it with no new argument. A
grid of 1 kpc cells names, for any point, the few components whose tidal sphere reaches it, so a cell far from every
component costs one dictionary lookup.

The field is immutable and compares and hashes by identity (the density model caches its terms by shape).
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import satellite_history
from planetgen.physics import component_profiles


KINDS_IN_THE_GATE = ("cluster", "satellite")
"""tuple: The component kinds that put systems in a cell (a stream, later, would too); each is a population of the mix."""

SATELLITE_POPULATIONS = ("satellite_young", "satellite_intermediate", "satellite")
"""tuple: The populations of a satellite's stars (GEN.526): the history's last 0.1 Gy, its last 3 Gy before that, and the
rest."""

POPULATIONS_IN_THE_GATE = ("cluster", *SATELLITE_POPULATIONS)
"""tuple: The populations of the mix that components put in a cell, which the scatters leave to the component pass and
the sector's own draw counts as none (a satellite's three are one kind, `satellite`)."""

KIND_OF_POPULATION = {"cluster": "cluster", "satellite": "satellite", "satellite_young": "satellite",
                      "satellite_intermediate": "satellite"}
"""dict: Population -> the component kind whose stars it holds."""


def population_shares(component):
    """
    The shares of a component's stars by population, `{population: share}` adding up to 1: a cluster's are all `cluster`;
    a satellite's are its history's (`galaxy.satellite_history`), young, intermediate and the rest.
    """
    if component.kind != "satellite":
        return {component.kind: 1.0}
    found = satellite_history.history(component.mass_solar, component.half_mass_radius_pc, component.age_gy)
    return {"satellite_young": found.young_share, "satellite_intermediate": found.intermediate_share,
            "satellite": found.old_share}



class ComponentField:
    """
    The clusters and satellites of a galaxy, for the sector gate.

    Args:
        found (iterable[galaxy.components.Component]): The galaxy's components; those of a kind in
            `KINDS_IN_THE_GATE` count here.
        edge_pc (float): The sector edge, parsecs: the cube a cell's density is integrated over.
        expected_at_density_1 (float): The systems a cell holds at `relative_density` 1.
    """

    __slots__ = ("edge_pc", "expected_at_density_1", "names", "kinds", "fe_h", "shares", "_centers", "_reach",
                 "_scale_pc", "_central", "_tables", "_index")

    def __init__(self, found, edge_pc, expected_at_density_1):
        held = [c for c in found if c.kind in KINDS_IN_THE_GATE]
        self.edge_pc = float(edge_pc)
        self.expected_at_density_1 = float(expected_at_density_1)
        self.names = tuple(c.name for c in held)
        self.kinds = tuple(c.kind for c in held)
        self.fe_h = tuple(c.metallicity_fe_h for c in held)
        self.shares = tuple(population_shares(c) for c in held)
        self._centers = np.array([c.center_pc for c in held], dtype=float).reshape(-1, 3)
        self._reach = np.array([c.tidal_radius_pc for c in held], dtype=float)
        self._scale_pc = np.array([component_profiles.scale_radius_pc(c) for c in held], dtype=float)
        # Systems per pc3 at the centre: the central density of the mass over the mean system mass.
        self._central = np.array([component_profiles.central_density_per_pc3(c) for c in held],
                                 dtype=float) / tuning.CLUSTER_MEAN_SYSTEM_MASS_SOL
        self._tables = tuple(component_profiles.table_for(c) for c in held)
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
        """Grid cell -> the components whose tidal sphere, grown by the cell cube's reach, touches it."""
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
        """Component `i`'s systems per cubic parsec at `radius_pc` from its centre (0 beyond the tidal radius)."""
        table = self._tables[i]
        x = radius_pc / self._scale_pc[i]
        if x >= table.tidal_x:
            return 0.0
        return float(self._central[i] * np.interp(x, table.x, table.density))

    def _cube_systems(self, i, position):
        """The systems of component `i` inside the cube of the sector edge about `position`."""
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

    def _contributions(self, position_pc):
        """`[(component index, systems in the cell's cube)]` for the components that reach `position_pc`."""
        near = self._index.get(self._grid_cell(position_pc))
        if not near:
            return []
        position = np.asarray(position_pc, dtype=float)
        return [(i, self._cube_systems(i, position)) for i in near]

    def densities(self, position_pc):
        """
        The components' part of the density at `position_pc` by population, in the units of `relative_density`: the
        systems each population puts in the cell's cube over the systems a cell holds at density 1. A cluster's stars are
        `cluster`; a satellite's are split by its history into `satellite_young`, `satellite_intermediate` and
        `satellite` (GEN.526). Every population in `POPULATIONS_IN_THE_GATE` is there (`0.0` away from its components).
        """
        total = dict.fromkeys(POPULATIONS_IN_THE_GATE, 0.0)
        for i, systems in self._contributions(position_pc):
            for population, share in self.shares[i].items():
                total[population] += systems * share
        return {population: systems / self.expected_at_density_1 for population, systems in total.items()}

    def relative_density(self, position_pc, kind=None):
        """
        The components' part of `relative_density` at `position_pc` (of one `kind`, or all): the systems they put in
        the cell's cube over the systems a cell holds at density 1. `0.0` away from every component.
        """
        return sum(systems for i, systems in self._contributions(position_pc)
                   if kind is None or self.kinds[i] == kind) / self.expected_at_density_1

    def fe_h_peaks(self, position_pc, kind):
        """
        The [Fe/H] distribution of the `kind` stars (a population of a satellite's three is the kind `satellite`) at
        `position_pc`, `((weight, mean, sigma), ...)`: one peak per
        component of the kind that reaches the point, weighted by the systems it puts in the cell and centred on its own
        [Fe/H] with `tuning.COMPONENT_FE_H_SIGMA[kind]`. Empty where no component of the kind reaches.
        """
        kind = KIND_OF_POPULATION.get(kind, kind)
        sigma = tuning.COMPONENT_FE_H_SIGMA[kind]
        return tuple((systems, self.fe_h[i], sigma) for i, systems in self._contributions(position_pc)
                     if self.kinds[i] == kind and systems > 0.0)

    def bound(self, r_cyl, z):
        """
        An upper bound of `relative_density` over every point at cylindrical radius `r_cyl` and height `z`, whatever
        its angle: each component at the nearest the ring comes to it, less the cell's reach (`CLUSTER_BOUND_SLACK_EDGES`
        edges), the components added. `0.0` for a field without components.
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
        The highest `|band|` of shell `shell` any component reaches (grown by the bound's slack), or `-1` when none
        does: a component far from the plane is no part of the run of bands the disk and bulge qualify, so the outline
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
