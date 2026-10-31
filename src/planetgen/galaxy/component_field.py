# planetgen/galaxy/component_field.py

"""
The components' density in the sector gate (GEN.160, GEN.237; docs/design/galaxy/globular-clusters.md sections 3 and 4).

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
from planetgen.galaxy import components, satellite_history
from planetgen.physics import component_profiles


KINDS_IN_THE_GATE = ("cluster", "satellite")
"""tuple: The component kinds that put systems in a cell (a stream, later, would too); each is a population of the mix."""

SATELLITE_POPULATIONS = ("satellite_young", "satellite_intermediate", "satellite")
"""tuple: The populations of a satellite's stars (GEN.526): the history's last 0.1 Gy, its last 3 Gy before that, and the
rest."""

OPEN_CLUSTER_POPULATIONS = tuning.OPEN_CLUSTER_POPULATIONS
"""tuple: The populations of a young open cluster's stars (GEN.244), one per quarter-decade of age, 1 Myr to 1 Gyr: the
cluster is coeval, so its stars belong to the population whose age window holds its age."""

POPULATIONS_IN_THE_GATE = ("cluster", *SATELLITE_POPULATIONS, *OPEN_CLUSTER_POPULATIONS)
"""tuple: The populations of the mix that components put in a cell, which the scatters leave to the component pass and
the sector's own draw counts as none (a satellite's three are one family, `satellite`, a young cluster's twelve one,
`open_cluster`)."""

FAMILY_OF_POPULATION = {"cluster": "cluster", "satellite": "satellite", "satellite_young": "satellite",
                        "satellite_intermediate": "satellite",
                        **{name: "open_cluster" for name in OPEN_CLUSTER_POPULATIONS}}
"""dict: Population -> the family of components whose stars it holds: `cluster` (the globular clusters), `satellite` or
`open_cluster`."""

KIND_OF_POPULATION = FAMILY_OF_POPULATION


def gate_components(found):
    """The components of `found` that put systems in a cell (a kind of `KINDS_IN_THE_GATE`), in order."""
    return [c for c in found if c.kind in KINDS_IN_THE_GATE]


def age_key(component):
    """The age a component's mean system mass is looked up at (`components.Component.mean_system_mass_solar`): Gy, rounded
    to 6 places."""
    return round(component.age_gy, 6)


def family_of(component):
    """A component's family: its kind, except that a `cluster` with the `plummer` profile is a young `open_cluster`
    (a globular cluster is a `king`)."""
    return components.family_of(component)


def fe_h_peaks_of(component):
    """The [Fe/H] distribution of one component's stars, `((weight, mean, sigma), ...)`: its own value (a globular cluster's,
    a young cluster's or a satellite's) with the scatter of its family (`tuning.COMPONENT_FE_H_SIGMA`). A satellite's three
    populations share it (GEN.662)."""
    return ((1.0, component.metallicity_fe_h, tuning.COMPONENT_FE_H_SIGMA[family_of(component)]),)


def population_shares(component):
    """
    The shares of a component's stars by population, `{population: share}` adding up to 1: a globular cluster's are all
    `cluster`; a young open cluster's all the population of its age (`components.open_cluster_population`); a
    satellite's are its history's (`galaxy.satellite_history`), young, intermediate and the rest.
    """
    if components.is_open_cluster(component):
        return {components.open_cluster_population(component.age_gy): 1.0}
    if component.kind != "satellite":
        return {component.kind: 1.0}
    found = satellite_history.history(component.mass_solar, component.half_mass_radius_pc, component.age_gy)
    return {"satellite_young": found.young_share, "satellite_intermediate": found.intermediate_share,
            "satellite": found.old_share}



_PREFILTER_ABOVE = 24
"""int: Components in one grid cell above which `_contributions` drops the far ones with one vector test."""


class ComponentField:
    """
    The clusters and satellites of a galaxy, for the sector gate.

    Args:
        found (iterable[galaxy.components.Component]): The galaxy's components; those of a kind in
            `KINDS_IN_THE_GATE` count here.
        edge_pc (float): The sector edge, parsecs: the cube a cell's density is integrated over.
        expected_at_density_1 (float): The systems a cell holds at `relative_density` 1.
        mean_by_age (dict, optional): Age (Gy, rounded to 6 places) -> the mean system mass there, solar masses
            (`components.Component.mean_system_mass_solar`), for every component's age. Left out, each is worked out
            here; that is almost all the time a field takes (about half a millisecond per distinct age, so tens of seconds
            for a galaxy of many young clusters), which is why `galaxy.field_job` works them out ahead, in small parts,
            and hands them in.
    """

    __slots__ = ("edge_pc", "expected_at_density_1", "names", "kinds", "families", "fe_h", "shares", "_centers", "_reach",
                 "_scale_pc", "_central", "_tables", "_index", "_big", "_small", "_small_by_cyl", "_small_cyl",
                 "_small_by_dist", "_small_dist", "_small_window", "_density_reach", "_density_window", "_encloses_origin")

    def __init__(self, found, edge_pc, expected_at_density_1, mean_by_age=None):
        held = gate_components(found)
        self.edge_pc = float(edge_pc)
        self.expected_at_density_1 = float(expected_at_density_1)
        self.names = tuple(c.name for c in held)
        self.kinds = tuple(c.kind for c in held)
        self.families = tuple(components.family_of(c) for c in held)
        self.fe_h = tuple(c.metallicity_fe_h for c in held)
        self.shares = tuple(population_shares(c) for c in held)
        self._centers = np.array([c.center_pc for c in held], dtype=float).reshape(-1, 3)
        self._reach = np.array([c.tidal_radius_pc for c in held], dtype=float)
        self._scale_pc = np.array([component_profiles.scale_radius_pc(c) for c in held], dtype=float)
        means = [c.mean_system_mass_solar if mean_by_age is None else mean_by_age[age_key(c)] for c in held]
        # Systems per pc3 at the centre: the central density of the mass over the component's mean system mass.
        self._central = np.array([component_profiles.central_density_per_pc3(c) / mean for c, mean in zip(held, means)],
                                 dtype=float)
        self._tables = tuple(component_profiles.table_for(c) for c in held)
        self._density_reach = self._qualifying_reach([c.mass_solar / mean for c, mean in zip(held, means)])
        self._index = self._build_index()
        self._group_by_reach()

    __hash__ = object.__hash__

    def __eq__(self, other):
        return self is other

    def __len__(self):
        return len(self.names)

    def _qualifying_reach(self, total_systems):
        """
        How far from its centre each component can put a cell at or over one expected system (GEN.627): the radius where
        its density times the cell's volume falls to one system, grown by the cell's half-diagonal, and never beyond its
        tidal radius. A component that holds under one system in all, or whose centre cell would hold under one, never
        qualifies a cell, and its reach is `0`. The outline takes a component in by this reach and not by its tidal one:
        a dwarf galaxy's tidal sphere (15.6 kpc for the Sagittarius remnant) covers the whole host, though its density
        never reaches one system in a cell.
        """
        one_cell = self.edge_pc ** 3
        reach = np.zeros(len(self._reach))
        for i, table in enumerate(self._tables):
            if total_systems[i] < 1.0 or self._central[i] * one_cell < 1.0:
                continue
            needed = 1.0 / (self._central[i] * one_cell)
            x = float(np.interp(needed, table.density[::-1], table.x[::-1]))
            reach[i] = min(self._reach[i], x * self._scale_pc[i] + 0.87 * self.edge_pc)
        return reach

    def qualifying_volume_pc3(self):
        """
        The volume, cubic parsecs, of the sectors the components qualify on their own (GEN.627): each component's sphere
        out to where its density reaches one system in a cell, at least one cell; the sum, so overlap counts twice.
        """
        radius = np.maximum(self._density_reach[self._density_reach > 0.0] - 0.87 * self.edge_pc, 0.0)
        return float(np.maximum(4.0 / 3.0 * math.pi * radius ** 3, self.edge_pc ** 3).sum())

    def _group_by_reach(self):
        """
        Splits the components into the few with a wide reach (globular clusters, satellites), which `bound` and `high_band`
        test one by one, and the many of small reach (the young open clusters of GEN.244, tens of thousands in a galaxy),
        which they find by a window in a sorted radius instead: a component of reach under
        `tuning.CLUSTER_SMALL_REACH_PC` can touch a ring or a shell only if its radius is within its reach of it.
        """
        reach = self._reach
        small = reach <= tuning.CLUSTER_SMALL_REACH_PC
        self._big = np.nonzero(~small)[0]
        self._small = np.nonzero(small)[0]
        cylindrical = np.hypot(self._centers[:, 0], self._centers[:, 1])[self._small]
        distance = np.linalg.norm(self._centers, axis=1)[self._small]
        order = np.argsort(cylindrical, kind="stable")
        self._small_by_cyl, self._small_cyl = self._small[order], cylindrical[order]
        order = np.argsort(distance, kind="stable")
        self._small_by_dist, self._small_dist = self._small[order], distance[order]
        slack = tuning.CLUSTER_BOUND_SLACK_EDGES * self.edge_pc
        self._small_window = (float(reach[self._small].max()) if len(self._small) else 0.0) + slack
        density_reach = self._density_reach
        self._density_window = (float(density_reach[self._small].max()) if len(self._small) else 0.0) + slack
        qualifying = self._small[density_reach[self._small] > 0.0]
        self._encloses_origin = bool(np.any(
            np.linalg.norm(self._centers[qualifying], axis=1) <= density_reach[qualifying] + slack)) if len(qualifying) else False

    def _grid_cell(self, point):
        grid = tuning.CLUSTER_GRID_PC
        return (math.floor(point[0] / grid), math.floor(point[1] / grid), math.floor(point[2] / grid))

    def _build_index(self):
        """Grid cell -> the components whose tidal sphere, grown by the cell cube's reach, touches it."""
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
        if len(near) > _PREFILTER_ABOVE:
            # Hundreds of young clusters share a grid cell: one vector test drops those too far to reach the cube.
            near = np.asarray(near)
            apart = np.linalg.norm(self._centers[near] - position, axis=1)
            near = near[apart - 0.87 * self.edge_pc < self._reach[near]]
        return [(int(i), self._cube_systems(int(i), position)) for i in near]

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

    def fe_h_peaks(self, position_pc, population):
        """
        The [Fe/H] distribution of the stars of `population` at `position_pc` (a family's populations share one:
        a satellite's three, a young cluster's twelve), `((weight, mean, sigma), ...)`: one peak per
        component of the family that reaches the point, weighted by the systems it puts in the cell and centred on its own
        [Fe/H] with `tuning.COMPONENT_FE_H_SIGMA[family]`. For a young cluster, only the clusters whose age is that
        population's. Empty where no component of the family reaches.
        """
        family = FAMILY_OF_POPULATION.get(population, population)
        sigma = tuning.COMPONENT_FE_H_SIGMA[family]
        exact = population in OPEN_CLUSTER_POPULATIONS
        return tuple((systems, self.fe_h[i], sigma) for i, systems in self._contributions(position_pc)
                     if self.families[i] == family and systems > 0.0
                     and (not exact or population in self.shares[i]))

    def bound(self, r_cyl, z):
        """
        An upper bound of `relative_density` over every point at cylindrical radius `r_cyl` and height `z`, whatever
        its angle: each component at the nearest the ring comes to it, less the cell's reach (`CLUSTER_BOUND_SLACK_EDGES`
        edges), the components added. `0.0` for a field without components.
        """
        if not len(self):
            return 0.0
        slack = tuning.CLUSTER_BOUND_SLACK_EDGES * self.edge_pc
        near = self._big
        if len(self._small):
            low, high = np.searchsorted(self._small_cyl, (r_cyl - self._small_window, r_cyl + self._small_window))
            near = np.concatenate((self._big, self._small_by_cyl[low:high]))
        if not len(near):
            return 0.0
        centers = self._centers[near]
        nearest = np.hypot(r_cyl - np.hypot(centers[:, 0], centers[:, 1]), z - centers[:, 2]) - slack
        total = 0.0
        for k in np.nonzero(nearest < self._reach[near])[0]:
            total += self._systems_per_pc3(int(near[k]), max(float(nearest[k]), 0.0)) * self.edge_pc ** 3
        return total / self.expected_at_density_1

    def high_band(self, shell, band_at, edge_pc):
        """
        The highest `|band|` of shell `shell` any component reaches (grown by the bound's slack), or `-1` when none
        does: a component far from the plane is no part of the run of bands the disk and bulge qualify, so the outline
        takes it in by hand. A component reaches as far as its density does (`_qualifying_reach`), not its tidal radius,
        and one whose density never reaches one system in a cell adds nothing (GEN.627). `band_at` is
        `galaxy.spherical.band_at`.
        """
        slack = tuning.CLUSTER_BOUND_SLACK_EDGES * edge_pc
        radius = shell * edge_pc
        if self._encloses_origin:
            return band_at(shell, 0.0) if shell else 0   # the sphere holds the pole: the shell's top band
        high = -1
        for i in self._big:
            center, reach = self._centers[i], self._density_reach[i]
            if reach <= 0.0:
                continue
            distance = float(np.linalg.norm(center))
            grown = reach + slack
            if distance <= grown:
                return band_at(shell, 0.0) if shell else 0
            if abs(distance - radius) > grown:
                continue
            latitude = math.asin(max(-1.0, min(1.0, center[2] / distance)))
            half_angle = math.asin(min(1.0, grown / distance))
            highest = min(math.pi / 2.0, abs(latitude) + half_angle)
            high = max(high, abs(band_at(shell, math.pi / 2.0 - highest)))
        if len(self._small):
            low, top = np.searchsorted(self._small_dist, (radius - self._density_window, radius + self._density_window))
            near = self._small_by_dist[low:top]
            near = near[self._density_reach[near] > 0.0]
            if len(near):
                centers = self._centers[near]
                distance = np.linalg.norm(centers, axis=1)
                grown = self._density_reach[near] + slack
                inside = np.abs(distance - radius) <= grown
                if inside.any():
                    latitude = np.arcsin(np.clip(centers[inside, 2] / distance[inside], -1.0, 1.0))
                    half_angle = np.arcsin(np.minimum(1.0, grown[inside] / distance[inside]))
                    highest = float(np.minimum(math.pi / 2.0, np.abs(latitude) + half_angle).max())
                    high = max(high, abs(band_at(shell, math.pi / 2.0 - highest)))
        return high
