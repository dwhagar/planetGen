# planetgen/galaxy/unresolved_crossings.py

"""
How many unresolved keep-out spheres a straight course crosses (NAV.74; docs/design/galaxy/realistic-galactic-core.md section 11)
=================================================================================================================================

A crowded cell is filled only with its bright stars (`galaxy.resolution`, GEN.281), so the rest of its stars are never made and a
course cannot bend round them. Their pull is already in the smooth galactic field, and a course does not bend to a smooth field,
so what is left to say is how many of their keep-out spheres a course passes through. A star's keep-out is its galactic Hill radius
(`galaxy.keepout`), worked out here with the same rule, from the mean mass of the cell's unresolved systems and the galaxy's own
mass inside the point's radius (`galaxy.enclosed_mass`).

For a straight line through a cloud of spheres of cross-section `pi R^2` at `n` to the cubic parsec, the mean number crossed is
`n pi R^2 L`. Here `n` is the cell's unresolved systems over its volume, `R` the Hill radius at the point, and the line is summed
a short step at a time so the cells and the radius change along it.

Nothing is stored: the answer is a function of the segment and the galaxy's skeleton. A course that passes farther from the centre
than the capped region reaches (`capped_reach_pc`) costs one distance and reads no cell; a course that does touch it reads the cells
it passes through, each worked out once and kept (`UnresolvedField` keeps them for as long as the caller keeps it).
"""

import math
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np

from planetgen.galaxy import density, enclosed_mass, keepout, mass_chain, resolution, spherical
from planetgen.physics import unit_registry as units

_km_to_pc = units.converter("km", "pc")
_msun_to_kg = units.converter("Msun", "kg")

REACH_DIRECTIONS = 128
"""int: The directions (a Fibonacci sphere) the capped region's reach is searched along."""
REACH_MARGIN = 0.5
"""float: A place counts toward the reach when it expects at least this share of the ceiling, so the reach errs wide (a cell's volume is not an edge cube's)."""
REACH_RADIUS_STEP = 1.05
"""float: The ratio between neighbouring radii of the reach search, once the radii are wider apart than a sector edge."""
REACH_LIMIT_PC = 2.0e4
"""float: The farthest radius the reach search looks at, parsecs."""
STEPS_PER_EDGE = 4
"""int: The steps a segment is summed in for each sector edge it spans."""
NEAR_CENTER_STEP = 0.125
"""float: The step is also no longer than this share of the distance from the center, where the density and the Hill radius change fastest."""
NEAR_CENTER_LEAST_DIVISOR = 32
"""int: ... but never shorter than a sector edge over this."""
INNER_RADIUS_PC = 0.05
"""float: The smallest radius a Hill radius is worked out at, parsecs (the galactic Hill radius is zero at the center; `enclosed_mass.INNER_RADIUS_PC`)."""


@dataclass(frozen=True)
class Crossings:
    """
    What a straight segment crosses.

    Attributes:
        crossings (float): The mean number of unresolved keep-out spheres the segment passes through.
        length_pc (float): The length of the segment inside capped cells, parsecs.
        unresolved_systems (float): The unresolved systems of the capped cells the segment passes through, each cell counted once.
        cells (int): The capped cells the segment passes through.
    """
    crossings: float
    length_pc: float
    unresolved_systems: float
    cells: int


NONE_CROSSED = Crossings(0.0, 0.0, 0.0, 0)
"""Crossings: A segment that touches no capped cell."""


def _fibonacci_directions(count):
    """`count` unit vectors spread over the sphere, as three arrays."""
    index = np.arange(count) + 0.5
    z = 1.0 - 2.0 * index / count
    phi = math.pi * (1.0 + math.sqrt(5.0)) * index
    ring = np.sqrt(1.0 - z * z)
    return ring * np.cos(phi), ring * np.sin(phi), z


def capped_reach_pc(skeleton):
    """
    The distance from the center beyond which no cell of the galaxy expects more systems than the cell star ceiling, so a course that stays
    that far out has nothing unresolved to count. Found by looking at the density model along `REACH_DIRECTIONS` directions at radii a sector
    edge apart (then `REACH_RADIUS_STEP` apart), vectorised; the farthest place that expects `REACH_MARGIN` of the ceiling sets it. A very thin
    feature between the directions could be missed (the core terms are round or discs, which the directions resolve). Zero when no cell is capped.
    Worked out once per galaxy: it is the one costly part of a query, a few thousand points.
    """
    edge = skeleton.edge_pc
    per_density_one = skeleton.expected_system_count_at_density_1
    threshold = REACH_MARGIN * skeleton.cell_star_ceiling
    ux, uy, uz = _fibonacci_directions(REACH_DIRECTIONS)
    radii, radius = [], 0.0
    while radius <= REACH_LIMIT_PC:
        radii.append(radius)
        radius = radius + edge if radius < 20.0 * edge else radius * REACH_RADIUS_STEP
    farthest = 0.0
    for radius in radii:
        values = density.total_density_array(radius * ux, radius * uy, radius * uz, skeleton.shape)
        if float(np.max(values)) * per_density_one >= threshold:
            farthest = radius
    return farthest + 2.0 * edge if farthest > 0.0 else 0.0


def _clip_to_ball(start, end, radius):
    """The parameters `(t0, t1)` of the part of the segment inside the ball of `radius` about the origin, or `None`."""
    direction = [e - s for s, e in zip(start, end)]
    a = sum(d * d for d in direction)
    if a <= 0.0:
        return None
    b = 2.0 * sum(s * d for s, d in zip(start, direction))
    c = sum(s * s for s in start) - radius * radius
    discriminant = b * b - 4.0 * a * c
    if discriminant < 0.0:
        return None
    root = math.sqrt(discriminant)
    t0, t1 = max(0.0, (-b - root) / (2.0 * a)), min(1.0, (-b + root) / (2.0 * a))
    return (t0, t1) if t1 > t0 else None


class UnresolvedField:
    """
    The unresolved stars of one galaxy, as far as a course needs them.

    Args:
        skeleton (store.GalaxySkeletonInfo): The galaxy's stored skeleton.
        mass_limit_at (callable, optional): `(shell, band) -> mass limit or None`, the scatter's mass limit
            (`ScatterLimits.mass_limit_at`) the fill resolved a cell under; `None` when no scatter ran.
        reach_pc (float, optional): The capped region's reach; worked out by `capped_reach_pc` on first use when not given.
    """

    def __init__(self, skeleton, mass_limit_at=None, reach_pc=None):
        self.skeleton = skeleton
        self._mass_limit_at = mass_limit_at or (lambda shell, band: None)
        self._reach = reach_pc
        self._cells = {}
        self._seed = skeleton.galaxy_seed or mass_chain.FALLBACK_SEED

    @property
    def reach_pc(self):
        """The capped region's reach, parsecs."""
        if self._reach is None:
            self._reach = capped_reach_pc(self.skeleton)
        return self._reach

    def cell(self, address):
        """`(unresolved_systems, volume_pc3, mean_unresolved_mass_kg)` of the capped cell `address`, or `None` when it is filled whole."""
        if address not in self._cells:
            skeleton = self.skeleton
            shell, band, wedge = address
            position = spherical.sector_position_pc(shell, band, wedge, skeleton.edge_pc)
            volume = spherical.cell_volume(shell, band, wedge, skeleton.edge_pc)
            mix = {name: max(value, 0.0) for name, value in density.population_densities(position, skeleton.shape).items()}
            light = SimpleNamespace(
                expected_systems=skeleton.expected_system_count_at_density_1 * max(density.relative_density(position, skeleton.shape), 0.0),
                mix=mix, volume_pc3=volume, edge_pc=skeleton.edge_pc, shape=skeleton.shape)
            found = resolution.cell_resolution(light, skeleton.cell_star_ceiling, skeleton.fill_luminosity_cap_sol,
                                               self._mass_limit_at(shell, band))
            if not found.capped or found.unresolved_systems <= 0.0:
                self._cells[address] = None
            else:
                mean_kg = _msun_to_kg(found.unresolved_mass_solar / found.unresolved_systems)
                self._cells[address] = (found.unresolved_systems, volume, mean_kg)
        return self._cells[address]

    def hill_radius_pc(self, mean_mass_kg, radius_pc):
        """A star's keep-out radius (its galactic Hill radius, `keepout`) at `radius_pc` from the center, parsecs."""
        radius_pc = max(radius_pc, INNER_RADIUS_PC)
        central = enclosed_mass.enclosed_mass_kg(self.skeleton.shape, self._seed, radius_pc)
        return _km_to_pc(keepout.galactic_hill_radius_km(mean_mass_kg, radius_pc, central))

    def crossings(self, start_pc, end_pc, step_pc=None):
        """
        The mean number of unresolved keep-out spheres the straight segment `start_pc` to `end_pc` (galaxy-frame parsecs) passes through,
        with the length it spends in capped cells and the unresolved systems of those cells (`Crossings`). `step_pc` is the step it is summed
        in (a quarter of a sector edge when not given). The same segment in either direction gives the same answer.
        """
        start, end = tuple(float(c) for c in start_pc), tuple(float(c) for c in end_pc)
        clipped = _clip_to_ball(start, end, self.reach_pc) if self.reach_pc > 0.0 else None
        if clipped is None:
            return NONE_CROSSED
        t0, t1 = clipped
        edge = self.skeleton.edge_pc
        step = edge / STEPS_PER_EDGE if step_pc is None else float(step_pc)
        span = math.dist(start, end)
        total = inside = 0.0
        seen = set()
        t = t0
        while t < t1 - 1e-12:
            here = tuple(s + (e - s) * t for s, e in zip(start, end))
            near = max(math.sqrt(sum(c * c for c in here)) * NEAR_CENTER_STEP, edge / NEAR_CENTER_LEAST_DIVISOR)
            ds = min(step, near, (t1 - t) * span)
            middle = tuple(s + (e - s) * (t + 0.5 * ds / span) for s, e in zip(start, end))
            t += ds / span
            address = spherical.sector_address_at(middle, edge)
            found = self.cell(address)
            if found is None:
                continue
            systems, volume, mean_kg = found
            seen.add(address)
            hill = self.hill_radius_pc(mean_kg, math.sqrt(sum(c * c for c in middle)))
            total += systems / volume * math.pi * hill * hill * ds
            inside += ds
        return Crossings(total, inside, sum(self.cell(a)[0] for a in seen), len(seen))
