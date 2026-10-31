# planetgen/galaxy/influence.py

"""
The influence set of an object: the masses that bend its path (GEN.341)
=======================================================================

The orbital update changes an object's motion vector from the point masses around it and the galaxy's smooth pull
(`docs/design/orbits/orbital-updates.md` section 4; Boss's rule of 2026-10-18, card F1: "the largest Hill radius present in the
galaxy becomes the search radius for gravitational objects"). This module is the search and nothing else: it moves no
object and reads no database.

- **The radius** is one number for the whole galaxy, `largest_hill_radius_pc`: the galactic Hill radius of the most
  extended object, `r (m / 3 (M(<r) + m))^(1/3)` with `M(<r)` the galaxy's own mass inside the object's radius
  (`galaxy.enclosed_mass`, GEN.370). It is worked out from each sector's heaviest mass (`sectors.max_mass_kg`, DB.55), as
  a sector's most extended object is its heaviest, and the central black hole is left out of it (it is in the smooth
  potential, GEN.285, and counted again it would stretch the radius to hundreds of parsecs).
- **The sectors** are those that can hold part of the sphere, found by sector bounds (`spherical.cells_touching_sphere`:
  the radius widened by the cells' reach, 1.3 edges, the margin without which a fifth of the neighbours is missed), not
  by a ring-and-slot stencil.
- **The members** are the point masses of those sectors within the radius of the object (`SectorMasses`, what a sector's
  point-mass table holds: a hierarchy is one point mass at its barycentre, a star cluster is its members, a lone rogue
  has none within reach, so its set is empty). The caller reads a sector through `read_sector(address)`.
- **The cap**: the `max_exact_members` nearest (default 64) are returned one by one, and the rest of those inside the
  radius are folded into one point mass per sector, at that part's centre of mass, so the work stays bounded in a cluster
  or beside a heavy object. `InfluenceSet.truncated` says so.

The order is fixed: sectors are read in address order and members within a sector by key, then everything is sorted by
distance with a stable sort, so the same catalogue gives the same set whatever order its sectors were stored or its
members listed in.
"""

import dataclasses
import math

import numpy as np

from planetgen.galaxy import enclosed_mass, spherical
from planetgen.physics import constants
from planetgen.physics import unit_registry as units

_msun_to_kg = units.converter("Msun", "kg")

_ly_to_m = units.converter("ly", "m")
_m_to_pc = units.converter("m", "pc")

DEFAULT_MAX_EXACT_MEMBERS = 64
"""int: The nearest members returned one by one; the rest are one point mass per sector (Boss, 2026-10-18, card F1)."""


@dataclasses.dataclass(frozen=True)
class SectorMasses:
    """
    The point masses one sector holds: `address` (shell, band, wedge), `positions_pc` (`n` by 3, galaxy frame),
    `masses_kg` (`n`) and `keys` (`n` hashable, mutually orderable names, say object IDs). Built through `of`, which puts
    the members in key order.
    """

    address: tuple
    positions_pc: np.ndarray
    masses_kg: np.ndarray
    keys: tuple

    @staticmethod
    def of(address, positions_pc, masses_kg, keys):
        """The `SectorMasses` of a sector's members in any order."""
        positions = np.asarray(positions_pc, dtype=float).reshape(-1, 3)
        masses = np.asarray(masses_kg, dtype=float).reshape(-1)
        keys = tuple(keys)
        if not len(positions) == len(masses) == len(keys):
            raise ValueError("positions, masses and keys must be the same length")
        order = sorted(range(len(keys)), key=lambda i: keys[i])
        return SectorMasses(tuple(address), positions[order], masses[order], tuple(keys[i] for i in order))


@dataclasses.dataclass(frozen=True)
class Monopole:
    """The members of one sector beyond the exact ones, as one point mass at their centre of mass."""

    address: tuple
    position_pc: tuple
    mass_kg: float
    count: int


@dataclasses.dataclass(frozen=True)
class InfluenceSet:
    """
    What bends one object's path: the `max_exact_members` nearest masses within `radius_pc` (`keys`, `positions_pc`,
    `masses_kg`, nearest first) and, for the rest inside the radius, one `Monopole` per sector (address order).
    `sectors_read` counts the sectors the search read; `truncated` is true when any member was folded into a monopole.
    """

    radius_pc: float
    keys: tuple
    positions_pc: np.ndarray
    masses_kg: np.ndarray
    monopoles: tuple
    sectors_read: int

    @property
    def truncated(self):
        """Whether any member inside the radius was folded into a monopole."""
        return bool(self.monopoles)

    @property
    def member_count(self):
        """Every mass inside the radius, the exact ones and the folded ones."""
        return len(self.keys) + sum(m.count for m in self.monopoles)

    @property
    def total_mass_kg(self):
        """The mass inside the radius, exact and folded."""
        return float(self.masses_kg.sum()) + sum(m.mass_kg for m in self.monopoles)


def galactic_hill_radius_pc(masses_kg, radii_pc, shape, galaxy_seed):
    """
    The galactic Hill radius, parsecs, of objects of `masses_kg` at galactocentric `radii_pc` (arrays or numbers) in the
    galaxy of `shape` and `galaxy_seed`: `r (m / 3 (M(<r) + m))^(1/3)`, with `M(<r)` the mass the galaxy holds inside
    `r` (`enclosed_mass`) and `r` no smaller than the nearest a system's Hill sphere is worked out at
    (`constants.MIN_TIDAL_DISTANCE_LY`), as for a star's own `Star.calculate_system_perimeter`.
    """
    floor_pc = _m_to_pc(_ly_to_m(constants.MIN_TIDAL_DISTANCE_LY))
    radii = np.maximum(np.asarray(radii_pc, dtype=float), floor_pc)
    masses = np.asarray(masses_kg, dtype=float)
    central = _msun_to_kg(enclosed_mass.mass_model(shape, galaxy_seed).enclosed_solar(radii))
    return radii * np.cbrt(masses / (3.0 * (central + masses)))


def largest_hill_radius_pc(masses_kg, radii_pc, shape, galaxy_seed, skip=None):
    """
    The one search radius for the whole galaxy: the largest galactic Hill radius over the given masses, parsecs.

    Args:
        masses_kg (array): The heaviest mass of each sector (`sectors.max_mass_kg`), or every object's.
        radii_pc (array): Each one's distance from the galaxy's centre (a sector's centre is near enough: a sector is
            4 pc across and the radius moves the enclosed mass slowly).
        shape (GalaxyShape): The galaxy's stored shape.
        galaxy_seed (bytes or None): Its seed.
        skip (array of bool, optional): True for the entries left out of the maximum: the central black hole, which the
            smooth potential holds (the caller names it; it is the object at the centre).

    Returns:
        float: The radius, or `0.0` when nothing is left to count.
    """
    masses = np.asarray(masses_kg, dtype=float)
    radii = np.asarray(radii_pc, dtype=float)
    if skip is not None:
        keep = ~np.asarray(skip, dtype=bool)
        masses, radii = masses[keep], radii[keep]
    if masses.size == 0:
        return 0.0
    return float(np.max(galactic_hill_radius_pc(masses, radii, shape, galaxy_seed)))


def influence_set(center_pc, radius_pc, edge_pc, read_sector, exclude_keys=(),
                  max_exact_members=DEFAULT_MAX_EXACT_MEMBERS, reach_edges=spherical.CELL_REACH_EDGES):
    """
    The masses inside `radius_pc` of `center_pc`, read from the sectors that can hold part of that sphere.

    Args:
        center_pc (tuple): The object's `(x, y, z)`, galaxy frame, parsecs.
        radius_pc (float): The search radius (`largest_hill_radius_pc`).
        edge_pc (float): The sector edge, parsecs.
        read_sector (callable): `read_sector(address)` returns the sector's `SectorMasses`, or `None` for a sector
            holding no mass; called once per sector, in address order.
        exclude_keys (iterable): Keys left out: the object itself and anything it is part of (its own system's members
            have their own part in the update).
        max_exact_members (int): The nearest members kept one by one; the rest, per sector, become one `Monopole`.
        reach_edges (float): The margin the sector search is widened by, in edges. The default holds every cell's
            farthest point (`spherical.CELL_REACH_EDGES`); a test sets 0 to show what the margin is for.

    Returns:
        InfluenceSet: The nearest members and the folded rest.

    Raises:
        ValueError: For a negative radius or margin, an edge that is not positive, or `max_exact_members` below 1.
    """
    if radius_pc < 0 or reach_edges < 0:
        raise ValueError(f"radius_pc and reach_edges must be >= 0, got {radius_pc} and {reach_edges}")
    if edge_pc <= 0:
        raise ValueError(f"edge_pc must be > 0, got {edge_pc}")
    if max_exact_members < 1:
        raise ValueError(f"max_exact_members must be >= 1, got {max_exact_members}")
    excluded = frozenset(exclude_keys)
    center = np.asarray(center_pc, dtype=float)
    cells = sorted(spherical.cells_touching_sphere(tuple(center), radius_pc, edge_pc, reach_edges), key=lambda cell: cell[:3])
    keys, positions, masses, sectors, distances = [], [], [], [], []
    for cell in cells:
        sector = read_sector(cell[:3])
        if sector is None or not len(sector.keys):
            continue
        distance = np.sqrt(((sector.positions_pc - center) ** 2).sum(axis=1))
        for i in np.flatnonzero(distance <= radius_pc):
            if sector.keys[i] in excluded:
                continue
            keys.append(sector.keys[i])
            positions.append(sector.positions_pc[i])
            masses.append(sector.masses_kg[i])
            sectors.append(sector.address)
            distances.append(distance[i])
    if not keys:
        return InfluenceSet(float(radius_pc), (), np.zeros((0, 3)), np.zeros(0), (), len(cells))
    order = np.argsort(np.asarray(distances), kind="stable")
    positions = np.asarray(positions)[order]
    masses = np.asarray(masses)[order]
    keys = [keys[i] for i in order]
    sectors = [sectors[i] for i in order]
    exact = slice(0, max_exact_members)
    return InfluenceSet(float(radius_pc), tuple(keys[exact]), positions[exact], masses[exact],
                        _monopoles(sectors[max_exact_members:], positions[max_exact_members:],
                                   masses[max_exact_members:]), len(cells))


def _monopoles(sectors, positions, masses):
    """One `Monopole` per sector for the members beyond the exact ones, in address order."""
    folded = {}
    for sector, position, mass in zip(sectors, positions, masses):
        folded.setdefault(sector, []).append((position, mass))
    result = []
    for address in sorted(folded):
        members = folded[address]
        mass = math.fsum(m for _p, m in members)
        weighted = np.array([p * m for p, m in members]).sum(axis=0) / mass
        result.append(Monopole(address, tuple(float(c) for c in weighted), mass, len(members)))
    return tuple(result)
