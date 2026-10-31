# planetgen/galaxy/enclosed_mass.py

"""
The mass a galaxy holds inside a radius, for the galaxy's tides (GEN.370)
==========================================================================

A Hill sphere is `r (m / 3 M)^(1/3)`; the `M` is the mass of the galaxy inside the object's own radius `r`, from the
galaxy's own density model, not the whole halo's. The whole halo (`mass_chain`'s virial mass) is one number for every
radius and overstates what pulls on a star that sits well inside it (most of the halo lies beyond), so a star near the
middle got a Hill sphere several times too small.

`enclosed_mass_solar(shape, galaxy_seed, radius_pc)` adds three parts, all from the one galaxy:

1. **The stars.** The smooth density model (`galaxy.density.total_density_array`: thin and thick disk, bar bulge,
   stellar halo, nuclear cluster and disc, the arms averaged round the circle) integrated over the cells of a grid in
   cylindrical radius and height, each cell's mass given to its radius from the centre, then summed in order of
   radius. The model's units become solar masses through `model_terms["mass_per_solar"]`. The star clusters and
   satellites the shape's `components` add are left out: they are stars of the smooth model already, or lie in the
   halo where their mass is small beside the dark halo's.
2. **The central black hole** (`density.central_black_hole_mass_solar`), a point at the centre.
3. **The dark halo.** An NFW halo (Navarro, Frenk and White 1996) of the mass chain's halo mass less its stars, with the
   chain's own virial radius and concentration, so the halo's share inside `r` is
   `M_dark (ln(1 + x) - x / (1 + x)) / (ln(1 + c) - c / (1 + c))`, `x = r c / r200`, held at `M_dark` beyond `r200`.

The gas is left out (a stand-in: GEN.338's potential gives the circular speed and the enclosed mass from all the
components and replaces this module's total). The grid is built once per galaxy and cached.
"""

import collections
import functools
import math

import numpy as np

from planetgen.galaxy import density, disk_history, mass_chain, seed as galaxySeed, shape_draw
from planetgen.physics import unit_registry as units

_msun_to_kg = units.converter("Msun", "kg")

INNER_RADIUS_PC = 0.05
"""float: The radius the grid starts at, parsecs; inside it the stars' mass falls as the cube of the radius."""
RADIAL_CELLS = 200
HEIGHT_CELLS = 200
AZIMUTH_STEPS = 16
"""int: Points round each circle (the arms, bar and boxy bulge carry harmonics of up to 4 round it)."""
NO_HALO_OUTER_PC = 100_000.0
"""float: The outer radius of the grid for a galaxy too light for a halo (no virial radius), parsecs."""
TABLE_POINTS = 512
"""int: Points of the table `enclosed_mass_solar` interpolates, evenly spaced in the log of the radius."""


def _stellar_table(shape):
    """
    `(log10 radius_pc, enclosed solar masses)` arrays of the stars inside each radius of the model `shape` describes.
    """
    terms = density.model_terms(shape)
    stellar_mass, halo_mass = disk_history.stellar_and_halo_masses(shape)
    outer = (NO_HALO_OUTER_PC if halo_mass is None else mass_chain.virial_radius_kpc(halo_mass) * 1000.0)
    radial_edges = np.geomspace(INNER_RADIUS_PC, outer, RADIAL_CELLS + 1)
    height_edges = np.concatenate(([0.0], np.geomspace(INNER_RADIUS_PC, outer, HEIGHT_CELLS)))
    # The area-weighted middle of each ring and the middle of each height step.
    ring = np.sqrt((radial_edges[:-1] ** 2 + radial_edges[1:] ** 2) / 2.0)
    height = (height_edges[:-1] + height_edges[1:]) / 2.0
    volume = (math.pi * (radial_edges[1:] ** 2 - radial_edges[:-1] ** 2))[:, None] * (2.0 * np.diff(height_edges))[None, :]
    angles = (np.arange(AZIMUTH_STEPS) + 0.5) * (2.0 * math.pi / AZIMUTH_STEPS)
    x = ring[:, None, None] * np.cos(angles)[None, None, :]
    y = ring[:, None, None] * np.sin(angles)[None, None, :]
    z = height[None, :, None]
    # Per solar mass of the model's units: the density grid is `k_norm` times the raw density the masses are in.
    raw = density.total_density_array(x, y, z, shape).mean(axis=2) / shape.k_norm
    mass = raw * volume / terms["mass_per_solar"]
    radius = np.hypot(ring[:, None], height[None, :])
    order = np.argsort(radius, axis=None)
    cumulative = np.cumsum(mass.ravel()[order])
    sorted_radius = radius.ravel()[order]
    grid = np.linspace(math.log10(INNER_RADIUS_PC), math.log10(outer), TABLE_POINTS)
    enclosed = np.interp(grid, np.log10(sorted_radius), cumulative, left=0.0, right=float(cumulative[-1]))
    return grid, enclosed


DarkHalo = collections.namedtuple("DarkHalo", ("dark_solar", "concentration", "virial_radius_pc"))
"""The dark NFW halo of a galaxy: its mass (the halo's less its stars), concentration and virial radius in parsecs."""


def dark_halo(shape, galaxy_seed):
    """
    The `DarkHalo` of the galaxy of `shape`: the halo mass of the mass chain less the galaxy's stars, with the chain's own virial radius
    and the concentration the seed's scatter gives at the stored (possibly dial-pinned) masses, or the median one for `galaxy_seed` `None`.
    A galaxy too light for a halo has none (`dark_solar` 0, `concentration` and `virial_radius_pc` 1).
    """
    stellar_mass, halo_mass = disk_history.stellar_and_halo_masses(shape)
    if halo_mass is None:
        return DarkHalo(0.0, 1.0, 1.0)
    if galaxy_seed is None:
        concentration = mass_chain.median_concentration(halo_mass)
    else:
        concentration = mass_chain.draw_mass_chain(galaxy_seed, stellar_mass, halo_mass).concentration
    return DarkHalo(max(halo_mass - stellar_mass, 0.0), concentration, mass_chain.virial_radius_kpc(halo_mass) * 1000.0)


class MassModel:
    """The mass a galaxy holds inside a radius: its stars, its central black hole and its dark halo (module doc)."""

    def __init__(self, shape, galaxy_seed):
        self.log_radius, self.stellar_solar = _stellar_table(shape)
        self.black_hole_solar = density.central_black_hole_mass_solar(shape)
        self.dark_solar, self.concentration, self.virial_radius_pc = dark_halo(shape, galaxy_seed)

    def stars_solar(self, radius_pc):
        """The solar masses of stars inside `radius_pc`, from the cached table; below its first radius the mass falls as
        the cube of the radius."""
        radius = np.asarray(radius_pc, dtype=float)
        inner = 10.0 ** self.log_radius[0]
        inside = np.interp(np.log10(np.maximum(radius, 1e-300)), self.log_radius, self.stellar_solar)
        small = self.stellar_solar[0] * (radius / inner) ** 3
        return np.where(radius < inner, small, inside)

    def dark_halo_solar(self, radius_pc):
        """The solar masses of dark halo inside `radius_pc` (the NFW profile, held at its total beyond the virial
        radius)."""
        c = self.concentration
        x = np.minimum(np.asarray(radius_pc, dtype=float), self.virial_radius_pc) * c / self.virial_radius_pc
        return self.dark_solar * (np.log1p(x) - x / (1.0 + x)) / (math.log1p(c) - c / (1.0 + c))

    def enclosed_solar(self, radius_pc):
        """Everything inside `radius_pc`, solar masses: the stars, the central black hole and the dark halo."""
        return self.stars_solar(radius_pc) + self.black_hole_solar + self.dark_halo_solar(radius_pc)


@functools.lru_cache(maxsize=16)
def _model(shape, galaxy_seed):
    return MassModel(shape, galaxy_seed)


def mass_model(shape, galaxy_seed):
    """
    The `MassModel` of the galaxy of `shape` and `galaxy_seed`, built once and kept.

    Args:
        shape (GalaxyShape): The galaxy's shape; its star clusters and satellites (`components`) are not read.
        galaxy_seed (bytes or None): The galaxy's seed, for its drawn NFW concentration; `None` takes the median
            concentration of the halo mass.
    """
    return _model(shape._replace(components=None), None if galaxy_seed is None else bytes(galaxy_seed))


def enclosed_mass_solar(shape, galaxy_seed, radius_pc):
    """The mass, solar masses, inside `radius_pc` of the galaxy of `shape` and `galaxy_seed` (`MassModel.enclosed_solar`);
    a float for one radius."""
    result = mass_model(shape, galaxy_seed).enclosed_solar(radius_pc)
    return float(result) if np.ndim(result) == 0 else result


def enclosed_mass_kg(shape, galaxy_seed, radius_pc):
    """`enclosed_mass_solar` in kilograms."""
    return _msun_to_kg(enclosed_mass_solar(shape, galaxy_seed, radius_pc))


def current_enclosed_mass_kg(radius_pc):
    """
    The mass, kg, inside `radius_pc` of the galaxy the running unit belongs to: its stored shape and seed when
    `galaxy.seed.seeded` was given them, else the galaxy the seed draws with every shape option at its default (the
    fallback seed's outside any unit).
    """
    seed = galaxySeed.current_galaxy_seed() or mass_chain.FALLBACK_SEED
    shape = galaxySeed.current_galaxy_shape()
    if shape is None:
        shape = shape_draw.default_shape(seed)
    return enclosed_mass_kg(shape, seed, radius_pc)
