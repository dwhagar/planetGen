# planetgen/galaxy/core_gas.py

"""
The gas of a galaxy's centre and the young stars it forms (GEN.282; docs/design/galaxy/realistic-galactic-core.md 3.4 and 6).

A bar carries the gas of its gap inward and piles it into a ring near the centre, the central molecular zone, which forms
stars at a small share of what comes in. Nothing here is typed in for a galaxy:

1. **The inflow.** The disk's gas inside the bar's gap (`galaxy.disk_history`, GEN.469: each ring's gas at the end of
   its history, inside `galaxy.star_formation.gap_radius_pc`) leaves on `tuning.CORE_GAS_INFLOW_TIME_GY`. The Milky
   Way's model gives 0.7 solar masses a year against the 0.8 measured (Hatchfield et al. 2021, in Henshaw et al. 2023).
   A galaxy without a bar has no inflow, so no ring.
2. **The rate.** The ring forms stars at `tuning.CORE_STAR_FORMATION_PER_INFLOW` of the inflow (0.07 over 0.8 in the Milky
   Way: the review's continuity equation, inflow = star formation + outflow + growth).
3. **The ring.** Radius, width and height are the bulge's scale radius times ratios (`tuning.CORE_RING_*`), as the
   nuclear disc's are. The rate over its area is a star formation rate density; the cluster formation efficiency of that
   (`physics.open_cluster_evolution.cluster_formation_efficiency`, Kruijssen 2012) is the share of the new stars that
   form in bound clusters, which `generation/open_clusters.py` draws as components; the rest are the field's.
4. **The field's young stars.** The rate times the age of the `young` population (0.1 Gy) times the unbound share is
   the mass of field stars the ring holds younger than that. They are part of the nuclear disc, not on top of it: the
   density model takes them out of the `nuclear_disc` population and puts them in `young`, so the total density, the
   skeleton's bound and the map prisms do not change, while the places the profile calls young, with their supernova
   remnants, nebulae and massive stars, follow the ring.
"""

import dataclasses
import functools
import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import disk_history, star_formation
from planetgen.physics import open_cluster_evolution

_GY_YR = 1.0e9
_PC2_PER_KPC2 = 1.0e6


@dataclasses.dataclass(frozen=True)
class CoreGas:
    """
    A galaxy's central ring as `core_gas` works it out.

    Attributes:
        inflow_msun_yr (float): The gas the bar carries to the ring, solar masses per year.
        rate_msun_yr (float): The ring's star formation rate, solar masses per year.
        radius_pc, sigma_pc, height_pc (float): The ring's radius, the standard deviation of its Gaussian width and its
            `sech^2` scale height.
        rate_density (float): The rate over the ring's area, solar masses per year per kpc2.
        bound_share (float): The share of the new stars that form in bound clusters.
        young_density (float): The amplitude of the field's young stars at the ring, in the model's unnormalized
            units (like `nuclear_disc_density`).
        gas_surface_density_msun_pc2 (float): The ring's gas over its area: the inflow times
            `tuning.CORE_RING_GAS_RESIDENCE_MYR`, over `2 pi R FWHM` (GEN.1195; the most massive cluster the ring forms
            reads it, `physics.cluster_maximum_mass`).
    """
    inflow_msun_yr: float
    rate_msun_yr: float
    radius_pc: float
    sigma_pc: float
    height_pc: float
    rate_density: float
    bound_share: float
    young_density: float
    gas_surface_density_msun_pc2: float = 0.0


NO_RING = CoreGas(0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0)
"""CoreGas: What a galaxy without a bar or without a nucleus has."""


def _radial_integral(radius_pc, sigma_pc):
    """`integral_0^inf r exp(-(r - R)^2 / (2 sigma^2)) dr`, the ring's radial weight for its area element."""
    z = radius_pc / (sigma_pc * math.sqrt(2.0))
    return (sigma_pc ** 2 * math.exp(-z * z) + radius_pc * sigma_pc * math.sqrt(math.pi / 2.0) * (1.0 + math.erf(z)))


@functools.lru_cache(maxsize=64)
def core_gas(shape):
    """
    The central ring of the galaxy of this shape (see the module docstring); `NO_RING` for one without a bar.

    Args:
        shape (GalaxyShape): The galaxy's shape.

    Returns:
        CoreGas
    """
    if not shape.structure.bar_present:
        return NO_RING
    from planetgen.galaxy import density
    found = disk_history.history_of(shape)
    radii = found.radius_scale_lengths * shape.disk_scale_length_pc
    width = float(radii[1] - radii[0]) if len(radii) > 1 else float(radii[0])
    area_pc2 = 2.0 * math.pi * radii * width
    inside = radii < star_formation.gap_radius_pc(shape)
    gap_gas = float((found.gas_now[inside] * area_pc2[inside]).sum())
    inflow = gap_gas / (tuning.CORE_GAS_INFLOW_TIME_GY * _GY_YR)
    rate = tuning.CORE_STAR_FORMATION_PER_INFLOW * inflow
    if rate <= 0.0:
        return NO_RING
    bulge = shape.bulge_scale_radius_pc
    radius = tuning.CORE_RING_RADIUS_BULGE_RATIO * bulge
    fwhm = tuning.CORE_RING_WIDTH_BULGE_RATIO * bulge
    height = tuning.CORE_RING_HEIGHT_BULGE_RATIO * bulge
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    ring_area_kpc2 = 2.0 * math.pi * radius * fwhm / _PC2_PER_KPC2
    rate_density = rate / ring_area_kpc2
    bound = float(open_cluster_evolution.cluster_formation_efficiency(rate_density))
    young_window_yr = tuning.STELLAR_POPULATION_AGE_RANGES_GY["young"][1] * _GY_YR
    field_young_mass = rate * young_window_yr * (1.0 - bound)
    stellar_mass, _halo = disk_history.stellar_and_halo_masses(shape)
    model_units_per_solar = sum(density.integrated_masses(shape).values()) / stellar_mass
    amplitude = (field_young_mass * model_units_per_solar
                 / (2.0 * math.pi * _radial_integral(radius, sigma) * 4.0 * height))
    return CoreGas(inflow_msun_yr=inflow, rate_msun_yr=rate, radius_pc=radius, sigma_pc=sigma, height_pc=height,
                   rate_density=rate_density, bound_share=bound, young_density=amplitude,
                   gas_surface_density_msun_pc2=(inflow * tuning.CORE_RING_GAS_RESIDENCE_MYR * 1.0e6
                                                 / (ring_area_kpc2 * _PC2_PER_KPC2)))


def _sech_squared(value):
    return 1.0 / math.cosh(value) ** 2 if abs(value) < 350.0 else 0.0


def young_density(r_cyl, z, core):
    """The field's young stars of the ring at cylindrical radius `r_cyl` and height `z`, unnormalized (`0.0` for no ring)."""
    if not core.young_density:
        return 0.0
    return (core.young_density * math.exp(-0.5 * ((r_cyl - core.radius_pc) / core.sigma_pc) ** 2)
            * _sech_squared(z / (2.0 * core.height_pc)))


def young_density_array(r_cyl, z, core):
    """`young_density` on arrays (`z` one height or an array that broadcasts against `r_cyl`)."""
    r_cyl = np.asarray(r_cyl, dtype=float)
    if not core.young_density:
        return np.zeros(np.broadcast(r_cyl, np.asarray(z)).shape)
    e = np.exp(-2.0 * np.minimum(np.abs(np.asarray(z, dtype=float)) / (2.0 * core.height_pc), 350.0))
    vertical = np.minimum(1.0, 4.0 * e / (1.0 + e) ** 2)
    return core.young_density * np.exp(-0.5 * ((r_cyl - core.radius_pc) / core.sigma_pc) ** 2) * vertical


def young_bound(inner, outer, z_near, core):
    """
    An upper bound on `young_density` in a region whose cylindrical radius runs from `inner` to `outer` and whose nearest
    `|z|` is `z_near` (arrays): the Gaussian ring's peak where the region holds the ring's radius, else its value at
    the nearer edge, times the vertical profile at the nearest height.
    """
    inner, outer, z_near = (np.asarray(value, dtype=float) for value in (inner, outer, z_near))
    if not core.young_density:
        return np.zeros(np.broadcast(inner, outer, z_near).shape)
    gap = np.maximum(np.maximum(inner - core.radius_pc, core.radius_pc - outer), 0.0)
    e = np.exp(-2.0 * np.minimum(z_near / (2.0 * core.height_pc), 350.0))
    vertical = np.minimum(1.0, 4.0 * e / (1.0 + e) ** 2)
    return core.young_density * np.exp(-0.5 * (gap / core.sigma_pc) ** 2) * vertical
