# planetgen/physics/open_cluster_evolution.py

"""
How a young open cluster forms, loses mass and sizes itself (GEN.244; docs/design/galaxy/intergalactic-and-halo-contents.md
5.14). Pure functions on numpy arrays, each from a paper whose number is in `tuning` ("a generated galaxy's young open
clusters"): `generation/open_clusters.py` draws the clusters with them.

* **Formation.** The share of the star formation that ends up in bound clusters, the cluster formation efficiency
  (`cluster_formation_efficiency`, Kruijssen 2012), and the initial masses, a power law of slope -2
  (`initial_masses`, Lamers et al. 2005).
* **Loss.** Stellar evolution (`mass_lost_by_evolution`, Lamers et al. 2005 Table 1) and tidal disruption on a timescale
  `t_dis = t0 (M / Msun) ** gamma` (`disruption_t0_yr`, `present_masses`), the mass a cluster holds at an age.
* **Size.** A half-mass radius from a half-mass density (`half_mass_radii_pc`, Portegies Zwart et al. 2010), grown by a
  factor of 5 to 10 over the first 100 Myr and held inside the cluster's Jacobi radius (`jacobi_radii_pc`).
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.physics import constants

_YEAR_PER_MYR = 1.0e6


def cluster_formation_efficiency(sigma_sfr_msun_yr_kpc2):
    """
    The share of star formation that ends up in bound clusters (Kruijssen 2012, equation 45, a fit to his model for a
    typical parameter set): `1 / (1.15 + 0.6 S ** -0.4 + 0.05 S ** -1)` for the star formation rate density `S` in
    solar masses per year per kpc2. Rises from 1% where the gas is thin to about 0.7 in the densest galaxies; 0 where
    there is no star formation.
    """
    sigma = np.maximum(np.asarray(sigma_sfr_msun_yr_kpc2, dtype=float), 0.0)
    a, b, power, c, minus_one = tuning.OPEN_CLUSTER_CFE_COEFFICIENTS
    with np.errstate(divide="ignore", invalid="ignore"):
        denominator = a + b * sigma ** power + c * sigma ** minus_one
    return np.where(sigma > 0.0, 1.0 / denominator, 0.0)


def mean_initial_mass_sol():
    """The mean initial mass of a cluster, solar masses, for the power-law mass function between
    `tuning.OPEN_CLUSTER_MASS_RANGE_SOL`: how many clusters a mass of bound star formation makes."""
    low, high = tuning.OPEN_CLUSTER_MASS_RANGE_SOL
    slope = tuning.OPEN_CLUSTER_MASS_FUNCTION_SLOPE
    if abs(slope - 2.0) < 1.0e-12:
        return math.log(high / low) / (1.0 / low - 1.0 / high)
    one = 1.0 - slope
    two = 2.0 - slope
    return ((high ** two - low ** two) / two) / ((high ** one - low ** one) / one)


def initial_masses(uniform, low_mass_sol=None):
    """Initial masses (solar masses) by inverse transform of the power-law mass function, from uniform draws in 0..1,
    between `low_mass_sol` (the range's low end when not given) and its high end."""
    low, high = tuning.OPEN_CLUSTER_MASS_RANGE_SOL
    if low_mass_sol is not None:
        low = low_mass_sol
    one = 1.0 - tuning.OPEN_CLUSTER_MASS_FUNCTION_SLOPE
    u = np.asarray(uniform, dtype=float)
    return (low ** one + u * (high ** one - low ** one)) ** (1.0 / one)


def low_mass_keeping_share(share):
    """The low end of the initial masses (solar masses) that leaves `share` of the clusters of the full range
    (`tuning.OPEN_CLUSTER_MASS_RANGE_SOL`): the number above a mass falls as `m ** (1 - slope)`, so it is the mass whose
    count is that share of the whole."""
    low, high = tuning.OPEN_CLUSTER_MASS_RANGE_SOL
    one = 1.0 - tuning.OPEN_CLUSTER_MASS_FUNCTION_SLOPE
    return (high ** one + share * (low ** one - high ** one)) ** (1.0 / one)


def mass_lost_by_evolution(age_yr, z):
    """
    The share of a cluster's initial mass that stellar evolution has returned to the gas by `age_yr` (Lamers et al.
    2005, equation 2): `log q = (log t - a) ** b + c`, with `(a, b, c)` interpolated between the metallicities of their
    Table 1 in `log Z`. Nothing is lost before `tuning.OPEN_CLUSTER_EVOLUTION_LOSS_START_MYR`.

    Args:
        age_yr (array): Ages, years.
        z (float or array): The mass fraction of metals (`tuning.OPEN_CLUSTER_SOLAR_Z` times 10 ** [Fe/H]).
    """
    rows = np.array(tuning.OPEN_CLUSTER_EVOLUTION_LOSS, dtype=float)
    log_z = np.log10(np.clip(np.asarray(z, dtype=float), rows[0, 0], rows[-1, 0]))
    a = np.interp(log_z, np.log10(rows[:, 0]), rows[:, 1])
    b = np.interp(log_z, np.log10(rows[:, 0]), rows[:, 2])
    c = np.interp(log_z, np.log10(rows[:, 0]), rows[:, 3])
    age = np.asarray(age_yr, dtype=float)
    start = tuning.OPEN_CLUSTER_EVOLUTION_LOSS_START_MYR * _YEAR_PER_MYR
    log_t = np.log10(np.maximum(age, start))
    lost = 10.0 ** (np.maximum(log_t - a, 0.0) ** b + c)
    return np.where(age > start, np.minimum(lost, 1.0), 0.0)


def disruption_t0_yr(rho_amb_msun_pc3, c_env_myr):
    """
    The constant `t0` of the disruption time `t_dis = t0 (M / Msun) ** gamma`, years (Lamers, Gieles and Portegies Zwart
    2005): `Cenv 1e-4 gamma (rho_amb / (Msun / pc3)) ** -0.5`, the cluster's `rho_amb` the mean density of its place.
    """
    gamma = tuning.OPEN_CLUSTER_DISRUPTION_GAMMA
    density = np.maximum(np.asarray(rho_amb_msun_pc3, dtype=float), 1.0e-6)
    return c_env_myr * _YEAR_PER_MYR * 10.0 ** (-4.0 * gamma) * density ** -0.5


def present_masses(initial_sol, age_yr, t0_yr, z):
    """
    The mass (solar masses) each cluster holds at `age_yr`, 0 for one disrupted (Lamers et al. 2005, equation 6):
    `mu ** gamma = mu_ev ** gamma - gamma t Msun ** gamma / (t0 Mi ** gamma)`, `mu_ev` what stellar evolution alone
    leaves of the initial mass.
    """
    gamma = tuning.OPEN_CLUSTER_DISRUPTION_GAMMA
    initial = np.asarray(initial_sol, dtype=float)
    kept_by_evolution = 1.0 - mass_lost_by_evolution(age_yr, z)
    inside = kept_by_evolution ** gamma - gamma * np.asarray(age_yr, dtype=float) / (t0_yr * initial ** gamma)
    return np.where(inside > 0.0, initial * np.maximum(inside, 0.0) ** (1.0 / gamma), 0.0)


def half_mass_radii_pc(mass_sol, age_gy, density_z, expansion_unit):
    """
    The half-mass radius (pc) of clusters of `mass_sol` at `age_gy`. At birth `(3 M / (8 pi rho)) ** (1/3)` with the
    half-mass density `rho` lognormal about `tuning.OPEN_CLUSTER_HALF_MASS_DENSITY_LOG10` (`density_z` a standard normal
    draw); it grows by a factor from `tuning.OPEN_CLUSTER_EXPANSION_FACTOR_RANGE` (`expansion_unit`, uniform in 0..1)
    over `tuning.OPEN_CLUSTER_EXPANSION_AGE_GY`, linearly in age, and stays there (Portegies Zwart et al. 2010).
    """
    mean, sigma = tuning.OPEN_CLUSTER_HALF_MASS_DENSITY_LOG10["mean"], tuning.OPEN_CLUSTER_HALF_MASS_DENSITY_LOG10["sigma"]
    density = 10.0 ** (mean + sigma * np.asarray(density_z, dtype=float))
    born = (3.0 * np.asarray(mass_sol, dtype=float) / (8.0 * math.pi * density)) ** (1.0 / 3.0)
    low, high = tuning.OPEN_CLUSTER_EXPANSION_FACTOR_RANGE
    factor = low + (high - low) * np.asarray(expansion_unit, dtype=float)
    grown = np.minimum(np.asarray(age_gy, dtype=float) / tuning.OPEN_CLUSTER_EXPANSION_AGE_GY, 1.0)
    return born * (1.0 + (factor - 1.0) * grown)


def jacobi_radii_pc(mass_sol, galactocentric_pc):
    """`physics.king.jacobi_radius_pc` on arrays: the largest radius a cluster can keep, `R (M / (2 M_gal(<R))) ** (1/3)`
    in an isothermal halo of circular speed `tuning.CLUSTER_HALO_CIRCULAR_KMS`."""
    speed = tuning.CLUSTER_HALO_CIRCULAR_KMS * 1.0e3
    radius_pc = np.maximum(np.asarray(galactocentric_pc, dtype=float), 1.0)
    galaxy_mass = speed * speed * radius_pc * constants.PARSEC_M / constants.G / constants.SOLAR_MASS_TO_KG
    return radius_pc * (np.asarray(mass_sol, dtype=float) / (2.0 * galaxy_mass)) ** (1.0 / 3.0)

