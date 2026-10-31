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
from planetgen.physics import constants, unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")

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


def mean_initial_mass_sol(high_mass_sol=None):
    """The mean initial mass of a cluster, solar masses, for the power-law mass function between
    `tuning.OPEN_CLUSTER_MASS_RANGE_SOL` (its high end held under `high_mass_sol`, a maximum cluster mass,
    `physics.cluster_maximum_mass`; an array gives an array): how many clusters a mass of bound star formation makes."""
    low, high = _mass_range(high_mass_sol)
    slope = tuning.OPEN_CLUSTER_MASS_FUNCTION_SLOPE
    if abs(slope - 2.0) < 1.0e-12:
        return np.log(high / low) / (1.0 / low - 1.0 / high)
    one = 1.0 - slope
    two = 2.0 - slope
    return ((high ** two - low ** two) / two) / ((high ** one - low ** one) / one)


def _mass_range(high_mass_sol, low_mass_sol=None):
    """`tuning.OPEN_CLUSTER_MASS_RANGE_SOL` with its high end held under `high_mass_sol` (a number or an array) and over twice
    its low end (`low_mass_sol` when given)."""
    low, high = tuning.OPEN_CLUSTER_MASS_RANGE_SOL
    if low_mass_sol is not None:
        low = low_mass_sol
    if high_mass_sol is not None:
        high = np.minimum(high, np.maximum(high_mass_sol, 2.0 * low))
    return low, high


def initial_masses(uniform, low_mass_sol=None, high_mass_sol=None):
    """Initial masses (solar masses) by inverse transform of the power-law mass function, from uniform draws in 0..1,
    between `low_mass_sol` (the range's low end when not given) and its high end (held under `high_mass_sol`, a maximum
    cluster mass: a number, or an array with a value for each draw)."""
    low, high = _mass_range(high_mass_sol, low_mass_sol)
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


def birth_half_mass_radii_pc(mass_sol, density_z):
    """The half-mass radius (pc) of clusters of `mass_sol` at birth: `(3 M / (8 pi rho)) ** (1/3)` with the half-mass density
    `rho` lognormal about `tuning.OPEN_CLUSTER_HALF_MASS_DENSITY_LOG10` (`density_z` a standard normal draw)."""
    mean, sigma = tuning.OPEN_CLUSTER_HALF_MASS_DENSITY_LOG10["mean"], tuning.OPEN_CLUSTER_HALF_MASS_DENSITY_LOG10["sigma"]
    density = 10.0 ** (mean + sigma * np.asarray(density_z, dtype=float))
    return (3.0 * np.asarray(mass_sol, dtype=float) / (8.0 * math.pi * density)) ** (1.0 / 3.0)


def expansion_factors(expansion_unit):
    """The factor a cluster's half-mass radius grows by, from `tuning.OPEN_CLUSTER_EXPANSION_FACTOR_RANGE` (`expansion_unit`
    uniform in 0..1)."""
    low, high = tuning.OPEN_CLUSTER_EXPANSION_FACTOR_RANGE
    return low + (high - low) * np.asarray(expansion_unit, dtype=float)


def half_mass_radii_pc(mass_sol, age_gy, density_z, expansion_unit):
    """
    The half-mass radius (pc) of clusters of `mass_sol` at `age_gy`. At birth `birth_half_mass_radii_pc`; it grows by
    `expansion_factors` over `tuning.OPEN_CLUSTER_EXPANSION_AGE_GY`, linearly in age, and stays there (Portegies Zwart et al.
    2010).
    """
    born = birth_half_mass_radii_pc(mass_sol, density_z)
    factor = expansion_factors(expansion_unit)
    grown = np.minimum(np.asarray(age_gy, dtype=float) / tuning.OPEN_CLUSTER_EXPANSION_AGE_GY, 1.0)
    return born * (1.0 + (factor - 1.0) * grown)


def cloud_shock_mass_lost_sol(born_radius_pc, factor, age_gy, cloud_surface_density_msun_pc2, cloud_density_msun_pc3):
    """
    The mass (solar masses) passing molecular clouds have taken from clusters by `age_gy` (Gieles et al. 2006, Eq 24, as
    Kruijssen et al. 2014 put it in the central ring): the disruption time `t_dis = 2.0 S (M / 1e4) (3.75 pc / r_h) ** 3` Gyr
    with `S = 5.1 / (Sigma_n rho_n)` makes the loss rate `M / t_dis = 1e4 (r_h / 3.75 pc) ** 3 Sigma_n rho_n / (2.0 x 5.1)` per
    Gyr, the same for every mass, so it is integrated along the half-mass radius `r_h` each cluster has had: `born_radius_pc`
    growing by `factor` over `tuning.OPEN_CLUSTER_EXPANSION_AGE_GY` and staying (`half_mass_radii_pc`). Loss on top of the tidal
    and evolutionary loss of `present_masses`, taken from the mass left.
    """
    born = np.asarray(born_radius_pc, dtype=float)
    growth = np.asarray(factor, dtype=float) - 1.0
    age = np.asarray(age_gy, dtype=float)
    span = tuning.OPEN_CLUSTER_EXPANSION_AGE_GY
    slope = born * growth / span                   # pc per Gy while the radius grows
    grown_for = np.minimum(age, span)
    with np.errstate(divide="ignore", invalid="ignore"):
        growing = np.where(slope > 0.0, ((born + slope * grown_for) ** 4 - born ** 4) / (4.0 * slope), born ** 3 * grown_for)
    steady = (born * (1.0 + growth)) ** 3 * np.maximum(age - span, 0.0)
    rate = (tuning.CLUSTER_SHOCK_REFERENCE_MASS_SOL * cloud_surface_density_msun_pc2 * cloud_density_msun_pc3
            / (tuning.CLUSTER_SHOCK_GIELES_COEFFICIENT_GYR * tuning.CLUSTER_SHOCK_GIELES_NUMERATOR
               * tuning.CLUSTER_SHOCK_REFERENCE_RADIUS_PC ** 3))
    return rate * (growing + steady)


def jacobi_radii_pc(mass_sol, galactocentric_pc):
    """`physics.king.jacobi_radius_pc` on arrays: the largest radius a cluster can keep, `R (M / (2 M_gal(<R))) ** (1/3)`
    in an isothermal halo of circular speed `tuning.CLUSTER_HALO_CIRCULAR_KMS`."""
    speed = tuning.CLUSTER_HALO_CIRCULAR_KMS * 1.0e3
    radius_pc = np.maximum(np.asarray(galactocentric_pc, dtype=float), 1.0)
    galaxy_mass = _kg_to_msun(speed * speed * radius_pc * constants.PARSEC_M / constants.G)
    return radius_pc * (np.asarray(mass_sol, dtype=float) / (2.0 * galaxy_mass)) ** (1.0 / 3.0)

