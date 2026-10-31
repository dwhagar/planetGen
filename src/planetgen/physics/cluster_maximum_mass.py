# planetgen/physics/cluster_maximum_mass.py

"""
The most massive cluster a ring of gas can form (GEN.1195; Reina-Campos and Kruijssen 2017, MNRAS 469, 1282, arXiv:1704.00732).
Pure functions of the gas surface density `Sigma_g`, the epicyclic frequency `kappa` and the one-dimensional velocity
dispersion `sigma` of the gas (the Toomre parameter is `Q = kappa sigma / (pi G Sigma_g)`, Eq 5):

* the Toomre mass `M_T = 4 pi^5 G^2 Sigma_g^3 / kappa^4` (Eq 2), the largest cloud shear lets collapse;
* the two-dimensional free-fall time `t_ff,2D = sqrt(2 pi) / kappa` (Eq 3) and the feedback time `t_fb` (Eqs 4 to 7), after
  which feedback stops a cloud's collapse: the clouds collapse only the fraction `min(1, t_fb / t_ff,2D) ** 4` of the
  Toomre mass (Eqs 8 and 9), shear-limited where feedback is slower and feedback-limited where it is faster;
* the cluster, `eps Gamma M_GMC,max` (Eq 10), a share `eps` of the cloud into stars and `Gamma` of those bound.

For the Milky Way's central ring (`Sigma_g` 1e3 solar masses per pc2, a circular speed of 180 km/s at 100 pc) it gives
3e4, against the 2e4 of the Arches and the Quintuplet. The mass goes as `Sigma_g ** 3` (a factor of two in the gas is a factor
of eight in the mass), and the paper's uncertainty is about 0.4 dex. The ceiling of `tuning.CLUSTER_MAX_MASS_CEILING_SOL`
is the paper's own footnote 9.
"""

import math

from planetgen import tuning
from planetgen.physics import constants
from planetgen.physics import unit_registry as units

_PC_M = constants.PARSEC_M
_MYR_S = constants.SECONDS_PER_YEAR * 1.0e6
_G_PC3_MSUN_MYR2 = constants.G * units.converter("Msun", "kg")(1.0) * _MYR_S ** 2 / _PC_M ** 3
"""float: The gravitational constant in parsecs cubed per solar mass per Myr squared (4.5e-3)."""
_CM2_S3_TO_PC2_MYR3 = 1.0e-4 * _MYR_S ** 3 / _PC_M ** 2


def epicyclic_frequency_per_myr(circular_speed_kms, radius_pc):
    """`kappa` for a flat rotation curve (Eq 11): `sqrt(2) V / R`, per Myr."""
    return math.sqrt(2.0) * circular_speed_kms * 1.0e3 * _MYR_S / (_PC_M * radius_pc)


def toomre_mass_sol(gas_surface_density_msun_pc2, kappa_per_myr):
    """The Toomre mass (Eq 2), solar masses: `4 pi^5 G^2 Sigma_g^3 / kappa^4`."""
    return (4.0 * math.pi ** 5 * _G_PC3_MSUN_MYR2 ** 2 * gas_surface_density_msun_pc2 ** 3 / kappa_per_myr ** 4)


def free_fall_time_2d_myr(kappa_per_myr):
    """The free-fall time of a collapsing sheet (Eq 3), Myr: `sqrt(2 pi) / kappa`."""
    return math.sqrt(2.0 * math.pi) / kappa_per_myr


def toomre_q(gas_surface_density_msun_pc2, kappa_per_myr, dispersion_kms):
    """The Toomre parameter of the gas (Eq 5): `kappa sigma / (pi G Sigma_g)`."""
    sigma_pc_myr = dispersion_kms * 1.0e3 * _MYR_S / _PC_M
    return kappa_per_myr * sigma_pc_myr / (math.pi * _G_PC3_MSUN_MYR2 * gas_surface_density_msun_pc2)


def feedback_time_myr(gas_surface_density_msun_pc2, kappa_per_myr, dispersion_kms):
    """
    The time after which feedback stops a cloud's collapse (Eqs 4 to 7), Myr: `t_sn / 2 (1 + sqrt(1 + 4 pi^2 G^2 t_ff,g Q^2
    Sigma_g^2 / (phi_fb eps_ff t_sn^2 kappa^2)))`, with the mid-plane density `phi_P kappa^2 / (2 pi Q^2 G)` and its free-fall
    time `t_ff,g = sqrt(3 pi / (32 G rho_g))`, `Q` from `toomre_q`. Never under the time of the first supernova.
    """
    q = toomre_q(gas_surface_density_msun_pc2, kappa_per_myr, dispersion_kms)
    t_sn = tuning.CLUSTER_FEEDBACK_FIRST_SUPERNOVA_MYR
    density = tuning.CLUSTER_MIDPLANE_GRAVITY_FACTOR * kappa_per_myr ** 2 / (2.0 * math.pi * q * q * _G_PC3_MSUN_MYR2)
    free_fall_g = math.sqrt(3.0 * math.pi / (32.0 * _G_PC3_MSUN_MYR2 * density))
    phi_fb = tuning.CLUSTER_FEEDBACK_ENERGY_RATE_CM2_S3 * _CM2_S3_TO_PC2_MYR3
    inside = (4.0 * math.pi ** 2 * _G_PC3_MSUN_MYR2 ** 2 * free_fall_g * q * q * gas_surface_density_msun_pc2 ** 2
              / (phi_fb * tuning.CLUSTER_STAR_FORMATION_PER_FREE_FALL * t_sn ** 2 * kappa_per_myr ** 2))
    return 0.5 * t_sn * (1.0 + math.sqrt(1.0 + inside))


def collapsed_fraction(gas_surface_density_msun_pc2, kappa_per_myr, dispersion_kms):
    """The share of the Toomre mass that collapses before feedback (Eq 8): `min(1, t_fb / t_ff,2D) ** 4`."""
    ratio = (feedback_time_myr(gas_surface_density_msun_pc2, kappa_per_myr, dispersion_kms)
             / free_fall_time_2d_myr(kappa_per_myr))
    return min(1.0, ratio) ** 4


def maximum_cloud_mass_sol(gas_surface_density_msun_pc2, kappa_per_myr, dispersion_kms):
    """The most massive cloud that collapses (Eq 9), solar masses: the Toomre mass times the collapsed fraction."""
    if gas_surface_density_msun_pc2 <= 0.0 or kappa_per_myr <= 0.0:
        return 0.0
    return (toomre_mass_sol(gas_surface_density_msun_pc2, kappa_per_myr)
            * collapsed_fraction(gas_surface_density_msun_pc2, kappa_per_myr, dispersion_kms))


def maximum_cluster_mass_sol(gas_surface_density_msun_pc2, kappa_per_myr, bound_share, dispersion_kms):
    """
    The most massive cluster the ring forms (Eq 10), solar masses: `eps Gamma M_GMC,max` for `bound_share` the share of
    the new stars that are bound (`Gamma`), held under `tuning.CLUSTER_MAX_MASS_CEILING_SOL`.
    """
    mass = (tuning.CLUSTER_MAX_MASS_EFFICIENCY * bound_share
            * maximum_cloud_mass_sol(gas_surface_density_msun_pc2, kappa_per_myr, dispersion_kms))
    return min(mass, tuning.CLUSTER_MAX_MASS_CEILING_SOL)
