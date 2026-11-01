# tests/test_cluster_maximum_mass.py

"""
GEN.1195: the most massive cluster a ring of gas forms (`physics/cluster_maximum_mass.py`, [R3720]) and
the clouds that disrupt a ring's clusters (`open_cluster_evolution.cloud_shock_mass_lost_sol`, [R1555] as
[R2505] put it in the central molecular zone).
"""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.physics import cluster_maximum_mass as cmm
from planetgen.physics import open_cluster_evolution as evolution

CMZ_KAPPA = cmm.epicyclic_frequency_per_myr(180.5, 100.0)      # 2.61 per Myr
SIGMA = 5.0         # km/s, the gas of the central ring [R1911]


def test_the_central_rings_toomre_mass_and_cluster_maximum_are_the_papers():
    assert CMZ_KAPPA == pytest.approx(2.61, abs=0.01)
    assert cmm.toomre_mass_sol(1.0e3, CMZ_KAPPA) == pytest.approx(5.3e5, rel=0.02)
    # their central molecular zone: a shear-limited maximum cluster mass of about 3e4 (Section 5), Arches and Quintuplet 2e4
    mass = cmm.maximum_cluster_mass_sol(1.0e3, CMZ_KAPPA, 0.6, SIGMA)
    assert mass == pytest.approx(3.2e4, rel=0.05)
    assert cmm.collapsed_fraction(1.0e3, CMZ_KAPPA, SIGMA) == 1.0           # shear-limited


def test_the_toomre_parameter_is_kappa_sigma_over_pi_g_sigma_g_and_the_discs_solar_neighbourhood_is_feedback_limited():
    # the paper's solar neighbourhood: Sigma_g 13 solar masses per pc2, Omega 0.029 per Myr (kappa 0.041), sigma 10 km/s
    kappa = math.sqrt(2.0) * 0.029
    assert cmm.toomre_q(13.0, kappa, 10.0) == pytest.approx(2.3, abs=0.1)
    assert cmm.toomre_q(1.0e3, CMZ_KAPPA, 5.0) == pytest.approx(0.94, abs=0.05)
    # beyond 4.3 kpc the disc is feedback-limited with a cloud maximum near 1e6 (their section 4.1), not the Toomre mass
    assert cmm.collapsed_fraction(13.0, kappa, 10.0) < 1.0
    assert 1.0e5 < cmm.maximum_cloud_mass_sol(13.0, kappa, 10.0) < 1.0e7


def test_the_toomre_mass_goes_as_the_gas_cubed_and_the_epicyclic_frequency_to_the_minus_fourth():
    base = cmm.toomre_mass_sol(500.0, 2.0)
    assert cmm.toomre_mass_sol(1000.0, 2.0) == pytest.approx(8.0 * base)
    assert cmm.toomre_mass_sol(500.0, 4.0) == pytest.approx(base / 16.0)


def test_feedback_never_acts_before_the_first_supernova_and_limits_a_thin_slow_ring():
    assert cmm.feedback_time_myr(1.0, 5.0, 10.0) >= tuning.CLUSTER_FEEDBACK_FIRST_SUPERNOVA_MYR
    # a slow rotator (low kappa) has a long free-fall time, so feedback ends the collapse first: feedback-limited
    thin = cmm.collapsed_fraction(10.0, 0.05, 10.0)
    assert 0.0 < thin < 1.0
    assert cmm.maximum_cloud_mass_sol(10.0, 0.05, 10.0) == pytest.approx(cmm.toomre_mass_sol(10.0, 0.05) * thin)
    # the same gas on a fast rotator is shear-limited
    assert cmm.collapsed_fraction(1.0e3, 3.0, SIGMA) == 1.0


def test_the_cluster_maximum_is_a_tenth_of_the_bound_share_of_the_cloud_and_held_under_the_ceiling():
    cloud = cmm.maximum_cloud_mass_sol(600.0, CMZ_KAPPA, SIGMA)
    assert cmm.maximum_cluster_mass_sol(600.0, CMZ_KAPPA, 0.5, SIGMA) == pytest.approx(0.1 * 0.5 * cloud)
    assert cmm.maximum_cluster_mass_sol(5.0e3, CMZ_KAPPA, 0.9, SIGMA) == tuning.CLUSTER_MAX_MASS_CEILING_SOL
    assert cmm.maximum_cluster_mass_sol(0.0, CMZ_KAPPA, 0.6, SIGMA) == 0.0


def test_a_cloud_takes_five_point_seven_myr_to_take_ten_thousand_solar_masses_in_the_central_ring():
    # [R2504, R2505]: a cluster of 1e4 solar masses and half-mass radius 0.5 pc lives 5.7 Myr among the ring's clouds
    lost = evolution.cloud_shock_mass_lost_sol(0.5, 1.0, 5.7e-3, tuning.CORE_RING_CLOUD_SURFACE_DENSITY_MSUN_PC2,
                                               tuning.CORE_RING_CLOUD_DENSITY_MSUN_PC3)
    assert lost == pytest.approx(1.0e4, rel=0.03)


def test_the_cloud_loss_follows_the_growing_radius_and_is_zero_at_birth():
    born, factor = np.array([0.4, 0.8]), np.array([5.0, 10.0])
    args = (tuning.CORE_RING_CLOUD_SURFACE_DENSITY_MSUN_PC2, tuning.CORE_RING_CLOUD_DENSITY_MSUN_PC3)
    assert np.all(evolution.cloud_shock_mass_lost_sol(born, factor, 0.0, *args) == 0.0)
    for age in (0.03, 0.1, 0.4):
        # the integral of r(t)^3 along the radius the cluster has had, by quadrature
        t = np.linspace(0.0, age, 20001)
        radius = born[:, None] * (1.0 + (factor[:, None] - 1.0) * np.minimum(t / tuning.OPEN_CLUSTER_EXPANSION_AGE_GY, 1.0))
        rate = (tuning.CLUSTER_SHOCK_REFERENCE_MASS_SOL * args[0] * args[1]
                / (tuning.CLUSTER_SHOCK_GIELES_COEFFICIENT_GYR * tuning.CLUSTER_SHOCK_GIELES_NUMERATOR
                   * tuning.CLUSTER_SHOCK_REFERENCE_RADIUS_PC ** 3))
        expected = rate * np.trapezoid(radius ** 3, t, axis=1)
        assert evolution.cloud_shock_mass_lost_sol(born, factor, age, *args) == pytest.approx(expected, rel=1e-4)
    older = evolution.cloud_shock_mass_lost_sol(born, factor, 0.4, *args)
    assert np.all(older > evolution.cloud_shock_mass_lost_sol(born, factor, 0.1, *args))


def test_a_capped_mass_function_draws_nothing_above_the_cap_and_makes_more_lighter_clusters():
    uniform = (np.arange(2000) + 0.5) / 2000
    assert evolution.initial_masses(uniform, high_mass_sol=2.0e4).max() <= 2.0e4
    assert evolution.initial_masses(uniform).max() > 1.0e5
    assert evolution.mean_initial_mass_sol(2.0e4) < evolution.mean_initial_mass_sol()
    assert evolution.mean_initial_mass_sol(1.0e9) == pytest.approx(evolution.mean_initial_mass_sol())
