# tests/test_cluster_unbinding.py

"""GEN.1112 RP7: the planets a birth cluster strips from their stars (`planetgen.physics.cluster_unbinding`, Spurzem et al. 2009 eqs. 36 and 37)."""

import math

import pytest

from planetgen.physics import cluster_unbinding as cu


def test_the_unbinding_time_is_the_papers_formula():
    # tau = 3 V / (40 pi G rho a), G = 4.30091e-3 pc (km/s)^2 / Msun; one pc/(km/s) is 0.97779 Myr
    a_pc = 100.0 * 4.848136811e-6
    expected = 3.0 * 2.0 / (40.0 * math.pi * 4.30091e-3 * 100.0 * a_pc) * 0.97779
    assert cu.unbinding_time_myr(100.0, 100.0, 2.0) == pytest.approx(expected, rel=1e-4)
    assert cu.unbinding_time_myr(1000.0, 100.0, 2.0) == pytest.approx(22.4, abs=0.1)


def test_the_time_is_one_over_the_density_and_the_orbit_and_grows_with_the_speed():
    base = cu.unbinding_time_myr(300.0, 500.0, 1.5)
    assert cu.unbinding_time_myr(300.0, 5000.0, 1.5) == pytest.approx(base / 10.0)
    assert cu.unbinding_time_myr(3000.0, 500.0, 1.5) == pytest.approx(base / 10.0)
    assert cu.unbinding_time_myr(300.0, 500.0, 3.0) == pytest.approx(2.0 * base)


def test_the_cross_section_and_the_time_agree_through_n_sigma_v():
    star_mass, density, speed, axis = 1.0, 250.0, 2.0, 400.0
    sigma_pc2 = cu.cross_section_au2(axis, star_mass, speed) * (4.848136811e-6) ** 2
    rate_per_myr = (density / star_mass) * sigma_pc2 * speed / 0.97779               # n sigma V in 1/Myr
    assert 1.0 / rate_per_myr == pytest.approx(cu.unbinding_time_myr(axis, density, speed), rel=1e-3)


def test_the_cross_section_grows_with_the_orbit_and_falls_with_the_speed():
    assert cu.cross_section_au2(200.0, 1.0, 2.0) == pytest.approx(2.0 * cu.cross_section_au2(100.0, 1.0, 2.0))
    assert cu.cross_section_au2(100.0, 1.0, 4.0) == pytest.approx(cu.cross_section_au2(100.0, 1.0, 2.0) / 4.0)


def test_the_survival_over_a_history_is_a_product_of_its_pieces():
    one = cu.survival_fraction(1000.0, [(10.0, 100.0, 2.0)])
    assert one == pytest.approx(math.exp(-10.0 / 22.39), rel=1e-3)
    two = cu.survival_fraction(1000.0, [(10.0, 100.0, 2.0), (10.0, 100.0, 2.0)])
    assert two == pytest.approx(one ** 2)
    assert cu.survival_fraction(1000.0, [(10.0, 100.0, 2.0), (500.0, 0.0, 2.0)]) == pytest.approx(one)
    assert cu.unbound_fraction(1000.0, [(10.0, 100.0, 2.0)]) == pytest.approx(1.0 - one)
    assert cu.survival_fraction(1000.0, []) == 1.0


def test_a_denser_longer_cluster_strips_more_and_wider_orbits_first():
    short = cu.unbound_fraction(500.0, [(5.0, 100.0, 2.0)])
    long = cu.unbound_fraction(500.0, [(50.0, 100.0, 2.0)])
    dense = cu.unbound_fraction(500.0, [(5.0, 1000.0, 2.0)])
    assert short < long and short < dense
    assert cu.unbound_fraction(5000.0, [(5.0, 100.0, 2.0)]) > short > cu.unbound_fraction(50.0, [(5.0, 100.0, 2.0)])


def test_the_relative_speed_and_the_half_mass_density():
    assert cu.relative_speed_km_s(1.0) == pytest.approx(math.sqrt(6.0))
    assert cu.half_mass_density(1000.0, 1.0) == pytest.approx(3000.0 / (8.0 * math.pi), rel=1e-9)
    assert cu.half_mass_density(1000.0, 2.0) == pytest.approx(cu.half_mass_density(1000.0, 1.0) / 8.0)


def test_soft_planets_are_the_ones_the_formula_was_derived_for():
    assert cu.soft(1000.0, 1.0, 2.0)                                   # circular speed 0.94 km/s under 2 km/s
    assert not cu.soft(5.0, 1.0, 2.0)                                  # Jupiter's orbit is hard against a young cluster


def test_bad_inputs_are_refused():
    with pytest.raises(cu.ClusterUnbindingError):
        cu.unbinding_time_myr(0.0, 100.0, 2.0)
    with pytest.raises(cu.ClusterUnbindingError):
        cu.survival_fraction(100.0, [(-1.0, 100.0, 2.0)])
    with pytest.raises(cu.ClusterUnbindingError):
        cu.relative_speed_km_s(0.0)
    with pytest.raises(cu.ClusterUnbindingError):
        cu.half_mass_density(1000.0, 0.0)
