# tests/test_core_gas.py

"""
GEN.282: a barred galaxy's centre has a gas ring, the central molecular zone, fed by the gas of the bar's gap; its young
field stars are part of the nuclear disc (moved into the `young` population), and its bound clusters are components.
"""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import components, core_gas, nebula_field, profile
from planetgen.galaxy import density as galaxyDensity
from planetgen.generation import open_clusters

BASE = dict(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
            bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
EDGE = 4.0
E1 = 6.3


def _shape(nuclear=1.0, barred=True):
    shape = galaxyDensity.build_galaxy_shape(**BASE, nuclear_amplitude=nuclear)
    if not barred:
        shape = shape._replace(structure=shape.structure._replace(bar_present=False))
    return shape


def _without_ring(monkeypatch):
    monkeypatch.setattr(core_gas, "core_gas", lambda shape: core_gas.NO_RING)


def test_the_milky_ways_ring_matches_the_measured_inflow_and_star_formation_rate():
    ring = core_gas.core_gas(_shape())
    assert ring.inflow_msun_yr == pytest.approx(0.8, rel=0.3)          # 0.8 +- 0.6
    assert ring.rate_msun_yr == pytest.approx(0.07, rel=0.3)           # 0.07 (Henshaw et al. 2023, Table 1)
    assert ring.rate_msun_yr == pytest.approx(tuning.CORE_STAR_FORMATION_PER_INFLOW * ring.inflow_msun_yr)
    assert 0.0 < ring.bound_share < 1.0
    assert ring.radius_pc == pytest.approx(100.0)


def test_a_galaxy_without_a_bar_has_no_ring():
    assert core_gas.core_gas(_shape(barred=False)) is core_gas.NO_RING


def test_the_ring_scales_with_the_bulge_and_the_gas_of_the_gap():
    small = galaxyDensity.build_galaxy_shape(**{**BASE, "bulge_scale_radius_pc": 158.0}, nuclear_amplitude=1.0)
    assert core_gas.core_gas(small).radius_pc == pytest.approx(core_gas.core_gas(_shape()).radius_pc / 10.0)


def test_the_rings_young_stars_come_out_of_the_nuclear_disc_and_the_total_does_not_move(monkeypatch):
    shape = _shape()
    with_ring = galaxyDensity.population_densities((100.0, 0.0, 0.0), shape)
    total_with = galaxyDensity.relative_density((100.0, 0.0, 0.0), shape)
    _without_ring(monkeypatch)
    plain = galaxyDensity.population_densities((100.0, 0.0, 0.0), shape)
    assert with_ring["young"] > plain["young"] + 1.0
    assert with_ring["nuclear_disc"] == pytest.approx(plain["nuclear_disc"] - (with_ring["young"] - plain["young"]))
    assert sum(with_ring.values()) == pytest.approx(sum(plain.values()), rel=1e-12)
    assert sum(with_ring.values()) == pytest.approx(total_with, rel=1e-9)
    assert galaxyDensity.relative_density((100.0, 0.0, 0.0), shape) == pytest.approx(total_with)


def test_the_young_share_of_the_nuclear_disc_is_about_a_percent_at_the_ring_and_negligible_far_from_it():
    shape = _shape()
    at_ring = galaxyDensity.population_densities((100.0, 0.0, 0.0), shape)
    share = core_gas.young_density(100.0, 0.0, core_gas.core_gas(shape)) * shape.k_norm / (at_ring["nuclear_disc"] + at_ring["young"])
    assert 0.003 < share < 0.05
    ring = core_gas.core_gas(shape)
    assert core_gas.young_density(1000.0, 0.0, ring) < 1.0e-6 * ring.young_density
    assert core_gas.young_density(100.0, 600.0, ring) < 1.0e-6 * ring.young_density


def test_a_galaxy_with_no_nucleus_is_unchanged(monkeypatch):
    shape = _shape(nuclear=0.0)
    with_ring = galaxyDensity.population_densities((100.0, 0.0, 0.0), shape)
    _without_ring(monkeypatch)
    assert galaxyDensity.population_densities((100.0, 0.0, 0.0), shape) == pytest.approx(with_ring)


def test_the_array_form_agrees_with_the_point_form():
    shape = _shape()
    xs = np.array([0.0, 60.0, 100.0, 130.0, 400.0])
    zs = np.array([0.0, 10.0, -20.0, 5.0, 0.0])
    found = galaxyDensity.population_densities_array(xs, np.zeros(5), zs, shape)
    for i in range(5):
        one = galaxyDensity.population_densities((float(xs[i]), 0.0, float(zs[i])), shape)
        for name in ("young", "nuclear_disc"):
            assert found[name][i] == pytest.approx(one[name], rel=1e-9)


def test_the_ring_raises_the_gas_the_profile_and_the_nebulae_read(monkeypatch):
    shape = _shape()
    at = (100.0, 0.0, 0.0)
    with_ring = nebula_field.gas_factor(at, shape)
    where = profile.profile_at(at, shape, EDGE, E1)
    _without_ring(monkeypatch)
    assert with_ring >= nebula_field.gas_factor(at, shape)
    assert where.populations["young"].share > profile.profile_at(at, shape, EDGE, E1).populations["young"].share


def test_the_bound_share_of_the_rings_stars_forms_clusters_where_the_ring_is():
    shape = _shape()
    found = open_clusters.generate_open_clusters(b"core-seed", shape, EDGE, E1)
    ring = core_gas.core_gas(shape)
    nearby = [c for c in found if math.hypot(c.center_pc[0], c.center_pc[1]) < ring.radius_pc + 4.0 * ring.sigma_pc]
    assert 5 < len(nearby) < 3000         # the clouds that pass the ring disrupt most of what it forms (GEN.1195)
    for cluster in nearby:
        assert components.is_open_cluster(cluster) and cluster.age_gy <= tuning.CORE_STAR_FORMATION_SPAN_GY
        assert cluster.half_mass_radius_pc < 3.0          # a tide of the centre holds a cluster to a pc or so
    radius = np.array([math.hypot(c.center_pc[0], c.center_pc[1]) for c in nearby])
    assert abs(np.median(radius) - ring.radius_pc) < 2.0 * ring.sigma_pc
    assert open_clusters.generate_open_clusters(b"core-seed", shape, EDGE, E1) == found


def test_without_the_ring_the_galaxy_draws_no_cluster_near_its_centre(monkeypatch):
    _without_ring(monkeypatch)
    found = open_clusters.generate_open_clusters(b"core-seed", _shape(), EDGE, E1)
    assert not [c for c in found if math.hypot(c.center_pc[0], c.center_pc[1]) < 300.0 and c.mass_solar > 5.0e3]


def test_the_circular_speed_of_the_centre_is_that_of_the_nucleus_inside_the_radius():
    speed = open_clusters._circular_speed_kms(np.array([10.0, 100.0, 200.0]), _shape())
    assert np.all(speed > 30.0) and np.all(speed < 250.0)
    assert speed[1] == pytest.approx(150.0, rel=0.4)


def test_the_young_stars_bound_covers_the_rings_young_stars():
    """The scatter's majorant for the young population has to hold where the ring's stars were moved into it."""
    from planetgen.generation import object_first
    shape = _shape()
    ring = core_gas.core_gas(shape)
    fractions = {name: 0.0 for name in profile.POPULATIONS}
    fractions["young"] = 1.0
    for radius in (0.5 * ring.radius_pc, ring.radius_pc, 1.5 * ring.radius_pc):
        point = (radius, 0.0, 0.0)
        young = galaxyDensity.population_densities(point, shape)["young"]
        bound = object_first.extent_majorants(shape, np.array([0.9 * radius]), np.array([0.0]), np.array([5.0]),
                                              fractions, np.array([1.1 * radius]))
        assert young <= bound[0]


def test_the_rings_gas_is_the_inflow_over_its_residence_and_its_cluster_maximum_is_a_few_ten_thousand_solar_masses():
    """GEN.1195: the central molecular zone's gas is about 3e7 to 5e7 solar masses at about 1e3 per pc2, and its most massive
    cluster is a few 1e4 (the Arches and the Quintuplet 2e4), certainly under 1e5 (Reina-Campos and Kruijssen 2017)."""
    shape = _shape()
    ring = core_gas.core_gas(shape)
    area = 2.0 * math.pi * ring.radius_pc * 2.0 * math.sqrt(2.0 * math.log(2.0)) * ring.sigma_pc
    assert ring.gas_surface_density_msun_pc2 * area == pytest.approx(ring.inflow_msun_yr * tuning.CORE_RING_GAS_RESIDENCE_MYR * 1.0e6)
    assert 2.0e7 < ring.gas_surface_density_msun_pc2 * area < 6.0e7 and 300.0 < ring.gas_surface_density_msun_pc2 < 2000.0
    assert 5.0e3 < open_clusters._ring_maximum_mass_sol(shape, ring) < tuning.CLUSTER_MAX_MASS_CEILING_SOL
    assert open_clusters._ring_maximum_mass_sol(shape, core_gas.NO_RING) == tuning.OPEN_CLUSTER_MASS_RANGE_SOL[1]


def test_the_rings_young_clusters_are_a_few_not_dozens_and_none_of_the_old_ones_is_heavy():
    """GEN.1195: the ring held 5.8 young clusters above 5e3 solar masses and 12 older ones above 5e4 (the Milky Way has the
    Arches and the Quintuplet and no old one); with the maximum from its gas and the clouds' disruption it holds a few young
    ones and no heavy old one. The observed count is the plausibility target, not a preset."""
    shape = _shape()
    young, old = [], []
    for seed in range(8):
        found = open_clusters.generate_open_clusters(f"ring-count-{seed}".encode(), shape, EDGE, E1)
        ring = core_gas.core_gas(shape)
        near = [c for c in found if math.hypot(c.center_pc[0], c.center_pc[1]) < ring.radius_pc + 4.0 * ring.sigma_pc]
        young.append(sum(c.age_gy < 0.01 and c.mass_solar > 5.0e3 for c in near))
        old.append(sum(c.age_gy > 0.1 and c.mass_solar > 5.0e4 for c in near))
    assert np.mean(young) <= 3.0 and max(young) <= 6
    assert max(old) == 0


def test_the_milky_way_shapes_disc_rings_change_from_shear_to_feedback_limited_at_about_four_kiloparsecs():
    """GEN.1196: Reina-Campos and Kruijssen 2017 (section 4.1) find the Milky Way's disc shear-limited inside 4.3 (3.5 to 6.3)
    kpc and feedback-limited beyond, the cloud maximum there about 1e6; the rings of the Milky Way's shape, with their own
    gas and the flat curve, change over in the same band, and the cluster maximum stays under the observed 1e4 to 4e4."""
    from planetgen.galaxy import disk_history
    from planetgen.physics import cluster_maximum_mass as cmm

    shape = _shape()
    _ages, _step, radii, _width, _formed, efficiency, maximum = open_clusters._formation_table(shape)
    gas = disk_history.history_of(shape).gas_now
    kappa = np.array([cmm.epicyclic_frequency_per_myr(tuning.CLUSTER_HALO_CIRCULAR_KMS, r) for r in radii])
    fraction = np.array([cmm.collapsed_fraction(g, k, tuning.DISC_GAS_DISPERSION_KMS) for g, k in zip(gas, kappa)])
    assert np.all(fraction[(radii > 1000.0) & (radii < 3500.0)] == 1.0)
    far = (radii > 6300.0) & (radii < 15000.0)
    assert np.all(fraction[far] < 1.0) and np.all(fraction[far] > 0.0)
    cloud = np.array([cmm.maximum_cloud_mass_sol(g, k, tuning.DISC_GAS_DISPERSION_KMS) for g, k in zip(gas, kappa)])
    assert 1.0e5 < cloud[(radii > 4000.0) & (radii < 10000.0)].max() < 1.0e7
    assert maximum[(radii > 3000.0) & (radii < 12000.0)].max() < 1.0e5
    assert maximum[(radii > 3000.0) & (radii < 12000.0)].max() > 5.0e2
