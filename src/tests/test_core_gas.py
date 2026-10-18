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
from planetgen.physics import open_cluster_evolution as evolution

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
    assert ring.inflow_msun_yr == pytest.approx(0.8, rel=0.3)          # 0.8 +- 0.6 (Hatchfield et al. 2021)
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
    nearby = [c for c in found if math.hypot(c.center_pc[0], c.center_pc[1]) < 400.0]
    assert 50 < len(nearby) < 3000
    ring = core_gas.core_gas(shape)
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
