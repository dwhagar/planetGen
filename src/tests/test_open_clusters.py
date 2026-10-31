# tests/test_open_clusters.py

"""
GEN.244: a generated galaxy's young open clusters. The draw (`generation/open_clusters.py`) is deterministic and follows
the papers' laws (`physics/open_cluster_evolution.py`); the clusters reach the density model, the profile and the
component pass through the twelve age-bin populations.
"""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import component_field, components, core_gas, potential, profile
from planetgen.galaxy import density as galaxyDensity
from planetgen.generation import component_pass, open_clusters
from planetgen.physics import open_cluster_evolution as evolution

EDGE = 4.0
E1 = 6.3


def _shape(found=()):
    shape = galaxyDensity.build_galaxy_shape(
        disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=700.0, bulge_amplitude=2.0,
        arm_count=2, pitch_angle_rad=math.radians(12.0), arm_amplitude=0.3)
    if found:
        shape = shape._replace(components=component_field.ComponentField(found, EDGE, E1))
    return shape


@pytest.fixture(scope="module")
def drawn():
    shape = _shape()
    return shape, open_clusters.generate_open_clusters(b"gen244-test-seed", shape, EDGE, E1)


def test_the_same_seed_and_shape_give_the_same_clusters(drawn):
    shape, found = drawn
    assert found == open_clusters.generate_open_clusters(b"gen244-test-seed", shape, EDGE, E1)
    assert found != open_clusters.generate_open_clusters(b"another-seed", shape, EDGE, E1)


def test_the_clusters_are_bound_young_and_inside_the_mass_function(drawn):
    _shape_, found = drawn
    assert len(found) > 100
    low, high = tuning.OPEN_CLUSTER_MASS_RANGE_SOL
    for cluster in found:
        assert cluster.kind == "cluster" and cluster.profile == "plummer" and components.is_open_cluster(cluster)
        assert low <= cluster.mass_solar <= high
        assert tuning.OPEN_CLUSTER_AGE_RANGE_MYR[0] / 1000.0 <= cluster.age_gy <= tuning.OPEN_CLUSTER_AGE_RANGE_MYR[1] / 1000.0
        assert cluster.tidal_radius_pc > cluster.half_mass_radius_pc > 0.0
    assert [c.name for c in found[:2]] == ["Open cluster 00001", "Open cluster 00002"]


def test_most_clusters_are_young_because_the_old_ones_are_disrupted(drawn):
    _shape_, found = drawn
    ages = np.array([c.age_gy for c in found])
    assert np.median(ages) < 0.3
    assert (ages < 0.1).sum() > (ages > 0.5).sum()


def test_the_clusters_lie_in_the_disk_on_circular_orbits(drawn):
    shape, found = drawn
    radius = np.array([math.hypot(c.center_pc[0], c.center_pc[1]) for c in found])
    height = np.array([abs(c.center_pc[2]) for c in found])
    assert np.median(height) < shape.disk_scale_height_pc
    assert radius.max() < 10.0 * shape.disk_scale_length_pc
    speed = np.array([math.hypot(*c.velocity_kms[:2]) for c in found])
    circular = potential.circular_speed_kms(potential.for_galaxy(shape, b"gen244-test-seed"), radius)
    assert abs(np.median(speed / circular) - 1.0) < 0.1


def test_the_formation_efficiency_rises_with_the_rate_density():
    sigma = np.array([1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0, 10.0, 100.0, 0.0])
    efficiency = evolution.cluster_formation_efficiency(sigma)
    assert np.all(np.diff(efficiency[:-1]) > 0.0)
    assert efficiency[-1] == 0.0 and efficiency[-2] < 1.0
    assert efficiency[0] < 0.05      # Kruijssen 2012: a few per cent where the gas is thin


def test_the_initial_masses_follow_a_power_law_of_slope_minus_two():
    masses = evolution.initial_masses(np.linspace(0.0005, 0.9995, 200000))
    low, high = tuning.OPEN_CLUSTER_MASS_RANGE_SOL
    assert masses.min() >= low and masses.max() <= high
    grid = np.geomspace(low, high, 400001)
    pdf = grid ** -tuning.OPEN_CLUSTER_MASS_FUNCTION_SLOPE
    mean = np.trapezoid(grid * pdf * grid, np.log(grid)) / np.trapezoid(pdf * grid, np.log(grid))
    assert evolution.mean_initial_mass_sol() == pytest.approx(mean, rel=1e-4)
    # dN/dM ~ M^-2: equal numbers per factor-of-two in the cumulative of 1/M
    assert (masses < 2.0 * low).mean() == pytest.approx(0.5, abs=0.01)


def test_stellar_evolution_takes_nothing_before_12_myr_and_more_at_older_ages():
    z = tuning.OPEN_CLUSTER_SOLAR_Z
    age = np.array([1.0e6, 1.2e7, 1.0e8, 1.0e9, 1.0e10])
    lost = evolution.mass_lost_by_evolution(age, z)
    assert lost[0] == 0.0 and lost[1] == 0.0
    assert np.all(np.diff(lost[1:]) > 0.0) and lost[3] == pytest.approx(0.245, abs=0.05)


def test_a_massive_cluster_outlives_a_light_one_and_a_dense_place_disrupts_faster():
    t0 = evolution.disruption_t0_yr(0.1, 500.0)
    kept = evolution.present_masses(np.array([1.0e2, 1.0e3, 1.0e4]), 3.0e8, t0, 0.02)
    assert kept[0] == 0.0 and kept[2] > kept[1] > 0.0
    assert evolution.disruption_t0_yr(1.0, 500.0) < t0
    assert evolution.present_masses(np.array([1.0e3]), 3.0e8, evolution.disruption_t0_yr(1.0, 500.0), 0.02)[0] < kept[1]


def test_the_half_mass_radius_grows_over_the_first_100_myr_and_stays_inside_the_tidal_limit():
    young = evolution.half_mass_radii_pc(1.0e3, 0.0, 0.0, 0.5)
    grown = evolution.half_mass_radii_pc(1.0e3, 0.1, 0.0, 0.5)
    assert grown == pytest.approx(young * 7.5, rel=1e-9)
    assert evolution.half_mass_radii_pc(1.0e3, 1.0, 0.0, 0.5) == pytest.approx(grown)
    assert evolution.jacobi_radii_pc(1.0e3, 8000.0, 9.0e10) == pytest.approx(14.17, abs=0.3)


def test_the_population_of_an_age_is_its_quarter_decade_bin():
    names = tuning.OPEN_CLUSTER_POPULATIONS
    assert len(names) == 12 and names[0] == "open_cluster_00"
    assert components.open_cluster_population(0.0015) == "open_cluster_00"
    assert components.open_cluster_population(0.9) == "open_cluster_11"
    for k, name in enumerate(names):
        low, high = tuning.STELLAR_POPULATION_AGE_RANGES_GY[name]
        assert components.open_cluster_population(math.sqrt(low * high)) == name
        assert low == pytest.approx(tuning.OPEN_CLUSTER_AGE_EDGES_MYR[k] / 1000.0)


def _one(age_gy=0.05, mass=800.0, center=(3000.0, 0.0, 10.0)):
    return components.Component(
        kind="cluster", name="Test open", source="generated", center_pc=center, velocity_kms=(0.0, 220.0, 0.0),
        epoch_unix=0.0, mass_solar=mass, half_mass_radius_pc=3.0, profile="plummer", age_gy=age_gy,
        metallicity_fe_h=-0.05, seed=7, tidal_radius_pc=15.0)


def test_a_cluster_adds_its_density_to_the_population_of_its_age():
    found = [_one()]
    shape = _shape(found)
    plain = _shape()
    population = components.open_cluster_population(0.05)
    inside = (3000.0, 0.0, 10.0)
    mix = galaxyDensity.population_densities(inside, shape)
    assert mix[population] > 0.0
    assert all(mix[name] == 0.0 for name in tuning.OPEN_CLUSTER_POPULATIONS if name != population)
    assert galaxyDensity.relative_density(inside, shape) > galaxyDensity.relative_density(inside, plain)
    assert sum(mix.values()) == pytest.approx(galaxyDensity.relative_density(inside, shape), rel=1e-9)
    far = (-6000.0, 4000.0, 0.0)
    assert galaxyDensity.population_densities(far, shape)[population] == 0.0


def test_the_array_form_leaves_every_population_of_the_gate_to_the_component_pass():
    shape = _shape([_one()])
    xs, ys, zs = np.array([3000.0, -6000.0]), np.array([0.0, 4000.0]), np.array([10.0, 0.0])
    found = galaxyDensity.population_densities_array(xs, ys, zs, shape)
    assert len(component_field.POPULATIONS_IN_THE_GATE) == 16
    for name in component_field.POPULATIONS_IN_THE_GATE:
        assert not found[name].any()


def test_the_profile_names_the_open_cluster_stars_of_the_place():
    where = profile.profile_at((3000.0, 0.0, 10.0), _shape([_one()]), EDGE, E1)
    name = components.open_cluster_population(0.05)
    assert where.populations[name].share > 0.0
    assert where.populations[name].age_gy == tuple(tuning.STELLAR_POPULATION_AGE_RANGES_GY[name])
    assert all(where.populations[other].share == 0.0 for other in tuning.OPEN_CLUSTER_POPULATIONS if other != name)


def test_the_component_pass_places_no_phenomena_for_a_young_cluster():
    assert component_pass.expected_phenomena(_one(), 8.0) == 0.0


def _brute_bound(field, r_cyl, z):
    cylindrical = np.hypot(field._centers[:, 0], field._centers[:, 1])
    nearest = np.hypot(r_cyl - cylindrical, z - field._centers[:, 2]) - tuning.CLUSTER_BOUND_SLACK_EDGES * field.edge_pc
    total = sum(field._systems_per_pc3(int(i), max(float(nearest[i]), 0.0)) * field.edge_pc ** 3
                for i in np.nonzero(nearest < field._reach)[0])
    return total / field.expected_at_density_1


def _brute_high_band(field, shell, band_at):
    high, radius = -1, shell * field.edge_pc
    slack = tuning.CLUSTER_BOUND_SLACK_EDGES * field.edge_pc
    for center, reach in zip(field._centers, field._density_reach):
        if reach <= 0.0:
            continue
        distance, grown = float(np.linalg.norm(center)), reach + slack
        if distance <= grown:
            return band_at(shell, 0.0) if shell else 0
        if abs(distance - radius) > grown:
            continue
        highest = min(math.pi / 2.0, abs(math.asin(max(-1.0, min(1.0, center[2] / distance)))) + math.asin(min(1.0, grown / distance)))
        high = max(high, abs(band_at(shell, math.pi / 2.0 - highest)))
    return high


def test_the_window_search_of_the_bound_and_the_top_band_equals_testing_every_component(drawn):
    from planetgen.galaxy import spherical
    _shape_, found = drawn
    wide = components.Component(
        kind="satellite", name="Wide", source="generated", center_pc=(20000.0, 3000.0, 9000.0), velocity_kms=(0.0, 0.0, 0.0),
        epoch_unix=0.0, mass_solar=1.0e8, half_mass_radius_pc=1500.0, profile="exponential", age_gy=8.0,
        metallicity_fe_h=-1.5, seed=3, tidal_radius_pc=9000.0)
    field = component_field.ComponentField(list(found) + [wide], EDGE, E1)
    assert len(field._big) == 1 and len(field._small) == len(found)
    rng = np.random.default_rng(5)
    for r_cyl, z in zip(rng.uniform(0.0, 12000.0, 150), rng.uniform(-1200.0, 1200.0, 150)):
        assert field.bound(r_cyl, z) == pytest.approx(_brute_bound(field, r_cyl, z), rel=1e-9, abs=1e-12)
    for shell in (0, 1, 7, 100, 400, 900, 1500, 2200, 2800, 5000):
        assert field.high_band(shell, spherical.band_at, EDGE) == _brute_high_band(field, shell, spherical.band_at)


def test_a_grid_cell_crowded_with_clusters_gives_the_same_densities_as_testing_each(drawn):
    shape, found = drawn
    field = component_field.ComponentField(found, EDGE, E1)
    assert max(len(v) for v in field._index.values()) > 24
    rng = np.random.default_rng(9)
    for x, y in zip(rng.uniform(-9000.0, 9000.0, 60), rng.uniform(-9000.0, 9000.0, 60)):
        position = (float(x), float(y), 5.0)
        total = field.relative_density(position)
        near = field._index.get(field._grid_cell(position), ())
        direct = sum(field._cube_systems(i, np.asarray(position)) for i in near) / field.expected_at_density_1
        assert total == pytest.approx(direct, rel=1e-12, abs=0.0)


def test_a_galaxy_that_would_form_millions_draws_the_cap_and_keeps_the_heavier_ones(monkeypatch):
    """A program limit (`OPEN_CLUSTER_MAX_COUNT`): the low end of the masses rises until that many are expected."""
    low, high = tuning.OPEN_CLUSTER_MASS_RANGE_SOL
    assert evolution.low_mass_keeping_share(1.0) == pytest.approx(low)
    floor = evolution.low_mass_keeping_share(0.1)
    assert floor == pytest.approx(1.0 / (0.1 / low + 0.9 / high), rel=1e-9) and floor > 8.0 * low
    assert evolution.initial_masses(np.array([0.0, 0.5]), floor)[0] == pytest.approx(floor)
    monkeypatch.setattr(tuning, "OPEN_CLUSTER_MAX_COUNT", 2_000)
    monkeypatch.setattr(core_gas, "core_gas", lambda shape: core_gas.NO_RING)      # the ring's clusters are not capped
    found = open_clusters.generate_open_clusters(b"capped-seed", _shape(), EDGE, E1)
    assert 0 < len(found) <= 2_000 + 5 * 2_000 ** 0.5


# --- GEN.1196: the disc's clusters are held to the most massive cluster each ring's gas can form ---------------------------


def test_no_disc_cluster_is_heavier_than_the_most_massive_its_rings_can_form(drawn):
    shape, found = drawn
    table = open_clusters._formation_table(shape)
    maximum, efficiency = table[6], table[5]
    assert maximum.shape == efficiency.shape and np.all(maximum >= 0.0)
    assert np.all(maximum <= tuning.CLUSTER_MAX_MASS_CEILING_SOL)
    disc = [c for c in found if np.hypot(c.center_pc[0], c.center_pc[1]) > 1000.0]
    assert disc and max(c.mass_solar for c in disc) <= float(maximum.max()) * (1.0 + 1.0e-9)
    assert max(c.mass_solar for c in disc) < 1.0e5      # not the mass function's 1e6: the Milky Way's heaviest are near 1e4
