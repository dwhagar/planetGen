# tests/test_gen_nuclear_core.py

"""
GEN.280: the nuclear star cluster and the nuclear stellar disc are density terms of the model, scaled to the bulge,
and the Gaussian core of ADM.49 is gone. The ball of shell 0 holds exactly the cluster's mass inside it.
"""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import density, resolution
from planetgen.galaxy.profile import profile_at
from planetgen.galaxy.skeleton import expected_system_count_at_density_1
from planetgen.galaxy.spherical import cell_volume
from planetgen.generation import object_first
from planetgen.physics import constants as physicalConstants
from planetgen.physics.units import ly_to_pc

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
BASE = dict(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
            bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)


def _shape(nuclear=1.0, **changes):
    return density.build_galaxy_shape(**{**BASE, **changes}, nuclear_amplitude=nuclear)


def test_the_bulges_unit_integral_matches_the_numerical_one():
    shape = _shape(nuclear=0.0)
    numerical = density.integrated_masses(shape, steps=60)["bulge"]
    assert density.model_terms(shape)["bulge_mass"] == pytest.approx(numerical, rel=2e-3)


def test_the_masses_follow_the_bulge_and_the_amplitude():
    shape = _shape()
    terms = density.model_terms(shape)
    structure = shape.structure
    # The cluster is the drawn mass in solar masses, in the model's units (GEN.404: those hold the stellar mass).
    assert terms["nuclear_cluster_mass"] == pytest.approx(structure.nuclear_cluster_mass_solar * terms["mass_per_solar"])
    systems_per_unit = shape.k_norm * physicalConstants.LOCAL_STELLAR_DENSITY_LY3 / ly_to_pc(1.0) ** 3
    # (the cluster and the nuclear disc add 3% to the stars the stellar mass makes, so the solar masses come out so far under)
    assert terms["nuclear_cluster_mass"] * systems_per_unit * density.mean_system_mass_solar(shape) == pytest.approx(
        structure.nuclear_cluster_mass_solar, rel=0.05)
    assert terms["nuclear_disc_mass"] == pytest.approx(structure.nuclear_disc_mass_solar * terms["mass_per_solar"])
    half = density.model_terms(_shape(nuclear=0.5))
    assert half["nuclear_cluster_mass"] == pytest.approx(0.5 * terms["nuclear_cluster_mass"])
    for amplitude in (0.0, -1.0):
        none = density.model_terms(_shape(nuclear=amplitude))
        assert none["nuclear_cluster_mass"] == 0.0 and none["nuclear_disc_mass"] == 0.0
        assert density.nuclear_bound(0.0, 0.0, _shape(nuclear=amplitude)) == 0.0


def test_a_smaller_galaxy_gets_a_proportionately_smaller_centre():
    small = density.model_terms(_shape(bulge_scale_radius_pc=158.0, disk_scale_length_pc=260.0, disk_scale_height_pc=30.0))
    large = density.model_terms(_shape())
    assert small["nuclear_cluster_scale_pc"] == pytest.approx(large["nuclear_cluster_scale_pc"] / 10)
    assert small["nuclear_disc_length_pc"] == pytest.approx(large["nuclear_disc_length_pc"] / 10)
    assert small["nuclear_disc_height_pc"] == pytest.approx(large["nuclear_disc_height_pc"] / 10)


def test_outside_the_softening_the_cluster_is_a_hernquist_profile():
    terms = density.model_terms(_shape())
    mass, scale = terms["nuclear_cluster_mass"], terms["nuclear_cluster_scale_pc"]
    for radius in (3.0, 10.0, 100.0, 1000.0):
        assert density._nuclear_cluster(radius, terms) == pytest.approx(mass * scale / (2 * math.pi * radius * (radius + scale) ** 3))
    # The mass inside r is M r^2 / (r + a)^2, and the half-mass radius is (1 + sqrt 2) a.
    assert (1 + math.sqrt(2)) * scale == pytest.approx(5.0, rel=0.01)


def test_the_cusp_is_held_flat_inside_the_softening_and_continuous_at_it():
    terms = density.model_terms(_shape())
    soft = terms["nuclear_cluster_softening_pc"]
    assert 0.0 < soft < tuning.NUCLEAR_RESOLUTION_PC
    inside = [density._nuclear_cluster(radius, terms) for radius in (0.0, soft / 3, soft)]
    assert inside == pytest.approx([inside[0]] * 3)
    just_out = density._nuclear_cluster(soft * (1 + 1e-9), terms)
    assert just_out == pytest.approx(inside[0], rel=1e-6)


def test_the_core_ball_holds_exactly_the_clusters_mass_inside_it():
    for shape in (_shape(), _shape(bulge_scale_radius_pc=40.0, disk_scale_length_pc=40.0, disk_scale_height_pc=12.0)):
        terms = density.model_terms(shape)
        ball = tuning.NUCLEAR_RESOLUTION_PC
        mass_inside = terms["nuclear_cluster_mass"] * ball ** 2 / (ball + terms["nuclear_cluster_scale_pc"]) ** 2
        at_centre = density._nuclear_cluster(0.0, terms)
        assert at_centre * 4.0 / 3.0 * math.pi * ball ** 3 == pytest.approx(mass_inside, rel=1e-9)


def test_the_cell_at_the_centre_resolves_about_the_clusters_systems():
    shape = _shape()
    expected_per_cell_at_1 = expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)
    profile = profile_at((0.0, 0.0, 0.0), shape, EDGE_PC, expected_per_cell_at_1, address=(0, 0, 0),
                         volume_pc3=cell_volume(0, 0, 0, EDGE_PC))
    assert profile.region == "nuclear"
    systems = resolution.expected_in_cell(profile)
    # The reference cluster (2.5e7 solar masses, at the galaxy's mean system mass a system) holds the share `r^2 / (r + a)^2` of
    # them inside the resolution ball, 24% for the 2 pc ball and the 2.07 pc scale radius; the softened cusp and the
    # bulge keep the cell within 5% of it.
    terms = density.model_terms(shape)
    ball = tuning.NUCLEAR_RESOLUTION_PC
    inside = ball ** 2 / (ball + terms["nuclear_cluster_scale_pc"]) ** 2
    assert systems == pytest.approx(shape.structure.nuclear_cluster_mass_solar / density.mean_system_mass_solar(shape) * inside,
                                    rel=0.05)
    assert profile.mix["nuclear_cluster"] > 0.95 * profile.density


def test_the_nuclear_disc_has_the_observed_run_of_density():
    terms = density.model_terms(_shape())
    per_pc3 = terms["nuclear_disc_density"]
    assert density._nuclear_disc(0.0, 0.0, terms) == pytest.approx(per_pc3)
    assert density._nuclear_disc(terms["nuclear_disc_length_pc"], 0.0, terms) == pytest.approx(per_pc3 / math.e)
    assert density._nuclear_disc(0.0, 10 * terms["nuclear_disc_height_pc"], terms) < 1e-3 * per_pc3


def test_the_array_forms_match_the_scalar_ones():
    shape = _shape()
    terms = density.model_terms(shape)
    radii = np.array([0.0, 0.5, 2.0, 7.0, 90.0, 700.0])
    assert density._nuclear_cluster_array(radii, terms) == pytest.approx([density._nuclear_cluster(r, terms) for r in radii])
    rs, zs = np.array([0.0, 30.0, 150.0, 900.0]), np.array([0.0, 10.0, 80.0, 400.0])
    assert density._nuclear_disc_array(rs, zs, terms) == pytest.approx([density._nuclear_disc(r, z, terms) for r, z in zip(rs, zs)])
    assert density._nuclear_disc_array(rs, 25.0, terms) == pytest.approx([density._nuclear_disc(r, 25.0, terms) for r in rs])


def test_the_populations_sum_to_the_density_at_the_centre_too():
    shape = _shape()
    for point in ((0.0, 0.0, 0.0), (3.0, 1.0, 0.5), (60.0, 20.0, 5.0), (400.0, -200.0, 50.0)):
        parts = density.population_densities(point, shape)
        assert sum(parts.values()) == pytest.approx(density.relative_density(point, shape), rel=1e-9)
        arrays = density.population_densities_array(np.array([point[0]]), np.array([point[1]]), point[2], shape)
        for name, value in parts.items():
            if name != "cluster":
                assert arrays[name][0] == pytest.approx(value, rel=1e-9, abs=1e-12), name


def test_the_outline_bound_and_the_band_majorants_bound_the_nuclear_density():
    from planetgen.galaxy import skeleton
    shape = _shape()
    for radius, z in ((0.0, 0.0), (5.0, 2.0), (40.0, 30.0), (300.0, 100.0)):
        bound = skeleton.bound_relative_density_at(shape, radius, z)
        for step in range(12):
            theta = step * math.pi / 6
            point = (radius * math.cos(theta), radius * math.sin(theta), z)
            assert density.relative_density(point, shape) <= bound * (1 + 1e-9) + density.halo_density(point, shape)
    fractions = {name: 1.0 for name in tuning.STELLAR_POPULATION_AGE_RANGES_GY}
    for shell in (0, 1, 2, 5, 12):
        high = 0 if shell < 2 else 2
        majorants = object_first.band_majorants(shape, shell, high, EDGE_PC, fractions)
        assert (majorants > 0).all()
        for index, band in enumerate(range(-high, high + 1)):
            from planetgen.galaxy.spherical import sector_position_pc, wedge_count
            for wedge in range(0, wedge_count(shell, band), max(1, wedge_count(shell, band) // 5)):
                point = sector_position_pc(shell, band, wedge, EDGE_PC)
                total = sum(density.population_densities(point, shape).values())
                assert total <= majorants[index] * (1 + 1e-9) + 1e-9, (shell, band, wedge)
