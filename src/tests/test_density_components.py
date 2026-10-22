# tests/test_density_components.py

"""
GEN.379: the density model as a list of components, each with an amplitude that may be 0 and a mass. The masses are in
closed form, agree with the numerical integral (`component_masses`), sum to the galaxy's stellar mass, and a component with
amplitude 0 has mass exactly 0 (the potential, GEN.338, reads the list).
"""

import math

import pytest

from planetgen.galaxy import density
from planetgen.galaxy import structure as galaxyStructure
from planetgen.physics import mathcheck

SHAPE = mathcheck._reference_galaxy_shape()
BASE = dict(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=700.0, bulge_amplitude=3.11,
            arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
NAMES = ("thin_disk", "thick_disk", "bulge", "halo", "nuclear_cluster", "nuclear_disc")


def _by_name(shape):
    return {component.name: component for component in density.density_components(shape)}


def test_the_components_come_in_a_fixed_order_with_a_population_and_shape_each():
    components = density.density_components(SHAPE)
    assert tuple(c.name for c in components) == NAMES
    assert {c.population for c in components} == {"disk", "old", "bulge", "halo", "nuclear_cluster", "nuclear_disc"}
    assert all(c.parameters and c.mass >= 0.0 for c in components)


def test_the_masses_sum_to_the_total_and_the_solar_masses_to_the_stellar_mass():
    for shape in (SHAPE, density.build_galaxy_shape(**BASE, nuclear_amplitude=1.0)):
        components = density.density_components(shape)
        assert sum(c.mass for c in components) == pytest.approx(density.total_mass(shape), rel=1e-12)
        assert sum(c.mass_solar for c in components) == pytest.approx(shape.structure.stellar_mass_solar, rel=1e-12)


def test_the_closed_form_masses_agree_with_the_numerical_integral():
    by_name = _by_name(SHAPE)
    numeric = density.component_masses(SHAPE, steps=80)
    assert by_name["thin_disk"].mass == pytest.approx(numeric["thin_disk"], rel=0.01)
    assert by_name["thick_disk"].mass == pytest.approx(numeric["thick_disk"], rel=0.01)
    assert by_name["bulge"].mass == pytest.approx(numeric["bulge"], rel=0.03)
    assert by_name["halo"].mass == pytest.approx(density.halo_mass(SHAPE), rel=1e-12)


def test_the_halo_is_the_stellar_halo_share_of_the_galaxys_stars():
    by_name = _by_name(SHAPE)
    ratio = SHAPE.structure.halo_stellar_ratio
    # the nuclear cluster and disc sit on top of the stars the visible components hold, so the halo is a little under its share
    assert by_name["halo"].mass_solar / SHAPE.structure.stellar_mass_solar == pytest.approx(ratio, rel=0.05)


def test_a_component_with_amplitude_zero_has_mass_zero_and_the_rest_are_unchanged():
    with_nucleus = density.build_galaxy_shape(**BASE, nuclear_amplitude=1.0)
    without = density.build_galaxy_shape(**BASE, nuclear_amplitude=0.0)
    for name in ("nuclear_cluster", "nuclear_disc"):
        assert _by_name(without)[name].amplitude == 0.0 and _by_name(without)[name].mass == 0.0
        assert _by_name(with_nucleus)[name].mass > 0.0
    no_thick = density.build_galaxy_shape(**BASE, structure=galaxyStructure.REFERENCE_STRUCTURE._replace(thick_mass_ratio=0.0))
    assert _by_name(no_thick)["thick_disk"].amplitude == 0.0 and _by_name(no_thick)["thick_disk"].mass == 0.0
    assert _by_name(no_thick)["thin_disk"].mass > 0.0
    assert density._thick_disk(1000.0, 10.0, no_thick, density.model_terms(no_thick)) == 0.0
    assert sum(c.mass_solar for c in density.density_components(no_thick)) == pytest.approx(no_thick.structure.stellar_mass_solar)


def test_a_galaxy_with_no_bulge_has_a_zero_mass_bulge():
    flat = density.build_galaxy_shape(**{**BASE, "bulge_amplitude": 0.0})
    assert _by_name(flat)["bulge"].mass == 0.0
    assert sum(c.mass_solar for c in density.density_components(flat)) == pytest.approx(flat.structure.stellar_mass_solar)
    assert density.relative_density((0.0, 0.0, 0.0), flat) > 0.0
