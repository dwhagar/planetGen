# tests/test_density_components.py

"""
GEN.379: the density model as a list of components, each with an amplitude that may be 0 and a mass. The masses are in
closed form, agree with an independent quadrature of each profile (to 1e-6, the item asks 1e-4), sum to the galaxy's stellar mass, and a
component with amplitude 0 has mass exactly 0 (the potential, GEN.338, reads the list). Today's spiral keeps its density
to the last digit (`fixtures/density_identity.json`, made from the code before the rewrite), and the stored family and
component numbers (DB.66) are read with NULL meaning the structure draw's own value.
"""

import json
import math
import os

import numpy as np
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


def _gauss(function, low, high, nodes):
    """The integral of `function` (an array function) over `[low, high]` by Gauss-Legendre quadrature."""
    x, w = np.polynomial.legendre.leggauss(nodes)
    return 0.5 * (high - low) * float(np.sum(w * function(0.5 * (high - low) * x + 0.5 * (high + low))))


def _quadrature_masses(shape):
    """Each component's mass by a quadrature of its profile as the model evaluates it, in the model's units (independent
    of the closed forms): the disks in cylinders, the bulge in the bar's scaled spherical coordinates, the halo in
    spheroidal shells (`m` from a hundredth of its core to 1e12 of its break), the nuclear cluster in shells and the
    nuclear disc in cylinders."""
    terms, halo = density.model_terms(shape), density.halo_terms(shape)
    out = {}

    def cylinder(scale_length_pc, scale_height_pc, value, full_height):
        radii, rw = np.polynomial.legendre.leggauss(160)
        radii = 0.5 * 50.0 * scale_length_pc * (radii + 1.0)
        rw = rw * 0.5 * 50.0 * scale_length_pc
        heights, zw = np.polynomial.legendre.leggauss(160)
        heights = 0.5 * full_height * scale_height_pc * (heights + 1.0)
        zw = zw * 0.5 * full_height * scale_height_pc
        r, z = np.meshgrid(radii, heights, indexing="ij")
        found = value(density.component_arrays(r, np.zeros_like(r), z, shape))
        return float(2.0 * 2.0 * math.pi * np.einsum("i,j,ij->", rw * radii, zw, found))

    out["thin_disk"] = cylinder(shape.disk_scale_length_pc, shape.disk_scale_height_pc,
                                lambda v: v.thin_profile * shape.arm_level, 40.0)
    out["thick_disk"] = cylinder(terms["thick_disk_scale_length_pc"], terms["thick_disk_scale_height_pc"],
                                 lambda v: v.thick_disk, 40.0)
    out["nuclear_disc"] = cylinder(terms["nuclear_disc_length_pc"], terms["nuclear_disc_height_pc"],
                                   lambda v: v.nuclear_disc, 40.0)

    # the bulge, in the bar's frame scaled by its three axes: `rho` out to 14, theta over a half, phi over a quarter
    rho, rw = np.polynomial.legendre.leggauss(240)
    rho, rw = 7.0 * (rho + 1.0), rw * 7.0
    theta, tw = np.polynomial.legendre.leggauss(60)
    theta, tw = 0.25 * math.pi * (theta + 1.0), tw * 0.25 * math.pi
    phi, pw = np.polynomial.legendre.leggauss(8)
    phi, pw = 0.25 * math.pi * (phi + 1.0), pw * 0.25 * math.pi
    r, t, f = np.meshgrid(rho, theta, phi, indexing="ij")
    sx, sy, sz = shape.bulge_scale_radius_pc, terms["bulge_scale_y_pc"], terms["bulge_scale_z_pc"]
    bar_x, bar_y = sx * r * np.sin(t) * np.cos(f), sy * r * np.sin(t) * np.sin(f)
    cos_a, sin_a = terms["bar_cos"], terms["bar_sin"]
    found = density.component_arrays(bar_x * cos_a - bar_y * sin_a, bar_x * sin_a + bar_y * cos_a, sz * r * np.cos(t), shape).bulge
    out["bulge"] = float(8.0 * sx * sy * sz * np.einsum("i,j,k,ijk->", rw, tw, pw, found * r ** 2 * np.sin(t)))

    # the halo, in spheroidal shells of volume 4 pi q m^2 dm (its density depends on `m` alone), m = e^u
    low, high = math.log(0.01 * halo["halo_core_pc"]), math.log(1e12 * halo["halo_break_pc"])
    shells = lambda u: np.array([
        float(density.component_arrays(np.array([m]), np.array([0.0]), np.array([0.0]), shape).halo[0]) * m ** 3
        for m in np.exp(u)])
    out["halo"] = 4.0 * math.pi * halo["halo_axis_ratio"] * (
        _gauss(shells, low, math.log(halo["halo_core_pc"]), 200)
        + _gauss(shells, math.log(halo["halo_core_pc"]), math.log(halo["halo_break_pc"]), 400)
        + _gauss(shells, math.log(halo["halo_break_pc"]), high, 2000))

    # the nuclear cluster, softened: shells out to 1e9 pc
    soft = terms["nuclear_cluster_softening_pc"]
    if terms["nuclear_cluster_mass"]:
        cluster = lambda r: 4.0 * math.pi * r ** 2 * density._nuclear_cluster_array(r, terms)
        out["nuclear_cluster"] = (_gauss(cluster, 0.0, soft, 100)
                                  + _gauss(lambda u: cluster(np.exp(u)) * np.exp(u), math.log(soft), math.log(1e9), 1500))
    else:
        out["nuclear_cluster"] = 0.0
    return out


def test_the_closed_form_masses_agree_with_a_quadrature_of_each_profile_to_1e_6():
    shape = density.build_galaxy_shape(**BASE, nuclear_amplitude=1.0)
    numeric = _quadrature_masses(shape)
    for component in density.density_components(shape):
        wanted = component.mass * (component.parameters["softened_mass_fraction"] if component.name == "nuclear_cluster" else 1.0)
        assert numeric[component.name] == pytest.approx(wanted, rel=1e-6), component.name


def test_the_solar_masses_are_the_components_shares_of_the_stellar_mass():
    masses = density.component_masses(SHAPE)
    assert tuple(masses) == NAMES
    assert sum(masses.values()) == pytest.approx(SHAPE.structure.stellar_mass_solar, rel=1e-12)
    assert masses["thin_disk"] > masses["bulge"] > masses["nuclear_cluster"] > 0.0


def test_the_integrated_masses_the_old_checks_use_still_agree_with_the_closed_forms():
    by_name = _by_name(SHAPE)
    numeric = density.integrated_masses(SHAPE, steps=80)
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


IDENTITY = os.path.join(os.path.dirname(__file__), "fixtures", "density_identity.json")
REFERENCE = galaxyStructure.REFERENCE_STRUCTURE
IDENTITY_SHAPES = {
    "reference": lambda: mathcheck._reference_galaxy_shape(),
    "nuclear": lambda: density.build_galaxy_shape(**BASE, nuclear_amplitude=1.0),
    "round_bulge": lambda: density.build_galaxy_shape(
        **BASE, nuclear_amplitude=0.5,
        structure=REFERENCE._replace(bar_present=False, bar_angle_deg=63.0, bar_axis_ratio_y=0.8, bar_axis_ratio_z=0.6)),
    "no_thick": lambda: density.build_galaxy_shape(**BASE, structure=REFERENCE._replace(thick_mass_ratio=0.0)),
    "no_bulge": lambda: density.build_galaxy_shape(**{**BASE, "bulge_amplitude": 0.0}, arm_level=0.8),
    "four_arms": lambda: density.build_galaxy_shape(
        **{**BASE, "arm_count": 4, "pitch_angle_rad": math.radians(22.0), "disk_scale_length_pc": 3800.0,
           "disk_scale_height_pc": 450.0, "bulge_scale_radius_pc": 500.0}, nuclear_amplitude=1.0),
}


def _close(got, wanted):
    got, wanted = np.asarray(got, dtype=float), np.asarray(wanted, dtype=float)
    assert np.all((got == wanted) | (np.abs(got - wanted) <= 1e-12 * np.abs(wanted))), (got, wanted)


@pytest.mark.parametrize("name", sorted(IDENTITY_SHAPES))
def test_the_spiral_written_as_components_gives_the_density_it_gave_before(name):
    """The fixture holds the density, the population split, the array forms and the halo of six shapes at 24 points, made
    by the code before GEN.379's rewrite (`relative_density`, `population_densities`, `population_densities_array`,
    `total_density_array`, `halo_density`, `total_mass`): the same to 1e-12 (it is the same to the last digit). The fixture was made again with GEN.929
    and GEN.931, which move two things of the old code's output and nothing else: every density is the old one times one
    number per shape (0.9515 to 0.9540: the model's units per solar mass follow the primary mass function's systems per solar
    mass, GEN.929; checked point by point before the fixture was rewritten, and the total mass is unchanged), and the split of
    a place into populations follows the arm contrast of each (GEN.931; it moves the mix of a place, never its total)."""
    with open(IDENTITY, encoding="utf-8") as handle:
        data = json.load(handle)
    shape, wanted = IDENTITY_SHAPES[name](), data["shapes"][name]
    points = [tuple(point) for point in data["points"]]
    x, y, z = (np.array([point[axis] for point in points]) for axis in range(3))
    _close([density.relative_density(point, shape) for point in points], wanted["relative_density"])
    for point, split in zip(points, wanted["population_densities"]):
        got = density.population_densities(point, shape)
        for population, value in split.items():
            _close(got[population], value)
    arrays = density.population_densities_array(x, y, z, shape)
    for population, values in wanted["population_densities_array"].items():
        _close(arrays[population], values)
    _close(density.total_density_array(x, y, z, shape), wanted["total_density_array"])
    _close([density.halo_density(point, shape) for point in points], wanted["halo_density"])
    _close(density.total_mass(shape), wanted["total_mass"])


def test_the_density_is_the_sum_of_the_components_and_a_zero_component_adds_nothing():
    shape = density.build_galaxy_shape(**BASE, nuclear_amplitude=1.0)
    for point in ((0.0, 0.0, 0.0), (500.0, 200.0, 30.0), (9000.0, -4000.0, 700.0), (2.0, -1.0, 0.5)):
        values = density.component_values(point, shape)
        total = sum(getattr(values, name) for name in density.COMPONENT_NAMES)
        assert shape.k_norm * total == pytest.approx(density.relative_density(point, shape), rel=1e-12)
    none = density.build_galaxy_shape(**{**BASE, "bulge_amplitude": 0.0}, structure=REFERENCE._replace(thick_mass_ratio=0.0))
    values = density.component_values((600.0, 100.0, 20.0), none)
    assert values.bulge == 0.0 and values.thick_disk == 0.0 and values.nuclear_cluster == 0.0 and values.nuclear_disc == 0.0
    arrays = density.component_arrays(np.array([600.0, 5000.0]), np.array([100.0, 0.0]), 20.0, none)
    assert not arrays.bulge.any() and not arrays.thick_disk.any()


def _row(**stored):
    row = dict(family="spiral", bar_present=None, bar_axis_ratio_b=None, bar_axis_ratio_c=None, bar_angle_rad=None,
               thick_disk_amplitude=None, spheroid_amplitude=0.0, spheroid_sersic_n=0.0, spheroid_effective_radius_pc=0.0,
               spheroid_axis_ratio_b=1.0, spheroid_axis_ratio_c=1.0, spheroid_orientation_1_rad=0.0,
               spheroid_orientation_2_rad=0.0, spheroid_orientation_3_rad=0.0, body_clump_share=0.0)
    return {**row, **stored}


def test_a_null_bar_and_thick_disk_are_the_structure_draws_own_and_a_number_pins_them():
    shape = density.build_galaxy_shape(**BASE)
    read = density.with_family_columns(shape, _row())
    assert read.structure == shape.structure and read.family == "spiral" and read.bar_present is None
    assert density.thick_disk_amplitude(read) == density.thick_disk_amplitude(shape)
    assert density.component_values((900.0, 100.0, 40.0), read) == density.component_values((900.0, 100.0, 40.0), shape)

    pinned = density.with_family_columns(shape, _row(
        bar_present=0, bar_axis_ratio_b=0.5, bar_axis_ratio_c=0.25, bar_angle_rad=math.radians(40.0),
        thick_disk_amplitude=0.07))
    assert pinned.bar_present is False and pinned.structure.bar_present is False
    assert (pinned.structure.bar_axis_ratio_y, pinned.structure.bar_axis_ratio_z) == (0.5, 0.25)
    assert pinned.structure.bar_angle_deg == pytest.approx(40.0)
    assert density.thick_disk_amplitude(pinned) == 0.07
    assert pinned.structure.thick_mass_ratio == pytest.approx(
        0.07 * shape.structure.thick_length_ratio ** 2 * shape.structure.thick_height_ratio / shape.arm_level)
    assert density.model_terms(pinned)["thick_disk_amplitude"] == 0.07
    assert density.model_terms(pinned)["bulge_scale_y_pc"] == 0.5 * shape.bulge_scale_radius_pc

    removed = density.with_family_columns(shape, _row(thick_disk_amplitude=0.0))
    assert _by_name(removed)["thick_disk"].mass == 0.0


def test_a_stored_spheroid_is_refused_until_the_spheroid_component_exists():
    shape = density.with_family_columns(density.build_galaxy_shape(**BASE), _row(
        family="elliptical", spheroid_amplitude=2.0, spheroid_sersic_n=4.0, spheroid_effective_radius_pc=3000.0,
        spheroid_orientation_2_rad=0.3, body_clump_share=0.0))
    assert shape.family == "elliptical" and shape.spheroid_orientation_rad == (0.0, 0.3, 0.0)
    with pytest.raises(ValueError, match="Sersic spheroid"):
        density.model_terms(shape)


def test_the_galaxy_maps_density_shape_is_the_one_it_was_before_the_family_numbers():
    """The prisms' shape (`shape_with_terms`) holds the spiral's own numbers and the model's terms, none of the stored
    family numbers: `bar_angle_rad` and `thick_disk_amplitude` in it are the terms the model works out."""
    shape = density.build_galaxy_shape(**BASE)
    terms = density.shape_with_terms(shape)
    for name in ("family", "bar_present", "spheroid_amplitude", "spheroid_orientation_rad", "body_clump_share"):
        assert name not in terms
    assert terms["bar_angle_rad"] == density.model_terms(shape)["bar_angle_rad"]
    assert terms["thick_disk_amplitude"] == density.model_terms(shape)["thick_disk_amplitude"]
    assert terms["k_norm"] == shape.k_norm and terms["disk_scale_length_pc"] == 2600.0
