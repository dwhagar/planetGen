# tests/test_enclosed_mass.py

"""GEN.370: Hill spheres take the galaxy's own mass inside the object's own radius (`galaxy.enclosed_mass`), not the
whole halo's mass, and the stored shape the unit runs in."""

import math

import numpy as np
import pytest

from planetgen.galaxy import density, disk_history, enclosed_mass, keepout, mass_chain, seed as galaxySeed, shape_draw
from planetgen.generation.binary import BinaryStarProxy
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import constants
from planetgen.physics.orbits import calculate_hill_sphere
from planetgen.physics.units import pc_to_ly

SEED = bytes(range(1, 17))
RADII_PC = np.array([2.0, 30.0, 400.0, 2000.0, 8000.0, 30000.0, 2.0e5])


@pytest.fixture(scope="module")
def shape():
    return shape_draw.default_shape(SEED)


@pytest.fixture(scope="module")
def model(shape):
    return enclosed_mass.mass_model(shape, SEED)


def _doubled(shape):
    """The same stars spread over twice the size (every length of the shape doubled)."""
    return shape._replace(
        disk_scale_length_pc=2.0 * shape.disk_scale_length_pc, disk_scale_height_pc=2.0 * shape.disk_scale_height_pc,
        bulge_scale_radius_pc=2.0 * shape.bulge_scale_radius_pc,
        spiral_reference_radius_pc=2.0 * shape.spiral_reference_radius_pc)


def test_the_mass_grows_with_the_radius_and_ends_at_the_galaxys_whole_mass(shape, model):
    enclosed = model.enclosed_solar(np.geomspace(0.1, 1.0e6, 400))
    assert np.all(np.diff(enclosed) >= 0.0)
    stellar_mass, halo_mass = disk_history.stellar_and_halo_masses(shape)
    whole = model.enclosed_solar(1.0e7)
    assert whole == pytest.approx(halo_mass + model.black_hole_solar, rel=0.05)
    assert model.stars_solar(1.0e7) == pytest.approx(stellar_mass, rel=0.05)


def test_at_the_centre_the_black_hole_is_most_of_what_is_inside(shape, model):
    black_hole = density.central_black_hole_mass_solar(shape)
    assert model.enclosed_solar(0.5) > black_hole
    assert model.enclosed_solar(0.5) == pytest.approx(black_hole, rel=0.5)


def test_inside_the_galaxy_less_than_the_halo_mass_is_enclosed(shape, model):
    _stellar, halo_mass = disk_history.stellar_and_halo_masses(shape)
    inside = model.enclosed_solar(RADII_PC[:5])
    assert np.all(inside < 0.5 * halo_mass)
    # the old central mass, the halo's whole mass, overstated it by this factor at the disk's own scale length
    at_scale_length = model.enclosed_solar(shape.disk_scale_length_pc)
    assert halo_mass / at_scale_length > 5.0


def test_the_dark_halo_follows_nfw_and_the_stars_follow_the_model(shape, model):
    c = model.concentration
    assert model.dark_halo_solar(0.0) == 0.0
    assert model.dark_halo_solar(model.virial_radius_pc) == pytest.approx(model.dark_solar)
    assert model.dark_halo_solar(10.0 * model.virial_radius_pc) == pytest.approx(model.dark_solar)
    half = model.dark_halo_solar(model.virial_radius_pc / 2.0)
    expected = model.dark_solar * (math.log(1.0 + c / 2.0) - (c / 2.0) / (1.0 + c / 2.0)) / (math.log1p(c) - c / (1.0 + c))
    assert half == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("radius_pc", [300.0, 2500.0, 12000.0])
def test_the_stars_inside_a_sphere_match_a_monte_carlo_of_the_density_model(shape, model, radius_pc):
    rng = np.random.default_rng(370)
    points = rng.uniform(-1.0, 1.0, size=(1_500_000, 3))
    points = points[np.einsum("ij,ij->i", points, points) <= 1.0] * radius_pc
    volume = 4.0 / 3.0 * math.pi * radius_pc ** 3
    raw = density.total_density_array(points[:, 0], points[:, 1], points[:, 2], shape) / shape.k_norm
    expected = raw.mean() * volume / density.model_terms(shape)["mass_per_solar"]
    assert model.stars_solar(radius_pc) == pytest.approx(expected, rel=0.06)


def test_a_galaxy_twice_the_size_changes_the_mass_and_so_the_hill_radius(shape):
    radius = 3.0 * shape.disk_scale_length_pc
    small = enclosed_mass.enclosed_mass_kg(shape, SEED, radius)
    large = enclosed_mass.enclosed_mass_kg(_doubled(shape), SEED, radius)
    assert large < small
    star_mass = constants.SOLAR_MASS_TO_KG
    hill_small = keepout.galactic_hill_radius_km(star_mass, radius, small)
    hill_large = keepout.galactic_hill_radius_km(star_mass, radius, large)
    assert hill_large > hill_small
    assert hill_large / hill_small == pytest.approx((small / large) ** (1.0 / 3.0), rel=1e-6)


def test_a_star_and_a_pair_take_their_hill_spheres_from_the_mass_inside_their_own_radius(shape):
    star = Star(SystemConfig())
    with galaxySeed.seeded(SEED, "sector", "1/1/1", shape=shape):
        results = {}
        for distance_pc in (50.0, 700.0, 9000.0):
            distance_ly = pc_to_ly(distance_pc)
            perimeter = star.calculate_system_perimeter(distance_ly)
            pair = BinaryStarProxy._calculate_system_perimeter_static(star.mass, distance_ly)
            central = enclosed_mass.enclosed_mass_kg(shape, SEED, distance_pc)
            expected = calculate_hill_sphere(distance_pc * constants.PARSEC_M, star.mass, central) / constants.AU_TO_M
            assert perimeter == pytest.approx(expected, rel=1e-9)
            assert pair == pytest.approx(perimeter, rel=1e-9)
            results[distance_pc] = central
    assert results[50.0] < results[700.0] < results[9000.0]


def test_the_shape_the_unit_is_given_is_the_one_read(shape):
    radius = 2000.0
    drawn = shape_draw.default_shape(SEED)
    with galaxySeed.seeded(SEED, "sector", "1/1/1"):
        assert enclosed_mass.current_enclosed_mass_kg(radius) == enclosed_mass.enclosed_mass_kg(drawn, SEED, radius)
    other = _doubled(shape)
    with galaxySeed.seeded(SEED, "sector", "1/1/1", shape=other):
        assert galaxySeed.current_galaxy_shape() is other
        assert enclosed_mass.current_enclosed_mass_kg(radius) == enclosed_mass.enclosed_mass_kg(other, SEED, radius)
        assert enclosed_mass.current_enclosed_mass_kg(radius) != enclosed_mass.enclosed_mass_kg(drawn, SEED, radius)
    assert galaxySeed.current_galaxy_shape() is None


def test_outside_any_unit_the_fallback_seeds_galaxy_stands_in():
    radius = 1500.0
    fallback = mass_chain.FALLBACK_SEED
    assert enclosed_mass.current_enclosed_mass_kg(radius) == enclosed_mass.enclosed_mass_kg(
        shape_draw.default_shape(fallback), fallback, radius)


def test_two_galaxies_enclose_different_masses():
    other = bytes(range(21, 37))
    first = enclosed_mass.enclosed_mass_kg(shape_draw.default_shape(SEED), SEED, 2000.0)
    second = enclosed_mass.enclosed_mass_kg(shape_draw.default_shape(other), other, 2000.0)
    assert first != second


def test_the_star_clusters_and_satellites_of_a_shape_are_not_part_of_the_cache_key(shape):
    with_components = shape._replace(components=object())
    assert enclosed_mass.mass_model(with_components, SEED) is enclosed_mass.mass_model(shape, SEED)


def test_the_old_whole_halo_central_mass_is_gone():
    assert not hasattr(mass_chain, "central_mass_kg")
    assert not hasattr(mass_chain, "current_central_mass_kg")
    assert not hasattr(mass_chain.MassChain, "halo_mass_kg")
