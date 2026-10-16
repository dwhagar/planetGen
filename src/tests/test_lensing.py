"""GEN.154: the Einstein radius of a compact object."""

import math

import pytest

from planetgen.physics import constants, lensing

AU_M = constants.AU_TO_M


def test_a_ten_solar_mass_hole_at_a_kiloparsec_has_an_einstein_radius_of_about_eight_au():
    mass = 10 * constants.SOLAR_MASS_TO_KG
    source = 8.0 * 1000 * constants.PARSEC_M
    radius = lensing.einstein_radius_m(mass, 1000 * constants.PARSEC_M, source)
    expected = math.sqrt(4 * constants.G * mass / constants.SPEED_OF_LIGHT_M_S ** 2
                         * 1000 * constants.PARSEC_M * (1 - 1000 / 8000))
    assert radius == pytest.approx(expected)
    assert 7 < radius / AU_M < 9  # the design note's 1 to 10 AU for 100 pc to 1 kpc


def test_the_radius_grows_as_the_square_root_of_the_mass():
    near, far = 1000 * constants.PARSEC_M, 8000 * constants.PARSEC_M
    one = lensing.einstein_radius_m(constants.SOLAR_MASS_TO_KG, near, far)
    ten = lensing.einstein_radius_m(10 * constants.SOLAR_MASS_TO_KG, near, far)
    assert ten / one == pytest.approx(math.sqrt(10))


def test_a_lens_not_between_the_observer_and_the_source_has_none():
    mass = constants.SOLAR_MASS_TO_KG
    assert lensing.einstein_radius_m(mass, 0.0, 10.0) is None
    assert lensing.einstein_radius_m(mass, 10.0, 10.0) is None
    assert lensing.einstein_radius_m(mass, 20.0, 10.0) is None
    assert lensing.einstein_radius_m(0.0, 5.0, 10.0) is None


def test_against_the_galactic_centre_the_lens_distance_is_sols_less_the_radius():
    mass = constants.SOLAR_MASS_TO_KG
    source = constants.GALACTIC_CENTER_DISTANCE_LY * constants.LY_TO_M
    radius_pc = source / constants.PARSEC_M - 1000.0  # 1 kpc from Sol
    assert lensing.einstein_radius_against_galactic_centre_m(mass, radius_pc) == pytest.approx(
        lensing.einstein_radius_m(mass, 1000 * constants.PARSEC_M, source))
    assert lensing.einstein_radius_against_galactic_centre_m(mass, source / constants.PARSEC_M + 1.0) is None
