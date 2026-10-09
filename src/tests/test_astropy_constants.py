"""
Constants from astropy (GEN.66).

`planetgen.physics.constants` takes its measured and defined values from
`astropy.constants` and `astropy.units`. These tests hold them against the
hand-set values the program used before (copied below), each within a stated
relative tolerance, so a change to astropy's tables that moves a figure the
program depends on is noticed. The tolerance of each is the size of the
difference the move to astropy was known to make.
"""
import math

import pytest

from planetgen.physics import constants

# (name, value used before GEN.66, relative tolerance)
TODAYS_VALUES = [
    ("EARTH_GRAVITY", 9.807, 1e-4),
    ("AU_M", 149_597_870_700, 0),
    ("LIGHTYEAR_M", 9_460_730_472_580_800, 0),
    ("PARSEC_M", 3.085677581491367e16, 1e-15),
    ("G", 6.6743e-11, 1e-12),
    ("SPEED_OF_LIGHT_M_S", 299_792_458, 0),
    ("SPEED_OF_LIGHT_KMS", 299_792.458, 1e-15),
    ("R", 8.314, 1e-4),
    ("BOLTZMANN", 1.381e-23, 5e-4),
    ("REDUCED_PLANCK", 1.054571817e-34, 1e-8),
    ("STEFAN_BOLTZMANN_CONSTANT", 5.67e-8, 1e-3),
    ("SOLAR_MASS_TO_KG", 1.989e30, 5e-4),
    ("SOLAR_LUMINOSITY", 3.82e26, 3e-3),
    ("EARTH_MASS_TO_KG", 5.972e24, 1e-4),
    ("JUPITER_MASS_TO_KG", 1.898e27, 1e-3),
    ("JUPITER_RADIUS_KM", 71492, 0),
    ("HYDROGEN_ATOM_MASS_KG", 1.6735575e-27, 1e-4),
    ("SOLAR_RADIUS_M", 6.957e8, 0),
    ("SOLAR_ESCAPE_VELOCITY", 617.7 * 1000, 5e-4),
    ("SECONDS_PER_YEAR", 365.25 * 24 * 3600, 0),
]


@pytest.mark.parametrize("name, before, tolerance", TODAYS_VALUES)
def test_constant_agrees_with_the_value_used_before(name, before, tolerance):
    assert getattr(constants, name) == pytest.approx(before, rel=tolerance, abs=0)


def test_exact_definitions_stay_integers():
    assert isinstance(constants.AU_M, int)
    assert isinstance(constants.LIGHTYEAR_M, int)
    assert isinstance(constants.SPEED_OF_LIGHT_M_S, int)


def test_derived_figures_follow_their_sources():
    assert constants.LIGHTYEAR_M == constants.SPEED_OF_LIGHT_M_S * 365.25 * 86400
    assert constants.PARSEC_M == pytest.approx(constants.AU_M * 648000 / math.pi, rel=1e-14)
    assert constants.AU_TO_KM == constants.AU_M / 1000
    assert constants.SOLAR_ESCAPE_VELOCITY == pytest.approx(
        math.sqrt(2 * constants.G * constants.SOLAR_MASS_TO_KG / constants.SOLAR_RADIUS_M), rel=1e-15)
