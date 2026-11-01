"""
Constants table, held against astropy (GEN.66, GEN.211).

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
from planetgen.physics import constants as constants_module

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


ASTROPY_VALUES = [
    # (name, astropy accessor, relative tolerance)
    ("EARTH_GRAVITY", lambda c, u: c.g0.value, 0),
    ("AU_M", lambda c, u: c.au.value, 0),
    ("LIGHTYEAR_M", lambda c, u: u.lyr.to(u.m), 0),
    ("PARSEC_M", lambda c, u: c.pc.value, 1e-15),
    ("G", lambda c, u: c.G.value, 1e-12),
    ("SPEED_OF_LIGHT_M_S", lambda c, u: c.c.value, 0),
    ("R", lambda c, u: c.R.value, 1e-12),
    ("BOLTZMANN", lambda c, u: c.k_B.value, 1e-12),
    ("REDUCED_PLANCK", lambda c, u: c.hbar.value, 1e-12),
    ("STEFAN_BOLTZMANN_CONSTANT", lambda c, u: c.sigma_sb.value, 1e-12),
    ("SOLAR_MASS_TO_KG", lambda c, u: c.M_sun.value, 1e-12),
    ("SOLAR_LUMINOSITY", lambda c, u: c.L_sun.value, 1e-12),
    ("EARTH_MASS_TO_KG", lambda c, u: c.M_earth.value, 1e-12),
    ("JUPITER_MASS_TO_KG", lambda c, u: c.M_jup.value, 1e-12),
    ("JUPITER_RADIUS_KM", lambda c, u: c.R_jup.to_value(u.km), 0),
    # [R4361] (astropy 8 defaults to [R3117], which moves this by 1.35e-9; GEN.211).
    ("HYDROGEN_ATOM_MASS_KG", lambda c, u: c.m_p.value + c.m_e.value, 2e-9),
    ("SOLAR_RADIUS_M", lambda c, u: c.R_sun.value, 0),
    ("SECONDS_PER_YEAR", lambda c, u: u.year.to(u.s), 0),
    ("ELECTRON_VOLT_J", lambda c, u: c.e.value, 0),
]


@pytest.mark.parametrize("name, accessor, tolerance", ASTROPY_VALUES)
def test_constant_table_matches_astropy(name, accessor, tolerance):
    """GEN.211: the table in `physics/constants.py` is held against astropy (a test-only dependency)."""
    pytest.importorskip("astropy")
    from astropy import constants, units
    assert getattr(constants_module, name) == pytest.approx(float(accessor(constants, units)), rel=tolerance, abs=0)


def test_importing_the_constants_does_not_import_astropy():
    import subprocess
    import sys
    out = subprocess.run([sys.executable, "-c", "import sys, planetgen.physics.constants; print('astropy' in sys.modules)"],
                         capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"
