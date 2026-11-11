"""
GEN.1525 (row 9 of the GEN.945 audit): the slope of the giant planet mass function of `tuning` held against the occurrence fit the repository holds
(`physics/data/planet_occurrence.json`, band `cold_giants_cumming`, [R0960], READ).
"""

import re

import pytest

from planetgen import tuning
from planetgen.physics import planet_occurrence

_BAND = "cold_giants_cumming"


def _held_slope():
    definition = planet_occurrence.BANDS[_BAND].definition
    found = re.search(r"dN/dlnM proportional to M to the (-?[0-9.]+) \+- ([0-9.]+)", definition)
    assert found, "the mass function of the Cumming band is no longer stated in its definition"
    return -float(found.group(1)), float(found.group(2))     # the band writes the exponent (M to the -0.31); the tuning writes the slope of M^-slope


def test_the_giant_mass_function_slope_is_the_one_the_occurrence_fit_holds():
    slope, error = _held_slope()
    assert error == pytest.approx(0.2)
    assert tuning.GIANT_MASS_FUNCTION_SLOPE == pytest.approx(slope, abs=1e-9)        # 0.31; dN/dlnM and dN/dlogM have the same exponent
