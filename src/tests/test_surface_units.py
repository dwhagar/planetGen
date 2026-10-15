# tests/test_surface_units.py

"""
Surface conditions in customary units alongside metric (Boss, 2026-10-01):
`planetgen.util.format.format_temperature_k` shows K with °C and °F, and
`format_pressure_pa` a pressure on the nPa to TPa ladder with atm and psi, bar
from 10 atm to 10,000 atm and mbar from 1 Pa to 0.01 atm; `format_gravity_g`
shows m/s² with g.
"""

import pytest

from planetgen.physics import constants as pc
from planetgen.util.format import (
    format_gravity_g, format_gravity_g_short, format_pressure_atm, format_pressure_pa, format_temperature_k,
)


@pytest.mark.parametrize("kelvin,expected", [
    (288.0, "288 K (15 °C, 59 °F)"),
    (737.0, "737 K (464 °C, 867 °F)"),
    (273.15, "273 K (0 °C, 32 °F)"),
    (276.0, "276 K (2.9 °C, 37 °F)"),
    (250.0, "250 K (-23 °C, -9.7 °F)"),
    (233.15, "233 K (-40 °C, -40 °F)"),
    (165.0, "165 K (-108 °C, -163 °F)"),
    (1500.0, "1,500 K (1,227 °C, 2,240 °F)"),
    (2.7, "2.7 K (-270 °C, -455 °F)"),
])
def test_temperature_in_k_c_and_f(kelvin, expected):
    assert format_temperature_k(kelvin) == expected


def test_temperature_none_is_an_en_dash():
    assert format_temperature_k(None) == "–"


@pytest.mark.parametrize("pascals,expected", [
    (0.0, "0 Pa (0 atm, 0 psi)"),
    (2.5e-13, "2.50 × 10⁻⁴ nPa (2.47 × 10⁻¹⁸ atm, 3.63 × 10⁻¹⁷ psi)"),
    (0.5, "500 mPa (4.93 × 10⁻⁶ atm, 7.25 × 10⁻⁵ psi)"),
    (610.0, "610 Pa (0.00602 atm, 0.0885 psi, 6.1 mbar)"),
    (999.0, "999 Pa (0.00986 atm, 0.145 psi, 9.99 mbar)"),
    (1000.0, "1 kPa (0.00987 atm, 0.145 psi, 10 mbar)"),
    (101_325.0, "101 kPa (1 atm, 14.7 psi)"),
    (9.2e6, "9.2 MPa (90.8 atm, 1,334 psi, 92 bar)"),
    (9.2e8, "920 MPa (9,080 atm, 133,435 psi, 9,200 bar)"),
    (3e11, "300 GPa (2.96 × 10⁶ atm, 4.35 × 10⁷ psi)"),
])
def test_pressure_ladder_with_atm_and_psi(pascals, expected):
    assert format_pressure_pa(pascals) == expected


def test_pressure_wrappers_and_none():
    assert format_pressure_atm(1.0) == "101 kPa (1 atm, 14.7 psi)"
    assert format_pressure_pa(None) == "–"
    assert format_pressure_atm(None) == "–"


def test_customary_constants():
    assert pc.STANDARD_ATMOSPHERE_PA == 101_325.0
    assert pc.PSI_PA == pytest.approx(4.4482216152605 / 0.0254 ** 2, rel=1e-12)
    assert pc.CELSIUS_ZERO_K == 273.15


@pytest.mark.parametrize("kelvin,expected", [
    (0.0, "0 K (-273 °C, -460 °F)"),
    (0.5, "500 mK (-273 °C, -459 °F)"),
    (1e-8, "10 nK (-273 °C, -460 °F)"),
])
def test_temperature_under_a_kelvin_takes_an_si_prefix(kelvin, expected):
    assert format_temperature_k(kelvin) == expected


def test_a_star_temperature_is_kelvin_alone():
    assert format_temperature_k(5778.0, companions=False) == "5,778 K"
    assert format_temperature_k(None, companions=False) == "–"


@pytest.mark.parametrize("gravity_g,expected,short", [
    (1.0, "9.81 m/s² (1 g)", "1 g"),
    (0.378, "3.71 m/s² (0.378 g)", "0.378 g"),
    (12.7, "125 m/s² (12.7 g)", "12.7 g"),
    (0.0, "0 m/s² (0 g)", "0 g"),
])
def test_gravity_in_m_s2_with_g(gravity_g, expected, short):
    assert format_gravity_g(gravity_g) == expected
    assert format_gravity_g_short(gravity_g) == short


def test_gravity_none_is_an_en_dash():
    assert format_gravity_g(None) == "–" and format_gravity_g_short(None) == "–"
