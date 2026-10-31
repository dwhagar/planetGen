# tests/test_nh3_activity_factor.py

"""GEN.794: the arithmetic of scripts/nh3_activity_factor.py on synthetic spectra (the real files are never in the repository)."""

import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import nh3_activity_factor as script  # noqa: E402

WAVELENGTH = np.arange(120.5, 300.0, 1.0)
SIGMA_WAVELENGTH = np.arange(100.0, 240.0, 0.1)
SIGMA = np.full_like(SIGMA_WAVELENGTH, 1.0e-17)


def _flat(irradiance):
    return np.full_like(WAVELENGTH, irradiance)


def test_the_rate_is_the_cross_section_times_the_photon_flux_over_120_to_230_nm():
    rate = script.photodissociation_rate_per_s(WAVELENGTH, _flat(1.0), SIGMA_WAVELENGTH, SIGMA)
    kept = WAVELENGTH[(WAVELENGTH >= 120.0) & (WAVELENGTH < 230.0)]
    flux = 1.0e-3 * kept * 1.0e-9 / (script.PLANCK_J_S * script.LIGHT_M_S) / 1.0e4
    assert rate == pytest.approx(1.0e-17 * flux.sum())


def test_doubling_the_light_doubles_the_rate_and_a_dark_bin_adds_nothing():
    rate = script.photodissociation_rate_per_s(WAVELENGTH, _flat(1.0), SIGMA_WAVELENGTH, SIGMA)
    assert script.photodissociation_rate_per_s(WAVELENGTH, _flat(2.0), SIGMA_WAVELENGTH, SIGMA) == pytest.approx(2.0 * rate)
    dark = _flat(1.0)
    dark[(WAVELENGTH > 200.0) & (WAVELENGTH < 217.0)] = 0.0
    assert script.photodissociation_rate_per_s(WAVELENGTH, dark, SIGMA_WAVELENGTH, SIGMA) < rate


def test_the_factor_pairs_the_strongest_year_with_the_weakest_of_the_same_record():
    years = np.repeat([2000, 2001, 2002], 40)
    scale = np.repeat([1.0, 1.08, 1.04], 40)
    irradiance = np.outer(np.ones_like(WAVELENGTH), scale)
    flags = np.zeros_like(irradiance)
    rates = script.yearly_rates(years, WAVELENGTH, irradiance, flags, SIGMA_WAVELENGTH, SIGMA)
    factor, quiet, active = script.activity_factor(rates)
    assert (quiet, active) == (2000, 2001)
    assert factor == pytest.approx(1.08)


def test_flagged_days_and_thin_years_are_left_out():
    years = np.repeat([2000, 2001], [60, 10])
    irradiance = np.ones((len(WAVELENGTH), 70))
    irradiance[:, 40:60] = 9.0
    flags = np.zeros_like(irradiance)
    flags[:, 40:60] = 1
    rates = script.yearly_rates(years, WAVELENGTH, irradiance, flags, SIGMA_WAVELENGTH, SIGMA)
    assert list(rates) == [2000]
    assert rates[2000] == pytest.approx(script.photodissociation_rate_per_s(WAVELENGTH, _flat(1.0), SIGMA_WAVELENGTH, SIGMA))
