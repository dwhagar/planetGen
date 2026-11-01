#!/usr/bin/env python3
"""
The active-Sun factor of the NH3 photodissociation rate (GEN.794), from a measured cross section and measured solar spectra.

    python scripts/nh3_activity_factor.py /path/to/folder

The folder holds four files (downloaded by hand, never committed; `h5py` is needed to read them):

* `NH3_partial_cross_sections.h5`: the photodissociation cross section of NH3 (from [R3889], Wu 2007,), https://home.strw.leidenuniv.nl/~ewine/photo/data/photo_data/NH3_partial_cross_sections.h5
* `sorcesols_v16_1nm.h5`, `susim_v22_1nm.h5`, `uarssols_v18_1nm.h5`: the daily 1 nm solar irradiance of SORCE SOLSTICE (2003 to 2018), UARS SUSIM
  (1991 to 2005) and UARS SOLSTICE (1991 to 2001), M. DeLand, https://sbuv.gsfc.nasa.gov/solar/compositev2/

For each instrument the rate J is the integral of the cross section times the photon flux over 120 to 230 nm, for every year of its record (the
median of the days of the year that carry no fill flag). The quiet Sun of an instrument is its weakest year and the active Sun its strongest, and
the factor is the strongest year's rate over the weakest year's. Each instrument is paired with itself, never with another: the three instruments
differ in absolute scale by several percent at 200 nm (Yeo and others 2015), which would otherwise pass for solar activity.

Limits: below 120 nm none of the instruments measure (Lyman-alpha at 121.6 nm is inside the range; the EUV below 120 nm, which varies more, carries
a few percent of the rate), so the factor is a little low; the instruments' cycle amplitudes near 200 nm disagree (SORCE SOLSTICE more than the
SUSIM), which is why the three are all kept; the cross section is that of 295 K; the [R2040] "active Sun" is stronger than any
cycle maximum measured (for H2O 1.68 against 1.33 by this method), so a factor on their scale would be higher.
"""

import argparse
import os
import sys

import numpy as np

PLANCK_J_S = 6.62607015e-34
LIGHT_M_S = 2.99792458e8
LOW_NM = 120.0
HIGH_NM = 230.0
MINIMUM_DAYS = 30
FILES = {
    "cross section": "NH3_partial_cross_sections.h5",
    "SORCE SOLSTICE": "sorcesols_v16_1nm.h5",
    "UARS SUSIM": "susim_v22_1nm.h5",
    "UARS SOLSTICE": "uarssols_v18_1nm.h5",
}


def photodissociation_rate_per_s(wavelength_nm, irradiance_mw_m2_nm, sigma_wavelength_nm, sigma_cm2):
    """
    The photodissociation rate at 1 AU, per second: the cross section averaged over each 1 nm bin times the photon flux of the bin.

    Args:
        wavelength_nm (ndarray): The bin centres of the spectrum, nm, 1 nm apart.
        irradiance_mw_m2_nm (ndarray): The irradiance of each bin, mW per m2 per nm.
        sigma_wavelength_nm, sigma_cm2 (ndarray): The cross section and its wavelengths, nm.
    """
    wavelength = np.asarray(wavelength_nm, dtype=float)
    keep = (wavelength >= LOW_NM) & (wavelength < HIGH_NM) & np.isfinite(irradiance_mw_m2_nm)
    photons_cm2_s_nm = np.asarray(irradiance_mw_m2_nm, dtype=float)[keep] * 1.0e-3 * (wavelength[keep] * 1.0e-9) / (PLANCK_J_S * LIGHT_M_S) / 1.0e4
    sigma_wavelength = np.asarray(sigma_wavelength_nm, dtype=float)
    sigma = np.asarray(sigma_cm2, dtype=float)
    mean_sigma = np.zeros(int(keep.sum()))
    for index, centre in enumerate(wavelength[keep]):
        inside = (sigma_wavelength >= centre - 0.5) & (sigma_wavelength < centre + 0.5)
        mean_sigma[index] = sigma[inside].mean() if inside.any() else 0.0
    return float(np.sum(mean_sigma * photons_cm2_s_nm))


def yearly_rates(years, wavelength_nm, irradiance, fill_flag, sigma_wavelength_nm, sigma_cm2):
    """`{year: rate}` from the median of the unflagged days of each year that has `MINIMUM_DAYS` of them."""
    years = np.asarray(years)
    rates = {}
    for year in sorted(set(int(y) for y in years)):
        days = np.nonzero(years == year)[0]
        good = np.where(fill_flag[:, days] == 0, irradiance[:, days], np.nan)
        if np.sum(np.any(np.isfinite(good), axis=0)) < MINIMUM_DAYS:
            continue
        spectrum = np.nanmedian(good, axis=1)
        if not np.isfinite(spectrum[(wavelength_nm >= LOW_NM) & (wavelength_nm < HIGH_NM)]).all():
            continue
        rates[year] = photodissociation_rate_per_s(wavelength_nm, spectrum, sigma_wavelength_nm, sigma_cm2)
    return rates


def activity_factor(rates):
    """`(factor, quiet_year, active_year)` of one instrument: the strongest year's rate over the weakest year's."""
    quiet = min(rates, key=rates.get)
    active = max(rates, key=rates.get)
    return rates[active] / rates[quiet], quiet, active


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("folder", help="the folder holding the four files named in the module docstring")
    folder = parser.parse_args().folder
    try:
        import h5py
    except ImportError:
        sys.exit("h5py is needed to read the files (pip install h5py)")
    cross = h5py.File(os.path.join(folder, FILES["cross section"]), "r")
    sigma_wavelength, sigma = cross["wavelength"][:], cross["photodissociation"][:]
    factors = {}
    for name in ("SORCE SOLSTICE", "UARS SUSIM", "UARS SOLSTICE"):
        data = h5py.File(os.path.join(folder, FILES[name]), "r")
        rates = yearly_rates(data["Year"][:], data["WavelengthScale"][:], data["Irradiance"][:, :len(data["Year"])], data["FillFlag"][:, :len(data["Year"])],
                             sigma_wavelength, sigma)
        factor, quiet, active = activity_factor(rates)
        factors[name] = factor
        print(f"{name}: quiet {quiet} {rates[quiet]:.4e} per s, active {active} {rates[active]:.4e} per s, factor {factor:.4f}")
    print(f"mean of the instruments {np.mean(list(factors.values())):.4f}, range {min(factors.values()):.4f} to {max(factors.values()):.4f}")


if __name__ == "__main__":
    main()
