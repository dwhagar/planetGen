# planetgen/physics/stellar_library.py

"""
The light of a star by band from model atmospheres and measured chromospheres (GEN.970, parent GEN.745).

`physics.stellar_spectrum` (GEN.443) gives a star a blackbody of its effective temperature with its ultraviolet cut by a seeded
share. This module answers from the spectra instead (`datasets.stellar_spectra`: ATLAS9 and PHOENIX-ACES photospheres) and adds
what the photosphere lacks, the chromosphere and corona, from the MUSCLES Treasury Survey. `physics.activity.uv_shares_of(star)` (which calls `uv_shares`) is what generation
calls (GEN.971): `physics.activity` stores its UV-C share as the star's `uv_c_share`, and `physics.body_radiation` takes the three
shares of a lit body's star, in place of the seeded share and the family built on it.

What the library does against the sky (the tests, `test_stellar_library_real.py`):

- The Sun (5777 K, log g 4.44): the photosphere gives the shares of the blackbody's power in UV-C, UV-B and UV-A of 0.29, 0.71 and 0.77;
  the measured values (ATLAS-1, Thuillier et al. 2004, `stellar_spectrum.UV_B_EXPONENT`) are 0.256, 0.625 and 0.787.
- The eleven M and K stars of MUSCLES (France et al. 2016, ApJ 820, 89, Tables 1 and 3, `MUSCLES_STARS`): the photosphere gives
  the near-ultraviolet (170 to 320 nm) to within 0.3 dex of the measured for all eleven, and it gives almost nothing in the far
  ultraviolet (91.2 to 170 nm; 1e-8 to 1e-15 of the bolometric power for the M dwarfs against the measured 1e-4.2 to 1e-4.8):
  that light is the chromosphere's and transition region's, and it follows the star's X-ray and extreme-ultraviolet output.
  Across the eleven stars the far-ultraviolet power is the XUV power (5 to 91.2 nm) within a factor of 3 (`FUV_OVER_XUV`, the
  median ratio 1.12).

So a star's band shares are its photosphere's plus an excess in the far ultraviolet that is `FUV_OVER_XUV` times its XUV power over
its luminosity (`physics.activity` gives `l_xuv_w`) and an excess XUV in the ionising band. The excess is zero for a star with no
corona (a white dwarf, a hot star). The near-ultraviolet excess is left out: the library already matches the measured near-ultraviolet.

Where the chromospheric far ultraviolet lies in the generator's bands (`stellar_spectrum.BAND_EDGES_NM`): all of it from 100 to 170 nm,
in UV-C (the 91.2 to 100 nm sliver goes with the ionising light, `FUV_IN_UV_C`).
"""

import functools
import math

import numpy as np

from planetgen.datasets import fetch, stellar_spectra
from planetgen.physics import constants, photometry, stellar_spectrum

MUSCLES_STARS = (
    # name, T_eff K (France et al. 2016 Table 1), then log10 of the share of the bolometric luminosity (Table 3) in
    # XUV (5 to 911 A), FUV (912 to 1700 A, with Lyman alpha) and NUV (1700 to 3200 A)
    ("GJ 1214", 2935, -4.50, -4.45, -4.20),
    ("GJ 876", 3062, -4.04, -4.53, -4.33),
    ("GJ 581", 3295, -4.50, -4.84, -4.18),
    ("GJ 436", 3281, -4.41, -4.52, -4.02),
    ("GJ 176", 3416, -4.46, -4.40, -3.81),
    ("GJ 667C", 3327, -4.21, -4.17, -3.86),
    ("GJ 832", 3816, -4.30, -4.29, -3.60),
    ("HD 85512", 4305, -4.82, -4.45, -3.01),
    ("HD 40307", 4783, -4.57, -4.31, -2.41),
    ("eps Eri", 5162, -4.15, -4.10, -1.97),
    ("HD 97658", 5156, -4.61, -4.36, -2.05),
)
"""tuple: The eleven stars of the MUSCLES Treasury Survey (France et al. 2016, ApJ 820, 89), read from the paper's tables."""

FUV_OVER_XUV = 1.12
"""float: The far-ultraviolet (912 to 1700 A) power of a star's chromosphere over its XUV (5 to 911 A) power: the median of
10 ** (f_FUV - f_XUV) over `MUSCLES_STARS` (0.32 to 1.15 for the 7 M dwarfs, 1.1 to 2.3 for the 4 K stars). [Stand-in SL2, tier T2,
tweak L0: eleven stars, one ratio, no dependence on temperature or age.]"""

BANDS = ("ionising", "uv_c", "uv_b", "uv_a", "visible", "near_infrared", "thermal_infrared")
"""tuple: The bands of `stellar_spectrum.BAND_EDGES_NM`, in order."""

FUV_NM = (91.2, 170.0)
XUV_NM = (0.5, 91.2)
FUV_IN_UV_C = (170.0 - 100.0) / (170.0 - 91.2)
"""float: The share of the chromospheric far ultraviolet (91.2 to 170 nm) that lies in the generator's UV-C (100 to 280 nm), taken as
the share of the band's width (the light is spread evenly in wavelength; Lyman alpha at 121.6 nm sits inside UV-C either way)."""

DEFAULT_LOGG = 4.5


def photosphere_shares(teff_k, logg=DEFAULT_LOGG, lib=None):
    """The share of the bolometric power of a model-atmosphere star in each of `BANDS`, as a dict, or None when the
    effective temperature is outside the library's 2300 to 50000 K."""
    edges = [stellar_spectrum.BAND_EDGES_NM[band] for band in BANDS]
    shares = stellar_spectra.band_fractions(teff_k, logg, edges, lib)
    return None if shares is None else dict(zip(BANDS, shares))


def band_shares(teff_k, logg=DEFAULT_LOGG, xuv_over_lbol=0.0, lib=None):
    """
    The share of a star's bolometric power in each of `BANDS`: its photosphere (`photosphere_shares`) plus the chromospheric
    and coronal far ultraviolet and the XUV of a star whose XUV power over its luminosity is `xuv_over_lbol` (0 for a star with no
    corona). None when the effective temperature is outside the library.
    """
    shares = photosphere_shares(teff_k, logg, lib)
    if shares is None:
        return None
    if xuv_over_lbol > 0.0:
        fuv = FUV_OVER_XUV * xuv_over_lbol
        shares["ionising"] += xuv_over_lbol + fuv * (1.0 - FUV_IN_UV_C)
        shares["uv_c"] += fuv * FUV_IN_UV_C
    return shares


def shares_of_blackbody(teff_k, logg=DEFAULT_LOGG, xuv_over_lbol=0.0, lib=None):
    """
    The UV-C, UV-B and UV-A power of the star as a multiple of its blackbody's in the same band (the `stellar_spectrum.uv_shares`
    quantities, 1 meaning the blackbody's value), or None outside the library. The seeded `uv_c_share` of GEN.443 is replaced by these.
    """
    shares = band_shares(teff_k, logg, xuv_over_lbol, lib)
    if shares is None:
        return None
    result = {}
    for band in ("uv_c", "uv_b", "uv_a"):
        low, high = stellar_spectrum.BAND_EDGES_NM[band]
        result[band] = shares[band] / float(stellar_spectrum.blackbody_fraction(low, high, teff_k))
    return result


FAR_UV_NM = (91.2, 200.0)
"""tuple: The band the photolysis rates of `physics.photolysis` are per unit flux of (91.2 to 200 nm)."""


def far_uv_power_w(teff_k, logg, luminosity_w, xuv_over_lbol=0.0):
    """
    The power a star puts into 91.2 to 200 nm, W (GEN.982): its photosphere's share of the luminosity in the band (the model atmosphere, or the
    blackbody when the library is missing and the run disregards that) plus the chromosphere's far ultraviolet, `FUV_OVER_XUV` times the
    star's XUV power over its luminosity `xuv_over_lbol` (all of it lies under 170 nm, inside the band). Zero for a star with no luminosity.
    """
    share = None
    if not fetch.locate(stellar_spectra.DATASET).use_stand_in("the blackbody share of the star's far ultraviolet"):
        shares = stellar_spectra.band_fractions(teff_k, logg, [FAR_UV_NM], stellar_spectra.library())
        share = None if shares is None else shares[0]
    if share is None:
        share = float(stellar_spectrum.blackbody_fraction(FAR_UV_NM[0], FAR_UV_NM[1], teff_k))
    return luminosity_w * (share + FUV_OVER_XUV * xuv_over_lbol)


def surface_logg(mass_kg, radius_km):
    """log10 of the surface gravity, cm/s^2, of a star of this mass and radius."""
    return math.log10(constants.G * mass_kg / (radius_km * 1.0e3) ** 2 * 100.0)


@functools.lru_cache(maxsize=8192)
def _shares_of_blackbody_cached(teff_k, logg, xuv_over_lbol):
    return shares_of_blackbody(teff_k, logg, xuv_over_lbol, stellar_spectra.library())


def uv_shares(mass_kg, radius_km, teff_k, xuv_over_lbol):
    """
    The UV-C, UV-B and UV-A power of a star of this mass, radius and effective temperature as multiples of its blackbody's
    (`shares_of_blackbody`; `xuv_over_lbol` is the corona's XUV power over the luminosity). A star outside the library's
    temperatures (2300 to 50000 K) takes the blackbody (1 in every band). When the library is not downloaded and the run
    disregards missing data, None: the caller takes the seeded shares (`fetch.stand_ins_used` lists it); any other run raises
    `fetch.DataMissing`.
    """
    if fetch.locate(stellar_spectra.DATASET).use_stand_in("the blackbody with the seeded UV shares of stellar_spectrum"):
        return None
    shares = _shares_of_blackbody_cached(round(float(teff_k), 1), round(surface_logg(mass_kg, radius_km), 3), float(xuv_over_lbol))
    return {"uv_c": 1.0, "uv_b": 1.0, "uv_a": 1.0} if shares is None else shares


def logg_of(star):
    """log10 of the surface gravity, cm/s^2, of an object with `mass` (kg) and `radius` (km), `DEFAULT_LOGG` when it lacks either."""
    mass, radius = getattr(star, "mass", None), getattr(star, "radius", None)
    return surface_logg(mass, radius) if mass and radius else DEFAULT_LOGG


_MISSING = "the blackbody with the BT-Settl offsets of physics/photometry.py"


def visible_spectrum(teff_k, logg, wavelengths_nm):
    """
    The relative spectral power per nm of a star at `wavelengths_nm` (nm, ascending, inside 380 to 780), read off its model atmosphere
    (`stellar_spectra.spectrum`: the share of the bolometric power in each log-spaced bin, divided by the bin's width), or None
    outside the library's temperatures. None too when the library is not downloaded and the run disregards missing data
    (`fetch.stand_ins_used` lists it); any other run raises `fetch.DataMissing`.
    """
    if fetch.locate(stellar_spectra.DATASET).use_stand_in(_MISSING):
        return None
    return _visible_cached(round(float(teff_k), 1), round(float(logg), 3), tuple(float(nm) for nm in wavelengths_nm))


@functools.lru_cache(maxsize=2048)
def _visible_cached(teff_k, logg, wavelengths_nm):
    shares = stellar_spectra.spectrum(teff_k, logg, stellar_spectra.library())
    if shares is None:
        return None
    edges = stellar_spectra.EDGES_NM
    density = shares / np.diff(edges)
    centres = np.sqrt(edges[:-1] * edges[1:])
    return np.interp(np.asarray(wavelengths_nm), centres, density) / density.max()


@functools.lru_cache(maxsize=2048)
def _colour(teff_k, logg):
    grid = tuple(float(nm) for nm in np.arange(photometry.WAVELENGTH_MIN_NM, photometry.WAVELENGTH_MAX_NM + 1.0))
    density = visible_spectrum(teff_k, logg, grid)
    if density is None:
        return None
    return photometry.colour_from_xyz(photometry.spectrum_to_xyz([float(v) for v in density], list(grid)))


def star_colour(teff_k, logg=DEFAULT_LOGG):
    """The `photometry.Colour` of a star as seen from space (absolute, no white adaptation) from its model atmosphere, or, outside the
    library's temperatures or for a run that disregards a missing library, `photometry.star_colour`'s blackbody with the BT-Settl offsets."""
    found = _colour(round(float(teff_k), 1), round(float(logg), 3))
    return photometry.star_colour(teff_k) if found is None else found
