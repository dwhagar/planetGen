# planetgen/physics/stellar_spectrum.py

"""
The light a star puts on a planet, band by band (GEN.443; docs/design/radiation-model.md R6 and R17).

The star's spectrum is a blackbody of its effective temperature scaled to the bolometric flux at the planet, with the
ultraviolet cut by a per-star share. The blackbody overpredicts the Sun's measured ultraviolet by 3.9 times in UV-C, 1.6 in
UV-B and 1.27 in UV-A (radiation-model.md section 12; ATLAS-1, [R4354]), and the research has no source for the
ultraviolet against effective temperature and activity. So the shares are an **open, seeded bracket (R17, Boss 9:54 PM EDT
2026-10-17: "UV-C share 0.1 to 1 of the blackbody")**: `uv_c_share` is drawn once per star in `[UV_C_SHARE_RANGE]` and the UV-B and UV-A
shares follow it on one family, `share_b = share_c ** UV_B_EXPONENT`, `share_a = share_c ** UV_A_EXPONENT`, fitted to the Sun's
measured shares (0.256, 0.625, 0.787) as a labelled stand-in: the Sun's values are a **target, never a preset**. Replace the
family when a source for the ultraviolet against T_eff and activity is read (BT-Settl and MUSCLES spectra are reachable).
**Since GEN.971 the seeded share and the family are the fallback only**: a star's shares come from `physics.stellar_library` (model atmospheres
and the MUSCLES corona), and this module's seeded share serves a run that disregards a missing library and a body whose star lacks the
quantities the library needs.

Bands (R6, nm): UV-C 100 to 280, UV-B 280 to 315, UV-A 315 to 400, visible 400 to 700, near infrared 700 to 4000, and the
thermal infrared beyond 4000 nm; below 100 nm is ionising and never reaches the ground (R18).

Units: wavelengths in nm, fluxes in W/m2 (a band) or W/m2/nm (per wavelength).
"""

import math

import numpy as np

from planetgen.physics import constants

BAND_EDGES_NM = {
    "ionising": (0.0, 100.0), "uv_c": (100.0, 280.0), "uv_b": (280.0, 315.0), "uv_a": (315.0, 400.0),
    "visible": (400.0, 700.0), "near_infrared": (700.0, 4000.0), "thermal_infrared": (4000.0, math.inf),
}
"""dict: The band edges of R6 (partly verified: the UV names are second-hand, 700 nm is a convention, 4 micron leaks under 1 percent
for Sun-like and cool-Earth blackbodies)."""

UV_C_SHARE_RANGE = (0.1, 1.0)
"""tuple: R17's bracket for a star's UV-C power as a share of its blackbody's (log-uniform). Open: no source against T_eff."""

UV_B_EXPONENT = 0.345
UV_A_EXPONENT = 0.176
"""float: The family's exponents, `ln(0.625)/ln(0.256)` and `ln(0.787)/ln(0.256)`, from the Sun's measured-to-blackbody shares
in UV-B and UV-A against UV-C (radiation-model.md section 12). A stand-in until a law of T_eff and activity is sourced."""

SOLAR_CONSTANT_PER_LUMINOSITY = 1.0 / (4.0 * math.pi * constants.AU_M ** 2)
"""float: W/m2 at 1 AU per watt of luminosity (the Sun, 3.828e26 W, gives 1361 W/m2)."""

_H = 6.62607015e-34
_C = 299792458.0
_K = constants.BOLTZMANN

_X = np.concatenate([np.geomspace(1e-6, 1.0, 400, endpoint=False), np.geomspace(1.0, 120.0, 2400)])
_PLANCK = _X ** 3 / np.expm1(_X)
_CUMULATIVE = np.concatenate([[0.0], np.cumsum(0.5 * (_PLANCK[1:] + _PLANCK[:-1]) * np.diff(_X))])
_TOTAL = math.pi ** 4 / 15.0


def _fraction_longer_than(x):
    """The share of a blackbody's power at wavelengths longer than the one whose `x = h c / (lambda k T)` is `x`."""
    x = np.asarray(x, dtype=float)
    inside = np.interp(x, _X, _CUMULATIVE)                         # 0 below the first node, the whole above the last
    tiny = np.where(x < _X[0], x ** 3 / 3.0, 0.0)                  # Rayleigh-Jeans: the integrand is x^2 for small x
    return (inside + tiny) / _TOTAL


def blackbody_fraction(low_nm, high_nm, temperature_k):
    """The share of a blackbody's total power between `low_nm` and `high_nm` (scalars or arrays; `high_nm` may be `inf`)."""
    if temperature_k <= 0:
        raise ValueError("a blackbody needs a temperature above zero")
    scale = _H * _C / (_K * temperature_k) * 1e9           # nm
    low = np.asarray(low_nm, dtype=float)
    high = np.asarray(high_nm, dtype=float)
    x_low = np.where(low > 0.0, scale / np.maximum(low, 1e-12), np.inf)
    x_high = np.where(np.isfinite(high), scale / np.maximum(high, 1e-12), 0.0)
    return _fraction_longer_than(x_low) - _fraction_longer_than(x_high)


def uv_shares(uv_c_share):
    """
    The UV-C, UV-B and UV-A shares of the blackbody for a star whose UV-C share is `uv_c_share`.

    Returns:
        dict: `uv_c`, `uv_b`, `uv_a` (1 means the blackbody's value).
    """
    share = min(max(float(uv_c_share), 1e-6), 1.0)
    return {"uv_c": share, "uv_b": share ** UV_B_EXPONENT, "uv_a": share ** UV_A_EXPONENT}


def draw_uv_c_share(rng):
    """One star's UV-C share, log-uniform in `UV_C_SHARE_RANGE` (R17), from `rng` (an object with `random()`)."""
    low, high = UV_C_SHARE_RANGE
    return math.exp(math.log(low) + rng.random() * (math.log(high) - math.log(low)))


def irradiance_w_m2(edges_nm, temperature_k, flux_w_m2, uv_c_share=1.0, shares=None):
    """
    The star's flux in each cell of a wavelength grid at the planet.

    Args:
        edges_nm (ndarray): The `n + 1` cell edges, ascending, nm.
        temperature_k (float): The star's effective temperature.
        flux_w_m2 (float): The bolometric flux at the planet, W/m2 (`L / 4 pi d^2`).
        uv_c_share (float, optional): The star's UV-C share of its blackbody (`uv_shares` gives the other two); 1 is the blackbody.
        shares (dict, optional): The `uv_c`, `uv_b` and `uv_a` shares themselves (`stellar_library.uv_shares`), in place of `uv_c_share`'s family.

    Returns:
        ndarray: W/m2 in each cell, with the ultraviolet bands cut by their shares.
    """
    edges = np.asarray(edges_nm, dtype=float)
    low, high = edges[:-1], edges[1:]
    power = flux_w_m2 * blackbody_fraction(low, high, temperature_k)
    shares = uv_shares(uv_c_share) if shares is None else shares
    middle = 0.5 * (low + high)
    factor = np.ones_like(middle)
    for band, share in shares.items():
        lo, hi = BAND_EDGES_NM[band]
        factor = np.where((middle >= lo) & (middle < hi), share, factor)
    return power * factor


def bolometric_flux_w_m2(luminosity_w, distance_au):
    """The flux of a star of `luminosity_w` at `distance_au`, W/m2."""
    return luminosity_w * SOLAR_CONSTANT_PER_LUMINOSITY / distance_au ** 2


def band_flux_w_m2(edges_nm, cell_values, band):
    """The sum of per-cell fluxes `cell_values` (W/m2 per cell) over the cells whose centre lies in `band` (a key of `BAND_EDGES_NM`)."""
    edges = np.asarray(edges_nm, dtype=float)
    middle = 0.5 * (edges[:-1] + edges[1:])
    lo, hi = BAND_EDGES_NM[band]
    return float(np.sum(np.asarray(cell_values)[(middle >= lo) & (middle < hi)]))
