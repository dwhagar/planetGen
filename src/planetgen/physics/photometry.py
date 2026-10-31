# planetgen/physics/photometry.py

"""
The photometry library, first part (GEN.416): a star's colour from its temperature.

A star's light is a Planck spectrum of its effective temperature. This module takes a spectrum to a colour: the
spectrum into the CIE 1931 tristimulus values X, Y and Z, into linear sRGB (D65), the brightest channel scaled to 1
(a hue swatch, not a brightness), the sRGB transfer curve on top for the 8-bit colour and the hex string. The same
chain is what the atmosphere colour work (GEN.418, GEN.419) and the view from a planet (VIEW.11) run a spectrum
through, so every step takes a spectrum and a wavelength grid, not only a blackbody.

* **The CIE curves** are the Wyman, Sloan and Shirley (2013, J. Computer Graphics Techniques 2(2)) multi-lobe fit of
  the 1931 two-degree observer (docs/design/orbits/sky-view.md section 4). Against the tabulated curves the fit is within
  0.015, 0.007 and 0.024 of peaks 1.06, 1.00 and 1.78, and the chromaticity of a blackbody is within 0.001 of the
  table of docs/design/composition/atmosphere-composition-parts/02-gas-optics-uv-visible.md section 2.9 from 3000 K up (0.0025 at
  2000 K), the 8-bit colour within 2 of it.
* **Hotter than 4500 K the star is a blackbody.** Cooler stars are not: TiO and VO bands cut the red-optical flux.
  Section 2.9 gives the offset of the BT-Settl atmosphere's chromaticity from the blackbody's at nine temperatures
  (`COOL_OFFSETS`); a cooler star takes its blackbody chromaticity plus that offset, and a seeded spread of 0.02 in
  x and 0.01 in y that the caller draws. The offset is the 4000 K one tapering to nothing at 4500 K, so the colour
  does not step there, and below 2300 K it stays at the 2300 K one (the grid ends there).
* **Absolute colour.** No adaptation to the star: the white point is D65 and a cool star is orange.

* **The bolometric correction** (VIEW.9, `bolometric_correction_v`) is Flower's (1996) polynomial fit of BC_V against log Teff in three
  ranges, with the coefficients as Torres (2010, AJ 140, 1158, arXiv:1008.3913, Table 1) reprinted them (the original had misprints). Its
  zero point is the one the fit carries, BC_V of the Sun -0.08 and so M_bol of the Sun 4.74 (`SUN_ABSOLUTE_BOLOMETRIC_MAG`), which the
  caller must use with it (Torres's point: the two go together). It replaces the blackbody table of sky-view.md section 2.4.

The functions are pure and use only operations whose result is fixed (`math.fsum` for the sums), so a colour is the
same on every machine. Nothing here is stored: a colour is worked out from the stored temperature when it is shown.
"""

import math
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

from planetgen.physics import constants

PLANCK_J_S = 2 * math.pi * constants.REDUCED_PLANCK
"""float: The Planck constant, J s (exact in the SI)."""

SECOND_RADIATION_CONSTANT_M_K = PLANCK_J_S * constants.SPEED_OF_LIGHT_M_S / constants.BOLTZMANN
"""float: c2 = h c / k, metre kelvin (0.01438776877)."""

WAVELENGTH_MIN_NM = 380
WAVELENGTH_MAX_NM = 780
"""int: The visible range the tristimulus integral runs over, nm, at 1 nm."""

BLACKBODY_ABOVE_K = 4500.0
"""float: Hotter than this a star's colour is a blackbody's."""

COOL_SPREAD_X = 0.02
COOL_SPREAD_Y = 0.01
"""float: The half-width of the seeded spread of a cool star's chromaticity (section 2.9)."""

COOL_OFFSETS = (
    # Teff K, dx, dy: BT-Settl (solar metallicity, log g 5) minus blackbody, CIE 1931 2-degree, section 2.9.
    (2300.0, -0.0100, -0.0222),
    (2500.0, -0.0495, -0.0154),
    (2800.0, -0.0273, -0.0104),
    (3000.0, -0.0208, -0.0079),
    (3300.0, -0.0018, -0.0029),
    (3500.0, +0.0026, -0.0039),
    (4000.0, +0.0180, -0.0023),
)
"""tuple: (temperature K, dx, dy) rows of the BT-Settl chromaticity offset from a blackbody's."""

COOL_TAPER_FROM_K = 4000.0
"""float: The offset of the last `COOL_OFFSETS` row falls linearly to zero at `BLACKBODY_ABOVE_K` from here."""

_XYZ_TO_LINEAR_SRGB = (
    (3.2404542, -1.5371385, -0.4985314),
    (-0.9692660, 1.8760108, 0.0415560),
    (0.0556434, -0.2040259, 1.0572252),
)
"""tuple: The D65 matrix from XYZ to linear sRGB (IEC 61966-2-1)."""


class NoVisibleLight(ValueError):
    """A spectrum holds no light between 380 and 780 nm (a black hole's or a cold remnant's), so it has no colour."""


@dataclass(frozen=True)
class Colour:
    """
    A colour: the chromaticity, the linear and the encoded sRGB.

    Attributes:
        x, y (float): CIE 1931 chromaticity.
        linear (tuple): Linear sRGB `(r, g, b)` with the brightest channel 1 and negatives clipped to 0.
        rgb8 (tuple): The 8-bit sRGB `(r, g, b)`.
        hex (str): `#rrggbb` of `rgb8`.
    """

    x: float
    y: float
    linear: Tuple[float, float, float]
    rgb8: Tuple[int, int, int]
    hex: str


def planck_relative(wavelength_nm, temperature_k):
    """
    A blackbody's spectral radiance at `wavelength_nm`, on a scale that is the same at every wavelength (only the
    ratios matter to a colour): `1 / (lam^5 (exp(c2 / (lam T)) - 1))` with `lam` in micrometres.
    """
    if not temperature_k > 0 or not math.isfinite(temperature_k):
        raise ValueError(f"a temperature must be positive and finite, not {temperature_k!r}")
    lam_um = wavelength_nm * 1e-3
    x = SECOND_RADIATION_CONSTANT_M_K * 1e6 / (lam_um * temperature_k)
    if x > 700.0:
        return 0.0
    return 1.0 / (lam_um ** 5 * math.expm1(x))


def gives_visible_light(temperature_k):
    """Whether a blackbody this hot gives any light between 380 and 780 nm (a black hole's temperature gives none, and a spectrum
    with no light has no colour); the long end of the range is the last to go dark."""
    return temperature_k > 0 and math.isfinite(temperature_k) and planck_relative(780.0, temperature_k) > 0.0


def _lobe(lam, mu, sigma_low, sigma_high):
    t = (lam - mu) / (sigma_low if lam < mu else sigma_high)
    return math.exp(-0.5 * t * t)


def cie_xyz_bar(wavelength_nm):
    """The CIE 1931 colour matching functions `(xbar, ybar, zbar)` at `wavelength_nm` (Wyman, Sloan and Shirley 2013)."""
    lam = float(wavelength_nm)
    xbar = 1.056 * _lobe(lam, 599.8, 37.9, 31.0) + 0.362 * _lobe(lam, 442.0, 16.0, 26.7) - 0.065 * _lobe(lam, 501.1, 20.4, 26.2)
    ybar = 0.821 * _lobe(lam, 568.8, 46.9, 40.5) + 0.286 * _lobe(lam, 530.9, 16.3, 31.1)
    zbar = 1.217 * _lobe(lam, 437.0, 11.8, 36.0) + 0.681 * _lobe(lam, 459.0, 26.0, 13.8)
    return xbar, ybar, zbar


_GRID_NM = tuple(range(WAVELENGTH_MIN_NM, WAVELENGTH_MAX_NM + 1))
_BARS = tuple(cie_xyz_bar(lam) for lam in _GRID_NM)


def spectrum_to_xyz(spectrum, wavelengths_nm: Optional[Sequence[float]] = None):
    """
    The tristimulus values `(X, Y, Z)` of a spectrum, normalised so `Y` is 1.

    Args:
        spectrum (sequence or callable): Spectral power per wavelength: values on `wavelengths_nm`, or a function
            of the wavelength in nm called on the 1 nm grid from 380 to 780 nm when `wavelengths_nm` is `None`.
        wavelengths_nm (sequence): Ascending wavelengths of the values, nm; `None` for the 1 nm visible grid.

    Returns:
        tuple: `(X, Y, Z)` with `Y == 1`.

    Raises:
        ValueError: If the spectrum has no light in the visible range.
    """
    if wavelengths_nm is None:
        grid = _GRID_NM
        bars = _BARS
        values = [spectrum(lam) for lam in grid] if callable(spectrum) else list(spectrum)
        if len(values) != len(grid):
            raise ValueError(f"a spectrum without wavelengths must have {len(grid)} values (380 to 780 nm at 1 nm)")
    else:
        grid = list(wavelengths_nm)
        values = list(spectrum)
        if len(values) != len(grid) or len(grid) < 2:
            raise ValueError("a spectrum and its wavelengths must have the same length, at least 2")
        bars = tuple(cie_xyz_bar(lam) for lam in grid)
    totals = []
    for k in range(3):
        pieces = [
            0.5 * (grid[i + 1] - grid[i]) * (values[i] * bars[i][k] + values[i + 1] * bars[i + 1][k])
            for i in range(len(grid) - 1)
        ]
        totals.append(math.fsum(pieces))
    x, y, z = totals
    if not y > 0:
        raise NoVisibleLight("the spectrum has no light between 380 and 780 nm")
    return x / y, 1.0, z / y


def chromaticity(xyz):
    """The CIE 1931 `(x, y)` of tristimulus values."""
    x, y, z = xyz
    total = x + y + z
    return x / total, y / total


def xyz_from_chromaticity(x, y):
    """Tristimulus values with `Y == 1` for the chromaticity `(x, y)`."""
    return x / y, 1.0, (1.0 - x - y) / y


def xyz_to_linear_srgb(xyz):
    """Linear sRGB (D65) of tristimulus values; may lie outside 0 to 1."""
    return tuple(math.fsum(m * v for m, v in zip(row, xyz)) for row in _XYZ_TO_LINEAR_SRGB)


def srgb_encode(linear):
    """The sRGB transfer curve of a linear value between 0 and 1."""
    linear = min(max(linear, 0.0), 1.0)
    return 12.92 * linear if linear <= 0.0031308 else 1.055 * linear ** (1.0 / 2.4) - 0.055


def colour_from_xyz(xyz):
    """
    The `Colour` of tristimulus values: negatives clipped, the brightest channel scaled to 1 (a hue, not a
    brightness), then the sRGB transfer curve.
    """
    rgb = [max(c, 0.0) for c in xyz_to_linear_srgb(xyz)]
    peak = max(rgb)
    if not peak > 0:
        raise ValueError("the colour is outside the sRGB gamut on every channel")
    linear = tuple(c / peak for c in rgb)
    rgb8 = tuple(int(round(255.0 * srgb_encode(c))) for c in linear)
    x, y = chromaticity(xyz)
    return Colour(x, y, linear, rgb8, "#{:02x}{:02x}{:02x}".format(*rgb8))


def blackbody_xyz(temperature_k):
    """The tristimulus values (`Y == 1`) of a Planck spectrum of `temperature_k`."""
    return spectrum_to_xyz(lambda lam: planck_relative(lam, temperature_k))


def blackbody_colour(temperature_k):
    """The `Colour` of a blackbody of `temperature_k`."""
    return colour_from_xyz(blackbody_xyz(temperature_k))


def cool_offset(temperature_k):
    """
    The chromaticity offset `(dx, dy)` of a cool star's atmosphere from a blackbody of the same temperature, and the
    weight `(1 at 4000 K and below, 0 at 4500 K and above)` that the seeded spread takes too.
    """
    if temperature_k >= BLACKBODY_ABOVE_K:
        return 0.0, 0.0, 0.0
    weight = 1.0 if temperature_k <= COOL_TAPER_FROM_K else (BLACKBODY_ABOVE_K - temperature_k) / (BLACKBODY_ABOVE_K - COOL_TAPER_FROM_K)
    rows = COOL_OFFSETS
    if temperature_k <= rows[0][0]:
        return rows[0][1], rows[0][2], weight
    if temperature_k >= rows[-1][0]:
        return rows[-1][1] * weight, rows[-1][2] * weight, weight
    for (t0, dx0, dy0), (t1, dx1, dy1) in zip(rows, rows[1:]):
        if t0 <= temperature_k <= t1:
            f = (temperature_k - t0) / (t1 - t0)
            return dx0 + f * (dx1 - dx0), dy0 + f * (dy1 - dy0), weight
    raise AssertionError("unreachable")  # pragma: no cover


def star_colour(temperature_k, spread_x=0.0, spread_y=0.0):
    """
    The colour of a star of effective temperature `temperature_k`.

    A blackbody from 4500 K up; a cooler star takes the BT-Settl offset (`cool_offset`) and a spread of
    `spread_x * 0.02` in x and `spread_y * 0.01` in y, scaled by the same taper, where `spread_x` and `spread_y`
    are numbers from -1 to 1 the caller draws from the star's seeded stream (0 for the mean colour).

    Args:
        temperature_k (float): The effective temperature, kelvin.
        spread_x, spread_y (float): The star's seeded draws, -1 to 1; ignored at 4500 K and above.

    Returns:
        Colour: The star's colour as seen from space, absolute (no white adaptation).
    """
    xyz = blackbody_xyz(temperature_k)
    dx, dy, weight = cool_offset(temperature_k)
    if weight == 0.0:
        return colour_from_xyz(xyz)
    x, y = chromaticity(xyz)
    x += dx + weight * spread_x * COOL_SPREAD_X
    y += dy + weight * spread_y * COOL_SPREAD_Y
    return colour_from_xyz(xyz_from_chromaticity(x, y))


FLOWER_COEFFICIENTS = (
    (3.70, (-0.190537291496456e+05, 0.155144866764412e+05, -0.421278819301717e+04, 0.381476328422343e+03)),
    (3.90, (-0.370510203809015e+05, 0.385672629965804e+05, -0.150651486316025e+05, 0.261724637119416e+04, -0.170623810323864e+03)),
    (math.inf, (-0.118115450538963e+06, 0.137145973583929e+06, -0.636233812100225e+05, 0.147412923562646e+05, -0.170587278406872e+04,
                0.788731721804990e+02)),
)
"""tuple: `(upper log Teff, coefficients a, b, c, ...)` of Flower's (1996) fit `BC_V = a + b (log Teff) + c (log Teff)^2 + ...` for log Teff below
3.70, from 3.70 to 3.90 and above 3.90, as Torres (2010, Table 1) corrects them. [Verified, read from the paper.]"""

BOLOMETRIC_CORRECTION_RANGE_K = (3500.0, 50000.0)
"""tuple: The effective temperatures, K, `bolometric_correction_v` answers for. Below 3500 K the polynomial runs away (-24 at 2000 K) and Torres says
the fit "breaks down completely for the M dwarfs"; the hot end is where the fit's data stop being dense. The caller decides what a cooler
or hotter star gets; none is invented here."""

SUN_ABSOLUTE_BOLOMETRIC_MAG = 4.74
"""float: M_bol of the Sun that goes with Flower's zero point (Torres 2010: BC_V of the Sun -0.08, M_bol 4.74, which with V = -26.76 and the 1 au
distance modulus 31.572 gives BC_V of the Sun -0.068)."""


def bolometric_correction_v(temperature_k):
    """
    The V-band bolometric correction `BC_V` (magnitudes, `M_V = M_bol - BC_V`) of a star of effective temperature `temperature_k`: Flower's (1996)
    fit with Torres's (2010) coefficients. About -0.08 for the Sun, +0.03 near 7000 K, -1.9 at 20,000 K and -4.2 at 45,000 K.

    Raises:
        ValueError: If `temperature_k` is outside `BOLOMETRIC_CORRECTION_RANGE_K` (or not finite).
    """
    low, high = BOLOMETRIC_CORRECTION_RANGE_K
    if not (isinstance(temperature_k, (int, float)) and low <= temperature_k <= high):
        raise ValueError(f"the bolometric correction is fitted from {low:g} to {high:g} K, not {temperature_k!r}")
    log_t = math.log10(temperature_k)
    for upper, coefficients in FLOWER_COEFFICIENTS:
        if log_t < upper:
            return sum(coefficient * log_t ** power for power, coefficient in enumerate(coefficients))
