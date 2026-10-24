# planetgen/physics/atmosphere_colour.py

"""
The atmosphere colour of a world with a ground (GEN.418; docs/design/composition/atmosphere-colour.md sections 10 to 12 and 14.1).

The colour is the absolute colour of the light the air scatters down to the ground (sky and atmosphere colour are the same thing,
Boss 2026-10-16): the star's spectrum through the Rayleigh depth of the mixture, the ozone's Chappuis absorption, the haze, dust and
cloud, and the ground's bounce, run through the CIE 1931 curves and the sRGB matrix with the brightest channel at 1 (a hue, not a
brightness, and not adapted to the star: a red dwarf's air is cream, not blue). The result also carries the sky's luminance as a share
of the star's, since the swatch hides it, the colour of the sun disc, and a word.

The function is pure: every stored input (ozone column, aerosol depths, cloud cover, ground albedo) is an argument, so the callers that
read the columns of DB.68 and the ground of GEN.441 and GEN.462 come later without touching this module. An argument left out is the
clear, bare-ground case.

**One spectrum from two models** (design section 12, items 5 and 6). Thin air is the single-scatter zenith radiance of the design's
worked example, `L = F0 (w p) mu0 [exp(-tau/mu0) - exp(-tau/mu)] / (4 pi (mu0 - mu))`, with the sun `mu0` from the zenith, the view at
the zenith (the scattering angle is the sun's zenith angle), the Rayleigh phase 3/4 (1 + cos^2) mixed with a Henyey-Greenstein phase for
each aerosol by its scattering depth. Thick air is the delta-Eddington two-stream diffuse flux at the ground (`physics/two_stream.py`), which
is where single scatter saturates to red. Between a total depth of 0.2 and 0.4 at 550 nm the two are blended so no colour steps (smoothstep,
weight on the log of the ratio of the two). The ozone and the ground's bounce enter as the ratio of two-stream runs with them and without
them, so they apply the same way in both. The sky's luminance is always the two-stream diffuse flux over the star's flux on a horizontal
surface at the top.

**The cool star.** A star hotter than 4500 K is a blackbody. A cooler one is a blackbody bent by `exp(a t + b t^2)` (`t = (lambda - 580 nm)
/ 200`) until its chromaticity has the BT-Settl offset of `photometry.cool_offset` plus the caller's seeded spread, so a red dwarf's light
is the right colour with a smooth spectrum.

Limits: a plane-parallel atmosphere (use `mu0` of 0.2 or more), the single-scatter regime under-reflects the Monte Carlo flux by 10 to 21 percent
in brightness (the hue is what the blend keeps), the cloud is a layer at the ground's air (cloud optics are GEN.414's), and the giants are
GEN.419's (they have a reference level, not a ground).
"""

import functools
import json
import math
import os
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from planetgen.physics import air_column, atmospheric_radiation as radiation, elements, photometry, two_stream
from planetgen.util import detmath

WAVELENGTHS_NM = np.arange(380.0, 781.0, 5.0)
"""ndarray: The grid of the colour integral, nm (cell centres)."""

EDGES_NM = np.arange(377.5, 782.6, 5.0)
"""ndarray: The cell edges of `WAVELENGTHS_NM`, nm."""

DEFAULT_MU0 = math.cos(math.radians(30.0))
"""float: The design's worked-example sun, 30 degrees from the zenith."""

THIN_BELOW_TAU = 0.2
THICK_ABOVE_TAU = 0.4
"""float: The total depth at 550 nm below which the single-scatter hue is used alone, and above which the two-stream alone."""

AEROSOL_DOMINATES_SHIFT = 0.03
"""float: The distance in CIE xy that the aerosols move the colour by when the label calls it aerosol-dominated (about a word)."""

DOBSON_CM2 = 2.687e16
"""float: Molecules per square centimetre in one Dobson unit."""

CHAPPUIS_FILE = os.path.join(elements.DATA_DIRECTORY, "ozone_chappuis.json")
"""str: The ozone cross sections of the Chappuis band on this grid (`scripts/build_ozone_chappuis.py`)."""

_REFERENCE_NM = 550.0
_TILT_CENTRE_NM = 580.0
_TILT_WIDTH_NM = 200.0
_YBAR = np.array([photometry.cie_xyz_bar(nm)[1] for nm in WAVELENGTHS_NM])
_STEP = 5.0


@functools.lru_cache(maxsize=None)
def chappuis_cross_section_cm2():
    """The ozone cross section on `WAVELENGTHS_NM`, cm2 per molecule (295 K, Brion-team data of the MPI-Mainz atlas; DATA-LICENSE.md)."""
    with open(CHAPPUIS_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    return np.interp(WAVELENGTHS_NM, document["wavelength_nm"], document["sigma_cm2"])


@dataclass(frozen=True)
class Aerosol:
    """
    A haze, dust or cloud layer: its optical depth at 550 nm and how it scatters.

    Attributes:
        tau_550 (float): The vertical extinction optical depth at 550 nm.
        omega (float or tuple): The single-scattering albedo: one number, or `((nm, omega), ...)` ascending, held at the ends.
        g (float): The asymmetry factor of the scattering (0 isotropic, 0.85 cloud droplets, 0.67 dust).
        exponent (float): The depth goes as `lambda^-exponent` (0 grey, 1 to 2 aggregates, 4 for molecules).
    """

    tau_550: float
    omega: object = 1.0
    g: float = 0.0
    exponent: float = 0.0

    def __post_init__(self):
        if not self.tau_550 >= 0.0 or not math.isfinite(self.tau_550):
            raise ValueError(f"an aerosol's optical depth must be zero or more, not {self.tau_550!r}")
        if not -1.0 < self.g < 1.0:
            raise ValueError(f"an aerosol's asymmetry factor must be inside (-1, 1), not {self.g!r}")
        values = [v for _nm, v in self.omega] if isinstance(self.omega, tuple) else [self.omega]
        if not all(0.0 <= v <= 1.0 for v in values):
            raise ValueError("a single-scattering albedo must lie between 0 and 1")

    def optics(self, wavelengths_nm):
        """The depth, single-scattering albedo and asymmetry at each wavelength, as arrays."""
        tau = self.tau_550 * (_REFERENCE_NM / wavelengths_nm) ** self.exponent
        if isinstance(self.omega, tuple):
            nm = [row[0] for row in self.omega]
            omega = np.interp(wavelengths_nm, nm, [row[1] for row in self.omega])
        else:
            omega = np.full_like(tau, float(self.omega))
        return tau, omega, np.full_like(tau, float(self.g))


@dataclass(frozen=True)
class AtmosphereColour:
    """
    The colour of a world's air.

    Attributes:
        colour (Colour): The sky's hue, absolute, the brightest channel 1.
        word (str): A short word for it (`colour_word`).
        luminance_ratio (float): The sky's luminance as a share of the star's on a horizontal surface at the top of the air.
        sun (Colour or None): The colour of the sun disc through the air, `None` when no direct light reaches the ground.
        sun_luminance_ratio (float): The sun disc's luminance over the star's, 0 to 1.
        tau_550 (float): The vertical extinction depth of the air at 550 nm (gas and aerosols, clear sky).
        tau_rayleigh_550 (float): The Rayleigh part of it.
        regime (str): `thin` (single scatter), `blend` or `thick` (two-stream).
        aerosol_shift (float): The distance in CIE xy between the colour and the colour of the same air with no aerosol.
        aerosol_dominated (bool): The aerosols move the colour by `AEROSOL_DOMINATES_SHIFT` or more.
    """

    colour: photometry.Colour
    word: str
    luminance_ratio: float
    sun: Optional[photometry.Colour]
    sun_luminance_ratio: float
    tau_550: float
    tau_rayleigh_550: float
    regime: str
    aerosol_shift: float
    aerosol_dominated: bool

    @property
    def hex(self):
        return self.colour.hex

    @property
    def label(self):
        """The word, with a note when aerosols decide the colour."""
        return f"{self.word} (aerosol-dominated)" if self.aerosol_dominated else self.word


def _hue_saturation(rgb8):
    r, g, b = (c / 255.0 for c in rgb8)
    top, bottom = max(r, g, b), min(r, g, b)
    span = top - bottom
    if span == 0.0:
        return 0.0, 0.0
    if top == r:
        hue = ((g - b) / span) % 6.0
    elif top == g:
        hue = (b - r) / span + 2.0
    else:
        hue = (r - g) / span + 4.0
    return hue * 60.0, span / top


_HUES = ((15.0, "red"), (45.0, "orange"), (70.0, "yellow"), (165.0, "green"), (200.0, "cyan"), (255.0, "blue"), (290.0, "violet"),
         (340.0, "magenta"), (360.0, "red"))
_PALE_WORDS = {"red": "pink", "orange": "peach", "yellow": "cream"}


def colour_word(rgb8):
    """
    A short word for a hue swatch (the brightest channel at 1, so only the hue and saturation count): white, a creamy or near white,
    a pale tint (`pink`, `peach`, `cream`, `pale blue`), the hue itself, or a `deep` one.
    """
    hue, saturation = _hue_saturation(rgb8)
    name = next(word for limit, word in _HUES if hue < limit)
    if saturation < 0.04:
        return "white"
    if saturation < 0.12:
        return "cream" if name in ("red", "orange", "yellow") else "near white"
    if saturation < 0.5:
        return _PALE_WORDS.get(name, f"pale {name}")
    return f"deep {name}" if saturation >= 0.8 else name


def _xyz(spectrum):
    return photometry.spectrum_to_xyz([float(v) for v in spectrum], [float(nm) for nm in WAVELENGTHS_NM])


def _luminance(spectrum):
    """The integral of a spectrum against the luminance curve (nm times the spectrum's unit)."""
    products = [float(v) * float(y) for v, y in zip(spectrum, _YBAR)]
    return math.fsum(0.5 * _STEP * (a + b) for a, b in zip(products, products[1:]))


def star_spectrum(temperature_k, spread_x=0.0, spread_y=0.0):
    """
    The star's relative spectral power on `WAVELENGTHS_NM`: Planck from 4500 K up; cooler, the Planck spectrum bent by `exp(a t + b t^2)`
    until its chromaticity is the blackbody's plus the BT-Settl offset and the seeded spread (`spread_x`, `spread_y` from -1 to 1).
    """
    base = np.array([photometry.planck_relative(float(nm), temperature_k) for nm in WAVELENGTHS_NM])
    dx, dy, weight = photometry.cool_offset(temperature_k)
    if weight == 0.0:
        return base
    x0, y0 = photometry.chromaticity(_xyz(base))
    target = np.array([x0 + dx + weight * spread_x * photometry.COOL_SPREAD_X, y0 + dy + weight * spread_y * photometry.COOL_SPREAD_Y])
    t = (WAVELENGTHS_NM - _TILT_CENTRE_NM) / _TILT_WIDTH_NM

    def chromaticity_of(params):
        return np.array(photometry.chromaticity(_xyz(base * np.exp(params[0] * t + params[1] * t * t))))

    params = np.zeros(2)
    for _ in range(20):
        residual = chromaticity_of(params) - target
        if float(np.max(np.abs(residual))) < 1e-7:
            break
        step = 1e-4
        jacobian = np.column_stack([(chromaticity_of(params + step * unit) - chromaticity_of(params)) / step for unit in np.eye(2)])
        params = params - np.clip(np.linalg.solve(jacobian, residual), -1.0, 1.0)
    return base * np.exp(params[0] * t + params[1] * t * t)


def _smoothstep(tau):
    u = min(max((tau - THIN_BELOW_TAU) / (THICK_ABOVE_TAU - THIN_BELOW_TAU), 0.0), 1.0)
    return u * u * (3.0 - 2.0 * u)


def _single_scatter(tau_rayleigh, tau_extinction, aerosol_optics, mu0):
    """The zenith radiance per unit star flux of a homogeneous layer, with the phase mixed by the scattering depths (the sun is
    `acos(mu0)` from the zenith, so is the scattering angle)."""
    mu = 1.0
    cos_scatter = mu0
    scattering_phase = tau_rayleigh * 0.75 * (1.0 + cos_scatter ** 2)
    for tau, omega, g in aerosol_optics:
        scattering_phase = scattering_phase + tau * omega * (1.0 - g * g) / (1.0 + g * g - 2.0 * g * cos_scatter) ** 1.5
    albedo_phase = scattering_phase / tau_extinction
    if abs(mu0 - mu) < 1e-6:
        f = tau_extinction / (mu * mu) * np.exp(-tau_extinction / mu)
    else:
        f = (np.exp(-tau_extinction / mu0) - np.exp(-tau_extinction / mu)) / (mu0 - mu)
    return albedo_phase * mu0 * f / (4.0 * math.pi)


def _case(tau_rayleigh, tau_absorption, tau_ozone, above, aerosols, ground_albedo, mu0):
    """One sky (clear or cloudy): its diffuse flux, hue correction and direct beam, as fractions of the star's horizontal flux."""
    optics = [aerosol.optics(WAVELENGTHS_NM) for aerosol in aerosols]
    tau_aerosol = sum((tau for tau, _omega, _g in optics), np.zeros_like(WAVELENGTHS_NM))
    tau_extinction = tau_rayleigh + tau_absorption + tau_aerosol

    def gas_slab(share):
        tau = share * (tau_rayleigh + tau_absorption)
        omega = np.where(tau > 0.0, share * tau_rayleigh / np.where(tau > 0.0, tau, 1.0), 1.0)
        return tau, omega, np.zeros_like(tau)

    lower_slab = gas_slab(1.0 - above)
    for tau, omega, g in optics:
        lower_slab = radiation.mix_layers(lower_slab, (tau, omega, g))
    upper = two_stream.scattering_layer(*gas_slab(above), mu0)
    lower = two_stream.scattering_layer(*lower_slab, mu0)
    ozone = two_stream.absorbing_layer(tau_ozone, mu0)

    plain = two_stream.stack_layers([upper, lower])
    arriving_plain, _out = two_stream.ground_flux(plain, 0.0)
    diffuse_plain = np.maximum(arriving_plain - plain.beam_transmission, 0.0)
    full = two_stream.stack_layers([upper, ozone, lower])
    arriving, _out = two_stream.ground_flux(full, ground_albedo)
    diffuse = np.maximum(arriving - full.beam_transmission, 0.0)

    weight = _smoothstep(float(np.interp(_REFERENCE_NM, WAVELENGTHS_NM, tau_extinction)))
    correction = np.ones_like(diffuse)
    if weight < 1.0:
        radiance = _single_scatter(tau_rayleigh, tau_extinction, optics, mu0)
        ratio = np.where(diffuse_plain > 0.0, radiance / np.where(diffuse_plain > 0.0, diffuse_plain, 1.0), 0.0)
        reference = float(np.interp(_REFERENCE_NM, WAVELENGTHS_NM, ratio))
        usable = (ratio > 0.0) & (reference > 0.0)
        correction = np.where(usable, (np.where(usable, ratio, 1.0) / (reference if reference > 0.0 else 1.0)) ** (1.0 - weight), 1.0)
    return diffuse, diffuse * correction, full.beam_transmission, weight


def _colour_of(column, star, ozone_du, ground_albedo, mu0, haze, dust, cloud, cloud_cover):
    tau_rayleigh, tau_absorption = radiation.gas_layers(column, WAVELENGTHS_NM, EDGES_NM)
    tau_ozone = np.zeros_like(WAVELENGTHS_NM)
    above = 0.0
    if ozone_du > 0.0:
        tau_ozone = ozone_du * DOBSON_CM2 * chappuis_cross_section_cm2()
        above = min(max(radiation.ozone_peak_pressure_pa(column, ozone_du) / column.pressure_pa, 0.0), 1.0)
    always = [layer for layer in (haze, dust) if layer is not None and layer.tau_550 > 0.0]
    clear = _case(tau_rayleigh, tau_absorption, tau_ozone, above, always, ground_albedo, mu0)
    cases = [(1.0, clear)]
    if cloud is not None and cloud.tau_550 > 0.0 and cloud_cover > 0.0:
        cover = min(max(cloud_cover, 0.0), 1.0)
        cases = [(1.0 - cover, clear), (cover, _case(tau_rayleigh, tau_absorption, tau_ozone, above, always + [cloud], ground_albedo, mu0))]
    diffuse = sum((share * case[0] for share, case in cases), np.zeros_like(WAVELENGTHS_NM))
    hue = sum((share * case[1] for share, case in cases), np.zeros_like(WAVELENGTHS_NM))
    beam = sum((share * case[2] for share, case in cases), np.zeros_like(WAVELENGTHS_NM))
    weights = [case[3] for _share, case in cases if _share > 0.0]
    regime = "thin" if max(weights) == 0.0 else "thick" if min(weights) == 1.0 else "blend"
    tau_550 = float(np.interp(_REFERENCE_NM, WAVELENGTHS_NM, tau_rayleigh + tau_absorption + sum(
        (layer.optics(WAVELENGTHS_NM)[0] for layer in always), np.zeros_like(WAVELENGTHS_NM))))
    tau_rayleigh_550 = float(np.interp(_REFERENCE_NM, WAVELENGTHS_NM, tau_rayleigh))
    sky = photometry.colour_from_xyz(_xyz(star * hue))
    star_luminance = _luminance(star)
    sun_luminance = _luminance(star * beam)
    sun = None
    if sun_luminance > 0.0 and np.all(star * beam >= 0.0):
        try:
            sun = photometry.colour_from_xyz(_xyz(star * beam))
        except ValueError:
            sun = None
    return sky, _luminance(star * diffuse) / star_luminance, sun, sun_luminance / star_luminance, tau_550, tau_rayleigh_550, regime


def atmosphere_colour(column, *, temperature_k, ozone_du=0.0, ground_albedo=0.0, mu0=DEFAULT_MU0,
                      haze=None, dust=None, cloud=None, cloud_cover=1.0, spread_x=0.0, spread_y=0.0):
    """
    The atmosphere colour of a world with a ground.

    Args:
        column (AirColumn): The air above the place.
        temperature_k (float): The star's effective temperature, K.
        ozone_du (float, optional): The ozone column, Dobson units (GEN.413).
        ground_albedo (float or ndarray, optional): The ground's reflectance, one number or one per wavelength of `WAVELENGTHS_NM`
            (GEN.441, GEN.462); 0 is a black ground.
        mu0 (float, optional): The cosine of the sun's zenith angle, 0.2 to 1 (the view is at the zenith).
        haze (Aerosol, optional): Aerosol mixed with the air (GEN.414).
        dust (Aerosol, optional): Dust mixed with the lower air (GEN.414).
        cloud (Aerosol, optional): Cloud near the ground (GEN.414), covering `cloud_cover` of the sky.
        cloud_cover (float, optional): The share of the sky the cloud covers; the rest is clear.
        spread_x, spread_y (float, optional): A cool star's seeded chromaticity draws, -1 to 1 (`photometry.star_colour`).

    Returns:
        AtmosphereColour

    Raises:
        ValueError: `mu0` is outside 0.2 to 1 or the temperature is not positive.
    """
    if not 0.2 <= mu0 <= 1.0:
        raise ValueError("the two-stream model is made for the star between 0.2 and 1 in cosine of the zenith angle")
    star = star_spectrum(temperature_k, spread_x, spread_y)
    albedo = np.broadcast_to(np.asarray(ground_albedo, dtype=float), WAVELENGTHS_NM.shape)
    sky, ratio, sun, sun_ratio, tau_550, tau_rayleigh_550, regime = _colour_of(
        column, star, ozone_du, albedo, mu0, haze, dust, cloud, cloud_cover)
    shift = 0.0
    if any(layer is not None and layer.tau_550 > 0.0 for layer in (haze, dust, cloud)):
        bare = _colour_of(column, star, ozone_du, albedo, mu0, None, None, None, 1.0)[0]
        shift = detmath.hypot(sky.x - bare.x, sky.y - bare.y)
    return AtmosphereColour(sky, colour_word(sky.rgb8), ratio, sun, sun_ratio, tau_550, tau_rayleigh_550, regime, shift,
                            shift >= AEROSOL_DOMINATES_SHIFT)


STAND_IN_GROUND_ALBEDO = 0.15
"""float: The ground's reflectance for a body whose ground is not stored (GEN.441 and GEN.462 replace it); the radiation model's default."""


def body_atmosphere_colour(body):
    """
    The `AtmosphereColour` of a planet or moon from what it stores: its air (`air_column.column_of`), its star's temperature and its
    ozone column, over a ground of `STAND_IN_GROUND_ALBEDO` and with no aerosol (haze, dust and cloud are not stored yet: DB.68).

    Returns `None` for an airless body, a body with no star, and a gas or ice giant (a giant needs its cloud deck, which is not stored yet,
    to have a colour worth showing).
    """
    if getattr(body, "body_type", "t") == "g":
        return None
    temperature = getattr(getattr(body, "star", None), "temperature", None)
    if not temperature or temperature <= 0.0:
        return None
    column = air_column.column_of(body)
    if column is None:
        return None
    return atmosphere_colour(column, temperature_k=temperature, ozone_du=getattr(body, "ozone_column_du", None) or 0.0,
                             ground_albedo=STAND_IN_GROUND_ALBEDO)
