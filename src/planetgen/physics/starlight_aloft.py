# planetgen/physics/starlight_aloft.py

"""
The starlight an air takes up before it reaches its reference level (GEN.1717; docs/design/composition/infrared-opacity.md)
==========================================================================================================================

A body with no ground (a giant, or a planet under a deep hydrogen envelope, GEN.1334) has no surface to heat: the temperature it stores is the one
at its 1 bar reference level (`surface_air.REFERENCE_PRESSURE_BAR`). `radiative_convective` deposits the sunlight the caller says goes aloft and
puts all of the rest into that level, as if it were ground. For an air thin in absorbers that is right (Jupiter's 1 bar air is hydrogen and helium,
the water, methane and ammonia of its sunlight bands are far below or far above it). For an envelope with tens of times the solar metals it is not:
the near-infrared bands of its water, methane and ammonia take the starlight high up, where the air radiates it back to space, and the 1 bar level
gets only what reaches it (GEN.1717: the table stored 2,600 K at 0.5 AU for a 5 Earth-mass core under 4 percent hydrogen, the model with the sunlight
at the 1 bar level 1,200 K).

This module finds how much of the starlight the gases take above the reference level and where, from the same correlated-k tables and collision-induced
pairs the model reads for the infrared (`line_opacity`), on the bands below 12,000 cm-1 (0.83 micron and longer), weighted by the star's own
blackbody. What it returns is `(share, pressure_pa)` in the form `radiative_convective.solve` takes (`aloft_share`, `aloft_pressure_pa`).

Limits (stand-ins, registry rows SA1 to SA3, docs/design/research/gap-filling-plan.md section 6.293):

* Starlight shortward of 12,000 cm-1 (the visible and the ultraviolet, about 60 percent of a Sun's) is taken to reach the reference level (the tables
  stop there; Rayleigh scattering and the visible bands of methane and ammonia are left to the albedo and to the haze).
* The air is isothermal at the temperature the caller gives (the step's start), because the line strengths change little over the few tens of kelvin
  the fixed point moves, and the profile itself is not known before the solve.
* The beam crosses the air at a fixed slant, `SLANT` (the global mean of the sunlight a sphere takes).

Use:

    deposit = starlight_aloft.of_air(partial_pressures_pa, gravity_m_s2, star_temperature_k, temperature_k)
    deposit.share, deposit.pressure_pa
"""

import functools
import math
from collections import namedtuple

import numpy as np

from planetgen.physics import constants, line_opacity, radiative_convective

SLANT = 2.0
"""float: STAND-IN (SA1, tier T7, L0): the secant of the angle the beam crosses the air at, 2 (a cosine of one half, the one a global mean over day and
night takes for the flux a sphere absorbs). The shortwave in a 1D model is conventionally taken with such a fixed angle."""
COVERED_BELOW_CM1 = float(line_opacity.EDGES_CM1[-1])
"""float: The highest wavenumber the tables cover, cm-1 (12,000: 0.83 micron). The starlight above it is not absorbed here (SA2)."""
MAX_SHARE = 0.9
"""float: STAND-IN (SA3, tier T7, L0): The most the gases take aloft, as a share of the absorbed starlight: the model needs a share under 1, and the
table that stands for it was fitted up to 0.9 (`greenhouse_table`)."""
MIN_SHARE = 1.0e-3
"""float: A share under this is no sunlight aloft at all (the bands the air has no gas for, or too little of)."""
SUB_STEPS = 12
"""int: Points of the Planck function across one band, for the star's share in it."""

Deposit = namedtuple("Deposit", "share pressure_pa")
"""The share of the absorbed starlight the gases take above the reference level, and the pressure (Pa) it is taken at (the mean of ln P, weighted by
the heating)."""
NONE = Deposit(0.0, 0.0)


@functools.lru_cache(maxsize=256)
def band_shares(star_temperature_k):
    """The share of a blackbody's flux (`star_temperature_k`) that falls in each of the tables' bands, over the whole spectrum (the bands cover its part below
    `COVERED_BELOW_CM1`)."""
    edges = line_opacity.EDGES_CM1
    c2 = line_opacity.SECOND_RADIATION_CM_K
    shares = np.empty(len(edges) - 1)
    for index in range(len(shares)):
        nu = np.linspace(edges[index], edges[index + 1], SUB_STEPS + 1)
        nu = 0.5 * (nu[1:] + nu[:-1])
        x = c2 * nu / star_temperature_k
        shares[index] = np.sum(nu ** 3 / np.expm1(np.minimum(x, 700.0))) * (edges[index + 1] - edges[index]) / SUB_STEPS
    total = (math.pi ** 4 / 15.0) * (star_temperature_k / c2) ** 4          # the integral of nu^3 / (exp(c2 nu / T) - 1) over all wavenumbers
    return shares / total


def of_air(partial_pressures_pa, gravity_m_s2, star_temperature_k, temperature_k):
    """
    The starlight the air takes above its bottom level: `(share, pressure_pa)`, `NONE` when the air takes none or the body has no star.

    Args:
        partial_pressures_pa (dict): Formula to partial pressure at the reference level, Pa.
        gravity_m_s2 (float): Surface gravity.
        star_temperature_k (float or None): The effective temperature of the star, K.
        temperature_k (float): The temperature the air is taken at (isothermal), K.
    """
    if not star_temperature_k or star_temperature_k <= 0.0 or not temperature_k or temperature_k <= 0.0:
        return NONE
    column = radiative_convective.Column(partial_pressures_pa, gravity_m_s2)
    if not column.prepared.absorbers:
        return NONE
    depth, _ = column.optical_depth(np.full(column.layers + 1, float(temperature_k)))        # (layer, band, g)
    weights = line_opacity.G_WEIGHTS
    # transmission at the top of each layer and at the bottom of the last, band-mean over the g points
    above = np.concatenate([np.zeros((1,) + depth.shape[1:]), np.cumsum(depth, axis=0)])
    transmitted = np.einsum("lbg,g->lb", np.exp(-SLANT * above), weights)             # (layer + 1, band)
    taken = (transmitted[:-1] - transmitted[1:]) * band_shares(float(star_temperature_k))[None, :]        # (layer, band)
    heating = taken.sum(axis=1)
    share = float(heating.sum())
    if share < MIN_SHARE:
        return NONE
    pressure = math.exp(float(np.sum(heating * np.log(column.layer_pressure)) / share))
    return Deposit(min(share, MAX_SHARE), pressure)


def combine(haze_share, haze_pressure_pa, gas):
    """
    The haze's share and the gases' as one: the haze takes its share first (it sits above), the gases take their own share of what is left, at their
    pressure, and the pressure of the two is the mean of ln P weighted by the share each takes.

    Args:
        haze_share (float): The share the haze takes (`aerosols.Haze.absorbed_share`).
        haze_pressure_pa (float or None): The pressure it is taken at.
        gas (Deposit): `of_air`.

    Returns:
        Deposit: `(share, pressure_pa)`, `NONE` if neither takes any.
    """
    gas_share = (1.0 - haze_share) * gas.share
    total = haze_share + gas_share
    if total <= 0.0:
        return NONE
    if haze_share <= 0.0 or not haze_pressure_pa:
        return Deposit(min(total, MAX_SHARE), gas.pressure_pa)
    if gas_share <= 0.0:
        return Deposit(min(total, MAX_SHARE), haze_pressure_pa)
    ln_p = (haze_share * math.log(haze_pressure_pa) + gas_share * math.log(gas.pressure_pa)) / total
    return Deposit(min(total, MAX_SHARE), math.exp(ln_p))
