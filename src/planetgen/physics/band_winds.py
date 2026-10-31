# planetgen/physics/band_winds.py

"""
Average surface winds by day and by night in each latitude band (GEN.1110)
==========================================================================

Boss asked (4:53 PM EDT 2026-10-30) whether the wind speed by day and by night across the latitude bands can be estimated. As an order of magnitude it can, from the
temperatures the bands already have (`regional_temperature`, in seven bands once GEN.1108 lands), the air's pressure and molar mass, and the planet's radius and spin.
This module is that estimate as a pure function; nothing calls it yet (GEN.1110's wiring and storage follow GEN.1108).

**The heat-engine speed.** A temperature difference `dT` at a mean `T` can drive air no faster than the kinetic energy a heat engine of Carnot efficiency `dT / T` gets from
the work `R_s dT` the pressure difference does: `u0 = sqrt((dT / T) R_s dT)`, `R_s` the air's specific gas constant. It is the scale of a sea breeze, of the overturning
of a slowly turning world and of the flow from the day side to the night side of a locked one. It needs no rotation and it is an upper bound on the wind aloft.

**The thermal wind.** Where the world turns fast enough that the Coriolis force balances the pressure gradient, the equator-to-pole difference `dT_ep` sets a shear
`u_g = (R_s / f) (dT_ep / L) ln(P_s / P_top)` with `f = 2 Omega sin(latitude)`, `L` the distance from the equator to the pole `pi a / 2` and `ln(P_s / P_top)`
`THERMAL_WIND_LOG_PRESSURE` (the surface to the tropopause). The band takes the smooth minimum of the two, `u0 u_g / sqrt(u0^2 + u_g^2)`: the thermal wind where `u_g < u0`
(rapid rotation, a small thermal Rossby number `R_s dT_ep / (Omega a)^2`), the heat-engine speed where it is not (the equator, where `f` goes to zero, and the slow rotators).
Each band is named for the regime it is in: **overturning** for a world with a thermal Rossby number of 1 or more (the Hadley cell reaches the poles, Held and Hou 1980);
otherwise **trades** inside the Hadley cell's edge `sqrt(5/3 Ro_T)` radians, **westerlies** out to 60 degrees and **polar** beyond.

**The surface.** Winds aloft reach the ground at a share of their speed (`SURFACE_COUPLING`, the boundary layer's drag), and at night, when the layer at the ground is
stable and cut off from the flow above, at `NIGHT_COUPLING` of that. The day-to-night flow of a band (`u0` of the band's own day and night temperatures, cut at high
latitude by the inertial limit `1 / sqrt(1 + (2 sin(latitude))^2)`: a diurnal flow cannot keep a direction where the Coriolis period is shorter than the day) is added in
quadrature to the zonal flow. A locked world has no day and night to cycle through: the day side and the night side each take the flow from one to the other.

**Thin air.** Below about 1 kPa the air has too little mass to carry a wind to the ground (the stress it can exert falls with its density): every speed is multiplied by
`1 / (1 + (P0 / P)^2)`, a half at `PRESSURE_FADE_PA` (Boss's "about 1 kPa" in the item) and a fade to nothing below.

**Test targets and limits.** The Earth's surface winds are 5 to 8 m/s in the mid-latitudes and 4 to 6 in the trades, Mars's 3 to 10, Venus's under 1 m/s and Titan's under 1;
the module is checked to the same order, not to the values (the two coupling stand-ins are tuned to the Earth, tweak level L3). It is a mean speed, not a gust, a direction
or a storm; a world without a surface (a giant) has no surface wind and its speed here is the flow at the level of the temperatures given. Mountains, coasts and the seasons
are not in it. Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.98.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants
from planetgen.util import detmath
from planetgen.util.checks import finite_domain

PRESSURE_FADE_PA = 1000.0
"""float: The pressure at which the surface wind is half what the air aloft allows (the item's "about 1 kPa"). Source tier T9, tweak level L2."""
PRESSURE_FADE_EXPONENT = 2.0
"""float: How sharply the fade falls below `PRESSURE_FADE_PA`: `1 / (1 + (P0 / P)^n)`. Source tier T9, tweak level L3."""
THERMAL_WIND_LOG_PRESSURE = 1.6
"""float: `ln(P_surface / P_tropopause)` of the thermal wind, `ln 5` (a tropopause at a fifth of the surface pressure; the Earth's is 0.2). Source tier T6, tweak level L1."""
SURFACE_COUPLING = 0.3
"""float: The share of the wind aloft that reaches the surface by day. Source tier T9, tweak level L3 (tuned to the Earth's mid-latitude surface winds, 5 to 8 m/s)."""
NIGHT_COUPLING = 0.5
"""float: The share of the daytime surface coupling that is left at night, when the layer at the ground is stable. Source tier T9, tweak level L3 (the Earth's night winds are about half the day's)."""
TROPICS_LIMIT_DEG = 60.0
"""float: The latitude where the westerlies give way to the polar winds, degrees. Source tier T6, tweak level L0."""

OVERTURNING = "overturning"
TRADES = "trades"
WESTERLIES = "westerlies"
POLAR = "polar"
LOCKED = "day to night"
REGIMES = (OVERTURNING, TRADES, WESTERLIES, POLAR, LOCKED)


class BandWindError(ValueError):
    """An input the model cannot use."""


@dataclass(frozen=True)
class Band:
    """A latitude band's centre and its mean temperatures by day and by night (K)."""

    latitude_deg: float
    day_k: float
    night_k: float


@dataclass(frozen=True)
class BandWind:
    """
    The mean surface wind of one band.

    Attributes:
        latitude_deg (float): The band's.
        day_m_s, night_m_s (float): The mean speed by day and by night, m/s (registry unit `m_s`).
        regime (str): One of `REGIMES`.
        zonal_aloft_m_s (float): The zonal flow aloft that the surface winds are a share of.
        day_night_m_s (float): The band's day-to-night flow aloft, before the surface share and the fade.
    """

    latitude_deg: float
    day_m_s: float
    night_m_s: float
    regime: str
    zonal_aloft_m_s: float
    day_night_m_s: float


def specific_gas_constant(mean_molar_mass_kg_mol):
    """`R_s = R / M` of the air, J/(kg K)."""
    if not mean_molar_mass_kg_mol > 0.0:
        raise BandWindError("the air's molar mass must be positive")
    return constants.R / mean_molar_mass_kg_mol


def pressure_fade(pressure_pa):
    """The share of the wind the air's mass lets reach the ground: a half at `PRESSURE_FADE_PA`, 1 in thick air, 0 with no air."""
    if pressure_pa < 0.0:
        raise BandWindError("the pressure must not be negative")
    if pressure_pa == 0.0:
        return 0.0
    return 1.0 / (1.0 + (PRESSURE_FADE_PA / pressure_pa) ** PRESSURE_FADE_EXPONENT)


def heat_engine_speed_m_s(delta_k, mean_k, gas_constant_j_kg_k):
    """`u0 = sqrt((dT / T) R_s dT)` for a temperature difference `delta_k` at a mean `mean_k`, m/s."""
    if not mean_k > 0.0 or not gas_constant_j_kg_k > 0.0:
        raise BandWindError("the mean temperature and the gas constant must be positive")
    return math.sqrt(abs(delta_k) / mean_k * gas_constant_j_kg_k * abs(delta_k))


def thermal_rossby(delta_ep_k, gas_constant_j_kg_k, radius_m, rotation_rate_rad_s):
    """`R_s dT_ep / (Omega a)^2`: infinite for a world that does not turn."""
    omega_a = abs(rotation_rate_rad_s) * radius_m
    return math.inf if omega_a == 0.0 else gas_constant_j_kg_k * abs(delta_ep_k) / omega_a ** 2


def thermal_wind_m_s(latitude_deg, delta_ep_k, gas_constant_j_kg_k, radius_m, rotation_rate_rad_s):
    """The thermal-wind shear speed `(R_s / f)(dT_ep / L) ln(P_s / P_top)` at a latitude; infinite where `f` is zero."""
    f = 2.0 * abs(rotation_rate_rad_s) * abs(math.sin(math.radians(latitude_deg)))
    if f == 0.0:
        return math.inf
    return gas_constant_j_kg_k / f * abs(delta_ep_k) / (0.5 * math.pi * radius_m) * THERMAL_WIND_LOG_PRESSURE


def _smooth_minimum(a, b):
    if math.isinf(a):
        return b
    if math.isinf(b):
        return a
    return a * b / detmath.hypot(a, b) if a > 0.0 and b > 0.0 else 0.0


def hadley_edge_deg(thermal_rossby_number):
    """The latitude, degrees, where the Hadley cell ends: `sqrt(5/3 Ro_T)` radians (Held and Hou 1980), the pole for a thermal Rossby number of 1 or more."""
    if thermal_rossby_number >= 1.0:
        return 90.0
    return min(90.0, math.degrees(math.sqrt(5.0 / 3.0 * thermal_rossby_number)))


@finite_domain(allow_inf=("rotation_rate_rad_s",))
def band_winds(bands, pressure_pa, mean_molar_mass_kg_mol, radius_m, rotation_rate_rad_s, locked=False):
    """
    The mean surface wind by day and by night in each band.

    Args:
        bands (sequence): `Band`s, at least two with different `|latitude|`; the equator-to-pole contrast is between the bands nearest and farthest from the equator.
        pressure_pa (float): The surface pressure.
        mean_molar_mass_kg_mol (float): The air's mean molar mass.
        radius_m (float): The planet's radius.
        rotation_rate_rad_s (float): The spin rate, rad/s (0 for a world that does not turn).
        locked (bool): Tidally locked: the day side and the night side take the flow from one to the other, and the bands' latitudes only name where.

    Returns:
        tuple: A `BandWind` for each band, in the order given.
    """
    bands = tuple(bands)
    if len(bands) < 2 or len({abs(b.latitude_deg) for b in bands}) < 2:
        raise BandWindError("winds need at least two bands at different latitudes")
    if any(not (b.day_k > 0.0 and b.night_k > 0.0) or not -90.0 <= b.latitude_deg <= 90.0 for b in bands):
        raise BandWindError("a band needs temperatures above zero and a latitude between -90 and 90")
    if not radius_m > 0.0:
        raise BandWindError("the radius must be positive")
    gas = specific_gas_constant(mean_molar_mass_kg_mol)
    fade = pressure_fade(pressure_pa)
    mean = {id(b): 0.5 * (b.day_k + b.night_k) for b in bands}
    equatorial = min(bands, key=lambda b: abs(b.latitude_deg))
    polar = max(bands, key=lambda b: abs(b.latitude_deg))
    delta_ep = abs(mean[id(equatorial)] - mean[id(polar)])
    mean_world = math.fsum(mean.values()) / len(bands)
    ro = thermal_rossby(delta_ep, gas, radius_m, rotation_rate_rad_s)
    edge = hadley_edge_deg(ro)
    u0_ep = heat_engine_speed_m_s(delta_ep, mean_world, gas)
    out = []
    for band in bands:
        latitude = abs(band.latitude_deg)
        t_mean = mean[id(band)]
        if locked:
            day_side = heat_engine_speed_m_s(band.day_k - band.night_k, t_mean, gas) if band.day_k != band.night_k else 0.0
            flow = max(day_side, heat_engine_speed_m_s(delta_ep, mean_world, gas))
            out.append(BandWind(band.latitude_deg, SURFACE_COUPLING * flow * fade, SURFACE_COUPLING * NIGHT_COUPLING * flow * fade, LOCKED, 0.0, flow))
            continue
        aloft = _smooth_minimum(u0_ep, thermal_wind_m_s(latitude, delta_ep, gas, radius_m, rotation_rate_rad_s))
        inertial = 1.0 / math.sqrt(1.0 + (2.0 * math.sin(math.radians(latitude))) ** 2)
        tide = heat_engine_speed_m_s(band.day_k - band.night_k, t_mean, gas) * inertial
        day = detmath.hypot(SURFACE_COUPLING * aloft, SURFACE_COUPLING * tide) * fade
        night = detmath.hypot(SURFACE_COUPLING * NIGHT_COUPLING * aloft, SURFACE_COUPLING * tide) * fade
        if ro >= 1.0:
            regime = OVERTURNING
        elif latitude < edge:
            regime = TRADES
        elif latitude < TROPICS_LIMIT_DEG:
            regime = WESTERLIES
        else:
            regime = POLAR
        out.append(BandWind(band.latitude_deg, day, night, regime, aloft, tide))
    return tuple(out)
