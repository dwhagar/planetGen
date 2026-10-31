# planetgen/physics/latitude_bands.py

"""
The climate of a body in seven bands (GEN.1108; Boss, 4:52 and 4:53 PM EDT 2026-10-30)
=======================================================================================

Every body says how much day and night it has, and how hot and cold, at various places from the equator. The places are seven bands of latitude of equal width
(`BAND_WIDTH_DEG`, 90 over 3.5): the equator band (12.9 degrees each side of the equator), and three bands on each side of it, the last of them taking its pole. This
module is a pure function of the stored columns (no draws, no epoch); `regional_temperature` stores them (GEN.1108).

For a body that turns on its own, each band holds, from a circular orbit round a star of angular radius `alpha`:

* **The hours of day and night** in a solar day: the share of one turn the upper limb of the star is up at the band's latitudes (`daylight.day_share_at`), averaged over
  the year and over the band's surface, times the solar day. A band inside the polar circle averages the long summer day and the long winter night over the year.
* **The longest polar day and the longest polar night** of the year, in days, at the band's poleward edge: the star is up all day while its declination `d` is at least
  `90 - |latitude| - alpha`, and `sin d = sin(tilt) sin(year angle)`, so the share of the year is `(pi - 2 asin(cos(|latitude| + alpha) / sin(tilt))) / 2 pi` (0 when
  the ratio is 1 or more); the night is the same with `cos(|latitude| - alpha)`. Earth's pole: 184 days of day and 179 of night (refraction would add a few days of day).
* **The temperature** as `regional_temperature` works it out (the annual-mean insolation of the band, the fourth-root law, the air carrying heat poleward), the **day and the night
  temperature** (the mean plus and minus half the day and night swing, which falls with the cosine of the latitude), and the **state** word: `day and night`, `polar day and
  night` (the band has a day or a night of a day or more at its pole side), `non-rotating` (the solar day is the year).

A world that turns once an orbit has no latitude day: its seven bands are rings round the point under the star, the same width, from the day side to the night side. Their
states are `permanent day side`, `terminator band` (the ring that holds the line between them) and `permanent night side`; the temperatures follow the cosine of the angle
from the substellar point, shared by the air with `regional_temperature.LOCKED_TRANSPORT_X0`.

Limits. The orbit is circular, so a north band and its south twin are equal (the bands stay separate for the eccentric orbit the model will gain). The star is a disc without
refraction. The day and night swing is the damped-oscillator estimate of `regional_temperature` for the equator times the cosine of the band's latitude.
"""

import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import constants, daylight, regional_temperature as rt
from planetgen.physics import unit_registry as units

_yr_to_h = units.converter("yr", "h")

BAND_COUNT = 7
BAND_WIDTH_DEG = 90.0 / 3.5
"""float: The width of each band in latitude, degrees: the equator band is this wide centred on the equator and each of the three on a side is this wide. [Stand-in LB1, tier T9, tweak L0.]"""

NAMES = ("north_polar", "north_middle", "north_tropic", "equator", "south_tropic", "south_middle", "south_polar")
"""tuple: The bands from the north pole to the south pole."""

RING_NAMES = ("day_centre", "day_middle", "day_edge", "terminator", "night_edge", "night_middle", "night_centre")
"""tuple: The rings of a world that turns once an orbit, from the point under the star to the point away from it."""

CYCLE = "day and night"
POLAR = "polar day and night"
NON_ROTATING = "non-rotating"
DAY_SIDE, TERMINATOR, NIGHT_SIDE = "permanent day side", "terminator band", "permanent night side"

RING_STATES = (DAY_SIDE, TERMINATOR, NIGHT_SIDE)
"""tuple: The states of a world's rings, which tell a ring from a latitude band."""

_STEPS = 12
"""int: The latitudes sampled across a band."""

_SEASONS = 72
"""int: The points of the year the day share is averaged over."""


@dataclass(frozen=True)
class Band:
    """One band: its place (degrees of latitude from the north edge to the south edge, or for a locked world's ring the angle from the point under the star), its share
    of the surface and its climate. A field is `None` where the band has none (a locked world's rings have no hours)."""

    name: str
    from_deg: float
    to_deg: float
    area_share: float
    mean_k: float | None
    day_k: float | None
    night_k: float | None
    swing_k: float | None
    day_hours: float | None
    night_hours: float | None
    longest_day_days: float | None
    longest_night_days: float | None
    state: str | None


def edges():
    """The latitude edges, north to south: `[90, 64.3, 38.6, 12.9, -12.9, -38.6, -64.3, -90]`."""
    half = BAND_WIDTH_DEG / 2.0
    return [90.0, 90.0 - BAND_WIDTH_DEG, 90.0 - 2.0 * BAND_WIDTH_DEG, half, -half, -(half + BAND_WIDTH_DEG), -(half + 2.0 * BAND_WIDTH_DEG), -90.0]


def _bounds():
    e = edges()
    return [(e[i], e[i + 1]) for i in range(BAND_COUNT)]


def _area_share(high, low):
    return (math.sin(math.radians(high)) - math.sin(math.radians(low))) / 2.0


def _latitudes(high, low):
    return low + (np.arange(_STEPS) + 0.5) / _STEPS * (high - low)


def _weights(latitudes):
    return np.cos(np.radians(latitudes))


def _mean_cosine(high, low):
    """The area-weighted mean cosine of the latitude over a band: the day and night swing of a band over the equator band's goes as this (see `regional_temperature`)."""
    weight = _weights(_latitudes(high, low))
    return float(np.sum(weight * weight) / np.sum(weight))


def polar_day_share(latitude_deg, tilt_deg, alpha):
    """The share of a circular year the star stays up all day at a latitude (the north's summer is the south's winter): `0` when the tilt is too small to reach it."""
    sin_tilt = math.sin(math.radians(daylight.folded_tilt_deg(tilt_deg)))
    if sin_tilt <= 0.0:
        return 0.0
    ratio = math.cos(math.radians(abs(latitude_deg)) + alpha) / sin_tilt
    return 0.0 if ratio >= 1.0 else (math.pi - 2.0 * math.asin(max(ratio, -1.0))) / (2.0 * math.pi)


def polar_night_share(latitude_deg, tilt_deg, alpha):
    """The share of a circular year the star stays down all day at a latitude."""
    sin_tilt = math.sin(math.radians(daylight.folded_tilt_deg(tilt_deg)))
    if sin_tilt <= 0.0:
        return 0.0
    ratio = math.cos(math.radians(abs(latitude_deg)) - alpha) / sin_tilt
    return 0.0 if ratio >= 1.0 else (math.pi - 2.0 * math.asin(min(ratio, 1.0))) / (2.0 * math.pi)


def _mean_day_share(latitudes, weights, tilt_deg, alpha):
    """The share of the time the star is up at a band, averaged over its latitudes (by area) and over the year."""
    sin_tilt = math.sin(math.radians(daylight.folded_tilt_deg(tilt_deg)))
    total = 0.0
    for k in range(_SEASONS):
        declination = math.degrees(math.asin(sin_tilt * math.sin(2.0 * math.pi * (k + 0.5) / _SEASONS)))
        total += sum(w * daylight.day_share_at(lat, declination, alpha) for lat, w in zip(latitudes, weights))
    return total / _SEASONS / float(np.sum(weights))


def _band_insolation(tilt_deg):
    """Each band's annual-mean insolation over the planet's mean, area-weighted over its latitudes (a sphere averages 1)."""
    tilt = math.radians(min(max(tilt_deg or 0.0, 0.0), 180.0))
    out = []
    for high, low in _bounds():
        latitudes = np.radians(_latitudes(high, low))
        weight = np.cos(latitudes)
        out.append(float(4.0 * np.sum(weight * rt._annual_insolation(np.abs(latitudes), tilt)) / np.sum(weight)))
    return out


def _ring_insolation():
    """Each ring's mean insolation over the planet's mean for a world locked to its star: a ring `i` of the angle `theta` from the point under the star gets `max(cos theta, 0)`."""
    out = []
    for i in range(BAND_COUNT):
        low, high = i * BAND_WIDTH_DEG, (i + 1) * BAND_WIDTH_DEG
        theta = np.radians(low + (np.arange(_STEPS) + 0.5) / _STEPS * (high - low))
        weight = np.sin(theta)
        out.append(float(4.0 * np.sum(weight * np.maximum(np.cos(theta), 0.0)) / np.sum(weight)))
    return out


def _ring_state(i):
    centre = (i + 0.5) * BAND_WIDTH_DEG
    if centre < 90.0 - BAND_WIDTH_DEG / 2.0:
        return DAY_SIDE
    if centre > 90.0 + BAND_WIDTH_DEG / 2.0:
        return NIGHT_SIDE
    return TERMINATOR


def _temperatures(radiative, areas, mean, share, floor=False):
    """`radiative` ratios (fourth root of the insolation) turned into band temperatures: scaled so the area-weighted mean of the bands is the world's mean, then pulled toward it by the air."""
    scale = sum(a * r for a, r in zip(areas, radiative)) / sum(areas)
    out = []
    for r in radiative:
        temperature = mean + (r / scale * mean - mean) * (1.0 - share)
        out.append(max(temperature, rt.NIGHT_FLOOR_K) if floor else temperature)
    return out


def alpha_deg_of(values):
    """The angular radius of the star as the body sees it, degrees, from the stored mean share of the time the star is up (`daylight.mean_day_share`: `(1 + sin alpha) / 2`); 0 for none."""
    share = values.get("day_share_mean")
    if share is None:
        return 0.0
    return math.degrees(math.asin(min(max(2.0 * share - 1.0, 0.0), 1.0)))


def compute(values, alpha_deg=0.0):
    """
    The seven `Band`s from north to south (rings from the day side to the night side for a world that turns once an orbit) for a rocky body from its column-named `values`:
    `surface_temperature_k`, `atmospheric_pressure_pa`, `gravity_g`, `rotation_period_hours`, `axial_tilt_deg` and `period_years`, and the angular radius of the star as it sees
    it, `alpha_deg` (0 for a point star). `[]` for a body without a temperature or a rotation period.
    """
    mean = values.get("surface_temperature_k")
    rotation = values.get("rotation_period_hours")
    if not mean or mean <= 0.0 or not rotation or rotation <= 0.0:
        return []
    pressure = values.get("atmospheric_pressure_pa") or 0.0
    gravity = (values.get("gravity_g") or 0.0) * constants.EARTH_GRAVITY
    period_years = values.get("period_years")
    tilt = values.get("axial_tilt_deg") or 0.0
    alpha = math.radians(alpha_deg or 0.0)
    column = rt.column_in_earths(pressure, gravity)
    locked = rt.is_locked(rotation, period_years)
    share = rt.transport_share(column, locked)
    areas = [_area_share(high, low) for high, low in _bounds()]
    bands = []
    if locked:
        areas = [(math.cos(math.radians(i * BAND_WIDTH_DEG)) - math.cos(math.radians((i + 1) * BAND_WIDTH_DEG))) / 2.0 for i in range(BAND_COUNT)]
        insolation = _ring_insolation()
        temperatures = _temperatures([max(i, 0.0) ** 0.25 for i in insolation], areas, mean, share, floor=True)
        for i, (high, low) in enumerate(_bounds()):
            bands.append(Band(RING_NAMES[i], i * BAND_WIDTH_DEG, (i + 1) * BAND_WIDTH_DEG, areas[i], temperatures[i], None, None, 0.0, None, None, None, None, _ring_state(i)))
        return bands
    orbital_hours = _yr_to_h(period_years) if period_years else None
    solar_day = daylight.solar_day_hours(rotation, orbital_hours, tilt <= 90.0) if orbital_hours else None
    insolation = _band_insolation(tilt)
    temperatures = _temperatures([i ** 0.25 for i in insolation], areas, mean, share)
    equator_swing = rt.swing_k(mean, pressure, gravity, rotation) * (1.0 - share)
    cosines = [_mean_cosine(high, low) for high, low in _bounds()]
    equator_cosine = cosines[BAND_COUNT // 2]
    for i, (high, low) in enumerate(_bounds()):
        swing = equator_swing * cosines[i] / equator_cosine
        latitudes = _latitudes(high, low)
        edge = max(abs(high), abs(low))
        day_days = polar_day_share(edge, tilt, alpha) * (orbital_hours or 0.0) / 24.0
        night_days = polar_night_share(edge, tilt, alpha) * (orbital_hours or 0.0) / 24.0
        day_hours = night_hours = None
        if solar_day:
            day_share = _mean_day_share(latitudes, _weights(latitudes), tilt, alpha)
            day_hours, night_hours = day_share * solar_day, (1.0 - day_share) * solar_day
        state = NON_ROTATING if solar_day and orbital_hours and abs(solar_day - orbital_hours) <= 0.01 * orbital_hours else (
            POLAR if day_days >= 1.0 or night_days >= 1.0 else CYCLE)
        bands.append(Band(NAMES[i], high, low, areas[i], temperatures[i], temperatures[i] + swing / 2.0, temperatures[i] - swing / 2.0, swing,
                          day_hours, night_hours, day_days if orbital_hours else None, night_days if orbital_hours else None, state))
    return bands
