# planetgen/physics/temperature_coverage.py

"""
What is read off the temperature of every place (GEN.1264)
===========================================================

`seasonal_swing.temperature_series` gives the temperature of each cell of a body at each step of its year and its turns. This module reads from it what the Temperature domain of
the habitability score and the pages need, as pure functions of that series (no draw, no stored column; the score, the columns and the migration are GEN.1265 and GEN.1266). The plan
is `reports/temperature-coverage-2026-11-04.md`, with the decisions of Boss of 7:32 PM EST 11-04.

**The bands.** The PHI-6 bands of `habitability.TEMPERATURE_BAND`: Blue 0 to 31 C, Green -20 to 45 C, Yellow -50 to 122 C, Red beyond. A sample is *livable* when it is in the Green
band or better (no suit needed). The wet bulb rule is applied to each sample: where the humidity is given, a sample whose wet bulb temperature (`habitability.wet_bulb_c`) is at or
above `habitability.WET_BULB_CRIT_C` is Yellow, not Blue or Green (cooling garments). The wet bulb rises with the dry bulb, so the rule is a dry-bulb ceiling found once. Shares are
over area (a cell of a row weighs `grid.row_weight / columns`: equal area) and time (each step weighs the same), so `blue + green + yellow + red = 1`, where `green` is Green and
not Blue, and `cold` (below -20 C) and `hot` (above 45 C) are the same samples counted by which end they fall off.

**The livable share.** `livable_share` L is the area-time mean of the livable indicator, the one number the score reads. The two pieces to show beside it are the *livable area*
(the share of the world livable at least `MOST_OF_THE_TIME` of the time) and the *livable time* a typical place of that area has (the area-weighted mean of its livable
time). `area_always`, `area_most` (90 percent) and `area_ever` are the shares livable all the time, at least `NEARLY_ALWAYS` of it and at some time.

**Extremes.** The coldest and hottest sample anywhere and where and when it is, and the same by latitude band (the seven of `latitude_bands`), by column of longitude, by quarter
of the orbit from the equinox, and, for a world that turns once an orbit, by sector of the sky of its star (`DAY_SIDE_DEG`, `NIGHT_SIDE_DEG` from the point under the star:
the day side, the terminator belt and the night side), with the livable share of each sector so that a world livable only in its terminator belt can be told.

**The best site** is the cell livable the largest share of the time; cells within `TIE_TOLERANCE` of the best share tie, and the tie goes to the more moderate one (the least
difference between its hottest and coldest sample). It carries its own coldest and hottest sample and the longest unbroken spell out of the livable band there: how long a person
would have to be inside (the series repeats, so a spell that runs over the end joins the start). A world with no livable cell has no best site; the cell with the most livable
time is still reported, with `livable_time` 0.

Limits. The series is a circular orbit with an evenly spread sea (`seasonal_swing`'s limits); an hour is a step of the series, not an instant, so a spell is a whole number of
steps; the seasons are quarters of the orbit and not the seasons of a hemisphere.
"""

import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import habitability, latitude_bands

KELVIN = 273.15
MOST_OF_THE_TIME = 0.5
"""float: A place counts in the livable area when it is livable at least this share of the time ("livable most of the time", Boss 7:32 PM EST 11-04). Source tier T9, tweak level L0 (a definition)."""
NEARLY_ALWAYS = 0.9
"""float: The share of the time a place is livable for `area_most` (the research note, "at least 90 percent of the time"). Source tier T9, tweak level L0 (a definition)."""
TIE_TOLERANCE = 1.0e-3
"""float: Cells whose livable time is within this of the best tie for the best site (a hundredth of a percent is a step or two of a year). Source tier T6, tweak level L0 (numerical)."""
DAY_SIDE_DEG = 60.0
"""float: A cell of a world that turns once an orbit is on the day side within this angle of the point under the star. Source tier T9, tweak level L0 (a definition: a quarter of the surface)."""
NIGHT_SIDE_DEG = 120.0
"""float: ... and on the night side beyond this angle, with the terminator belt (half the surface) between. Source tier T9, tweak level L0 (a definition)."""
QUARTERS = 4
"""int: The parts of the orbit the extremes by season are given for, from the equinox."""
SECTORS = ("day side", "terminator belt", "night side")
"""tuple: The sectors of the sky of the star for a world that turns once an orbit."""


@dataclass(frozen=True)
class Extreme:
    """One sample: the temperature, kelvin, the cell it is in (row, column, latitude and longitude in degrees) and the time (days from an equinox)."""

    temperature_k: float
    row: int
    column: int
    latitude_deg: float
    longitude_deg: float
    time_days: float


@dataclass(frozen=True)
class BestSite:
    """
    The cell livable the largest share of the time.

    Attributes:
        row, column, latitude_deg, longitude_deg: The cell.
        livable_time (float): The share of the time its temperature is livable.
        coldest_k, hottest_k (float): Its coldest and hottest sample, kelvin.
        longest_out_days (float): The longest unbroken spell out of the livable band, days (0 when it is always livable, the whole span when never).
        mean_k (float): Its mean temperature.
    """

    row: int
    column: int
    latitude_deg: float
    longitude_deg: float
    livable_time: float
    coldest_k: float
    hottest_k: float
    longest_out_days: float
    mean_k: float


@dataclass(frozen=True)
class Coverage:
    """
    What is read off a `TemperatureSeries` (the module text).

    Attributes:
        blue, green, yellow, red (float): The shares of area and time in each band (Green not counting Blue); they add to 1. The wet bulb rule moves samples to Yellow.
        cold, hot (float): The shares below -20 C and above 45 C.
        humid (float): The share the wet bulb rule moved to Yellow (0 when no humidity was given).
        livable_share (float): L, `blue + green`.
        livable_area, livable_time (float): The two pieces of L: the share of the world livable at least `MOST_OF_THE_TIME` of the time, and the mean livable time of it.
        area_always, area_most, area_ever (float): The shares of the world livable all the time, at least `NEARLY_ALWAYS` of it, and at some time.
        place_time (ndarray): `[row, column]`, the share of the time each cell is livable.
        coldest, hottest (Extreme): The coldest and hottest sample anywhere.
        by_band (tuple): `(coldest_k, hottest_k)` for each of `latitude_bands.NAMES`, north to south.
        by_column (tuple): `(coldest_k, hottest_k)` for each column of longitude.
        by_quarter (tuple): `(coldest_k, hottest_k)` for each quarter of the orbit from the equinox.
        sectors (dict or None): For a world that turns once an orbit, `{sector name: (livable share, coldest_k, hottest_k)}`; `None` for any other.
        best_site (BestSite): The best cell.
        ceiling_c (float): The dry bulb ceiling of the wet bulb rule used (`math.inf` when none).
    """

    blue: float
    green: float
    yellow: float
    red: float
    cold: float
    hot: float
    humid: float
    livable_share: float
    livable_area: float
    livable_time: float
    area_always: float
    area_most: float
    area_ever: float
    place_time: np.ndarray
    coldest: Extreme
    hottest: Extreme
    by_band: tuple
    by_column: tuple
    by_quarter: tuple
    sectors: object
    best_site: BestSite
    ceiling_c: float


def wet_bulb_ceiling_c(relative_humidity, limit_c=None):
    """
    The dry bulb temperature at which the wet bulb reaches the limit (`habitability.WET_BULB_CRIT_C`) at this humidity, deg C: a sample at or above it is Yellow. `math.inf` when no
    humidity is given or when the wet bulb stays under the limit up to 50 deg C, the top of the fit of `habitability.wet_bulb_c`.
    """
    if relative_humidity is None:
        return math.inf
    limit_c = habitability.WET_BULB_CRIT_C if limit_c is None else limit_c
    if habitability.wet_bulb_c(50.0, relative_humidity) < limit_c:
        return math.inf
    low, high = -20.0, 50.0
    for _ in range(60):
        middle = 0.5 * (low + high)
        if habitability.wet_bulb_c(middle, relative_humidity) < limit_c:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def band_of_latitude(latitude_deg):
    """The index into `latitude_bands.NAMES` of a latitude (north is 0)."""
    return min(max(int((90.0 - latitude_deg) // latitude_bands.BAND_WIDTH_DEG), 0), latitude_bands.BAND_COUNT - 1)


def longest_out_days(livable, step_days):
    """The longest unbroken run of samples that are not livable, days, in a series that repeats (a run that reaches the end joins the start). 0 when all are livable, the whole
    span when none is."""
    flags = np.asarray(livable, dtype=bool)
    out = ~flags
    if not out.any():
        return 0.0
    if out.all():
        return float(out.size * step_days)
    start = int(np.argmax(flags))
    rolled = np.roll(out, -start)
    change = np.diff(np.concatenate(([0], rolled.astype(np.int8), [0])))
    begins, ends = np.flatnonzero(change == 1), np.flatnonzero(change == -1)
    return float(np.max(ends - begins) * step_days)


def _extreme(series, mask, hottest):
    """The coldest (or hottest) sample of the series among the cells and steps where `mask` is true, as an `Extreme`."""
    temperature = np.where(mask, series.temperature_k, np.inf if not hottest else -np.inf)
    flat = int(np.argmax(temperature) if hottest else np.argmin(temperature))
    row, column, step = np.unravel_index(flat, temperature.shape)
    return Extreme(float(temperature[row, column, step]), int(row), int(column), float(series.grid.latitude_deg[row]), float(series.grid.longitude_deg[column]),
                   float(series.time_days[step]))


def _range(temperature, selector):
    """`(coldest, hottest)` over the cells of `selector` (a boolean `[row, column]` mask) and all steps."""
    chosen = temperature[selector]
    return (float(chosen.min()), float(chosen.max())) if chosen.size else (math.nan, math.nan)


def sector_masks(series):
    """The `[row, column]` masks of the sectors of the sky of the star, `{name: mask}`, by the angle from the point under the star (latitude 0, longitude 0)."""
    latitude = np.radians(series.grid.latitude_deg)[:, None]
    longitude = np.radians(series.grid.longitude_deg)[None, :]
    angle = np.degrees(np.arccos(np.clip(np.cos(latitude) * np.cos(longitude), -1.0, 1.0)))
    return {SECTORS[0]: angle < DAY_SIDE_DEG, SECTORS[1]: (angle >= DAY_SIDE_DEG) & (angle <= NIGHT_SIDE_DEG), SECTORS[2]: angle > NIGHT_SIDE_DEG}


def coverage(series, relative_humidity=None):
    """
    The `Coverage` of a `TemperatureSeries`; `relative_humidity` (0 to 1) switches the wet bulb rule on.
    """
    temperature = series.temperature_k
    rows, columns, steps = temperature.shape
    weights = np.repeat(series.grid.row_weight / columns, columns).reshape(rows, columns)
    celsius = temperature - KELVIN
    band = habitability.TEMPERATURE_BAND
    ceiling = wet_bulb_ceiling_c(relative_humidity)

    def share(mask):
        return float(np.sum(mask.mean(axis=2) * weights))

    humid = (celsius >= band[2]) & (celsius <= band[5]) & (celsius >= ceiling)
    cold = celsius < band[2]
    hot = celsius > band[5]
    blue = (celsius >= band[3]) & (celsius <= band[4]) & ~humid
    green = (celsius >= band[2]) & (celsius <= band[5]) & ~blue & ~humid
    red = (celsius < band[1]) | (celsius > band[6])
    yellow = ~(blue | green | red)
    livable = blue | green
    place_time = livable.mean(axis=2)
    most = place_time >= MOST_OF_THE_TIME
    area_most_of_time = float(np.sum(most * weights))
    livable_time = float(np.sum(place_time * most * weights) / area_most_of_time) if area_most_of_time > 0.0 else 0.0

    by_band = []
    for index in range(latitude_bands.BAND_COUNT):
        selector = np.repeat(np.array([band_of_latitude(lat) == index for lat in series.grid.latitude_deg])[:, None], columns, axis=1)
        by_band.append(_range(temperature, selector))
    by_column = tuple(_range(temperature[:, [c], :], np.ones((rows, 1), dtype=bool)) for c in range(columns))
    phase = (series.time_days % series.orbital_days) / series.orbital_days
    quarter = np.minimum((phase * QUARTERS).astype(int), QUARTERS - 1)
    by_quarter = tuple((float(temperature[:, :, quarter == q].min()), float(temperature[:, :, quarter == q].max())) if np.any(quarter == q) else (math.nan, math.nan)
                       for q in range(QUARTERS))
    sectors = None
    if not math.isfinite(series.solar_days):
        sectors = {}
        for name, selector in sector_masks(series).items():
            area = float(np.sum(weights[selector]))
            low, high = _range(temperature, selector)
            sectors[name] = (float(np.sum(place_time[selector] * weights[selector]) / area) if area > 0.0 else 0.0, low, high)

    best = float(place_time.max())
    tied = np.argwhere(place_time >= best - TIE_TOLERANCE)
    swing = temperature.max(axis=2) - temperature.min(axis=2)
    row, column = min((tuple(cell) for cell in tied), key=lambda cell: (swing[cell], cell))
    cell_series = temperature[row, column]
    site = BestSite(int(row), int(column), float(series.grid.latitude_deg[row]), float(series.grid.longitude_deg[column]), float(place_time[row, column]),
                    float(cell_series.min()), float(cell_series.max()), longest_out_days(livable[row, column], series.step_days), float(cell_series.mean()))
    everywhere = np.ones_like(temperature, dtype=bool)
    return Coverage(
        blue=share(blue), green=share(green), yellow=share(yellow), red=share(red), cold=share(cold), hot=share(hot), humid=share(humid),
        livable_share=share(livable), livable_area=area_most_of_time, livable_time=livable_time,
        area_always=float(np.sum((place_time >= 1.0 - 1.0e-9) * weights)), area_most=float(np.sum((place_time >= NEARLY_ALWAYS) * weights)),
        area_ever=float(np.sum((place_time > 0.0) * weights)), place_time=place_time,
        coldest=_extreme(series, everywhere, False), hottest=_extreme(series, everywhere, True),
        by_band=tuple(by_band), by_column=by_column, by_quarter=by_quarter, sectors=sectors, best_site=site, ceiling_c=ceiling)
