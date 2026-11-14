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

import functools
import math
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np

from planetgen.physics import albedo_color, climate_branches, climate_draws, constants, daylight, habitability, latitude_bands, regional_temperature, ring_darkening, seasonal_swing, sky_suns, surface_cycle

KELVIN = constants.CELSIUS_ZERO_K
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


WET_BULB_TOP_C = 50.0
"""float: The top of the fit of `habitability.wet_bulb_c` (-20 to 50 C): the wet bulb rule is not applied above it, as `habitability` did not."""


def wet_bulb_ceiling_c(relative_humidity, limit_c=None):
    """
    The dry bulb temperature at which the wet bulb reaches the limit (`habitability.WET_BULB_CRIT_C`) at this humidity, deg C: a sample at or above it is Yellow. `math.inf` when no
    humidity is given or when the wet bulb stays under the limit up to 50 deg C, the top of the fit of `habitability.wet_bulb_c`.
    """
    if relative_humidity is None:
        return math.inf
    limit_c = habitability.WET_BULB_CRIT_C if limit_c is None else limit_c
    if habitability.wet_bulb_c(WET_BULB_TOP_C, relative_humidity) < limit_c:
        return math.inf
    low, high = -20.0, WET_BULB_TOP_C
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


def _substellar_angle_deg(series):
    """The `[row, column]` angle of each cell from the point under the star (latitude 0, longitude 0), degrees."""
    latitude = np.radians(series.grid.latitude_deg)[:, None]
    longitude = np.radians(series.grid.longitude_deg)[None, :]
    return np.degrees(np.arccos(np.clip(np.cos(latitude) * np.cos(longitude), -1.0, 1.0)))


def _cell_weights(series):
    """The `[row, column]` area weight of each cell: `grid.row_weight / columns` (equal area along a row)."""
    rows, columns = series.temperature_k.shape[:2]
    return np.repeat(series.grid.row_weight / columns, columns).reshape(rows, columns)


def _latitude_band_masks(series):
    """`[row, column]` masks of the seven bands of latitude (`latitude_bands.NAMES`), north to south."""
    columns = series.temperature_k.shape[1]
    index = np.repeat(np.array([band_of_latitude(lat) for lat in series.grid.latitude_deg])[:, None], columns, axis=1)
    return [index == number for number in range(latitude_bands.BAND_COUNT)]


def sector_masks(series):
    """The `[row, column]` masks of the sectors of the sky of the star, `{name: mask}`, by the angle from the point under the star (latitude 0, longitude 0)."""
    angle = _substellar_angle_deg(series)
    return {SECTORS[0]: angle < DAY_SIDE_DEG, SECTORS[1]: (angle >= DAY_SIDE_DEG) & (angle <= NIGHT_SIDE_DEG), SECTORS[2]: angle > NIGHT_SIDE_DEG}


def coverage(series, relative_humidity=None, band=None, wet_bulb_limit=None):
    """
    The `Coverage` of a `TemperatureSeries`; `relative_humidity` (0 to 1) switches the wet bulb rule on. `band` is the eight edges of a visitor's temperature bands
    (`habitability.Bands.temperature`; the human ones when `None`) and `wet_bulb_limit` its uncompensable wet bulb (`habitability.Bands.wet_bulb`).
    """
    temperature = series.temperature_k
    rows, columns, steps = temperature.shape
    weights = _cell_weights(series)
    celsius = temperature - KELVIN
    band = habitability.TEMPERATURE_BAND if band is None else band
    ceiling = wet_bulb_ceiling_c(relative_humidity, wet_bulb_limit)

    total = float(np.sum(weights))

    def share(mask):
        return float(np.sum(mask.mean(axis=2) * weights) / total)

    humid = (celsius >= band[2]) & (celsius <= band[5]) & (celsius <= WET_BULB_TOP_C) & (celsius >= ceiling)
    cold = celsius < band[2]
    hot = celsius > band[5]
    blue = (celsius >= band[3]) & (celsius <= band[4]) & ~humid
    green = (celsius >= band[2]) & (celsius <= band[5]) & ~blue & ~humid
    red = (celsius < band[1]) | (celsius > band[6])
    yellow = ~(blue | green | red)
    livable = blue | green
    place_time = livable.mean(axis=2)
    most = place_time >= MOST_OF_THE_TIME
    area_most_of_time = float(np.sum(most * weights) / total)
    livable_time = float(np.sum(place_time * most * weights) / (area_most_of_time * total)) if area_most_of_time > 0.0 else 0.0

    by_band = [_range(temperature, selector) for selector in _latitude_band_masks(series)]
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
        area_always=float(np.sum((place_time >= 1.0 - 1.0e-9) * weights) / total), area_most=float(np.sum((place_time >= NEARLY_ALWAYS) * weights) / total),
        area_ever=float(np.sum((place_time > 0.0) * weights) / total), place_time=place_time,
        coldest=_extreme(series, everywhere, False), hottest=_extreme(series, everywhere, True),
        by_band=tuple(by_band), by_column=by_column, by_quarter=by_quarter, sectors=sectors, best_site=site, ceiling_c=ceiling)


# --- From a body's stored columns, and what is stored (GEN.1265) ----------------------------------------------------------------------------------------------------------

FIELDS = (
    "temp_livable_share", "temp_livable_area", "temp_livable_time", "temp_blue_share", "temp_green_share", "temp_cold_share", "temp_hot_share", "temp_area_always", "temp_area_most",
    "temp_livable_belt_only",
    "temp_coldest_k", "temp_coldest_lat_deg", "temp_coldest_lon_deg", "temp_coldest_day", "temp_hottest_k", "temp_hottest_lat_deg", "temp_hottest_lon_deg", "temp_hottest_day",
    "temp_site_lat_deg", "temp_site_lon_deg", "temp_site_livable_time", "temp_site_coldest_k", "temp_site_hottest_k", "temp_site_mean_k", "temp_site_longest_out_days",
)
"""tuple: The columns stored on `planets` and `moons` (GEN.1265), all for the human visitor: the livable share L and its two pieces, the Blue, Green, cold and hot shares, the area livable
all the time and at least `NEARLY_ALWAYS` of it, whether the livable ground is only the terminator belt (0 or 1), the coldest and hottest sample anywhere with where (degrees) and
when (days from an equinox), and the best site (its place, the share of the time it is livable, its coldest, hottest and mean temperature, kelvin, and its longest spell out of the
livable band, days). All `None` for a body with no temperature or no habitability score."""

BELT_ONLY_SHARE = 0.1
"""float: A world that turns once an orbit is livable only in its terminator belt when the livable share of the day side and of the night side are each under this ("livable only in
the terminator belt", Boss 7:32 PM EST 11-04). Source tier T9, tweak level L3 (a stand-in, set on the sweep of GEN.1265)."""

SHEET_ROWS = 12
"""int: The rows of the grid a still body (no spin or starlight to work from) is laid on; it is the cell grid of `surface_cycle`, so its shares read as any body's."""


def still_series(temperature_k):
    """
    The series of a body that has one temperature everywhere and always (no star to light it, or no spin or orbit stored): every cell of the grid at `temperature_k` for four steps. A
    hand-built `habitability.World` is read this way too.
    """
    grid = surface_cycle.cell_grid()
    rows, columns = grid.latitude_deg.size, grid.longitude_deg.size
    steps = 4
    return SimpleNamespace(grid=grid, temperature_k=np.full((rows, columns, steps), float(temperature_k)), time_days=np.arange(steps) + 0.5, step_days=1.0,
                           orbital_days=float(steps), solar_days=1.0)


def starlight_w_m2(values):
    """
    The starlight a stored body takes, W/m2, and the Bond albedo, from the columns that fix its temperature: `T_s = ((1 - A) S (1 + G) / 4 sigma)^(1/4)`
    (`climate_draws.surface_temperature_k`), so `S = 4 sigma T_s^4 / ((1 - A)(1 + G))`. `(None, None)` for a body without the stored albedo and greenhouse quantity, with no air
    model, or on the 2.725 K floor (no light to work from).
    """
    albedo, greenhouse, mean = values.get("bond_albedo"), values.get("greenhouse_factor"), values.get("surface_temperature_k")
    if albedo is None or greenhouse is None or not mean or not 0.0 <= albedo < 1.0:
        return None, None
    if mean <= constants.COSMIC_BACKGROUND_TEMPERATURE_K * (1.0 + climate_draws.FLOOR_TOLERANCE):
        return None, None
    return 4.0 * constants.STEFAN_BOLTZMANN_CONSTANT * mean ** 4 / ((1.0 - albedo) * (1.0 + greenhouse)), albedo


def ice_law_of(values, albedo):
    """
    The ice-albedo law (`climate_branches.IceAlbedoLaw`, GEN.972) of a stored body with surface water, or `None`: the ice-covered albedo is the stand-in of `climate_branches`
    turned to the colour of the star; the ice-free one is the albedo that gives the stored Bond albedo `albedo` with the stored share of ground under ice (`ground_ice_fraction`), so the
    area mean of the cells' albedo is the stored one.
    """
    water, mass = values.get("water_mass_fraction"), values.get("mass_kg")
    if not water or not mass or water <= 0.0 or mass <= 0.0:
        return None
    covered = max(albedo_color.shifted_range((climate_branches.ICE_COVERED_ALBEDO,) * 2, "frozen", values.get("star_teff_k"))[0], albedo)
    ice = min(max(values.get("ground_ice_fraction") or 0.0, 0.0), 0.95)
    free = min(max((albedo - covered * ice) / (1.0 - ice), 0.0), albedo)
    return climate_branches.IceAlbedoLaw(ice_free=free, ice_covered=max(covered, free))


def ring_transmission(values):
    """
    `transmission(latitudes_deg, declinations_rad)` of `surface_cycle.series` for a ringed planet that turns on its own (GEN.1055 R8, `ring_darkening.shadow_transmission`), or
    `None` for a body with no ring that shades. The ring lies in the equatorial plane, so the Sun's declination is the angle the shadow is worked out at; a world with no tilt sees
    the ring edge on. The share is read at declinations rounded to a quarter degree.
    """
    inner_km, outer_km = values.get("ring_inner_km"), values.get("ring_outer_km")
    depth, radius_km = values.get("ring_optical_depth"), values.get("radius_km")
    if not (inner_km and outer_km and depth and radius_km) or depth <= 0.0 or radius_km <= 0.0 or outer_km < inner_km:
        return None
    if math.sin(math.radians(daylight.folded_tilt_deg(values.get("axial_tilt_deg") or 0.0))) < 1.0e-9:
        return None
    inner, outer = inner_km / radius_km, outer_km / radius_km
    cache = {}

    def transmission(latitudes_deg, declinations_rad):
        degrees = np.round(np.degrees(declinations_rad) * 4.0) / 4.0
        out = np.ones((len(latitudes_deg), len(degrees)))
        for column, declination in enumerate(degrees):
            for row, latitude in enumerate(latitudes_deg):
                key = (row, float(declination))
                if key not in cache:
                    cache[key] = ring_darkening.shadow_transmission(float(latitude), float(declination), inner, outer, depth, RING_LONGITUDES)
                out[row, column] = cache[key]
        return out

    return transmission


RING_LONGITUDES = 24
"""int: The longitudes that sample a turn of the planet for the ring's shade at one latitude (`latitude_bands` uses the same)."""


def _inputs(values):
    """The hashable key of what the series of `values` is made from, or `None` when a series cannot be made (no spin, no orbit, or no stored albedo to work the starlight from)."""
    rotation, period, mean = values.get("rotation_period_hours"), values.get("period_years"), values.get("surface_temperature_k")
    if not rotation or rotation <= 0.0 or not period or period <= 0.0 or not mean or mean <= 0.0:
        return None
    flux, albedo = starlight_w_m2(values)
    if flux is None or flux <= 0.0:
        return None
    ring = tuple(values.get(name) for name in ("ring_inner_km", "ring_outer_km", "ring_optical_depth", "radius_km"))
    water = (values.get("water_mass_fraction"), values.get("mass_kg"), values.get("star_teff_k"), values.get("ground_ice_fraction"))
    solar_day, days = values.get("solar_day_hours"), values.get("days_per_year")
    year_h = solar_day * days if solar_day and days else period * _YEAR_H       # a moon is lit through its planet's year, not its own orbit of the planet
    return (rotation, year_h, values.get("axial_tilt_deg") or 0.0, mean, flux, albedo, values.get("atmospheric_pressure_pa") or 0.0,
            (values.get("gravity_g") or 0.0) * constants.EARTH_GRAVITY, values.get("ground_ocean_fraction") or 0.0, ring, water, values.get("sky_suns") or None)


_YEAR_H = constants.HOURS_PER_YEAR
"""float: Hours in a year of 365.25 days, the unit of the orbital period `temperature_series` takes (`period_years` is stored in years)."""


@functools.lru_cache(maxsize=64)
def _series(key):
    rotation, orbital, tilt, mean, flux, albedo, pressure, gravity, ocean, ring, water, suns_text = key
    values = {"ring_inner_km": ring[0], "ring_outer_km": ring[1], "ring_optical_depth": ring[2], "radius_km": ring[3], "axial_tilt_deg": tilt,
              "water_mass_fraction": water[0], "mass_kg": water[1], "star_teff_k": water[2], "ground_ice_fraction": water[3]}
    law = ice_law_of(values, albedo)
    locked = regional_temperature.is_locked(rotation, orbital / _YEAR_H)
    return seasonal_swing.temperature_series(
        rotation, orbital, tilt, mean, flux, albedo, pressure, gravity, ocean, ice_fraction=water[3], ice_law=law,
        transmission=None if locked else ring_transmission(values), suns=sky_suns.decode(suns_text))        # GEN.1458: the other suns of the sky


def series_of_values(values):
    """
    The `TemperatureSeries` of a rocky body from its column-named `values` (`habitability_world.values_from_body`, or a stored row): the rotation period, the orbit (`period_years`),
    the axial tilt, the mean temperature, the air (`atmospheric_pressure_pa`, `gravity_g`), the sea and the ice (`ground_ocean_fraction`, `ground_ice_fraction`, `water_mass_fraction`),
    the albedo and greenhouse quantity that fix the starlight (`starlight_w_m2`) and the ring (`ring_*`). `None` for a body that cannot be worked (no spin or orbit stored, or no
    starlight): its temperature is then one everywhere (`still_series`). Kept for the last few bodies, so a re-score does not work the same series twice.
    """
    key = _inputs(values)
    return None if key is None else _series(key)


@functools.lru_cache(maxsize=64)
def _coverage(key, humidity, band, wet_bulb_limit):
    return coverage(_series(key), humidity, band, wet_bulb_limit)


def coverage_of_values(values, relative_humidity=None, band=None, wet_bulb_limit=None):
    """
    The `Coverage` of a body from its stored columns for a visitor whose temperature bands are `band` (the human ones when `None`), or of its one temperature (`still_series`) when
    no series can be made. Kept for the last few bodies.
    """
    key = _inputs(values)
    band = tuple(habitability.TEMPERATURE_BAND if band is None else band)
    if key is None:
        return coverage(still_series(values["surface_temperature_k"]), relative_humidity, band, wet_bulb_limit)
    return _coverage(key, relative_humidity, band, wet_bulb_limit)


def belt_only(read):
    """Whether the livable ground of `read` (a `Coverage`) is only the terminator belt of a world that turns once an orbit: some of it is livable and each of the day side and the night
    side is under `BELT_ONLY_SHARE`."""
    if read.sectors is None or read.livable_share <= 0.0:
        return False
    return read.sectors[SECTORS[0]][0] < BELT_ONLY_SHARE and read.sectors[SECTORS[2]][0] < BELT_ONLY_SHARE


def fields(read):
    """`{FIELDS name: value}` of a `Coverage`; every field `None` for `None`."""
    if read is None:
        return dict.fromkeys(FIELDS)
    site = read.best_site
    return {
        "temp_livable_share": read.livable_share, "temp_livable_area": read.livable_area, "temp_livable_time": read.livable_time, "temp_blue_share": read.blue,
        "temp_green_share": read.green, "temp_cold_share": read.cold, "temp_hot_share": read.hot, "temp_area_always": read.area_always, "temp_area_most": read.area_most,
        "temp_livable_belt_only": int(belt_only(read)),
        "temp_coldest_k": read.coldest.temperature_k, "temp_coldest_lat_deg": read.coldest.latitude_deg, "temp_coldest_lon_deg": read.coldest.longitude_deg,
        "temp_coldest_day": read.coldest.time_days,
        "temp_hottest_k": read.hottest.temperature_k, "temp_hottest_lat_deg": read.hottest.latitude_deg, "temp_hottest_lon_deg": read.hottest.longitude_deg,
        "temp_hottest_day": read.hottest.time_days,
        "temp_site_lat_deg": site.latitude_deg, "temp_site_lon_deg": site.longitude_deg, "temp_site_livable_time": site.livable_time, "temp_site_coldest_k": site.coldest_k,
        "temp_site_hottest_k": site.hottest_k, "temp_site_mean_k": site.mean_k, "temp_site_longest_out_days": site.longest_out_days,
    }


# --- The seven bands of latitude from the same series (GEN.1265) ------------------------------------------------------------------------------------------------------

def _band_masks(series):
    """`[row, column]` masks of the seven bands, north to south: by latitude, or for a world that turns once an orbit the rings by the angle from the point under the star
    (`latitude_bands.RING_NAMES`)."""
    if math.isfinite(series.solar_days):
        return _latitude_band_masks(series)
    index = np.minimum((_substellar_angle_deg(series) // latitude_bands.BAND_WIDTH_DEG).astype(int), latitude_bands.BAND_COUNT - 1)
    return [index == number for number in range(latitude_bands.BAND_COUNT)]


def band_climate(series):
    """
    `[(mean_k, swing_k)]` for each of the seven bands of `latitude_bands` (north to south, or day center to night center for a world that turns once an orbit), read off the series: the
    mean of the band is the mean over its cells (by area) and the whole span, the swing the day and night swing, peak to peak, of the band's cells (the mean over its cells and its solar
    days of the hottest minus the coldest sample of the day). For a world that turns once an orbit there is no day, and a ring's swing is the range of its cells over the span, hottest minus
    coldest (the ring's real range).
    """
    temperature = series.temperature_k
    rows, columns, steps = temperature.shape
    weights = _cell_weights(series)
    mean = temperature.mean(axis=2)
    if math.isfinite(series.solar_days):
        # PERF.249: a world that turns many times per sample (a hypergiant's planet) has more solar days than samples; the windows past the last sample are empty and were
        # dropped, so the count is held at the number of samples (it was split into hundreds of millions of empty windows) and the windows are reduced in one pass
        days = min(max(int(round(steps * series.step_days / series.solar_days)), 1), steps)
        size, extra = divmod(steps, days)
        lengths = np.full(days, size)
        lengths[:extra] += 1
        starts = np.concatenate(([0], np.cumsum(lengths)[:-1]))
        spans = np.maximum.reduceat(temperature, starts, axis=2) - np.minimum.reduceat(temperature, starts, axis=2)
        daily = np.mean(np.ascontiguousarray(np.moveaxis(spans, 2, 0)), axis=0)
    out = []
    for mask in _band_masks(series):
        if not mask.any():
            out.append((math.nan, math.nan))
            continue
        share = weights[mask]
        if math.isfinite(series.solar_days):
            swing = float(np.sum(daily[mask] * share) / np.sum(share))
        else:
            swing = float(temperature[mask].max() - temperature[mask].min())
        out.append((float(np.sum(mean[mask] * share) / np.sum(share)), swing))
    return out


def summary(values):
    """
    The `temperature_coverage` object the API serves (`None` for a body that has none): `livable_share`, `livable_area`, `livable_time`, `blue_share`, `green_share`, `cold_share`,
    `hot_share`, `area_always`, `area_most`, `livable_belt_only` (a bool), `coldest` and `hottest` (`{"temperature_k", "latitude_deg", "longitude_deg", "day"}`: the worst conditions
    anywhere, with where and when) and `best_site` (`{"latitude_deg", "longitude_deg", "livable_time", "coldest_k", "hottest_k", "mean_k", "longest_out_days"}`), from the stored columns
    (`FIELDS`).
    """
    if values.get("temp_livable_share") is None:
        return None
    return {
        "livable_share": values["temp_livable_share"], "livable_area": values.get("temp_livable_area"), "livable_time": values.get("temp_livable_time"),
        "blue_share": values.get("temp_blue_share"), "green_share": values.get("temp_green_share"), "cold_share": values.get("temp_cold_share"),
        "hot_share": values.get("temp_hot_share"), "area_always": values.get("temp_area_always"), "area_most": values.get("temp_area_most"),
        "livable_belt_only": bool(values.get("temp_livable_belt_only")),
        "coldest": {"temperature_k": values.get("temp_coldest_k"), "latitude_deg": values.get("temp_coldest_lat_deg"), "longitude_deg": values.get("temp_coldest_lon_deg"),
                    "day": values.get("temp_coldest_day")},
        "hottest": {"temperature_k": values.get("temp_hottest_k"), "latitude_deg": values.get("temp_hottest_lat_deg"), "longitude_deg": values.get("temp_hottest_lon_deg"),
                    "day": values.get("temp_hottest_day")},
        "best_site": {"latitude_deg": values.get("temp_site_lat_deg"), "longitude_deg": values.get("temp_site_lon_deg"), "livable_time": values.get("temp_site_livable_time"),
                      "coldest_k": values.get("temp_site_coldest_k"), "hottest_k": values.get("temp_site_hottest_k"), "mean_k": values.get("temp_site_mean_k"),
                      "longest_out_days": values.get("temp_site_longest_out_days")},
    }
