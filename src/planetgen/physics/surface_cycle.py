# planetgen/physics/surface_cycle.py

"""
The sunlight on every place of a body through its year and its turns (GEN.1262)
===============================================================================

One calculation for every body, whatever it does: the cosine of the angle between the local vertical and the direction of the star, for each cell of a grid of
latitudes and longitudes and each step of time. A tidally locked world, a world that turns backwards, a 28-day spinner near resonance and Earth all come out of the same
formula; there is no case for any of them.

**The geometry.** The star's declination `d` (its latitude on the sky of the planet) follows the planet round a circular orbit, `sin d = sin(tilt) sin(L)`, `L = 2 pi t / P` the
position in the orbit counted from an equinox, `P` the orbital period; a tilt past 90 degrees is a retrograde spin with the same seasons (`daylight.folded_tilt_deg`). The star's
hour angle `h` at a longitude `lambda` is `2 pi t / D + lambda`, `D` the solar day: the time from one noon to the next, which is `1 / |1 / R - 1 / P|` for a body that turns
prograde in `R` and `1 / (1 / R + 1 / P)` for one that turns backwards, and endless for a body that turns once an orbit, whose star then stands still in the sky of every
cell. The cosine of the zenith angle is

    cos z = sin(phi) sin(d) + cos(phi) cos(d) cos(h)

and the star's upper limb is up while `cos z > -sin(alpha)`, `alpha` the star's angular radius (`daylight.day_share_at` is the same law integrated over a turn).

**The grid.** `CELL_ROWS` rows of equal area (equal steps in the sine of the latitude) by `CELL_COLUMNS` columns of equal longitude, so every cell of a row has one weight and
a world average is the plain mean of the rows. The seven latitude bands (`latitude_bands`) are a coarse subset of it.

**The time.** One orbit (or one solar day, whichever is longer) in steps of at most a `SAMPLES_PER_SOLAR_DAY`-th of a solar day, a `STEPS_PER_TIME_CONSTANT`-th of the thermal time
constant of the surface when the caller has one, and a `STEPS_PER_YEAR`-th of the year. The longest of the three limits that fits in `MAX_STEPS` is taken.

**Shade.** A ring (GEN.1055) takes sunlight off the cells behind it; `series` takes an optional `transmission(latitudes_deg, declinations_rad)` that returns the light left for each
latitude and each step, and multiplies it into the cosine where the star is up.

Limits. The orbit is circular, so the declination is a sine of the orbital angle; the axis does not precess; the star is a disc without refraction; the cosine is the instantaneous
one, not the day's mean.
"""

import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import daylight

CELL_ROWS = 12
"""int: The rows of latitude of the grid, each of the same area. [GEN.1262; resolution of the grid, a numerical choice from the plan of Boss, 7:32 PM EST 2026-11-04.]"""

CELL_COLUMNS = 12
"""int: The columns of longitude of the grid. Same source as `CELL_ROWS`."""

SAMPLES_PER_SOLAR_DAY = 48
"""int: The least time steps of one solar day when `MAX_STEPS` allows them, so a day's rise and set are found to a 48th of it (30 minutes of an Earth day)."""

STEPS_PER_TIME_CONSTANT = 4
"""int: The least time steps of the thermal time constant of the surface, when the caller gives it."""

STEPS_PER_YEAR = 3000
"""int: The least time steps of one orbit."""

MAX_STEPS = 20000
"""int: The most time steps of the series, which bounds the memory (a cell of the grid is `MAX_STEPS` floats) and the time of a body."""

MAX_PERIODIC_ORBITS = 12
"""int: The most orbits a periodic series may span to make a whole number of solar days of them. [GEN.1263; numerical.]"""

PERIODIC_TOLERANCE = 0.02
"""float: How far the solar day may be stretched to fit a whole number of them in the span, as a share of the day (and, for a body whose star hardly moves in its sky, the share of a
turn the star may move in the span). [GEN.1263; numerical.]"""


@dataclass(frozen=True)
class Grid:
    """The cells: the latitude of each row and the longitude of each column, degrees, and the share of the surface each row holds (one cell holds `row_weight / columns`)."""

    latitude_deg: np.ndarray
    longitude_deg: np.ndarray
    row_weight: np.ndarray


@dataclass(frozen=True)
class Series:
    """The sunlight of a body: `cosine[row, column, step]` is the cosine of the zenith angle (negative when the star is below the horizon, scaled by the transmission where
    there is one), `time_days` the time of each step from an equinox, `step_days` the step, and `declination_rad` the star's declination at each step."""

    grid: Grid
    time_days: np.ndarray
    step_days: float
    declination_rad: np.ndarray
    cosine: np.ndarray
    orbital_days: float
    solar_days: float


def cell_grid(rows=CELL_ROWS, columns=CELL_COLUMNS):
    """The grid: `rows` rows of equal area and `columns` columns of equal longitude."""
    sine = -1.0 + (np.arange(rows) + 0.5) * 2.0 / rows
    latitude = np.degrees(np.arcsin(sine))
    longitude = (np.arange(columns) + 0.5) * 360.0 / columns
    return Grid(latitude, longitude, np.full(rows, 1.0 / rows))


def solar_day_days(rotation_hours, orbital_hours, tilt_deg):
    """The time from one noon to the next, days: `1 / |1 / R - 1 / P|` for a body that turns prograde (tilt up to 90 degrees), `1 / (1 / R + 1 / P)` for one that turns backwards,
    and `inf` for a body that turns once an orbit (or has no orbit): the star stands still in its sky."""
    if not rotation_hours or rotation_hours <= 0.0:
        return math.inf
    if not orbital_hours or orbital_hours <= 0.0:
        return rotation_hours / 24.0
    prograde = abs(tilt_deg or 0.0) % 360.0 <= 90.0 or abs(tilt_deg or 0.0) % 360.0 >= 270.0
    rate = 1.0 / rotation_hours + (-1.0 if prograde else 1.0) / orbital_hours
    return 1.0 / (24.0 * abs(rate)) if abs(rate) > 1.0e-12 else math.inf


def periodic_span(orbital_days, solar_days):
    """
    The span and the solar day that make the sunlight repeat: `(span_days, solar_days)`, the span a whole number of orbits (the declination repeats) holding a whole number of
    solar days (the hour angle repeats), the solar day stretched by at most `PERIODIC_TOLERANCE` of itself. A series cut at an orbit when a solar day is 14 days and the orbit 28.5
    would end in a jump; one that spans the orbit and holds exactly two days does not. A star that stands still (`inf`) or moves less than `PERIODIC_TOLERANCE` of a turn in
    the span stays still. When no span of up to `MAX_PERIODIC_ORBITS` fits, the one that fits best is taken.
    """
    if not math.isfinite(solar_days):
        return orbital_days, math.inf
    turns_per_orbit = orbital_days / solar_days
    best = None
    for orbits in range(1, MAX_PERIODIC_ORBITS + 1):
        turns = round(orbits * turns_per_orbit)
        error = abs(orbits * turns_per_orbit - turns) / max(turns, 1)
        if best is None or error < best[0] - 1.0e-12:
            best = (error, orbits, turns)
        if error <= PERIODIC_TOLERANCE:
            best = (error, orbits, turns)
            break
    _, orbits, turns = best
    span = orbits * orbital_days
    return span, (span / turns if turns > 0 else math.inf)


def time_step_days(orbital_days, solar_days, thermal_days=None, samples_per_day=SAMPLES_PER_SOLAR_DAY, max_steps=MAX_STEPS, span=None, steps_per_year=STEPS_PER_YEAR):
    """The step and the number of steps of the series: the span is the longer of the orbit and a finite solar day, the step the least of a `samples_per_day`-th of the solar day,
    a quarter of the thermal time constant (when given) and a 3000th of the orbit, but never so small the span needs more than `max_steps`."""
    span = span or max(orbital_days, solar_days if math.isfinite(solar_days) else 0.0)
    wanted = min(solar_days / samples_per_day if math.isfinite(solar_days) else math.inf, orbital_days / steps_per_year)
    if thermal_days:
        wanted = min(wanted, thermal_days / STEPS_PER_TIME_CONSTANT)
    count = int(min(max(math.ceil(span / wanted), 2), max_steps))
    return span / count, count


def declination_rad(time_days, orbital_days, tilt_deg):
    """The star's declination at each time, radians: `sin d = sin(tilt) sin(2 pi t / P)`, from an equinox."""
    return np.arcsin(math.sin(math.radians(daylight.folded_tilt_deg(tilt_deg))) * np.sin(2.0 * math.pi * np.asarray(time_days) / orbital_days))


def cosine_zenith(latitude_deg, longitude_deg, time_days, orbital_days, solar_days, tilt_deg):
    """`cos z` for each latitude, longitude and time: an array of shape (latitudes, longitudes, times). The star's hour angle is `2 pi t / D + longitude`, so a longitude is a
    starting angle of the turn (the phase), and `solar_days` may be `inf` for a body that turns once an orbit."""
    latitude = np.radians(np.atleast_1d(np.asarray(latitude_deg, dtype=float)))[:, None, None]
    longitude = np.radians(np.atleast_1d(np.asarray(longitude_deg, dtype=float)))[None, :, None]
    time = np.asarray(time_days, dtype=float)[None, None, :]
    declination = declination_rad(time, orbital_days, tilt_deg)
    hour = 2.0 * math.pi * time / solar_days + longitude
    return np.sin(latitude) * np.sin(declination) + np.cos(latitude) * np.cos(declination) * np.cos(hour)


def series(rotation_hours, orbital_hours, tilt_deg, thermal_days=None, transmission=None, rows=CELL_ROWS, columns=CELL_COLUMNS,
           samples_per_day=SAMPLES_PER_SOLAR_DAY, max_steps=MAX_STEPS, periodic=False, steps_per_year=STEPS_PER_YEAR):
    """The `Series` of a body from its stored rotation period, orbital period (both hours) and tilt (degrees). `thermal_days` is the thermal time constant of the surface if the
    caller has one; `transmission(latitudes_deg, declinations_rad)`, if given, returns the share of the light left, shape (rows, steps), and multiplies the cosine of the cells
    where the star is up. `periodic` stretches the solar day a little and spans several orbits where it must (`periodic_span`) so the series repeats exactly: the form a Fourier
    filter (GEN.1263) needs."""
    orbital_days = orbital_hours / 24.0
    solar_days = solar_day_days(rotation_hours, orbital_hours, tilt_deg)
    span = None
    if periodic:
        span, solar_days = periodic_span(orbital_days, solar_days)
    step, count = time_step_days(orbital_days, solar_days, thermal_days, samples_per_day, max_steps, span, steps_per_year)
    grid = cell_grid(rows, columns)
    time = (np.arange(count) + 0.5) * step
    cosine = cosine_zenith(grid.latitude_deg, grid.longitude_deg, time, orbital_days, solar_days, tilt_deg)
    declination = declination_rad(time, orbital_days, tilt_deg)
    if transmission is not None:
        left = np.asarray(transmission(grid.latitude_deg, declination), dtype=float)[:, None, :]
        cosine = np.where(cosine > 0.0, cosine * left, cosine)
    return Series(grid, time, step, declination, cosine, orbital_days, solar_days)
