# planetgen/physics/seasonal_swing.py

"""
The seasonal swing of the temperature at a latitude, from the planet's tilt and orbit (GEN.1039)
================================================================================================

`regional_temperature` spreads a world's mean temperature over latitude bands by the annual-mean light of each band, so the stored tilt already moves the
*means* (a pole is warmer than the equator over the year above 54.8 degrees of tilt), but nothing says how far a place swings between its summer and its winter. A
tilted planet's seasons are its weather: Uranus-like tilts give a pole in daylight at noon-equator brightness for decades, Mars's eccentric orbit adds a second
season, and the Earth without its Moon is the middle case ([R4708]:254; [R4709], 2003, Int. J. Astrobiol. 1:61 and
2:1; [R4171]:596, Icarus; [R0149]:277; recalled, none opened). This module is that
calculation as a pure function of stored quantities. Nothing calls it yet (GEN.1048 wires it into the regional temperatures).

**The light.** The day's mean insolation at latitude `phi` for the Sun at declination `delta` is the closed form of `regional_temperature._daily_insolation`
(`Q = (h sin phi sin delta + cos phi cos delta sin h) / pi`, `h` the half day length). The declination follows the planet round its orbit: `sin delta = sin(eps)
sin(L)`, with `L = v + varpi` the true longitude (`v` the true anomaly by Kepler's equation, `varpi` the longitude of periapsis from the equinox) and the light goes
as `(a / r)^2`. The year is cut into `STEPS` equal times, so an eccentric orbit spends longer at apoapsis. The annual mean over the sphere is `S / (4 sqrt(1 - e^2))`
(checked in the tests).

**The response.** About the place's annual-mean state, one energy balance per latitude (the `energy balance model` of [R3313], J. Atmos. Sci. 36:1189, local
form, recalled):

    C dT/dt = (1 - A) S [Q(t) - <Q>]  -  eps_out sigma (T^4 - Tm^4)  -  B_t (T - Tm)

with `Tm` the place's annual mean (from `regional_temperature`'s law: the fourth root of its annual-mean light, the air carrying the share `transport_share` of the
contrast poleward), `eps_out` the planet's outgoing emissivity (its outgoing long-wave over `sigma Tm^4`, closed by the global mean: absorbed `(1 - A) S / 4 sqrt(1 - e^2)`
over `sigma T^4`, so a greenhouse is a smaller `eps_out`), and `B_t` the extra damping by heat the air carries away, which gives the means' own share
(`regional_temperature.transport_share`): `B_t = 4 eps_out sigma Tm^3 s / (1 - s)`. Solved for the periodic year by implicit steps, starting from the exact
solution of the linearised equation (harmonic by harmonic, `T_n = F_n / (B + i n omega C)`) and repeating until the year closes to a hundredth of a kelvin, so the
nonlinear emission of a cold dry planet (Mars) is kept. The mean of the year can come out a little off `Tm` in a strong nonlinearity; the swing is reported about the
solution's own mean and shifted to `Tm`.

**The heat capacity** `C` of the surface layer the season reaches (per square metre): the ocean's mixed layer (`OCEAN_MIXED_LAYER_M` of water, the share
`ocean_fraction` of the place), the ground's skin on the rest (`regional_temperature.GROUND_INERTIA` times `sqrt(period / pi)`) and the air column
(`regional_temperature.AIR_CP_J_KG_K` times pressure over gravity). The Earth's ocean makes it two orders larger than Mars's ground: a deep sea damps the seasons to a few
kelvin, a dry planet follows the light.

**What it does not do.** No ice-albedo feedback or snow (GEN.972 owns it: the albedo is the caller's, one number), no latent heat, no land and sea laid out by longitude, no
precession of the axis, and the tilt is taken as constant. Days are taken short against the year, so the daily mean applies; a planet whose day is a large share of its year
(a slow rotator) should be given the daily mean of its own day by the caller. The response is a zonal mean: a continental interior swings more than the value here.
Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.91.

The temperature of every place (GEN.1263)
=========================================

`temperature_series` is the same energy balance with the longitude and the time of day added: the temperature of each cell of `surface_cycle`'s grid (12 rows of equal area by 12
columns) at each step of a span that holds a whole number of orbits and of solar days, from the stored rotation period, orbital period, tilt, mean temperature, flux, albedo, air,
gravity and the shares of ocean and ice. A tidally locked world, a world that turns backwards, a 28-day spinner near resonance and Earth come out of the one calculation.

    C(w) dT/dt = a (1 - A) I(t) - eps sigma T^4 - h (T - T_world)

`I` is the sunlight of the cell (`surface_cycle`, the cosine of the zenith angle where the star is up), `A` the Bond albedo of the cell, `a` the scale that makes the area-and-time
mean of the absorbed light the world's own `(1 - A_world) S / 4`, and `h = 4 eps sigma Tm^3 s / (1 - s)` the heat the air carries toward the world's mean temperature at the same
step, with `s = regional_temperature.transport_share` (the locked share, `LOCKED_TRANSPORT_X0`, for a body that turns once an orbit). The heat capacity depends on the cycle: the
layer a cycle reaches is thin for a day's and thick for a year's, so `C(w)` is the air column, the sea to the depth `sqrt(kappa P / pi)` (at most `OCEAN_MIXED_LAYER_M`) on the
share `ocean_fraction` of the cell, and the ground's skin `GROUND_INERTIA sqrt(P / pi)` on the rest, for each harmonic of period `P`. That is a filter in the frequency domain:

1. Each cell's steady temperature balances its mean light with the air's share (Newton's method on `eps sigma T^4 + h (T - Tm) = a (1 - A) <I>`), shifted to the stored mean.
2. The ice is laid by place: the GEN.972 law (`climate_branches.IceAlbedoLaw`, its width and exponent) is moved along the temperature axis until the area mean of the cover
   is the stored `ground_ice_fraction`, then read at the cell's steady temperature; the albedo of the cell is the law's. Twice, since the ice changes the light and the light the ice.
3. The response is the filter `1 / (B + h + i w C(w))` on the Fourier series of the light (`B = 4 eps sigma T_cell^3` the slope of the emission at the cell's mean, the transport
   solved together with the world's mean at each harmonic). What the linear filter leaves out of the fourth-power emission is put right by passes of the same filter on the
   remaining imbalance, accelerated by Anderson mixing (the filter is a preconditioner: a poor slope costs passes, not accuracy). A cold dry night needs 15 to 30 passes, a mild
   world 3 to 5.
4. The emission `eps` is scaled until the mean of the series is the stored mean (the hot and cold places make the mean of the fourth powers more than the fourth power of the
   mean), and what is left is added as a shift, so the area-and-time mean is the stored mean exactly. Against a time-stepped implicit integration of the same equation with one
   heat capacity, the filter agrees to 0.3 K rms for a 0.3 bar Nitrikin-like world and 0.8 K for an Earth-like one (the tests hold the first); for a Mars-like one the two differ by
   about 3 K rms, the implicit stepping's own error at a few steps a day.

**Cost (PERF.238).** The default resolution is 12 rows by 4 columns (12 for a body that turns fewer than `MANY_TURNS` times an orbit), 8 steps a solar day, 120 steps an orbit
   and at most 1500 steps: on 310 rocky planets and moons drawn from four star types and three zones the series and `temperature_coverage.coverage` took a mean of 66 ms,
   a median of 26 ms, 158 ms at the 90th percentile and 1.0 s at the most (40 ms with air, 72 ms for the airless ones, whose cold dry nights need many passes), against a mean
   of 1.1 s at `FINE`; the livable share is within 0.005 of the fine one on average and 0.033 at the 95th percentile (0.21 at the worst, a body on the edge of a band).

Limits. The orbit is circular (the eccentricity and `flux_variation` are not in the series yet); the sea is spread evenly over the longitudes and the latitudes (no land and sea
laid out); the ice follows the cell's mean temperature, not its season, and carries no latent heat or snow; the air shares the heat toward the world's mean at each step, one
share for all hours (a locked world and a slow spinner near resonance use different shares, `regional_temperature.is_locked`); no place is colder than `FLOOR_K`.
"""

import math
from dataclasses import dataclass, replace

import numpy as np

from planetgen.physics import climate_branches, constants, habitability, regional_temperature, surface_cycle

STEPS = 96
"""int: Equal steps of time (mean anomaly) in the year. Source tier T6, tweak level L0 (numerical)."""
WATER_HEAT_CAPACITY_J_M3_K = 4.1e6
"""float: The volumetric heat capacity of liquid water (about 4,186 J/kg/K at 1,000 to 1,025 kg/m3). Source tier T1, tweak level L0."""
OCEAN_MIXED_LAYER_M = 50.0
"""float: The depth of the ocean that follows the seasons, metres (the 50 m mixed layer of the simple climate models of [R4708] and [R4171, R4172],
recalled; the Earth's seasonal mixed layer is 20 to 100 m). Source tier T4, tweak level L0."""
LATITUDE_SAMPLES = 36
"""int: Latitudes (equal in sine of latitude, so equal in area) used to find the annual means' scale. Source tier T6, tweak level L0 (numerical)."""
MAX_SHARE = 0.999
"""float: The transport share is cut here so the damping stays finite. Source tier T6, tweak level L0 (numerical)."""
CONVERGENCE_K = 0.01
"""float: The year is repeated until no temperature moves more than this, kelvin. Source tier T6, tweak level L0 (numerical)."""
MAX_YEARS = 40
"""int: The most years the periodic solution is iterated for."""
FREEZING_K = 273.15
"""float: The freezing point of water, kelvin."""


class SeasonalSwingError(ValueError):
    """An input the model cannot use."""


@dataclass(frozen=True)
class Season:
    """
    The year at one latitude.

    Attributes:
        latitude_deg (float): The latitude (north positive).
        mean_k (float): The annual mean temperature there.
        max_k, min_k (float): The warmest and coldest temperature of the year.
        swing_k (float): `max_k - min_k`, peak to peak.
        series_k (tuple): The temperature at each of `STEPS` equal times of the year, starting at periapsis.
        warm_fraction (float): The share of the year above freezing.
        warmest_fraction (float): The time of the warmest step as a share of the year since periapsis.
    """

    latitude_deg: float
    mean_k: float
    max_k: float
    min_k: float
    swing_k: float
    series_k: tuple
    warm_fraction: float
    warmest_fraction: float


def _check(value, name, low=None, high=None):
    if value is None or not math.isfinite(value) or (low is not None and value < low) or (high is not None and value > high):
        raise SeasonalSwingError(f"{name} is out of range: {value}")


def orbit_year(eccentricity, perihelion_longitude_deg=0.0, steps=STEPS):
    """
    The orbit at `steps` equal times of the year from periapsis: `(sin L, (a / r)^2)`, `L` the true longitude from the equinox and `r` the distance (Kepler's equation
    by Newton's method; the eccentricity must be under 1).
    """
    _check(eccentricity, "eccentricity", 0.0, 0.999)
    mean_anomaly = (np.arange(steps) + 0.5) / steps * 2.0 * math.pi
    anomaly = mean_anomaly + eccentricity * np.sin(mean_anomaly)
    for _ in range(30):
        step = (anomaly - eccentricity * np.sin(anomaly) - mean_anomaly) / (1.0 - eccentricity * np.cos(anomaly))
        anomaly = anomaly - step
        if float(np.max(np.abs(step))) < 1e-13:
            break
    true = 2.0 * np.arctan2(math.sqrt(1.0 + eccentricity) * np.sin(anomaly / 2.0), math.sqrt(1.0 - eccentricity) * np.cos(anomaly / 2.0))
    distance = 1.0 - eccentricity * np.cos(anomaly)
    return np.sin(true + math.radians(perihelion_longitude_deg)), 1.0 / distance ** 2


def daily_insolation_year(latitudes_deg, obliquity_deg, eccentricity, perihelion_longitude_deg=0.0, steps=STEPS):
    """
    The day's mean insolation at each latitude through the year, in units of the flux at the orbit's semi-major axis: an array of shape `(latitudes, steps)`. The tilt is
    folded about 90 degrees (a retrograde spin is the same seasons the other way round).
    """
    _check(obliquity_deg, "obliquity", 0.0, 180.0)
    tilt = math.radians(obliquity_deg)
    sin_longitude, inverse_square = orbit_year(eccentricity, perihelion_longitude_deg, steps)
    declination = np.arcsin(math.sin(tilt) * sin_longitude)
    latitude = np.radians(np.asarray(latitudes_deg, dtype=float))[:, None]
    flux = regional_temperature._daily_insolation(latitude, declination[None, :])
    return flux * inverse_square[None, :]


def sphere_mean_insolation(eccentricity):
    """The annual mean over the whole sphere in units of the flux at the semi-major axis: `1 / (4 sqrt(1 - e^2))`."""
    return 1.0 / (4.0 * math.sqrt(1.0 - eccentricity ** 2))


def heat_capacity_j_m2_k(pressure_pa, gravity_ms2, ocean_fraction, period_years, ocean_heat_capacity_j_m3_k=WATER_HEAT_CAPACITY_J_M3_K,
                         mixed_layer_m=OCEAN_MIXED_LAYER_M):
    """
    The heat capacity per square metre of the layer the seasons reach: the sea's mixed layer on its share, the ground's skin on the rest, and the air column.
    """
    _check(ocean_fraction, "ocean fraction", 0.0, 1.0)
    _check(period_years, "period", 1e-6)
    period_s = period_years * 365.25 * 86400.0
    ground = regional_temperature.GROUND_INERTIA * math.sqrt(period_s / math.pi)
    sea = ocean_heat_capacity_j_m3_k * mixed_layer_m
    column = regional_temperature.AIR_CP_J_KG_K * pressure_pa / gravity_ms2 if pressure_pa and gravity_ms2 and pressure_pa > 0.0 else 0.0
    return ocean_fraction * sea + (1.0 - ocean_fraction) * ground + column


def _sample_latitudes():
    """Latitudes of equal area from pole to pole, degrees."""
    return np.degrees(np.arcsin(-1.0 + (np.arange(LATITUDE_SAMPLES) + 0.5) / LATITUDE_SAMPLES * 2.0))


def annual_mean_temperatures(latitudes_deg, obliquity_deg, eccentricity, mean_k, share, perihelion_longitude_deg=0.0):
    """
    The annual mean surface temperature at each latitude, kelvin: the world's mean spread by the fourth root of each place's annual-mean light, scaled so the area
    mean of the temperatures is `mean_k`, with the air carrying the share `share` of the contrast (the law of `regional_temperature.compute`).
    """
    grid = _sample_latitudes()
    ratio_grid = daily_insolation_year(grid, obliquity_deg, eccentricity, perihelion_longitude_deg).mean(axis=1) ** 0.25
    scale = float(np.mean(ratio_grid))
    ratio = daily_insolation_year(latitudes_deg, obliquity_deg, eccentricity, perihelion_longitude_deg).mean(axis=1) ** 0.25 / scale
    return mean_k + (ratio * mean_k - mean_k) * (1.0 - share)


def seasonal_swing(latitudes_deg, obliquity_deg, eccentricity, mean_k, stellar_flux_w_m2, bond_albedo, pressure_pa, gravity_ms2, ocean_fraction,
                   period_years, perihelion_longitude_deg=0.0, ocean_heat_capacity_j_m3_k=WATER_HEAT_CAPACITY_J_M3_K, mixed_layer_m=OCEAN_MIXED_LAYER_M):
    """
    The year at each of the latitudes (the module docstring).

    Args:
        latitudes_deg (sequence): Latitudes, degrees, north positive.
        obliquity_deg, eccentricity (float): The tilt and the orbit's eccentricity.
        mean_k (float): The world's annual mean surface temperature.
        stellar_flux_w_m2 (float): The star's flux at the semi-major axis.
        bond_albedo (float): The planet's Bond albedo, 0 to 1.
        pressure_pa, gravity_ms2 (float): The air (pressure 0 for none) and the surface gravity.
        ocean_fraction (float): The share of the surface that is liquid ocean, 0 to 1.
        period_years (float): The orbital period.
        perihelion_longitude_deg (float): The longitude of periapsis from the equinox (the phase of the seasons against the orbit).
        ocean_heat_capacity_j_m3_k, mixed_layer_m (float): The sea's volumetric heat capacity and mixed-layer depth.

    Returns:
        tuple: A `Season` for each latitude.
    """
    for value, name, low, high in ((mean_k, "mean temperature", 1.0, None), (stellar_flux_w_m2, "flux", 1e-9, None), (bond_albedo, "albedo", 0.0, 0.999),
                                   (gravity_ms2, "gravity", 1e-6, None)):
        _check(value, name, low, high)
    _check(eccentricity, "eccentricity", 0.0, 0.999)
    _check(obliquity_deg, "obliquity", 0.0, 180.0)
    latitudes = np.asarray(latitudes_deg, dtype=float)
    if latitudes.size == 0 or float(np.max(np.abs(latitudes))) > 90.0:
        raise SeasonalSwingError("latitudes must be given, each -90 to 90")
    column = regional_temperature.column_in_earths(pressure_pa, gravity_ms2)
    share = min(regional_temperature.transport_share(column), MAX_SHARE)
    capacity = heat_capacity_j_m2_k(pressure_pa, gravity_ms2, ocean_fraction, period_years, ocean_heat_capacity_j_m3_k, mixed_layer_m)
    sigma = constants.STEFAN_BOLTZMANN_CONSTANT
    absorbed_mean = (1.0 - bond_albedo) * stellar_flux_w_m2 * sphere_mean_insolation(eccentricity)
    emissivity = min(absorbed_mean / (sigma * mean_k ** 4), 1.0)
    local_mean = annual_mean_temperatures(latitudes, obliquity_deg, eccentricity, mean_k, share, perihelion_longitude_deg)
    insolation = daily_insolation_year(latitudes, obliquity_deg, eccentricity, perihelion_longitude_deg)
    forcing = (1.0 - bond_albedo) * stellar_flux_w_m2 * (insolation - insolation.mean(axis=1, keepdims=True))
    linear_slope = 4.0 * emissivity * sigma * local_mean ** 3
    damping = linear_slope * (1.0 + share / (1.0 - share))
    year_s = period_years * 365.25 * 86400.0
    omega = 2.0 * math.pi / year_s
    spectrum = np.fft.rfft(forcing, axis=1)
    response = spectrum / (damping[:, None] + 1j * capacity * omega * np.arange(spectrum.shape[1])[None, :])
    response[:, 0] = 0.0
    temperature = local_mean[:, None] + np.fft.irfft(response, n=STEPS, axis=1)
    step_s = year_s / STEPS
    transport_slope = linear_slope * share / (1.0 - share)
    state = temperature[:, -1].copy()
    for _ in range(MAX_YEARS):
        start = state.copy()
        series = np.empty_like(temperature)
        for k in range(STEPS):
            guess = state
            for _newton in range(5):
                residual = (capacity * (guess - state) / step_s - forcing[:, k] + emissivity * sigma * (guess ** 4 - local_mean ** 4)
                            + transport_slope * (guess - local_mean))
                slope = capacity / step_s + 4.0 * emissivity * sigma * guess ** 3 + transport_slope
                guess = guess - residual / slope
            state = guess
            series[:, k] = state
        temperature = series
        if float(np.max(np.abs(state - start))) < CONVERGENCE_K:
            break
    out = []
    for index, latitude in enumerate(latitudes):
        row = temperature[index]
        shifted = row - row.mean() + local_mean[index]
        out.append(Season(float(latitude), float(local_mean[index]), float(shifted.max()), float(shifted.min()), float(shifted.max() - shifted.min()),
                          tuple(float(v) for v in shifted), float(np.mean(shifted > FREEZING_K)), float((int(np.argmax(shifted)) + 0.5) / STEPS)))
    return tuple(out)


def peak_wet_bulb_c(season, relative_humidity):
    """The wettest-bulb hour of the year at the latitude, deg C: the warmest temperature of the year at a relative humidity (a fraction, 0.05 to 0.99), by Stull's fit
    (`habitability.wet_bulb_c`, valid -20 to 50 deg C; a hotter day is held at the fit's edge)."""
    hottest_c = min(season.max_k - FREEZING_K, 50.0)
    return habitability.wet_bulb_c(hottest_c, relative_humidity)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# The temperature of every place at every time (GEN.1263)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------

OCEAN_EDDY_DIFFUSIVITY_M2_S = 2.5e-4
"""float: How fast the sea's heat mixes down, m2/s: the depth a cycle of period `P` reaches is `sqrt(kappa P / pi)`, 50 m (`OCEAN_MIXED_LAYER_M`) for a year and about 2.6 m for a day, the
order of the observed mixed layers (a few metres by day, tens of metres by season). Source tier Confident educated guess from physics, tweak level L1 (set to give the 50 m of a
year; GEN.725)."""
FLOOR_K = regional_temperature.NIGHT_FLOOR_K
"""float: No place of the series is colder than this, kelvin (the heat from the ground, `regional_temperature.NIGHT_FLOOR_K`)."""
NONLINEAR_PASSES = 30
"""int: The most passes that put back the part of the emission the linear response leaves out. Source tier T6, tweak level L0 (numerical)."""
ANDERSON_DEPTH = 6
"""int: The passes the mixing remembers. Source tier T6, tweak level L0 (numerical)."""
NONLINEAR_TOLERANCE_K = 1.0
"""float: The passes stop when no temperature moves more than this, kelvin. Source tier T6, tweak level L0 (numerical)."""
EMISSION_PASSES = 5
"""int: The most times the outgoing emission is scaled to close the stored mean. Source tier T6, tweak level L0 (numerical)."""
EMISSION_TOLERANCE_K = 1.0
"""float: The mean of the series is taken as the stored one within this, kelvin, before the shift that makes it exact. Source tier T6, tweak level L0 (numerical)."""
MIN_EMISSIVITY = 0.02
"""float: The least outgoing emissivity the closure takes (a greenhouse of 98 percent). Source tier T6, tweak level L0 (numerical)."""
TEMPERATURE_MAX_STEPS = 1500
"""int: The most time steps of a temperature series (PERF.238: the median body costs about 26 ms; `FINE` is the resolution the tests compare it with). Source tier T6, tweak level L0 (numerical)."""
SAMPLES_PER_DAY = 8
"""int: The least time steps of a solar day when `TEMPERATURE_MAX_STEPS` allows them. Source tier T6, tweak level L0 (numerical)."""
STEPS_PER_YEAR = 120
"""int: The least time steps of an orbit (the declination moves slowly; `surface_cycle.STEPS_PER_YEAR` is 3000). Source tier T6, tweak level L0 (numerical)."""
FEW_COLUMNS = 4
"""int: The columns of longitude of a body that turns many times in an orbit (`MANY_TURNS`): every column then sees every hour of the day in each season, so its statistics are the same
and more columns only repeat them. Source tier T6, tweak level L0 (numerical)."""
MANY_TURNS = 24
"""float: The solar days in an orbit from which `FEW_COLUMNS` are used; a body with fewer (a slow spinner, a near-resonant or a locked body, where longitude is the place of the cell
in the day) has `surface_cycle.CELL_COLUMNS`. Source tier T6, tweak level L0 (numerical)."""
FINE = {"columns": surface_cycle.CELL_COLUMNS, "samples_per_day": 24, "steps_per_year": surface_cycle.STEPS_PER_YEAR, "max_steps": 6000}
"""dict: The keyword arguments of `temperature_series` for the fine resolution, which the default is checked against (it costs three to ten times as much)."""
CEILING_MARGIN = 1.05
"""float: No place of the series is hotter than the fourth root of the hottest light a cell gets over its emission, times this (the passes of a hot airless body overshoot and overflow
without it). Source tier T6, tweak level L0 (numerical)."""
ICE_OUTER_PASSES = 2
"""int: The times the ice cover is laid again from the cell means the temperatures give. Source tier T6, tweak level L0 (numerical)."""


@dataclass(frozen=True)
class TemperatureSeries:
    """
    The temperature of every place of a body through its year and its turns.

    Attributes:
        grid (surface_cycle.Grid): The cells (rows of equal area by columns of equal longitude).
        time_days (ndarray): The time of each step, days from an equinox.
        step_days, span_days (float): The step and the span (a whole number of orbits that holds a whole number of solar days).
        temperature_k (ndarray): `[row, column, step]`, kelvin, shifted so the area-and-time mean is the stored mean exactly.
        mean_k (ndarray): `[row, column]`, the time mean of each cell.
        albedo (ndarray): `[row, column]`, the Bond albedo of each cell (the ice law laid over the cells).
        ice_cover (ndarray): `[row, column]`, the share of each cell under ice, 0 to 1 (zeros when no ice law was given).
        share (float): The share of the contrast the air carries (`regional_temperature.transport_share`).
        solar_days, orbital_days (float): The solar day as the series holds it, and the orbit.
        shift_k (float): The shift that was added to make the mean the stored one exactly (what the closure of the emission left, a fraction of a kelvin unless its limits were reached).
        emissivity (float): The outgoing emissivity the mean was closed with.
        passes (int): The passes of the nonlinear emission that were made.
        converged (bool): Whether the last pass moved no temperature more than `NONLINEAR_TOLERANCE_K`.
    """

    grid: object
    time_days: np.ndarray
    step_days: float
    span_days: float
    temperature_k: np.ndarray
    mean_k: np.ndarray
    albedo: np.ndarray
    ice_cover: np.ndarray
    share: float
    solar_days: float
    orbital_days: float
    shift_k: float
    emissivity: float
    passes: int
    converged: bool

    def area_weights(self):
        """The weight of each row of cells, `[rows, 1]`, adding to 1 over the rows (one cell holds `weight / columns`)."""
        return self.grid.row_weight[:, None]

    def area_time_mean_k(self):
        """The mean over the cells (by area) and the time."""
        return float(np.sum(self.temperature_k.mean(axis=(1, 2)) * self.grid.row_weight))


def columns_for(rotation_hours, orbital_hours, tilt_deg):
    """The columns of longitude a body needs: `FEW_COLUMNS` when it turns `MANY_TURNS` times or more in an orbit, else `surface_cycle.CELL_COLUMNS`."""
    solar_days = surface_cycle.solar_day_days(rotation_hours, orbital_hours, tilt_deg)
    if not math.isfinite(solar_days) or solar_days <= 0.0:
        return surface_cycle.CELL_COLUMNS
    return FEW_COLUMNS if orbital_hours / 24.0 / solar_days >= MANY_TURNS else surface_cycle.CELL_COLUMNS


def _ocean_depth_m(period_s, mixed_layer_m, diffusivity):
    return np.minimum(mixed_layer_m, np.sqrt(diffusivity * period_s / math.pi))


def capacity_by_harmonic(span_s, harmonics, pressure_pa, gravity_ms2, ocean_fraction, ocean_heat_capacity_j_m3_k=WATER_HEAT_CAPACITY_J_M3_K,
                         mixed_layer_m=OCEAN_MIXED_LAYER_M, diffusivity=OCEAN_EDDY_DIFFUSIVITY_M2_S):
    """
    The heat capacity per square metre of the layer a cycle of each harmonic reaches, J/(m2 K), for the harmonics `1 ... harmonics - 1` of a span (index 0 is not used): the air
    column, the sea to the depth the cycle mixes (at most `mixed_layer_m`) on its share, and the ground's skin (`regional_temperature.GROUND_INERTIA` times `sqrt(period / pi)`) on the rest.
    A day's cycle reaches a thin layer and a year's a thick one, which one capacity for all cycles cannot say.
    """
    _check(ocean_fraction, "ocean fraction", 0.0, 1.0)
    period = span_s / np.maximum(np.arange(harmonics), 1)
    ground = regional_temperature.GROUND_INERTIA * np.sqrt(period / math.pi)
    sea = ocean_heat_capacity_j_m3_k * _ocean_depth_m(period, mixed_layer_m, diffusivity)
    column = regional_temperature.AIR_CP_J_KG_K * pressure_pa / gravity_ms2 if pressure_pa and gravity_ms2 and pressure_pa > 0.0 else 0.0
    return ocean_fraction * sea + (1.0 - ocean_fraction) * ground + column


def _lay_ice(cell_mean_k, weights, law, ice_fraction):
    """The share of each cell under ice, 0 to 1: the GEN.972 law's shape (how wide the ice line is) with its window moved until the area mean of the cover is the stored share of the
    ice (`ground_ice_fraction`); the law as it is when no share is stored."""
    if ice_fraction is None:
        return np.vectorize(law.ice_fraction)(cell_mean_k)
    if ice_fraction <= 0.0:
        return np.zeros_like(cell_mean_k)
    if ice_fraction >= 1.0:
        return np.ones_like(cell_mean_k)

    def cover(offset):
        moved = replace(law, free_above_k=law.free_above_k + offset, covered_below_k=law.covered_below_k + offset)
        return np.vectorize(moved.ice_fraction)(cell_mean_k)

    low, high = -400.0, 400.0
    for _ in range(50):
        middle = 0.5 * (low + high)
        if float(np.sum(cover(middle) * weights)) < ice_fraction:
            low = middle
        else:
            high = middle
    return cover(0.5 * (low + high))


def _anderson(step, start):
    """
    The fixed point of `step` by Anderson mixing (depth `ANDERSON_DEPTH`): the series where one more pass moves no temperature more than `NONLINEAR_TOLERANCE_K`. Returns the series, the
    passes made and whether it settled. The passes of a plain iteration converge slowly where the emission is steep and not at all where a night is cold and dry.
    """
    state = start
    updated = step(state)
    shape = updated.shape
    size = updated.size
    delta_f = np.empty((ANDERSON_DEPTH, size))
    delta_g = np.empty((ANDERSON_DEPTH, size))
    count = 0
    previous_f = previous_g = None
    passes = 1
    for passes in range(1, NONLINEAR_PASSES + 1):
        change = (updated - state).ravel()
        if float(np.max(np.abs(change))) < NONLINEAR_TOLERANCE_K:
            return updated, passes, True
        current_g = updated.ravel()
        if previous_f is not None:
            slot = (count) % ANDERSON_DEPTH
            np.subtract(change, previous_f, out=delta_f[slot])
            np.subtract(current_g, previous_g, out=delta_g[slot])
            count += 1
        previous_f, previous_g = change, current_g
        used = min(count, ANDERSON_DEPTH)
        if used:
            matrix = delta_f[:used]
            gram = matrix @ matrix.T
            gamma = np.linalg.solve(gram + (1.0e-10 * np.trace(gram) + 1.0e-30) * np.eye(used), matrix @ change)
            mixed = current_g - gamma @ delta_g[:used]
        else:
            mixed = current_g
        state = np.maximum(mixed.reshape(shape), FLOOR_K)
        updated = step(state)
    return updated, passes, False


def _mean_balance(absorbed_w_m2, emissivity_sigma, transport, world_mean_k, weights):
    """Each cell's steady temperature: `absorbed = eps sigma T^4 + h (T - Tm)`, by Newton from the fourth root of the absorbed light; `h` the transport (W/m2/K)."""
    state = np.maximum((absorbed_w_m2 / emissivity_sigma) ** 0.25, FLOOR_K)
    state = 0.5 * state + 0.5 * world_mean_k
    for _ in range(40):
        residual = emissivity_sigma * state ** 4 + transport * (state - world_mean_k) - absorbed_w_m2
        slope = 4.0 * emissivity_sigma * state ** 3 + transport
        step = residual / slope
        state = np.maximum(state - step, FLOOR_K)
        if float(np.max(np.abs(step))) < 1.0e-6:
            break
    return state


def temperature_series(rotation_hours, orbital_hours, tilt_deg, mean_k, stellar_flux_w_m2, bond_albedo, pressure_pa, gravity_ms2, ocean_fraction,
                       ice_fraction=None, ice_law=None, transmission=None, rows=surface_cycle.CELL_ROWS, columns=None, samples_per_day=SAMPLES_PER_DAY, steps_per_year=STEPS_PER_YEAR, max_steps=TEMPERATURE_MAX_STEPS,
                       ocean_heat_capacity_j_m3_k=WATER_HEAT_CAPACITY_J_M3_K, mixed_layer_m=OCEAN_MIXED_LAYER_M, diffusivity=OCEAN_EDDY_DIFFUSIVITY_M2_S):
    """
    The temperature of every cell of the grid at every step of the year (the module section "The temperature of every place"): `TemperatureSeries`.

    Args:
        rotation_hours, orbital_hours (float): The stored rotation period and orbital period.
        tilt_deg (float): The axial tilt (past 90 the spin is backwards).
        mean_k (float): The stored mean surface temperature, which the area-and-time mean of the series is exactly.
        stellar_flux_w_m2, bond_albedo (float): The flux at the planet and its mean Bond albedo (these give the outgoing emissivity as `seasonal_swing` does).
        pressure_pa, gravity_ms2 (float): The air (0 for none) and the gravity.
        ocean_fraction (float): The share of the surface that is liquid ocean (the stored `ground_ocean_fraction`).
        ice_fraction (float or None): The share of the surface under ice (the stored `ground_ice_fraction`).
        ice_law (climate_branches.IceAlbedoLaw or None): The ice-albedo law of GEN.972; `None` keeps one albedo for every cell.
        transmission (callable or None): `surface_cycle.series`'s ring shade.
        samples_per_day, max_steps (int): The time resolution (`surface_cycle.time_step_days`).
    """
    for value, name, low, high in ((mean_k, "mean temperature", 1.0, None), (stellar_flux_w_m2, "flux", 1e-9, None), (bond_albedo, "albedo", 0.0, 0.999),
                                   (gravity_ms2, "gravity", 1e-6, None), (rotation_hours, "rotation period", 1e-6, None), (orbital_hours, "orbital period", 1e-6, None)):
        _check(value, name, low, high)
    _check(tilt_deg, "tilt", 0.0, 360.0)
    if columns is None:
        columns = columns_for(rotation_hours, orbital_hours, tilt_deg)
    light = surface_cycle.series(rotation_hours, orbital_hours, tilt_deg, None, transmission, rows=rows, columns=columns, samples_per_day=samples_per_day, max_steps=max_steps,
                                 periodic=True, steps_per_year=steps_per_year)
    grid = light.grid
    rows, columns, steps = light.cosine.shape
    weights = np.repeat(grid.row_weight / columns, columns).reshape(rows, columns)
    span_days = steps * light.step_days
    span_s = span_days * 86400.0
    sunlight = np.maximum(light.cosine, 0.0)
    sunlight_mean = sunlight.mean(axis=2)

    sigma = constants.STEFAN_BOLTZMANN_CONSTANT
    absorbed_world = (1.0 - bond_albedo) * stellar_flux_w_m2 / 4.0
    emissivity = min(absorbed_world / (sigma * mean_k ** 4), 1.0)
    eps_sigma = emissivity * sigma
    column = regional_temperature.column_in_earths(pressure_pa, gravity_ms2)
    locked = regional_temperature.is_locked(rotation_hours, orbital_hours / (365.25 * 24.0))
    share = min(regional_temperature.transport_share(column, locked), MAX_SHARE)
    transport = 4.0 * eps_sigma * mean_k ** 3 * share / (1.0 - share)

    albedo = np.full((rows, columns), bond_albedo)
    cover = np.zeros((rows, columns))
    passes_outer = ICE_OUTER_PASSES if ice_law is not None else 1
    for _ in range(passes_outer):
        scale = eps_sigma * mean_k ** 4 / float(np.sum((1.0 - albedo) * sunlight_mean * weights))
        cell_mean = _mean_balance(scale * (1.0 - albedo) * sunlight_mean, eps_sigma, transport, mean_k, weights)
        cell_mean += mean_k - float(np.sum(cell_mean * weights))
        if ice_law is not None:
            cover = _lay_ice(cell_mean, weights, ice_law, ice_fraction)
            albedo = ice_law.ice_free + (ice_law.ice_covered - ice_law.ice_free) * cover
    scale = eps_sigma * mean_k ** 4 / float(np.sum((1.0 - albedo) * sunlight_mean * weights))
    forcing = scale * (1.0 - albedo)[:, :, None] * sunlight

    harmonics = steps // 2 + 1
    omega = 2.0 * math.pi * np.arange(harmonics) / span_s
    capacity = capacity_by_harmonic(span_s, harmonics, pressure_pa, gravity_ms2, ocean_fraction, ocean_heat_capacity_j_m3_k, mixed_layer_m, diffusivity)
    spectrum_forcing = np.fft.rfft(forcing, axis=2)
    capacity_term = 1j * (omega * capacity)[None, None, :]
    weight = weights[:, :, None]
    absorbed_mean = scale * (1.0 - albedo) * sunlight_mean

    def solve(emission, start):
        """The periodic temperatures for an outgoing emission `emission = eps sigma`, from `start` (or from the cells' steady temperatures)."""
        carried = 4.0 * emission * mean_k ** 3 * share / (1.0 - share)
        steady = _mean_balance(absorbed_mean, emission, carried, mean_k, weights)
        steady += mean_k - float(np.sum(steady * weights))
        state = np.maximum(steady[:, :, None] + np.zeros_like(forcing), FLOOR_K) if start is None else start
        # The linear operator that is inverted at every pass: the emission linearised about the cell's own steady temperature (a slope that is the damping of a cycle), the air's
        # transport towards the world's mean, and the capacity of the layer each cycle reaches. It is a preconditioner: the passes put right whatever it leaves out, so a poor
        # slope costs passes, not accuracy.
        inverse = 1.0 / ((4.0 * emission * steady ** 3 + carried)[:, :, None] + capacity_term)
        coupling_denominator = 1.0 - carried * np.sum(weight * inverse, axis=(0, 1))
        # No place can be hotter than the one that sends the hottest light it ever gets straight back out: a hot airless body with a short day (GEN.1265) overshoots otherwise,
        # the fourth power overflows and the whole series is not a number.
        ceiling = CEILING_MARGIN * (float(np.max(forcing)) / emission) ** 0.25

        def step(temperature):
            spectrum = np.fft.rfft(temperature, axis=2)
            world = np.sum(weight * spectrum, axis=(0, 1))
            residual = spectrum_forcing - np.fft.rfft(emission * temperature ** 4, axis=2) - carried * (spectrum - world[None, None, :]) - capacity_term * spectrum
            coupling = np.sum(weight * residual * inverse, axis=(0, 1)) / coupling_denominator
            change = np.fft.irfft((residual + carried * coupling[None, None, :]) * inverse, n=steps, axis=2)
            return np.clip(temperature + change, FLOOR_K, ceiling)

        return _anderson(step, state)

    # The outgoing emission is closed by the stored mean: the part of the world that runs hot or cold makes the mean of the fourth powers more than the fourth power of the mean, so the
    # emission is scaled until the mean of the series is the stored one, and the shift that is left is rounding.
    emission = eps_sigma
    temperature, passes, converged = solve(emission, None)
    for _ in range(EMISSION_PASSES):
        mean = float(np.sum(temperature.mean(axis=2) * weights))
        if abs(mean - mean_k) < EMISSION_TOLERANCE_K:
            break
        emission = min(max(emission * (mean / mean_k) ** 4, MIN_EMISSIVITY * sigma), sigma)
        temperature, more, converged = solve(emission, temperature)
        passes += more
    shift = mean_k - float(np.sum(temperature.mean(axis=2) * weights))
    temperature = temperature + shift
    return TemperatureSeries(grid, light.time_days, light.step_days, span_days, temperature, temperature.mean(axis=2), albedo, cover, share, light.solar_days,
                             light.orbital_days, float(shift), emission / sigma, passes, converged)
