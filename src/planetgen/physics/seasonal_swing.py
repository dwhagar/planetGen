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
Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.81.
"""

import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import constants, habitability, regional_temperature

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
