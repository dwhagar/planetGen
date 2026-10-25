# planetgen/physics/radius_valley.py

"""
The radius valley of a simulated Kepler survey (GEN.509, docs/design/composition/radius-valley-tuning.md)
=========================================================================================================

Rows e1 and e2 of the atmosphere recipes say whether a close-in planet of 1 to 30 Earth masses keeps its hydrogen-helium envelope
(`weight = 1 - f_lost/f_env0`, with `f_lost` the energy-limited XUV loss). Their two knobs, the starting envelope fraction `f_env0`
and the loss efficiency `eta`, cannot be read from any measurement of a single planet: they are degenerate with the core-mass
distribution. They are tuned to the one thing the population shows, the valley between bare cores and envelope planets in the
period-radius plane of the planets Kepler found.

A survey is `PLANETS` planets drawn by the one random wrapper (`draw`):

1. a host of 0.8 to 1.3 solar masses (the California-Kepler Survey range) with `L = M^3.5` and `R = M^0.9` (Kepler magnitude 14);
2. a period log-uniform from 2 to 100 days;
3. a core mass log-normal around `CORE_MASS_MEDIAN_EARTH` with `CORE_MASS_SIGMA` (the peak of Rogers and Owen 2021) and a core
   radius `M^0.25` (Lopez and Fortney 2014 Eq. 1);
4. a starting envelope `f_env0 * exp(ENVELOPE_SIGMA z)` of the core mass;
5. the loss of the recipe row: `f_lost = LOSS_PER_ENVELOPE f_env0 (v_esc/15 km/s)^-2 pi R_core^3 E_xuv F / (G M_core^2)`, with `E_xuv = 3.1e16 F` joules
   per square metre and `F` the insolation (Owen and Wu 2017's variable efficiency, decision D2);
6. what is left is an envelope of `f_env0 - f_lost` (bare below `BARE_BELOW`) whose radius is the Lopez and Fortney 2014 Eq. 3 law
   at the planet's flux and 5 Gyr (`bulk_structure.envelope_radius_earth`);
7. the Kepler transit-plus-recovery chance of the planet (`kepler_recovery.detection_probability`, GEN.504) weights it.

The valley of a period bin is the minimum of the smoothed histogram of log10 R between 1.2 and 3.2 Earth radii, weighted by that
chance. A line through the valleys of the period bins gives the slope d log10 R / d log10 P and the radius at 10 days. They are
compared with the published -0.11 +- 0.02 (Petigura et al. 2022) and 1.84 +- 0.03 Earth radii at 10 days (the target table).

Stellar-mass slopes (GEN.722): `survey(..., host_range=WIDE_HOST_MASS_RANGE_SOLAR)` draws hosts of 0.5 to 1.4 solar masses and
`mass_slopes_of` fits `log R = c + alpha log F + beta log M*` to the valley of cells of flux and host mass, the form Berger et al. 2023
fit (`alpha` 0.069, `beta` -0.046 +- 0.12). The loss as written gives core-powered mass loss's pair (alpha 0.07, beta 0.00); the data
cannot tell it from photoevaporation's (beta -0.17), so the mechanism stays an open stand-in (the tuning doc, section 4).

The simulated survey stands for the generator's planets only in its statistics; it draws no planet the generator keeps.
"""

import math
from typing import NamedTuple

from planetgen.physics import bulk_structure, kepler_recovery
from planetgen.physics import constants as physical_constants
from planetgen.util import draw

PLANETS = 4000
"""int: Planets per survey."""
WIDE_SURVEY_PLANETS = 12000
"""int: Planets per survey of the stellar-mass slopes (eighteen cells of them, so three times as many)."""

CORE_MASS_MEDIAN_EARTH = 3.72
"""float: Median core mass, Earth masses (the peak of Rogers and Owen 2021, arXiv:2007.11006, 3.72 +0.45 -0.33)."""

CORE_MASS_SIGMA = 0.4
"""float: Log-normal width of the core mass (a stand-in: no published width is read)."""

ENVELOPE_SIGMA = 0.6
"""float: Log-normal width of the starting envelope fraction (a stand-in)."""

MAX_ENVELOPE = 0.5
"""float: An envelope fraction is held at this (Lopez and Fortney's law covers a few percent to tens of percent; beyond it the planet is a giant)."""

BARE_BELOW = 0.002
"""float: An envelope under this fraction of the core mass is stripped."""

PERIOD_RANGE_D = (2.0, 100.0)
HOST_MASS_RANGE_SOLAR = (0.8, 1.3)
WIDE_HOST_MASS_RANGE_SOLAR = (0.5, 1.4)
"""tuple: The hosts of the survey that measures the valley against stellar mass (GEN.722): the 0.5 to 1.4 solar masses of Petigura et al.
2022's sample, where the stellar-mass slopes are measured; `HOST_MASS_RANGE_SOLAR` stays the range the period slope and the radius at 10 days were
tuned on."""
MASS_BINS_SOLAR = ((0.5, 0.7), (0.7, 1.0), (1.0, 1.4))
FLUX_BINS_EARTH = ((1.0, 3.0), (3.0, 10.0), (10.0, 30.0), (30.0, 100.0), (100.0, 300.0), (300.0, 1000.0))
"""tuple: The stellar-mass bins and the insolation bins (Earth units) in which `slopes_of` finds a valley each, as Berger et al. 2023 bin the
Kepler planets."""
KEPMAG = 14.0
AGE_GY = 5.0
XUV_ENERGY_J_M2 = 3.1e16
"""float: The XUV energy per square metre a planet at one Earth insolation receives over its life (the row's constant)."""

BIN_EDGES_D = (2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
RADIUS_RANGE_EARTH = (1.2, 3.2)
HISTOGRAM_BINS = 20
MIN_PLANETS_PER_BIN = 40

G = physical_constants.G
M_EARTH = physical_constants.EARTH_MASS_TO_KG
R_EARTH = physical_constants.EARTH_RADIUS_KM * 1000.0
M_SUN = physical_constants.SOLAR_MASS_TO_KG
R_SUN = physical_constants.SOLAR_RADIUS_M
MAX_RADIUS_RATIO = 0.07
"""float: The largest radius ratio the recovery model covers (kepler_recovery validity); a larger planet is not counted."""

LOSS_PER_ENVELOPE = 10.0
"""float: The efficiency `eta` along the degenerate ridge of the tuning, `eta = LOSS_PER_ENVELOPE * f_env0` (see the note, section 2)."""
AU = 1.496e11


class Candidate(NamedTuple):
    """One planet of a survey: its period (d), host mass (solar), insolation (Earth), radius (Earth radii) and the chance Kepler finds it."""

    period_d: float
    host_mass_solar: float
    insolation_earth: float
    radius_earth: float
    weight: float


def planet(f_env0, eta, host_range=HOST_MASS_RANGE_SOLAR):
    """One planet of the survey, drawn: a `Candidate`."""
    host = draw.uniform(*host_range)
    period = math.exp(draw.uniform(math.log(PERIOD_RANGE_D[0]), math.log(PERIOD_RANGE_D[1])))
    core_mass = CORE_MASS_MEDIAN_EARTH * math.exp(CORE_MASS_SIGMA * draw.gauss())
    start = f_env0 * math.exp(ENVELOPE_SIGMA * draw.gauss())
    a_au = (G * host * M_SUN * (period * 86400.0) ** 2 / (4.0 * math.pi ** 2)) ** (1.0 / 3.0) / AU
    insolation = host ** 3.5 / a_au ** 2
    core_radius = bulk_structure.core_radius_earth(core_mass)
    v_esc = math.sqrt(2.0 * G * core_mass * M_EARTH / (core_radius * R_EARTH)) / 1000.0
    efficiency = eta * (v_esc / 15.0) ** -2
    lost = efficiency * math.pi * (core_radius * R_EARTH) ** 3 * XUV_ENERGY_J_M2 * insolation / (G * (core_mass * M_EARTH) ** 2)
    left = min(start - lost, MAX_ENVELOPE)
    if left < BARE_BELOW:
        radius = core_radius
    else:
        radius = bulk_structure.envelope_radius_earth(core_mass / (1.0 - left), left, insolation, AGE_GY)
    host_radius = host ** 0.9
    if radius * R_EARTH / (host_radius * R_SUN) > MAX_RADIUS_RATIO:
        return Candidate(period, host, insolation, radius, 0.0)
    weight = kepler_recovery.detection_probability({
        "period_d": period, "radius_earth": radius, "host_radius_solar": host_radius, "host_mass_solar": host,
        "kepmag": KEPMAG, "cdpp_z": 0.0})
    return Candidate(period, host, insolation, radius, weight)


def _valley_log_radius(radii_weights):
    """The log10 radius of the minimum of the smoothed weighted histogram of `(radius_earth, weight)`, or None with fewer than
    `MIN_PLANETS_PER_BIN` planets."""
    if len(radii_weights) < MIN_PLANETS_PER_BIN:
        return None
    low, high = math.log10(RADIUS_RANGE_EARTH[0]), math.log10(RADIUS_RANGE_EARTH[1])
    step = (high - low) / HISTOGRAM_BINS
    counts = [0.0] * HISTOGRAM_BINS
    for r, w in radii_weights:
        k = int((math.log10(r) - low) / step)
        if 0 <= k < HISTOGRAM_BINS:
            counts[k] += w
    smooth = [counts[max(i - 1, 0)] + 2.0 * counts[i] + counts[min(i + 1, HISTOGRAM_BINS - 1)] for i in range(HISTOGRAM_BINS)]
    edge = 3
    k = min(range(edge, HISTOGRAM_BINS - edge), key=lambda i: smooth[i])
    return low + (k + 0.5) * step


def valley_of(planets):
    """
    The valley of a list of `(period_d, radius_earth, weight)`: `(slope, radius at 10 days)`, or None when fewer than three period
    bins hold `MIN_PLANETS_PER_BIN` planets.
    """
    points = []
    for lo, hi in zip(BIN_EDGES_D, BIN_EDGES_D[1:]):
        valley = _valley_log_radius([(r, w) for p, r, w in planets if lo <= p < hi])
        if valley is not None:
            points.append((math.log10(math.sqrt(lo * hi)), valley))
    if len(points) < 3:
        return None
    n = len(points)
    mean_x = sum(x for x, _ in points) / n
    mean_y = sum(y for _, y in points) / n
    slope = sum((x - mean_x) * (y - mean_y) for x, y in points) / sum((x - mean_x) ** 2 for x, _ in points)
    return slope, 10.0 ** (mean_y + slope * (1.0 - mean_x))


def survey(f_env0, eta, planets=PLANETS, host_range=HOST_MASS_RANGE_SOLAR):
    """A survey drawn on the bound stream: the list of `Candidate`s."""
    return [planet(f_env0, eta, host_range) for _ in range(planets)]


def triples(candidates):
    """The `(period_d, radius_earth, weight)` of each candidate, what `valley_of` reads."""
    return [(c.period_d, c.radius_earth, c.weight) for c in candidates]


def _fit_plane(rows):
    """Least squares `y = c + a x1 + b x2` over `(x1, x2, y)`: `(a, b)`."""
    import numpy as np
    matrix = np.array([[1.0, x1, x2] for x1, x2, _ in rows])
    solution = np.linalg.lstsq(matrix, np.array([y for _, _, y in rows]), rcond=None)[0]
    return float(solution[1]), float(solution[2])


def mass_slopes_of(candidates):
    """
    The valley against insolation and stellar mass, as Berger et al. 2023 fit it: a valley in each cell of `MASS_BINS_SOLAR` by
    `FLUX_BINS_EARTH` (the minimum of the smoothed weighted histogram, as `valley_of`), then `log10 R_valley = c + alpha log10 F + beta log10 M`
    through the cells. Returns `(alpha, beta)`, or None when fewer than four cells hold `MIN_PLANETS_PER_BIN` planets.
    """
    rows = []
    for mass_lo, mass_hi in MASS_BINS_SOLAR:
        for flux_lo, flux_hi in FLUX_BINS_EARTH:
            valley = _valley_log_radius([(c.radius_earth, c.weight) for c in candidates
                                         if mass_lo <= c.host_mass_solar < mass_hi and flux_lo <= c.insolation_earth < flux_hi])
            if valley is not None:
                rows.append((math.log10(math.sqrt(flux_lo * flux_hi)), math.log10(math.sqrt(mass_lo * mass_hi)), valley))
    if len(rows) < 4:
        return None
    return _fit_plane(rows)


def valley_sampler(f_env0, eta, field, planets=PLANETS):
    """A sampler for `tuning_harness.simulate`: one survey per call, its `field` ("slope" or "radius_10d") or None when too sparse."""
    def sample(_i):
        found = valley_of(triples(survey(f_env0, eta, planets)))
        if found is None:
            return None
        return found[0] if field == "slope" else found[1]
    return sample


def mass_slope_sampler(f_env0, eta, field, planets=WIDE_SURVEY_PLANETS):
    """A sampler for `tuning_harness.simulate`: one survey of hosts of 0.5 to 1.4 solar masses per call, its `field` ("flux_slope", the
    valley's slope against insolation, or "mass_slope", against stellar mass at a given insolation) or None when too sparse."""
    def sample(_i):
        found = mass_slopes_of(survey(f_env0, eta, planets, WIDE_HOST_MASS_RANGE_SOLAR))
        if found is None:
            return None
        return found[0] if field == "flux_slope" else found[1]
    return sample
