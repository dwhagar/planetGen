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

The simulated survey stands for the generator's planets only in its statistics; it draws no planet the generator keeps.
"""

import math

from planetgen.physics import bulk_structure, kepler_recovery
from planetgen.util import draw

PLANETS = 4000
"""int: Planets per survey."""

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
KEPMAG = 14.0
AGE_GY = 5.0
XUV_ENERGY_J_M2 = 3.1e16
"""float: The XUV energy per square metre a planet at one Earth insolation receives over its life (the row's constant)."""

BIN_EDGES_D = (2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
RADIUS_RANGE_EARTH = (1.2, 3.2)
HISTOGRAM_BINS = 20
MIN_PLANETS_PER_BIN = 40

G = 6.674e-11
M_EARTH = 5.972e24
R_EARTH = 6.371e6
M_SUN = 1.989e30
R_SUN = 6.957e8
MAX_RADIUS_RATIO = 0.07
"""float: The largest radius ratio the recovery model covers (kepler_recovery validity); a larger planet is not counted."""

LOSS_PER_ENVELOPE = 10.0
"""float: The efficiency `eta` along the degenerate ridge of the tuning, `eta = LOSS_PER_ENVELOPE * f_env0` (see the note, section 2)."""
AU = 1.496e11


def planet(f_env0, eta):
    """One planet of the survey, drawn: `(period_d, radius_earth, weight)`."""
    host = draw.uniform(*HOST_MASS_RANGE_SOLAR)
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
        return period, radius, 0.0
    weight = kepler_recovery.detection_probability({
        "period_d": period, "radius_earth": radius, "host_radius_solar": host_radius, "host_mass_solar": host,
        "kepmag": KEPMAG, "cdpp_z": 0.0})
    return period, radius, weight


def valley_of(planets):
    """
    The valley of a list of `(period_d, radius_earth, weight)`: `(slope, radius at 10 days)`, or None when fewer than three period
    bins hold `MIN_PLANETS_PER_BIN` planets.
    """
    low, high = math.log10(RADIUS_RANGE_EARTH[0]), math.log10(RADIUS_RANGE_EARTH[1])
    step = (high - low) / HISTOGRAM_BINS
    points = []
    for lo, hi in zip(BIN_EDGES_D, BIN_EDGES_D[1:]):
        inside = [(r, w) for p, r, w in planets if lo <= p < hi]
        if len(inside) < MIN_PLANETS_PER_BIN:
            continue
        counts = [0.0] * HISTOGRAM_BINS
        for r, w in inside:
            k = int((math.log10(r) - low) / step)
            if 0 <= k < HISTOGRAM_BINS:
                counts[k] += w
        smooth = [counts[max(i - 1, 0)] + 2.0 * counts[i] + counts[min(i + 1, HISTOGRAM_BINS - 1)] for i in range(HISTOGRAM_BINS)]
        edge = 3
        k = min(range(edge, HISTOGRAM_BINS - edge), key=lambda i: smooth[i])
        points.append((math.log10(math.sqrt(lo * hi)), low + (k + 0.5) * step))
    if len(points) < 3:
        return None
    n = len(points)
    mean_x = sum(x for x, _ in points) / n
    mean_y = sum(y for _, y in points) / n
    slope = sum((x - mean_x) * (y - mean_y) for x, y in points) / sum((x - mean_x) ** 2 for x, _ in points)
    return slope, 10.0 ** (mean_y + slope * (1.0 - mean_x))


def survey(f_env0, eta, planets=PLANETS):
    """A survey drawn on the bound stream: the list of planets."""
    return [planet(f_env0, eta) for _ in range(planets)]


def valley_sampler(f_env0, eta, field, planets=PLANETS):
    """A sampler for `tuning_harness.simulate`: one survey per call, its `field` ("slope" or "radius_10d") or None when too sparse."""
    def sample(_i):
        found = valley_of(survey(f_env0, eta, planets))
        if found is None:
            return None
        return found[0] if field == "slope" else found[1]
    return sample
