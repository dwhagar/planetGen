# planetgen/physics/orbit_elements.py

"""
Eccentricity and inclination of the orbits of planets and moons (GEN.727)
=========================================================================

A planet or moon used to go round on a circle with a plane tilt drawn flat between two bounds. Real orbits are ellipses, and the
ellipse carries the light a planet gets (the mean flux rises as `1 / sqrt(1 - e^2)`), the seasons (the flux swings between periapsis
and apoapsis), the room the neighbours need (apoapsis against periapsis) and a body's tide (`tidal_heating`, GEN.737).
This module draws the elements and holds the few functions of them the rest of the generator uses.

**Planets.** The eccentricity of a planet of a system with two or more planets follows the positive half of a zero-mean Gaussian of
dispersion 0.083 ([R1307], the Kepler systems with several transiting planets: 0.083 +0.015 -0.020). A
system with one planet is a mixture of that law and a broader one: 0.76 of single planets draw from a Rayleigh law peaking at 0.26,
the rest from the multi-planet law (the same paper's mixture model: "a fraction of 0.76 +0.21 -0.12 ... a Rayleigh distribution that
peaks at 0.26 +0.04 -0.06"). Both laws are read from the paper's abstract (source tier T1, tweak L0), but they describe small
transiting planets close to their stars; here they are used for every planet at every distance (an extrapolation the paper does not
test: the laws are the only measured ones, so they stand in for the rest, tier T2). Eccentricity is cut at the largest value that keeps
the periapsis outside the star's Roche limit (`max_planet_eccentricity`; the truncation is in the draw, nothing is thrown away) and
the tides of the orbit take what they can (`circularised`): the laws were fitted to planets that tides had already acted on, so the
planets closest in are circularised twice over, a small double count (the fit has no period dependence to remove).

**Inclination.** A planet's plane lies above or below the system's reference plane (the star's equator) by an angle that follows a
Rayleigh law. The mutual inclination of the planets of the Kepler multiples is Rayleigh with a dispersion of 1.0 to 2.2 degrees
[R1311]; the mutual angle of two planets drawn independently from one law is `sqrt(2)` times the
dispersion of that law, so the planet-to-plane dispersion is the mutual dispersion over `sqrt(2)`, here the middle of the range
(source tier T2, tweak L0: the paper was read in an earlier session, not again).

**Moons.** A regular moon (disk, impact, ring or binary partner) has the small forced eccentricity of a resonance or of the neighbours
that keep it from circularising, drawn log-uniform from `REGULAR_MOON_ECCENTRICITY_RANGE` (source tier T3, tweak L0: the regular
satellites of the giants lie between 0.001 and 0.03, the Moon at 0.055; the range is the middle of that). A captured moon starts on
an orbit with the thermal distribution of eccentricities, `f(e) = 2 e`, cut at the largest value that keeps its apoapsis inside the
stable part of its parent's Hill sphere (`MOON_PROGRADE_STABLE_HILL_FRACTION`, tier T5, tweak L0).

**Tides.** Raised in the body by what it orbits, they circularise the orbit of a planet or a captured moon (a regular moon's forced
eccentricity is renewed by its neighbours and stays, so the tide it raises keeps heating it, GEN.737) on `circularisation_time_yr`,
`tau_e = (2 / 21) (Q / k2) (m / M) (a / R)^5 / n` ([R1612], the standard result for a synchronously rotating body;
`e` falls as `exp(-t / tau_e)`). Q and k2 are the body's own (`interior_heat.tidal_response`, GEN.731 and GEN.955: the Andrade k2 and Q of its mass, radius and
mantle temperature for rock, the icy and the gas-giant values for the rest), the same that despin it (`tidal_spin`).

Use:

    e = draw_planet_eccentricity(multiple=True, e_max=max_planet_eccentricity(0.4, star, planet))
    e = circularised(e, age_yr, tau_yr)
    swing = flux_swing(e)                      # peak to peak over the mean, of the light the orbit alone gives
"""

import math

from scipy import special

from planetgen import tuning
from planetgen.physics import constants, interior_heat, tidal_rheology
from planetgen.util import draw

MULTI_SIGMA = 0.083
"""float: Dispersion of the half-Gaussian eccentricity law of a system with several planets ([R1307], abstract; T1, L0)."""
SINGLE_MODERATE_FRACTION = 0.76
"""float: The share of single planets that draw from the broader Rayleigh law ([R1307], abstract; T1, L0)."""
SINGLE_RAYLEIGH_MODE = 0.26
"""float: Where the broader Rayleigh law peaks (its dispersion; [R1307], abstract; T1, L0)."""
MAX_ECCENTRICITY = 0.9
"""float: No planet or moon is drawn above this (a stand-in: the laws put almost nothing there; the Roche cut comes first for a close orbit; T4, L0)."""
ROCHE_COEFFICIENT = 2.44
"""float: The fluid Roche limit is this times the primary's radius times the cube root of the density ratio (Chandrasekhar 1969; T5, L0)."""

MUTUAL_INCLINATION_SIGMA_DEG = 1.6
"""float: Rayleigh dispersion of the mutual inclination of two planets, degrees: the middle of the 1.0 to 2.2 of [R1311] (T2, L0)."""
PLANET_INCLINATION_SIGMA_DEG = MUTUAL_INCLINATION_SIGMA_DEG / math.sqrt(2.0)
"""float: Rayleigh dispersion of a planet's inclination to the system's reference plane, degrees."""

REGULAR_MOON_ECCENTRICITY_RANGE = (0.001, 0.01)
"""tuple: Log-uniform range of the forced eccentricity of a regular moon (the range `hydrosphere` drew its tide from before GEN.727; T3, L0)."""


def _half_gaussian(sigma, e_max):
    """The positive half of a zero-mean Gaussian of dispersion `sigma`, cut at `e_max`, by its inverse CDF (one draw)."""
    reach = special.erf(e_max / (sigma * math.sqrt(2.0)))
    return float(sigma * math.sqrt(2.0) * special.erfinv(draw.random() * reach))


def _rayleigh(sigma, e_max):
    """A Rayleigh law of dispersion `sigma`, cut at `e_max`, by its inverse CDF (one draw)."""
    reach = 1.0 - math.exp(-e_max * e_max / (2.0 * sigma * sigma))
    return sigma * math.sqrt(-2.0 * math.log(1.0 - draw.random() * reach))


def draw_planet_eccentricity(multiple, e_max):
    """
    The eccentricity of a planet, from the law of the system it is in.

    Args:
        multiple (bool): Whether the system has two or more planets (`None`, a system whose planets are not counted, is taken as one that does).
        e_max (float): The largest eccentricity its orbit allows (`max_planet_eccentricity`); 0 gives a circle.

    Returns:
        float: In [0, e_max].
    """
    e_max = min(max(e_max, 0.0), MAX_ECCENTRICITY)
    if e_max <= 0.0:
        return 0.0
    if multiple is False and draw.random() < SINGLE_MODERATE_FRACTION:
        return _rayleigh(SINGLE_RAYLEIGH_MODE, e_max)
    return _half_gaussian(MULTI_SIGMA, e_max)


def draw_planet_inclination_deg():
    """A planet's inclination to the system's reference plane, degrees: a Rayleigh law of `PLANET_INCLINATION_SIGMA_DEG` (one draw)."""
    return PLANET_INCLINATION_SIGMA_DEG * math.sqrt(-2.0 * math.log(1.0 - draw.random()))


def draw_moon_eccentricity(origin, e_max):
    """
    The eccentricity of a moon before the tides act, from how it came to be there (`moon_composition.ORIGINS`).

    A regular moon takes a forced eccentricity (log-uniform over `REGULAR_MOON_ECCENTRICITY_RANGE`); a captured one the thermal law
    `f(e) = 2 e` on `[0, e_max]`, `e = e_max sqrt(u)`. One draw either way.

    Args:
        origin (str or None): The moon's origin; `None` is a regular moon.
        e_max (float): The largest eccentricity its orbit allows (apoapsis inside the parent's stable Hill sphere).
    """
    e_max = min(max(e_max, 0.0), MAX_ECCENTRICITY)
    point = draw.random()
    if origin == "captured":
        return e_max * math.sqrt(point)
    low, high = REGULAR_MOON_ECCENTRICITY_RANGE
    return min(low * (high / low) ** point, e_max)


def draw_argument_of_periapsis_deg():
    """Where periapsis lies from the ascending node, degrees: uniform (one draw)."""
    return draw.uniform(0.0, 360.0)


def periapsis_au(semi_major_axis_au, eccentricity):
    """The closest distance of an orbit, AU."""
    return semi_major_axis_au * (1.0 - (eccentricity or 0.0))


def apoapsis_au(semi_major_axis_au, eccentricity):
    """The farthest distance of an orbit, AU."""
    return semi_major_axis_au * (1.0 + (eccentricity or 0.0))


def mean_flux_factor(eccentricity):
    """The time-average of the light over an orbit against the light at the semi-major axis, `1 / sqrt(1 - e^2)` (Kepler's second law)."""
    return 1.0 / math.sqrt(1.0 - (eccentricity or 0.0) ** 2)


def flux_distance_au(semi_major_axis_au, eccentricity):
    """The radius of the circle that gives the same mean light as the orbit, `a (1 - e^2)^(1/4)`: the light goes as `1 / r^2`."""
    return semi_major_axis_au * (1.0 - (eccentricity or 0.0) ** 2) ** 0.25


def flux_swing(eccentricity):
    """
    The peak-to-peak swing of the light an orbit gives against its mean, `4 e / (1 - e^2)^(3/2)`: the light goes as `1 / r^2`, from
    `1 / (a (1 + e))^2` at apoapsis to `1 / (a (1 - e))^2` at periapsis, and the mean is `1 / (a^2 sqrt(1 - e^2))`. Capped at 1, the
    largest swing `flux_variation` holds (GEN.597).
    """
    e = eccentricity or 0.0
    return min(4.0 * e / (1.0 - e * e) ** 1.5, 1.0)


def roche_limit_au(primary_radius_km, primary_density, body_density):
    """The fluid Roche limit of a body of `body_density` round a primary of `primary_density` and `primary_radius_km`, AU (any one density unit)."""
    return ROCHE_COEFFICIENT * primary_radius_km / constants.AU_TO_KM * (primary_density / body_density) ** (1.0 / 3.0)


def max_planet_eccentricity(semi_major_axis_au, star, planet):
    """
    The largest eccentricity a planet at `semi_major_axis_au` can have and keep its periapsis outside the Roche limit of its star; 0 when
    the orbit is already inside it, `MAX_ECCENTRICITY` when the star or the planet has no size to compare.
    """
    radius_km = getattr(star, "radius", None)
    mass_kg = getattr(star, "mass", None)
    if not radius_km or not mass_kg or not planet.mass or not planet.radius:
        return MAX_ECCENTRICITY
    star_density = mass_kg / (4.0 / 3.0 * math.pi * (radius_km * 1000.0) ** 3)
    planet_density = planet.mass / (4.0 / 3.0 * math.pi * (planet.radius * 1000.0) ** 3)
    limit = roche_limit_au(radius_km, star_density, planet_density)
    return min(max(1.0 - limit / semi_major_axis_au, 0.0), MAX_ECCENTRICITY)


def max_moon_eccentricity(semi_major_axis_au, hill_radius_au):
    """The largest eccentricity that keeps a moon's apoapsis inside the stable part of the parent's Hill sphere; 0 for a moon already at the edge."""
    if hill_radius_au is None or semi_major_axis_au <= 0.0:
        return MAX_ECCENTRICITY
    reach = tuning.MOON_PROGRADE_STABLE_HILL_FRACTION * hill_radius_au
    return min(max(reach / semi_major_axis_au - 1.0, 0.0), MAX_ECCENTRICITY)


def circularisation_time_yr(body, primary_mass_kg, semi_major_axis_au):
    """
    How long the tide raised in `body` by the body it orbits takes to cut its eccentricity by a factor of e, years:
    `(2 / 21) (Q / k2) (m / M) (a / R)^5 / n` [R1612]. `inf` for a body with no size or mass.
    """
    if not body.mass or not body.radius or not primary_mass_kg:
        return math.inf
    a_m = semi_major_axis_au * constants.AU_TO_M
    n = math.sqrt(constants.G * primary_mass_kg / a_m ** 3)
    k2, q = interior_heat.tidal_response(body, tidal_rheology.tidal_frequency_rad_s(n))
    seconds = (2.0 / 21.0) * (q / k2) * (body.mass / primary_mass_kg) * (a_m / (body.radius * 1000.0)) ** 5 / n
    return seconds / constants.SECONDS_PER_YEAR


def circularised(eccentricity, age_yr, time_yr):
    """The eccentricity left after `age_yr` of tides that cut it by e every `time_yr` years."""
    if time_yr <= 0.0:
        return 0.0
    return eccentricity * math.exp(-age_yr / time_yr)
