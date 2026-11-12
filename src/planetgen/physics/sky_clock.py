# planetgen/physics/sky_clock.py

"""
The clock, the calendar of eclipses and the gravitational lens of a body (GEN.1434)
==================================================================================

Three figures a world's sky and a world's clocks carry, worked out on demand from the stored masses, radii and orbits with no draw and no stored column. They are the derived rows of the
object sheet (`describe/clock_rows.py`), the way the moon's longest eclipse in its planet's shadow is.

**The clock** (`clock_lag_fraction`). A clock on the surface of a body that circles a host runs slow against a clock at rest far from the whole nest, by the weak-field sum

    lag = G m / (R c^2)  +  sum over each host of 3 G M_h / (2 r_h c^2)

the first term the body's own surface potential, and for each host it circles (a planet: its star; a moon: its planet and the planet's star) the potential at the orbit, `G M_h / (r_h c^2)`, plus the
speed term `v^2 / (2 c^2)` with `v^2 = G M_h / r_h` of a circular orbit, which together are one and a half times the potential. The velocities of a nest add in quadrature (the cross terms average
out over the orbits), so the terms of the hosts simply add. `lag_per_year_s` turns the fraction into the seconds a year the clock falls behind. The Earth's is 0.49 s a year (of which 0.022 s is
the Earth's own surface), a star's own surface alone is 67 s a year for the Sun and 4.7e3 s for a white dwarf of 0.6 solar masses.

**The calendar of eclipses** (`moon_eclipse`, `planet_eclipse`). The blocker, seen from the observer, is a disk of angular radius `alpha_c = R_c / D` that crosses the star's disk of angular radius
`alpha_s = R_s / d` once in each synodic period `T`, at the angular rate `omega` the two directions part at (the same geometry `physics/eclipses` averages the covered light over, and the
same `beta = beta_max sin(u)` height of a crossing, `u` uniform). A crossing touches the disk only when the height is inside the reach `alpha_c + alpha_s`, which has the share

    P = (2 / pi) asin(min(1, (alpha_c + alpha_s) / beta_max))        of the crossings        (1 when the orbit is in the observer's plane)

so the mean time between eclipses is `T / P`. The longest an eclipse lasts is the central crossing, `2 (alpha_c + alpha_s) / omega` between the first and the last contact, and when the blocker is
the larger disk the central crossing covers the whole star for `2 (alpha_c - alpha_s) / omega` (a total eclipse; an annular one covers at most `(alpha_c / alpha_s)^2` of the star). The Moon
and the Sun give an eclipse about every 1.2 years of a point on the Earth's orbit, up to 2.1 hours between first and last contact; Venus crosses the Sun every 81 years against the 61 of the
real cycle (the circular orbits the model keeps; EC1 of section 6.249).

**The lens** (`focal_distance_m`). A ray that grazes a body of mass `M` and radius `R` bends by `4 G M / (R c^2)` and meets the axis at `R^2 c^2 / (4 G M)` from the body's center: 548 AU for the
Sun, the nearest a spacecraft sees the whole Einstein ring of a distant star. A body whose radius is not larger than four Schwarzschild radii bends by more than a weak field allows and has no
such focus (`None`). The Sun's corona bends light too and pushes the real focus to about 650 AU; the figure here is the vacuum one.

Stand-ins are the rows SK1 to SK3 of section 6.263 of the gap plan.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, eclipses

_C2 = float(constants.SPEED_OF_LIGHT_M_S) ** 2
_HOURS_PER_YEAR = constants.SECONDS_PER_YEAR / 3600.0
_TWO_PI = 2.0 * math.pi


def clock_lag_fraction(own_mass_kg, own_radius_km, hosts=()):
    """
    How much slower a clock on the surface of a body runs than one at rest far from the nest, a fraction (see the module docstring).

    Args:
        own_mass_kg (float): The body's mass.
        own_radius_km (float): The body's radius.
        hosts (iterable): Each body it circles, as `(mass_kg, orbit_radius_km)`: a planet gives its star, a moon its planet then the planet's star.
    """
    lag = 0.0
    if own_mass_kg and own_radius_km:
        lag += constants.G * own_mass_kg / (own_radius_km * 1.0e3 * _C2)
    for mass_kg, orbit_km in hosts:
        if mass_kg and orbit_km:
            lag += 1.5 * constants.G * mass_kg / (orbit_km * 1.0e3 * _C2)
    return lag


def lag_per_year_s(lag_fraction):
    """The seconds a year a clock that runs `lag_fraction` slow falls behind."""
    return lag_fraction * constants.SECONDS_PER_YEAR


def focal_distance_m(mass_kg, radius_km):
    """The distance from a body's center, meters, at which the rays that graze its surface meet (`R^2 c^2 / (4 G M)`), or `None` where the field is too strong for the weak-field figure."""
    if not mass_kg or not radius_km or mass_kg <= 0.0 or radius_km <= 0.0:
        return None
    radius_m = radius_km * 1.0e3
    distance = radius_m * radius_m * _C2 / (4.0 * constants.G * mass_kg)
    schwarzschild = 2.0 * constants.G * mass_kg / _C2
    return distance if radius_m > 4.0 * schwarzschild else None


@dataclass(frozen=True)
class Eclipse:
    """One blocker's eclipses of a star as seen from one body: the mean `interval_yr` between them, the `longest_contact_h` from first contact to last, the `longest_total_h` the whole star is
    covered (0 for a blocker that is the smaller disk) and the `deepest_cover`, the largest share of the star's disk hidden."""

    interval_yr: float
    longest_contact_h: float
    longest_total_h: float
    deepest_cover: float


def touching_share(alpha_c, alpha_s, beta_max):
    """The share of the crossings that touch the star's disk, for a crossing at the height `beta_max sin(u)` (`u` uniform) of the reach `alpha_c + alpha_s`."""
    if alpha_c <= 0.0 or alpha_s <= 0.0:
        return 0.0
    if beta_max <= 0.0:
        return 1.0
    return 2.0 / math.pi * math.asin(min(1.0, (alpha_c + alpha_s) / beta_max))


def _eclipse(alpha_c, alpha_s, beta_max, synodic_yr, omega_rad_yr):
    """The `Eclipse` of the disks `alpha_c` and `alpha_s` (radians) crossing every `synodic_yr` at the angular rate `omega_rad_yr`, or `None` when they never touch."""
    share = touching_share(alpha_c, alpha_s, beta_max)
    if share <= 0.0 or synodic_yr <= 0.0 or omega_rad_yr <= 0.0:
        return None
    reach = alpha_c + alpha_s
    total_h = 2.0 * max(alpha_c - alpha_s, 0.0) / omega_rad_yr * _HOURS_PER_YEAR
    deepest = 1.0 if alpha_c >= alpha_s else (alpha_c / alpha_s) ** 2
    return Eclipse(synodic_yr / share, 2.0 * reach / omega_rad_yr * _HOURS_PER_YEAR, total_h, deepest)


def moon_eclipse(star_radius_km, star_distance_km, planet_period_yr, moon_radius_km, moon_distance_km, moon_period_yr, inclination_deg):
    """The eclipses of a star by a moon for a body (its planet) at `star_distance_km` from the star: the moon `moon_distance_km` from the planet, its orbit tilted `inclination_deg` to the planet's plane.
    `None` for a moon that never touches the star's disk from there."""
    if not (star_radius_km and star_distance_km and moon_radius_km and moon_distance_km and moon_period_yr and moon_period_yr > 0.0):
        return None
    rate = abs(1.0 / moon_period_yr - (1.0 / planet_period_yr if planet_period_yr and planet_period_yr > 0.0 else 0.0))
    if rate <= 0.0:
        return None
    synodic = 1.0 / rate
    return _eclipse(moon_radius_km / moon_distance_km, star_radius_km / star_distance_km, math.radians(inclination_deg), synodic, _TWO_PI / synodic)


def planet_eclipse(star_radius_km, observer_km, observer_period_yr, blocker_radius_km, blocker_km, blocker_period_yr, inclination_deg):
    """The eclipses (transits) of a star seen from a planet at `observer_km` from it by an inner planet at `blocker_km`: `None` for a blocker that is not nearer the star or never touches its disk."""
    if not (star_radius_km and blocker_radius_km and observer_km > blocker_km > 0.0 and observer_period_yr and blocker_period_yr
            and observer_period_yr > 0.0 and blocker_period_yr > 0.0):
        return None
    gap_km = observer_km - blocker_km
    synodic = 1.0 / abs(1.0 / blocker_period_yr - 1.0 / observer_period_yr)
    beta_max = blocker_km * math.sin(math.radians(inclination_deg)) / gap_km
    return _eclipse(blocker_radius_km / gap_km, star_radius_km / observer_km, beta_max, synodic, eclipses.transit_rate(observer_km, observer_period_yr, blocker_km, blocker_period_yr))
