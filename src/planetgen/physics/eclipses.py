# planetgen/physics/eclipses.py

"""
Eclipses of a star by a body of the system, orbit-averaged (GEN.1462)
====================================================================

Every body of a system can stand in front of a star as a body sees it: a sibling planet nearer the star crosses the host (a transit), a moon crosses it for its planet or for a
sister moon, one star of a close pair crosses the other for a circumbinary planet. `daylight.parent_shadow_share` has the one case a moon always has, its own planet's shadow; this module
is the rest, for any blocker, with no draw.

**The covered share.** The blocker, seen from the observer, is a disk of angular radius `alpha_c = R_c / D` that crosses the star's disk of angular radius `alpha_s = R_s / d`
once in each synodic period `T`, at the angular rate `omega` the two directions part at. Along a track at the vertical offset `beta` it covers the area `A(phi, beta)` (the overlap of
two disks, `lens_area`) of the star's `pi alpha_s^2`, so over one crossing the star is covered for the time-integral `I(beta) / omega`, `I(beta) = integral A dphi`. The offset is the
height of the blocker above the line to the star at the crossing, `beta = beta_max sin(u)` with `u` uniform (a blocker whose orbit is tilted by `i` to the observer's swings through
`+-beta_max` the way a node does); a blocker in the observer's plane has `beta = 0`. The mean share of the star's light covered over a long time is

    mean_covered = <I(beta)> / (pi alpha_s^2 omega T)        and the share of the time the star is touched at all is    <J(beta)> / (omega T),

with `J(beta)` the length of the track along which `A > 0`.

**Three geometries.**

* `orbiting_the_star` (a sibling planet inside the observer's orbit; both circle the star): `D = a_obs - a_c`, `omega = |v_c - v_obs| / D` at conjunction, `T = 1 / |1/P_c - 1/P_obs|`,
  `beta_max = a_c sin(i) / D`. The star is at `a_obs` and subtends `R_s / a_obs`.
* `orbiting_the_observer` (a moon of the observer, or of a sister moon taken at the planet's distance): `D = a`, the blocker's direction goes round once in the synodic period
  `T = 1 / |1/P_m - 1/P_obs|`, `omega = 2 pi / T`, `beta_max = i`.
* `pair` (one star of a close pair in front of the other, for a circumbinary planet at `d` from the barycenter): the stars swing across each other with the pair's period and the
  planet, at the pair's plane within `+- (R1 + R2) / s`, sees each cross the other twice a period; the covered share of the star behind is the share of the pair period that the
  disks overlap times the mean area, which `mean_covered` gives with `omega = 2 pi / P / (d / s)` taken at the crossing (the stars' relative speed `2 pi s / P` seen from `d`).

The orbits are circular and the observer is taken as a point, as everywhere in the model (the disk of a world is a few arc seconds in the sky of the Moon-Earth pair). The textbook
cases are the Moon and the Sun for the Earth (the Moon touches the Sun's disk for 1.5e-4 of the time and covers 3.6e-5 of its light on average, from `5.1` degrees of inclination) and the transits of
Venus (eight in 243 years, about six hours each, 2.4e-5 of the time).

Stand-ins are the rows EC1 to EC4 of section 6.249 of the gap plan. GEN.1484 applies `light_factor_of` (a planet of a close pair) to the light the planet is generated with.
"""

import functools
import math

import numpy as np

from planetgen.physics import constants, spin

TRACK_SAMPLES = 400
"""int: Points along a crossing the area is integrated at."""

OFFSET_SAMPLES = 24
"""int: Points of `u` the vertical offset is averaged over."""

_TWO_PI = 2.0 * math.pi


def lens_area(radius_a, radius_b, separation):
    """The area where two disks of radii `radius_a` and `radius_b` overlap, their centers `separation` apart (all in one unit)."""
    if separation >= radius_a + radius_b:
        return 0.0
    small, large = min(radius_a, radius_b), max(radius_a, radius_b)
    if separation <= large - small:
        return math.pi * small * small
    a2, b2, d2 = radius_a * radius_a, radius_b * radius_b, separation * separation
    part_a = a2 * math.acos(max(-1.0, min(1.0, (d2 + a2 - b2) / (2.0 * separation * radius_a))))
    part_b = b2 * math.acos(max(-1.0, min(1.0, (d2 + b2 - a2) / (2.0 * separation * radius_b))))
    kite = 0.5 * math.sqrt(max(0.0, (-separation + radius_a + radius_b) * (separation + radius_a - radius_b)
                               * (separation - radius_a + radius_b) * (separation + radius_a + radius_b)))
    return part_a + part_b - kite


def _track(alpha_c, alpha_s, beta):
    """`(I, J)` along a crossing at the offset `beta` (radians): `I` the integral of the overlap area over the angle along the track, `J` the length of track where the disks touch."""
    reach = alpha_c + alpha_s
    if abs(beta) >= reach:
        return 0.0, 0.0
    half = math.sqrt(reach * reach - beta * beta)
    xs = np.linspace(-half, half, TRACK_SAMPLES)
    areas = np.array([lens_area(alpha_c, alpha_s, math.hypot(x, beta)) for x in xs])
    return float(np.trapezoid(areas, xs)), 2.0 * half


def crossing(alpha_c, alpha_s, beta_max):
    """`(I, J)` averaged over the vertical offset `beta = beta_max sin(u)`, `u` uniform: the covered area integrated over a crossing and the length of a crossing, radians. Only the offsets
    inside the reach of the two disks cross at all, so the average runs over those (`u` within `asin(reach / beta_max)`), weighted by their share of the arcsine law."""
    if alpha_c <= 0.0 or alpha_s <= 0.0:
        return 0.0, 0.0
    if beta_max <= 0.0:
        return _track(alpha_c, alpha_s, 0.0)
    limit = math.asin(min(1.0, (alpha_c + alpha_s) / beta_max))
    total_i = total_j = 0.0
    for k in range(OFFSET_SAMPLES):
        u = limit * (2.0 * (k + 0.5) / OFFSET_SAMPLES - 1.0)
        i, j = _track(alpha_c, alpha_s, beta_max * math.sin(u))
        total_i, total_j = total_i + i, total_j + j
    weight = 2.0 * limit / math.pi / OFFSET_SAMPLES
    return total_i * weight, total_j * weight


def _shares(alpha_c, alpha_s, beta_max, omega_times_period):
    """`(mean_covered, touched_share)` for a crossing: see the module docstring."""
    if omega_times_period <= 0.0:
        return 0.0, 0.0
    crossed_i, crossed_j = crossing(alpha_c, alpha_s, beta_max)
    return crossed_i / (math.pi * alpha_s * alpha_s * omega_times_period), crossed_j / omega_times_period


def orbiting_the_star(star_radius_km, observer_au, observer_period_yr, blocker_radius_km, blocker_au, blocker_period_yr, inclination_deg, km_per_au):
    """`(mean_covered, touched_share)` of a star's light at a body `observer_au` from it, by a blocker that circles the same star inside the observer's orbit (a sibling planet)."""
    if not (observer_au > blocker_au > 0.0) or observer_period_yr <= 0.0 or blocker_period_yr <= 0.0 or not star_radius_km or not blocker_radius_km:
        return 0.0, 0.0
    gap_au = observer_au - blocker_au
    alpha_c = blocker_radius_km / (gap_au * km_per_au)
    alpha_s = star_radius_km / (observer_au * km_per_au)
    synodic = 1.0 / abs(1.0 / blocker_period_yr - 1.0 / observer_period_yr)
    speed = abs(_TWO_PI * blocker_au / blocker_period_yr - _TWO_PI * observer_au / observer_period_yr)
    omega = speed / gap_au                                        # rad per yr at the conjunction
    beta_max = blocker_au * math.sin(math.radians(inclination_deg)) / gap_au
    return _shares(alpha_c, alpha_s, beta_max, omega * synodic)


def orbiting_the_observer(star_radius_km, star_au, observer_period_yr, blocker_radius_km, blocker_distance_km, blocker_period_yr, inclination_deg, km_per_au):
    """`(mean_covered, touched_share)` of a star's light at a body by a blocker that circles that body (its moon): `star_au` the star's distance, `inclination_deg` the blocker's orbit against the body's."""
    if blocker_distance_km <= 0.0 or not (star_radius_km and blocker_radius_km and star_au) or blocker_period_yr <= 0.0:
        return 0.0, 0.0
    alpha_c = blocker_radius_km / blocker_distance_km
    alpha_s = star_radius_km / (star_au * km_per_au)
    rate = abs(1.0 / blocker_period_yr - (1.0 / observer_period_yr if observer_period_yr and observer_period_yr > 0.0 else 0.0))
    if rate <= 0.0:
        return 0.0, 0.0
    return _shares(alpha_c, alpha_s, math.radians(inclination_deg), _TWO_PI)


def pair(front_radius_km, behind_radius_km, separation_au, pair_period_yr, barycenter_au, inclination_deg, km_per_au):
    """
    `(mean_covered_of_the_star_behind, touched_share)` for a circumbinary planet `barycenter_au` from the pair's barycenter, the pair at `separation_au` and `pair_period_yr`, the planet's
    orbit tilted `inclination_deg` to the pair's plane: the star in front covers the one behind twice a pair period. The seen separation is `s / d`, the crossing rate `2 pi s / (P d)`.
    """
    if not (front_radius_km and behind_radius_km and separation_au > 0.0 and pair_period_yr > 0.0 and barycenter_au > 0.0):
        return 0.0, 0.0
    alpha_c = front_radius_km / (barycenter_au * km_per_au)
    alpha_s = behind_radius_km / (barycenter_au * km_per_au)
    beta_max = (separation_au / barycenter_au) * math.sin(math.radians(inclination_deg))
    omega = _TWO_PI * separation_au / (pair_period_yr * barycenter_au)        # rad per yr
    covered, touched = _shares(alpha_c, alpha_s, beta_max, omega * pair_period_yr / 2.0)      # two crossings a period
    return covered, touched


@functools.lru_cache(maxsize=4096)
def _pair_light_lost(radius_a, luminosity_a, radius_b, luminosity_b, separation_au, pair_period_yr, distance_au, inclination_deg, km_per_au):
    """`pair_light_lost` on numbers (the cache key: a planet asks again for each light it takes, and the crossing integral costs tens of milliseconds)."""
    total = luminosity_a + luminosity_b
    if not total > 0.0:
        return 0.0
    (large_r, large_l), (small_r, small_l) = sorted(((radius_a, luminosity_a), (radius_b, luminosity_b)), key=lambda star: -star[0])
    small_covered, _ = pair(large_r, small_r, separation_au, pair_period_yr, distance_au, inclination_deg, km_per_au)
    large_covered, _ = pair(small_r, large_r, separation_au, pair_period_yr, distance_au, inclination_deg, km_per_au)
    return (small_covered * small_l + large_covered * large_l) / total / 2.0


def pair_light_lost(star_a, star_b, separation_au, pair_period_yr, distance_au, inclination_deg, km_per_au):
    """
    The share of the light of a close pair that a planet `distance_au` from its barycenter loses to the stars eclipsing each other, averaged over a long time (GEN.1484). Each period the
    larger star crosses the smaller one's disk once (covering a share of the smaller one's light, `pair(larger, smaller, ...)`) and the smaller crosses the larger once, so each
    crossing type is counted once a period, half of what `pair` gives (it counts two crossings of one kind a period). The lost light is the covered share of each star behind times its luminosity,
    over the pair's total.
    """
    return _pair_light_lost(float(star_a.radius), float(star_a.luminosity), float(star_b.radius), float(star_b.luminosity), float(separation_au), float(pair_period_yr),
                            float(distance_au), float(inclination_deg), float(km_per_au))


def light_factor_of(planet):
    """
    The share of its host's light a planet keeps after the eclipses of its sky (GEN.1484): 1 less `pair_light_lost` for a planet of a close pair (its host is the pair's proxy, which holds both stars, their
    separation and period, and the plane they swing in), `None` for every other planet (a lone star, a moon, a planet of one star of a wide pair), so their light is as it was. The sibling planets and
    moons that cross a star cover under 1e-4 of the light for 93 percent of planets (section 6.249 of the gap plan, rows EC1 to EC4) and are not applied.
    """
    host = getattr(planet, "star", None)
    stars = getattr(host, "stars", None)
    if getattr(planet, "is_moon", False) or not stars or len(stars) != 2 or not getattr(planet, "distance", None):
        return None
    period = getattr(host, "binary_mutual_orbital_period_years", None)
    plane = (getattr(host, "binary_mutual_orbital_inclination_deg", None), getattr(host, "binary_mutual_orbital_ascending_node_deg", None))
    if not period or None in plane:
        return None
    tilt = spin.obliquity_between(spin.orbit_normal(getattr(planet, "orbital_inclination_deg", None) or 0.0, getattr(planet, "orbital_ascending_node_deg", None) or 0.0),
                                  spin.orbit_normal(*plane))
    lost = pair_light_lost(stars[0], stars[1], host.binary_separation_au, period, float(planet.distance), tilt, constants.AU_TO_KM)
    return None if lost <= 0.0 else max(0.0, 1.0 - lost)
