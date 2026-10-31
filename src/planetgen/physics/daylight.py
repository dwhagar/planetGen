# planetgen/physics/daylight.py

"""
Day and night time shares, and the hours of daylight (GEN.767; Boss, 7:54 PM and 8:22 PM EDT 2026-10-24).

Every planet and moon turns under its star, so each has a share of the time its star is above the horizon. This module derives
that share and the length of the daylight from the geometry the generator already stores, with no draw and no stand-in:
the rotation period, the orbital period round the star, the axial tilt (the spin axis and the orbit normal), the star's
radius and its distance. For a moon it adds the extra night when its planet stands between it and the star.

* **The share at the equator** (`equator_day_share`): the star is up while its centre is above `-alpha`, `alpha` being the
  angular radius of the star (its upper limb counts, no refraction): the half day arc `H` has `cos H = -sin(alpha)/cos(delta)`
  and the share is `H / pi`, averaged over a circular year, the declination following `sin(delta) = sin(e) sin(lambda)` with
  `e` the angle between the spin axis and the normal of the orbit round the star. Earth gives 50.15 percent, 12.04 h of a
  24 h day.
* **The mean over the surface** (`mean_day_share`): the lit part of a sphere seen by a star of angular radius `alpha` is
  `(1 + sin(alpha)) / 2` of it at every instant, whatever the tilt, so the time share averaged over the surface is that number
  too (a point star gives one half). It holds for a locked world as well.
* **The solar day** (`solar_day_hours`): the star returns to the same place in `1 / |1/P_rot -/+ 1/P_orb|` (the sign for a
  prograde or a retrograde spin), `P_orb` the period round the star (a moon's star orbit is its parent's). It is `None` for a world
  that turns once an orbit: it has no day, only the day side, the terminator and the night side (the day side is lit all the time, the terminator
  half of it, the night side never). `daylight_hours` is the solar day times the equatorial share.
* **The polar circle** (`polar_circle_deg`): the latitude beyond which the star does not set for at least a day,
  `90 - e - alpha` degrees (zero or less becomes 0: the whole world has polar days and nights; a tilt past 90 degrees
  counts as `180 - tilt`).
* **The parent's shadow** (`parent_shadow_share`): a moon is in its planet's shadow while the centre of the star is behind the
  planet's disc, which the planet fills with the angular radius `theta = asin(R_planet / a)` as the moon sees it. With the
  angle `u` along the orbit from the point towards the star and `beta` the angle between the star's direction and the moon's
  orbit plane, the condition is `-cos(beta) cos(u) >= cos(theta)`, so a share `arccos(cos(theta)/cos(beta))/pi` of each orbit
  when `beta < theta`, none otherwise. The year turns the star round the planet's orbit plane, so `sin(beta) = sin(Theta)
  |sin(lambda)|`, `Theta` the angle between the moon's orbit normal and the planet's. The moon's day share (both of them) is
  cut by that share: `day_share * (1 - parent_shadow_share)`. For a moon locked to its planet the extra night falls on its
  near side only (the hemisphere facing the planet is the one that faces the star while the moon is in the shadow).

The orbits are circular, as everywhere in the model, so the libration of an eccentric orbit is not there to note; the star is
taken at the parent's distance for a moon (the moon's own distance from its planet is a few light seconds).
"""

import math

from planetgen.physics import constants, day_night, spin
from planetgen.physics import unit_registry as units

_yr_to_h = units.converter("yr", "h")

FIELDS = ("day_share_equator", "day_share_mean", "parent_shadow_share", "solar_day_hours", "daylight_hours", "polar_circle_deg",
          "day_night_state", "night_hours", "days_per_year")
"""tuple: The stored columns on `planets` and `moons`. `day_share_equator`, `day_share_mean` and `parent_shadow_share` are shares
of the time (0 to 1); `solar_day_hours` and `daylight_hours` are hours; `polar_circle_deg` is a latitude in degrees. A
world that turns once an orbit has no `day_share_equator`, `solar_day_hours` or `daylight_hours` (`None`); only a moon has a
`parent_shadow_share` (`None` for a planet). `day_night_state` is the word of `day_night.state`, `night_hours` the solar day less the daylight and `days_per_year` the
solar days in one orbit (GEN.1097); the last two are `None` where the solar day is."""

YEAR_SAMPLES = 720
"""int: The points of the circular year the declination and the shadow are averaged over."""

LOCK_TOLERANCE = 1e-3
"""float: How close the rotation period is to the orbital period to count as turning once an orbit (the generator sets them equal)."""


def angular_radius(star_radius_km, distance_au):
    """The angular radius of the star as the world sees it, radians; 0 for a star of unknown size (a point source)."""
    if not star_radius_km or not distance_au or distance_au <= 0.0:
        return 0.0
    return math.asin(min(1.0, star_radius_km / (distance_au * constants.AU_TO_KM)))


def folded_tilt_deg(tilt_deg):
    """The obliquity as an angle from 0 to 90 degrees (a tilt past 90 is a retrograde spin with the same seasons)."""
    tilt = abs(tilt_deg) % 360.0
    tilt = 360.0 - tilt if tilt > 180.0 else tilt
    return 180.0 - tilt if tilt > 90.0 else tilt


def day_share_at(latitude_deg, declination_deg, alpha):
    """The share of one turn the star's upper limb is above the horizon at a latitude, the star at a declination (a point of the
    year), for a star of angular radius `alpha` radians."""
    latitude, declination = math.radians(latitude_deg), math.radians(declination_deg)
    denominator = math.cos(latitude) * math.cos(declination)
    numerator = -math.sin(alpha) - math.sin(latitude) * math.sin(declination)
    if denominator <= 1e-12:
        return 1.0 if numerator < 0.0 else 0.0
    ratio = numerator / denominator
    if ratio >= 1.0:
        return 0.0
    if ratio <= -1.0:
        return 1.0
    return math.acos(ratio) / math.pi


def equator_day_share(tilt_deg, alpha):
    """The share of the time the star is up at the equator over a circular year (before any shadow of a parent)."""
    sin_tilt = math.sin(math.radians(folded_tilt_deg(tilt_deg)))
    total = 0.0
    for k in range(YEAR_SAMPLES):
        longitude = 2.0 * math.pi * (k + 0.5) / YEAR_SAMPLES
        declination = math.degrees(math.asin(sin_tilt * math.sin(longitude)))
        total += day_share_at(0.0, declination, alpha)
    return total / YEAR_SAMPLES


def mean_day_share(alpha):
    """The share of the time the star is up, averaged over the whole surface: the lit part of a sphere, `(1 + sin alpha) / 2`."""
    return 0.5 * (1.0 + math.sin(alpha))


def polar_circle_deg(tilt_deg, alpha):
    """The latitude beyond which the star stays above (or below) the horizon for at least a day, degrees; 0 when the whole world has
    them."""
    return max(0.0, 90.0 - folded_tilt_deg(tilt_deg) - math.degrees(alpha))


def is_locked(rotation_period_hours, orbital_period_hours):
    """Whether the body turns once an orbit round its star (its rotation period is the orbital period)."""
    if not rotation_period_hours or not orbital_period_hours:
        return False
    return abs(rotation_period_hours - orbital_period_hours) <= LOCK_TOLERANCE * orbital_period_hours


def solar_day_hours(rotation_period_hours, orbital_period_hours, prograde=True):
    """The time from one noon to the next, hours; `None` for a world that turns once an orbit (no day)."""
    if not rotation_period_hours or not orbital_period_hours or rotation_period_hours <= 0.0 or orbital_period_hours <= 0.0:
        return None
    if prograde and is_locked(rotation_period_hours, orbital_period_hours):
        return None
    rate = 1.0 / rotation_period_hours + (-1.0 if prograde else 1.0) / orbital_period_hours
    return 1.0 / abs(rate)


def shadow_share_at(orbit_radius_km, parent_radius_km, beta_deg):
    """The share of one orbit a moon spends hiding the centre of the star behind its planet while the star stands `beta_deg` degrees
    from the moon's orbit plane: `arccos(cos(theta) / cos(beta)) / pi`, none once `beta` passes the planet's angular radius."""
    if not orbit_radius_km or not parent_radius_km or orbit_radius_km <= parent_radius_km:
        return 0.0
    cos_theta = math.cos(math.asin(parent_radius_km / orbit_radius_km))
    cos_beta = math.cos(math.radians(beta_deg))
    ratio = cos_theta / cos_beta if cos_beta > 0.0 else 2.0
    return math.acos(ratio) / math.pi if ratio < 1.0 else 0.0


def longest_eclipse_hours(orbit_radius_km, parent_radius_km, orbit_period_hours):
    """The longest time a moon spends hiding the centre of the star behind its planet in one pass, hours (UX.205): the pass through the
    middle of the shadow (the star in the moon's orbit plane), which is `theta / pi` of one orbit, `theta` the angular radius of the
    planet as the moon sees it. Zero when the moon is not outside the planet or has no period."""
    if not orbit_period_hours or orbit_period_hours <= 0.0:
        return 0.0
    return shadow_share_at(orbit_radius_km, parent_radius_km, 0.0) * orbit_period_hours


def parent_shadow_share(orbit_radius_km, parent_radius_km, mutual_inclination_deg):
    """The share of a moon's time its planet hides the centre of the star over a circular year (see the module text): the orbit radius
    and the planet's radius in kilometres, and the angle between the moon's orbit normal and the planet's, degrees."""
    sin_incl = math.sin(math.radians(mutual_inclination_deg))
    total = 0.0
    for k in range(YEAR_SAMPLES):
        longitude = 2.0 * math.pi * (k + 0.5) / YEAR_SAMPLES
        total += shadow_share_at(orbit_radius_km, parent_radius_km, math.degrees(math.asin(abs(sin_incl * math.sin(longitude)))))
    return total / YEAR_SAMPLES


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _spin_axis(body):
    axis = tuple(getattr(body, name, None) for name in ("spin_axis_x", "spin_axis_y", "spin_axis_z"))
    return None if any(value is None for value in axis) else axis


def _orbit_normal(body):
    return spin.orbit_normal(getattr(body, "orbital_inclination_deg", None) or 0.0,
                             getattr(body, "orbital_ascending_node_deg", None) or 0.0)


def compute(body, parent=None):
    """
    The `FIELDS` of a planet or moon as a dict (`None` for what the body has not): `parent` is the planet a moon orbits. A body
    with no rotation period gets all `None`. Pure: no draws.
    """
    out = dict.fromkeys(FIELDS)
    rotation = getattr(body, "rotation_period_hours", None)
    star = getattr(body, "star", None)
    host = parent if parent is not None else body
    solid = getattr(body, "body_type", "t") != "g"
    if not solid:
        out["day_night_state"] = day_night.NO_SURFACE
    if not rotation or rotation <= 0.0 or star is None:
        return out
    distance_au = getattr(host, "distance", None)
    orbital_hours = _yr_to_h(getattr(host, "period", None) or 0.0)
    alpha = angular_radius(getattr(star, "radius", None), distance_au)
    if parent is None:
        tilt = getattr(body, "axial_tilt_deg", None)
    else:
        axis = _spin_axis(body)
        tilt = None if axis is None else spin.obliquity_between(axis, _orbit_normal(parent))
    if tilt is None:
        return out
    prograde = tilt <= 90.0
    out["day_share_mean"] = mean_day_share(alpha)
    out["polar_circle_deg"] = polar_circle_deg(tilt, alpha)
    day = solar_day_hours(rotation, orbital_hours, prograde)
    if day is not None:
        out["day_share_equator"] = equator_day_share(tilt, alpha)
        out["solar_day_hours"] = day
    if parent is not None:
        mutual = spin.obliquity_between(_orbit_normal(body), _orbit_normal(parent))
        shadow = parent_shadow_share(body.distance * constants.AU_TO_KM, getattr(parent, "radius", None), mutual)
        out["parent_shadow_share"] = shadow
        out["day_share_mean"] *= 1.0 - shadow
        if out["day_share_equator"] is not None:
            out["day_share_equator"] *= 1.0 - shadow
    if out["day_share_equator"] is not None:
        out["daylight_hours"] = out["solar_day_hours"] * out["day_share_equator"]
    out["day_night_state"] = day_night.state(rotation, getattr(host, "period", None), out["polar_circle_deg"], solid)
    out["night_hours"] = day_night.night_hours(out["solar_day_hours"], out["daylight_hours"])
    out["days_per_year"] = day_night.days_per_year(orbital_hours, out["solar_day_hours"])
    return out


def generate(body, parent=None):
    """Sets `FIELDS` on a planet or moon (`parent` the planet a moon orbits)."""
    for name, value in compute(body, parent).items():
        setattr(body, name, value)
