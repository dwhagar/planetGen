# planetgen/physics/daylight.py

"""
Day and night time shares, and the hours of daylight (GEN.767; Boss, 7:54 PM and 8:22 PM EDT 2026-10-24).

Every planet and moon turns under its star, so each has a share of the time its star is above the horizon. This module derives
that share and the length of the daylight from the geometry the generator already stores, with no draw and no stand-in:
the rotation period, the orbital period round the star, the axial tilt (the spin axis and the orbit normal), the star's
radius and its distance. For a moon it adds the extra night when its planet stands between it and the star.

* **The share at the equator** (`equator_day_share`): the star is up while its center is above `-alpha`, `alpha` being the
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
* **The parent's shadow** (`parent_shadow_share`): a moon is in its planet's shadow while the center of the star is behind the
  planet's disk, which the planet fills with the angular radius `theta = asin(R_planet / a)` as the moon sees it. With the
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

from planetgen.physics import constants, day_night, orbit_elements, planet_type, radiation_field, sky_suns, spin
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

DAY_SUN_W_M2 = 0.036
"""float: The light, W/m2, that makes a star a sun of the day (GEN.1461): civil twilight, 3.4 lux of the 1.3e5 lux of the Sun at Earth (2.6e-5 of its 1361 W/m2); a star that gives a body less than
this does not make a day of its own (a stand-in, section 6.243 of the gap plan)."""

SKY_YEAR_SAMPLES = 36
"""int: The points of the circular year the suns' places are worked at when more than one sun counts (GEN.1461)."""

SKY_LATITUDE_STEPS = 48
"""int: The strips of equal area the mean over the surface is worked in when more than one sun counts (GEN.1461)."""

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


def lit_arcs_share(arcs):
    """
    The share of one turn that is lit by at least one sun (GEN.1461): `arcs` is `[(center, half_width)]` in radians, the hour angle a sun is up around (its right ascension) and half the
    time it is up, as an angle (`pi` for a sun that does not set). A circle of arcs merged where they overlap: two suns that stand together light the arc of one, two that stand
    opposite light both.
    """
    spans = []
    for center, half in arcs:
        if half >= math.pi:
            return 1.0
        if half <= 0.0:
            continue
        low = (center - half) % (2.0 * math.pi)
        high = low + 2.0 * half
        if high > 2.0 * math.pi:
            spans.append((low, 2.0 * math.pi))
            spans.append((0.0, high - 2.0 * math.pi))
        else:
            spans.append((low, high))
    if not spans:
        return 0.0
    spans.sort()
    covered, end = 0.0, spans[0][0]
    for low, high in spans:
        if high <= end:
            continue
        covered += high - max(low, end)
        end = high
    return min(1.0, covered / (2.0 * math.pi))


def _sky_sources(tilt_rad, distance_au, counted, longitude):
    """`[(sin_declination, right_ascension, angular_radius)]` of the counted suns of a planet at a point of its year: the host stands at the longitude `longitude`, the others at the place
    their offset from the host and the planet's own place (`distance_au` on the other side of the host) give; `counted` is `[(star, offset_au or None)]`, the host's offset `None`."""
    sin_e, cos_e = math.sin(tilt_rad), math.cos(tilt_rad)
    host = (math.cos(longitude), math.sin(longitude), 0.0)
    out = []
    for source, offset in counted:
        if offset is None:
            vector, separation = host, distance_au
        else:
            vector = tuple(offset[axis] + distance_au * host[axis] for axis in range(3))
            separation = math.sqrt(sum(c * c for c in vector))
            vector = tuple(c / separation for c in vector)
        out.append((sin_e * vector[1] + cos_e * vector[2], math.atan2(cos_e * vector[1] - sin_e * vector[2], vector[0]),
                    angular_radius(getattr(source, "radius", None), separation)))
    return out


def sky_day_shares(tilt_rad, distance_au, counted):
    """
    `(equator_share, mean_share)` of the time some counted sun is up, for a planet that circles the first star of `counted` at `distance_au` (GEN.1461): the union of the arcs the suns
    light, `lit_arcs_share`, at the equator and averaged over the surface, over a circular year. The hour angle of a sun is its right ascension against the turning ground; the
    surface mean is taken in strips of equal area. A lone counted sun gives `equator_day_share` and `mean_day_share`.
    """
    equator = mean = 0.0
    for k in range(SKY_YEAR_SAMPLES):
        longitude = 2.0 * math.pi * (k + 0.5) / SKY_YEAR_SAMPLES
        suns = _sky_sources(tilt_rad, distance_au, counted, longitude)
        equator += lit_arcs_share([(ra, math.pi * day_share_at(0.0, math.degrees(math.asin(sin_dec)), alpha)) for sin_dec, ra, alpha in suns])
        strips = 0.0
        for j in range(SKY_LATITUDE_STEPS):
            latitude = math.degrees(math.asin(-1.0 + 2.0 * (j + 0.5) / SKY_LATITUDE_STEPS))
            strips += lit_arcs_share([(ra, math.pi * day_share_at(latitude, math.degrees(math.asin(sin_dec)), alpha)) for sin_dec, ra, alpha in suns])
        mean += strips / SKY_LATITUDE_STEPS
    return equator / SKY_YEAR_SAMPLES, mean / SKY_YEAR_SAMPLES


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
    """The share of one orbit a moon spends hiding the center of the star behind its planet while the star stands `beta_deg` degrees
    from the moon's orbit plane: `arccos(cos(theta) / cos(beta)) / pi`, none once `beta` passes the planet's angular radius."""
    if not orbit_radius_km or not parent_radius_km or orbit_radius_km <= parent_radius_km:
        return 0.0
    cos_theta = math.cos(math.asin(parent_radius_km / orbit_radius_km))
    cos_beta = math.cos(math.radians(beta_deg))
    ratio = cos_theta / cos_beta if cos_beta > 0.0 else 2.0
    return math.acos(ratio) / math.pi if ratio < 1.0 else 0.0


def longest_eclipse_hours(orbit_radius_km, parent_radius_km, orbit_period_hours):
    """The longest time a moon spends hiding the center of the star behind its planet in one pass, hours (UX.205): the pass through the
    middle of the shadow (the star in the moon's orbit plane), which is `theta / pi` of one orbit, `theta` the angular radius of the
    planet as the moon sees it. Zero when the moon is not outside the planet or has no period."""
    if not orbit_period_hours or orbit_period_hours <= 0.0:
        return 0.0
    return shadow_share_at(orbit_radius_km, parent_radius_km, 0.0) * orbit_period_hours


def parent_shadow_share(orbit_radius_km, parent_radius_km, mutual_inclination_deg):
    """The share of a moon's time its planet hides the center of the star over a circular year (see the module text): the orbit radius
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


def _orbit_normal(body, planet=None):
    """The normal of the orbit of `body` in the system's frame: a moon's elements are against its `planet`'s equator (GEN.1352)."""
    if planet is not None:
        return spin.moon_orbit_normal(body, planet)
    return spin.orbit_normal(getattr(body, "orbital_inclination_deg", None) or 0.0,
                             getattr(body, "orbital_ascending_node_deg", None) or 0.0)


def _counted_sky(body, parent, star, host):
    """
    `(tilt_rad, flux_distance_au, counted)` when more than the host's own light makes the day of this body (GEN.1461), else `None` (the host alone, or a lone star, is the day as it was):
    `counted` is `[(star, offset_au or None)]` of the host (offset `None`) and each companion that puts at least `DAY_SUN_W_M2` on the body, the offset in the equinox frame.
    """
    if not radiation_field.companions_of(star):
        return None
    axis = _spin_axis(body)
    distance_au = getattr(host, "distance", None)
    if axis is None or not distance_au:
        return None
    flux_au = orbit_elements.flux_distance_au(distance_au, getattr(host, "eccentricity", None) or 0.0)
    tilt, total, found = sky_suns.placed(star, axis, _orbit_normal(parent) if parent is not None else _orbit_normal(body), flux_au)
    if not found:
        return None
    counted = [(star, None)] if total - sum(source.light_w_m2 for source in found) >= DAY_SUN_W_M2 else []
    counted += [(source.star, source.offset_au) for source in found if source.light_w_m2 >= DAY_SUN_W_M2]
    if not counted or (len(counted) == 1 and counted[0][1] is None):
        return None
    return tilt, flux_au, counted


def compute(body, parent=None):
    """
    The `FIELDS` of a planet or moon as a dict (`None` for what the body has not): `parent` is the planet a moon orbits. A body
    with no rotation period gets all `None`. Pure: no draws.
    """
    out = dict.fromkeys(FIELDS)
    rotation = getattr(body, "rotation_period_hours", None)
    star = getattr(body, "star", None)
    host = parent if parent is not None else body
    solid = planet_type.has_ground(body)      # GEN.1334: a planet under a deep envelope has no surface for a day to fall on
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
    sky = _counted_sky(body, parent, star, host)       # GEN.1461: the day is the union of the arcs every counted sun lights
    day = solar_day_hours(rotation, orbital_hours, prograde)
    if sky is not None:
        equator_share, out["day_share_mean"] = sky_day_shares(*sky)
    if day is not None:
        out["day_share_equator"] = equator_day_share(tilt, alpha) if sky is None else equator_share
        out["solar_day_hours"] = day
    if parent is not None:
        mutual = spin.obliquity_between(_orbit_normal(body, parent), _orbit_normal(parent))
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
