# planetgen/describe/clock_rows.py

"""
The derived rows of the object sheet for a body's clock, the eclipses of its star and its gravitational lens (GEN.1434)
=======================================================================================================================

`rows(conn, kind, object_id)` reads the stored star, planet and moon and gives the rows `physics/sky_clock` works out: how far a clock on the body falls behind a distant one each year (a star, a
planet, a moon), how often and how long the moons and the inner planets eclipse its star for a planet, and where a star's gravity brings the light that grazes it to a focus. They are worked out when
the sheet is made and are not columns: each is a pure function of stored data, the way the moon's longest eclipse in its planet's shadow is. Each row is `(group, key, label, value, raw)`.
"""

import math

from planetgen.describe import day_rows, heavy_rows, host_star, model, sky_size_rows
from planetgen.physics import orbit_elements, sky_clock, sky_suns, sky_view, spin, vec3

CLOCK_KEY = "clock_lag_s_per_yr"
"""str: The key of the derived row for the seconds a year a clock on the body falls behind a distant clock."""

LENS_KEY = "lens_focus_au"
"""str: The key of the star's derived row for the distance of the focus of its gravitational lens."""

ECLIPSE_KEY_PREFIX = "star_eclipses_"
"""str: The start of the key of the planet's derived rows for the eclipses of its star, one per blocker, numbered from 1 in the order of the shortest interval."""

SUN_ECLIPSE_KEY_PREFIX = "sun_eclipses_"
"""str: The start of the key of a planet's derived rows for the eclipses of its other suns by its moons (GEN.1443), numbered from 1 in the order of the shortest interval."""

MAX_ECLIPSE_ROWS = 4
"""int: The most blockers a planet's sheet lists; the ones that eclipse its star most often."""

def is_derived_key(key):
    """Whether `key` is one of the rows this module gives (the clock, the lens, or an eclipse row)."""
    return key in (CLOCK_KEY, LENS_KEY) or key.startswith((ECLIPSE_KEY_PREFIX, SUN_ECLIPSE_KEY_PREFIX)) or heavy_rows.is_derived_key(key) or day_rows.is_derived_key(key)


def _star_of(conn, planet):
    return host_star.of(conn, planet, "mass_kg, radius_km")


def _clock_row(own, hosts):
    lag = sky_clock.lag_per_year_s(sky_clock.clock_lag_fraction(own["mass_kg"], own["radius_km"], hosts))
    if not lag:
        return None
    return ("physical", CLOCK_KEY, "Clock lag behind a distant clock, a year", [model.qty("time", lag, "s")], lag)


def _clock_rows(conn, kind, object_id):
    if kind == "star":
        star = conn.execute("SELECT mass_kg, radius_km FROM stars WHERE id = ?", (object_id,)).fetchone()
        row = _clock_row(star, ()) if star is not None else None
    elif kind == "planet":
        planet = conn.execute("SELECT star_system_id, star_id, mass_kg, radius_km, distance_km FROM planets WHERE id = ?", (object_id,)).fetchone()
        star = _star_of(conn, planet) if planet is not None else None
        row = _clock_row(planet, [(star["mass_kg"], planet["distance_km"])]) if star is not None and planet["distance_km"] else None
    elif kind == "moon":
        moon = conn.execute("SELECT parent_type, parent_id, mass_kg, radius_km, distance_km FROM moons WHERE id = ?", (object_id,)).fetchone()
        if moon is None or moon["parent_type"] != "planet" or not moon["distance_km"]:
            return []
        planet = conn.execute("SELECT star_system_id, star_id, mass_kg, distance_km FROM planets WHERE id = ?", (moon["parent_id"],)).fetchone()
        if planet is None:
            return []
        star = _star_of(conn, planet)
        hosts = [(planet["mass_kg"], moon["distance_km"])] + ([(star["mass_kg"], planet["distance_km"])] if star is not None and planet["distance_km"] else [])
        row = _clock_row(moon, hosts)
    else:
        return []
    return [row] if row else []


def _lens_rows(conn, kind, object_id):
    if kind != "star":
        return []
    star = conn.execute("SELECT mass_kg, radius_km FROM stars WHERE id = ?", (object_id,)).fetchone()
    focus = sky_clock.focal_distance_m(star["mass_kg"], star["radius_km"]) if star is not None else None
    if not focus:
        return []
    return [("physical", LENS_KEY, "Focus of its gravitational lens", [model.qty("length", focus, "m")], focus)]


def _eclipse_value(eclipse):
    """The inline text of one blocker's eclipses: how often, how long and how deep."""
    items = ["about every ", model.qty("time", eclipse.interval_yr, "yr"), ", up to ", model.qty("time", eclipse.longest_contact_h, "h"), " from first contact to last"]
    if eclipse.longest_total_h > 0.0:
        items += [", the whole star hidden for up to ", model.qty("time", eclipse.longest_total_h, "h")]
    else:
        items += [", at most ", model.qty("fraction", eclipse.deepest_cover, "pct"), " of it hidden"]
    return items


def _eclipse_rows(conn, kind, object_id):
    if kind != "planet":
        return []
    planet = conn.execute("SELECT id, star_system_id, star_id, distance_km, period_years FROM planets WHERE id = ?", (object_id,)).fetchone()
    star = _star_of(conn, planet) if planet is not None else None
    if star is None or not planet["distance_km"] or not planet["period_years"] or not star["radius_km"]:
        return []
    found = []
    for moon in conn.execute("SELECT name, radius_km, distance_km, period_years, orbital_inclination_deg FROM moons WHERE parent_type = 'planet' AND parent_id = ?", (object_id,)).fetchall():
        eclipse = sky_clock.moon_eclipse(star["radius_km"], planet["distance_km"], planet["period_years"], moon["radius_km"], moon["distance_km"], moon["period_years"],
                                         moon["orbital_inclination_deg"] or 0.0)
        if eclipse:
            found.append((moon["name"], eclipse))
    inner = conn.execute("SELECT name, radius_km, distance_km, period_years, orbital_inclination_deg FROM planets WHERE id <> ? AND distance_km < ? AND "
                         + ("star_id = ?" if planet["star_id"] else "star_system_id = ?"), (object_id, planet["distance_km"], planet["star_id"] or planet["star_system_id"])).fetchall()
    for other in inner:
        eclipse = sky_clock.planet_eclipse(star["radius_km"], planet["distance_km"], planet["period_years"], other["radius_km"], other["distance_km"], other["period_years"],
                                           other["orbital_inclination_deg"] or 0.0)
        if eclipse:
            found.append((other["name"], eclipse))
    found.sort(key=lambda pair: pair[1].interval_yr)
    return [("orbit", f"{ECLIPSE_KEY_PREFIX}{number}", f"Eclipses of its star by {name}", _eclipse_value(eclipse), eclipse.interval_yr)
            for number, (name, eclipse) in enumerate(found[:MAX_ECLIPSE_ROWS], 1)]


def _planet_eclipse_row(conn, kind, object_id):
    """A moon's sheet: the eclipses of its star by its planet, the planet's disk against the star's from the moon (the moon is the observer, the planet the blocker)."""
    if kind != "moon":
        return []
    moon = conn.execute("SELECT parent_type, parent_id, distance_km, period_years, orbital_inclination_deg FROM moons WHERE id = ?", (object_id,)).fetchone()
    if moon is None or moon["parent_type"] != "planet":
        return []
    planet = conn.execute("SELECT star_system_id, star_id, radius_km, distance_km, period_years FROM planets WHERE id = ?", (moon["parent_id"],)).fetchone()
    star = _star_of(conn, planet) if planet is not None else None
    if star is None or not star["radius_km"] or not planet["distance_km"] or not planet["radius_km"]:
        return []
    eclipse = sky_clock.moon_eclipse(star["radius_km"], planet["distance_km"], planet["period_years"], planet["radius_km"], moon["distance_km"], moon["period_years"],
                                     moon["orbital_inclination_deg"] or 0.0)
    if eclipse is None:
        return []
    return [("orbit", f"{ECLIPSE_KEY_PREFIX}1", "Eclipses of its star by its planet", _eclipse_value(eclipse), eclipse.interval_yr)]


def _sun_eclipse_rows(conn, kind, object_id):
    """A planet's sheet: the eclipses of each other sun of its sky that shows (`sky_view.bodies`) by each of its moons (`sky_clock.companion_eclipse`). Loads the system, so only a planet that has moons and stored other suns is asked."""
    if kind != "planet":
        return []
    planet = conn.execute("SELECT star_system_id, sky_suns FROM planets WHERE id = ?", (object_id,)).fetchone()
    if planet is None or not planet["sky_suns"]:
        return []
    moons = conn.execute("SELECT name, radius_km, distance_km, period_years, orbital_inclination_deg FROM moons WHERE parent_type = 'planet' AND parent_id = ?", (object_id,)).fetchall()
    if not moons:
        return []
    from planetgen.db import store
    system = store.load_star_system(conn, planet["star_system_id"])
    host = next((candidate for candidate in system.all_planets if getattr(candidate, "db_id", None) == object_id and getattr(candidate, "body_type", None) != "a"), None)
    star = getattr(host, "star", None)
    if host is None or star is None or not getattr(host, "distance", None) or host.spin_axis_x is None:
        return []
    light_au = orbit_elements.flux_distance_au(host.distance, getattr(host, "eccentricity", None) or 0.0)
    sky = sky_view.bodies(star, light_au)[1:sky_size_rows.MAX_OTHER_ROWS + 1]
    normal = spin.orbit_normal(host.orbital_inclination_deg, host.orbital_ascending_node_deg)
    _, _, placed = sky_suns.placed(star, (host.spin_axis_x, host.spin_axis_y, host.spin_axis_z), normal, light_au)
    found = []
    for number, body in enumerate(sky, 1):
        source = next((item for item in placed if item.star is body.star), None)
        if source is None:
            continue
        latitude = math.asin(max(-1.0, min(1.0, vec3.unit(source.offset_au)[2])))       # off the planet's orbit plane, the moons' reference plane
        for moon in moons:
            eclipse = sky_clock.companion_eclipse(math.radians(body.angular_deg) / 2.0, latitude, moon["radius_km"], moon["distance_km"], moon["period_years"], moon["orbital_inclination_deg"] or 0.0)
            if eclipse:
                found.append((f"Eclipses of its {sky_size_rows.source_phrase(number, body)} by {moon['name']}", eclipse))
    found.sort(key=lambda pair: pair[1].interval_yr)
    return [("orbit", f"{SUN_ECLIPSE_KEY_PREFIX}{number}", label, _eclipse_value(eclipse), eclipse.interval_yr) for number, (label, eclipse) in enumerate(found[:MAX_ECLIPSE_ROWS], 1)]


def rows(conn, kind, object_id):
    """The derived rows of the sheet of a star, planet or moon for its clock, its eclipses and its lens, `(group, key, label, value, raw)` each; none for any other kind."""
    return (_clock_rows(conn, kind, object_id) + _lens_rows(conn, kind, object_id) + _eclipse_rows(conn, kind, object_id)
            + _planet_eclipse_row(conn, kind, object_id) + _sun_eclipse_rows(conn, kind, object_id))
