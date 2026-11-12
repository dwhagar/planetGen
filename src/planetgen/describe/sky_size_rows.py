# planetgen/describe/sky_size_rows.py

"""
The derived rows of the object sheet for how wide the sources look in a body's sky (GEN.1434)
===========================================================================================

`rows(conn, kind, object_id)` reads the stored planet or moon and gives how wide its sun looks (a planet's star, a moon's too), how wide a moon's planet looks from the moon, and how wide each
other sun of a multiple system looks. They are worked out when the sheet is made and are not columns: each is a pure function of stored data, the way the clock rows of `clock_rows` are. Each
row is `(group, key, label, value, raw)`, the `raw` the angular diameter in degrees.

**How wide.** The angular diameter of a disk of radius `R` at distance `d` is `2 asin(R / d)` (`daylight.angular_radius`, the figure the sky text and the day-length laws already use). A planet's
star stands at the circle that gives the planet the same mean light as its orbit (`orbit_elements.flux_distance_au`, the root-mean-square distance, as `sky_view.bodies` takes it); a moon's star
is at its planet's. A moon's planet stands at the moon's orbit, `2 asin(R_planet / a_moon)`: the Moon sees the Earth 3.7 times as wide as the Earth sees the Moon.

**The other suns.** A body of a multiple system has a sky that holds every other source bright enough to show (`sky_view.bodies`, GEN.1443), each with the width its star has at its own distance.
They need the star objects of the whole system (the companions' radii and the orbit that gives their light), so they are read from the loaded system only for a body that has stored other suns
(`sky_suns` not empty), and a lone star's body never loads one. Up to `MAX_OTHER_ROWS` rows, the brightest first.

**What stays out.** A rogue planet has no sun and no row. A body's sibling planets and moons that stand in its sky as disks (not as suns) are not rows: their widths change by the hour, and the
eclipse and transit rows of the clock block (`clock_rows`) already give how often one crosses its sun and how long it takes.

Stand-ins are the rows SZ1 and SZ2 of section 6.280 of the gap plan.
"""

import math

from planetgen.physics import constants, daylight, orbit_elements, sky_view
from planetgen.util.format import three_figures_or_tiny

SUN_KEY = "sky_width_sun"
"""str: The key of the derived row for how wide a body's sun looks."""

PLANET_KEY = "sky_width_planet"
"""str: The key of a moon's derived row for how wide its planet looks."""

OTHER_KEY_PREFIX = "sky_width_other_"
"""str: The start of the key of a body's derived rows for its other suns, numbered from 1 by their light."""

MAX_OTHER_ROWS = 4
"""int: The most other suns a sheet lists; the brightest."""

_ORDINALS = ("second", "third", "fourth", "fifth", "sixth")


def is_derived_key(key):
    """Whether `key` is one of the rows this module gives."""
    return key in (SUN_KEY, PLANET_KEY) or key.startswith(OTHER_KEY_PREFIX)


def width_deg(radius_km, distance_km):
    """The angular diameter of a disk of `radius_km` at `distance_km`, degrees; 0 for no size or no distance."""
    if not radius_km or not distance_km or distance_km <= 0.0:
        return 0.0
    return math.degrees(2.0 * daylight.angular_radius(radius_km, distance_km / constants.AU_TO_KM))


def _row(key, label, degrees):
    return ("orbit", key, label, [f"{three_figures_or_tiny(degrees)} degrees"], degrees)


def _host_star(conn, planet):
    if planet["star_id"]:
        star = conn.execute("SELECT radius_km FROM stars WHERE id = ?", (planet["star_id"],)).fetchone()
        if star is not None:
            return star
    return conn.execute("SELECT radius_km FROM stars WHERE star_system_id = ? ORDER BY id LIMIT 1", (planet["star_system_id"],)).fetchone()


def _sun_row(conn, planet):
    """The row for the width of the star a planet (or the planet a moon circles) has in its sky."""
    star = _host_star(conn, planet)
    if star is None or not planet["distance_km"]:
        return []
    distance_km = orbit_elements.flux_distance_au(planet["distance_km"] / constants.AU_TO_KM, planet["eccentricity"] or 0.0) * constants.AU_TO_KM
    degrees = width_deg(star["radius_km"], distance_km)
    return [_row(SUN_KEY, "Width of its sun in the sky", degrees)] if degrees else []


def _other_rows(conn, kind, object_id, planet, parent):
    """The rows for the other suns of a multiple system's body: `parent` the planet a moon circles. Loads the system, so only a body with stored other suns is asked."""
    from planetgen.db import store
    system = store.load_star_system(conn, planet["star_system_id"])
    host = None
    for candidate in system.all_planets:
        if getattr(candidate, "db_id", None) == (parent["id"] if parent is not None else object_id) and getattr(candidate, "body_type", None) != "a":
            host = candidate
            break
    star = getattr(host, "star", None)
    if host is None or star is None or not getattr(host, "distance", None):
        return []
    sky = sky_view.bodies(star, orbit_elements.flux_distance_au(host.distance, getattr(host, "eccentricity", None) or 0.0))
    rows = []
    for number, body in enumerate(sky[1:MAX_OTHER_ROWS + 1], 1):
        word = "sun" if body.role == "sun" else "night light"
        rows.append(_row(f"{OTHER_KEY_PREFIX}{number}", f"Width of its {_ORDINALS[min(number - 1, len(_ORDINALS) - 1)]} {word} in the sky ({body.kind})", body.angular_deg))
    return rows


def rows(conn, kind, object_id):
    """The derived rows of the sheet of a planet or moon for how wide its sources look, `(group, key, label, value, raw)` each; none for any other kind and none for a rogue planet."""
    if kind == "planet":
        planet = conn.execute("SELECT star_system_id, star_id, distance_km, eccentricity, sky_suns FROM planets WHERE id = ?", (object_id,)).fetchone()
        if planet is None:
            return []
        found = _sun_row(conn, planet)
        return found + (_other_rows(conn, kind, object_id, planet, None) if found and planet["sky_suns"] else [])
    if kind == "moon":
        moon = conn.execute("SELECT parent_type, parent_id, distance_km, sky_suns FROM moons WHERE id = ?", (object_id,)).fetchone()
        if moon is None or moon["parent_type"] != "planet":
            return []
        planet = conn.execute("SELECT id, star_system_id, star_id, distance_km, eccentricity, radius_km FROM planets WHERE id = ?", (moon["parent_id"],)).fetchone()
        if planet is None:
            return []
        found = _sun_row(conn, planet)
        degrees = width_deg(planet["radius_km"], moon["distance_km"])
        if degrees:
            found = found + [_row(PLANET_KEY, "Width of its planet in the sky", degrees)]
        return found + (_other_rows(conn, kind, object_id, planet, planet) if found and moon["sky_suns"] else [])
    return []
