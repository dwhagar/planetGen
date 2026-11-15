# planetgen/generation/dust_supply.py

"""
The cosmic dust a rocky body takes in from its own system (GEN.1543, rewiring GEN.1542)
=========================================================================================

`sulfur_cloud_set` gives a body the iron-bearing carbon species V8 (hypothesis, off by default) only where cosmic dust reaches it. That dust is not a property of the body alone: it is what the belts and comets of
its star send inward (`physics/influx.py`: grain crossings, the share the body takes, what a giant in the way ejects). `apply` works it out once for every rocky planet and moon of a finished system, the belts
and comets all in place, and keeps it on the body as `_dust_influx_kg_yr` (not stored: a pure function of the belts, comets and giants that are). `sulfur_cloud_set.fill` reads it when the life pass scores the
body, so the dust the V8 gate sees is the system's own. A moon takes the dust at its planet's distance, with its own radius.

It is one call (`StarSystem.generate`, before the late accretion pass), so a reorder of the passes moves one line.
"""

from planetgen.generation import late_accretion
from planetgen.physics import constants, influx

ATTRIBUTE = "_dust_influx_kg_yr"
"""str: The attribute `apply` sets on each rocky planet and moon (kg a year); a body it never reached has none, and the V8 gate reads that as no dust."""

EARTH_RADII_KM = 6371.0


def reservoirs_of(belts, comets, star_age_gy):
    """The `influx.Reservoir`s of the belts (`AsteroidBelt`s with mass and ring) and comets (nucleus mass and perihelion) a star has, the belts aged at the star's age."""
    found = [influx.Reservoir("belt", belt.mass_kg, belt.lower_limit, belt.upper_limit, star_age_gy * 1.0e9)
             for belt in belts if getattr(belt, "mass_kg", None) and belt.lower_limit and belt.upper_limit]
    found += [influx.Reservoir("comet", comet.mass_kg, comet.perihelion_distance_au, comet.perihelion_distance_au)
              for comet in comets if getattr(comet, "mass_kg", None) and getattr(comet, "perihelion_distance_au", None)]
    return found


def giants_of(bodies, exclude=None):
    """The `influx.Giant`s among `bodies` (planets, not moons), leaving `exclude` out."""
    return [influx.Giant(body.mass / constants.EARTH_MASS_TO_KG, body.radius / EARTH_RADII_KM, body.distance)
            for body in bodies if body.body_type == "g" and body is not exclude and not getattr(body, "is_moon", False) and body.mass and body.radius and body.distance]


def influx_of(body, planet, siblings, comets):
    """The dust, kg a year, that reaches `body` (a planet, or a moon of `planet`) from the belts of `siblings` and the `comets` of its star. 0 where there are no reservoirs."""
    star = planet.star
    star_mass_solar = star.mass / constants.SOLAR_MASS_TO_KG
    belts = [b for b in siblings if b.body_type == "a"]
    reservoirs = reservoirs_of(belts, comets, star.age)
    if not reservoirs or not body.radius or not planet.distance:
        return 0.0
    luminosity = getattr(star, "luminosity", None)
    return influx.planet_influx_kg_per_yr(body.radius, planet.distance, star_mass_solar, reservoirs, giants_of(siblings, exclude=planet),
                                          luminosity / constants.SOLAR_LUMINOSITY if luminosity else None)[0]


def apply(system):
    """
    Gives every rocky planet and moon of `system` (a finished `StarSystem`, belts and comets in place) its `ATTRIBUTE`. A giant and a belt are left alone. Comets are the system's primary and, in a wide pair, secondary
    star's (`system.comets`, `system.secondary_comets`); a star of the tree beyond them has none.
    Returns how many bodies were given dust.
    """
    planets = system.all_planets
    comets_by_star = {id(system.star): system.comets}
    if getattr(system, "binary_type", None) == "wide":
        comets_by_star[id(system.secondary_star)] = system.secondary_comets
    given = 0
    for planet in planets:
        if planet.body_type in ("a", "g") or getattr(planet, "star", None) is None or getattr(planet.star, "starless", False):
            continue
        siblings = late_accretion.siblings_of(planet, planets)
        comets = comets_by_star.get(id(planet.star), [])
        for body in (planet, *planet.moons):
            if body.body_type in ("a", "g"):
                continue
            setattr(body, ATTRIBUTE, influx_of(body, planet, siblings, comets))
            given += 1
    return given
