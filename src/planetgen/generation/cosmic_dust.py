# planetgen/generation/cosmic_dust.py

"""
The cosmic dust a planet takes in from its system, given to the V8 layer of its sulfur cloud (GEN.1542 slice 1)
===============================================================================================================

V8 of the Venus sulfur-cloud set (`sulfur_cloud_set`) needs carbon from cosmic dust. The dust a rocky planet takes in each year is what the belts and comets of its own system send it
(`influx.planet_influx_kg_per_yr`, GEN.1543), and `dust_supply.v8_supplied` says whether that is enough carbon for the layer. This pass reads each rocky planet's supply once the system's belts are final
(after `late_accretion.apply`) and makes its cloud chemistry again with that answer (`sulfur_cloud_set.fill(..., carbon_dust=)`).

V8 is a single-source hypothesis and off by default (`sulfur_cloud_set.INCLUDE_HYPOTHESES`), and the supply decides nothing else, so the pass does nothing while the hypotheses are off: the stored cloud chemistry
of every body is what the first fill gave it.
"""

from planetgen.generation import late_accretion
from planetgen.physics import constants, dust_supply, habitability_world, influx, sulfur_cloud_set
from planetgen.physics import unit_registry as units

_kg_to_mearth = units.converter("kg", "Mearth")
_kg_to_msun = units.converter("kg", "Msun")
_km_to_rearth = units.converter("km", "Rearth")


def reservoirs_of(planet, planets, comets):
    """The belts of `planet`'s star among `planets` (as they stand now, their present mass and age) and its comets (`comets`: `(star, [Comet])` pairs, one for each star that has some), as `influx.Reservoir`."""
    star = planet.star
    age_yr = star.age * 1.0e9
    found = [influx.Reservoir("belt", belt.mass_kg, belt.lower_limit, belt.upper_limit, age_yr)
             for belt in planets if belt.body_type == "a" and late_accretion._same_star(belt, planet) and belt.mass_kg]
    found += [influx.Reservoir("comet", comet.mass_kg, comet.perihelion_distance_au, comet.perihelion_distance_au)
              for host, owned in comets if host is star for comet in owned if comet.mass_kg and comet.perihelion_distance_au]
    return found


def influx_of(planet, planets, comets):
    """The dust, kg a year, the planet takes in from the belts and comets of its star (`influx.planet_influx_kg_per_yr`), the giants of the star in the way."""
    star = planet.star
    giants = [influx.Giant(_kg_to_mearth(body.mass), _km_to_rearth(body.radius), body.distance)
              for body in planets if body.body_type == "g" and body.star is star and body.mass and body.radius]
    reservoirs = reservoirs_of(planet, planets, comets)
    if not reservoirs:
        return 0.0
    return influx.planet_influx_kg_per_yr(planet.radius, planet.mass, planet.distance, _kg_to_msun(star.mass),
                                          star.luminosity / constants.SOLAR_LUMINOSITY, reservoirs, giants)[0]


def apply(planets, comets):
    """
    Makes the cloud chemistry again, with the carbon its dust supplies, of every rocky planet of `planets` (`StarSystem.all_planets`) that has any (a body whose first fill gave it none has no acid cloud for V8 to
    live in). Does nothing while `sulfur_cloud_set.INCLUDE_HYPOTHESES` is off.

    Returns:
        int: The bodies made again.
    """
    if not sulfur_cloud_set.INCLUDE_HYPOTHESES:
        return 0
    made = 0
    for planet in planets:
        if planet.body_type == "a" or not late_accretion.is_rocky_of_a_star(planet) or getattr(planet, "cloud_chemistry", None) is None or not planet.star.luminosity:
            continue
        values = habitability_world.values_from_body(planet)
        values["energy_flux_w_m2"] = getattr(planet, "energy_flux_w_m2", None)
        sulfur_cloud_set.fill(planet, values, carbon_dust=dust_supply.v8_supplied(influx_of(planet, planets, comets), planet.radius))
        made += 1
    return made
