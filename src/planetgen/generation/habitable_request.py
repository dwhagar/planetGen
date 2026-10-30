# planetgen/generation/habitable_request.py

"""
What a request for a habitable world asks of a finished body (GEN.911, GEN.916)
=================================================================================

A request for a habitable world (`system_config.HABITABLE_WORLD`) asks for conditions, never for a class letter: a planet or moon is
what was asked for when its solved surface is temperate (`tuning.REQUEST_TEMPERATURE_K`), its air at least `tuning.REQUEST_MIN_PRESSURE_PA`
thick and its water a surface ocean (`tuning.REQUEST_HYDROSPHERES`). Its class is derived from the same conditions after it is built
(GEN.646), and the classes that "can carry life" (`tuning.HABITABLE_PLANET_CLASSES`) include a hot reducing world, a glaciated one and a
carbon dioxide one, so the letter is no test of them.

The orbit is aimed by the climate the draw solved (GEN.911). The zone's outer edge is where about 8 bar of carbon dioxide would hold 273 K,
and the carbonate-silicate thermostat leaves far less there, so a world drawn at a distance in the zone is most often frozen. A draw that
misses by its temperature moves the next draw to where the star's light would give `tuning.REQUEST_TARGET_K` (the surface temperature of a
world with the same air goes as the flux to the quarter, so the distance goes as the square of the ratio of temperatures), kept inside the
zone. The solve decides what each draw is; nothing is forced.
"""

from planetgen import tuning


def meets(body):
    """Whether the finished planet or moon `body` has the conditions a request for a habitable world asks for."""
    temperature = getattr(body, "surface_temperature", None)
    pressure = getattr(body, "atmospheric_pressure", None)
    low, high = tuning.REQUEST_TEMPERATURE_K
    return (getattr(body, "body_type", None) == "t" and temperature is not None and low <= temperature <= high
            and pressure is not None and pressure >= tuning.REQUEST_MIN_PRESSURE_PA
            and getattr(body, "hydrosphere", None) in tuning.REQUEST_HYDROSPHERES)


def meets_with_moons(planet):
    """Whether `planet` or one of its moons meets the request (`meets`)."""
    return meets(planet) or any(meets(moon) for moon in getattr(planet, "moons", ()))


def count(planets):
    """The planets and moons of `planets` (belts skipped) that meet the request (`meets`)."""
    return sum(1 for planet in planets if planet.body_type != "a" for body in (planet, *planet.moons) if meets(body))


def aimed_distance(distance_au, temperature_k, inner_au, outer_au):
    """
    Where the next draw of a request goes, in AU: the distance at which the light would give a world of this air the target temperature
    (`distance * (T / target)^2`), kept inside `[inner_au, outer_au]`. A draw with no solved temperature keeps its distance.
    """
    if not temperature_k or temperature_k <= 0.0:
        return distance_au
    aimed = distance_au * (temperature_k / tuning.REQUEST_TARGET_K) ** 2
    return min(max(aimed, inner_au), outer_au)
