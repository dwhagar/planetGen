# planetgen/generation/habitable_request.py

"""
What a request for a habitable world asks of a finished body (GEN.911, GEN.916, GEN.1130)
==========================================================================================

A request for a habitable world (`system_config.HABITABLE_WORLD`) asks for conditions, never for a class letter: a planet or moon is
what was asked for when its solved surface is temperate (`tuning.REQUEST_TEMPERATURE_K`), its air at least `tuning.REQUEST_MIN_PRESSURE_PA`
thick and its water a surface ocean (`tuning.REQUEST_HYDROSPHERES`). Its class is derived from the same conditions after it is built
(GEN.646), and the classes that "can carry life" (`tuning.HABITABLE_PLANET_CLASSES`) include a hot reducing world, a glaciated one and a
carbon dioxide one, so the letter is no test of them.

The request only raises the odds (GEN.1130, `util.request_odds`): the star, the planet in the ecosphere and the system's planets are
drawn as they would be, and a draw without the conditions is kept one time in `tuning.REQUEST_ODDS` (the host star, `physics.stellar_evolution`;
the planet, `StarSystem._planet_in_ecosphere`; a system with no rocky planet in the ecosphere, `StarSystem._habitable_request_stands`). Nothing is placed, aimed or adapted, and
a system the request did not reach says `"unmet"` (`StarSystem.habitable_request_status`).
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


def count(planets):
    """The planets and moons of `planets` (belts skipped) that meet the request (`meets`)."""
    return sum(1 for planet in planets if planet.body_type != "a" for body in (planet, *planet.moons) if meets(body))

