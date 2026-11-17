# planetgen/generation/habitable_request.py

"""
What a request for a habitable world asks of a finished body (GEN.911, GEN.916, GEN.1130)
==========================================================================================

A request for a habitable world (`system_config.HABITABLE_WORLD`) asks for conditions, never for a class letter: a planet or moon is
what was asked for when its solved surface is temperate (`tuning.REQUEST_TEMPERATURE_K`), its air between the pressures of `tuning.REQUEST_PRESSURE_PA`
(50 to 250 kPa), its water a surface ocean (`tuning.REQUEST_HYDROSPHERES`) and that ocean covers `tuning.REQUEST_WATER_COVER` of the surface (20 to 95 percent,
so there is land and sea; GEN.1702). Its class is derived from the same conditions after it is built
(GEN.646). Life is no more a matter of the letter: the life classes score the finished body's conditions (`physics.life_classes`), and a system that bars
habitable worlds keeps out a body that meets the conditions or has life on it (`carries_life`, GEN.1387).

The request only raises the odds (GEN.1130, `util.request_odds`): the star, the planet in the ecosphere and the system's planets are
drawn as they would be, and a draw without the conditions is kept one time in `tuning.REQUEST_ODDS` (the host star, `physics.stellar_evolution`;
the planet, `StarSystem._planet_in_ecosphere`; a system with no rocky planet in the ecosphere, `StarSystem._habitable_request_stands`). Nothing is placed, aimed or adapted, and
a system the request did not reach says `"unmet"` (`StarSystem.habitable_request_status`).

The request is reached when a body that meets the conditions also has life on it (`count_living`, GEN.1778): a temperate ocean world whose life never began (a star too young for the
first step of any class, or no class its air supports) is not a habitable world, and a system whose only such body is lifeless says `"unmet"` instead of `"aim"`, so a finished
system holds a habitable world (`StarSystem.habitable_class_count`) or says it did not.
"""

from planetgen import tuning


def meets(body):
    """Whether the finished planet or moon `body` has the conditions a request for a habitable world asks for."""
    temperature = getattr(body, "surface_temperature", None)
    pressure = getattr(body, "atmospheric_pressure", None)
    cover = getattr(body, "ocean_fraction", None)
    low, high = tuning.REQUEST_TEMPERATURE_K
    low_pa, high_pa = tuning.REQUEST_PRESSURE_PA
    low_cover, high_cover = tuning.REQUEST_WATER_COVER
    return (getattr(body, "body_type", None) == "t" and temperature is not None and low <= temperature <= high
            and pressure is not None and low_pa <= pressure <= high_pa
            and getattr(body, "hydrosphere", None) in tuning.REQUEST_HYDROSPHERES
            and cover is not None and low_cover <= cover <= high_cover)


def carries_life(body):
    """Whether life lives on the finished planet or moon `body`: the life pass gave it a life class (`planet.life_class`). Life comes from the conditions and the
    history of the body, never from its class letter (GEN.1387); a body that is not through the life pass has none."""
    return getattr(body, "life_class", None) is not None


def count(planets):
    """The planets and moons of `planets` (belts skipped) that meet the request (`meets`)."""
    return sum(1 for planet in planets if planet.body_type != "a" for body in (planet, *planet.moons) if meets(body))


def count_living(planets):
    """The planets and moons of `planets` (belts skipped) that meet the request (`meets`) and have life on them (`carries_life`): the habitable worlds a request for one is
    answered with (GEN.1778). Meaningful only once the life pass has run on every body."""
    return sum(1 for planet in planets if planet.body_type != "a" for body in (planet, *planet.moons) if meets(body) and carries_life(body))

