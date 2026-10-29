"""Planets of a named class for tests (GEN.646): the class is derived from conditions, so a test that needs a class draws bodies until one
comes out as it, from streams of its own."""

from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.util import draw

DEFAULT_TRIES = 300


def draw_planet(star, distance, seed=0, config=None, **kwargs):
    """One `Planet` at `distance`, drawn from the stream `seed`; `kwargs` are `Planet`'s own."""
    with draw.bound(f"class-support:{seed}"):
        return Planet(config or SystemConfig(), star, star.habitable_zone, distance, **kwargs)


def planet_of_class(planet_class, star, distance, config=None, tries=DEFAULT_TRIES, seed=0, **kwargs):
    """
    The first planet at `distance` whose derived class is `planet_class`, drawing from the streams `seed`, `seed + 1`, ...

    Raises:
        LookupError: When none of `tries` draws is that class (the class does not occur there; pick a place where it does).
    """
    for index in range(seed, seed + tries):
        planet = draw_planet(star, distance, index, config, **kwargs)
        if planet.planet_class == planet_class:
            return planet
    raise LookupError(f"no class {planet_class} planet in {tries} draws at {distance} AU")


def planets_where(star, distance, keep, count, config=None, tries=DEFAULT_TRIES, seed=0, **kwargs):
    """`count` planets at `distance` for which `keep(planet)` holds, drawn from successive streams."""
    found = []
    for index in range(seed, seed + tries):
        planet = draw_planet(star, distance, index, config, **kwargs)
        if keep(planet):
            found.append(planet)
            if len(found) == count:
                return found
    raise LookupError(f"only {len(found)} of {count} wanted planets in {tries} draws at {distance} AU")
