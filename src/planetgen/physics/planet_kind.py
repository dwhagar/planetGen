# planetgen/physics/planet_kind.py

"""
Planet kind from formation (GEN.681)
====================================

What sort of body the disc built at a formation distance, without a class letter. The mass law (`planet_mass`) gives the
mass and whether the core ran away into a gas giant; this module adds the kind the rest of the generator branches on and
the one fact about origin the composition draws need.

* **giant**: a runaway gas giant, a core of `GIANT_CORE_CRITICAL_MASS_EARTH` or more that formed beyond the snow line
  and kept no runaway envelope (an ice giant: ice-rich, with a hydrogen and helium envelope of a few percent), or any body
  heavier than the heaviest terrestrial one (20 Earth masses). Its radius and envelope come from the giant draws.
* **terrestrial**: every other core, rock-and-metal inside the snow line, ice-and-rock beyond it. Whether a heavy inner
  terrestrial body holds a hydrogen envelope (a gas dwarf) is the envelope law's call (GEN.822), not decided here.

`beyond_snow_line` is the origin flag: bodies formed beyond the star's snow line carry their water as ice (the water
draws read it), whatever they later become. A pure function of the star and the distance; it reads no class and writes none.
"""

import collections

from planetgen import tuning
from planetgen.physics import constants, planet_mass
from planetgen.util.random import sample_power_law

Formation = collections.namedtuple("Formation", ("kind", "mass_kg", "core_kg", "runaway", "beyond_snow_line"))
"""`kind` "terrestrial" or "giant" (`atmosphere_recipes.BODY_KINDS`); `mass_kg` the planet's mass; `core_kg` the assembled
solid mass it grew from; `runaway` whether it became a gas giant; `beyond_snow_line` whether it formed past the snow line."""

_EARTH_KG = constants.EARTH_MASS_TO_KG
_MAX_TERRESTRIAL_EARTH = 20.0


def draw_formation(distance_au, star_mass_kg, snow_line_au, fe_h=0.0, metal_ratio=1.0, stream=None):
    """
    The kind, mass and origin of the planet formed at `distance_au`. Draws exactly what `planet_mass.planet_mass` draws,
    so the same stream gives the same planet.

    Args:
        distance_au (float): The formation distance from the star, in AU.
        star_mass_kg (float): The star's mass, in kilograms.
        snow_line_au (float): The star's own snow line, in AU.
        fe_h (float): The star's [Fe/H] in dex.
        metal_ratio (float): The star's Z / Z0 (`bulk_draw.metal_ratio_to_solar`, GEN.464): the disc's solids scale with it.
        stream (draw.Stream or None): The body's own stream; `None` draws from the stream bound now.

    Returns:
        Formation: The kind, mass, core, runaway flag and origin.
    """
    mass = planet_mass.planet_mass(distance_au, star_mass_kg, snow_line_au, fe_h, metal_ratio, stream)
    beyond = distance_au >= snow_line_au
    critical = tuning.GIANT_CORE_CRITICAL_MASS_EARTH * _EARTH_KG
    giant = mass.runaway or mass.mass_kg > _MAX_TERRESTRIAL_EARTH * _EARTH_KG or (beyond and mass.core_kg >= critical)
    return Formation("giant" if giant else "terrestrial", mass.mass_kg, mass.core_kg, mass.runaway, beyond)


def kind_of_given(mass_kg=None, radius_km=None, beyond_snow_line=False):
    """
    The kind of a body whose caller gave its mass or its radius (the formation law is not asked): a mass above the heaviest
    terrestrial body (`bulk_draw.MAX_MASS_EARTH`, 20 Earth masses), or a core-critical mass beyond the snow line, is a giant; a radius
    alone names a giant from `GIANT_MIN_RADIUS_KM`; anything else is terrestrial. The mass wins when both are given.

    Returns:
        str: "terrestrial" or "giant".
    """
    if mass_kg is not None:
        mass_earth = mass_kg / _EARTH_KG
        if mass_earth > _MAX_TERRESTRIAL_EARTH:
            return "giant"
        return "giant" if beyond_snow_line and mass_earth >= tuning.GIANT_CORE_CRITICAL_MASS_EARTH else "terrestrial"
    if radius_km is not None and radius_km >= tuning.GIANT_MIN_RADIUS_KM:
        return "giant"
    return "terrestrial"


def moon_mass_ceiling_kg(parent_mass_kg, parent_is_giant):
    """The heaviest moon a parent holds: its share of the parent's mass (`MOON_GIANT_PARENT_MAX_MASS_RATIO` or
    `MOON_TERRESTRIAL_PARENT_MAX_MASS_RATIO`)."""
    ratio = (tuning.MOON_GIANT_PARENT_MAX_MASS_RATIO if parent_is_giant
             else tuning.MOON_TERRESTRIAL_PARENT_MAX_MASS_RATIO)
    return ratio * parent_mass_kg


def draw_moon_mass_kg(parent_mass_kg, parent_is_giant, limit_kg=None):
    """
    A moon's mass (GEN.646): `dN/dlogM ~ M^-MOON_MASS_FUNCTION_SLOPE` between `MOON_MIN_MASS_KG` and the parent's ceiling
    (`moon_mass_ceiling_kg`, or `limit_kg` where lower). One draw, from the stream bound now.

    Returns:
        float: The mass in kg, or `None` when the ceiling is under `MOON_MIN_MASS_KG` (the parent holds no moon that light).
    """
    ceiling = moon_mass_ceiling_kg(parent_mass_kg, parent_is_giant)
    if limit_kg is not None:
        ceiling = min(ceiling, limit_kg)
    if ceiling <= tuning.MOON_MIN_MASS_KG:
        return None
    return sample_power_law(tuning.MOON_MIN_MASS_KG, ceiling, tuning.MOON_MASS_FUNCTION_SLOPE)
