# planetgen/physics/planet_mass.py

"""
Planet mass (GEN.650)
=====================

How heavy a planet is when no class letter supplies the mass: the disc's solids at the distance it formed, the
oligarchs a feeding zone merges into one planet, and the core-accretion rule that turns a heavy enough core into a gas
giant. A pure law over the star and the formation distance; it reads no class and writes none.

1. **Solids.** The isolation mass at the formation distance (`formation.isolation_mass_kg`, Lissauer 1993; Kokubo & Ida
   2000) from the minimum-mass nebula scaled for the star's mass and [Fe/H] and boosted beyond its snow line.
2. **Assembly.** Giant impacts merge the oligarchs: a planet is `1 / GIANT_IMPACT_SURVIVAL_FRACTION` isolation masses on
   average (mass conservation against the survivors' share), spread log-normally by `PLANET_ASSEMBLY_SCATTER_DEX`.
3. **Runaway.** A core of `GIANT_CORE_CRITICAL_MASS_EARTH` or more, beyond the snow line, that finds gas left (a share
   scaled with the star's [Fe/H]) becomes a gas giant of the observed mass function `dN/dlogM ~ M^-0.31` (Cumming et al.
   2008) between `GIANT_RUNAWAY_MIN_MASS_EARTH` and the lighter of 13 Jupiter masses and `PLANET_MAX_STAR_MASS_RATIO` of
   the star. A core that finds no gas stays a solid or ice-giant core of its assembled mass.
"""

import collections
import math

from planetgen import tuning
from planetgen.physics import constants, formation, metallicity
from planetgen.util import draw

PlanetMass = collections.namedtuple("PlanetMass", ("mass_kg", "core_kg", "runaway"))
"""`mass_kg` the planet's mass; `core_kg` the assembled solid mass it grew from (the whole mass unless `runaway`);
`runaway` whether it became a gas giant."""

_EARTH_KG = constants.EARTH_MASS_TO_KG
_BROWN_DWARF_LIMIT_KG = 13.0 * constants.JUPITER_MASS_TO_KG


def assembled_core_kg(distance_au, star_mass_kg, snow_line_au, fe_h=0.0, stream=None):
    """
    The mass of the planet the disc's solids build at `distance_au` before any gas runaway: the isolation mass times the
    assembly factor times a log-normal scatter, never above `PLANET_MAX_STAR_MASS_RATIO` of the star.

    Args:
        distance_au (float): The formation distance from the star, in AU.
        star_mass_kg (float): The star's mass, in kilograms.
        snow_line_au (float): The star's own snow line, in AU.
        fe_h (float): The star's [Fe/H] in dex.
        stream (draw.Stream or None): The body's own stream; `None` draws from the stream bound now.

    Returns:
        float: The assembled solid mass, in kilograms.
    """
    source = stream or draw
    scale = formation.disk_surface_density_scale(star_mass_kg, fe_h) * tuning.PLANET_DISC_SOLIDS_ENHANCEMENT
    density = formation.mmsn_surface_density_gcm2(distance_au, snow_line_au, scale)
    isolation = formation.isolation_mass_kg(distance_au, density, star_mass_kg)
    scatter = 10.0 ** (tuning.PLANET_ASSEMBLY_SCATTER_DEX * source.gauss(0.0, 1.0))
    mass = isolation * scatter / tuning.GIANT_IMPACT_SURVIVAL_FRACTION
    return min(mass, tuning.PLANET_MAX_STAR_MASS_RATIO * star_mass_kg)


def runaway_mass_kg(star_mass_kg, stream=None):
    """A runaway giant's mass: `dN/dlogM ~ M^-slope` between `GIANT_RUNAWAY_MIN_MASS_EARTH` and the lighter of 13 Jupiter
    masses and `PLANET_MAX_STAR_MASS_RATIO` of the star (the floor wins when the star is too light for any giant)."""
    source = stream or draw
    low = tuning.GIANT_RUNAWAY_MIN_MASS_EARTH * _EARTH_KG
    high = max(low, min(_BROWN_DWARF_LIMIT_KG, tuning.PLANET_MAX_STAR_MASS_RATIO * star_mass_kg))
    slope = tuning.GIANT_MASS_FUNCTION_SLOPE
    # Inverse of the cumulative distribution of dN/dM ~ M^-(1 + slope) on [low, high].
    low_power, high_power = low ** -slope, high ** -slope
    return (low_power - source.random() * (low_power - high_power)) ** (-1.0 / slope)


def planet_mass(distance_au, star_mass_kg, snow_line_au, fe_h=0.0, stream=None):
    """
    A planet's mass at its formation distance from a star (module docstring). Draws, in this order, the assembly
    scatter, then (only for a core past the critical mass beyond the snow line) the gas-survival draw and the runaway
    mass, all from `stream`, so the same stream gives the same planet.

    Args:
        distance_au (float): The formation distance from the star, in AU.
        star_mass_kg (float): The star's mass, in kilograms.
        snow_line_au (float): The star's own snow line, in AU.
        fe_h (float): The star's [Fe/H] in dex.
        stream (draw.Stream or None): The body's own stream; `None` draws from the stream bound now.

    Returns:
        PlanetMass: The mass, the core it grew from and whether it ran away.
    """
    source = stream or draw
    core = assembled_core_kg(distance_au, star_mass_kg, snow_line_au, fe_h, source)
    critical = tuning.GIANT_CORE_CRITICAL_MASS_EARTH * _EARTH_KG
    if distance_au >= snow_line_au and core >= critical:
        share = tuning.GIANT_GAS_SURVIVAL_SHARE * metallicity.giant_occurrence_factor(fe_h)
        if source.random() < share:
            return PlanetMass(max(runaway_mass_kg(star_mass_kg, source), core), core, True)
    return PlanetMass(core, core, False)
