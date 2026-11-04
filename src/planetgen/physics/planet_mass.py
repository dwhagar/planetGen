# planetgen/physics/planet_mass.py

"""
Planet mass (GEN.650)
=====================

How heavy a planet is when no class letter supplies the mass: the disc's solids at the distance it formed, the
oligarchs a feeding zone merges into one planet, and the core-accretion rule that turns a heavy enough core into a gas
giant. A pure law over the star and the formation distance; it reads no class and writes none.

1. **Solids.** The isolation mass at the formation distance (`formation.isolation_mass_kg`, [R2709]; Kokubo & Ida
   2000) for a feeding zone `OLIGARCH_SEPARATION_MEDIAN` mutual Hill radii wide, from the minimum-mass nebula scaled for the star's
   mass (`PLANET_SOLIDS_STELLAR_MASS_EXPONENT`), [Fe/H] and the disc's own scatter (`DISC_SOLIDS_SCATTER_DEX`, [R0114]),
   boosted beyond the snow line and steeper inside 1 AU (`PLANET_INNER_DISC_STEEPENING`). Beyond `PLANET_OUTER_GROWTH_AU` the slow
   outer disc has built less: the mass falls as a power of the distance. `planet_architecture` builds a whole chain of these
   (`oligarch_core_kg`) and thins it; a planet asked for alone at a distance draws one.
2. **Runaway.** A core of `GIANT_CORE_CRITICAL_MASS_EARTH` or more, beyond the snow line, that finds gas left (a share
   scaled with the star's [Fe/H]) becomes a gas giant of the observed mass function `dN/dlogM ~ M^-0.31` (Cumming et al.
   2008) between `GIANT_RUNAWAY_MIN_MASS_EARTH` and the lighter of 13 Jupiter masses and `PLANET_MAX_STAR_MASS_RATIO` of
   the star. The share of cores that find gas falls beyond `GIANT_GAS_SURVIVAL_CUT_AU` as the square of the ratio. A core
   that finds no gas stays a solid or ice-giant core of its mass, at most `GIANT_NO_GAS_CORE_MAX_EARTH`.
"""

import collections
import math

from planetgen import tuning
from planetgen.physics import constants, formation, mass_boundaries, metallicity, unit_registry as units
from planetgen.util import draw

_kg_to_msun = units.converter("kg", "Msun")

PlanetMass = collections.namedtuple("PlanetMass", ("mass_kg", "core_kg", "runaway"))
"""`mass_kg` the planet's mass; `core_kg` the assembled solid mass it grew from (the whole mass unless `runaway`);
`runaway` whether it became a gas giant."""

_EARTH_KG = constants.EARTH_MASS_TO_KG
_BROWN_DWARF_LIMIT_KG = mass_boundaries.kg("deuterium_burning_limit")


def solids_surface_density_gcm2(distance_au, star_mass_kg, snow_line_au, metal_ratio=1.0, scatter=1.0):
    """
    The disc's solid surface density at `distance_au`, in g/cm^2: the minimum-mass nebula scaled for the star's mass
    (`PLANET_SOLIDS_STELLAR_MASS_EXPONENT`), its [Fe/H] (`metal_ratio`) and `PLANET_DISC_SOLIDS_ENHANCEMENT`, boosted beyond
    the snow line and steeper inside 1 AU (`PLANET_INNER_DISC_STEEPENING`); `scatter` multiplies it (a disc heavier or lighter than
    the median of its star's mass, 1 for the median).
    """
    solar_masses = _kg_to_msun(star_mass_kg)
    scale = solar_masses ** tuning.PLANET_SOLIDS_STELLAR_MASS_EXPONENT * metal_ratio * tuning.PLANET_DISC_SOLIDS_ENHANCEMENT * scatter
    density = formation.mmsn_surface_density_gcm2(distance_au, snow_line_au, scale)
    if distance_au < 1.0:
        density *= distance_au ** -tuning.PLANET_INNER_DISC_STEEPENING
    return density


def oligarch_core_kg(distance_au, star_mass_kg, snow_line_au, metal_ratio, scatter, separation):
    """
    The mass of the oligarch a feeding zone `separation` mutual Hill radii wide builds at `distance_au` in a disc scattered by `scatter`: the
    isolation mass, less beyond `PLANET_OUTER_GROWTH_AU` where the disc has built less, never above `PLANET_MAX_STAR_MASS_RATIO` of the star
    nor above the brown dwarf edge (GEN.1096: a metal-rich star of 1.5 to 2 solar masses has a 0.01 ratio of 15 to 21 Jupiter masses, and
    `gas_outcome` returns the larger of the runaway mass and the core, so the core was the planet's mass).
    """
    density = solids_surface_density_gcm2(distance_au, star_mass_kg, snow_line_au, metal_ratio, scatter)
    mass = formation.isolation_mass_kg(distance_au, density, star_mass_kg, separation)
    if distance_au > tuning.PLANET_OUTER_GROWTH_AU:
        mass *= (tuning.PLANET_OUTER_GROWTH_AU / distance_au) ** tuning.PLANET_OUTER_GROWTH_EXPONENT
    return min(mass, tuning.PLANET_MAX_STAR_MASS_RATIO * star_mass_kg, _BROWN_DWARF_LIMIT_KG)


def draw_disc_scatter(stream=None):
    """The factor a disc's solids are scattered by about the median of its star's mass: 10 to a normal of `DISC_SOLIDS_SCATTER_DEX` dex."""
    source = stream or draw
    return 10.0 ** (tuning.DISC_SOLIDS_SCATTER_DEX * source.gauss(0.0, 1.0))


def draw_separation(median, stream=None):
    """A spacing in mutual Hill radii: log-normal about `median` with `OLIGARCH_SEPARATION_SIGMA_LN`, floored at `OLIGARCH_SEPARATION_MIN`."""
    source = stream or draw
    return max(tuning.OLIGARCH_SEPARATION_MIN, median * math.exp(tuning.OLIGARCH_SEPARATION_SIGMA_LN * source.gauss(0.0, 1.0)))


def assembled_core_kg(distance_au, star_mass_kg, snow_line_au, metal_ratio=1.0, stream=None):
    """
    The mass of the planet the disc's solids build alone at `distance_au` before any gas runaway: one oligarch (`oligarch_core_kg`) of a disc
    scattered at random (`draw_disc_scatter`) and spaced at random (`draw_separation`).

    Args:
        distance_au (float): The formation distance from the star, in AU.
        star_mass_kg (float): The star's mass, in kilograms.
        snow_line_au (float): The star's own snow line, in AU.
        metal_ratio (float): The star's Z / Z0 (`bulk_draw.metal_ratio_to_solar`, GEN.464).
        stream (draw.Stream or None): The body's own stream; `None` draws from the stream bound now.

    Returns:
        float: The assembled solid mass, in kilograms.
    """
    source = stream or draw
    scatter = draw_disc_scatter(source)
    separation = draw_separation(tuning.OLIGARCH_SEPARATION_MEDIAN, source)
    return oligarch_core_kg(distance_au, star_mass_kg, snow_line_au, metal_ratio, scatter, separation)


def saturate_kg(mass_kg, cap_kg):
    """`mass_kg` pressed softly under `cap_kg`: close to itself well below the cap, never above it, with no pile-up at a hard edge
    (`mass / (1 + (mass / cap) ** k) ** (1 / k)`, k = `PLANET_CORE_SATURATION_SHARPNESS`)."""
    sharpness = tuning.PLANET_CORE_SATURATION_SHARPNESS
    return mass_kg / (1.0 + (mass_kg / cap_kg) ** sharpness) ** (1.0 / sharpness)


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


def planet_mass(distance_au, star_mass_kg, snow_line_au, fe_h=0.0, metal_ratio=1.0, stream=None):
    """
    A planet's mass at its formation distance from a star (module docstring). Draws, in this order, the disc scatter and
    the spacing, then (only for a core past the critical mass beyond the snow line) the gas-survival draw and the runaway
    mass, all from `stream`, so the same stream gives the same planet.

    Args:
        distance_au (float): The formation distance from the star, in AU.
        star_mass_kg (float): The star's mass, in kilograms.
        snow_line_au (float): The star's own snow line, in AU.
        fe_h (float): The star's [Fe/H] in dex: the giant-planet occurrence law (Fischer and Valenti) is stated in it.
        metal_ratio (float): The star's Z / Z0 (`bulk_draw.metal_ratio_to_solar`, GEN.464): the disc's solids law is stated in it.
        stream (draw.Stream or None): The body's own stream; `None` draws from the stream bound now.

    Returns:
        PlanetMass: The mass, the core it grew from and whether it ran away.
    """
    source = stream or draw
    core = assembled_core_kg(distance_au, star_mass_kg, snow_line_au, metal_ratio, source)
    return gas_outcome(core, distance_au, star_mass_kg, snow_line_au, fe_h, source)


def gas_outcome(core_kg, distance_au, star_mass_kg, snow_line_au, fe_h=0.0, stream=None):
    """
    What a solid core of `core_kg` formed at `distance_au` becomes: a gas giant when it reaches the critical mass beyond the snow
    line and finds gas left (one draw, then the runaway mass), else a solid or ice-giant core pressed under its cap. Draws, only
    for a core past the critical mass beyond the snow line, the gas-survival draw and then the runaway mass.

    Returns:
        PlanetMass: The mass, the core it grew from and whether it ran away.
    """
    source = stream or draw
    core = core_kg
    critical = tuning.GIANT_CORE_CRITICAL_MASS_EARTH * _EARTH_KG
    if distance_au >= snow_line_au and core >= critical:
        found = (tuning.GIANT_GAS_SURVIVAL_SHARE * 10.0 ** (tuning.GIANT_GAS_SURVIVAL_METALLICITY_EXPONENT * fe_h)
                 * metallicity.giant_occurrence_factor(fe_h))
        found *= min(1.0, (tuning.GIANT_GAS_SURVIVAL_CUT_AU / distance_au) ** 2)
        if source.random() < found:
            core = min(core, _BROWN_DWARF_LIMIT_KG)      # two giants' cores that merged past the deuterium limit leave a planet of at most the limit
            return PlanetMass(max(runaway_mass_kg(star_mass_kg, source), core), core, True)
    cap = tuning.GIANT_NO_GAS_CORE_MAX_EARTH if distance_au >= snow_line_au else tuning.PLANET_INNER_CORE_MAX_EARTH
    core = saturate_kg(core, cap * _EARTH_KG)
    return PlanetMass(core, core, False)
