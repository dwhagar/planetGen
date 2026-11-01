# planetgen/physics/body_strength.py

"""
What a body stores of its strength and its binding: the roundness ratio, the yield strength and the binding energy (GEN.1159, GEN.1178)
=========================================================================================================================================

A planet or a moon has its mass, radius, density and composition from the bulk draw (`bulk_draw`, `bulk_structure`). From them this
module works out, in one place and without drawing from the run's stream:

- `roundness_ratio`: the body's mean radius over the radius at which its own weight beats the strength of what it is made of
  (`roundness.roundness_ratio`); 1 or more is round, less is not (the classes and the Asteroid type of GEN.1160 read it).
- `yield_strength_pa`: the whole-body strength of the body's material mix, drawn from the range the research gives for each material
  (`roundness.drawn_strengths_pa`, GEN.1178), so the line between round and not round is fuzzy as it is in the Solar System. It is the
  cold strength: the engine's temperature is the one a body took its shape at (it creeps round when warm), not its surface temperature,
  and no body stores that yet, so none is passed.
- `binding_energy_j`: the energy that assembling the body from infinity releases (`binding_energy.binding_energy_j`): two uniform
  shells for a rocky body with a core, the n = 1 polytrope for a gas giant.

The strength quantile is drawn from a stream of the body's own (its mass, radius, distance and whether it is a moon), so the run's
stream does not move and the same body gets the same strength each time it is made again.

A gas giant has no strength (it is a fluid held up by its own pressure): its ratio and strength are `None` and its binding energy is
kept. A body without a mass or a radius has none of the three.

Use:

    body_strength.generate(planet)       # sets FIELDS on a planet or a moon
"""

import math

from planetgen.physics import binding_energy, planet_classes, roundness
from planetgen.physics import unit_registry as units
from planetgen.util import draw

FIELDS = ("roundness_ratio", "yield_strength_pa", "binding_energy_j")
"""tuple: What a planet or a moon stores of this module, as `planets` and `moons` columns (GEN.1159)."""

_km_to_m = units.converter("km", "m")
_g_cm3_to_kg_m3 = units.converter("g_cm3", "kg_m3")
_msun_to_kg = units.converter("Msun", "kg")
_kg_to_msun = units.converter("kg", "Msun")

STRENGTH_DRAW_STREAM = "body-strength"
"""str: The seed prefix of the stream a body's strength quantile is drawn from (GEN.1178)."""

CORE_FRACTION_LIMIT = 0.999
"""float: A body whose metal is this share or more of its dry mass is one uniform ball for the binding energy (no mantle to put around the core)."""


def ice_fraction(body):
    """The share of the body's mass that is ice of any kind: its water (every phase) and its other ices, 0 for a body with neither."""
    return (getattr(body, "water_mass_fraction", None) or 0.0) + (getattr(body, "ice_mass_fraction", None) or 0.0)


def _core(body, mass_kg, radius_m):
    """`(core mass, core radius)` in kg and m of a rocky body: the metal share of the mass in a core whose radius is the body's times the
    square root of the metal share of its dry mass (the convention of `interior_heat`), or `(0, 0)` where the body has no metal share."""
    metal = getattr(body, "metal_mass_fraction", None) or 0.0
    rock = getattr(body, "rock_mass_fraction", None) or 0.0
    if not (metal > 0.0 and rock > 0.0):
        return 0.0, 0.0
    share = metal / (metal + rock)
    if share >= CORE_FRACTION_LIMIT:
        return 0.0, 0.0
    return metal * mass_kg, radius_m * math.sqrt(share)


def _is_giant(body):
    """Whether the body is a giant: its body type says so, or its class is one of the giant classes (an ice giant of rock and water is a terrestrial body type with a giant class)."""
    letter = getattr(body, "planet_class", None)
    return getattr(body, "body_type", None) == planet_classes.GIANT or (letter in planet_classes.BY_LETTER and planet_classes.BY_LETTER[letter].kind == planet_classes.GIANT)


def solid(mass_kg, radius_m, density_kg_m3, ice, metal, porosity=0.0, temperature_k=None, quantile=0.5, core_kg=0.0, core_m=0.0):
    """
    The `FIELDS` of a solid body from its mass, mean radius, density and composition: `{name: value}`. `ice` and `metal` are the shares of its mass,
    `temperature_k` the temperature the body took its shape at (`None`: the cold whole-body strength, which every caller passes until a body stores its shape-setting temperature), `quantile` the strength draw from 0 to 1, and `core_kg` and
    `core_m` the core of a differentiated body (0 for none).
    """
    ice = min(max(ice or 0.0, 0.0), 1.0)
    metal = min(max(metal or 0.0, 0.0), 1.0 - ice)
    return {
        "roundness_ratio": roundness.roundness_ratio(radius_m, density_kg_m3, ice, metal, porosity or 0.0, temperature_k, draw=quantile),
        "yield_strength_pa": roundness.yield_strength_pa(ice, metal, porosity or 0.0, temperature_k, draw=quantile),
        "binding_energy_j": binding_energy.binding_energy_j("rocky", mass_kg, radius_m, core_kg, core_m),
    }


def values(body, quantile=None):
    """
    The `FIELDS` of a planet or a moon: `{name: value}`, every one `None` for a body that lacks a mass or a radius. `quantile` is the strength
    draw from 0 to 1 (`None` draws it from the body's own stream).
    """
    out = dict.fromkeys(FIELDS)
    mass_kg, radius_km = getattr(body, "mass", None), getattr(body, "radius", None)
    if not (mass_kg and radius_km and mass_kg > 0.0 and radius_km > 0.0):
        return out
    radius_m = _km_to_m(radius_km)
    if _is_giant(body):
        out["binding_energy_j"] = binding_energy.binding_energy_j("gas-giant", mass_kg, radius_m)
        return out
    core_kg, core_m = _core(body, mass_kg, radius_m)
    density = getattr(body, "density", None)
    density_kg_m3 = _g_cm3_to_kg_m3(density) if density else roundness.bulk_density_kg_m3(mass_kg, radius_m)
    if quantile is None:
        quantile = quantile_of(body)
    return solid(mass_kg, radius_m, density_kg_m3, ice_fraction(body), getattr(body, "metal_mass_fraction", None), getattr(body, "porosity", None),
                 None, quantile, core_kg, core_m)


def quantile_for(*parts):
    """A strength quantile, 0 to 1, from a stream seeded with `parts` (the numbers that make one body unlike another): the run's stream does not move."""
    with draw.bound(f"{STRENGTH_DRAW_STREAM}:" + ":".join(repr(part) for part in parts)):
        return draw.random()


def quantile_of(body):
    """The strength quantile of a planet or a moon, from a stream of its own (mass, radius, distance and whether it is a moon)."""
    return quantile_for(float(body.mass), float(body.radius), float(getattr(body, "distance", 0.0) or 0.0), int(bool(getattr(body, "is_moon", False))))


def generate(body):
    """Sets `FIELDS` on a planet or a moon (call once its mass, radius, composition and temperature are final)."""
    for name, value in values(body).items():
        setattr(body, name, value)


COMET_FIELDS = ("density_g_cm3", "ice_mass_fraction", *FIELDS)
"""tuple: What a comet stores of this module (GEN.1159): its bulk density, g/cm3, the share of its mass that is ice of any kind, and `FIELDS`."""


def generate_comet(comet):
    """Sets `COMET_FIELDS` on a comet or an interstellar comet (GEN.1162) from its makeup, nucleus diameter and mass (a sphere at the makeup's bulk density), at the cold strength."""
    makeup = comet.makeup
    radius_m = comet.nucleus_diameter_km * 500.0
    density = makeup.bulk_density_kg_m3 if makeup.bulk_density_kg_m3 is not None else makeup.density_at(comet.nucleus_diameter_km / 2.0)
    comet.density_g_cm3 = density / 1000.0
    comet.ice_mass_fraction = makeup.water + makeup.ice
    quantile = quantile_for(float(comet.mass_kg), float(comet.nucleus_diameter_km), float(getattr(comet, "perihelion_distance_au", 0.0)),
                            float(getattr(comet, "eccentricity", getattr(comet, "velocity_kms", 0.0))))
    for name, value in solid(comet.mass_kg, radius_m, density, comet.ice_mass_fraction, makeup.metal, makeup.porosity, None, quantile).items():
        setattr(comet, name, value)


POPULATION_FIELDS = ("density_g_cm3", "ice_mass_fraction")
"""tuple: What a belt or an asteroid field stores of its makeup (GEN.1159): the mean bulk density of its members, g/cm3, and the mean share of their mass that is ice
of any kind. The roundness ratio, strength and binding energy belong to a member, not to the population (the member rows, GEN.1070, GEN.1161)."""


def population_values(population):
    """`POPULATION_FIELDS` of a belt or an asteroid field from its drawn makeup, in that order: `(None, None)` for one read back without a makeup, and `None`
    for the density of a makeup that has none (an icy body's density depends on its size)."""
    makeup = getattr(population, "makeup", None)
    if makeup is None:
        return None, None
    density = None if makeup.bulk_density_kg_m3 is None else makeup.bulk_density_kg_m3 / 1000.0
    return density, makeup.water + makeup.ice


def grain_density_kg_m3(ice, metal=0.0):
    """The density of a body of this ice and metal share of its mass and no porosity: its grains mixed by volume (`roundness.GRAIN_DENSITIES_KG_M3`), the rest rock."""
    rock = 1.0 - ice - metal
    grains = roundness.GRAIN_DENSITIES_KG_M3
    return 1.0 / (ice / grains["ice"] + metal / grains["metal"] + rock / grains["rock"])


def round_mass_for(ice, metal=0.0, quantile=None):
    """The mass, kg, at which a body of this ice and metal share (no porosity) turns round: its round radius at the median strength (or at `quantile`), as a sphere of its grain density."""
    return roundness.round_mass_kg(grain_density_kg_m3(ice, metal), ice, metal, draw=quantile)


def ratio_for(mass_kg, ice, metal=0.0, quantile=None):
    """The roundness ratio of a body of this mass, ice and metal share (no porosity), at the median strength (or at `quantile`): 1 or more is round."""
    return roundness.roundness_ratio_from_mass(mass_kg, grain_density_kg_m3(ice, metal), ice, metal, draw=quantile)


STELLAR_FIELDS = ("density_g_cm3", "binding_energy_j")
"""tuple: What a star, a white dwarf, a neutron star or a black hole stores of its bulk (GEN.1162): the mean density, g/cm3, and the binding energy, J. A fluid or degenerate object has
no strength to measure against its weight, so it has no roundness ratio or yield strength."""

WHITE_DWARF_CLASSES = ("VII", "D")
"""tuple: The Yerkes classes of a white dwarf (`D` is the model's `VII` under its own name)."""

BROWN_DWARF_MAX_MASS_SOLAR = 0.08
"""float: A star under this mass (the hydrogen-burning limit, `mass_boundaries`) is a brown dwarf for the binding energy (the n = 3/2 polytrope)."""


def stellar_kind(yerkes_class, mass_kg):
    """The `binding_energy.KINDS` entry of a star row: `black-hole` for the marker `BH`, `neutron-star` for `NS`, `white-dwarf` for `VII` or `D`, `brown-dwarf` under the hydrogen-burning limit, else `star`."""
    if yerkes_class == "BH":
        return "black-hole"
    if yerkes_class == "NS":
        return "neutron-star"
    if yerkes_class in WHITE_DWARF_CLASSES:
        return "white-dwarf"
    if _kg_to_msun(mass_kg) < BROWN_DWARF_MAX_MASS_SOLAR:
        return "brown-dwarf"
    return "star"


def bulk_values(kind, mass_kg, radius_km):
    """`STELLAR_FIELDS` of a mass object of `kind` (a `binding_energy.KINDS` entry other than `rocky`) and this mass and radius: `(density g/cm3, binding energy J)`, `(None, None)` for one without
    a mass or a radius. A black hole's radius is its event horizon, so its density is the mean inside the horizon and its binding energy `binding_energy.INFINITE`."""
    if not (mass_kg and radius_km and mass_kg > 0.0 and radius_km > 0.0):
        return None, None
    radius_m = _km_to_m(radius_km)
    density = roundness.bulk_density_kg_m3(mass_kg, radius_m) / 1000.0
    return density, binding_energy.binding_energy_j(kind, mass_kg, radius_m)


def star_values(star):
    """`STELLAR_FIELDS` of a `Star` or a compact remnant standing in for one (its `mass`, `radius`, `yerkes_class`)."""
    mass_kg = getattr(star, "mass", None)
    return bulk_values(stellar_kind(getattr(star, "yerkes_class", None), mass_kg or 0.0), mass_kg, getattr(star, "radius", None))


def remnant_values(kind, mass_solar, radius_km):
    """`STELLAR_FIELDS` of a standalone black hole or neutron star from its mass, solar masses, and its horizon or radius, km (`kind` is `black-hole` or `neutron-star`)."""
    return bulk_values(kind, _msun_to_kg(mass_solar), radius_km)


ROGUE_FIELDS = ("density_g_cm3", *FIELDS)
"""tuple: What a rogue planet stores of this module (GEN.1162): the mean density and `FIELDS`."""


def rogue_values(planet, water_fraction=0.0):
    """
    `ROGUE_FIELDS` of a `RoguePlanet`: `{name: value}`. A giant or a brown dwarf keeps its density and binding energy (the n = 1 polytrope for a giant, n = 3/2 for a brown dwarf). A terrestrial one is a rocky body
    of its drawn density, its ice the `water_fraction` of its mass (its ice shell and ocean) and no metal share stored, so its core is not modelled and its binding energy is the uniform ball's.
    """
    out = dict.fromkeys(ROGUE_FIELDS)
    mass_kg, radius_km = planet.mass_kg, planet.radius_km
    if not (mass_kg and radius_km and mass_kg > 0.0 and radius_km > 0.0):
        return out
    radius_m = _km_to_m(radius_km)
    out["density_g_cm3"] = roundness.bulk_density_kg_m3(mass_kg, radius_m) / 1000.0
    if planet.mass_bin == "brown-dwarf":
        out["binding_energy_j"] = binding_energy.binding_energy_j("brown-dwarf", mass_kg, radius_m)
    elif planet.planet_type == "g":
        out["binding_energy_j"] = binding_energy.binding_energy_j("gas-giant", mass_kg, radius_m)
    else:
        quantile = quantile_for(float(mass_kg), float(radius_km), 0.0, 0)
        out.update(solid(mass_kg, radius_m, _g_cm3_to_kg_m3(out["density_g_cm3"]), water_fraction, 0.0, 0.0, None, quantile))
    return out
