# planetgen/physics/bulk_draw.py

"""
Drawing a planet's or moon's bulk composition and solving its mass from it (GEN.428, docs/design/bulk-composition.md,
the contracts of docs/design/composition-spine-contracts.md)
==========================================================================================================================

Route 2 of the composition model: the composition comes first and the mass, density and gravity follow from it
(`physics/bulk_structure.py`, GEN.429). The class still decides the radius (as it always did); the composition decides how
much mass that radius holds.

1. **The star sets the dry ratios.** The host star's own abundances (`stellar_composition.key_draws` and
   `stellar_abundances`, the same star key its [Fe/H] stream uses) give the number ratios Fe/Mg, Si/Mg, Al/Mg, Ca/Mg, Ni/Mg
   and C/O. A rocky body takes all of the star's Mg, Si, Ca and Al into its rock as oxides and the Fe that its mantle does
   not hold as FeO into its metal core (`dry_partition`). The Earth-like star gives a core mass fraction near 0.3 with no
   draw, and an iron-rich or iron-poor star moves it, which is the observed trend (Dorn et al. 2015; Santos et al. 2015).
2. **What is drawn** is only what the star cannot say, in this order: the kind of rock (`ROCK_KIND_TABLE`, a gated possibility
   table: silicate, or silicon carbide and graphite for a star with C/O of 0.8 or more, Bond et al. 2010; Unterborn et al.
   2014), the water share (log-uniform over its class's range in `hydrosphere.WATER_MASS_FRACTION`, as the hydrosphere drew it
   before), for an icy class the other ices (`comet_population`'s volatiles against the water), the mantle's FeO, the core's
   light elements, a small body's porosity, a carbon body's graphite. Every range is in `RANGES` with where it comes from.
3. **Mass from radius.** The class radius is fixed and the mass that the composition needs for it is found by bisection on
   `bulk_structure.structure`; a body too dense to reach its radius below 20 Earth masses is set to the radius of 20.
   A small body (below `SMALL_BODY_MASS_KG`) is a heap with a porosity and its mass is its volume times the mixing density.
4. **Class S and a body given both a mass and a radius** go the other way (`origin` is `inverse`): the mass and the radius
   (for class S, the radius and a density from the class's range) say the core mass fraction, as the observer's route does.

The draw order is the order above; it is part of the generator epoch.
"""

import functools
import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.galaxy import scales
from planetgen.generation import possibility_table as pt
from planetgen.physics import (bulk_composition, bulk_structure, comet_population, constants, elements, hydrosphere,
                               stellar_abundances, stellar_composition)
from planetgen.util import draw
from planetgen.util.random import log_uniform

SOURCES = {
    "bond2010": "Bond, O'Brien and Lauretta 2010, ApJ 715, 1050 (planets of stars with C/O above 0.8 are carbon-rich)",
    "unterborn2014": "Unterborn, Dismukes and Panero 2014, ApJ 793, 124 (carbon-rich interiors: silicon carbide, graphite)",
    "class": "planetGen class definition (tuning.PLANET_CLASSES): the class is named for a carbon-based surface",
}

RANGES = {
    "rock_feo": {"default": (0.02, 0.18), "B": (0.0, 0.01)},
    "core_light_elements": (0.03, 0.12),
    "graphite_over_rock": (0.05, 0.40),
    "porosity": (0.0, 0.30),
    "other_ice_comet_class": "jupiter_family",
}
"""dict: The drawn ranges.

- `rock_feo`: FeO as a mass fraction of the rock, uniform. Earth's bulk silicate has 0.081, the Moon about 0.10, Mars 0.18
  (Khan and Connolly 2008), Mercury under 0.01 (Nittler et al. 2011): class B takes the Mercury-like range.
- `core_light_elements`: light elements (S, Si, O, C, H) as a mass fraction of the metal, uniform; Earth's core has 0.05 to
  0.10 (Hirose et al. 2021).
- `graphite_over_rock`: a carbon body's graphite against its silicon carbide and oxide rock, uniform; a stand-in range until the
  carbon mantle's equation of state (GEN.428's open item) says what the interiors hold.
- `porosity`: the void share of a small icy body, uniform; a stand-in until viscous compaction (lane K) replaces it.
- `other_ice_comet_class`: whose volatile proportions an icy class takes for its ices other than water.
"""

DENSITY_BOUNDS_G_CM3 = (0.6, 25.0)
"""tuple: The lowest and highest mean density a terrestrial body can have, g/cm3: ice with the 30 percent porosity of `RANGES` (0.917 times 0.7) and
pure iron at 20 Earth masses (Seager et al. 2007: about 23 g/cm3). The hard bounds of the plausibility checks."""

MAX_MASS_EARTH = 20.0
"""float: The heaviest terrestrial body, Earth masses (the top of the closed forms)."""
MIN_MASS_KG = 1.0e12
"""float: The lightest body the mass solve searches."""
RADIUS_TOLERANCE = 1.0e-6
"""float: The mass solve ends when its bracket is this close, as a fraction of the mass."""
MAX_CMF = 0.9
"""float: The largest core mass fraction (of the dry part) the inverse route returns."""
CARBON_OVER_OXYGEN = 0.8
"""float: A star at or above this C/O makes carbon-rich rock (Bond et al. 2010)."""


class DrawError(ValueError):
    """A composition the draw cannot make."""


@dataclass(frozen=True)
class HostRatios:
    """
    The host star's number ratios against Mg and its C/O (`host_ratios`).

    Attributes:
        fe, si, al, ca, ni (float): Atoms per Mg atom.
        c_over_o (float): Carbon atoms per oxygen atom.
        y_envelope (float): The helium share of the star's hydrogen and helium by mass.
        z (float): The star's metal mass fraction.
    """

    fe: float
    si: float
    al: float
    ca: float
    ni: float
    c_over_o: float
    y_envelope: float
    z: float


@functools.lru_cache(maxsize=4096)
def _ratios(fe_h, age_gy, key):
    draws = stellar_composition.key_draws(key)
    branch = stellar_abundances.branch_from_uniform(fe_h, age_gy, draws.u_branch)
    offsets = stellar_abundances.abundance_offsets(fe_h, branch, draws.normals, draws.u_enhanced, draws.u_carbon)
    star = stellar_abundances.composition(fe_h, offsets)
    ratio = stellar_abundances.number_ratio
    return HostRatios(fe=ratio(star, "Fe", "Mg"), si=ratio(star, "Si", "Mg"), al=ratio(star, "Al", "Mg"),
                      ca=ratio(star, "Ca", "Mg"), ni=ratio(star, "Ni", "Mg"), c_over_o=ratio(star, "C", "O"),
                      y_envelope=star.y / (star.x + star.y), z=star.z)


def star_key(star):
    """The key of `star`'s abundance draws: the one its [Fe/H] stream is seeded with (mass, age and galactic distance)."""
    distance_ly = getattr(star, "galactic_center_dist_ly", None)
    if distance_ly is None:
        distance_ly = scales.default_galactocentric_ly()
    return stellar_composition.abundance_key(mass=star.mass, age_gy=star.age, distance_ly=distance_ly)


REMNANT_DISC_RATIOS = HostRatios(fe=0.8, si=0.9, al=0.07, ca=0.06, ni=0.045, c_over_o=2.0, y_envelope=0.25, z=1.0)
"""HostRatios: What the planets of a compact remnant are made of (the carbon-rich debris disc of a disrupted white dwarf, PSR B1257+12's;
exotic-environments-planets-and-compact-binaries.md section 2.1): refractories in the proportions of the Sun's and C/O of 2, which
makes carbon rock. A stand-in until the disc's own abundances are modelled; no star of that kind has an [Fe/H]."""


def host_ratios(star):
    """
    The number ratios of the elements a rocky body is made of, from `star`'s own abundances; a compact remnant (it has no [Fe/H]) gives
    `REMNANT_DISC_RATIOS`.
    """
    fe_h = getattr(star, "metallicity_fe_h", None)
    if fe_h is None:
        return REMNANT_DISC_RATIOS
    return _ratios(fe_h, star.age, star_key(star))


ROCK_KIND_TABLE = pt.REGISTRY.register(pt.Table(
    "bulk-rock-kind", 1,
    [
        pt.Row("silicate", tuple(tuning.PLANET_CLASSES.keys() - {"C"}), 1.0, SOURCES["bond2010"],
               conditions=(pt.Condition("c_over_o", maximum=CARBON_OVER_OXYGEN - 1e-12),),
               notes="a star under C/O 0.8 makes silicate rock"),
        pt.Row("carbide-graphite", (pt.ANY_CLASS,), 1.0, SOURCES["unterborn2014"],
               conditions=(pt.Condition("c_over_o", minimum=CARBON_OVER_OXYGEN),),
               notes="a star at C/O 0.8 or more makes silicon carbide and graphite"),
        pt.Row("carbide-graphite-by-class", ("C",), 1.0, SOURCES["class"], state=pt.PARTLY_VERIFIED,
               conditions=(pt.Condition("c_over_o", maximum=CARBON_OVER_OXYGEN - 1e-12),),
               limit="a class definition, not a measured population: a carbon surface needs no carbon-rich star here",
               notes="class C is carbon rock whatever its star"),
    ],
    description="Whether a rocky body's rock is silicate or silicon carbide and graphite, by the host star's C/O and the class."))
"""Table: The kind of rock, gated on the host's C/O (`c_over_o`) and the class; exactly one row passes."""


def _grams(formula):
    return elements.molar_mass(formula) * 1000.0


@dataclass(frozen=True)
class DryPartition:
    """The metal and rock of a rocky body per unit of its dry mass, and what the rock holds (`dry_partition`)."""

    metal: float
    rock: float
    rock_feo: float
    rock_al2o3: float
    rock_mg_si: float


def dry_partition(ratios, rock_feo, light_elements, carbon_rock, graphite_over_rock=0.0):
    """
    The metal and rock shares of a dry body from the star's ratios, per Mg atom of the body.

    The rock is MgO, SiO2, CaO and Al2O3 (SiC and graphite with `carbon_rock`) plus `rock_feo` of itself as FeO; the metal is
    the Fe the mantle does not hold, with the star's Ni per Fe, and `light_elements` of the metal's mass as S, Si, O and C.
    A mantle that would take more Fe than the star has takes all of it and the metal is zero.

    Args:
        ratios (HostRatios): The star's ratios.
        rock_feo (float): FeO as a mass fraction of the rock.
        light_elements (float): The light elements' share of the metal's mass.
        carbon_rock (bool): Silicon carbide and graphite in place of silicate.
        graphite_over_rock (float): For `carbon_rock`, the graphite mass against the rest of the rock.

    Returns:
        DryPartition: `metal + rock` is 1.
    """
    mgo = _grams("MgO")
    si_mass = ratios.si * (_grams("SiC") if carbon_rock else _grams("SiO2"))
    shell = mgo + si_mass + ratios.ca * _grams("CaO") + 0.5 * ratios.al * _grams("Al2O3")
    al2o3 = 0.5 * ratios.al * _grams("Al2O3")
    graphite = graphite_over_rock * shell if carbon_rock else 0.0
    rock_no_feo = shell + graphite
    rock = rock_no_feo / (1.0 - rock_feo)
    n_feo = rock_feo * rock / _grams("FeO")
    if n_feo > ratios.fe:
        n_feo = ratios.fe
        rock = rock_no_feo + n_feo * _grams("FeO")
        rock_feo = n_feo * _grams("FeO") / rock
    iron = (ratios.fe - n_feo) * _grams("Fe") + ratios.ni * _grams("Ni")
    metal = iron / (1.0 - light_elements)
    total = metal + rock
    return DryPartition(metal=metal / total, rock=rock / total, rock_feo=rock_feo, rock_al2o3=al2o3 / rock,
                        rock_mg_si=1.0 / ratios.si)


def _subject(body, ratios):
    return pt.Subject(body.planet_class, {"c_over_o": ratios.c_over_o})


def _feo_range(planet_class):
    ranges = RANGES["rock_feo"]
    return ranges.get(planet_class, ranges["default"])


def _water_and_ice(body):
    """Water as a share of the whole body from the class's range, and the other ices beside it for an icy class."""
    low, high = hydrosphere.WATER_MASS_FRACTION.get(body.planet_class, hydrosphere.DEFAULT_WATER_MASS_FRACTION)
    water = log_uniform(low, high)
    ice = 0.0
    if body.planet_class == "D":
        volatiles = comet_population.draw_volatiles(RANGES["other_ice_comet_class"])
        ice = water * sum(volatiles.percent.values()) / 100.0
        scale = min(1.0, 0.95 / (water + ice))
        water, ice = water * scale, ice * scale
    return water, ice


def _porosity(body):
    return draw.uniform(*RANGES["porosity"]) if body.planet_class == "D" else None


def mass_from_radius(bulk, radius_km, x_fe):
    """The mass in kg at which `bulk` has `radius_km`, by bisection on the structure; the radius it reaches if 20 Earth masses
    are not enough."""
    def radius_at(mass_kg):
        return bulk_structure.structure(bulk, mass_kg, x_fe=x_fe).radius_km

    top_kg = MAX_MASS_EARTH * constants.EARTH_MASS_TO_KG
    if radius_at(top_kg) < radius_km:
        return top_kg, radius_at(top_kg)
    low, high = math.log(MIN_MASS_KG), math.log(top_kg)
    if radius_at(MIN_MASS_KG) > radius_km:
        raise DrawError(f"a body of radius {radius_km!r} km is lighter than {MIN_MASS_KG!r} kg for its composition")
    while high - low > RADIUS_TOLERANCE:
        middle = 0.5 * (low + high)
        if radius_at(math.exp(middle)) > radius_km:
            high = middle
        else:
            low = middle
    return math.exp(0.5 * (low + high)), radius_km


def _normalise(metal, rock, water, ice):
    """The dry shares of `1 - water - ice`."""
    wet = water + ice
    return metal * (1.0 - wet), rock * (1.0 - wet)


def draw_terrestrial(body, radius_given, mass_given):
    """
    Draws `body`'s bulk composition and sets its mass, radius and density from it (see the module text).

    Args:
        body (Planet): A planet or moon of a terrestrial class with its class, radius and zone set; `body.star` is its star.
        radius_given, mass_given (bool): Whether the caller gave that value (it is kept).

    Sets:
        The attributes of `bulk_composition.BODY_FIELDS`, `mass`, `radius` and `density` (g/cm3).
    """
    ratios = host_ratios(body.star)
    kind = ROCK_KIND_TABLE.draw(body.planet_class, {"c_over_o": ratios.c_over_o}).name
    carbon_rock = kind != "silicate"
    water, ice = _water_and_ice(body)
    feo = draw.uniform(*_feo_range(body.planet_class))
    light = draw.uniform(*RANGES["core_light_elements"])
    porosity = _porosity(body)
    graphite = draw.uniform(*RANGES["graphite_over_rock"]) if carbon_rock else 0.0

    inverse = body.planet_class == "S" or (radius_given and mass_given)
    if inverse:
        _draw_inverse(body, ratios, water, ice, feo, light, carbon_rock, graphite, mass_given)
        return

    part = dry_partition(ratios, feo, light, carbon_rock, graphite)
    metal, rock = _normalise(part.metal, part.rock, water, ice)
    bulk = bulk_composition.BulkFractions(
        metal=metal, rock=rock, water=water, ice=ice, envelope=0.0, rock_feo=part.rock_feo, rock_al2o3=part.rock_al2o3,
        rock_mg_si=part.rock_mg_si, porosity=porosity, carbon_rock=carbon_rock, origin="drawn")
    _settle(body, bulk, mass_given)


def _x_fe(bulk):
    return bulk_structure.x_fe_from_feo_mass_fraction(bulk.rock_feo) if bulk.rock_feo else 0.0


def _settle(body, bulk, mass_given):
    """Sets the mass (or the radius, where the mass is given) of `body` from `bulk`, and stores the composition."""
    x_fe = _x_fe(bulk)
    if mass_given:
        body.radius = bulk_structure.structure(bulk, body.mass, x_fe=x_fe).radius_km
    else:
        body.mass, body.radius = mass_from_radius(bulk, body.radius, x_fe)
    _store(body, bulk)


def _store(body, bulk):
    for attribute, value in bulk.to_attributes().items():
        setattr(body, attribute, value)
    volume_m3 = 4.0 / 3.0 * math.pi * (body.radius * 1000.0) ** 3
    body.density = body.mass / volume_m3 / 1000.0


def _draw_inverse(body, ratios, water, ice, feo, light, carbon_rock, graphite, mass_given):
    """Class S and a body given both a mass and a radius: the core mass fraction from the mass and the radius."""
    if body.planet_class == "S" and not mass_given:
        low, high = tuning.PLANET_CLASSES["S"]["density_range"]
        density = draw.uniform(low, high)
        body.mass = 4.0 / 3.0 * math.pi * (body.radius * 1000.0) ** 3 * density * 1000.0
    part = dry_partition(ratios, feo, light, carbon_rock, graphite)
    x_fe = bulk_structure.x_fe_from_feo_mass_fraction(part.rock_feo) if part.rock_feo else 0.0
    mass_earth = body.mass / constants.EARTH_MASS_TO_KG
    radius_earth = body.radius / constants.EARTH_RADIUS_KM
    wet = water + ice
    cmf = bulk_structure.cmf_from_mass_radius(mass_earth, radius_earth, x_fe, 1.0, high=MAX_CMF)
    if cmf is None:
        cmf = MAX_CMF if radius_earth < bulk_structure.radius_fit_earth(mass_earth, MAX_CMF, x_fe, 1.0) else 0.0
    bulk = bulk_composition.BulkFractions(
        metal=cmf * (1.0 - wet), rock=(1.0 - cmf) * (1.0 - wet), water=water, ice=ice, envelope=0.0, rock_feo=part.rock_feo,
        rock_al2o3=part.rock_al2o3, rock_mg_si=part.rock_mg_si, porosity=None, carbon_rock=carbon_rock, origin="inverse")
    _store(body, bulk)
