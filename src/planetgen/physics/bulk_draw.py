# planetgen/physics/bulk_draw.py

"""
Drawing a planet's or moon's bulk composition and solving its mass from it (GEN.428, docs/design/bulk-composition.md,
the contracts of docs/design/composition-spine-contracts.md)
==========================================================================================================================

Route 2 of the composition model: the composition comes first and the radius, density and gravity follow from it
(`physics/bulk_structure.py`, GEN.429). The mass comes from the formation law (`planet_mass`, GEN.650) or from a moon's size law
(`planet_kind.draw_moon_mass_kg`); a caller that gives a radius instead gets the mass that radius holds. No class letter is read: the
body's origin (`origin_of`: formed beyond its star's snow line or not) and its star say what it is made of.

1. **The star sets the dry ratios.** The host star's own abundances (`stellar_composition.key_draws` and
   `stellar_abundances`, the same star key its [Fe/H] stream uses) give the number ratios Fe/Mg, Si/Mg, Al/Mg, Ca/Mg, Ni/Mg
   and C/O. A rocky body takes all of the star's Mg, Si, Ca and Al into its rock as oxides and the Fe that its mantle does
   not hold as FeO into its metal core (`dry_partition`). The Earth-like star gives a core mass fraction near 0.3 with no
   draw, and an iron-rich or iron-poor star moves it, which is the observed trend [R1146, R3916].
2. **What is drawn** is only what the star cannot say, in this order: the kind of rock (`ROCK_KIND_TABLE`, a gated possibility
   table: silicate, or silicon carbide and graphite for a star with C/O of 0.8 or more, [R0455]; Unterborn et al.
   2014), the water share (a rocky body's from the Bern delivery law by its mass and the snow lines it formed at, `delivered_water`, GEN.849; an
   icy body's log-uniform from `tuning.ICY_WATER_MASS_FRACTION_FLOOR` up to the water share of the solids at its formation temperature, GEN.950), for an icy body the other ices (`comet_population`'s volatiles against the water), the mantle's FeO, the core's
   light elements, a small icy body's porosity, a carbon body's graphite, and then the iron branch of a rocky planet (`IRON_BRANCH_TABLE`: its
   iron follows the star's with a scatter, or is enriched, GEN.512). Every range is in `RANGES` with where it comes from.
3. **Radius from mass.** The mass is fixed and `bulk_structure.structure` gives the radius the composition holds it at. A caller that
   gave a radius gets the mass that the composition needs for it, found by bisection; a body too dense to reach its radius below 20
   Earth masses is set to the radius of 20. A small body (below `SMALL_BODY_MASS_KG`) is a heap with a porosity and its mass is its
   volume times the mixing density.
4. **A body given both a mass and a radius** goes the other way (`origin` is `inverse`): the mass and the radius say the core mass
   fraction, as the observer's route does.

The draw order is the order above; it is part of the generator epoch.
"""

import dataclasses
import functools
import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.galaxy import scales
from planetgen.generation import possibility_table as pt
from planetgen.physics import (bulk_composition, bulk_structure, comet_population, compaction, composition_text, condensate_shares, constants, delivered_water, elements,
                               hydrosphere, ice_species, moon_composition, stellar_abundances, stellar_composition)
from planetgen.util import draw
from planetgen.util.random import log_uniform
from planetgen.physics import unit_registry as units

_kg_to_mearth = units.converter("kg", "Mearth")
_mearth_to_kg = units.converter("Mearth", "kg")

SOURCES = {
    "bond2010": "[R0455] (planets of stars with C/O above 0.8 are carbon-rich)",
    "unterborn2014": "Unterborn, ApJ 793, 124 (carbon-rich interiors: silicon carbide, graphite)",
}

RANGES = {
    "enriched_cmf": (0.5, 0.8),
    "rock_feo": (0.02, 0.18),
    "core_light_elements": (0.03, 0.12),
    "graphite_over_rock": (0.05, 0.40),
    "porosity": (0.0, 0.30),
    "other_ice_comet_class": ice_species.COMET_CLASS,
}
"""dict: The drawn ranges.

- `enriched_cmf`: the core mass fraction of the dry part of an iron-enriched planet (`IRON_BRANCH_TABLE`), uniform, floored at
  the star-matched value; its top is the cap of the observed iron-rich rocky planets (about 0.8, [R3545]).
- `rock_feo`: FeO as a mass fraction of the rock, uniform. Earth's bulk silicate has 0.081, the Moon about 0.10, Mars 0.18
  [R2334]; Mercury, under 0.01, is outside it: the Mercury-like range belonged to class B,
  a class the generator no longer draws, and no condition of the body says it is reduced (a stand-in, tier "from needed output").
- `core_light_elements`: light elements (S, Si, O, C, H) as a mass fraction of the metal, uniform; Earth's core has 0.05 to
  0.10 [R1950].
- `graphite_over_rock`: a carbon body's graphite against its silicon carbide and oxide rock, uniform; a stand-in range until the
  carbon mantle's equation of state (GEN.428's open item) says what the interiors hold.
- `porosity`: the void share of a small icy body, uniform; a stand-in until viscous compaction (lane K) replaces it.
- `other_ice_comet_class`: whose volatile proportions an icy body takes for its ices other than water.
"""

IRON_SCATTER_DEX = 0.05
"""float: The scatter of a star-matched planet's Fe/Mg about its star's, dex (a normal draw): GEN.512 (`docs/design/composition/iron-branch-tuning.md`
section 3 gives the tuning and the bracket)."""
IRON_ENRICHED_WEIGHT = 0.01
"""float: The weight of the iron-enriched branch against 1 minus it for the star-matched branch (GEN.512, tuned to 1 of 21 planets 2 sigma
iron-rich in a homogeneous sample, an upper limit)."""

DENSITY_BOUNDS_G_CM3 = (0.6, 25.0)
"""tuple: The lowest and highest mean density a terrestrial body can have, g/cm3: ice with the 30 percent porosity of `RANGES` (0.917 times 0.7) and
pure iron at 20 Earth masses ([R3992]: about 23 g/cm3). The hard bounds of the plausibility checks."""

MAX_MASS_EARTH = 20.0
"""float: The heaviest terrestrial body, Earth masses (the top of the closed forms)."""
MIN_MASS_KG = 1.0e12
"""float: The lightest body the mass solve searches."""
RADIUS_TOLERANCE = 1.0e-6
"""float: The mass solve ends when its bracket is this close, as a fraction of the mass."""
MAX_CMF = 0.9
"""float: The largest core mass fraction (of the dry part) the inverse route returns."""
CARBON_OVER_OXYGEN = 0.8
"""float: A star at or above this C/O makes carbon-rich rock [R0455]."""


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
def _star_composition(fe_h, age_gy, key):
    draws = stellar_composition.key_draws(key)
    branch = stellar_abundances.branch_of(fe_h, age_gy)
    offsets = stellar_abundances.abundance_offsets(fe_h, branch, draws.normals, draws.u_enhanced, draws.u_carbon)
    return stellar_abundances.composition(fe_h, offsets)


@functools.lru_cache(maxsize=4096)
def _ratios(fe_h, age_gy, key):
    star = _star_composition(fe_h, age_gy, key)
    ratio = stellar_abundances.number_ratio
    return HostRatios(fe=ratio(star, "Fe", "Mg"), si=ratio(star, "Si", "Mg"), al=ratio(star, "Al", "Mg"),
                      ca=ratio(star, "Ca", "Mg"), ni=ratio(star, "Ni", "Mg"), c_over_o=ratio(star, "C", "O"),
                      y_envelope=star.y / (star.x + star.y), z=star.z)


def star_key(star):
    """The key of `star`'s abundance draws: the star's mass and galactic distance, not its age (GEN.984)."""
    distance_ly = getattr(star, "galactic_center_dist_ly", None)
    if distance_ly is None:
        distance_ly = scales.default_galactocentric_ly()
    return stellar_composition.abundance_key(mass=star.mass, distance_ly=distance_ly)


def metal_ratio_to_solar(star):
    """
    Z / Z0 of `star` (GEN.464): its metal mass fraction over the reference Sun-like mix's, the factor the disc's solids scale with. 1 for a
    compact remnant, which has no [Fe/H] (its disc is `REMNANT_DISC_RATIOS`, not drawn from a metal fraction).
    """
    if getattr(star, "metallicity_fe_h", None) is None:
        return 1.0
    return stellar_abundances.metal_ratio_to_solar(_ratios(star.metallicity_fe_h, star.age, star_key(star)).z)


def ice_forming_ratios(star):
    """
    `(carbon, oxygen)` of `star`: its C and O per atom of the rock-forming metal Mg, over the Sun's, the amounts `comet_ice_budget` and `condensate_shares` scale the
    carbon oxides and the water with (GEN.595); `None` for a compact remnant, which has no [Fe/H].
    """
    fe_h = getattr(star, "metallicity_fe_h", None)
    if fe_h is None:
        return None
    composition = _star_composition(fe_h, star.age, star_key(star))
    ratio = stellar_abundances.number_ratio
    return tuple(ratio(composition, element, "Mg") / ratio({}, element, "Mg") for element in ("C", "O"))


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
        pt.Row("silicate", (pt.ANY_CLASS,), 1.0, SOURCES["bond2010"],
               conditions=(pt.Condition("c_over_o", maximum=CARBON_OVER_OXYGEN - 1e-12),),
               notes="a star under C/O 0.8 makes silicate rock"),
        pt.Row("carbide-graphite", (pt.ANY_CLASS,), 1.0, SOURCES["unterborn2014"],
               conditions=(pt.Condition("c_over_o", minimum=CARBON_OVER_OXYGEN),),
               notes="a star at C/O 0.8 or more makes silicon carbide and graphite"),
    ],
    description="Whether a rocky body's rock is silicate or silicon carbide and graphite, by the host star's C/O."))
"""Table: The kind of rock, gated on the host's C/O (`c_over_o`); exactly one row passes."""


ROCKY, ICY = "rocky", "icy"
"""str: The two origins a body's makeup is drawn for (`origin_of`)."""
IRON_BRANCH_ORIGINS = (ROCKY,)
"""tuple: The origins with an iron branch: a rocky planet, not an icy one (a moon inherits its star's, below)."""

IRON_BRANCH_TABLE = pt.REGISTRY.register(pt.Table(
    "bulk-iron-branch", 1,
    [
        pt.Row("star-matched", IRON_BRANCH_ORIGINS, 1.0 - IRON_ENRICHED_WEIGHT,
               "[R0541, R3546]: planet iron follows the host's, 1:1",
               state=pt.PARTLY_VERIFIED, tunable=True,
               limit="the slope of planet against star iron is 0.94 (+1.02/-1.07) and 1.3 to 5.6 by fit: consistent with 1:1, not measured; the scatter is an upper limit",
               notes="Fe/Mg of the planet = the star's times 10 to a normal of IRON_SCATTER_DEX; the slope stays 1"),
        pt.Row("iron-enriched", IRON_BRANCH_ORIGINS, IRON_ENRICHED_WEIGHT,
               "[R3546]: 1 of 21 planets 2 sigma more iron-rich than its star; [R3978]: 1 of 11",
               state=pt.PARTLY_VERIFIED, tunable=True,
               limit="an upper limit: 0 of 20 distinct [R0541] to 5 of 22 ([R0032], heterogeneous); mass-radius scatter alone makes about 15 percent of mock samples show one",
               notes="core mass fraction uniform 0.5 to 0.8 (RANGES), never below the star-matched value"),
    ],
    description="Whether a rocky planet's iron follows its star's or is enriched (an iron-rich core)."))
"""Table: The iron branch of a rocky planet (GEN.512). The coreless branch has no count and is not drawn (derive it from the water fraction)."""


def iron_branch_ratios(body, ratios, feo, light, carbon_rock, graphite):
    """
    The host ratios as the planet inherits them: the star's Fe/Mg with the scatter of a star-matched planet, or an Fe/Mg that makes the
    core mass fraction a draw from `RANGES["enriched_cmf"]`. A moon, and a body of an origin outside `IRON_BRANCH_ORIGINS`, inherit the star's.
    """
    origin = origin_of(body)
    if body.is_moon or origin not in IRON_BRANCH_ORIGINS:
        return ratios
    branch = IRON_BRANCH_TABLE.draw(origin, {}).name
    return branch_ratios(ratios, branch, IRON_SCATTER_DEX, feo, light, carbon_rock, graphite)


def branch_ratios(ratios, branch, scatter_dex, feo, light, carbon_rock, graphite):
    """`ratios` as a planet of `branch` (a row of `IRON_BRANCH_TABLE`) inherits them; the draws are the scatter's normal or the target's uniform."""
    if branch == "star-matched":
        scatter = draw.gauss(0.0, 1.0) * scatter_dex
        return dataclasses.replace(ratios, fe=ratios.fe * 10.0 ** scatter)
    target = max(draw.uniform(*RANGES["enriched_cmf"]), dry_partition(ratios, feo, light, carbon_rock, graphite).metal)
    low, high = ratios.fe, ratios.fe
    while dry_partition(dataclasses.replace(ratios, fe=high), feo, light, carbon_rock, graphite).metal < target:
        high *= 2.0
    for _ in range(60):
        middle = 0.5 * (low + high)
        if dry_partition(dataclasses.replace(ratios, fe=middle), feo, light, carbon_rock, graphite).metal < target:
            low = middle
        else:
            high = middle
    return dataclasses.replace(ratios, fe=0.5 * (low + high))


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


def origin_of(body):
    """
    `ICY` for a body that formed beyond its star's snow line, else `ROCKY`. A moon takes its parent's (the distance it formed at is
    `body.formation_distance_au`, the parent's for a moon); a body with no snow line to compare, or no distance, is rocky.
    """
    line = getattr(getattr(body, "star", None), "snow_line_au", None)
    distance = getattr(body, "formation_distance_au", None)
    if distance is None:
        distance = getattr(body, "distance", None)
    return ICY if line is not None and distance is not None and distance >= line else ROCKY


def _subject(body, ratios):
    return pt.Subject(origin_of(body), {"c_over_o": ratios.c_over_o})


def _snow_ratio(body):
    """The distance `body` formed at in its star's snow lines, or `None` when it has no star line or distance to say."""
    line = getattr(getattr(body, "star", None), "snow_line_au", None)
    distance = getattr(body, "formation_distance_au", None)
    if distance is None:
        distance = getattr(body, "distance", None)
    return None if not line or distance is None else distance / line


def _rocky_water(body):
    """A rocky body's water share from the Bern delivery law (`delivered_water`); a body with no mass yet is taken at the middle of the model's masses."""
    mass = getattr(body, "mass", None)
    mass_earth = _kg_to_mearth(mass) if mass else math.sqrt(math.prod(delivered_water.MASS_RANGE_EARTH))
    ratio = _snow_ratio(body)
    return delivered_water.sample_fraction(mass_earth, 0.5 if ratio is None else ratio)


def _icy_water_range(body):
    """
    The range of an icy body's water as a share of its mass: from `tuning.ICY_WATER_MASS_FRACTION_FLOOR` up to the share of the solids that
    is water at the disk temperature of its formation distance (`condensate_shares.water_share_of_solids`; the temperature falls as the
    square root of the distance, from the water condensation temperature at the star's own snow line).
    """
    distance = getattr(body, "formation_distance_au", None)
    if distance is None:
        distance = body.distance
    water_line = condensate_shares.ICE_CONDENSATION_K["H2O"]
    temperature = min(water_line * math.sqrt(body.star.snow_line_au / distance), math.nextafter(water_line, 0.0))
    top = condensate_shares.water_share_of_solids(temperature)
    return min(tuning.ICY_WATER_MASS_FRACTION_FLOOR, top), top


def _water_and_ice(body):
    """
    Water as a share of the whole body, the other ices beside it and those ices' shares by species (`ice_species`): the Bern delivery law for a rocky
    body, the icy range for one beyond the snow line, whose other ices follow the comet volatiles it draws. `(water, ice, shares)`; the shares are
    `None` for a body with no other ice.
    """
    if origin_of(body) == ICY:
        water = log_uniform(*_icy_water_range(body))
    else:
        water = _rocky_water(body)
    ice, shares = 0.0, None
    if origin_of(body) == ICY:
        volatiles = comet_population.draw_volatiles(RANGES["other_ice_comet_class"])
        ice = water * sum(volatiles.percent.values()) / 100.0
        shares = ice_species.shares_of_percent(volatiles.percent)
        scale = min(1.0, 0.95 / (water + ice))
        water, ice = water * scale, ice * scale
    return water, ice, shares


def _porosity(body):
    return draw.uniform(*RANGES["porosity"]) if origin_of(body) == ICY and _heap_mass(body) else None


SMALL_BODY_RADIUS_KM = 540.0
"""float: About the radius of a heap of `bulk_structure.SMALL_BODY_MASS_KG` of ice and rock, for a body given a radius and no mass."""


def _heap_mass(body):
    """Whether `body` is a heap of grains (lighter than `SMALL_BODY_MASS_KG`, or of a radius a heap that light has), so its porosity matters."""
    mass = getattr(body, "mass", None)
    if mass is not None:
        return mass < bulk_structure.SMALL_BODY_MASS_KG
    return body.radius is not None and body.radius < SMALL_BODY_RADIUS_KM


def mass_from_radius(bulk, radius_km, x_fe):
    """The mass in kg at which `bulk` has `radius_km`, by bisection on the structure; the radius it reaches if 20 Earth masses
    are not enough."""
    def radius_at(mass_kg):
        return bulk_structure.structure(bulk, mass_kg, x_fe=x_fe).radius_km

    top_kg = _mearth_to_kg(MAX_MASS_EARTH)
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
    kind = ROCK_KIND_TABLE.draw(origin_of(body), {"c_over_o": ratios.c_over_o}).name
    carbon_rock = kind != "silicate"
    inverse = radius_given and mass_given
    origin = None if inverse else getattr(body, "moon_origin", None)
    # A moon that holds its origin (GEN.466) takes its water, ice, rock and metal from it, not from the range of its parent's origin
    water, ice, ice_shares = (0.0, 0.0, None) if origin else _water_and_ice(body)
    feo = draw.uniform(*RANGES["rock_feo"])
    light = draw.uniform(*RANGES["core_light_elements"])
    porosity = None if origin else _porosity(body)
    graphite = draw.uniform(*RANGES["graphite_over_rock"]) if carbon_rock else 0.0

    if inverse:
        _draw_inverse(body, ratios, water, ice, ice_shares, feo, light, carbon_rock, graphite)
        return

    ratios = iron_branch_ratios(body, ratios, feo, light, carbon_rock, graphite)
    part = dry_partition(ratios, feo, light, carbon_rock, graphite)
    if origin:
        _draw_moon(body, origin, part, mass_given)
        return
    metal, rock = _normalise(part.metal, part.rock, water, ice)
    bulk = bulk_composition.BulkFractions(
        metal=metal, rock=rock, water=water, ice=ice, envelope=0.0, rock_feo=part.rock_feo, rock_al2o3=part.rock_al2o3,
        rock_mg_si=part.rock_mg_si, porosity=porosity, carbon_rock=carbon_rock, origin="drawn", ice_shares=ice_shares)
    _settle(body, bulk, mass_given)


def _draw_moon(body, origin, part, mass_given):
    """
    A moon's composition from its origin (GEN.466, `moon_composition`): the metal, rock, water and ice shares and the porosity of the
    origin's makeup, the rock's FeO, Al2O3 and Mg/Si from the star as for any rocky body. Below `compaction.TAPER_END_KG` the mass
    follows from the radius by the tapered brittle law of GEN.442; heavier, the structure solver decides, with no macroporosity.
    """
    system = getattr(body, "moon_disk_system", None)
    if origin == "disk" and system is None:    # a moon built on its own: the disk of a body of its parent's mass
        system = moon_composition.draw_disk_system(body.primary_mass_kg)
    makeup = moon_composition.draw_moon(origin, body.distance * constants.AU_TO_KM, system)
    bulk = bulk_composition.BulkFractions(
        metal=makeup.metal, rock=makeup.rock, water=makeup.water, ice=makeup.ice, envelope=0.0, rock_feo=part.rock_feo,
        rock_al2o3=part.rock_al2o3, rock_mg_si=part.rock_mg_si, porosity=makeup.porosity, carbon_rock=False, origin="drawn")
    if mass_given:
        # GEN.442: below the taper's end the radius follows from the mass by the same brittle law that gives the mass from the radius
        low, high = compaction.ROCK_FRACTION_RANGE
        try:
            if body.mass >= compaction.TAPER_END_KG:
                raise compaction.OutOfRange("the structure solver decides")
            radius_km = compaction.compacted_radius_km(body.mass, min(max(makeup.rock + makeup.metal, low), high))
        except compaction.OutOfRange:    # heavier than the taper's end (or lighter than its table): the solver
            _settle(body, dataclasses.replace(bulk, porosity=None), mass_given)
            return
        body.radius = radius_km
        volume_m3 = 4.0 / 3.0 * math.pi * (radius_km * 1000.0) ** 3
        grain = bulk_structure.mixing_density_kg_m3(bulk.metal, bulk.rock, bulk.water, bulk.ice, 0.0)   # what the stored shares hold with no pores
        _store(body, dataclasses.replace(bulk, porosity=max(0.0, 1.0 - body.mass / volume_m3 / grain)))
        return
    try:
        mass_kg = moon_composition.mass_from_radius_kg(makeup, body.radius)
    except compaction.OutOfRange:    # heavier than the taper's end: the solver, and the voids are closed
        _settle(body, dataclasses.replace(bulk, porosity=None), mass_given)
        return
    volume_m3 = 4.0 / 3.0 * math.pi * (body.radius * 1000.0) ** 3
    grain = makeup.grain_density_kg_m3
    body.mass = mass_kg
    _store(body, dataclasses.replace(bulk, porosity=max(0.0, 1.0 - mass_kg / volume_m3 / grain)))


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
    body.composition = composition_text.text(body)  # GEN.430: the text states the stored fractions
    volume_m3 = 4.0 / 3.0 * math.pi * (body.radius * 1000.0) ** 3
    body.density = body.mass / volume_m3 / 1000.0


def _draw_inverse(body, ratios, water, ice, ice_shares, feo, light, carbon_rock, graphite):
    """A body given both a mass and a radius: the core mass fraction from the mass and the radius."""
    part = dry_partition(ratios, feo, light, carbon_rock, graphite)
    x_fe = bulk_structure.x_fe_from_feo_mass_fraction(part.rock_feo) if part.rock_feo else 0.0
    mass_earth = _kg_to_mearth(body.mass)
    radius_earth = body.radius / constants.EARTH_RADIUS_KM
    wet = water + ice
    cmf = bulk_structure.cmf_from_mass_radius(mass_earth, radius_earth, x_fe, 1.0, high=MAX_CMF)
    if cmf is None:
        cmf = MAX_CMF if radius_earth < bulk_structure.radius_fit_earth(mass_earth, MAX_CMF, x_fe, 1.0) else 0.0
    bulk = bulk_composition.BulkFractions(
        metal=cmf * (1.0 - wet), rock=(1.0 - cmf) * (1.0 - wet), water=water, ice=ice, envelope=0.0, rock_feo=part.rock_feo,
        rock_al2o3=part.rock_al2o3, rock_mg_si=part.rock_mg_si, porosity=None, carbon_rock=carbon_rock, origin="inverse",
        ice_shares=ice_shares)
    _store(body, bulk)
