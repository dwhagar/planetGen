# planetgen/physics/surface_composition.py

"""
The surface and crust composition of a planet or moon (GEN.441, docs/design/composition/surface-composition.md)
==============================================================================================================

What lies on the ground, as shares of rows of the shared materials table (`physics.materials`, GEN.427): the minerals of the
crust, the ices that survive the surface temperature, the ocean and its salts. It reads the bulk fractions of GEN.428 (how
much ice, metal and carbon rock the body has), the hydrosphere of GEN.461 (what covers the surface: liquid, ice or nothing,
and how much of it), the surface temperature and pressure, the air's oxygen and the body's heat flow and age. Nothing here
is a Solar System preset: the ground is the one the body's own inputs allow.

**The draw, in order** (the order is part of the generator epoch):

1. *What covers the surface.* The hydrosphere gives a covered share `f`: liquid ocean (`liquid water` and its salts), or ice
   (`water ice`, with the other ices of the body that the temperature and the air let stand, `ice_stable`). A body that is
   ice-rich all through (water and other ices at least half the rock's mass) and frozen is ice from pole to pole.
2. *The dominant ground of the rest* (`DOMINANT_TABLE`, a gated possibility table): metal-rich where the mantle is stripped
   (a metal share of 0.7 or more, Hu et al. 2012), carbon where the rock is carbon rock, otherwise by the primary and
   secondary crust rules of Hu, Ehlmann and Seager 2012 (part E section 3.2): basaltic where decompression melting still
   resurfaces (the heat flow is over a floor `q_min`, log-uniform), else the primary crust of the magma ocean, ultramafic
   where the mantle overturned (more likely the more FeO) and anorthositic where it did not.
3. *Modifiers*, each an areal share of the land and drawn only when the inputs allow it: granitoid (`felsic`) with plate
   tectonics, clay where liquid water ever stood, iron-oxidised ground where the air holds free oxygen, evaporite where an
   ocean dried back, and sulfur where tidal heating beats the radiogenic flow.
4. *Grain state* (slab, crushed or powder): a surface under 0.1 Gyr old or molten is not yet regolith.

The assemblages are those of Hu et al. 2012 Table 2 (`GROUND_TYPES`); where the materials table has no row for one of their
minerals the nearest row stands in and the row says so. Every range that is a parametrisation rather than a measurement is
named in `RANGES` with what it stands in for; none is read from a Solar System body.

Use:

    composition = draw_surface(inputs)          # a SurfaceComposition, or None for a body with no solid ground
    composition.components                      # {materials row: mass share}, adds up to 1
    composition.mix                             # {ground type: areal share}, adds up to 1
    band_albedo(composition)                    # the albedo by band, area-weighted over the rows that have one
    crust_reservoir(composition)                # K, Th, U of the crust, for the radiogenic heat of GEN.459
"""

import math
from dataclasses import dataclass, field

from planetgen.generation import possibility_table as pt
from planetgen.physics import hydrosphere, materials, saturation
from planetgen.util import draw
from planetgen.util.random import log_uniform

SOURCES = {
    "hu2012": "Hu, Ehlmann and Seager 2012, ApJ 752, 7 (arXiv:1204.1544), Table 2: mineral assemblages of airless-planet crusts",
    "partE": "part E of the body composition research (docs/design/composition/body-composition-parts/E-ground-radiogenic-"
             "shielding.md) section 3.2, the ground-type decision rules",
    "palme2014": "Palme and O'Neill 2014 via pyrolite; Workman and Hart 2005; Gale et al. 2013; Rudnick and Gao 2014 "
                 "(part E table 2.8b, second-hand)",
    "millero": "mean ocean salinity 35 g/kg (Millero et al. 2008, Deep-Sea Research I 55, 50)",
}

GIANT_CLASSES = ("I", "J", "T")
"""tuple: The classes with no solid ground (every field empty)."""

RANGES = {
    "volcanism_floor_w_m2": (1.0e-4, 1.0e-3),
    "overturn_feo_scale": 0.20,
    "felsic_share": (0.0, 0.5),
    "clay_share": (0.0, 0.3),
    "oxidised_share": (0.05, 0.4),
    "evaporite_share": (0.0, 0.02),
    "sulfur_share": (0.05, 0.3),
    "young_surface_gy": 0.1,
    "ice_rich_over_rock": 0.5,
    "metal_stripped": 0.7,
    "free_oxygen_pa": 20.0,
    "molten_k": 1400.0,
    "ocean_salt": 0.035,
}
"""dict: The parametrisations of part E section 3.2 (none is a measurement; each is a placeholder for a fit against data).

- `volcanism_floor_w_m2`: the heat flow below which a body no longer melts and resurfaces, log-uniform (E 3.2 rule 4).
- `overturn_feo_scale`: the probability that the magma-ocean mantle overturned is the rock's FeO share over this (more FeO,
  denser late cumulates, E 3.2 rule 3).
- `felsic_share`: the granitoid share of the land of a body with plate tectonics, uniform (E 3.2 rule 5).
- `clay_share`, `oxidised_share`, `evaporite_share`, `sulfur_share`: the areal shares of the modifiers, uniform (E 3.2 rules 6, 7).
- `young_surface_gy`: a surface younger than this (the body's age) is not yet a regolith (E 3.2 rule 8).
- `ice_rich_over_rock`: water and other ices, over the rock's mass, from which a frozen body is ice at the surface.
- `metal_stripped`: the metal share of a body from which the mantle is taken to be stripped (Hu et al. 2012's metal-rich surface).
- `free_oxygen_pa`: the O2 partial pressure from which the ground oxidises, Pa (the ozone floor of 0.02 kPa, `physics.ozone`).
- `molten_k`: a surface at least this hot is molten (the lava-world threshold of `physics.mixture_draw`).
- `ocean_salt`: salts as a share of the ocean's mass (Earth's mean; no other world's is known).
"""

OCEAN, ICE = "ocean", "ice"
LIQUID_STATES = ("surface ocean", "hycean", "ice-covered ocean")
FROZEN_STATES = ("ice", "ice-covered ocean")
SLAB, CRUSHED, POWDER = "slab", "crushed", "powder"
GRAINS = (SLAB, CRUSHED, POWDER)

ICE_FORMULAS = {"water_ice": "H2O", "co2_ice": "CO2", "ch4_ice": "CH4", "nh3_ice": "NH3", "n2_ice": "N2", "so2_ice": "SO2",
                "co_ice": "CO", "hcn_ice": "HCN"}
"""dict: The ice rows of the materials table and the formula whose saturation curve says when each stands."""

VACUUM_ICE_PA = saturation.pressure_pa("H2O", hydrosphere.VACUUM_ICE_MAX_K)
"""float: The vapour pressure at which an ice survives in vacuum for billions of years: water's at the 220 K that
`hydrosphere` already uses, so every ice shares one criterion."""

# The assemblages of Hu et al. 2012 Table 2, on rows of the materials table. Hu's Mg-smectite has no row (kaolinite, the clay row
# with a measured albedo, stands in), its pyrite none (troilite and kamacite, the sulfide and the metal, stand in), its glass is
# the quenched melt (`basalt powder`) and its augite `diopside`; the others are the minerals Hu names.
GROUND_TYPES = pt.REGISTRY.register(pt.Table(
    "surface-ground-types", 1,
    [
        pt.Row("ultramafic", (pt.ANY_CLASS,), 1.0, SOURCES["hu2012"], fractions={"olivine": 0.60, "enstatite": 0.40},
               notes="primary crust after mantle overturn"),
        pt.Row("anorthositic", (pt.ANY_CLASS,), 1.0, SOURCES["hu2012"],
               fractions={"plagioclase": 0.97, "diopside": 0.03}, notes="primary flotation crust, no overturn"),
        pt.Row("basaltic", (pt.ANY_CLASS,), 1.0, SOURCES["hu2012"],
               fractions={"plagioclase": 0.76, "diopside": 0.08, "enstatite": 0.06, "basalt_powder": 0.05, "olivine": 0.01},
               notes="secondary crust of decompression melting; Hu's glass is `basalt powder`"),
        pt.Row("felsic", (pt.ANY_CLASS,), 1.0, SOURCES["hu2012"],
               fractions={"orthoclase": 0.40, "quartz": 0.35, "plagioclase": 0.20, "biotite": 0.05},
               notes="granitoid tertiary crust; Hu's K-feldspar is `orthoclase`"),
        pt.Row("clay", (pt.ANY_CLASS,), 1.0, SOURCES["hu2012"], state=pt.PARTLY_VERIFIED,
               limit="Hu's Mg-smectite has no row; kaolinite stands in",
               fractions={"kaolinite": 0.50, "serpentine": 0.50}, notes="aqueous alteration of basalt"),
        pt.Row("oxidised", (pt.ANY_CLASS,), 1.0, SOURCES["hu2012"],
               fractions={"hematite_fine": 0.50, "plagioclase": 0.38, "diopside": 0.04, "enstatite": 0.03,
                          "basalt_powder": 0.03, "olivine": 0.02}, notes="Hu's 50/50 nanophase hematite and basalt"),
        pt.Row("metal", (pt.ANY_CLASS,), 1.0, SOURCES["hu2012"], state=pt.PARTLY_VERIFIED,
               limit="Hu's pyrite has no row; troilite and kamacite stand in", fractions={"troilite": 0.5, "kamacite": 0.5},
               notes="a mantle stripped to metal and sulfide"),
        pt.Row("carbon", (pt.ANY_CLASS,), 1.0, SOURCES["partE"], state=pt.PARTLY_VERIFIED,
               limit="the carbon-rich crust rule was not reached in part E; graphite and the organics are the rows with albedo",
               fractions={"carbon_graphite": 0.40, "kerogen": 0.30, "anthracite": 0.30}, notes="carbon-rock ground"),
        pt.Row("evaporite", (pt.ANY_CLASS,), 1.0, SOURCES["partE"], fractions={"halite": 0.4, "gypsum": 0.3, "calcite": 0.3},
               state=pt.PARTLY_VERIFIED, limit="proportions are equal-ish placeholders until a brine evaporation model",
               notes="dried basins"),
        pt.Row("sulfur", (pt.ANY_CLASS,), 1.0, SOURCES["partE"], fractions={"sulfur": 1.0},
               notes="sulfur deposits of a tidally heated body"),
    ],
    components=pt.SHARED_COMPONENTS,
    description="The mineral assemblage of each kind of ground. Nothing is drawn from this table: it holds the recipes."))
"""Table: The recipes of the ground types, on rows of the shared materials table (weights are 1 and unused)."""

DOMINANT_TABLE = pt.REGISTRY.register(pt.Table(
    "surface-dominant-ground", 1,
    [
        pt.Row("metal", (pt.ANY_CLASS,), 1.0, SOURCES["hu2012"], state=pt.PARTLY_VERIFIED,
               limit="the stripped-mantle threshold is a parametrisation (RANGES)",
               conditions=(pt.Condition("metal_exposed", equals=1),), notes="a stripped mantle"),
        pt.Row("carbon", (pt.ANY_CLASS,), 1.0, SOURCES["partE"], state=pt.PARTLY_VERIFIED,
               limit="the carbon-rich crust rule was not reached in part E",
               conditions=(pt.Condition("metal_exposed", equals=0), pt.Condition("carbon_rock", equals=1))),
        pt.Row("basaltic", (pt.ANY_CLASS,), 1.0, SOURCES["partE"], state=pt.PARTLY_VERIFIED,
               limit="the heat-flow floor is a parametrisation",
               conditions=(pt.Condition("metal_exposed", equals=0), pt.Condition("carbon_rock", equals=0),
                           pt.Condition("volcanic", equals=1))),
        pt.Row("ultramafic", (pt.ANY_CLASS,), 1.0, SOURCES["partE"], state=pt.PARTLY_VERIFIED,
               limit="the overturn chance is a parametrisation",
               conditions=(pt.Condition("metal_exposed", equals=0), pt.Condition("carbon_rock", equals=0),
                           pt.Condition("volcanic", equals=0), pt.Condition("overturn", equals=1))),
        pt.Row("anorthositic", (pt.ANY_CLASS,), 1.0, SOURCES["partE"], state=pt.PARTLY_VERIFIED,
               limit="the overturn chance is a parametrisation",
               conditions=(pt.Condition("metal_exposed", equals=0), pt.Condition("carbon_rock", equals=0),
                           pt.Condition("volcanic", equals=0), pt.Condition("overturn", equals=0))),
    ],
    description="The ground type that covers a body's land, by whether its mantle is stripped, its rock is carbon rock, "
                "decompression melting still resurfaces it and the magma ocean overturned."))
"""Table: The dominant ground of the land; exactly one row passes."""

# K, Th and U (ppm by mass) of the reservoirs that model a ground type (part E table 2.8b; second-hand, via pyrolite).
RESERVOIRS = {
    "ultramafic": {"k_ppm": 51.0, "th_ppm": 0.0079, "u_ppm": 0.0032, "model": "depleted mantle (Workman and Hart 2005)"},
    "basaltic": {"k_ppm": 1162.0, "th_ppm": 0.252, "u_ppm": 0.083, "model": "N-MORB (Gale et al. 2013)"},
    "felsic": {"k_ppm": 23244.0, "th_ppm": 10.5, "u_ppm": 2.7, "model": "upper continental crust (Rudnick and Gao 2014)"},
}
"""dict: The K, Th and U of the ground types that have a measured reservoir model. The others (anorthositic, clay, oxidised,
metal, carbon, evaporite, sulfur, ice) have none in the sources read and are left to the caller's default."""


@dataclass(frozen=True)
class SurfaceInputs:
    """
    What the surface draw reads (the stored values of a planet or moon).

    Attributes:
        planet_class (str): The class letter.
        surface_temperature_k, pressure_pa (float): The surface temperature and the air pressure (0 for none).
        o2_partial_pa (float): The free oxygen's partial pressure.
        partial_pressures_pa (dict): Other gases' partial pressures by formula (`CO2`, `N2`, ...), for the ices that stand.
        water_mass_fraction, ice_mass_fraction, rock_mass_fraction, metal_mass_fraction (float): The bulk fractions (GEN.428).
        rock_feo_fraction (float or None): FeO as a share of the rock.
        carbon_rock (bool): The rock is silicon carbide and graphite.
        volatile_shares (dict): The other ices' mass shares by materials row (`co2 ice`, ...), adding up to the
            non-water ice; empty means water only.
        hydrosphere (str or None): The hydrosphere state (`hydrosphere.HYDROSPHERE_STATES`).
        ocean_fraction (float): The share of the surface the water covers, liquid or frozen.
        ocean_class (str or None): The ocean's chemistry class.
        internal_heat_flux_w_m2 (float): The body's own heat flow.
        tidal_heat_flux_w_m2 (float): The tidal heat flow (0 for none).
        age_gy (float): The body's age.
        mobile_lid_chance (float): The chance of plate tectonics (`magnetism.mobile_lid_chance`).
    """

    planet_class: str
    surface_temperature_k: float
    pressure_pa: float = 0.0
    o2_partial_pa: float = 0.0
    partial_pressures_pa: dict = field(default_factory=dict)
    water_mass_fraction: float = 0.0
    ice_mass_fraction: float = 0.0
    rock_mass_fraction: float = 1.0
    metal_mass_fraction: float = 0.0
    rock_feo_fraction: float | None = None
    carbon_rock: bool = False
    volatile_shares: dict = field(default_factory=dict)
    hydrosphere: str | None = None
    ocean_fraction: float = 0.0
    ocean_class: str | None = None
    internal_heat_flux_w_m2: float = 0.0
    tidal_heat_flux_w_m2: float = 0.0
    age_gy: float = 4.5
    mobile_lid_chance: float = 0.0


@dataclass(frozen=True)
class SurfaceComposition:
    """
    A body's surface composition.

    Attributes:
        ground_type (str): The ground type with the largest areal share (`ocean` and `ice` included).
        mix (dict): Areal share by ground type (`ocean`, `ice` and the rows of `GROUND_TYPES`), adding up to 1.
        components (dict): Mass share by materials row, adding up to 1 (a share of the surface is taken as a share of its mass).
        grain (str): `slab`, `crushed` or `powder`.
    """

    ground_type: str
    mix: dict
    components: dict
    grain: str


def ice_stable(material_name, temperature_k, pressure_pa, partial_pa=0.0):
    """
    Whether the ice `material_name` stands on a surface at `temperature_k`.

    Water ice stands under 273.15 K where the air holds at least water's triple-point pressure, and where it does not under the
    220 K that `hydrosphere` uses for exposed ice in vacuum; any other ice stands where its vapour pressure is no higher than its
    own partial pressure in the air or, with less than that, the vacuum floor `VACUUM_ICE_PA`.
    """
    formula = ICE_FORMULAS[material_name]
    if formula == "H2O":
        if pressure_pa >= saturation.WATER_TRIPLE_POINT_PA:
            return temperature_k < hydrosphere.WATER_MELTING_K
        return temperature_k < hydrosphere.VACUUM_ICE_MAX_K
    if not saturation.has(formula):
        return False
    return saturation.pressure_pa(formula, temperature_k) <= max(partial_pa, VACUUM_ICE_PA)


def _frozen_components(inputs):
    """The ices that stand at the surface, as mass shares adding up to 1 (water ice first); empty if none stands."""
    other = {m: s for m, s in inputs.volatile_shares.items() if s > 0.0}
    water = max(inputs.water_mass_fraction, 0.0)
    ices = {}
    if water > 0.0 and ice_stable("water_ice", inputs.surface_temperature_k, inputs.pressure_pa):
        ices["water_ice"] = water
    total_other = math.fsum(other.values())
    scale = inputs.ice_mass_fraction / total_other if total_other > 0.0 else 0.0
    for name, share in other.items():
        partial = inputs.partial_pressures_pa.get(ICE_FORMULAS[name], 0.0)
        if ice_stable(name, inputs.surface_temperature_k, inputs.pressure_pa, partial):
            ices[name] = share * scale
    total = math.fsum(ices.values())
    return {name: share / total for name, share in ices.items()} if total > 0.0 else {}


def _quantities(inputs, frozen):
    """The stored quantities the dominant-ground table's gates read, with the two draws that decide them."""
    floor = log_uniform(*RANGES["volcanism_floor_w_m2"])
    feo = inputs.rock_feo_fraction or 0.0
    overturn = draw.random() < min(1.0, feo / RANGES["overturn_feo_scale"])
    stripped = inputs.metal_mass_fraction >= RANGES["metal_stripped"] and not inputs.carbon_rock
    return {"metal_exposed": int(stripped), "carbon_rock": int(bool(inputs.carbon_rock)),
            "volcanic": int(inputs.internal_heat_flux_w_m2 >= floor), "overturn": int(overturn)}


def _covered(inputs, frozen):
    """The covered share of the surface and its kind: `(share, OCEAN | ICE | None)`."""
    share = min(max(inputs.ocean_fraction or 0.0, 0.0), 1.0)
    if share <= 0.0:
        return 0.0, None
    if inputs.hydrosphere in ("surface ocean", "hycean"):
        return share, OCEAN
    if inputs.hydrosphere in FROZEN_STATES and frozen:
        return share, ICE
    return 0.0, None


def _ice_rich(inputs, frozen):
    """Whether a frozen body is ice through and through, so its whole surface is ice."""
    if not frozen:
        return False
    rock = max(inputs.rock_mass_fraction, 1e-9)
    return (max(inputs.water_mass_fraction, 0.0) + max(inputs.ice_mass_fraction, 0.0)) >= RANGES["ice_rich_over_rock"] * rock


def _grain(inputs):
    if inputs.surface_temperature_k >= RANGES["molten_k"]:
        return SLAB
    if inputs.age_gy < RANGES["young_surface_gy"]:
        return draw.choice((SLAB, CRUSHED))
    return POWDER


def _modifiers(inputs, covered_kind):
    """The modifier shares of the land, drawn in a fixed order: `{ground type: share of the land}`."""
    shares = {}
    lid = inputs.mobile_lid_chance > 0.0 and draw.random() < inputs.mobile_lid_chance
    if lid:
        shares["felsic"] = draw.uniform(*RANGES["felsic_share"])
    water_ever = covered_kind is not None or inputs.hydrosphere in LIQUID_STATES
    if water_ever:
        shares["clay"] = draw.uniform(*RANGES["clay_share"])
    if inputs.o2_partial_pa >= RANGES["free_oxygen_pa"]:
        shares["oxidised"] = draw.uniform(*RANGES["oxidised_share"])
    if covered_kind == OCEAN:
        shares["evaporite"] = draw.uniform(*RANGES["evaporite_share"])
    if inputs.tidal_heat_flux_w_m2 > inputs.internal_heat_flux_w_m2 > 0.0 or \
            (inputs.tidal_heat_flux_w_m2 > 0.0 and inputs.internal_heat_flux_w_m2 == 0.0):
        shares["sulfur"] = draw.uniform(*RANGES["sulfur_share"])
    return shares


def draw_surface(inputs):
    """
    Draws the surface composition of a body (see the module text for the order).

    Args:
        inputs (SurfaceInputs): The body's stored values.

    Returns:
        SurfaceComposition or None: `None` for a giant, which has no solid ground.
    """
    if inputs.planet_class in GIANT_CLASSES:
        return None
    frozen = _frozen_components(inputs)
    quantities = _quantities(inputs, frozen)
    covered, kind = _covered(inputs, frozen)
    if _ice_rich(inputs, frozen):
        covered, kind = 1.0, ICE
    mix = {}
    if kind is not None:
        mix[kind] = covered
    land = 1.0 - covered
    if land > 1e-12:
        dominant = DOMINANT_TABLE.draw(inputs.planet_class, quantities).name
        mods = _modifiers(inputs, kind)
        if dominant in ("metal", "carbon"):
            mods.pop("felsic", None)
        total_mods = math.fsum(mods.values())
        if total_mods > 0.9:
            mods = {k: v * 0.9 / total_mods for k, v in mods.items()}
            total_mods = 0.9
        mix[dominant] = mix.get(dominant, 0.0) + land * (1.0 - total_mods)
        for name, share in mods.items():
            mix[name] = mix.get(name, 0.0) + land * share
    components = _components(mix, frozen, inputs)
    ground_type = max(mix, key=lambda name: (mix[name], name))
    return SurfaceComposition(ground_type=ground_type, mix=mix, components=components, grain=_grain(inputs))


def _components(mix, frozen, inputs):
    """The mass shares by materials row of the areal `mix`."""
    parts = {}

    def add(name, share):
        parts[name] = parts.get(name, 0.0) + share

    for kind, areal in mix.items():
        if kind == OCEAN:
            salt = RANGES["ocean_salt"]
            add("water_liquid", areal * (1.0 - salt))
            add("halite", areal * salt)
        elif kind == ICE:
            for name, share in frozen.items():
                add(name, areal * share)
        else:
            for name, share in GROUND_TYPES.row(kind).normalised_fractions().items():
                add(name, areal * share)
    total = math.fsum(parts.values())
    return {name: share / total for name, share in parts.items()}


ALBEDO_ROWS = {
    "ultramafic": {"ultramafic_dunite_powder": 0.5, "ultramafic_pyroxenite_powder": 0.5},
    "anorthositic": {"anorthosite_powder": 1.0},
    "basaltic": {"basalt_powder": 1.0},
    "felsic": {"felsic_granite_powder_a": 0.5, "felsic_granite_powder_b": 0.5},
    "clay": {"kaolinite": 0.5, "serpentine": 0.5},
    "oxidised": {"hematite_fine": 1.0},
    "carbon": {"kerogen": 0.5, "carbon_black": 0.5},
    "evaporite": {"halite": 0.4, "gypsum": 0.3, "calcite": 0.3},
    "sulfur": {"sulfur": 1.0},
    ICE: {"water_ice": 1.0},
}
"""dict: The lab-spectrum rows of the materials table that stand for a ground type's colour (part E table 2.6; the assemblage
minerals have no spectra of their own). A grain state of slab or crushed swaps the basalt powder for the slab. The ocean and the
metal have no measured visible spectrum in the sources read, so they are not covered (`band_albedo` reports the covered share)."""

BANDS = ("alb_040_050", "alb_050_060", "alb_060_070")


def band_albedo(composition):
    """
    The lab reflectance by band of the surface (`BANDS`), area-weighted over the ground types that have a measured spectrum
    (`ALBEDO_ROWS`), and the share of the surface they cover. It is a reflectance proxy, not a Bond albedo: the weathering and
    the ocean belong to the albedo of GEN.462.

    Returns:
        tuple: `({band: reflectance or None}, covered areal share)`.
    """
    covered = math.fsum(share for name, share in composition.mix.items() if name in ALBEDO_ROWS)
    result = {}
    for band in BANDS:
        total, weight = 0.0, 0.0
        for name, areal in composition.mix.items():
            rows = dict(ALBEDO_ROWS.get(name, {}))
            if name == "basaltic" and composition.grain != POWDER:
                rows = {"basalt_slab": 1.0}
            for row, part in rows.items():
                value = materials.material(row).prop(band).value
                if value is not None:
                    total += areal * part * value
                    weight += areal * part
        result[band] = total / weight if weight > 0.0 else None
    return result, covered


def crust_reservoir(composition):
    """
    The K, Th and U (ppm by mass) of the crust, area-weighted over the land types that have a measured reservoir
    (`RESERVOIRS`), and the share of the land they cover.

    Returns:
        tuple: `({k_ppm, th_ppm, u_ppm} or None, covered share of the land)`.
    """
    land = {name: share for name, share in composition.mix.items() if name not in (OCEAN, ICE)}
    total_land = math.fsum(land.values())
    known = {name: share for name, share in land.items() if name in RESERVOIRS}
    weight = math.fsum(known.values())
    if total_land <= 0.0 or weight <= 0.0:
        return None, 0.0
    result = {key: math.fsum(RESERVOIRS[name][key] * share for name, share in known.items()) / weight
              for key in ("k_ppm", "th_ppm", "u_ppm")}
    return result, weight / total_land
