# planetgen/physics/evaporite_sequence.py

"""
The salts a brine leaves when it dries, in the order its chemistry lets them go (GEN.1366)
==========================================================================================

`surface_composition` gave the evaporite ground a fixed recipe (four tenths halite, three tenths gypsum, three tenths calcite) and the ocean's dissolved salt one row,
halite. A brine is not that: what it leaves depends on what it holds, and the order is the one fact of evaporite chemistry that does not change from sea to lake. This
module allocates a brine's ions to the solids of the materials table, in that order, and picks the hydrate or polymorph of each by temperature and water activity.

**The ions** (`brine_ions`). The rock that weathers is the volcanic crust of the body (the melt of its mantle, `rock_types.partial_melt`); each of its cations dissolves in
proportion to its oxide moles times the mobility of the weathering series (`MOBILITY`: Ca, Na, Mg, K in the order of Polynov's series). The anions are the ocean class's:
`ANION_PARAMS` gives, for each class of `hydrosphere.OCEAN_CLASSES`, the carbonate per calcium and the chloride share of the rest, and the total anion charge
is the cation charge (a brine is neutral). Nothing is a named world: the class is a regime (a chloride brine, an acid sulfate water, a soda water, a neutral one).

**The sequence** (`allocate`). The chemical divides of Hardie and Eugster fix the order; every step takes what the ions left by the steps before it allow:

1. *Acid water only*: jarosite then alunite take potassium and sulfate (they need iron or aluminium from the acid-leached rock, which the rock oxides give).
2. *Calcite* takes calcium and carbonate in equal moles. What is left over decides the path: carbonate left over is the alkaline path (soda), calcium left over the calcium
   path.
3. *Alkaline path*: the carbonate left takes magnesium (hydromagnesite, magnesite from `MAGNESITE_MIN_K`) and then sodium, as nahcolite where the carbon dioxide is high
   (`NAHCOLITE_MIN_PCO2_PA`) and trona where it is low, or the hydrate of sodium carbonate the temperature allows (`sodium_carbonate_row`).
4. *Calcium path*: calcium sulfate (`calcium_sulfate_row`: gypsum or anhydrite from [R5143]), then, if calcium is still left, calcium chloride; sulfate left takes magnesium
   (`magnesium_sulfate_row`: epsomite, hexahydrite or kieserite), then sodium (`sodium_sulfate_row`: mirabilite or thenardite).
5. *Chlorides*: perchlorate (magnesium, then calcium) where the air supplies the oxidant and no ozone shield stands, halite, carnallite, sylvite, bischofite and last
   antarcticite.

A mole of a solid is its formula's mole; the mass share is mole times molar mass. The ions that no step can pair (an acid water's free acid, a cation with no anion left) are
reported as `unplaced` and carry no mass.

**The calcium sulfate** is the only boundary read in a source: [R5143] measured the gypsum and anhydrite equilibrium as a function of temperature and the activity of
water (`GYPSUM_ANHYDRITE_AW`), and a brine saturated with halite has an activity of water of 0.75, at which gypsum dehydrates above 18 degrees C. The other boundaries are
stand-ins (registry section 6.212): the mirabilite and thenardite transition, the epsomite, hexahydrite and kieserite ones, the sodium carbonate hydrates, the nahcolite and
trona boundary, and the anions of the classes.

Limits: the sequence is equilibrium on the complete desiccation of a brine, not the kinetics of a basin (a real sequence stops, mixes and redissolves); the double salts
(glauberite, polyhalite, langbeinite, kainite) and the magnesium chlorides below bischofite are not formed; bassanite, which forms in concentrated chloride brines
instead of anhydrite, is not formed; the activity of water is the halite one where halite forms and the gypsum stage's where it does not, not a curve along the evaporation.
"""

import math
from dataclasses import dataclass, field

from planetgen.physics import elements

OXIDE_CATIONS = {"CaO": ("Ca", 2, 1), "MgO": ("Mg", 2, 1), "Na2O": ("Na", 1, 2), "K2O": ("K", 1, 2)}
"""dict: Oxide -> (cation, its charge, cations per formula unit)."""

OXIDE_PER_ELEMENT = {oxide: elements.molar_mass(oxide) / (per_unit * elements.molar_mass(cation)) for oxide, (cation, _z, per_unit) in OXIDE_CATIONS.items()}
"""dict: Oxide -> its mass per unit mass of its cation, dimensionless (`Na2O` 1.348 per sodium, `K2O` 1.205 per potassium): how an element's weight percent is stated as its oxide's."""

MOBILITY = {"Ca": 3.0, "Na": 2.4, "Mg": 1.3, "K": 1.25}
"""dict: The relative mobility of each cation in weathering, dimensionless: Polynov's series as Perel'man tabulates it (Ca 3.0, Na 2.4, Mg 1.3, K 1.25; Si 0.20, Fe 0.04, Al 0.02).
Recalled from the series, not opened (registry row EV1)."""

ANION_PARAMS = {
    "chloride brine": (0.10, 0.95),
    "acid sulfate": (0.0, 0.20),
    "soda": (1.5, 0.60),
    "neutral": (1.0, 0.85),
}
"""dict: Ocean class -> `(carbonate moles per mole of calcium, chloride share of the anion charge that carbonate leaves)`, dimensionless (registry row EV2). A neutral water holds
carbonate enough to take its calcium (the carbon dioxide of a living-range air is the calcite sink), a soda water more (the alkaline path), a chloride brine hardly any (calcium is
left for the chloride), an acid water none; sulfate is the rest of the charge. `ice-sealed` and a body with no class take `neutral`'s."""
DEFAULT_CLASS = "neutral"

# --- the calcium sulfate, from [R5143] ----------------------------------------------------------------------------------------------

GYPSUM_ANHYDRITE_AW = (
    (285.15, 0.725), (291.15, 0.750), (296.15, 0.770), (303.15, 0.800), (312.15, 0.845), (318.15, 0.880), (323.15, 0.915), (328.15, 0.960), (331.15, 1.000),
)
"""tuple: The equilibrium of gypsum and anhydrite, `(temperature K, activity of water)`: Hardie's table of the gypsum-anhydrite equilibrium ([R5143: p. 181], 12 to 58 degrees C;
12, 18 and 58 are extrapolated in the source, 23, 39 and 55 measured). Gypsum is the stable form above the line (more water), anhydrite below it."""
HALITE_SATURATED_AW = 0.75
"""float: The activity of water of a brine saturated with halite, dimensionless: the lower end of the activity [R5143] gives for the calcium sulfate of natural brines (the assemblage gypsum or anhydrite with halite and solution)."""
GYPSUM_STAGE_AW = 0.94
"""float: The activity of water of the brine at which calcium sulfate first precipitates, before halite, dimensionless (registry row EV3)."""

# --- the other boundaries (stand-ins, registry section 6.212) -----------------------------------------------------------------------

SODIUM_SULFATE_K = 305.55
"""float: Below it the sodium sulfate that dries is mirabilite (the decahydrate), above it thenardite (the anhydrous salt), K (32.4 degrees C, registry row EV4)."""
MAGNESIUM_SULFATE_K = (321.15, 341.15)
"""tuple: The magnesium sulfate that dries is epsomite below the first temperature, hexahydrite below the second, kieserite above, K (48 and 68 degrees C, registry row EV5)."""
SODIUM_CARBONATE_K = (305.65, 382.15)
"""tuple: The sodium carbonate that dries is the decahydrate below the first temperature (natron), the monohydrate below the second (thermonatrite), anhydrous natrite above, K
(32.5 and 109 degrees C, registry row EV6)."""
NAHCOLITE_MIN_PCO2_PA = 5.0e4
"""float: Above this carbon dioxide partial pressure the sodium that is left with carbonate dries as the bicarbonate (nahcolite), below it as trona, Pa (registry row EV7)."""
MAGNESITE_MIN_K = 333.15
"""float: Below this temperature the magnesium carbonate that dries is the hydrated one (hydromagnesite), above it magnesite, K (60 degrees C, registry row EV8)."""
PERCHLORATE_SHARE = 0.01
"""float: The share of the chloride charge that is oxidised to perchlorate where the air supplies oxidants and no ozone shield stands, dimensionless (registry row EV9)."""
PERCHLORATE_MIN_PA = 1.0
"""float: The least surface pressure with oxidants for the share above, Pa."""

JAROSITE_SHARE = 0.5
"""float: The share of an acid water's potassium that jarosite takes (the rest goes to alunite), dimensionless: the acid-leached rock holds more iron than aluminium, but the aluminium of the
feldspars is the larger reservoir at low pH (registry row EV10)."""

# --- the solids ---------------------------------------------------------------------------------------------------------------------

FORMULAS = {
    "calcite": "CaCO3", "hydromagnesite": "Mg5(CO3)4(OH)2.4H2O", "magnesite": "MgCO3", "natron": "Na2CO3.10H2O", "thermonatrite": "Na2CO3.H2O", "natrite": "Na2CO3",
    "nahcolite": "NaHCO3", "trona": "Na3H(CO3)2.2H2O", "gypsum": "CaSO4.2H2O", "anhydrite": "CaSO4", "epsomite": "MgSO4.7H2O", "hexahydrite": "MgSO4.6H2O",
    "kieserite": "MgSO4.H2O", "mirabilite": "Na2SO4.10H2O", "thenardite": "Na2SO4", "halite": "NaCl", "sylvite": "KCl", "carnallite": "KMgCl3.6H2O",
    "bischofite": "MgCl2.6H2O", "antarcticite": "CaCl2.6H2O", "mg_perchlorate": "Mg(ClO4)2", "ca_perchlorate": "Ca(ClO4)2", "jarosite": "KFe3(SO4)2(OH)6",
    "alunite": "KAl3(SO4)2(OH)6",
}
"""dict: The materials rows this sequence makes, with the formula whose molar mass turns moles into mass."""

ROW_IDS = tuple(FORMULAS)
"""tuple: Every row the sequence can return."""

_MOLAR_MASS = {row: elements.molar_mass(formula) for row, formula in FORMULAS.items()}


def gypsum_anhydrite_aw(temperature_k):
    """The activity of water at which gypsum and anhydrite are in equilibrium at `temperature_k` [R5143]: the table's line, extended below 12 degrees C by its lowest
    slope and held at 1.0 above 58 (past it anhydrite is stable in water of any activity)."""
    points = GYPSUM_ANHYDRITE_AW
    if temperature_k >= points[-1][0]:
        return 1.0
    if temperature_k <= points[0][0]:
        slope = (points[1][1] - points[0][1]) / (points[1][0] - points[0][0])
        return max(0.0, points[0][1] + slope * (temperature_k - points[0][0]))
    for (t0, a0), (t1, a1) in zip(points, points[1:]):
        if t0 <= temperature_k <= t1:
            return a0 + (a1 - a0) * (temperature_k - t0) / (t1 - t0)
    raise AssertionError("unreachable")


def calcium_sulfate_row(temperature_k, water_activity):
    """`gypsum` where the brine's activity of water is above the equilibrium line at `temperature_k`, `anhydrite` where it is at or below it [R5143]."""
    return "gypsum" if water_activity > gypsum_anhydrite_aw(temperature_k) else "anhydrite"


def sodium_sulfate_row(temperature_k):
    """`mirabilite` below `SODIUM_SULFATE_K`, else `thenardite`."""
    return "mirabilite" if temperature_k < SODIUM_SULFATE_K else "thenardite"


def magnesium_sulfate_row(temperature_k):
    """`epsomite`, `hexahydrite` or `kieserite` by temperature (`MAGNESIUM_SULFATE_K`)."""
    low, high = MAGNESIUM_SULFATE_K
    return "epsomite" if temperature_k < low else "hexahydrite" if temperature_k < high else "kieserite"


def sodium_carbonate_row(temperature_k, pco2_pa):
    """`nahcolite` where the carbon dioxide is high, else `natron`, `thermonatrite` or `natrite` by temperature (`SODIUM_CARBONATE_K`) or, between them, `trona`."""
    if pco2_pa >= NAHCOLITE_MIN_PCO2_PA:
        return "nahcolite"
    low, high = SODIUM_CARBONATE_K
    if temperature_k < low:
        return "trona"
    return "thermonatrite" if temperature_k < high else "natrite"


def magnesium_carbonate_row(temperature_k):
    """`hydromagnesite` below `MAGNESITE_MIN_K`, else `magnesite`."""
    return "hydromagnesite" if temperature_k < MAGNESITE_MIN_K else "magnesite"


@dataclass(frozen=True)
class Sequence:
    """
    What a brine left.

    Attributes:
        moles (dict): Materials row -> moles of its formula, in the units of the ions given (relative; only the ratios matter).
        shares (dict): Materials row -> mass share, adding to 1 (empty where nothing precipitated).
        unplaced (dict): Ion -> the charge no step could pair, in the units of the ions given.
        calcium_sulfate (str or None): `gypsum` or `anhydrite`, whichever the brine left.
    """

    moles: dict
    shares: dict
    unplaced: dict = field(default_factory=dict)
    calcium_sulfate: str | None = None


def brine_ions(oxides, ocean_class):
    """
    The ions of the brine a weathered crust makes in an ocean of a class, relative moles of charge-balanced ions.

    Args:
        oxides (dict): The crust's oxides, weight percent (`CaO`, `MgO`, `Na2O`, `K2O` are read).
        ocean_class (str or None): One of `hydrosphere.OCEAN_CLASSES`; anything else takes `DEFAULT_CLASS`.

    Returns:
        dict: `Ca`, `Mg`, `Na`, `K`, `Cl`, `SO4`, `CO3`: the moles of each ion, with the anion charge equal to the cation charge; all zero for a rock with none of the four oxides.
    """
    cations, charge = {}, 0.0
    for oxide, (cation, z, per_unit) in OXIDE_CATIONS.items():
        moles = max(oxides.get(oxide, 0.0), 0.0) / (elements.molar_mass(oxide) * 1000.0) * per_unit * MOBILITY[cation]
        cations[cation] = moles
        charge += moles * z
    carbonate_per_calcium, chloride_share = ANION_PARAMS.get(ocean_class, ANION_PARAMS[DEFAULT_CLASS])
    carbonate = min(carbonate_per_calcium * cations["Ca"], charge / 2.0)
    rest = charge - 2.0 * carbonate
    return {**cations, "Cl": rest * chloride_share, "SO4": rest * (1.0 - chloride_share) / 2.0, "CO3": carbonate}


def allocate(ions, temperature_k, pco2_pa=0.0, pressure_pa=0.0, o2_shield=False, acid=False):
    """
    Allocates the ions of a brine to solids in the order of the module text.

    Args:
        ions (dict): `brine_ions`'s result.
        temperature_k (float): The temperature the salts dry at, K.
        pco2_pa (float): The carbon dioxide partial pressure, Pa.
        pressure_pa (float): The surface pressure, Pa (the perchlorate needs the air's oxidants).
        o2_shield (bool): The air holds enough free oxygen for an ozone shield (no perchlorate then).
        acid (bool): The water is an acid sulfate one (jarosite and alunite come first).

    Returns:
        Sequence: The moles, mass shares and unplaced charge.
    """
    left = {name: max(float(ions.get(name, 0.0)), 0.0) for name in ("Ca", "Mg", "Na", "K", "Cl", "SO4", "CO3")}
    solids = {}

    def take(row, limit=math.inf, **per_formula):
        """Precipitates the most of `row` the ions allow, at most `limit` formulas (each keyword the moles of that ion per formula), and returns the moles."""
        count = min(limit, *(left[ion] / need for ion, need in per_formula.items() if need > 0.0))
        if count <= 0.0:
            return 0.0
        for ion, need in per_formula.items():
            left[ion] -= need * count
        solids[row] = solids.get(row, 0.0) + count
        return count

    if acid:
        # iron and aluminium come from the acid-leached rock; potassium is the limit (one per formula), sulfate two per formula
        take("jarosite", limit=JAROSITE_SHARE * left["K"], K=1.0, SO4=2.0)
        take("alunite", K=1.0, SO4=2.0)
    take("calcite", Ca=1.0, CO3=1.0)
    if left["CO3"] > 0.0:                                                 # the alkaline path: carbonate outlasts calcium
        magnesium = magnesium_carbonate_row(temperature_k)
        if magnesium == "hydromagnesite":
            take("hydromagnesite", Mg=5.0, CO3=4.0)
        else:
            take("magnesite", Mg=1.0, CO3=1.0)
        sodium_carbonate = sodium_carbonate_row(temperature_k, pco2_pa)
        if sodium_carbonate == "nahcolite":
            take("nahcolite", Na=1.0, CO3=0.5)                            # one sodium and one bicarbonate (half a carbonate's charge) per formula
        elif sodium_carbonate == "trona":
            take("trona", Na=3.0, CO3=2.0)
        else:
            take(sodium_carbonate, Na=2.0, CO3=1.0)
    sulfate_form = None
    if left["SO4"] > 0.0 and left["Ca"] > 0.0:                            # the calcium path: calcium outlasts carbonate
        halite_forms = left["Na"] > 0.0 and left["Cl"] > 0.0
        water_activity = HALITE_SATURATED_AW if halite_forms else GYPSUM_STAGE_AW
        sulfate_form = calcium_sulfate_row(temperature_k, water_activity)
        take(sulfate_form, Ca=1.0, SO4=1.0)
    take(magnesium_sulfate_row(temperature_k), Mg=1.0, SO4=1.0)
    take(sodium_sulfate_row(temperature_k), Na=2.0, SO4=1.0)
    if pressure_pa >= PERCHLORATE_MIN_PA and not o2_shield and left["Cl"] > 0.0 and PERCHLORATE_SHARE > 0.0:
        oxidised = left["Cl"] * PERCHLORATE_SHARE
        left["Cl"] -= oxidised
        left["ClO4"] = oxidised
        take("mg_perchlorate", Mg=1.0, ClO4=2.0)
        take("ca_perchlorate", Ca=1.0, ClO4=2.0)
    take("halite", Na=1.0, Cl=1.0)
    take("carnallite", K=1.0, Mg=1.0, Cl=3.0)
    take("sylvite", K=1.0, Cl=1.0)
    take("bischofite", Mg=1.0, Cl=2.0)
    take("antarcticite", Ca=1.0, Cl=2.0)
    masses = {row: moles * _MOLAR_MASS[row] for row, moles in solids.items() if moles > 0.0}
    total = math.fsum(masses.values())
    shares = {row: mass / total for row, mass in masses.items()} if total > 0.0 else {}
    charge = {"Ca": 2.0, "Mg": 2.0, "Na": 1.0, "K": 1.0, "Cl": 1.0, "SO4": 2.0, "CO3": 2.0, "ClO4": 1.0}
    unplaced = {ion: moles * charge[ion] for ion, moles in left.items() if moles * charge[ion] > 1.0e-12 * max(1.0, sum(ions.values()))}
    return Sequence(moles=solids, shares=shares, unplaced=unplaced, calcium_sulfate=sulfate_form)


def salts_of(oxides, ocean_class, temperature_k, pco2_pa=0.0, pressure_pa=0.0, o2_shield=False):
    """
    The mass shares of the solids a body's weathered crust and ocean class leave on drying, `{materials row: share}` adding to 1, or `{}` where the crust has no cation to dissolve.

    Args:
        oxides (dict or None): The weathering crust's oxides, weight percent; `None` for a body with no silicate rock (the volcanic sodium chloride is then all there is, so halite).
        ocean_class (str or None): The ocean's class (`hydrosphere.OCEAN_CLASSES`).
        temperature_k (float): The surface temperature, K.
        pco2_pa, pressure_pa (float): The carbon dioxide partial pressure and the surface pressure, Pa.
        o2_shield (bool): Free oxygen enough for an ozone shield.
    """
    if oxides is None:
        return {"halite": 1.0}
    result = allocate(brine_ions(oxides, ocean_class), temperature_k, pco2_pa, pressure_pa, o2_shield, acid=ocean_class == "acid sulfate")
    return dict(result.shares)
