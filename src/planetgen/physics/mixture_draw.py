# planetgen/physics/mixture_draw.py

"""
The atmosphere mixture of a body, drawn from its regimes (GEN.451, GEN.479, GEN.454, GEN.499; docs/design/mixture-draw.md)
=======================================================================================================================

`physics/atmosphere_recipes.py` (GEN.478) draws which regime a body's air is in: the mantle redox, the volatile delivery, the
outgassed class, the retention class, the envelope, the giant, brown-dwarf and hot-Jupiter regimes, the lava-world vapour type, the
trace gases. Its rows say what the air is made of in words. This module turns the drawn regimes into mole fractions of the gases
of the table (`physics/gases.py`), by chemistry and not by a per-class template: the element budget of the body, the temperature
and pressure of its layer and the oxygen fugacity of its mantle go into the equilibrium solver (`physics/gas_equilibrium.py`,
over the NASA Glenn thermodynamics of `physics/species_thermo.py`), and the gas that comes out is the mixture. A class letter
is never an input to the chemistry: it only gates which rows of the recipe tables apply.

    draw = outgassed_mixture(Subject("A", {... the quantities required_inputs("A") names ...}))
    draw.mole_fractions       # {"H2O": 0.78, "CO2": 0.16, ...}, gas-table formulas, sums to 1
    draw.drawn                # {"a": "a6", "b": "b2", "c": "c3", ...}
    draw.redox_diw            # log10 fO2 above the iron-wuestite buffer (Hirschmann 2021 scale), None for a giant

What this module does not do: the total pressure of the air (the retention row `d` gives the regime; the pressure, the freeze-out
of what the surface temperature condenses and the sum of the partial pressures are the generator's, GEN.411 and GEN.422), the
clouds and hazes (GEN.414) and the ozone (GEN.413). Every draw goes through `util/draw.py`, in a fixed order of its own, after
the recipe draw: a change of the order changes every body and is a generator-epoch change.

Rocky bodies (this slice; the giant, brown-dwarf, hot-Jupiter, envelope and lava-world families join with GEN.454 and GEN.499)
-------------------------

The outgassed mixture is the gas in equilibrium with a mantle of a given oxygen fugacity (Schaefer and Fegley 2010; Gaillard and
Scaillet 2014; Ortenzi et al. 2020, part 05 section 3 of the atmosphere research), as the research recipe 3.5 states:

1. *Redox.* The magma-ocean stage's oxygen fugacity is the row `a1` of the recipes (Sossi et al. 2020: IW+0.5 for about an Earth
   mass, about IW-2 for a small body whose redox was set at core formation; Hirschmann 2021 scale, `physics/oxygen_fugacity.py`),
   and a body with long-lived geology moves to the later upper mantle `a2` over a seeded oxidation time, `OXIDATION_TIME_GY`. A
   body around a white dwarf takes `a4`, a small body `a3`. The mean redox is therefore a function of mass, age and the
   building blocks, never the modern Earth's constant.
2. *Inventory.* The delivered volatile class (`b1` dry, `b2` low, `b3` high, `b4` water-rich) and its continuous water fraction
   (`b5`) give the hydrogen and carbon in the outgassing stage in the proportions of Ortenzi et al. 2020 Table 1; the
   nitrogen, the sulfur and the radiogenic argon come from the design brackets of `INVENTORY`.
3. *The volcano's gas.* The gas in equilibrium with the mantle at the magmatic temperature (`MAGMATIC_TEMPERATURE_K`, the
   research's 1200 K illustration widened) and the pressure of the layer, with the oxygen fugacity held at the mantle's: H2, H2O,
   CO, CO2, CH4, N2, NH3, H2S, SO2, COS and the rest of the table.
4. *Cooling.* The gas then cools in a closed system (the oxygen it carries is all the oxygen it has) to the surface temperature
   when that is above the quench temperature (`QUENCH_TEMPERATURE_K`, seeded between 400 and 900 K: Zahnle and Marley 2014 put
   the methane to carbon monoxide quench of planets at 400 to 600 K, Liggins et al. 2022 compute at 800 K), and to the quench
   temperature when it is colder: the reactions freeze there. Water, sulfur and the other condensables that a cool surface
   takes out of the gas are the pressure bookkeeping's (GEN.422); this mixture is the gas before that.

What the research does not give is drawn from the brackets below and labelled as the design's: each is a stand-in for a research
item (GEN.481 volatile delivery, GEN.482 redox, GEN.492 radioisotopes) that replaces it.
"""

import math
from dataclasses import dataclass, field

from planetgen.generation.possibility_table import Subject
from planetgen.physics import atmosphere_recipes, gas_equilibrium, gases, species_thermo
from planetgen.physics import oxygen_fugacity as fo2
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.util import draw

# --- constants -------------------------------------------------------------------------------------------------------------

AMU_G_MOL = {"H2O": 18.01528, "CO2": 44.0095, "K": 39.0983}
"""dict: Molar masses (g/mol) of the compounds the inventory is given in."""

K40_ABUNDANCE = 1.17e-4
"""float: 40K per K (atom fraction, present day; NUBASE2020 via part 05 section 3.6, 2ND-HAND)."""

K40_DECAY_PER_GY = 0.5554
"""float: ln 2 / 1.248 Gyr (NUBASE2020)."""

K40_ARGON_BRANCH = 0.104
"""float: The share of 40K decays that gives 40Ar (electron capture; part 05 section 3.6)."""

# Design brackets: (low, high, shape); "log" is log-uniform, "uniform" is uniform. Each is a stand-in for the research item named.
INVENTORY = {
    # Ortenzi et al. 2020 Table 1, read in part 05 section 3.1: delivered H2O 0.05 and 0.5 wt%, CO2 0.01 and 0.1 wt% (CO2/H2O = 0.2 in
    # both), the dry class H2O 0 and CO2 0.6 wt%. The dry water is a bracket (their stored 0.005 wt% is the top of it).
    "co2_over_h2o_wet": (0.2, 0.2, "uniform"),
    "h2o_dry_wt_percent": (1e-4, 5e-3, "log"),
    "co2_dry_wt_percent": (0.6, 0.6, "uniform"),
    # The water-rich class (b4: at least 10 weight percent water): the carbon share is open.
    "h2o_rich_wt_percent": (10.0, 50.0, "log"),
    "co2_over_h2o_rich": (5e-3, 0.2, "log"),
    # The nitrogen and sulfur per carbon atom: the bulk silicate Earth (Kargel and Lewis 1993 via Schaefer, Lodders and Fegley
    # 2012 Table 2: N/C 0.013, S/C 1.7) and the CI chondrites (N/C 0.04, S/C 0.6) are inside; open, GEN.481.
    "nitrogen_per_carbon": (1e-3, 0.1, "log"),
    "sulfur_per_carbon": (0.1, 3.0, "log"),
    # Potassium of the rock (weight percent), for the 40Ar: the bulk silicate Earth 0.024, the continental crust 2.1: open, GEN.492.
    "potassium_wt_percent": (0.005, 0.1, "log"),
    # The share of the 40Ar produced that has reached the air (part 05 section 3.6: 0.05 to 0.7).
    "argon_degassed": (0.05, 0.7, "uniform"),
}
"""dict: name -> (low, high, shape), the draws of the inventory that no table of the recipes holds."""

MAGMATIC_TEMPERATURE_K = (1100.0, 1500.0, "uniform")
"""tuple: The temperature the volcano's gas equilibrates with the mantle at (part 05 section 3.4 computes at 1200 K; the
bracket is the design's: basalt erupts at 1450 to 1550 K, the gas re-equilibrates in the conduit and at the vent below that)."""

QUENCH_TEMPERATURE_K = (400.0, 900.0, "log")
"""tuple: The temperature a cool body's outgassed gas keeps (Zahnle and Marley 2014: 400 to 600 K for planets; Liggins 2022: 800 K)."""

OXIDATION_TIME_GY = (0.1, 1.5, "log")
"""tuple: The time a body with long-lived geology takes to move its mantle from the magma-ocean redox to the later upper mantle's:
the e-folding time of the oxidation (Aulbach and Stagno 2016 put the Earth's upper mantle at its present redox by about 3.5 Ga;
bracket, GEN.482)."""

MAX_EQUILIBRIUM_K = 5000.0
"""float: Hotter outgassing temperatures are held at this (the fits of the main species end at 6000 K)."""

UNLISTED_LIMIT = 1e-3
"""float: The most mole fraction that species the gas table does not hold may add up to before the draw says so (`MixtureError`)."""

ROCKY_CATEGORIES = "abc"
"""str: The recipe tables the rocky outgassing reads."""

ROCKY_REDOX_CLASSES = {"a5": "reduced", "a6": "intermediate", "a7": "oxidised"}


class MixtureError(RuntimeError):
    """A mixture that cannot be drawn: a quantity the tables read that the subject lacks, or a gas outside the gas table."""


@dataclass(frozen=True)
class MixtureDraw:
    """
    A drawn mixture.

    Attributes:
        mole_fractions (dict): Gas-table formula -> mole fraction, summing to 1 (trace gases below 1e-30 dropped).
        drawn (dict): Recipe group -> the row drawn (`{"a": "a6", "b": "b2", "c": "c3"}`).
        redox_diw (float or None): log10 fO2 above the iron-wuestite buffer of the outgassing mantle (Hirschmann 2021 scale).
        temperature_k (float): The temperature of the equilibrium.
        pressure_bar (float): The pressure of the equilibrium.
        inventory (dict): The element budget the chemistry read, atoms per 100 g of body (`H`, `C`, `N`, `S`, `Ar`).
        constants (dict): The design constants this draw drew, name -> value.
        unlisted (dict): Species the equilibrium holds that the gas table does not (with their mole fractions before the
            remainder was renormalised); empty unless a hot body makes atoms or radicals.
    """

    mole_fractions: dict
    drawn: dict
    redox_diw: float
    temperature_k: float
    pressure_bar: float
    inventory: dict
    constants: dict = field(default_factory=dict)
    unlisted: dict = field(default_factory=dict)


# --- drawing helpers -------------------------------------------------------------------------------------------------------

def _bracket(constants, name, bracket):
    """Draws one design constant of `bracket` = (low, high, shape) and records it under `name`."""
    low, high, shape = bracket
    if low == high:
        value = low
    elif shape == "log":
        value = 10.0 ** draw.uniform(math.log10(low), math.log10(high))
    else:
        value = draw.uniform(low, high)
    constants[name] = value
    return value


def _require(subject, names):
    missing = [name for name in names if subject.get(name) is None]
    if missing:
        raise MixtureError(f"the subject of class {subject.body_class!r} lacks {', '.join(missing)}")


def required_inputs(body_class, categories=ROCKY_CATEGORIES):
    """
    The quantities a subject of `body_class` must carry for `RECIPES.draw_all(subject, categories)` and the chemistry to run: every
    stored quantity a row of the tables reads (in a weight, an intermediate, a distribution, a condition), less the ones an
    earlier row of the same tables draws (`sets`) and the draws themselves, plus what this module reads.

    Returns:
        tuple: Sorted names from `atmosphere_recipes.INPUTS`.
    """
    names, drawn = set(), set()
    for letter in categories:
        for entry in RECIPES.in_category(letter):
            if atmosphere_recipes.pt.ANY_CLASS not in entry.classes and body_class not in entry.classes:
                continue
            expressions = list(entry.lets.values()) + list(RECIPES._dist_expressions(entry))
            if entry.weight is not None:
                expressions.append(entry.weight)
            for text in entry.spec.get("values", {}).values():
                expressions.append(atmosphere_recipes.Expression(text, f"row {entry.id}"))
            for expression in expressions:
                names |= expression.names
            names |= {condition.quantity for condition in entry.conditions}
            if entry.spec.get("sets"):
                drawn.add(entry.spec["sets"])
    wanted = (names - drawn) & set(atmosphere_recipes.INPUTS)
    return tuple(sorted(wanted | {"ts_k", "pressure_bar", "mass_earth", "age_gy"}))


# --- the redox ------------------------------------------------------------------------------------------------------------

def draw_redox(subject, recipe, constants):
    """
    The oxygen fugacity of the mantle that outgasses, log10 units above the iron-wuestite buffer (Hirschmann 2021 scale).

    The recipe draw holds the magma-ocean stage's value `a1` (Sossi 2020: about IW+0.5 at an Earth mass, about IW-2 for a small
    body), the later upper mantle's `a2`, the building blocks' `a3` and the white-dwarf debris `a4`; the body's own geology and age
    pick the one the air was outgassed over. A body with long-lived geology moves from `a1` to `a2` as
    `1 - exp(-age / tau)` with the oxidation time `tau` of `OXIDATION_TIME_GY`.

    Args:
        subject (Subject): With `age_gy` and, for a tectonic body, `long_lived_geology`.
        recipe (RecipeDraw): What `RECIPES.draw_all` drew for `subject`.
        constants (dict): The design constants drawn so far; this adds `oxidation_time_gy` when it is used.

    Returns:
        float: dIW, or `None` when no redox row applies (a giant).
    """
    values = recipe.values
    if "a4" in values:
        return values["a4"]
    if "a3" in values and "a1" not in values:
        return values["a3"]
    if "a1" not in values:
        return values.get("a2")
    first = values["a1"]
    if "a2" in values and subject.get("long_lived_geology"):
        tau = _bracket(constants, "oxidation_time_gy", OXIDATION_TIME_GY)
        progress = 1.0 - math.exp(-max(subject["age_gy"], 0.0) / tau)
        return first + (values["a2"] - first) * progress
    return first


# --- the rocky outgassing --------------------------------------------------------------------------------------------------

def _water_and_carbon_wt_percent(recipe, constants):
    """H2O and CO2 of the outgassing stage, weight percent of the body, from the delivered class and its continuous water fraction."""
    groups = recipe.groups
    row = groups["b"].name if "b" in groups else "b2"
    if row in ("b2", "b3"):
        water = 100.0 * recipe.values["b5"] if "b5" in recipe.values else 100.0 * 10.0 ** draw.uniform(math.log10(5e-4), math.log10(5e-3))
        return water, water * _bracket(constants, "co2_over_h2o_wet", INVENTORY["co2_over_h2o_wet"])
    if row == "b1":
        return (_bracket(constants, "h2o_dry_wt_percent", INVENTORY["h2o_dry_wt_percent"]),
                _bracket(constants, "co2_dry_wt_percent", INVENTORY["co2_dry_wt_percent"]))
    water = _bracket(constants, "h2o_rich_wt_percent", INVENTORY["h2o_rich_wt_percent"])
    return water, water * _bracket(constants, "co2_over_h2o_rich", INVENTORY["co2_over_h2o_rich"])


def rocky_inventory(subject, recipe, constants):
    """
    The element budget of the outgassing stage, atoms per 100 g of body: hydrogen and carbon from the delivered water and carbon
    dioxide, nitrogen and sulfur per carbon atom, and the argon the potassium of the rock has made by the body's age.
    """
    water, carbon_dioxide = _water_and_carbon_wt_percent(recipe, constants)
    hydrogen = 2.0 * water / AMU_G_MOL["H2O"]
    carbon = carbon_dioxide / AMU_G_MOL["CO2"]
    nitrogen = carbon * _bracket(constants, "nitrogen_per_carbon", INVENTORY["nitrogen_per_carbon"])
    sulfur = carbon * _bracket(constants, "sulfur_per_carbon", INVENTORY["sulfur_per_carbon"])
    potassium = _bracket(constants, "potassium_wt_percent", INVENTORY["potassium_wt_percent"]) / AMU_G_MOL["K"]
    degassed = _bracket(constants, "argon_degassed", INVENTORY["argon_degassed"])
    argon = degassed * K40_ARGON_BRANCH * K40_ABUNDANCE * potassium * math.expm1(K40_DECAY_PER_GY * max(subject["age_gy"], 0.0))
    return {"H": hydrogen, "C": carbon, "N": nitrogen, "S": sulfur, "Ar": argon}


def _atoms(equilibrium):
    """The moles of atoms of every element in the gas of `equilibrium`, oxygen included."""
    atoms = {}
    for name, moles in equilibrium.moles.items():
        for symbol, count in species_thermo.SPECIES[name].elements.items():
            atoms[symbol] = atoms.get(symbol, 0.0) + count * moles
    return atoms


def _to_gas_table(fractions):
    """The mole fractions of species names as the gas table's formulas; what the table does not hold is returned apart."""
    listed, unlisted = {}, {}
    for name, value in fractions.items():
        formula = species_thermo.SPECIES[name].formula
        if gases.has(formula):
            listed[formula] = listed.get(formula, 0.0) + value
        else:
            unlisted[formula] = unlisted.get(formula, 0.0) + value
    return listed, unlisted


def outgassed_mixture(subject, recipe=None):
    """
    The gas a rocky body's mantle outgasses, drawn for `subject` (see the module text).

    Args:
        subject (Subject): The body: its class and the quantities the tables of `ROCKY_CATEGORIES` read, with `ts_k` (surface
            temperature, K), `pressure_bar` (the pressure of the air's base), `mass_earth` and `age_gy`.
        recipe (RecipeDraw, optional): The recipe draw to use; drawn here (`RECIPES.draw_all`) when not given.

    Returns:
        MixtureDraw: The mixture.

    Raises:
        MixtureError: A quantity is missing, or the hot gas holds more than `UNLISTED_LIMIT` of species the gas table lacks.
    """
    _require(subject, ("ts_k", "pressure_bar", "mass_earth", "age_gy"))
    if recipe is None:
        try:
            recipe = RECIPES.draw_all(subject, categories=ROCKY_CATEGORIES)
        except KeyError as error:
            raise MixtureError(str(error)) from error
    constants = {}
    redox = draw_redox(subject, recipe, constants)
    if redox is None:
        raise MixtureError(f"no redox row applies to the class {subject.body_class!r}: it has no mantle to outgas")
    inventory = rocky_inventory(subject, recipe, constants)
    magmatic = _bracket(constants, "magmatic_temperature_k", MAGMATIC_TEMPERATURE_K)
    quench = _bracket(constants, "quench_temperature_k", QUENCH_TEMPERATURE_K)
    temperature = min(max(subject["ts_k"], quench), MAX_EQUILIBRIUM_K)
    pressure = subject["pressure_bar"]
    volcanic = gas_equilibrium.solve(inventory, magmatic, pressure, fixed_o2_log_fugacity=fo2.log_fo2(redox, magmatic))
    if temperature == magmatic:
        result = volcanic
    else:
        result = gas_equilibrium.solve(_atoms(volcanic), temperature, pressure)
    listed, unlisted = _to_gas_table({name: value for name, value in result.mole_fractions.items() if value > 1e-30})
    if sum(unlisted.values()) > UNLISTED_LIMIT:
        raise MixtureError(f"{sum(unlisted.values()):.3g} of the gas is species the gas table lacks: {sorted(unlisted)}")
    total = sum(listed.values())
    return MixtureDraw(
        mole_fractions={formula: value / total for formula, value in sorted(listed.items(), key=lambda item: -item[1])},
        drawn={group: drawn.name for group, drawn in recipe.groups.items()},
        redox_diw=redox,
        temperature_k=temperature,
        pressure_bar=subject["pressure_bar"],
        inventory=inventory,
        constants=constants,
        unlisted=unlisted,
    )
