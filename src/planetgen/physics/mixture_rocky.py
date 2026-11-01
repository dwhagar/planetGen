# planetgen/physics/mixture_rocky.py

"""
The gas a rocky body outgasses (GEN.451, GEN.479; docs/design/composition/mixture-draw.md)
==============================================================================

The rocky family of `mixture_draw.draw_mixture`. The outgassed mixture is the gas in equilibrium with a mantle of a given oxygen
fugacity (Schaefer and Fegley 2010; Gaillard and Scaillet 2014; Ortenzi et al. 2020, part 05 section 3 of the atmosphere
research), as the research recipe 3.5 states, with the chemistry of `physics/gas_equilibrium.py` over the NASA Glenn
thermodynamics of `physics/species_thermo.py` and no per-class template:

    draw = outgassed_mixture(Subject("A", {... the quantities required_inputs("A") names ...}))
    draw.mole_fractions       # {"H2O": 0.78, "CO2": 0.16, ...}, gas-table formulas, sums to 1

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

from planetgen.physics import dry_carbon, gas_equilibrium
from planetgen.physics import oxygen_fugacity as fo2
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.physics.mixture_base import MixtureError, NotApplicable, atoms_of, bracket, finish, require
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
    # both), the dry class H2O 0 and CO2 0.6 wt% delivered, 0.018 wt% stored. The dry water is a bracket (their stored 0.005 wt% is the
    # top of it). The dry carbon is the bracket from the stored 0.018 to the wet High case's 0.1 wt%, not the delivered 0.6: an Earth twin
    # draws the dry class nine times in ten, and 0.6 wt% CO2 is sixteen times the carbon of the bulk silicate Earth (100 ppm C, Bergin
    # et al. 2015 Table S1), where this bracket's log-middle (0.042 wt% CO2, 115 ppm C) is within a sixth of it (GEN.813).
    # GEN.1005: the dry carbon is now `dry_carbon.co2_wt_percent_bracket()`, 0.029 to 0.108 wt% CO2 from the carbon of non-carbonaceous building blocks,
    # the share already differentiated and the core's share (Blanchard et al. 2022); the bulk silicate Earth is a test of it.
    "co2_over_h2o_wet": (0.2, 0.2, "uniform"),
    "h2o_dry_wt_percent": (1e-4, 5e-3, "log"),
    "co2_dry_wt_percent": (*dry_carbon.co2_wt_percent_bracket(), "log"),
    # The water-rich class (b4: at least 10 weight percent water): the carbon share is open.
    # CO2 over H2O: Kite and Ford 2018 give c_C/f_W 1e-6 to about 0.1 (Sec 3.1) and 0.3 to 30 percent of the water (Sec 4.3), which disagree
    # inside the paper; the widest is kept, so the top of 0.2 is inside the second (GEN.824, AM14).
    "h2o_rich_wt_percent": (10.0, 50.0, "log"),
    "co2_over_h2o_rich": (5e-3, 0.2, "log"),
    # The nitrogen and sulfur per carbon atom: the bulk silicate Earth (Kargel and Lewis 1993 via Schaefer, Lodders and Fegley
    # 2012 Table 2: N/C 0.013, S/C 1.7) and the CI chondrites (N/C 0.04, S/C 0.6) are inside; open, GEN.481. The nitrogen bracket
    # runs from 0.005 to 0.1: Bergin et al. 2015 Table S1 gives N/C 0.020 for the bulk silicate Earth (C/N 49 plus or minus 9), 0.059
    # for the CI chondrites, 0.074 for the air of Venus and 0.055 for the air of Mars, and its log-middle 0.022 is the Earth's (GEN.813).
    "nitrogen_per_carbon": (5e-3, 0.1, "log"),
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
the e-folding time of the oxidation (studies put the Earth's upper mantle at its present redox by about 3.5 Ga;
bracket, GEN.482)."""

MIN_EQUILIBRIUM_K = 300.0
"""float: The least temperature that every species of `species_thermo` is fitted for (the sulfur and phosphorus gases start at 300 K)."""

MAX_EQUILIBRIUM_K = 5000.0
"""float: Hotter outgassing temperatures are held at this (the fits of the main species end at 6000 K)."""

ROCKY_CATEGORIES = "abc"
"""str: The recipe tables the rocky outgassing reads."""

# --- the redox ------------------------------------------------------------------------------------------------------------

def draw_redox(subject, recipe, constants):
    """
    The oxygen fugacity of the mantle that outgasses, log10 units above the iron-wuestite buffer (Hirschmann 2021 scale).

    The recipe draw holds the magma-ocean stage's value `a1` (Sossi 2020: about IW+0.5 at an Earth mass, about IW-2 for a small
    body), the later upper mantle's `a2`, the building blocks' `a3` and the white-dwarf debris `a4`; the body's own geology and age
    pick the one the air was outgassed over. A body with long-lived geology moves from `a1` to `a2` as
    `1 - exp(-age / tau)` with the oxidation time `tau` of `OXIDATION_TIME_GY`; the later mantle is never more reduced than the
    magma-ocean surface of the same body (Deng 2020, Hirschmann 2022: the later processes only oxidise; GEN.482b).

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
        tau = bracket(constants, "oxidation_time_gy", OXIDATION_TIME_GY)
        progress = 1.0 - math.exp(-max(subject["age_gy"], 0.0) / tau)
        return first + (max(values["a2"], first) - first) * progress
    return first


# --- the rocky outgassing --------------------------------------------------------------------------------------------------

def _water_and_carbon_wt_percent(recipe, constants):
    """H2O and CO2 of the outgassing stage, weight percent of the body, from the delivered class and its continuous water fraction."""
    groups = recipe.groups
    row = groups["b"].name if "b" in groups else "b2"
    if row not in ("b1", "b2", "b3", "b4"):
        raise MixtureError(f"the delivered-volatile table drew {row!r}, which has no inventory")
    if row in ("b2", "b3"):
        water = 100.0 * recipe.values["b5"] if "b5" in recipe.values else 100.0 * 10.0 ** draw.uniform(math.log10(5e-4), math.log10(5e-3))
        return water, water * bracket(constants, "co2_over_h2o_wet", INVENTORY["co2_over_h2o_wet"])
    if row == "b1":
        return (bracket(constants, "h2o_dry_wt_percent", INVENTORY["h2o_dry_wt_percent"]),
                bracket(constants, "co2_dry_wt_percent", INVENTORY["co2_dry_wt_percent"]))
    water = bracket(constants, "h2o_rich_wt_percent", INVENTORY["h2o_rich_wt_percent"])
    return water, water * bracket(constants, "co2_over_h2o_rich", INVENTORY["co2_over_h2o_rich"])


def rocky_inventory(subject, recipe, constants):
    """
    The element budget of the outgassing stage, atoms per 100 g of body: hydrogen and carbon from the delivered water and carbon
    dioxide, nitrogen and sulfur per carbon atom, and the argon the potassium of the rock has made by the body's age.
    """
    water, carbon_dioxide = _water_and_carbon_wt_percent(recipe, constants)
    hydrogen = 2.0 * water / AMU_G_MOL["H2O"]
    carbon = carbon_dioxide / AMU_G_MOL["CO2"]
    nitrogen = carbon * bracket(constants, "nitrogen_per_carbon", INVENTORY["nitrogen_per_carbon"])
    sulfur = carbon * bracket(constants, "sulfur_per_carbon", INVENTORY["sulfur_per_carbon"])
    potassium = bracket(constants, "potassium_wt_percent", INVENTORY["potassium_wt_percent"]) / AMU_G_MOL["K"]
    degassed = bracket(constants, "argon_degassed", INVENTORY["argon_degassed"])
    argon = degassed * K40_ARGON_BRANCH * K40_ABUNDANCE * potassium * math.expm1(K40_DECAY_PER_GY * max(subject["age_gy"], 0.0))
    return {"H": hydrogen, "C": carbon, "N": nitrogen, "S": sulfur, "Ar": argon}


def volcanic_gas(subject, recipe, constants, supply=None):
    """
    The first stage: the gas in equilibrium with the mantle at the magmatic temperature and the pressure of the layer, the
    oxygen fugacity held. Draws the redox, the inventory and the magmatic temperature (recorded in `constants`).

    Args:
        supply (callable, optional): `(inventory) -> {element: share}`, the share of the delivered inventory that the surface holds
            (`physics/surface_pressure.py`: the share the interior has degassed, less what the sinks took); an element it lacks is
            taken whole, and a body with no gas at the surface raises `NotApplicable`.

    Returns:
        tuple: `(Equilibrium, redox_diw, inventory)`; the inventory is the delivered one, before `supply`.

    Raises:
        MixtureError: A quantity is missing, or the class has no mantle (no redox row applies).
        NotApplicable: Nothing has reached the surface.
    """
    require(subject, ("ts_k", "pressure_bar", "mass_earth", "age_gy"))
    redox = draw_redox(subject, recipe, constants)
    if redox is None:
        raise MixtureError(f"no redox row applies to the class {subject.body_class!r}: it has no mantle to outgas")
    inventory = rocky_inventory(subject, recipe, constants)
    magmatic = bracket(constants, "magmatic_temperature_k", MAGMATIC_TEMPERATURE_K)
    budget = surface_budget(inventory, supply)
    volcanic = gas_equilibrium.solve(budget, magmatic, subject["pressure_bar"], fixed_o2_log_fugacity=fo2.log_fo2(redox, magmatic))
    return volcanic, redox, inventory


def surface_budget(inventory, supply):
    """The delivered `inventory` times the share of each element at the surface (`supply(inventory)`); elements with none are dropped."""
    shares = {} if supply is None else supply(inventory)
    budget = {element: atoms * shares.get(element, 1.0) for element, atoms in inventory.items()}
    budget = {element: atoms for element, atoms in budget.items() if atoms > 0.0}
    if not budget:
        raise NotApplicable("no gas has reached the surface")
    return budget


def cool(atoms, subject, constants, volcanic=None, held_diw=None, quenched=True):
    """
    The second stage: the gas, a closed system of `atoms` (moles of atoms, oxygen included), cooled to the surface temperature or,
    when the surface is colder, to the quench temperature (drawn here and recorded in `constants`).

    Args:
        volcanic (Equilibrium, optional): The first stage; where the result would be at the magmatic temperature already
            it is returned as it is.
        held_diw (float, optional): Hold the oxygen fugacity at this log10 units above iron-wuestite (a surface that buffers the
            air, as an evolved oxidised body's crust does); `atoms` then carry no oxygen.
        quenched (bool): Freeze the gas at the quench temperature when the surface is colder (a fresh outgassed air whose reactions
            stop, drawn here and recorded in `constants`); `False` for an old air that has had the time to settle to the surface
            temperature, with the least temperature the species of the table are all fitted for (`MIN_EQUILIBRIUM_K`).
    """
    if quenched:
        quench = bracket(constants, "quench_temperature_k", QUENCH_TEMPERATURE_K)
    else:
        quench = MIN_EQUILIBRIUM_K
    temperature = min(max(subject["ts_k"], quench), MAX_EQUILIBRIUM_K)
    if volcanic is not None and held_diw is None and temperature == volcanic.temperature_k:
        return volcanic
    held = None if held_diw is None else fo2.log_fo2(held_diw, temperature)
    return gas_equilibrium.solve(atoms, temperature, subject["pressure_bar"], fixed_o2_log_fugacity=held)


def outgassed_mixture(subject, recipe=None, supply=None):
    """
    The gas a rocky body's mantle outgasses, drawn for `subject` (see the module text).

    Args:
        subject (Subject): The body: its class and the quantities the tables of `ROCKY_CATEGORIES` read, with `ts_k` (surface
            temperature, K), `pressure_bar` (the pressure of the air's base), `mass_earth` and `age_gy`.
        recipe (RecipeDraw, optional): The recipe draw to use; drawn here (`RECIPES.draw_all`) when not given.
        supply (callable, optional): The share of each element of the delivered inventory that the surface holds (`volcanic_gas`).

    Returns:
        MixtureDraw: The mixture, with the `moles` of gas per 100 g of mantle.

    Raises:
        MixtureError: A quantity is missing, or the hot gas holds more than `UNLISTED_LIMIT` of species the gas table lacks.
    """
    if recipe is None:
        try:
            recipe = RECIPES.draw_all(subject, categories=ROCKY_CATEGORIES)
        except KeyError as error:
            raise MixtureError(str(error)) from error
    constants = {}
    volcanic, redox, inventory = volcanic_gas(subject, recipe, constants, supply)
    result = cool(atoms_of(volcanic), subject, constants, volcanic)
    return finish(result, family="rocky", drawn={group: drawn.name for group, drawn in recipe.groups.items()}, redox_diw=redox,
                  inventory=inventory, constants=constants, supplied=True)
