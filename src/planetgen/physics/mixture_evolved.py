# planetgen/physics/mixture_evolved.py

"""
Evolved secondary air and frost-equilibrium air (GEN.479; part 05 section 6 of the atmosphere research)
=======================================================================================================

The `j` rows of the recipes (Venus-like, Mars-like, Titan-like, Pluto-like and the abiotic-oxygen world) are class recipes whose
Solar System bodies are TEST rows of the research, never presets. This module draws each as chemistry on the same pieces as the
rocky family and gives the pieces their own seeded design brackets, each a stand-in for the research item that will replace it.

* **Oxidised terrestrial (`j1` Venus-like, `j2` Mars-like): the outgassed gas at high redox, minus its water and its sulfur.**
  Part 05 section 6.1: take the outgassed mixture, remove the water by escape or condensation, keep the nitrogen inventory and the
  radiogenic argon. The first stage is `mixture_rocky.volcanic_gas`; the water is removed down to the share `WATER_RESIDUAL` of the
  hydrogen, and the sulfur the surface sulfates and the cloud keep is
  taken out down to `SULFUR_SINK`; the rest cools to the surface or quench temperature (`mixture_rocky.cool`) with the oxygen
  fugacity held at the body's redox, because the crust's sulfates and ferric iron buffer an old dry air (it is CO2 with 3 percent N2 and
  a few ppm of CO and water, as Venus is). A Mars-like body's carbon dioxide ice and carbonate sink is the freeze-out of the
  pressure bookkeeping (GEN.452), applied to this gas.
* **Abiotic oxygen (`j5`).** The outgassed gas, its carbon, nitrogen and sulfur rescaled to a carbon of `CARBON_PER_OCEAN_WATER` per molecule of the
  surface water of one to ten oceans (`OCEAN_WATER_WT_PERCENT`) added to it; the
  hydrogen of the water escapes (the fraction `WATER_PHOTOLYSED`) and the oxygen it leaves behind stays in the closed gas, except for
  the share `OXYGEN_SINK` that the crust takes up ([R2769]: the O2 builds to tens or hundreds of bar around a late M
  dwarf). Nothing holds the oxygen fugacity here: the gas is free to have O2.
* **Titan-like (`j3`): nitrogen and methane at the methane's saturation.** The nitrogen is the balance; the methane is the
  surface's humidity `TITAN_HUMIDITY` of the saturation pressure at the surface temperature over the pressure of the air (the
  lakes and the damp ground are the reservoir), hydrogen a trace of the methane's photolysis (`TITAN_HYDROGEN`).
* **Frost equilibrium (`j4`, Pluto and Triton).** Each gas stands over the surface ice at the pressure of an ideal solution,
  `p_i = chi_i p_sat,i(T)` (Raoult), with nitrogen ice the balance of the ice; the nitrogen sets the total. The ice mixing
  ratios are open (`FROST_ICE`): Pluto's measured 3e-3 of methane in the nitrogen ice gives 1e-7 of methane in the gas at 40 K,
  where the air has 3e-3, so the methane that the air holds sits in methane-rich ice; the bracket spans both.
"""

import dataclasses

from planetgen.physics import saturation
from planetgen.physics import mixture_oxygen, mixture_rocky
from planetgen.physics.mixture_base import MixtureDraw, NotApplicable, atoms_of, bracket, finish, require

WATER_RESIDUAL = (1e-5, 1e-3, "log")
"""tuple: The share of the hydrogen a water-lost body keeps (Venus keeps 20 ppm of water in 96.5 percent carbon dioxide, TEST)."""

SULFUR_SINK = (1e-4, 1e-2, "log")
"""tuple: The share of the outgassed sulfur that stays in the air of a Venus-like or Mars-like body after the surface sulfates and the
cloud take the rest (Venus has 150 ppm of SO2 in CO2, a share of about 1e-4 of the BSE sulfur to carbon of 1.7; TEST); design."""

WATER_PHOTOLYSED = (0.3, 1.0, "uniform")
"""tuple: The share of the water an abiotic-oxygen world loses its hydrogen from (the tail of Luger and Barnes: up to ten oceans); design."""

OCEAN_WATER_WT_PERCENT = (0.0233, 0.233, "log")
"""tuple: The water an abiotic-oxygen world's surface holds before its hydrogen is lost, weight percent of the body: one to ten
Earth oceans (one ocean is 0.0233 percent of the Earth's mass; the ten are the top of the range [R2769] consider)."""

CARBON_PER_OCEAN_WATER = (0.01, 0.2, "log")
"""tuple: The carbon of an abiotic-oxygen world's gas per water molecule of its oceans, mole per mole: the Earth's atmosphere and carbonate
crust hold about 5e19 kg of carbon against 1.4e21 kg of ocean, 0.054 per molecule; a design bracket around it that GEN.481 replaces."""

OXYGEN_SINK = (0.0, 0.9, "uniform")
"""tuple: The share of the oxygen the lost water leaves that the crust oxidises away instead of the air keeping it; design."""

TITAN_HUMIDITY = (0.3, 1.0, "uniform")
"""tuple: The methane's relative humidity at the surface (Titan: 5.65 percent near the surface at 94 K against a saturation of 12 percent;
the measured humidity is about 0.5; [R3289])."""

TITAN_HYDROGEN = (3e-4, 3e-3, "log")
"""tuple: The hydrogen mole fraction of a methane atmosphere (Titan 1.0e-3 +/- 0.16e-3, [R3289], TEST); design bracket around it."""

FROST_ICE = {"CH4": (3e-3, 0.3, "log"), "CO": (1e-3, 0.1, "log")}
"""dict: The mole fraction of each minor gas in the surface ice (the rest is nitrogen), stand-ins for the research item on the ice
composition of cold bodies (part 05 section 6.2: the vapour-pressure relations of the ices are the Fray and Schmitt ones)."""

TITAN_MAX_METHANE = 0.5
"""float: The methane share the nitrogen-methane gas is held at when the saturation exceeds the pressure of the air."""

FAMILY_OF_ROW = {"j1": "venus", "j2": "mars", "j3": "titan", "j4": "frost", "j5": "abiotic-oxygen"}
"""dict: The `j` row -> the family the result records."""

FAMILIES = tuple(FAMILY_OF_ROW.values())


def _saturation_bar(formula, temperature_k):
    curve = saturation.saturation(formula, temperature_k)
    if curve.status == saturation.ABOVE_CRITICAL:
        raise NotApplicable(f"{formula} has no saturation pressure at {temperature_k} K (above its critical point)")
    return curve.pressure_pa / saturation.BAR_TO_PA


def _from_fractions(fractions, *, family, recipe, temperature_k, pressure_bar, inventory, constants):
    total = sum(fractions.values())
    return MixtureDraw(
        mole_fractions=dict(sorted(((formula, value / total) for formula, value in fractions.items() if value > 0.0),
                                   key=lambda item: -item[1])),
        family=family, drawn={group: drawn.name for group, drawn in recipe.groups.items()}, redox_diw=None,
        temperature_k=temperature_k, pressure_bar=pressure_bar, inventory=inventory, constants=constants)


def evolved_oxidised(subject, recipe, row, supply=None):
    """
    The gas of a Venus-like (`j1`), Mars-like (`j2`) or abiotic-oxygen (`j5`) body.

    Args:
        subject (Subject): As `mixture_rocky.outgassed_mixture` reads.
        recipe (RecipeDraw): With the evolved row drawn in group `j`.
        row (str): `j1`, `j2` or `j5`.
        supply (callable, optional): The share of each element of the delivered inventory that the surface holds
            (`mixture_rocky.volcanic_gas`).
    """
    constants = {}
    volcanic, redox, inventory = mixture_rocky.volcanic_gas(subject, recipe, constants, supply)
    atoms = atoms_of(volcanic)
    held = None
    if row in ("j1", "j2"):
        atoms["H"] = atoms.get("H", 0.0) * bracket(constants, "water_residual", WATER_RESIDUAL)
        atoms["S"] = atoms.get("S", 0.0) * bracket(constants, "sulfur_sink", SULFUR_SINK)
        atoms.pop("O", None)
        held = redox       # the crust buffers the oxygen: its sulfates and ferric iron take what the lost water's oxygen leaves
    else:
        ocean = 2.0 * bracket(constants, "ocean_water_wt_percent", OCEAN_WATER_WT_PERCENT) / mixture_rocky.AMU_G_MOL["H2O"]
        carbon = ocean / 2.0 * bracket(constants, "carbon_per_ocean_water", CARBON_PER_OCEAN_WATER)
        for symbol in ("N", "S"):
            atoms[symbol] = carbon * inventory[symbol] / inventory["C"]        # nitrogen and sulfur stayed per carbon atom
        atoms["C"] = carbon
        atoms["H"] = atoms.get("H", 0.0) + ocean
        atoms["O"] += ocean / 2.0
        lost = atoms["H"] * bracket(constants, "water_photolysed", WATER_PHOTOLYSED)
        atoms["H"] -= lost
        atoms["O"] -= (lost / 2.0) * bracket(constants, "oxygen_sink", OXYGEN_SINK)
    atoms = {symbol: value for symbol, value in atoms.items() if value > 0.0}
    result = mixture_rocky.cool(atoms, subject, constants, held_diw=held, quenched=False)
    draw_result = finish(result, family=FAMILY_OF_ROW[row], drawn={group: drawn.name for group, drawn in recipe.groups.items()},
                         redox_diw=redox, inventory=inventory, constants=constants, supplied=True)
    if row == "j5":
        capped = mixture_oxygen.cap_oxygen(draw_result.mole_fractions, mixture_oxygen.ABIOTIC_O2_MAX)
        # The oxygen over the cap is the crust's: the other gases keep their moles, so the air is what they make at the new shares.
        others = draw_result.moles * (1.0 - draw_result.mole_fractions.get("O2", 0.0))
        draw_result = dataclasses.replace(draw_result, mole_fractions=capped, moles=others / (1.0 - capped.get("O2", 0.0)))
    return draw_result


def titan_like(subject, recipe, supply=None):
    """
    The nitrogen-methane gas of a cold body with a thick air (`j3`), from the surface temperature and the pressure of the air; its
    moles are the nitrogen the surface holds (`supply` of the delivered inventory, as for the rocky family) over the nitrogen's share.
    """
    require(subject, ("ts_k", "pressure_bar", "age_gy"))
    constants = {}
    inventory = mixture_rocky.rocky_inventory(subject, recipe, constants)
    nitrogen = mixture_rocky.surface_budget(inventory, supply).get("N", 0.0)
    if nitrogen <= 0.0:
        raise NotApplicable("no nitrogen has reached the surface")
    temperature, pressure = subject["ts_k"], subject["pressure_bar"]
    humidity = bracket(constants, "methane_humidity", TITAN_HUMIDITY)
    methane = min(humidity * _saturation_bar("CH4", temperature) / pressure, TITAN_MAX_METHANE)
    hydrogen = bracket(constants, "hydrogen_fraction", TITAN_HYDROGEN)
    constants["methane_saturation_bar"] = _saturation_bar("CH4", temperature)
    result = _from_fractions({"N2": 1.0 - methane - hydrogen, "CH4": methane, "H2": hydrogen}, family="titan", recipe=recipe,
                             temperature_k=temperature, pressure_bar=pressure, inventory=inventory, constants=constants)
    return dataclasses.replace(result, moles=0.5 * nitrogen / result.mole_fractions["N2"])


def frost_equilibrium(subject, recipe):
    """
    The gas over a surface of nitrogen, methane and carbon monoxide ice (`j4`): the partial pressure of each is its ice mole
    fraction times its saturation pressure at the surface temperature, and the nitrogen sets the pressure the air has.
    """
    require(subject, ("ts_k",))
    constants = {}
    temperature = subject["ts_k"]
    ice = {formula: bracket(constants, f"ice_{formula.lower()}", spec) for formula, spec in FROST_ICE.items()}
    ice["N2"] = 1.0 - sum(ice.values())
    if ice["N2"] <= 0.0:
        raise NotApplicable("the minor ices exceed the nitrogen ice")
    partial = {formula: share * _saturation_bar(formula, temperature) for formula, share in ice.items()}
    total = sum(partial.values())
    return _from_fractions(partial, family="frost", recipe=recipe, temperature_k=temperature, pressure_bar=total,
                           inventory=ice, constants=constants)


def draw(subject, recipe, family, supply=None):
    """The mixture of a body whose recipe drew an evolved or cold-body row (`family` is `FAMILY_OF_ROW`'s value for it)."""
    if family == "titan":
        return titan_like(subject, recipe, supply)
    if family == "frost":
        return frost_equilibrium(subject, recipe)
    row = recipe.groups["j"].name
    return evolved_oxidised(subject, recipe, row, supply)
