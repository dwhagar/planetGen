# planetgen/physics/mixture_coma.py

"""
The gas of a comet or an icy asteroid's coma (GEN.479; part 05 section 9 of the atmosphere research)
=====================================================================================================

A coma is not an equilibrium: the ices sublimate into a free expansion, so the gas is the ice's composition. The recipe tables
(`k` rows) draw it as abundances relative to water: CO2/H2O and CO/H2O log-normal (k4, k5), methanol, ammonia and hydrogen cyanide
log-uniform over the cross-comet ranges (k6), and the species only one comet (67P) measured log-normal about that comet's value
(k7). The class row (k1 CO2-dominated, k2 CO-dominated, k3 comparable) orders the two: a CO2-dominated coma has at least as much CO2 as
CO, a CO-dominated one at least as much CO as CO2, and a comparable one has both at their geometric mean. The mole fractions are the
ratios over their sum (water 1). The pressure is the outflow's and not a balance: `pressure_bar` is 0 and the generator's coma
code (not the mixture) says how thick the coma is.
"""

import math

from planetgen.physics import gases
from planetgen.physics.mixture_base import MixtureDraw, MixtureError

COMA_CATEGORIES = "k"

FORMULA_OF = {"H2CO": "CH2O"}
"""dict: The recipe tables' name of a species -> the gas table's formula."""

DOMINANCE = {"k1": "CO2", "k2": "CO", "k3": None}
"""dict: The class row -> the species that is the larger of CO2 and CO (`None`: comparable)."""


def coma_mixture(subject, recipe):
    """
    The coma gas of `subject`: its class row, the ratios the recipe drew, normalised over water.

    Args:
        subject (Subject): The comet or icy asteroid (only for its temperature, if it carries `ts_k`).
        recipe (RecipeDraw): With the `k` group and the values `k4` to `k7` drawn.

    Raises:
        MixtureError: The recipe drew no coma ratios, or a species is not in the gas table.
    """
    try:
        group = recipe.groups["k"].name
        ratios = {"H2O": 1.0, "CO2": recipe.values["k4"], "CO": recipe.values["k5"]}
        ratios.update(recipe.values["k6"])
        ratios.update({FORMULA_OF.get(name, name): value for name, value in recipe.values["k7"].items()})
    except KeyError as error:
        raise MixtureError(f"the recipe drew no coma ({error})") from error
    dominant = DOMINANCE[group]
    if dominant is None:
        ratios["CO2"] = ratios["CO"] = math.sqrt(ratios["CO2"] * ratios["CO"])
    else:
        other = "CO" if dominant == "CO2" else "CO2"
        if ratios[dominant] < ratios[other]:
            ratios["CO2"], ratios["CO"] = ratios["CO"], ratios["CO2"]
    for formula in ratios:
        if not gases.has(formula):
            raise MixtureError(f"the coma species {formula} is not in the gas table")
    total = sum(ratios.values())
    return MixtureDraw(
        mole_fractions=dict(sorted(((formula, value / total) for formula, value in ratios.items()), key=lambda item: -item[1])),
        family="coma", drawn={name: drawn.name for name, drawn in recipe.groups.items()}, redox_diw=None,
        temperature_k=float(subject.get("ts_k") or 0.0), pressure_bar=0.0, inventory=dict(ratios),
        constants={"co2_to_h2o": ratios["CO2"], "co_to_h2o": ratios["CO"]})
