# planetgen/physics/mixture_trace.py

"""
The trace and minor gases of an air, drawn after the life stage (GEN.479; table m of the possibility tables)
=============================================================================================================

The major gases of `mixture_draw.draw_mixture` come from chemistry before anything lives on the body. The trace gases of table `m`
read the life stage (`o2_pal`, `life_chemical`, `lightning`, whether the air is CO2-rich), so they are drawn after it, from the same
seeded stream, on top of the mixture: `draw_trace_gases(subject, draw)` returns the mixture with the drawn gases added at their
mole fractions and everything renormalised. What the rows decide is presence; the amounts are design brackets, each a stand-in with the
measured value as its TEST anchor:

| row | gas | mole fraction bracket | TEST anchor |
|---|---|---|---|
| m3 | N2O | 1e-8 to 1e-6 | Earth, about 3.3e-7 ([R0052]: ~400 ppbv) |
| m4 | PH3 | 1e-10 to 2e-8 | the 20 ppb Venus claim [R1672], a lifetime of hours [R4163] |
| m5 | HCl | 1e-6 to 1e-2 | 1.4 percent at 1000 K and 100 bar over crust composition [R3945] |
| m6 | HF | 1e-7 to 1e-3 | 4.5 percent at 1000 K and 1e-6 bar [R3945] |
| m11 | NO2, CH3Cl | 1e-11 to 1e-8, 1e-11 to 1e-9 | none on any exoplanet |

The rows that are not abundances of this module: the ozone rows m1 and m2 are the ozone column's (GEN.413: the draw is in the result,
`TraceDraw.recipe.binaries`), the alkali rows m7 and m8 decide the sodium and potassium gas of a giant and are read by `mixture_giant`, m9 (OCS)
is a formation equilibrium the solver already gives, and m10 (radon) is the ground source of the radiation chain.
"""

import dataclasses

from planetgen.generation.possibility_table import Subject
from planetgen.physics import atmosphere_recipes
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.physics.mixture_base import MixtureError, bracket

ABUNDANCES = {
    "m3": {"N2O": (1e-8, 1e-6, "log")},
    "m4": {"PH3": (1e-10, 2e-8, "log")},
    "m5": {"HCl": (1e-6, 1e-2, "log")},
    "m6": {"HF": (1e-7, 1e-3, "log")},
    "m11": {"NO2": (1e-11, 1e-8, "log"), "CH3Cl": (1e-11, 1e-9, "log")},
}
"""dict: Row -> gas -> (low, high, shape) of the mole fraction the gas has when the row draws it present."""


def required_trace_inputs(body_class):
    """The stored quantities the rows of table m read for `body_class` (the life stage and the air's character among them)."""
    names = set()
    for entry in RECIPES.in_category("m"):
        if atmosphere_recipes.pt.ANY_CLASS not in entry.classes and body_class not in entry.classes:
            continue
        expressions = list(entry.lets.values()) + list(RECIPES._dist_expressions(entry))
        if entry.weight is not None:
            expressions.append(entry.weight)
        for expression in expressions:
            names |= expression.names
        names |= {condition.quantity for condition in entry.conditions if not condition.quantity.startswith("drawn.")}
    return tuple(sorted(names & set(atmosphere_recipes.INPUTS)))


def draw_trace_gases(subject, mixture):
    """
    The mixture with the trace gases of table m added.

    Args:
        subject (Subject): The body, with the quantities `required_trace_inputs(subject.body_class)` names.
        mixture (MixtureDraw): From `draw_mixture`; its `recipe` is the draw the rows of table m follow.

    Returns:
        MixtureDraw: The same mixture with the present trace gases added (mole fractions renormalised to 1), `trace` in `constants`
            holding the rows drawn and the amounts, and `recipe` the table-m draw's `RecipeDraw` merged into the first (so
            `recipe.binaries` has the ozone rows for the ozone code).

    Raises:
        MixtureError: A quantity the rows read is missing (named) or the mixture has no recipe.
    """
    if mixture.recipe is None:
        raise MixtureError("the mixture carries no recipe draw")
    quantities = {**mixture.recipe.subject.quantities, **subject.quantities}
    try:
        trace_recipe = RECIPES.draw_all(Subject(subject.body_class, quantities), categories="m")
    except KeyError as error:
        raise MixtureError(str(error)) from error
    constants = dict(mixture.constants)
    added = {}
    for row, gases_of_row in ABUNDANCES.items():
        if trace_recipe.binaries.get(row, (False, None))[0]:
            for formula, spec in gases_of_row.items():
                added[formula] = bracket(constants, f"trace_{formula}", spec)
    fractions = dict(mixture.mole_fractions)
    if fractions:
        for formula, value in added.items():
            fractions[formula] = fractions.get(formula, 0.0) + value
        total = sum(fractions.values())
        fractions = dict(sorted(((formula, value / total) for formula, value in fractions.items()), key=lambda item: -item[1]))
    constants["trace_rows"] = sorted(row for row, (present, _) in trace_recipe.binaries.items() if present)
    merged = dataclasses.replace(mixture.recipe, binaries={**mixture.recipe.binaries, **trace_recipe.binaries})
    return dataclasses.replace(mixture, mole_fractions=fractions, constants=constants, recipe=merged)
