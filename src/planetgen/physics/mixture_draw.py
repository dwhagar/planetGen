# planetgen/physics/mixture_draw.py

"""
The atmosphere mixture of a body, drawn from its regimes (GEN.479; docs/design/composition/mixture-draw.md)
===============================================================================================

`physics/atmosphere_recipes.py` (GEN.478) draws which regime a body's air is in: the mantle redox, the volatile delivery, the
outgassed class, the retention class, the envelope, the giant, brown-dwarf and hot-Jupiter regimes, the lava-world vapour type, the
trace gases. Its rows say what the air is made of in words. This module turns the drawn regimes into mole fractions of the gases
of the table (`physics/gases.py`) by chemistry and not by a per-class template: the element budget of the body, the temperature
and pressure of its layer and the oxygen fugacity of its mantle go into the equilibrium solver (`physics/gas_equilibrium.py`, over
the NASA Glenn thermodynamics of `physics/species_thermo.py`), and the gas that comes out is the mixture. A class letter is never an
input to the chemistry: it only gates which rows of the recipe tables apply.

    draw = draw_mixture(Subject("M", {... the quantities required_inputs("M") names ...}))
    draw.mole_fractions       # {"H2O": 0.78, "CO2": 0.16, ...}, gas-table formulas, sums to 1
    draw.family               # "rocky", "giant", "brown-dwarf", "hot-jupiter", "envelope", ...
    draw.drawn                # {"a": "a6", "b": "b2", "c": "c3", ...}
    draw.redox_diw            # log10 fO2 above the iron-wuestite buffer (Hirschmann 2021 scale); None for a giant

It is a pure function of the subject and the seeded stream (`util/draw.py`): no pressure logic and no storage. The total pressure of
the air (the retention row `d`), the freeze-out of what the surface temperature condenses, the sum of the partial pressures and the
mean molar mass are the generator's (GEN.411, GEN.422); the clouds and hazes are GEN.414's and the ozone GEN.413's. The draws go
after the recipe draw in a fixed order: a change of the order changes every body and is a generator-epoch change.

The families, by what the recipe drew:

| the recipe drew | family | module |
|---|---|---|
| a brown-dwarf row (`g`) | `brown-dwarf` | `mixture_giant` |
| a hot-Jupiter row (`h`) | `hot-jupiter` | `mixture_giant` |
| a giant row (`f`: the classes I, J and T, the class being the authority on what is a giant; a stripped envelope `e1` of such a body is for the class reconciliation, GEN.645) | `giant` | `mixture_giant` |
| the envelope kept (`e2`) of a body without giant rows (a super-Earth) | `envelope` | `mixture_giant` (the metallicity from the drawn Z) |
| a coma row (`k`: a comet or an icy asteroid) | `coma` | `mixture_coma` |
| an evolved-terrestrial row (`j1` Venus-like, `j2` Mars-like, `j5` abiotic oxygen) | `venus`, `mars`, `abiotic-oxygen` | `mixture_evolved` |
| a cold-body row (`j3` Titan-like, `j4` Pluto-like) | `titan`, `frost` | `mixture_evolved` |
| a molten surface (a hot-zone lava row `i` with the dayside temperature `tp_k` at or above `mixture_lava.MELT_MIN_K`; the silicate-vapour exosphere `l1` of a white-dwarf or pulsar planet at the surface temperature) | `lava` | `mixture_lava` |
| airless white-dwarf planet (`l2`, or `l1` of a surface below `MELT_MIN_K`) | `airless` (no gas: empty `mole_fractions`) | |
| anything else rocky (`a` to `c`, `e1`, `e3`, `e4`) | `rocky` | `mixture_rocky` |
"""

import dataclasses

from planetgen.physics import atmosphere_recipes, mixture_coma, mixture_escape, mixture_evolved, mixture_giant, mixture_lava, mixture_rocky
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.physics.mixture_base import MixtureDraw, MixtureError, NotApplicable  # noqa: F401  (the public names of the result)


CATEGORIES = "".join(letter for letter in atmosphere_recipes.DRAW_ORDER if letter <= "l")
"""str: The recipe tables that decide the major gases, in the recipes' draw order (a to l). The trace and minor gases `m` read the
life stage (`o2_pal`), so they are drawn after it (`mixture_trace.draw_trace_gases`), and the clouds `n`, the radioisotopes `o` and the
ultraviolet regime `p` by their own lanes, all from the same stream."""


ROCKY_INPUTS = ("ts_k", "pressure_bar", "mass_earth", "age_gy", "v_esc_kms", "teq_k", "xuv_index")
"""tuple: What the rocky chemistry reads besides what the recipe rows read."""

GIANT_INPUTS = ("radius_km", "t_layer_k")
"""tuple: What the giant chemistry reads besides what the recipe rows read (the mass is `mass_jupiter` or `mass_earth`; `fe_h`, the
host's [Fe/H], is optional and 0 without it)."""


def required_inputs(body_class):
    """
    The quantities a subject of `body_class` may have to carry for `draw_mixture`: every stored quantity a row of the recipe tables
    reads for the class (in a weight, an intermediate, a distribution or a condition), less the ones an earlier row draws itself
    (`sets`) and the draws themselves, plus what the chemistry of the families the class can draw reads (`ROCKY_INPUTS`,
    `GIANT_INPUTS`). A row whose conditions fail does not read its quantities, so a given subject may need fewer; one that lacks a
    quantity a row does read is refused by name.

    Returns:
        tuple: Sorted names, all in `atmosphere_recipes.INPUTS`.
    """
    names, drawn, categories = set(), set(), set()
    for letter in CATEGORIES:
        for entry in RECIPES.in_category(letter):
            if atmosphere_recipes.pt.ANY_CLASS not in entry.classes and body_class not in entry.classes:
                continue
            categories.add(letter)
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
    chemistry = set()
    if "a" in categories or "e" in categories:
        chemistry |= set(ROCKY_INPUTS)
    if categories & set("efg"):
        chemistry |= set(GIANT_INPUTS)
    return tuple(sorted(((names - drawn) & set(atmosphere_recipes.INPUTS)) | chemistry))


def _surface_temperature(subject, name):
    value = subject.get(name)
    return subject.get("ts_k") if value is None else value


def family_of(recipe, subject):
    """The family the recipe's draw and the subject's surface temperature call for (see the module text)."""
    groups = recipe.groups
    if "g" in groups:
        return "brown-dwarf"
    if "h" in groups:
        return "hot-jupiter"
    if "f" in groups:
        return "giant"
    if "k" in groups:
        if groups["k"].name != "none":
            return "coma"
        if subject.body_class == "comet":
            return "airless"        # a comet beyond the coma's perihelion has no gas; an asteroid goes on to the rocky rows (GEN.674)
    if "l" in groups:
        if groups["l"].name == "l2":
            return "airless"
        if groups["l"].name == "l1":
            temperature = _surface_temperature(subject, "ts_k")
            return "lava" if temperature is not None and temperature >= mixture_lava.MELT_MIN_K else "airless"
    if "i" in groups:
        temperature = _surface_temperature(subject, "tp_k")
        if temperature is not None and temperature >= mixture_lava.MELT_MIN_K:
            return "lava"
    if "j" in groups and groups["j"].name in mixture_evolved.FAMILY_OF_ROW:
        return mixture_evolved.FAMILY_OF_ROW[groups["j"].name]
    if "e" in groups and groups["e"].name == "e2":
        return "envelope"
    return "rocky"


def draw_mixture(subject, recipe_subject=None, supply=None, escape=True, escape_column_bar=None):
    """
    The mixture of `subject`'s atmosphere.

    Args:
        subject (Subject): The body: its class and the quantities of `required_inputs(subject.body_class)`.
        recipe_subject (Subject, optional): The subject the recipe tables are drawn for, when it is not `subject`: a caller that solves the
            mixture again at another temperature and pressure (the fixed point of `surface_air`) keeps the regimes of one step by
            drawing them for that step's subject, so that a body on the edge of a gate does not change regime between two steps
            while the stream is consumed exactly as before. `subject` itself when not given.
        supply (callable, optional): `(inventory) -> {element: share}`, the share of each element of the delivered inventory
            that the surface holds (`physics/surface_pressure.py`), for the families that outgas (`mixture_rocky.volcanic_gas`);
            the whole inventory when not given.
        escape (bool): Whether the light-gas escape (`mixture_escape.apply`) takes its share from the mixture. A caller that iterates on
            the pressure (`surface_air`) leaves it off until the air has settled without it.
        escape_column_bar (float, optional): The pressure, bar, the escape compares the loss to when it is not the subject's.

    Returns:
        MixtureDraw: The mixture. Its `recipe` is the `RecipeDraw` of `CATEGORIES` that decided it (the retention regime `d`, the
            trace gases `m`, the values drawn), so that the pressure, the cloud and the colour code read the same draw and never
            draw the tables again.

    Raises:
        MixtureError: A quantity a row reads is missing (named), a family has no chemistry yet, or the gas holds more than
            `mixture_base.UNLISTED_LIMIT` of species the gas table lacks.
    """
    recipe = draw_recipe(subject if recipe_subject is None else recipe_subject)
    family = family_of(recipe, subject)
    drawn = {group: row.name for group, row in recipe.groups.items()}
    try:
        result = _draw_family(subject, recipe, family, supply)
    except NotApplicable:
        result = MixtureDraw(mole_fractions={}, family="airless", drawn=drawn, redox_diw=None,
                             temperature_k=float(subject.get("ts_k") or 0.0), pressure_bar=0.0, inventory={})
    if escape:
        result = mixture_escape.apply(subject, result, escape_column_bar)
    return dataclasses.replace(result, recipe=recipe)


def draw_recipe(subject):
    """The recipe draw (`RECIPES.draw_all` over `CATEGORIES`) that decides the family and the regimes of `subject`'s air."""
    try:
        return RECIPES.draw_all(subject, categories=CATEGORIES)
    except KeyError as error:
        raise MixtureError(str(error)) from error


def _draw_family(subject, recipe, family, supply):
    """The family's mixture; `NotApplicable` for a regime the body's temperature or supply does not allow, with the evolved families
    falling back to the rocky outgassing."""
    if family == "rocky":
        result = mixture_rocky.outgassed_mixture(subject, recipe, supply)
    elif family == "coma":
        result = mixture_coma.coma_mixture(subject, recipe)
    elif family == "airless":
        result = MixtureDraw(mole_fractions={}, family="airless", drawn={g: d.name for g, d in recipe.groups.items()}, redox_diw=None,
                             temperature_k=float(subject.get("ts_k") or 0.0), pressure_bar=0.0, inventory={})
    elif family == "lava":
        result = mixture_lava.vapour_mixture(subject, recipe, temperature_name="ts_k" if "l" in recipe.groups else "tp_k")
    elif family in mixture_evolved.FAMILIES:
        try:
            result = mixture_evolved.draw(subject, recipe, family, supply)
        except NotApplicable:
            result = mixture_rocky.outgassed_mixture(subject, recipe, supply)
    elif family == "envelope":
        z = recipe.values.get("e5")
        if z is None:
            raise MixtureError("the envelope row drew no metal mass fraction (e5)")
        result = mixture_giant.giant_mixture(subject, recipe, family="envelope",
                                             log_enrichment=mixture_giant.enrichment_of_mass_fraction(z))
    else:
        result = mixture_giant.giant_mixture(subject, recipe, family=family)
    return result
