# planetgen/physics/adaptation_draw.py

"""
The draw of a body's adaptation traits (GEN.876, `physics/adaptation_traits.py`): from the stress the body has, how likely each trait is and how far it goes.

The stress score of a trait is `s = clip(log10(stress / onset) / log10(full / onset), 0, 1)` (a straight line for a shortfall, below): 0 where the ordinary limit of the index starts to cost the
life something (`ONSET`), 1 where the extreme grade is needed (`FULL`, the onset times the extreme grade's factor). A trait is drawn with chance
`P_CAP * s` and is extreme with chance `GRADE3_SHARE * s^2`, otherwise tolerant. Two uniform draws per stress, always, so ten per body whatever it holds:
the caller binds the stream of the body's own (`draw.Stream`), so no other draw of any body moves when a trait's numbers change.

The two stresses of GEN.1725 (`cold`, `heat`) are a **shortfall** in kelvin and not a ratio: the stress is how far the temperature of the water the life lives in (`life_supply.habitat_temperature_k`: the surface's under an open ocean or on a dry surface, the melting point of ice under an ice lid, the rise with depth in an aquifer; TEST.800, the surface mean before it drew a cold trait for most of the chemical-class worlds, whose crust is warm under a frozen surface) is past the ordinary limit of
an oxygen-making life (`life_adaptations`, axes `t_cold` and `t_hot`), the onset is 0 (the ordinary limit itself) and the full stress is the whole width the table gives (from the ordinary
limit to the extreme grade's). They come last in the draw order, so the three earlier traits of every body are the ones they were.
"""

import math

from planetgen.physics import adaptation_traits, atmosphere, habitability, habitability_world, life_adaptations, life_supply

ONSET = {"radiation": habitability.DOSE_THRESH_GY, "sulfide": habitability.H2S_TOX_KPA, "carbon_dioxide": habitability.CO2_TOX_KPA, "cold": 0.0, "heat": 0.0}
"""dict: The stress at which each ordinary limit starts to cost the life: the dose an animal takes (Gy/yr), the hydrogen sulfide and the carbon dioxide
scales of the toxins of the air (kPa), and 0 K of shortfall for cold and heat. Tier T7, tweak level L3 (registry LF69, and AD1 for cold and heat)."""

FULL = {"radiation": habitability.MICROBE_DOSE_THRESH_GY * adaptation_traits.MICROBE_FACTOR["radiation"][3],
        "sulfide": habitability.H2S_TOX_KPA * adaptation_traits.ANIMAL_FACTOR["sulfide"][3],
        "carbon_dioxide": habitability.CO2_TOX_KPA * 10.0,
        "cold": abs(life_adaptations.ordinary("t_cold") - life_adaptations.widest("t_cold", class_id="1a")),
        "heat": abs(life_adaptations.widest("t_hot", class_id="1a") - life_adaptations.ordinary("t_hot"))}
"""dict: The stress at which the extreme grade is needed: the extreme microbe threshold for the dose, the extreme animal limit for sulfide, and ten times the
carbon dioxide scale (the extreme grade of the design's table; an animal's tops out at 5 times); for cold and heat, the width in K from the ordinary limit to the table's extreme grade
for the oxygen-making class 1a (20 K and 33 K). Tier T7, tweak level L3 (registry LF69, and AD1 for cold and heat)."""

SHORTFALL_STRESSES = ("cold", "heat")
"""tuple: The stresses whose amount is a shortfall in K past the ordinary limit, scored on a straight line from the onset to the full stress, not on a logarithm."""


def stresses(values):
    """`{stress: the body's stress}` from a body's column-named `values` (`habitability_world.values_from_body`): the absorbed dose where the life is (Gy/yr)
    and the partial pressures of hydrogen sulfide and carbon dioxide (kPa)."""
    temperature = life_supply.habitat_temperature_k(values)
    return {"radiation": habitability_world.habitat_dose_sv_yr(values) / habitability.SV_PER_GY,
            "sulfide": atmosphere.partial_kpa(values, "h2s"),
            "carbon_dioxide": atmosphere.partial_kpa(values, "co2"),
            "cold": max(0.0, life_adaptations.ordinary("t_cold") - temperature) if temperature else 0.0,
            "heat": max(0.0, temperature - life_adaptations.ordinary("t_hot")) if temperature else 0.0}


def stress_score(stress, amount):
    """0 to 1: how hard `amount` of `stress` presses (see the module docstring)."""
    if amount <= ONSET[stress]:
        return 0.0
    if stress in SHORTFALL_STRESSES:
        return min(1.0, (amount - ONSET[stress]) / (FULL[stress] - ONSET[stress]))
    return min(1.0, math.log10(amount / ONSET[stress]) / math.log10(FULL[stress] / ONSET[stress]))


def chance(stress, amount):
    """The chance that a body with `amount` of `stress` has its trait; never over `adaptation_traits.P_CAP`, and never falls as the stress rises."""
    return adaptation_traits.P_CAP * stress_score(stress, amount)


def draw(values, stream):
    """The traits (`{stress: grade}`) of a body with `values`, drawn from `stream` (a `draw.Stream`): two draws for each stress, in the order of
    `adaptation_traits.STRESSES`."""
    amounts = stresses(values)
    found = {}
    for stress in adaptation_traits.STRESSES:
        present, grade = stream.random(), stream.random()
        score = stress_score(stress, amounts[stress])
        if present < adaptation_traits.P_CAP * score:
            found[stress] = 3 if grade < adaptation_traits.GRADE3_SHARE * score * score else 2
    return found
