# planetgen/physics/adaptation_draw.py

"""
The draw of a body's adaptation traits (GEN.876, `physics/adaptation_traits.py`): from the stress the body has, how likely each trait is and how far it goes.

The stress score of a trait is `s = clip(log10(stress / onset) / log10(full / onset), 0, 1)`: 0 where the ordinary limit of the index starts to cost the
life something (`ONSET`), 1 where the extreme grade is needed (`FULL`, the onset times the extreme grade's factor). A trait is drawn with chance
`P_CAP * s` and is extreme with chance `GRADE3_SHARE * s^2`, otherwise tolerant. Two uniform draws per stress, always, so six per body whatever it holds:
the caller binds the stream of the body's own (`draw.Stream`), so no other draw of any body moves when a trait's numbers change.
"""

import math

from planetgen.physics import adaptation_traits, atmosphere, habitability, habitability_world

ONSET = {"radiation": habitability.DOSE_THRESH_GY, "sulfide": habitability.H2S_TOX_KPA, "carbon_dioxide": habitability.CO2_TOX_KPA}
"""dict: The stress at which each ordinary limit starts to cost the life: the dose an animal takes (Gy/yr), the hydrogen sulfide and the carbon dioxide
scales of the toxins of the air (kPa). Tier T7, tweak level L3 (registry LF69)."""

FULL = {"radiation": habitability.MICROBE_DOSE_THRESH_GY * adaptation_traits.MICROBE_FACTOR["radiation"][3],
        "sulfide": habitability.H2S_TOX_KPA * adaptation_traits.ANIMAL_FACTOR["sulfide"][3],
        "carbon_dioxide": habitability.CO2_TOX_KPA * 10.0}
"""dict: The stress at which the extreme grade is needed: the extreme microbe threshold for the dose, the extreme animal limit for sulfide, and ten times the
carbon dioxide scale (the extreme grade of the design's table; an animal's tops out at 5 times). Tier T7, tweak level L3 (registry LF69)."""


def stresses(values):
    """`{stress: the body's stress}` from a body's column-named `values` (`habitability_world.values_from_body`): the absorbed dose where the life is (Gy/yr)
    and the partial pressures of hydrogen sulfide and carbon dioxide (kPa)."""
    return {"radiation": habitability_world.habitat_dose_sv_yr(values) / habitability.SV_PER_GY,
            "sulfide": atmosphere.partial_kpa(values, "h2s"),
            "carbon_dioxide": atmosphere.partial_kpa(values, "co2")}


def stress_score(stress, amount):
    """0 to 1: how hard `amount` of `stress` presses (see the module docstring)."""
    if amount <= ONSET[stress]:
        return 0.0
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
