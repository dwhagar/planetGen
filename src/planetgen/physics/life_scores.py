# planetgen/physics/life_scores.py

"""
The four life scores of a planet or moon (GEN.649, Boss 2026-10-23): Microbial, Flora, Fauna and Sentient, each a whole
number from 0 (none) to `tuning.LIFE_SCORE_MAX` (abundant), worked out from the body's evolutionary path (the latest
milestone its timeline reached, `life_stage`) and the life class that lives on it. They replace the stored `phi_bio`, `phi_cpx` and `phi_tech`.

A kind of life needs a stage before it can score: nothing without a milestone, microbial from abiogenesis, flora from
the producers of a light-eating class, fauna from multicellularity and an animal grade, a civilization at the technological milestone. How
much of it there is comes from the world's life class (`physics/life_classes.py`, `life_complexity.py`): the likelihood (0 to 1) of
the kind of life, stepped to 1..9. Flora before multicellular life is mats and algae and stays at `tuning.LIFE_SCORE_SIMPLE_FLORA_MAX`; fauna never
exceeds its class's producer base; a civilization fills its surface over `tuning.LIFE_SCORE_SENTIENT_FULL_GYR` after its milestone. A
likelihood under `tuning.HABITABILITY_LIFE_SCORE_MIN` is none. The thresholds are DESIGN values (docs/design/habitability-index.md 9).
"""

import math

from planetgen import tuning

FIELDS = ("life_microbial", "life_flora", "life_fauna", "life_sentient")
"""tuple: The stored columns, in the order the page shows them."""

NAMES = ("Microbial", "Flora", "Fauna", "Sentient")
"""tuple: What each of `FIELDS` is called."""


def _from_likelihood(likelihood):
    """1..`LIFE_SCORE_MAX` from a likelihood in [0, 1], 0 below the line `tuning.HABITABILITY_LIFE_SCORE_MIN`."""
    if likelihood is None or likelihood < tuning.HABITABILITY_LIFE_SCORE_MIN:
        return 0
    return max(1, min(tuning.LIFE_SCORE_MAX, math.ceil(tuning.LIFE_SCORE_MAX * likelihood)))


def stage_index(stage):
    """The rank of a milestone key (0 abiogenesis to 4 technological civilization), -1 for no milestone."""
    return tuning.LIFE_MILESTONES.index(stage) if stage in tuning.LIFE_MILESTONES else -1


NATIVE_LIFE_MILESTONE = {"microbial": "abiogenesis", "macroscopic": "multicellularity",
                         "biosphere": "technological_civilization"}
"""dict: The milestone a `habitability_world.native_life` stage stands for when a body has no timeline of its own to
say (a hand-built or edited body): microbes are abiogenesis, macroscopic life multicellularity, a biosphere a
civilization."""


def stage_of(native_life, milestone, has_timeline):
    """
    The milestone `scores` rates a body at, from its `habitability_world.native_life` (`"none"` for no life), the
    latest milestone its timeline reached (`life_stage`) and whether it has a timeline at all. Every body with a life
    chemistry has one (Boss, 2026-10-23); one without stands at the milestone of its native life.
    """
    if native_life == "none":
        return None
    return milestone if has_timeline else NATIVE_LIFE_MILESTONE.get(native_life, "abiogenesis")


def food_stage(microbial, flora, fauna):
    """
    The `habitability.LIFE_STAGES` rung a person's food comes from (Boss, 2026-10-23): flora is always good food, flora and
    fauna together are better, and fauna is all animal life. No life scores is `"none"`, microbes alone `"microbial"`
    (chemosynthetic mats), flora or fauna alone `"macroscopic"` (a Green food rung), both `"biosphere"` (Blue).
    """
    if flora and fauna:
        return "biosphere"
    if flora or fauna:
        return "macroscopic"
    return "microbial" if microbial else "none"


def scores(microbial_likelihood, complex_likelihood, stage, milestone_age_gyr=None, system_age_gyr=None, light_producer=True, grade=0):
    """
    `(microbial, flora, fauna, sentient)` for a body whose life class gives microbial life `microbial_likelihood` and complex life
    `complex_likelihood` (both 0 to 1; `habitability_world.compute`), whose timeline reached `stage` (a `tuning.LIFE_MILESTONES` key or
    `None`) at `milestone_age_gyr` of a system `system_age_gyr` old.

    Flora is the producer's own and needs a class that eats light (`light_producer`): the chemical classes' mats, in the dark or under ice,
    are Microbial only (GEN.871). Fauna needs an animal grade (`grade`, 0 to 3, `life_complexity`) and never exceeds the producer base of its
    class (the Flora a light class feeds on, the Microbial score of a chemical one: a vent world has Fauna with no Flora) or what its grade
    shows (`tuning.LIFE_COMPLEXITY["fauna_cap"]`).
    """
    rank = stage_index(stage)
    if rank < 0:
        return 0, 0, 0, 0
    microbial = _from_likelihood(microbial_likelihood)
    if microbial == 0:
        return 0, 0, 0, 0
    flora = fauna = sentient = 0
    if light_producer and rank >= stage_index("photosynthesis"):
        flora = microbial
        if rank < stage_index("multicellularity"):
            flora = min(flora, tuning.LIFE_SCORE_SIMPLE_FLORA_MAX)
    if rank >= stage_index("multicellularity"):
        complex_score = _from_likelihood(complex_likelihood)
        if light_producer:
            flora = max(flora, complex_score)
        if grade >= 1:
            fauna = min(complex_score, flora if light_producer else microbial, tuning.LIFE_COMPLEXITY["fauna_cap"][min(grade, 3)])
    if rank >= stage_index("technological_civilization") and fauna > 0:
        spread = 1.0
        if milestone_age_gyr is not None and system_age_gyr is not None:
            spread = min(1.0, max(0.0, system_age_gyr - milestone_age_gyr) / tuning.LIFE_SCORE_SENTIENT_FULL_GYR)
        sentient = max(1, round(fauna * spread))
    return microbial, flora, fauna, sentient
