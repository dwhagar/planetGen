# planetgen/physics/life_scores.py

"""
The four life scores of a planet or moon (GEN.649, Boss 2026-10-23): Microbial, Flora, Fauna and Sentient, each a whole
number from 0 (none) to `tuning.LIFE_SCORE_MAX` (abundant), worked out from the body's evolutionary path (the latest
milestone its timeline reached, `life_stage`) and its surface conditions (the habitability index's `PHI_bio` and
`PHI_cpx`). They replace the stored `phi_bio`, `phi_cpx` and `phi_tech`.

A kind of life needs a stage before it can score: nothing without a milestone, microbial from abiogenesis, flora from
photosynthesis, fauna from multicellularity, a civilization at the technological milestone. How much of it there is
comes from the surface: the likelihood (0 to 1) of the kind of life, stepped to 1..9. Flora before multicellular
life is mats and algae and stays at `tuning.LIFE_SCORE_SIMPLE_FLORA_MAX`; fauna never exceeds the flora it feeds on; a
civilization fills its surface over `tuning.LIFE_SCORE_SENTIENT_FULL_GYR` after its milestone. A likelihood under
`tuning.HABITABILITY_LIFE_SCORE_MIN` is none. The thresholds are DESIGN values (docs/design/habitability-index.md 9).
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


def scores(phi_bio, phi_cpx, stage, milestone_age_gyr=None, system_age_gyr=None):
    """
    `(microbial, flora, fauna, sentient)` for a body whose surface gives microbial life `phi_bio` and complex life
    `phi_cpx` (both 0 to 1), whose timeline reached `stage` (a `tuning.LIFE_MILESTONES` key or `None`) at
    `milestone_age_gyr` of a system `system_age_gyr` old.
    """
    rank = stage_index(stage)
    if rank < 0:
        return 0, 0, 0, 0
    microbial = _from_likelihood(phi_bio)
    if microbial == 0:
        return 0, 0, 0, 0
    flora = fauna = sentient = 0
    if rank >= stage_index("photosynthesis"):
        flora = microbial
        if rank < stage_index("multicellularity"):
            flora = min(flora, tuning.LIFE_SCORE_SIMPLE_FLORA_MAX)
    if rank >= stage_index("multicellularity"):
        complex_score = _from_likelihood(phi_cpx)
        flora = max(flora, complex_score)
        fauna = min(complex_score, flora)
    if rank >= stage_index("technological_civilization") and fauna > 0:
        spread = 1.0
        if milestone_age_gyr is not None and system_age_gyr is not None:
            spread = min(1.0, max(0.0, system_age_gyr - milestone_age_gyr) / tuning.LIFE_SCORE_SENTIENT_FULL_GYR)
        sentient = max(1, round(fauna * spread))
    return microbial, flora, fauna, sentient
