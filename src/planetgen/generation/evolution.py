# planetgen/generation/evolution.py

"""
Evolutionary Timeline
=====================

What a world's life has reached at its star's age: the stage whose age on the world's schedule is the latest not past it, capped by what
its life class can reach (`physics.life_complexity`) and by a forbidden civilization (`-intelligent_life`; a request for one only raises the odds, GEN.1130). The schedule itself, a draw of five hard steps (GEN.750), is
`physics.life_timeline`; the oxygen the oxygenic ones make is `physics.biogenic_oxygen`; `generation.life.apply_life_data` joins them.
"""


from planetgen import tuning
from planetgen.physics import life_class_text, life_classes

MILESTONE_KEYS = tuning.LIFE_MILESTONES
"""
tuple[str]: The milestone keys, in ascending order of life complexity -- shared by
`get_evolutionary_state` for both display-name lookup and comparing a
milestone against a cap (`physics.life_complexity.top_stage`).
"""

RUNG_OF_MILESTONE = {"abiogenesis": 1, "photosynthesis": 2, "complex_cells": 3, "multicellularity": 4, "technological_civilization": 6}
"""dict: The rung of the stage ladder (`physics.life_class_text`, life-classes.md 5.1) each milestone stands for; a world whose multicellular life has
animals (an animal grade of 1 to 3) is at rung 5."""

SCHEDULE_FIELDS = tuning.LIFE_SCHEDULE_FIELDS
"""tuple[str]: The body columns holding a world's milestone schedule (DB.100), one age in Gyr for each `MILESTONE_KEYS`
milestone in order."""


def get_evolutionary_state(star, ages, max_stage=None):
    """
    The evolutionary state of a planet at its host star's age, on the schedule it was drawn (`physics.life_timeline.schedule`): the
    latest milestone reached and the ages (DB.51). This is the structure the stored columns hold; `describe/build.py` writes it as
    prose. It draws nothing.

    Args:
        star (Star): The host star.
        ages (tuple): The age, Gyr, of each `MILESTONE_KEYS` milestone (`None` for one the world never reaches).
        max_stage (str, optional): A `MILESTONE_KEYS` key capping the milestone, whether from the natural age-based
            roll below: the complexity of its life class allows a world only so
            far (`physics.life_complexity.top_stage`). The planet class letter caps nothing (GEN.647).

    Returns:
        dict: `stage` (a `MILESTONE_KEYS` key, `None` when no milestone has been reached), `milestone_age_gyr` (`None` with no stage)
            and `system_age_gyr`.
    """
    timeline = dict(zip(MILESTONE_KEYS, ages))

    # The current system age is always the star's actual, already-finalized age
    # (see Star.adjust_age_for_planets) -- never inflated past it, so this sentence
    # can never contradict the star's own stated age/lifespan elsewhere in the output.
    current_system_age = star.age

    # `max_stage` (the habitability index's ceiling) and INTELLIGENT_LIFE=False's pre-existing
    # "never report a civilization when one was explicitly disallowed" rule are both just caps on the
    # same milestone index -- combined here so the "most recent milestone"
    # search below only ever considers what BOTH allow, rather than
    # rolling unrestricted and correcting the result after the fact.
    max_index = len(MILESTONE_KEYS) - 1
    if max_stage:
        max_index = min(max_index, MILESTONE_KEYS.index(max_stage))
    if star.system_config.INTELLIGENT_LIFE is False:
        max_index = min(max_index, MILESTONE_KEYS.index("multicellularity"))
    allowed_keys = MILESTONE_KEYS[:max_index + 1]

    most_recent_key = None
    most_recent_milestone_age = -1.0

    # Find the most recent milestone prior to or at the current system age
    for key in allowed_keys:
        age = timeline[key]
        if age is not None and age <= current_system_age and age > most_recent_milestone_age:
            most_recent_milestone_age = age
            most_recent_key = key

    # Hard invariant, not just a silent clamp: whatever got selected above
    # must never exceed max_index. Structurally impossible given the
    # candidate pool was already truncated to allowed_keys.
    if most_recent_key is not None:
        assert MILESTONE_KEYS.index(most_recent_key) <= max_index, (
            f"get_evolutionary_state picked {most_recent_key!r}, exceeding its cap {MILESTONE_KEYS[max_index]!r}"
        )

    return {"stage": most_recent_key,
            "milestone_age_gyr": most_recent_milestone_age if most_recent_key is not None else None,
            "system_age_gyr": current_system_age}


def _rung(stage, grade):
    """`(rung, grade label or None)` of a milestone on a world whose animals are of grade `grade` (0 for none)."""
    if stage == "multicellularity" and grade:
        return 5, life_class_text.GRADES[min(int(grade), 3) - 1]
    return RUNG_OF_MILESTONE[stage], None


def stage_name(stage, class_id, grade=0):
    """The display name of a milestone on a world of life class `class_id` (`physics.life_class_text`): the class's own words for the rung,
    the animal grade's when its multicellular life has animals (`grade`, 1 to 3)."""
    rung, label = _rung(stage, grade)
    return life_class_text.stage_name(class_id, rung, label)


def stage_text(stage, class_id, grade=0):
    """What the biosphere is like at a milestone on a world of life class `class_id` (`physics.life_class_text`)."""
    rung, label = _rung(stage, grade)
    return life_class_text.stage_description(class_id, rung, label)


def releases_oxygen(class_id):
    """Whether a life class splits water for its electrons and so makes free oxygen (`life_classes.LifeClass.releases_oxygen`;
    `physics.biogenic_oxygen`)."""
    return bool(class_id) and life_classes.describe(class_id).releases_oxygen
