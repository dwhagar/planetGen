# planetgen/generation/evolution.py

"""
Evolutionary Timeline
=====================

What a world's life has reached at its star's age: the stage whose age on the world's schedule is the latest not past it, capped by what
the habitability index allows and by a forced or forbidden civilization. The schedule itself, a draw of five hard steps (GEN.750), is
`physics.life_timeline`; the oxygen the oxygenic ones make is `physics.biogenic_oxygen`; `generation.life.apply_life_data` joins them.
"""


from planetgen import tuning

MILESTONE_KEYS = tuning.LIFE_MILESTONES
"""
tuple[str]: `EVOLUTIONARY_TEXT`'s own milestone
keys, in ascending order of life complexity -- shared by
`get_evolutionary_state` for both display-name lookup and comparing a
milestone against a cap (`life.stage_ceiling`).
"""

MILESTONE_DISPLAY_NAMES = {
    "abiogenesis": "Abiogenesis",
    "photosynthesis": "Photosynthesis",
    "complex_cells": "Complex Cells",
    "multicellularity": "Multicellularity",
    "technological_civilization": "Technological Civilization",
}


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
            roll below or a forced `INTELLIGENT_LIFE=True`: the habitability index allows a world only microbial life
            (GEN.206), so its timeline stops at "photosynthesis" (`life.stage_ceiling`). The class letter caps nothing
            (GEN.647).

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

    # Apply INTELLIGENT_LIFE=True: force the HIGHEST milestone `max_index`
    # above actually allows -- "Technological Civilization" unless the
    # habitability index caps it lower.
    if star.system_config.INTELLIGENT_LIFE is True:
        most_recent_key = allowed_keys[-1]
        # Capped at the current system age: a timeline that could not be fitted into the world's age (it never
        # makes the oxygen) still can't claim a milestone that predates "now".
        forced_age = timeline[most_recent_key]
        most_recent_milestone_age = current_system_age if forced_age is None else min(forced_age, current_system_age)

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


def stage_name(stage, life_chemical):
    """The display name of a milestone (`MILESTONE_DISPLAY_NAMES`), the "photosynthesis" one as the chemistry's own
    way of harvesting energy calls it (`tuning.LIFE_METABOLISMS`)."""
    if stage == "photosynthesis":
        return tuning.LIFE_METABOLISMS[metabolism(life_chemical)]["stage_name"].title()
    return MILESTONE_DISPLAY_NAMES[stage]


def stage_text(stage, life_chemical):
    """What the biosphere is like at a milestone: `tuning.EVOLUTIONARY_TEXT`, the "photosynthesis" one by the
    chemistry's metabolism (`tuning.LIFE_METABOLISMS`)."""
    if stage == "photosynthesis":
        return tuning.LIFE_METABOLISMS[metabolism(life_chemical)]["text"]
    return tuning.EVOLUTIONARY_TEXT[stage]


def metabolism(life_chemical):
    """How a life chemistry harvests its energy (`tuning.LIFE_CHEMICALS[...]["metabolism"]`)."""
    return tuning.LIFE_CHEMICALS[life_chemical]["metabolism"]


def releases_oxygen(life_chemical):
    """Whether a life chemistry splits water for its electrons and so makes free oxygen (`metabolism` is `oxygenic`;
    `physics.biogenic_oxygen`)."""
    return bool(life_chemical) and metabolism(life_chemical) == "oxygenic"
