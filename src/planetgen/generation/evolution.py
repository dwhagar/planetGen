# planetgen/generation/evolution.py

"""
Evolutionary Timeline Generation
================================

This module contains functions for generating a speculative evolutionary
timeline for a planet based on its host star's spectral class. The timelines
are based on a set of predefined milestones and are adjusted based on the
star's type, with hotter, more massive stars having accelerated timelines and
cooler, less massive stars having decelerated timelines.

The timelines are purely speculative and are intended for creative world-building
purposes. They are not based on any established scientific models of astrobiology.
"""


from planetgen import tuning
from planetgen.util import draw
from planetgen.util import log
from planetgen.generation.stellar_object import get_star_evolutionary_profile

MILESTONE_KEYS = ("abiogenesis", "photosynthesis", "complex_cells", "multicellularity", "technological_civilization")
"""
tuple[str]: `EVOLUTIONARY_TIMELINES`/`EVOLUTIONARY_TEXT`'s own milestone
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


def draw_scale(star, chemical_scales=None):
    """
    Draws a world's evolutionary pace once, as a `tuning.EVOLUTIONARY_TIMELINES` key. Everything the world says
    about its pace (the life sentence, the timeline, the stored column) is this one value.

    The pace is one the star's lifespan supports and, when the world has a life chemistry, one that chemistry
    needs (`chemical_scales`, `LIFE_CHEMICALS[...]["evolutionary_time_scale"]`; a star that supports none of them
    decides alone). A forced technological civilization (`INTELLIGENT_LIFE`) draws among the paces whose
    civilization fits inside the star's actual, already-finalized age (the fastest one when none does), so the
    milestone never lies in the future.

    Args:
        star (Star): The host star.
        chemical_scales (list[str], optional): The paces the life chemistry allows.

    Returns:
        str: `"fast"`, `"normal"` or `"slow"`.
    """
    # Uses get_star_evolutionary_profile rather than a raw STAR_EVOLUTION[spectral_class]
    # lookup so giants/supergiants/white dwarfs etc. (whose spectral letter reflects
    # only current temperature, not a main-sequence lifespan) are handled correctly.
    star_info = get_star_evolutionary_profile(star)

    # "normal" is the fallback if the star's supported_evolutionary_scales is empty or missing
    supported_scales = star_info.get("supported_evolutionary_scales")
    if supported_scales and chemical_scales:
        supported_scales = [scale for scale in supported_scales if scale in chemical_scales] or supported_scales
    if not supported_scales:
        evolutionary_scale = "normal"
    elif star.system_config.INTELLIGENT_LIFE is True:
        # Pick randomly among whichever supported scales are actually reachable at that age (so
        # "fast" isn't always favored just because it's usually reachable --
        # a G-type star forcing intelligent life should still often read as
        # "Standard" pace, not always "Hyper-Accelerated"). If none
        # fit, fall back to the fastest supported scale, and the milestone
        # age is capped to the star's current age when reporting it.
        scale_speed_order = [s for s in ["fast", "normal", "slow"] if s in supported_scales]
        reachable_scales = [
            s for s in scale_speed_order
            if tuning.EVOLUTIONARY_TIMELINES[s]['technological_civilization'] <= star.age
        ]
        if reachable_scales:
            evolutionary_scale = draw.choice(reachable_scales)
            log.choice("Evolutionary scale", evolutionary_scale,
                       f"forced intelligent life: uniform draw among reachable scales {reachable_scales} "
                       f"(fit within star age {star.age:.4g})")
        else:
            evolutionary_scale = scale_speed_order[0]
            log.choice("Evolutionary scale", evolutionary_scale,
                       f"forced intelligent life: no scale in {scale_speed_order} was reachable by star "
                       f"age {star.age:.4g}, falling back to the fastest supported one")
    else:
        evolutionary_scale = draw.choice(supported_scales)
        log.choice("Evolutionary scale", evolutionary_scale,
                   f"uniform draw among supported scales {supported_scales}")
    return evolutionary_scale


def get_evolutionary_state(star, scale, max_stage=None):
    """
    The speculative evolutionary state of a planet at its host star's age, at the pace `draw_scale` drew for it:
    the latest milestone reached and the ages (DB.51). This is the structure the stored columns hold;
    `describe/build.py` writes it as prose. It draws nothing.

    Args:
        star (Star): The host star.
        scale (str): The pace, a `tuning.EVOLUTIONARY_TIMELINES` key (`draw_scale`).
        max_stage (str, optional): A `MILESTONE_KEYS` key capping the milestone, whether from the natural age-based
            roll below or a forced `INTELLIGENT_LIFE=True`: the habitability index allows a world only microbial life
            (GEN.206), so its timeline stops at "photosynthesis" (`life.stage_ceiling`). The class letter caps nothing
            (GEN.647).

    Returns:
        dict: `scale` (the one given), `stage` (a `MILESTONE_KEYS` key, `None` when no milestone has been reached),
            `milestone_age_gyr` (`None` with no stage) and `system_age_gyr`.
    """
    timeline = tuning.EVOLUTIONARY_TIMELINES[scale]

    # The current system age is always the star's actual, already-finalized age
    # (see Star.adjust_age_for_planets) — never inflated past it, so this sentence
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
        if age <= current_system_age and age > most_recent_milestone_age:
            most_recent_milestone_age = age
            most_recent_key = key

    # Apply INTELLIGENT_LIFE=True: force the HIGHEST milestone `max_index`
    # above actually allows -- "Technological Civilization" unless the
    # habitability index caps it lower.
    if star.system_config.INTELLIGENT_LIFE is True:
        most_recent_key = allowed_keys[-1]
        # Capped at the current system age: even on the fastest supported scale, a
        # very young star may still be younger than that scale's own tech-civ
        # threshold, and the milestone can never be claimed to predate "now".
        most_recent_milestone_age = min(timeline[most_recent_key], current_system_age)

    # Hard invariant, not just a silent clamp: whatever got selected above
    # must never exceed max_index. Structurally impossible given the
    # candidate pool was already truncated to allowed_keys.
    if most_recent_key is not None:
        assert MILESTONE_KEYS.index(most_recent_key) <= max_index, (
            f"get_evolutionary_state picked {most_recent_key!r}, exceeding its cap {MILESTONE_KEYS[max_index]!r}"
        )

    return {"scale": scale, "stage": most_recent_key,
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
