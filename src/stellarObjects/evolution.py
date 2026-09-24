# stellarObjects/evolution.py

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

import random

from . import log, program_constants
from .utils import format_age_string, get_star_evolutionary_profile, to_paragraph
# Removed: from . import config # Import the config module

MILESTONE_KEYS = ("abiogenesis", "photosynthesis", "complex_cells", "multicellularity", "technological_civilization")
"""
tuple[str]: `EVOLUTIONARY_TIMELINES`/`EVOLUTIONARY_TEXT`'s own milestone
keys, in ascending order of life complexity -- shared by
`get_evolutionary_timeline` for both display-name lookup and comparing a
milestone against a `program_constants.PLANET_CLASS_MAX_LIFE_STAGE` cap.
"""

_MILESTONE_DISPLAY_NAMES = {
    "abiogenesis": "Abiogenesis",
    "photosynthesis": "Photosynthesis",
    "complex_cells": "Complex Cells",
    "multicellularity": "Multicellularity",
    "technological_civilization": "Technological Civilization",
}


def life_stage_from_paragraphs(paragraphs):
    """
    Recovers the milestone `get_evolutionary_timeline` chose from the
    paragraph it wrote -- the only place it's recorded, stored as-is in
    `planet_evolutionary_paragraphs`/`moon_evolutionary_paragraphs`.

    Args:
        paragraphs (list[str]): A body's evolutionary paragraphs.

    Returns:
        str or None: One of `MILESTONE_KEYS`, or `None` when no milestone
            has been reached (or the body has no timeline at all).
    """
    for paragraph in paragraphs:
        for key, display_name in _MILESTONE_DISPLAY_NAMES.items():
            if f"would have been {display_name} at " in paragraph:
                return key
    return None


def get_evolutionary_timeline(star, planet_class=None):
    """
    Generates a speculative evolutionary timeline for a planet based on its
    host star's spectral class.

    This function takes a `Star` object as input and uses its spectral class
    to determine the appropriate evolutionary timeline. The timelines are
    defined in the `EVOLUTIONARY_TIMELINES` constant and are categorized as
    "Fast," "Earth Norm," or "Slow." The function retrieves the appropriate
    timeline and formats it into a human-readable string.

    The function also includes a summary of the star's characteristics and
    its likely impact on the evolution of life. The output is formatted as a
    simple paragraph string, suitable for wikitext or markdown.

    Args:
        star (Star): The star for which to generate an evolutionary timeline.
        planet_class (str, optional): The planet's own class code (e.g.
            "G"). When it appears in `program_constants.
            PLANET_CLASS_MAX_LIFE_STAGE`, the milestone this function can
            report -- whether from the natural age-based roll below or a
            forced `INTELLIGENT_LIFE=True` -- is capped at that class's own
            ceiling, so the returned narrative never contradicts what the
            class's own generated description already committed to (e.g. a
            Class G planet, "simple life", can never report a
            "Technological Civilization"). `None`/an unlisted class stays
            uncapped, matching this function's previous behavior exactly.

    Returns:
        list: A list of strings, where each string is a paragraph describing
              the speculative evolutionary timeline for a planet orbiting the star.
    """
    # Uses get_star_evolutionary_profile rather than a raw STAR_EVOLUTION[spectral_class]
    # lookup so giants/supergiants/white dwarfs etc. (whose spectral letter reflects
    # only current temperature, not a main-sequence lifespan) are handled correctly.
    star_info = get_star_evolutionary_profile(star)

    # Ensure "normal" is a fallback if the star's supported_evolutionary_scales is empty or missing
    supported_scales = star_info.get("supported_evolutionary_scales")
    if not supported_scales:
        evolutionary_scale = "normal"
    elif star.system_config.INTELLIGENT_LIFE is True:
        # A forced technological civilization needs to fit inside the star's
        # actual (already-finalized, lifespan-capped) age. Pick randomly among
        # whichever supported scales are actually reachable at that age (so
        # "fast" isn't always favored just because it's usually reachable —
        # a G-type star forcing intelligent life should still often read as
        # "Standard" pace, not always "Hyper-Accelerated"), so the milestone
        # below never has to be claimed at a time later than "now". If none
        # fit, fall back to the fastest supported scale, and the milestone
        # age is capped to the star's current age when reporting it.
        scale_speed_order = [s for s in ["fast", "normal", "slow"] if s in supported_scales]
        reachable_scales = [
            s for s in scale_speed_order
            if program_constants.EVOLUTIONARY_TIMELINES[s]['technological_civilization'] <= star.age
        ]
        if reachable_scales:
            evolutionary_scale = random.choice(reachable_scales)
            log.choice("Evolutionary scale", evolutionary_scale,
                       f"forced intelligent life: uniform draw among reachable scales {reachable_scales} "
                       f"(fit within star age {star.age:.4g})")
        else:
            evolutionary_scale = scale_speed_order[0]
            log.choice("Evolutionary scale", evolutionary_scale,
                       f"forced intelligent life: no scale in {scale_speed_order} was reachable by star "
                       f"age {star.age:.4g}, falling back to the fastest supported one")
    else:
        # If there are multiple supported scales, pick one randomly
        evolutionary_scale = random.choice(supported_scales)
        log.choice("Evolutionary scale", evolutionary_scale,
                   f"uniform draw among supported scales {supported_scales}")

    timeline = program_constants.EVOLUTIONARY_TIMELINES[evolutionary_scale]

    # The current system age is always the star's actual, already-finalized age
    # (see Star.adjust_age_for_planets) — never inflated past it, so this sentence
    # can never contradict the star's own stated age/lifespan elsewhere in the output.
    current_system_age = star.age

    # A planet class's own ceiling (program_constants.PLANET_CLASS_MAX_LIFE_STAGE,
    # e.g. Class G is capped at "photosynthesis" -- "simple life"), and
    # INTELLIGENT_LIFE=False's pre-existing "never report a civilization
    # when one was explicitly disallowed" rule, are both just caps on the
    # same milestone index -- combined here so the "most recent milestone"
    # search below only ever considers what BOTH allow, rather than
    # rolling unrestricted and correcting the result after the fact.
    class_max_key = program_constants.PLANET_CLASS_MAX_LIFE_STAGE.get(planet_class)
    max_index = MILESTONE_KEYS.index(class_max_key) if class_max_key else len(MILESTONE_KEYS) - 1
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

    # Apply INTELLIGENT_LIFE=True: force the HIGHEST milestone this planet
    # class (and `max_index` above) actually allows -- for an uncapped
    # class this is still "Technological Civilization" exactly as before;
    # for a class capped below that (e.g. Class G, "simple life"), this
    # forces its own ceiling instead of overriding the class's own
    # generated description.
    if star.system_config.INTELLIGENT_LIFE is True:
        most_recent_key = allowed_keys[-1]
        # Capped at the current system age: even on the fastest supported scale, a
        # very young star may still be younger than that scale's own tech-civ
        # threshold, and the milestone can never be claimed to predate "now".
        most_recent_milestone_age = min(timeline[most_recent_key], current_system_age)

    # Hard invariant, not just a silent clamp: whatever got selected above
    # must never exceed max_index. Structurally impossible given the
    # candidate pool was already truncated to allowed_keys -- a violation
    # here means a future change reintroduced an unrestricted milestone
    # path somewhere above, and that's a real bug worth failing loudly on
    # rather than silently shipping a planet's life narrative that
    # contradicts its own class description again.
    if most_recent_key is not None:
        assert MILESTONE_KEYS.index(most_recent_key) <= max_index, (
            f"get_evolutionary_timeline picked {most_recent_key!r} for planet_class={planet_class!r}, "
            f"exceeding its own cap {MILESTONE_KEYS[max_index]!r}"
        )

    most_recent_milestone_name = (
        _MILESTONE_DISPLAY_NAMES[most_recent_key] if most_recent_key is not None
        else "No significant evolutionary milestone"
    )

    # Constructing the paragraph output
    output_sentences = []

    output_sentences.append(
        f"A speculative evolutionary timeline for a planet orbiting this star indicates a {timeline['evolutionary_pace'].lower()} evolutionary pace. "
        f"The current estimated age of the system is {format_age_string(current_system_age)}."
    )

    if most_recent_milestone_age > -1.0:
        output_sentences.append(
            f"The most recent significant evolutionary milestone prior to this age would have been "
            f"{most_recent_milestone_name} at {format_age_string(most_recent_milestone_age)}. {program_constants.EVOLUTIONARY_TEXT[most_recent_milestone_name.lower().replace(' ', '_')]}"
        )
    else:
        output_sentences.append(
            "No significant evolutionary milestones are predicted to have occurred yet at this system's age."
        )

    # Join the sentences into a single paragraph
    return [to_paragraph(output_sentences)]