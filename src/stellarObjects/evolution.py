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

from . import program_constants
from .utils import format_age_string, get_star_evolutionary_profile, to_paragraph
# Removed: from . import config # Import the config module

def get_evolutionary_timeline(star):
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
        evolutionary_scale = random.choice(reachable_scales) if reachable_scales else scale_speed_order[0]
    else:
        # If there are multiple supported scales, pick one randomly
        evolutionary_scale = random.choice(supported_scales)

    timeline = program_constants.EVOLUTIONARY_TIMELINES[evolutionary_scale]

    # The current system age is always the star's actual, already-finalized age
    # (see Star.adjust_age_for_planets) — never inflated past it, so this sentence
    # can never contradict the star's own stated age/lifespan elsewhere in the output.
    current_system_age = star.age

    milestones = {
        "Abiogenesis": timeline['abiogenesis'],
        "Photosynthesis": timeline['photosynthesis'],
        "Complex Cells": timeline['complex_cells'],
        "Multicellularity": timeline['multicellularity'],
        "Technological Civilization": timeline['technological_civilization'],
    }

    most_recent_milestone_name = "No significant evolutionary milestone"
    most_recent_milestone_age = -1.0

    # Find the most recent milestone prior to or at the current system age
    for name, age in milestones.items():
        if age <= current_system_age and age > most_recent_milestone_age:
            most_recent_milestone_age = age
            most_recent_milestone_name = name

    # Apply INTELLIGENT_LIFE tri-state logic
    if star.system_config.INTELLIGENT_LIFE is True:
        most_recent_milestone_name = "Technological Civilization"
        # Capped at the current system age: even on the fastest supported scale, a
        # very young star may still be younger than that scale's own tech-civ
        # threshold, and the milestone can never be claimed to predate "now".
        most_recent_milestone_age = min(milestones["Technological Civilization"], current_system_age)
    elif star.system_config.INTELLIGENT_LIFE is False:
        if most_recent_milestone_name == "Technological Civilization":
            most_recent_milestone_name = "Multicellularity"
            most_recent_milestone_age = milestones["Multicellularity"]

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
            f"No significant evolutionary milestones are predicted to have occurred yet at this system's age."
        )

    # Join the sentences into a single paragraph
    return [to_paragraph(output_sentences)]