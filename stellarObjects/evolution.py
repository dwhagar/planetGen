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
from .constants import EVOLUTIONARY_TIMELINES, STAR_EVOLUTION
from .utils import to_paragraph, _format_age_string # Import the new utility function

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
    spectral_class = star.type[0]
    star_info = STAR_EVOLUTION.get(spectral_class, {})
    
    # Ensure "normal" is a fallback if the star's supported_evolutionary_scales is empty or missing
    supported_scales = star_info.get("supported_evolutionary_scales")
    if not supported_scales:
        evolutionary_scale = "normal"
    else:
        # If there are multiple supported scales, pick one randomly
        evolutionary_scale = random.choice(supported_scales)

    timeline = EVOLUTIONARY_TIMELINES[evolutionary_scale]

    # star_lifespan is now directly a float in billions of years
    star_lifespan_value = timeline['star_lifespan']

    # The current system age should be the star's actual age, not a randomly generated one
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

    # Constructing the paragraph output
    output_parts = []

    output_parts.append(
        f"A speculative evolutionary timeline for a planet orbiting this star indicates a star lifespan of "
        f"{_format_age_string(star_lifespan_value)} and an evolutionary pace described as {timeline['evolutionary_pace'].lower()}. "
        f"The current estimated age of the system is {_format_age_string(current_system_age)}."
    )

    if most_recent_milestone_age > -1.0:
        output_parts.append(
            f"The most recent significant evolutionary milestone prior to this age would have been "
            f"{most_recent_milestone_name} at {_format_age_string(most_recent_milestone_age)}."
        )
    else:
        output_parts.append(
            f"No significant evolutionary milestones are predicted to have occurred yet at this system's age."
        )

    summary_text = (
        f"The star's spectral class of '{spectral_class}' suggests a {timeline['evolutionary_pace'].lower()} "
        f"evolutionary timeline. The star's lifespan of {_format_age_string(star_lifespan_value)} would likely have a "
        f"significant impact on the development of life on any orbiting planets. "
    )
    output_parts.append(summary_text)

    constraint_notes = star_info.get('evolutionary_constraint_notes', '')
    if constraint_notes:
        output_parts.append(constraint_notes) # Append the constraint notes directly

    return output_parts