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
from .constants import EVOLUTIONARY_TIMELINES, STAR_EVOLUTION, EVOLUTIONARY_TEXT
from .utils import to_paragraph, _format_age_string
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

    # The current system age should be the star's actual age, not a randomly generated one
    current_system_age = star.age

    # If FORCE_INT is true, ensure the star's age is sufficient for technological civilization
    if star.system_config.FORCE_INT: # Use star.system_config
        tech_civ_age = timeline['technological_civilization']
        if current_system_age < tech_civ_age:
            # Nudge the star's age to be just enough for technological civilization
            # Add a small random increment to avoid exact boundary issues
            current_system_age = tech_civ_age + (random.random() * 0.1) # Add up to 0.1 billion years

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

    # Apply FORCE_INT and NO_INT logic
    if star.system_config.FORCE_INT: # Use star.system_config
        most_recent_milestone_name = "Technological Civilization"
        most_recent_milestone_age = milestones["Technological Civilization"]
    elif star.system_config.NO_INT: # Use star.system_config
        if most_recent_milestone_name == "Technological Civilization":
            most_recent_milestone_name = "Multicellularity"
            most_recent_milestone_age = milestones["Multicellularity"]

    # Constructing the paragraph output
    output_sentences = []

    output_sentences.append(
        f"A speculative evolutionary timeline for a planet orbiting this star indicates an evolutionary pace described as {timeline['evolutionary_pace'].lower()}. "
        f"The current estimated age of the system is {_format_age_string(current_system_age)}."
    )

    if most_recent_milestone_age > -1.0:
        output_sentences.append(
            f"The most recent significant evolutionary milestone prior to this age would have been "
            f"{most_recent_milestone_name} at {_format_age_string(most_recent_milestone_age)}. {EVOLUTIONARY_TEXT[most_recent_milestone_name.lower().replace(' ', '_')]}"
        )
    else:
        output_sentences.append(
            f"No significant evolutionary milestones are predicted to have occurred yet at this system's age."
        )

    # Join the sentences into a single paragraph
    return [to_paragraph(output_sentences)]