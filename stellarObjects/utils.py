# stellarObjects/utils.py

"""
Utilities
=========

This module contains utility functions for the planetGen package, including
mathematical calculations, string formatting, and name generation.
"""

import math
import random
from .constants import SOLAR_LUMINOSITY, STEFAN_BOLTZMANN_CONSTANT
from .names import VOWELS, BAD_CONSONANTS, AVOIDED_NAMES, NSFW_WORDS, WORD_SIZE_MEAN

def to_scientific_notation(number, precision=2):
    """
    Converts a number to scientific notation with the specified precision.

    Args:
        number (float): The number to convert.
        precision (int, optional): The number of decimal places to show. Defaults to 2.

    Returns:
        str: The number in scientific notation (e.g., "1.23e+04").
    """
    if number == 0:
        return "0"
    exponent = int(math.floor(math.log10(abs(number))))
    coefficient = number / 10**exponent
    output = f"Exp|{coefficient:.{precision}f}|{exponent:d}"
    return "{{" + output + "}}"


def years_to_time_string(years):
    """
    Converts a decimal number of years into a human-readable string
    like "x years y days z hours m minutes" (omitting any components with zero values).
    """
    total_minutes = round(years * 365.25 * 24 * 60)
    years = total_minutes // (365 * 24 * 60)
    remaining_minutes = total_minutes % (365 * 24 * 60)
    days = remaining_minutes // (24 * 60)
    remaining_minutes %= 24 * 60
    hours = remaining_minutes // 60
    minutes = remaining_minutes % 60

    time_parts = []
    if years > 0:
        time_parts.append(f"{years} year{'s' if years > 1 else ''}")
    if days > 0:
        time_parts.append(f"{days} day{'s' if days > 1 else ''}")
    if hours > 0:
        time_parts.append(f"{hours} hour{'s' if hours > 1 else ''}")
    if minutes > 0:
        time_parts.append(f"{minutes} minute{'s' if minutes > 1 else ''}")

    if len(time_parts) > 1:
        time_parts[-1] = f"and {time_parts[-1]}"
    return " ".join(time_parts)


def calc_object_mass(object_class, object_radius, PLANET_CLASSES, PLANET_DENSITY, object_density=None):
    """
    Calculates the mass of an object in kg.
    """
    if object_density is None:
        min_density, max_density = PLANET_DENSITY[PLANET_CLASSES[object_class]['type']]
        p_density = random.uniform(min_density, max_density)
    else:
        p_density = object_density

    volume = (4 / 3) * math.pi * (object_radius * 1000) ** 3
    mass = volume * p_density * 1000
    return volume, mass


def calculate_habitable_zone(luminosity):
    """
    Calculates the inner and outer boundaries of the habitable zone for a star.
    """
    solar_lum = luminosity / SOLAR_LUMINOSITY
    inner_radius = math.sqrt(solar_lum / 1.1)
    outer_radius = math.sqrt(solar_lum / 0.53)
    return (inner_radius, outer_radius)


def calculate_stellar_radius(luminosity, temperature):
    """
    Calculates the radius of a star in meters.
    """
    return math.sqrt(luminosity / (4 * math.pi * STEFAN_BOLTZMANN_CONSTANT * temperature ** 4))


def split_into_syllables(name):
    """
    Splits a word into a list of syllables.
    """
    syllables = []
    current_syllable = ""
    for i, char in enumerate(name):
        current_syllable += char
        if char in VOWELS and i < len(name) - 1 and name[i+1] not in VOWELS:
            syllables.append(current_syllable)
            current_syllable = ""
    if current_syllable:
        syllables.append(current_syllable)
    return syllables


def is_name_valid(name):
    """
    Checks if a generated name is valid.
    """
    name_lower = name.lower()
    if name_lower in AVOIDED_NAMES:
        return False
    for word in NSFW_WORDS:
        if word in name_lower:
            return False
    if len(name) > 1 and name[-1] in VOWELS and name[-2] in VOWELS:
        return False
    return True


def split_long_word(name):
    """
    Splits a long word into two, capitalizing the second word.
    """
    if len(name) > WORD_SIZE_MEAN:
        split_point = len(name) // 2
        return name[:split_point] + " " + name[split_point:].capitalize()
    return name


def generate_phoneme_salad_name(name_list, prefix_list, suffix_list):
    """
    Generates a unique, phonetically pleasing name from a list of base names.
    """
    while True:
        name = random.choice(name_list)
        
        syllables = split_into_syllables(name)
        if len(syllables) > 1:
            random.shuffle(syllables)
            name = "".join(syllables)

        prefix = random.choice(prefix_list)
        if prefix[-1] in VOWELS and name[0].lower() in VOWELS:
            name = prefix + name[1:]
        elif prefix[-1] not in VOWELS and name[0].lower() not in VOWELS:
            if (prefix[-1] + name[0].lower()) in BAD_CONSONANTS:
                name = prefix + "'" + name
            else:
                name = prefix + name
        else:
            name = prefix + name

        suffix = random.choice(suffix_list)
        if name[-1] in VOWELS and suffix[0].lower() in VOWELS:
            name = name + suffix[1:]
        elif name[-1] not in VOWELS and suffix[0].lower() not in VOWELS:
            if (name[-1] + suffix[0].lower()) in BAD_CONSONANTS:
                name = name + "'" + suffix
            else:
                name = name + suffix
        else:
            name = name + suffix
            
        name = name.lower()
        
        if is_name_valid(name):
            name = split_long_word(name)
            name = name[0].upper() + name[1:]
            if "'" in name:
                parts = name.split("'")
                name = "'".join([part[0].upper() + part[1:] for part in parts])
            return name


def to_paragraph(sentences):
    """
    Converts a list of sentences into a single, well-formed paragraph.
    """
    return " ".join(sentences)