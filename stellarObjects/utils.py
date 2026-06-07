# stellarObjects/utils.py

"""
Utilities
=========

This module contains utility functions for the planetGen package, including
mathematical calculations, string formatting, and name generation. These
functions support various aspects of the celestial body generation process,
from calculating physical properties to creating unique and plausible names.
"""

import math
import random
from .constants import SOLAR_LUMINOSITY, STEFAN_BOLTZMANN_CONSTANT
from .names import VOWELS, BAD_CONSONANTS, DICTIONARY_WORDS, NSFW_WORDS, WORD_SIZE_MEAN
from . import config

def to_scientific_notation(number, precision=2):
    """
    Converts a number to scientific notation with the specified precision.

    This function is used for formatting large numbers in a compact and
    standardized way, suitable for data templates and displays. It checks the
    global `config.MARKDOWN` flag to determine the output format.

    Args:
        number (float): The number to convert.
        precision (int, optional): The number of decimal places to show. Defaults to 2.

    Returns:
        str: The number in scientific notation, either as a wikitext template
             or as HTML for markdown.
    """
    if number == 0:
        return "0"
    exponent = int(math.floor(math.log10(abs(number))))
    coefficient = number / 10**exponent
    if config.MARKDOWN:
        return f"{coefficient:.{precision}f} × 10<sup>{exponent}</sup>"
    else:
        output = f"Exp|{coefficient:.{precision}f}|{exponent:d}"
        return "{{" + output + "}}"


def years_to_time_string(years):
    """
    Converts a decimal number of years into a human-readable string.

    The output format is "x years y days z hours m minutes", with any zero-value
    components omitted for brevity. This provides a more intuitive representation
    of orbital periods.

    Args:
        years (float): The number of years to convert.

    Returns:
        str: A human-readable string representing the time duration.
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
    Calculates the mass of a celestial object in kilograms.

    This function computes the mass based on the object's radius and density.
    If the density is not provided, it is randomly determined based on the
    object's class and type.

    Args:
        object_class (str): The class of the object (e.g., 'M', 'N').
        object_radius (float): The radius of the object in kilometers.
        PLANET_CLASSES (dict): A dictionary defining the properties of planet classes.
        PLANET_DENSITY (dict): A dictionary of density ranges for planet types.
        object_density (float, optional): The density of the object in g/cm³.

    Returns:
        tuple: A tuple containing the volume in km³ and the mass in kg.
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

    The habitable zone is defined as the region around a star where liquid
    water could exist on a planet's surface. This calculation is based on the
    star's luminosity.

    Args:
        luminosity (float): The luminosity of the star in Watts.

    Returns:
        tuple: A tuple containing the inner and outer radii of the habitable
               zone in AU.
    """
    solar_lum = luminosity / SOLAR_LUMINOSITY
    inner_radius = math.sqrt(solar_lum / 1.1)
    outer_radius = math.sqrt(solar_lum / 0.53)
    return (inner_radius, outer_radius)


def calculate_stellar_radius(luminosity, temperature):
    """
    Calculates the radius of a star in meters.

    This function uses the Stefan-Boltzmann law to calculate the star's radius
    from its luminosity and surface temperature.

    Args:
        luminosity (float): The luminosity of the star in Watts.
        temperature (float): The surface temperature of the star in Kelvin.

    Returns:
        float: The radius of the star in meters.
    """
    return math.sqrt(luminosity / (4 * math.pi * STEFAN_BOLTZMANN_CONSTANT * temperature ** 4))


def calculate_hill_sphere(distance_m, body_mass_kg, central_mass_kg):
    """
    Calculates the Hill sphere radius for a celestial body.

    The Hill sphere is the region around a celestial body where its own
    gravity is the dominant force for attracting satellites. This function
    calculates the radius of this sphere.

    Args:
        distance_m (float): The distance (semi-major axis) between the smaller
                            body and the larger central body, in meters.
        body_mass_kg (float): The mass of the smaller body (e.g., a planet) in kilograms.
        central_mass_kg (float): The mass of the larger central body (e.g., a star) in kilograms.

    Returns:
        float: The radius of the Hill sphere in meters.
    """
    return distance_m * (body_mass_kg / (3 * central_mass_kg)) ** (1 / 3)


def split_into_syllables(name):
    """
    Splits a word into a list of syllables.

    This is a basic syllable splitting function that helps in the process of
    generating new names by rearranging syllables from existing names.

    Args:
        name (str): The word to be split into syllables.

    Returns:
        list: A list of strings, where each string is a syllable.
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
    Checks if a generated name is valid based on multiple criteria.

    This function ensures that the generated name is not a common English word,
    does not contain any offensive terms, and follows basic phonetic rules.
    The validation checks are as follows:
    - The name should not exist in the NLTK dictionary of words.
    - The name should not contain any substring from the NSFW (Not Safe For Work) word list.
    - The name should not contain more than two consecutive vowels.
    - The name should not contain any of the defined bad consonant clusters.

    These checks help in generating names that are unique, appropriate, and sound plausible.

    Args:
        name (str): The name to validate. The function expects a lowercase string.

    Returns:
        bool: Returns `True` if the name is valid according to all the rules,
              otherwise returns `False`.
    """
    name_lower = name.lower()
    if name_lower in DICTIONARY_WORDS:
        return False
    for word in NSFW_WORDS:
        if word in name_lower:
            return False
    
    vowel_count = 0
    for char in name_lower:
        if char in VOWELS:
            vowel_count += 1
        else:
            vowel_count = 0
        if vowel_count > 2:
            return False

    for cluster in BAD_CONSONANTS:
        if cluster in name_lower:
            return False

    return True


def split_long_word(name):
    """
    Splits a long word into two, capitalizing the second word.

    This function improves the readability of long generated names by splitting
    them into two parts, creating a compound name effect.

    Args:
        name (str): The long word to be split.

    Returns:
        str: The split and capitalized name, or the original name if not long enough.
    """
    if len(name) > WORD_SIZE_MEAN:
        split_point = len(name) // 2
        return name[:split_point] + " " + name[split_point:].capitalize()
    return name


def generate_phoneme_salad_name(name_list, prefix_list, suffix_list):
    """
    Generates a unique, phonetically pleasing name from a list of base names.

    This function creates new names by taking a base name, shuffling its
    syllables, and adding a prefix and suffix. It includes logic to ensure
    the resulting name is phonetically plausible and passes validation checks.

    Args:
        name_list (list): A list of base names to choose from.
        prefix_list (list): A list of possible prefixes.
        suffix_list (list): A list of possible suffixes.

    Returns:
        str: A newly generated, unique name.
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

    This function is a simple utility to combine a list of strings into a
    single paragraph, which is used for generating descriptive text for
    celestial bodies.

    Args:
        sentences (list): A list of strings, where each string is a sentence.

    Returns:
        str: A single string representing the combined paragraph.
    """
    return " ".join(sentences)