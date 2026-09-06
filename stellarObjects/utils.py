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
import secrets

from .config import SystemConfig
from .names import BAD_CONSONANTS, DICTIONARY_WORDS, NSFW_WORDS, VOWELS, WORD_SIZE_MEAN
from . import physical_constants, program_constants

def reseed_rng():
    """
    Re-seeds the `random` module with a high-entropy seed from `secrets`.

    This is called at the start of most generation methods to avoid
    correlated/repeating sequences across successive calls.
    """
    random.seed(secrets.randbits(128))


def get_star_spectral_class(star):
    """
    Returns the uppercase spectral class character (e.g. 'G') for a star.

    Works for both a plain `Star` and a `BinaryStarProxy`, without importing
    `BinaryStarProxy` here (which would create a circular import) — a proxy
    is detected by duck-typing on its `_primary` attribute, and its primary
    star's type is used as the representative spectral class.

    Args:
        star: A `Star` or `BinaryStarProxy` instance.

    Returns:
        str: The uppercase spectral class character.
    """
    reference_star = star._primary if hasattr(star, '_primary') else star
    return reference_star.type[0].upper()


def get_star_evolutionary_profile(star):
    """
    Returns the `STAR_EVOLUTION`-shaped profile describing which life
    chemicals and evolutionary paces are plausible for planets orbiting
    `star`, correctly accounting for its Yerkes luminosity class.

    `STAR_EVOLUTION` is keyed by spectral letter (e.g. 'G'), which for a
    main-sequence (Yerkes 'V') star fully determines both its current
    temperature and its total lifespan -- so that class's fixed entry
    applies directly.

    For any evolved or remnant star (giants, supergiants, subgiants, bright
    giants, hypergiants, subdwarfs, white dwarfs), the spectral letter only
    reflects the star's *current* temperature/color, not a lifespan -- e.g.
    an "F5VII" white dwarf merely glows at F-like temperature; it did not
    live and die as an F-type main-sequence star. Using STAR_EVOLUTION['F']
    directly for such a star would misapply a main-sequence lifespan/scale
    table to an unrelated evolutionary history. `potentially_viable_chemicals`
    (driven by the star's current emission spectrum) is still taken from the
    current-temperature letter, since that part genuinely is about current
    color. `supported_evolutionary_scales`, which is really a proxy for how
    much time was available for a biosphere to develop, is instead derived
    from the star's own already-computed age/lifespan (already correctly
    yerkes-class-aware -- see `Star._calculate_initial_star_age_and_lifespan`),
    by checking which evolutionary paces could reach their technological
    civilization milestone within that time budget.

    Args:
        star: A `Star` or `BinaryStarProxy` instance.

    Returns:
        dict: A `STAR_EVOLUTION`-entry-shaped dict (with at least
              `potentially_viable_chemicals` and `supported_evolutionary_scales`
              keys), or `{}` if the spectral letter has no entry.
    """
    reference_star = star._primary if hasattr(star, '_primary') else star
    spectral_class_char = reference_star.type[0].upper()
    base_info = program_constants.STAR_EVOLUTION.get(spectral_class_char, {})
    if not base_info:
        return {}

    if reference_star.yerkes_class == "V":
        return base_info

    # A white dwarf's lifespan is infinite (it just cools forever), so "how
    # much time has been available so far" is better represented by its age.
    # Every other evolved class has a finite lifespan, used as-is.
    time_budget = reference_star.age if reference_star.lifespan == float('inf') else reference_star.lifespan

    reachable_scales = [
        scale for scale in ["fast", "normal", "slow"]
        if program_constants.EVOLUTIONARY_TIMELINES[scale]['technological_civilization'] <= time_budget
    ]
    if not reachable_scales:
        # Even the fastest pace doesn't fit -- still return it so callers have
        # something to work with, mirroring the fallback in get_evolutionary_timeline.
        reachable_scales = ["fast"]

    return {**base_info, "supported_evolutionary_scales": reachable_scales}


def format_length_km(system_config: SystemConfig, value_km, threshold, round_digits, scientific_precision=None):
    """
    Formats a length in kilometers, switching between a comma-grouped plain
    number and scientific notation based on a threshold.

    Args:
        system_config (SystemConfig): The shared SystemConfig object.
        value_km (float): The length to format, in kilometers.
        threshold (float): The value above which scientific notation is used.
        round_digits (int): Decimal places for the plain-number form (and the
                            default precision for the scientific-notation form).
        scientific_precision (int, optional): Decimal places for the
                                              scientific-notation form, if
                                              different from `round_digits`.

    Returns:
        str: The formatted length string, including the " km" unit.
    """
    if value_km <= threshold:
        return f"{round(value_km, round_digits):,} km"
    precision = scientific_precision if scientific_precision is not None else round_digits
    return f"{to_scientific_notation(system_config, value_km, precision)} km"


def format_relative_to_sol(system_config: SystemConfig, value, sol_constant, unit, low_percent_precision=2):
    """
    Formats a physical quantity (mass, luminosity, ...) with its value in
    scientific notation alongside a comparison to the Sol reference value:
    a percentage when the ratio is small, or a "×" multiplier when it's large.

    Args:
        system_config (SystemConfig): The shared SystemConfig object.
        value (float): The quantity's value, in SI units.
        sol_constant (float): The Sol-reference value for the same quantity
                              (e.g. `physical_constants.SOLAR_MASS_TO_KG`).
        unit (str): The unit label for the SI value (e.g. "kg", "W").
        low_percent_precision (int, optional): Decimal places used for the
                                               percentage when the ratio is
                                               below `PERCENT_SOL_THRESHOLD_LOW`.
                                               Defaults to 2.

    Returns:
        str: The formatted string, e.g. "1.23 × 10^30 kg (1.00× Sol)".
    """
    sol_val = value / sol_constant
    sci_notation = to_scientific_notation(system_config, value)
    if sol_val < program_constants.PERCENT_SOL_THRESHOLD_LOW:
        return f"{sci_notation} {unit} ({sol_val * program_constants.PERCENT_MULTIPLIER:.{low_percent_precision}f}% of Sol)"
    elif sol_val < program_constants.PERCENT_SOL_THRESHOLD_HIGH:
        return f"{sci_notation} {unit} ({sol_val * program_constants.PERCENT_MULTIPLIER:.1f}% of Sol)"
    else:
        return f"{sci_notation} {unit} ({sol_val:.1f}× Sol)"


def ly_to_milliparsecs(ly):
    """
    Converts a distance in light-years to milliparsecs (mpc).

    For the database persistence layer only -- sector-scale position/
    geometry columns (star_systems.position_x/y/z_mpc, sectors.edge_mpc;
    see stellarObjects/schema.sql) are stored in milliparsecs specifically,
    distinct from every other distance-shaped column in the schema (which
    is kilometers). Generation/physics code keeps its own native light-year
    units for sector geometry throughout and never calls this.

    Args:
        ly (float): The distance in light-years.

    Returns:
        float: The distance in milliparsecs.
    """
    au = ly * physical_constants.LY_TO_AU
    return au / physical_constants.AU_PER_MILLIPARSEC


def milliparsecs_to_ly(mpc):
    """
    Converts a distance in milliparsecs (mpc) back to light-years -- the
    inverse of `ly_to_milliparsecs`, for reconstructing live objects
    (native light-year sector geometry) from database rows.

    Args:
        mpc (float): The distance in milliparsecs.

    Returns:
        float: The distance in light-years.
    """
    au = mpc * physical_constants.AU_PER_MILLIPARSEC
    return au * physical_constants.AU_TO_LY


def to_scientific_notation(system_config: SystemConfig, number, precision=2):
    """
    Converts a number to scientific notation with the specified precision.

    This function is used for formatting large numbers in a compact and
    standardized way, suitable for data templates and displays. It checks the
    global `config.MARKDOWN` flag to determine the output format.

    Args:
        system_config (SystemConfig): The shared SystemConfig object.
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
    if system_config.MARKDOWN:
        return f"{coefficient:.{precision}f} × 10<sup>{exponent}</sup>"
    else:
        output = f"Exp|{coefficient:.{precision}f}|{exponent:d}"
        return "{{" + output + "}}"

def format_age_string(age_gy, precision=2):
    """
    Formats an age in billions of years (GY) into a human-readable string,
    dynamically choosing between "Million Years" and "Billion Years".

    Args:
        age_gy (float): The age in billions of years.
        precision (int, optional): The number of decimal places to show. Defaults to 2.

    Returns:
        str: A formatted string (e.g., "12.50 Billion Years" or "500 Million Years").
    """
    if age_gy >= 1.0:
        return f"{age_gy:.{precision}f} Billion Years"
    else:
        return f"{age_gy * 1000:.{precision}f} Million Years"

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


def calculate_object_mass(object_class, object_radius, planet_classes, planet_density, object_density=None):
    """
    Calculates the mass of a celestial object in kilograms.

    This function computes the mass based on the object's radius and density.
    If the density is not provided, it is randomly determined based on the
    object's class and type.

    Args:
        object_class (str): The class of the object (e.g., 'M', 'N').
        object_radius (float): The radius of the object in kilometers.
        planet_classes (dict): A dictionary defining the properties of planet classes.
        planet_density (dict): A dictionary of density ranges for planet types.
        object_density (float, optional): The density of the object in g/cm³.

    Returns:
        tuple: A tuple containing the volume in km³ and the mass in kg.
    """
    reseed_rng()
    if object_density is None:
        min_density, max_density = planet_density[planet_classes[object_class]['type']]
        p_density = random.uniform(min_density, max_density)
    else:
        p_density = object_density

    volume_km3 = (4 / 3) * math.pi * object_radius ** 3
    volume_m3 = volume_km3 * physical_constants.KM_TO_M_FACTOR ** 3
    mass = volume_m3 * p_density * 1000
    return volume_km3, mass


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
    solar_lum = luminosity / physical_constants.SOLAR_LUMINOSITY
    inner_radius = math.sqrt(solar_lum / 1.1)
    outer_radius = math.sqrt(solar_lum / 0.53)
    return (inner_radius, outer_radius)


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
    reseed_rng()
    while True:
        name = secrets.choice(name_list)
        
        syllables = split_into_syllables(name)
        if len(syllables) > 1:
            random.shuffle(syllables)
            name = "".join(syllables)

        prefix = secrets.choice(prefix_list)
        if prefix[-1] in VOWELS and name[0].lower() in VOWELS:
            name = prefix + name[1:]
        elif prefix[-1] not in VOWELS and name[0].lower() not in VOWELS:
            if (prefix[-1] + name[0].lower()) in BAD_CONSONANTS:
                name = prefix + "'" + name
            else:
                name = prefix + name
        else:
            name = prefix + name

        suffix = secrets.choice(suffix_list)
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


def properties_to_string(system_config: SystemConfig, properties, template_name, markdown_header=None, markdown_key_map=None):
    """
    Converts a dictionary of properties into either a Markdown table or a Wiki template block,
    based on the `system_config.MARKDOWN` flag.

    Args:
        system_config (SystemConfig): The shared SystemConfig object.
        properties (dict): A dictionary where keys are property names (str) in Wikitext format
                           and values are their corresponding values (str or any type convertible to str).
        template_name (str): The name of the template to use for Wiki format (e.g., "Planet Data").
        markdown_header (str, optional): An optional header to prepend to the Markdown table.
                                         Defaults to None.
        markdown_key_map (dict, optional): A dictionary mapping Wikitext property keys to their
                                            desired Markdown table header names. If a key is not
                                            found in this map, it will be converted from
                                            lower_snake_case to Title Case for Markdown.

    Returns:
        str: A string formatted as either a Markdown table or a Wiki template block.
    """
    output_lines = []

    if system_config.MARKDOWN:
        if markdown_header:
            output_lines.append(markdown_header)
        output_lines.extend([
            "| Property | Value |",
            "|---|---|"
        ])
        for prop_key, value in properties.items():
            # Determine the Markdown header name
            if markdown_key_map and prop_key in markdown_key_map:
                markdown_prop_name = markdown_key_map[prop_key]
            else:
                # Default to converting lower_snake_case to Title Case
                markdown_prop_name = prop_key.replace('_', ' ').title()
            output_lines.append(f"| {markdown_prop_name} | {value} |")
        return "\n".join(output_lines)
    else:
        output_lines.append(f"{{{{{template_name}")
        for prop_key, value in properties.items():
            # For Wikitext, the keys are already in the correct format
            output_lines.append(f"|{prop_key}={value}")
        output_lines.append("}}")
        return "\n".join(output_lines)