# stellarObjects/planetLife.py

"""
Planet Life Chemistry and Evolutionary Data
=============================================

This module determines whether life could plausibly arise on a given planet
or moon, which biochemistry it would be based on, how fast its biosphere
would evolve, and (for habitable, non-moon worlds) its evolutionary
timeline narrative.

Life data is intentionally NOT computed during `Planet` construction (see
`planetPhysics.generate_planet_properties`) — every planet and moon is
generated in a "no life data yet" state. Instead, `StarSystem.__init__`
calls `apply_life_data` on every planet and moon in the system, after all
of them have been generated and after `Star.adjust_age_for_planets` has
finalized the star's age, so evolutionary timelines are computed against
the star's final age rather than a provisional pre-adjustment one.
"""

import random
import secrets
from stellarObjects import constants
from .utils import get_star_spectral_class, reseed_rng
from .evolution import get_evolutionary_timeline


def get_viable_life_chemicals(planet, spectral_class=None):
    """
    Determines the viable life chemicals for a planet by finding the
    intersection of chemicals supported by both the planet's class and
    its star's spectral type. Normalizes the probabilities to sum to 1.0.

    Args:
        planet (Planet): The planet or moon to evaluate.
        spectral_class (str, optional): The star's spectral class, if the
                                        caller already computed it (e.g.
                                        `apply_life_data`, which also calls
                                        `get_evolutionary_speed` for the same
                                        star). Computed from `planet.star`
                                        if not given.

    Returns:
        dict: A dictionary mapping viable life chemical strings to their
              normalized float probabilities (e.g., {"Melanin": 0.6}).
    """
    # Retrieve the base list of possible chemicals for this specific planet class
    planet_data = constants.PLANET_CLASSES.get(planet.planet_class)
    if not planet_data or not planet_data.get("life_chemical"):
        return {}

    planet_chems = planet_data["life_chemical"]

    if not planet.star_type:
        return {}

    # Determine the star's spectral class (e.g., 'G' from 'G2V')
    spectral_class = spectral_class or get_star_spectral_class(planet.star)

    # Retrieve the potentially viable chemicals for the star's evolutionary scale
    star_data = constants.STAR_EVOLUTION.get(spectral_class)
    if not star_data or not star_data.get("supported_evolutionary_scales"):
        return {}

    star_chems = star_data["potentially_viable_chemicals"]

    raw_probabilities = {}
    total_raw_prob = 0

    # Intersect the lists using substring matching and collect raw probabilities
    for p_chem in planet_chems:
        if any(p_chem in s_chem for s_chem in star_chems):
            chem_prob = constants.LIFE_CHEMICALS.get(p_chem, {}).get(
                "star_spectra_probabilities", {}).get(spectral_class, 0)

            if chem_prob > 0:
                raw_probabilities[p_chem] = chem_prob
                total_raw_prob += chem_prob

    # Normalize the probabilities proportionally so they sum to 1.0
    normalized_chemicals = {}
    if total_raw_prob > 0:
        for chem, prob in raw_probabilities.items():
            normalized_chemicals[chem] = prob / total_raw_prob

    return normalized_chemicals


def get_evolutionary_speed(planet, spectral_class=None):
    """
    Determines the evolutionary speed for a planet's biosphere.
    Finds the intersection of speeds supported by the star's lifespan
    and the speed required by the planet's specific life chemical, then
    chooses one randomly.

    Args:
        planet (Planet): The planet or moon to evaluate. Its `life_chemical`
                         should already be set (by `apply_life_data`) if a
                         chemical-specific speed is to be considered.
        spectral_class (str, optional): The star's spectral class, if the
                                        caller already computed it. Computed
                                        from `planet.star` if not given.

    Returns:
        str: The chosen evolutionary speed (e.g., 'fast', 'normal', 'slow'),
             or None if data is missing.
    """
    if not planet.star_type:
        return None

    # Determine the star's spectral class (e.g., 'G' from 'G2V')
    spectral_class = spectral_class or get_star_spectral_class(planet.star)

    # 1. Get the speeds supported by the star
    star_data = constants.STAR_EVOLUTION.get(spectral_class)
    if not star_data or not star_data.get("supported_evolutionary_scales"):
        return None

    star_speeds = star_data["supported_evolutionary_scales"]

    # 2. Get the speeds supported by the chemical (if one is assigned)
    if planet.life_chemical:
        chem_data = constants.LIFE_CHEMICALS.get(planet.life_chemical)
        if chem_data and chem_data.get("evolutionary_time_scale"):
            chem_speeds = chem_data["evolutionary_time_scale"]

            # Ensure it's a list for intersection logic, even though
            # constants.py currently stores it as a single string
            if isinstance(chem_speeds, str):
                chem_speeds = [chem_speeds]

            # Find the intersection
            valid_speeds = [speed for speed in star_speeds if speed in chem_speeds]

            if valid_speeds:
                return secrets.choice(valid_speeds)

    # 3. Fallback: If no chemical is set (or if there's somehow no overlap),
    # just pick a random speed supported by the star.
    if star_speeds:
        return secrets.choice(star_speeds)

    return None


def apply_life_data(planet):
    """
    Determines and assigns a planet's (or moon's) life chemistry,
    evolutionary speed, and — for habitable, non-moon worlds — its
    evolutionary timeline narrative.

    Called once per planet/moon by `StarSystem.__init__`, after all planets
    and moons in the system have been generated and the star's age has been
    finalized by `Star.adjust_age_for_planets`.

    Args:
        planet (Planet): The planet or moon to apply life data to. Its
                         `life_chemical`, `reflection_spectrum_visible`,
                         `reflection_spectrum_non_visible`, `evolution`, and
                         (if applicable) `evolutionary_data` attributes are
                         set in place.
    """
    reseed_rng()
    spectral_class = get_star_spectral_class(planet.star) if planet.star_type else None

    # Get the dictionary of viable chemicals and their weights
    viable_chems = get_viable_life_chemicals(planet, spectral_class)

    if viable_chems:
        # random.choices returns a list, so we grab the first [0] element
        planet.life_chemical = random.choices(
            population=list(viable_chems.keys()),
            weights=list(viable_chems.values()),
            k=1
        )[0]
        life_chem_data = constants.LIFE_CHEMICALS.get(planet.life_chemical, {})
        planet.reflection_spectrum_visible = life_chem_data.get("reflection_spectrum_visible")
        planet.reflection_spectrum_non_visible = life_chem_data.get("reflection_spectrum_non_visible")
    else:
        planet.life_chemical = None

    # Determine the evolutionary speed based on the star and chosen chemical
    planet.evolution = get_evolutionary_speed(planet, spectral_class)

    # Generate evolutionary timeline data if the planet is habitable and not a moon
    if planet.zone == 'e' and not planet.is_moon: # Only for habitable planets
        planet.evolutionary_data = get_evolutionary_timeline(planet.star)
