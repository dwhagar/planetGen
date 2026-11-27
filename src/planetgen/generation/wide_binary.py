# planetgen/generation/wide_binary.py

"""
Wide (S-type) binary draws
==========================

The draws for the pair of an **S-type (wide) binary**: the widest separation its place keeps (`wide_pair_limit_au`, GEN.103), the separation
(`sample_wide_binary_separation_au`, Moe and Di Stefano's companions per decade of period, GEN.448) and the eccentricity
(`sample_wide_binary_eccentricity`, a thermal distribution capped at `tuning.WIDE_BINARY_ECCENTRICITY_MAX`). The pair's orbit once it is drawn is a node of
the system's orbit tree, and its read-only view is `generation.pair_orbit` (GEN.640 slice b replaced the `WideBinaryPair` class that carried a copy of it).
"""

import math

from planetgen import tuning
from planetgen.generation import companion_classes, place_terms
from planetgen.util import draw


def wide_pair_limit_au(system_config, primary, secondary):
    """
    The widest axis a wide pair of these two stars keeps in its place (GEN.103,
    `place_terms.wide_axis_limit_au`): the place is `system_config.PLACE_TERMS`, or the Sun's neighbourhood for a
    system with none; the masses are the stars' at birth and the age the primary's.
    """
    terms = getattr(system_config, "PLACE_TERMS", None) or place_terms.SOLAR_NEIGHBOURHOOD
    total_mass_solar = place_terms.born_mass_solar(primary) + place_terms.born_mass_solar(secondary)
    limit = place_terms.wide_axis_limit_au(terms, total_mass_solar, primary.age)
    return min(limit, tuning.WIDE_BINARY_SEPARATION_CEILING_AU)


def sample_wide_binary_separation_au(max_au, born_mass_solar=1.0, total_mass_solar=2.0, inner_boost=1.0):
    """
    Draws an S-type (wide) binary's separation (semi-major axis) from Moe and Di Stefano's companions per decade of
    period for a primary of this mass (`generation.companion_classes`, GEN.448), between
    `tuning.WIDE_BINARY_SEPARATION_MIN_AU` and `max_au` (the place's survival limit, `wide_pair_limit_au`; a limit under
    the lower bound gives the lower bound, which only a pair forced wide reaches). One uniform draw.

    Args:
        max_au (float): The widest axis the pair keeps, AU.
        born_mass_solar (float): The primary's birth mass, solar masses (it picks the distribution).
        total_mass_solar (float): The pair's birth mass in total (it sets the period of an axis).
        inner_boost (float): The weight of the 0.25 to 10 AU range after the metallicity trend (1 with none).

    Returns:
        float: A separation, in AU.
    """
    return companion_classes.sample_separation_au(born_mass_solar, total_mass_solar, max_au, inner_boost)


def sample_wide_binary_eccentricity():
    """
    Draws an S-type (wide) binary's orbital eccentricity from a "thermal"
    distribution, `f(e) = 2e`, capped at
    `tuning.WIDE_BINARY_ECCENTRICITY_MAX` -- see that constant's
    own docstring for why wide pairs (unlike the close/P-type pair) keep a
    realistic, generally non-zero eccentricity.

    Derivation: the thermal PDF `f(e) = 2e` restricted to `[0, e_max]` and
    renormalized is still exactly proportional to `e` (just rescaled), so
    its CDF is `F(e) = (e / e_max)^2` and the closed-form inverse-CDF
    sample is `e = e_max * sqrt(u)`, `u ~ Uniform(0, 1)` -- no rejection
    sampling needed.

    Returns:
        float: An eccentricity, in [0, `WIDE_BINARY_ECCENTRICITY_MAX`).
    """
    return tuning.WIDE_BINARY_ECCENTRICITY_MAX * math.sqrt(draw.random())
