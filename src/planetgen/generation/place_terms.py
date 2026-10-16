# generation/place_terms.py

"""
What a place does to the stars drawn in it (GEN.103, docs/design/placement-by-population.md 3.1 and 3.3):
the close-companion factor that follows a star's own [Fe/H] and the widest axis a wide pair can keep.

Both are functions of numbers the sector profile already holds (the density, the populations' velocity
dispersions and the star's own [Fe/H] and age), so nothing here is a catalogue by place. `PlaceTerms` carries the
profile's part to the system draw on `SystemConfig.PLACE_TERMS`; a system drawn with no galaxy placement uses
`SOLAR_NEIGHBOURHOOD`.
"""

import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import constants


@dataclass(frozen=True)
class PlaceTerms:
    """
    The part of a sector profile the pair draws read.

    Attributes:
        mass_density_solar_pc3 (float): Stellar mass per cubic parsec at the place.
        dispersion_kms (float): The velocity dispersion of the field there, km/s.
        perturber_mass_solar (float): The mean mass of the stars and remnants that pass a pair.
    """
    mass_density_solar_pc3: float
    dispersion_kms: float
    perturber_mass_solar: float = tuning.MEAN_SYSTEM_MASS_SOLAR


SOLAR_NEIGHBOURHOOD = PlaceTerms(
    mass_density_solar_pc3=tuning.SOLAR_NEIGHBOURHOOD_MASS_DENSITY_SOLAR_PC3,
    dispersion_kms=tuning.SOLAR_NEIGHBOURHOOD_DISPERSION_KMS)
"""PlaceTerms: The terms of a system drawn with no galaxy placement (a single system from the command line)."""


def from_profile(profile):
    """
    The terms of the place a `galaxy.profile.SectorProfile` describes: the mass density is the cell's expected
    systems over its volume times the mean system mass, the dispersion the root mean square of the populations'
    dispersions weighted by their shares (the heavier the share, the more of the stars that move that way).
    """
    systems_per_pc3 = profile.expected_systems / profile.volume_pc3 if profile.volume_pc3 > 0.0 else 0.0
    weighted = sum(p.share * p.dispersion_kms ** 2 for p in profile.populations.values())
    dispersion = math.sqrt(weighted) if weighted > 0.0 else tuning.SOLAR_NEIGHBOURHOOD_DISPERSION_KMS
    return PlaceTerms(mass_density_solar_pc3=systems_per_pc3 * tuning.MEAN_SYSTEM_MASS_SOLAR,
                      dispersion_kms=dispersion)


def close_companion_factor(fe_h):
    """
    `F(fe_h) / F(0)`, the factor on a star's chance of a companion inside 10 AU: `F` is the straight line through
    `tuning.CLOSE_BINARY_FRACTION_BY_FE_H`, held flat past its ends. 1 for a star with no [Fe/H].
    """
    if fe_h is None:
        return 1.0
    return _fraction_at(fe_h) / _fraction_at(0.0)


def _fraction_at(fe_h):
    points = sorted(tuning.CLOSE_BINARY_FRACTION_BY_FE_H.items())
    if fe_h <= points[0][0]:
        return points[0][1]
    for (low_x, low_y), (high_x, high_y) in zip(points, points[1:]):
        if fe_h <= high_x:
            return low_y + (high_y - low_y) * (fe_h - low_x) / (high_x - low_x)
    return points[-1][1]


def encounter_limit_au(terms, total_mass_solar, age_gy):
    """
    The semi-major axis at which half the pairs of this age survive passing stars (Bahcall, Hut and Tremaine
    1985, as Jiang and Tremaine 2010 restate it): `0.002 (M1 + M2) sigma / (G rho_2 t)`, with `rho_2` the mass
    density times the mean perturber mass. Infinite where nothing passes (a zero density) or the pair has no age.
    """
    rho2 = terms.mass_density_solar_pc3 * terms.perturber_mass_solar
    if rho2 <= 0.0 or age_gy <= 0.0:
        return math.inf
    rho2_si = rho2 * constants.SOLAR_MASS_TO_KG ** 2 / constants.PARSEC_M ** 3
    seconds = age_gy * 1.0e9 * constants.SECONDS_PER_YEAR
    metres = (tuning.ENCOUNTER_LIMIT_COEFFICIENT * total_mass_solar * constants.SOLAR_MASS_TO_KG
              * terms.dispersion_kms * 1000.0) / (constants.G * rho2_si * seconds)
    return metres / constants.AU_M


def tidal_limit_au(terms, total_mass_solar):
    """The Jacobi radius written with the local mass density, `(M / (3 rho))^(1/3)`, in AU: infinite at zero density."""
    if terms.mass_density_solar_pc3 <= 0.0:
        return math.inf
    parsecs = (total_mass_solar / (3.0 * terms.mass_density_solar_pc3)) ** (1.0 / 3.0)
    return parsecs * constants.PARSEC_M / constants.AU_M


def wide_axis_limit_au(terms, total_mass_solar, age_gy):
    """The widest axis a wide pair keeps (a sharp cut, GEN.103): the smaller of the encounter and tidal limits."""
    return min(encounter_limit_au(terms, total_mass_solar, age_gy), tidal_limit_au(terms, total_mass_solar))


def born_mass_solar(star):
    """A star's mass at birth in solar masses (its present mass for a star with no model initial mass)."""
    initial = getattr(star, "initial_mass_sol", None)
    return initial if initial is not None else star.mass / constants.SOLAR_MASS_TO_KG
