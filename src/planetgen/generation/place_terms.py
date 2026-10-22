# generation/place_terms.py

"""
What a place does to the stars drawn in it (GEN.103, docs/design/placement-by-population.md 3.1 and 3.3):
the close-companion odds factor that follows a star's own [Fe/H] and the widest axis a wide pair can keep.

Both are functions of numbers the sector profile already holds (the density, the populations' velocity
dispersions and the star's own [Fe/H] and age), so nothing here is a catalogue by place. `PlaceTerms` carries the
profile's part to the system draw on `SystemConfig.PLACE_TERMS`; a system drawn with no galaxy placement uses
`SOLAR_NEIGHBOURHOOD`.
"""

import functools
import math
from dataclasses import dataclass

import numpy as np

from planetgen import tuning
from planetgen.physics import constants, ms_lifetime, remnant_mass, remnant_progenitors


RELATIVE_VELOCITY_FACTOR = math.sqrt(2.0)
"""float: A pair's perturbers move past it with the *relative* velocity dispersion, root 2 times a population's own
one-dimensional dispersion for two stars drawn from the same distribution (Jiang and Tremaine 2010, section 5: "the
relative velocity dispersion to be 2^(1/2) times this", 29 km/s becoming 40 km/s). The population dispersions of
`tuning.STELLAR_POPULATION_DISPERSION_KMS` are one-dimensional, so `encounter_limit_au` multiplies by this (GEN.447)."""

_COHORT_AGES = 24
"""int: Ages of the midpoint rule that averages a population's cohorts, which formed at a steady rate over its range."""

_REMNANT_GRID_POINTS = 400
"""int: Points of the log grid the white dwarf masses are integrated on (the initial-final relation is clipped)."""


def _imf_segments():
    """`((low, high, alpha, k), ...)`: the generator's initial mass function (`tuning.IMF_BREAKS_SOL`, `IMF_SLOPES`),
    each segment `k M^-alpha`, joined continuously."""
    breaks, slopes = tuning.IMF_BREAKS_SOL, tuning.IMF_SLOPES
    norms = [1.0]
    for i in range(1, len(slopes)):
        norms.append(norms[-1] * breaks[i] ** (slopes[i] - slopes[i - 1]))
    return tuple((breaks[i], breaks[i + 1], slopes[i], norms[i]) for i in range(len(slopes)))


def _imf_integral(power, low, high):
    """`integral of M^power xi(M) dM` from `low` to `high` (initial masses, solar masses; `xi` the IMF)."""
    total = 0.0
    for seg_low, seg_high, alpha, k in _imf_segments():
        a, b = max(low, seg_low), min(high, seg_high)
        if a >= b:
            continue
        exponent = power + 1.0 - alpha
        total += k * (math.log(b / a) if exponent == 0.0 else (b ** exponent - a ** exponent) / exponent)
    return total


def _present_day_moments(age_gy, remnants=True):
    """
    `(sum of M, sum of M squared)` over everything a cohort of one age holds now, per unit of the IMF's number: the
    stars still on the main sequence (initial mass up to the turn-off, where the main-sequence lifetime reaches the
    age) and the remnants of the rest, by the progenitor ranges and windows the generator draws remnants from
    (`tuning.WHITE_DWARF_PROGENITOR_SOL` and the others) and their mass rules: white dwarfs by the initial-final
    relation, neutron stars and black holes uniform over their mass ranges. Neutron stars and black holes take the
    progenitors of their windows above the white dwarf range's top (`remnant_progenitors`).
    """
    top = tuning.IMF_BREAKS_SOL[-1]
    turn_off = min(top, ms_lifetime.mass_for_lifetime_gy(age_gy))
    low_edge = tuning.IMF_BREAKS_SOL[0]
    mass = _imf_integral(1.0, low_edge, turn_off)
    square = _imf_integral(2.0, low_edge, turn_off)
    if not remnants:
        return mass, square

    wd_low, wd_high = tuning.WHITE_DWARF_PROGENITOR_SOL
    wd_low = max(wd_low, turn_off)
    if wd_low < wd_high:
        grid = np.geomspace(wd_low, wd_high, _REMNANT_GRID_POINTS)
        middle = 0.5 * (grid[:-1] + grid[1:])
        weight = np.array([_imf_integral(0.0, a, b) for a, b in zip(grid[:-1], grid[1:])])
        remnant = remnant_mass.white_dwarf_mass_sol(middle)
        mass += float(np.sum(weight * remnant))
        square += float(np.sum(weight * remnant * remnant))

    for kind in remnant_progenitors.KINDS:
        number = sum(_imf_integral(0.0, max(low, wd_high, turn_off), high)
                     for low, high in remnant_progenitors.kind_windows(kind) if max(low, wd_high, turn_off) < high)
        first, second = remnant_mass.uniform_moments(kind)
        mass += number * first
        square += number * second
    return mass, square


@functools.lru_cache(maxsize=256)
def _cohort_moments(age_low_gy, age_high_gy, remnants=True):
    """`_present_day_moments` averaged over a population whose stars formed at a steady rate between two ages."""
    low, high = sorted((max(age_low_gy, 1.0e-4), max(age_high_gy, 1.0e-4)))
    ages = low + (high - low) * (np.arange(_COHORT_AGES) + 0.5) / _COHORT_AGES
    moments = [_present_day_moments(float(age), remnants) for age in ages]
    return (sum(m for m, _s in moments) / _COHORT_AGES, sum(s for _m, s in moments) / _COHORT_AGES)


def mass_weighted_perturber_mass_solar(populations, remnants=True):
    """
    The perturber mass of Jiang and Tremaine 2010's `rho_2` (their eq. 47): the second moment of the mass function over
    its first, `<M^2> / <M>`, so that `rho_2 = rho * <M^2> / <M>` for a mass density `rho` (GEN.447). The mass function is
    the generator's own: its initial mass function, aged to each population's age range, with the remnants of what has
    left the main sequence.

    Args:
        populations (iterable): `(weight, (youngest_age_gy, oldest_age_gy))` pairs, the weight each population's share
            of the stars formed.
        remnants (bool): Whether the white dwarfs, neutron stars and black holes count (they do; `False` is the stars
            alone, the comparison with Jiang and Tremaine's observed mass function).

    Returns:
        float: Solar masses.
    """
    first = second = 0.0
    for weight, (age_low, age_high) in populations:
        mass, square = _cohort_moments(round(float(age_low), 6), round(float(age_high), 6), remnants)
        first += weight * mass
        second += weight * square
    return second / first if first > 0.0 else tuning.MEAN_SYSTEM_MASS_SOLAR


@dataclass(frozen=True)
class PlaceTerms:
    """
    The part of a sector profile the pair draws read.

    Attributes:
        mass_density_solar_pc3 (float): Stellar mass per cubic parsec at the place.
        dispersion_kms (float): The velocity dispersion of the field there, km/s.
        perturber_mass_solar (float): The mass-weighted mean mass, `<M^2> / <M>`, of the stars and remnants that pass a
            pair (`mass_weighted_perturber_mass_solar`); the default is the whole disk's star formation history.
    """
    mass_density_solar_pc3: float
    dispersion_kms: float
    perturber_mass_solar: float = mass_weighted_perturber_mass_solar([(1.0, tuning.STAR_FORMATION_AGE_RANGE_GY)])


SOLAR_NEIGHBOURHOOD = PlaceTerms(
    mass_density_solar_pc3=tuning.SOLAR_NEIGHBOURHOOD_MASS_DENSITY_SOLAR_PC3,
    dispersion_kms=tuning.SOLAR_NEIGHBOURHOOD_DISPERSION_KMS)
"""PlaceTerms: The terms of a system drawn with no galaxy placement (a single system from the command line)."""


def from_profile(profile):
    """
    The terms of the place a `galaxy.profile.SectorProfile` describes: the mass density is the cell's expected
    systems over its volume times the mean system mass, the dispersion the root mean square of the populations'
    dispersions weighted by their shares (the heavier the share, the more of the stars that move that way), the
    perturber mass the mass-weighted mean of the mass function of the populations the cell holds (GEN.447).
    """
    systems_per_pc3 = profile.expected_systems / profile.volume_pc3 if profile.volume_pc3 > 0.0 else 0.0
    weighted = sum(p.share * p.dispersion_kms ** 2 for p in profile.populations.values())
    dispersion = math.sqrt(weighted) if weighted > 0.0 else tuning.SOLAR_NEIGHBOURHOOD_DISPERSION_KMS
    perturber = mass_weighted_perturber_mass_solar(
        (p.share, p.age_gy) for p in profile.populations.values() if p.share > 0.0)
    return PlaceTerms(mass_density_solar_pc3=systems_per_pc3 * tuning.MEAN_SYSTEM_MASS_SOLAR,
                      dispersion_kms=dispersion, perturber_mass_solar=perturber)


def close_companion_factor(fe_h):
    """
    `F(fe_h) / F(0)`, the factor on the odds of a star's companion inside 10 AU (`scaled_close_chance`): `F` is the
    straight line through `tuning.CLOSE_BINARY_FRACTION_BY_FE_H`, held flat past its ends. 1 for a star with no [Fe/H].
    """
    if fe_h is None:
        return 1.0
    return _fraction_at(fe_h) / _fraction_at(0.0)


def scaled_close_chance(chance, fe_h):
    """
    A close-companion chance moved by the star's [Fe/H] (Boss, 2026-10-16, 7:19 PM EDT): the factor
    `close_companion_factor` multiplies the *odds* `p / (1 - p)`, not `p`, for every primary mass. The result never passes 1,
    and a star already likely to have a companion rises less in absolute terms: a G dwarf (0.22) goes 0.43, 0.36, 0.25
    and 0.12 at [Fe/H] -3.0, -1.0, -0.2 and +0.5, inside Moe, Kratter and Badenes 2019's ranges, and an O star (0.45)
    goes to 0.58 at -0.7, inside the "no change within about 30%" the O and early-B surveys find.
    """
    if chance <= 0.0 or chance >= 1.0:
        return min(max(chance, 0.0), 1.0)
    odds = chance / (1.0 - chance) * close_companion_factor(fe_h)
    return odds / (1.0 + odds)


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
    density times the mass-weighted perturber mass (`<M^2> / <M>`) and `sigma` the dispersion of the relative velocity,
    `RELATIVE_VELOCITY_FACTOR` times the field's one-dimensional dispersion (their eq. 1 and section 5; GEN.447).
    Infinite where nothing passes (a zero density) or the pair has no age.
    """
    rho2 = terms.mass_density_solar_pc3 * terms.perturber_mass_solar
    if rho2 <= 0.0 or age_gy <= 0.0:
        return math.inf
    rho2_si = rho2 * constants.SOLAR_MASS_TO_KG ** 2 / constants.PARSEC_M ** 3
    seconds = age_gy * 1.0e9 * constants.SECONDS_PER_YEAR
    metres = (tuning.ENCOUNTER_LIMIT_COEFFICIENT * total_mass_solar * constants.SOLAR_MASS_TO_KG
              * RELATIVE_VELOCITY_FACTOR * terms.dispersion_kms * 1000.0) / (constants.G * rho2_si * seconds)
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
