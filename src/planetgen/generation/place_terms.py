# generation/place_terms.py

"""
What a place does to the stars drawn in it (GEN.103, docs/design/generation/placement-by-population.md 3.1 and 3.3):
the widest axis a wide pair can keep and the widest orbit a planet keeps (GEN.325). The close-companion odds factor that
follows a star's own [Fe/H] is `generation.close_binary_trend`.

All are functions of numbers the sector profile already holds (the density, the populations' velocity
dispersions and the star's own [Fe/H] and age), so nothing here is a catalogue by place. `PlaceTerms` carries the
profile's part to the system draw on `SystemConfig.PLACE_TERMS`; a system drawn with no galaxy placement uses
`SOLAR_NEIGHBOURHOOD`.
"""

import math
from dataclasses import dataclass, field

from planetgen import tuning
from planetgen.galaxy import system_mass as galaxy_system_mass
from planetgen.generation import system_mass
from planetgen.physics import constants, planet_survival
from planetgen.physics import unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")

_yr_to_s = units.converter("yr", "s")
_gyr_to_yr = units.converter("Gyr", "yr")
_m_to_au = units.converter("m", "au")
_pc_to_m = units.converter("pc", "m")


RELATIVE_VELOCITY_FACTOR = math.sqrt(2.0)
"""float: A pair's perturbers move past it with the *relative* velocity dispersion, root 2 times a population's own
one-dimensional dispersion for two stars drawn from the same distribution ([R2166], section 5: "the
relative velocity dispersion to be 2^(1/2) times this", 29 km/s becoming 40 km/s). The population dispersions of
`tuning.STELLAR_POPULATION_DISPERSION_KMS` are one-dimensional, so `encounter_limit_au` multiplies by this (GEN.447)."""


def mass_weighted_perturber_mass_solar(populations, remnants=True):
    """
    The perturber mass of [R2166]'s `rho_2` (their eq. 47): the second moment of the mass function over
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
        mass, square = system_mass.cohort_moments(round(float(age_low), 6), round(float(age_high), 6), remnants)
        first += weight * mass
        second += weight * square
    return second / first if first > 0.0 else system_mass.window_mean_system_mass_solar(*tuning.STAR_FORMATION_AGE_RANGE_GY)


@dataclass(frozen=True)
class PlaceTerms:
    """
    The part of a sector profile the pair draws read.

    Attributes:
        mass_density_solar_pc3 (float): Stellar mass per cubic parsec at the place.
        dispersion_kms (float): The velocity dispersion of the field there, km/s.
        perturber_mass_solar (float): The mass-weighted mean mass, `<M^2> / <M>`, of the stars and remnants that pass a
            pair (`mass_weighted_perturber_mass_solar`); the default is the whole disk's star formation history.
        mean_system_mass_solar (float): The mean mass of a system there (`galaxy.system_mass.mean_over_solar`), which
            turns the mass density into the systems per volume that strip planets (`systems_per_pc3`); the default is
            the whole disk's star formation history.
    """
    mass_density_solar_pc3: float
    dispersion_kms: float
    perturber_mass_solar: float = field(default_factory=lambda: mass_weighted_perturber_mass_solar(
        [(1.0, tuning.STAR_FORMATION_AGE_RANGE_GY)]))
    mean_system_mass_solar: float = field(default_factory=lambda: system_mass.window_mean_system_mass_solar(
        *tuning.STAR_FORMATION_AGE_RANGE_GY))

    @property
    def systems_per_pc3(self):
        """Systems per cubic parsec: the mass density over the mean system mass."""
        return self.mass_density_solar_pc3 / self.mean_system_mass_solar


def __getattr__(name):
    """
    `SOLAR_NEIGHBOURHOOD`, the terms of a system drawn with no galaxy placement (a single system from the command line),
    built on first use: its defaults read the post-main-sequence tracks, which setup downloads, so importing this module
    (and the setup command itself) must not need them.
    """
    if name == "SOLAR_NEIGHBOURHOOD":
        terms = PlaceTerms(mass_density_solar_pc3=tuning.SOLAR_NEIGHBOURHOOD_MASS_DENSITY_SOLAR_PC3,
                           dispersion_kms=tuning.SOLAR_NEIGHBOURHOOD_DISPERSION_KMS)
        globals()[name] = terms
        return terms
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def from_profile(profile):
    """
    The terms of the place a `galaxy.profile.SectorProfile` describes: the mass density is the cell's expected
    systems over its volume times the mean system mass, the dispersion the root mean square of the populations'
    dispersions (the profile's three-axis totals, so divided by `constants.ISOTROPIC_TOTAL_OVER_ONE_AXIS` to the one-axis
    value the limits take) weighted by their shares (the heavier the share, the more of the stars that move that way), the
    perturber mass the mass-weighted mean of the mass function of the populations the cell holds (GEN.447).
    """
    systems_per_pc3 = profile.expected_systems / profile.volume_pc3 if profile.volume_pc3 > 0.0 else 0.0
    weighted = sum(p.share * p.dispersion_kms ** 2 for p in profile.populations.values())
    dispersion = (math.sqrt(weighted) / constants.ISOTROPIC_TOTAL_OVER_ONE_AXIS if weighted > 0.0
                  else tuning.SOLAR_NEIGHBOURHOOD_DISPERSION_KMS)
    perturber = mass_weighted_perturber_mass_solar(
        (p.share, p.age_gy) for p in profile.populations.values() if p.share > 0.0)
    mean_mass = galaxy_system_mass.mean_over_solar({n: p.share for n, p in profile.populations.items()}, profile.shape)
    return PlaceTerms(mass_density_solar_pc3=systems_per_pc3 * mean_mass, dispersion_kms=dispersion,
                      perturber_mass_solar=perturber, mean_system_mass_solar=mean_mass)


def encounter_limit_au(terms, total_mass_solar, age_gy):
    """
    The semi-major axis at which half the pairs of this age survive passing stars (Bahcall, Hut and Tremaine
    1985, as [R2166] restate it): `0.002 (M1 + M2) sigma / (G rho_2 t)`, with `rho_2` the mass
    density times the mass-weighted perturber mass (`<M^2> / <M>`) and `sigma` the dispersion of the relative velocity,
    `RELATIVE_VELOCITY_FACTOR` times the field's one-dimensional dispersion (their eq. 1 and section 5; GEN.447).
    Infinite where nothing passes (a zero density) or the pair has no age.
    """
    rho2 = terms.mass_density_solar_pc3 * terms.perturber_mass_solar
    if rho2 <= 0.0 or age_gy <= 0.0:
        return math.inf
    rho2_si = rho2 * constants.SOLAR_MASS_TO_KG ** 2 / constants.PARSEC_M ** 3
    seconds = _yr_to_s(_gyr_to_yr(age_gy))
    metres = (tuning.ENCOUNTER_LIMIT_COEFFICIENT * total_mass_solar * constants.SOLAR_MASS_TO_KG
              * RELATIVE_VELOCITY_FACTOR * terms.dispersion_kms * 1000.0) / (constants.G * rho2_si * seconds)
    return _m_to_au(metres)


def planet_limit_au(terms, age_gy=None):
    """
    The largest semi-major axis at which a planet of a star of `age_gy` outlives the passing systems of the place
    (`physics.planet_survival.planet_survival_limit_au` on the place's systems per volume and dispersion). A star
    younger than `tuning.PLANET_SURVIVAL_AGE_GY` has been exposed to its own age only; one with no age (or an older
    one) to that age. Infinite where nothing passes (GEN.325).
    """
    age = tuning.PLANET_SURVIVAL_AGE_GY if age_gy is None else min(age_gy, tuning.PLANET_SURVIVAL_AGE_GY)
    return planet_survival.planet_survival_limit_au(terms.systems_per_pc3, terms.dispersion_kms, age)


def tidal_limit_au(terms, total_mass_solar):
    """The Jacobi radius written with the local mass density, `(M / (3 rho))^(1/3)`, in AU: infinite at zero density."""
    if terms.mass_density_solar_pc3 <= 0.0:
        return math.inf
    parsecs = (total_mass_solar / (3.0 * terms.mass_density_solar_pc3)) ** (1.0 / 3.0)
    return _m_to_au(_pc_to_m(parsecs))


def wide_axis_limit_au(terms, total_mass_solar, age_gy):
    """The widest axis a wide pair keeps (a sharp cut, GEN.103): the smaller of the encounter and tidal limits."""
    return min(encounter_limit_au(terms, total_mass_solar, age_gy), tidal_limit_au(terms, total_mass_solar))


def born_mass_solar(star):
    """A star's mass at birth in solar masses (its present mass for a star with no model initial mass)."""
    initial = getattr(star, "initial_mass_sol", None)
    return initial if initial is not None else _kg_to_msun(star.mass)
