# planetgen/physics/black_hole_formation.py

"""
How many stellar-mass black holes a stellar population makes (GEN.563; docs/design/systems/black-hole-density.md,
docs/design/systems/black-hole-formation-sources.md). Pure functions; `galaxy.remnant_density` puts them in a place.

A black hole is the remnant of a massive star that fails to explode. Of the stars born with 8 to 150 solar masses
(`imf_massive_share`) a share `f` fail, and `f` falls as the star's metallicity rises: the failed-supernova share of a
metal-poor population is the larger. Only the range of `f` at two metallicities is known (`tuning.BLACK_HOLE_*`), so
a galaxy draws its two ends from its seed (`draw_failed_shares`), and `failed_share` runs linearly in the metal mass
fraction Z between them (the share saturates toward low Z, as Ebinger et al. 2020's grid does) and is held at the solar
end above solar. A population whose stars are older than the lives of its massive stars has made all it will
(`remnant_formed_share`); one still forming has made the share of them whose lives are over.
"""

import functools
import math

from planetgen import tuning
from planetgen.physics import imf, ms_lifetime, stellar_evolution


@functools.lru_cache(maxsize=1)
def imf_massive_share():
    """The share of the stars born that have `tuning.CORE_COLLAPSE_MASS_RANGE_SOL` solar masses (8 to 150): 0.713% of
    the generator's own IMF (the massive stars that can end as a neutron star or a black hole; 0.637% on the Kroupa function
    before GEN.927)."""
    low, high = tuning.CORE_COLLAPSE_MASS_RANGE_SOL
    return imf.number(low, high) / imf.number(imf.lowest_mass_sol(), imf.highest_mass_sol())


@functools.lru_cache(maxsize=1)
def _lifetime_table():
    """`(masses, remnant ages in Gy, IMF weights)` on a log grid over the massive stars: the age a progenitor of each
    mass has by the time its remnant exists (`remnant_age_fraction` times the main-sequence life), and its weight."""
    low, high = tuning.CORE_COLLAPSE_MASS_RANGE_SOL
    steps = tuning.BLACK_HOLE_FORMED_MASS_STEPS
    edges = [low * (high / low) ** (i / steps) for i in range(steps + 1)]
    masses = [math.sqrt(a * b) for a, b in zip(edges, edges[1:])]
    ages = [ms_lifetime.lifetime_gy(m) * stellar_evolution.remnant_age_fraction(m) for m in masses]
    weights = [imf.number(a, b) for a, b in zip(edges, edges[1:])]
    return tuple(ages), tuple(weights)


def _formed_at(age_gy):
    """The share of the massive stars (by number, IMF weighted) whose remnant exists in a population of this age."""
    ages, weights = _lifetime_table()
    return sum(w for a, w in zip(ages, weights) if a <= age_gy) / sum(weights)


@functools.lru_cache(maxsize=64)
def remnant_formed_share(age_low_gy, age_high_gy):
    """
    The share of a population's massive stars whose remnants have formed, for a population whose stars have ages
    uniform in `[age_low_gy, age_high_gy]`: the mean over the window of the IMF share of progenitors that lived
    shorter than the age. 1 for a window that starts after the longest life of a massive star (about 45 Myr to the
    end of the last giant phase of an 8 solar mass star, `stellar_evolution.remnant_age_fraction`), 0.79 for the young
    window (0 to 0.1 Gy).
    """
    if age_high_gy <= age_low_gy:
        return _formed_at(age_low_gy)
    samples = tuning.BLACK_HOLE_FORMED_AGE_STEPS
    return sum(_formed_at(age_low_gy + (age_high_gy - age_low_gy) * (i + 0.5) / samples)
               for i in range(samples)) / samples


def failed_share(fe_h, solar_share, metal_free_share):
    """
    The share of the massive stars of metallicity `fe_h` ([Fe/H], dex) that fail to explode and leave a black hole:
    `f_low + (f_sun - f_low) min(10^[Fe/H], 1)`, linear in the metal mass fraction between the metal-free share
    `f_low` and the solar share `f_sun`, and `f_sun` above solar. Works on numpy arrays too.
    """
    ratio = 10.0 ** fe_h
    try:
        ratio = ratio.clip(max=1.0)
    except AttributeError:
        ratio = min(ratio, 1.0)
    return metal_free_share + (solar_share - metal_free_share) * ratio


def share_cap(neutron_star_unbound_share):
    """
    The most `f` a galaxy can have so that the neutron stars and black holes formed fit among the massive stars born:
    `1 - NS_formed / IMF_SHARE`, `NS_formed` the neutron star rate per star (the bound count; the unbound share ejected
    on top is added) times `1 + unbound`. 21.4% to 32.7% for the unbound shares of 0.4 to 0.2 (12.1% to 24.7% on the Kroupa function before GEN.927; Boss, GEN.563 card 5).
    """
    formed = tuning.phenomenon_rate_per_star("neutron-star") * (1.0 + neutron_star_unbound_share)
    return max(0.0, 1.0 - formed / imf_massive_share())


def draw_failed_shares(uniform_solar, uniform_metal_free, cap):
    """
    A galaxy's two ends of `failed_share`, `(f_sun, f_low)`, from two uniform draws in 0 to 1 (the order matters: the
    solar one first). `f_sun` is uniform over 0 to the lesser of the widest verified range's top and `cap`; `f_low` is
    uniform over the metal-free range, cut at `cap`, and not below `f_sun` (no source gives fewer holes at low metallicity
    than at solar).
    """
    sun_low, sun_high = tuning.BLACK_HOLE_FAILED_SHARE_SOLAR_RANGE
    free_low, free_high = tuning.BLACK_HOLE_FAILED_SHARE_METAL_FREE_RANGE
    solar = sun_low + uniform_solar * (min(sun_high, cap) - sun_low)
    free = min(free_low + uniform_metal_free * (free_high - free_low), cap)
    return solar, max(free, solar)


def reference_failed_shares(neutron_star_unbound_share):
    """`draw_failed_shares` for the draws of a galaxy built with no seed (`tuning.REFERENCE_BLACK_HOLE_DRAWS`) and a neutron star
    unbound share: the shares of the reference structure (`galaxy.structure.REFERENCE_STRUCTURE`, whose unbound share is
    `galaxy.natal_kick`'s for the Milky Way's masses)."""
    solar_draw, metal_free_draw = tuning.REFERENCE_BLACK_HOLE_DRAWS
    return draw_failed_shares(solar_draw, metal_free_draw, share_cap(neutron_star_unbound_share))
