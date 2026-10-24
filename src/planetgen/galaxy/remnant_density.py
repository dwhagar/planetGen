# planetgen/galaxy/remnant_density.py

"""
How many stellar-mass black holes and neutron stars a place holds (GEN.563 and GEN.554; docs/design/systems/black-hole-density.md,
docs/design/systems/neutron-star-density.md).

A black hole is what a massive star leaves when it fails to explode, and the share that fails falls as the star's metals
rise (`physics.black_hole_formation`). So the black holes per star at a place are

    r_BH = IMF_SHARE x RETAIN x sum over populations p [ w_p x C_p x F_p ]  /  (1 - IMF_SHARE)

with `w_p` the population's share of the place's stars, `C_p` the share of its massive stars whose remnants have formed
(1 for every population older than 100 Myr, 0.79 for the young one), `F_p` the mean over its [Fe/H] distribution at the
place of the galaxy's failed-supernova share (`Structure.black_hole_solar_share` and `black_hole_metal_free_share`), and
`RETAIN` the share that stays bound. The first factor is `base_rate_per_star`, the rest `place_factor`; the scatter
weights its bins by the sum alone (`weights_array`), the sector fill and the profile multiply the base rate by it.

The populations' [Fe/H] distributions are the generator's own (`physics.metallicity`): the disc's gradient, the bulge's
and the halo's means, and for a cluster or a satellite the generic peaks of `tuning` (a percent of the systems; the
sector profile uses the exact peaks of the components that reach the cell, `place_factor_from_peaks`).

Every star of 8 to 150 solar masses leaves one or the other (GEN.554), so the neutron stars are what the black holes leave
and the same sum with `1 - F_p` for `F_p` gives them: `r_NS = IMF_SHARE x (1 - u) x sum [ w_p x C_p x (1 - F_p) ] / (1 - IMF_SHARE)`,
`u` the share of the neutron stars formed that the galaxy ejects (`Structure.neutron_star_unbound_share`). Each function
takes the `kind` ("black-hole", the default, or "neutron-star").

The intermediate-mass black holes are not remnants of failed supernovae and keep a rate of their own.
"""

import functools
import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import density as galaxyDensity
from planetgen.physics import black_hole_formation, metallicity
from planetgen.physics.units import pc_to_ly

BLACK_HOLE = "black-hole"
NEUTRON_STAR = "neutron-star"
KINDS = (BLACK_HOLE, NEUTRON_STAR)
"""tuple: The kinds of remnant a massive star leaves, whose rates follow the metals (`kind` of every function here)."""

_GENERIC_PEAKS = {
    "cluster": tuning.METALLICITY_CLUSTER_PEAKS,
    "satellite": tuning.METALLICITY_SATELLITE_PEAKS,
    "satellite_young": tuning.METALLICITY_SATELLITE_PEAKS,
    "satellite_intermediate": tuning.METALLICITY_SATELLITE_PEAKS,
}
"""dict: The populations that carry their components' [Fe/H], and the peaks a place with no component of its own takes."""


def retained_share():
    """The share of the black holes formed that stay bound to the galaxy (the rest are ejected, `ejecta_field`)."""
    return 1.0 - tuning.BLACK_HOLE_UNBOUND_FRACTION


@functools.lru_cache(maxsize=2)
def base_rate_per_star(kind=BLACK_HOLE):
    """
    The bound remnants per star of a population in which every massive star has made its remnant and all of them are
    this `kind`: the IMF's 8 to 150 solar mass share (times `retained_share()` for the black holes; the neutron stars'
    bound share is in their place factor, which the galaxy's own unbound share sets). The scatter's rate for the class.
    """
    _check_kind(kind)
    share = black_hole_formation.imf_massive_share()
    return share * retained_share() if kind == BLACK_HOLE else share


def follows_metals(kind, subtype):
    """Whether a scatter class's rate follows the metals: the stellar-mass black holes and the neutron stars that are not
    magnetars (the magnetars follow the young stars, the intermediate-mass black holes have a rate of their own)."""
    return (kind == BLACK_HOLE and subtype == "stellar") or (kind == NEUTRON_STAR and subtype is None)


def _check_kind(kind):
    if kind not in KINDS:
        raise ValueError(f"{kind!r} is not a remnant kind whose rate follows the metals ({', '.join(KINDS)})")


@functools.lru_cache(maxsize=1)
def _per_generator_star():
    """Divides a population's share by the stars the generator keeps: its star counts leave out the stars that have
    already collapsed (`1 - IMF_SHARE`)."""
    return 1.0 / (1.0 - black_hole_formation.imf_massive_share())


def _formed(name):
    """`C_p`, the share of population `name`'s massive stars whose remnants exist, over the generator's stars."""
    low, high = tuning.STELLAR_POPULATION_AGE_RANGES_GY[name]
    return black_hole_formation.remnant_formed_share(float(low), float(high)) * _per_generator_star()


@functools.lru_cache(maxsize=256)
def _failed_mean_table(solar, free, sigma):
    """`(means, F)`: `F_p` on a grid of mean [Fe/H] from four scatters under the metallicity range to four over it. The
    mean over a normal peak is a smooth function of its mean alone, so the scatter's many points read it from here."""
    low, high = tuning.METALLICITY_RANGE_FE_H
    means = np.linspace(low - 4.0 * sigma, high + 4.0 * sigma, tuning.BLACK_HOLE_FAILED_MEAN_TABLE_STEPS)
    values = metallicity.expected_factor_array(lambda fe_h: black_hole_formation.failed_share(fe_h, solar, free), means,
                                               sigma)
    return means, values


def _failed_mean(structure, means, sigma):
    """`F_p` for a normal [Fe/H] distribution of mean `means` (array) and scatter `sigma`."""
    table_means, table = _failed_mean_table(float(structure.black_hole_solar_share),
                                            float(structure.black_hole_metal_free_share), float(sigma))
    return np.interp(means, table_means, table)


def _failed_mean_of_peaks(structure, peaks):
    """`F_p` for peaks `((weight, mean, sigma), ...)`."""
    total = sum(weight for weight, _mean, _sigma in peaks)
    return sum(weight * _failed_mean(structure, np.asarray(mean, dtype=float), sigma)
               for weight, mean, sigma in peaks) / total


def _own_share(kind, structure, failed, bound):
    """A population's mean share of the massive stars that leave this `kind`, from its mean failed share `failed`: the
    failed share for a black hole; what is left, less the unbound share when `bound`, for a neutron star."""
    if kind == BLACK_HOLE:
        return failed
    left = 1.0 - failed
    return left * (1.0 - structure.neutron_star_unbound_share) if bound else left


def population_factors(structure, distance_ly, kind=BLACK_HOLE, bound=True):
    """
    `{population: C_p F_p}` at arrays of distances from the centre, in light years (or one): the factor each population's
    density is weighed by (`C_p (1 - F_p)`, less the unbound share, for the neutron stars: `bound` False leaves that out,
    for the count formed). The generic peaks stand for the cluster and satellite populations.
    """
    _check_kind(kind)
    distance = np.asarray(distance_ly, dtype=float)
    factors = {}
    by_disc = {}
    for name in tuning.STELLAR_POPULATION_AGE_RANGES_GY:
        if name in _GENERIC_PEAKS:
            mean = _failed_mean_of_peaks(structure, _GENERIC_PEAKS[name]) + np.zeros_like(distance)
        else:
            key = metallicity.open_cluster_disk_population(name) if name in tuning.OPEN_CLUSTER_POPULATIONS else name
            if key not in by_disc:
                means, sigma = metallicity.mean_fe_h_array(key, distance)
                by_disc[key] = _failed_mean(structure, means, sigma)
            mean = by_disc[key]
        factors[name] = _formed(name) * _own_share(kind, structure, mean, bound)
    return factors


def weights_array(x, y, z, shape, kind=BLACK_HOLE):
    """
    The `kind`'s weight per unit of `base_rate_per_star` at arrays of points, parsecs: the populations'
    densities (`galaxy.density.population_densities_array`, the model's relative density) each times `C_p F_p` at the
    point. What the scatter weighs a longitude bin by in place of the total density.
    """
    mix = galaxyDensity.population_densities_array(x, y, z, shape)
    distance = pc_to_ly(np.sqrt(np.asarray(x, dtype=float) ** 2 + np.asarray(y, dtype=float) ** 2 + np.asarray(z) ** 2))
    factors = population_factors(shape.structure, distance, kind)
    return sum(np.maximum(density, 0.0) * factors[name] for name, density in mix.items())


def place_factor(position_pc, shape, kind=BLACK_HOLE):
    """
    The place's `kind` per star over `base_rate_per_star`, `sum w_p C_p F_p`, at one point (the one-point form of
    `weights_array` over the total density, with the clusters' and satellites' densities in the mix). 0 where nothing is.
    """
    mix = {name: max(value, 0.0) for name, value in galaxyDensity.population_densities(position_pc, shape).items()}
    total = sum(mix.values())
    if total <= 0.0:
        return 0.0
    distance = pc_to_ly(math.sqrt(sum(c * c for c in position_pc)))
    factors = population_factors(shape.structure, distance, kind)
    return sum(value * float(factors[name]) for name, value in mix.items()) / total


def place_factor_from_peaks(populations, structure, kind=BLACK_HOLE):
    """
    `place_factor` for the sector profile's `populations` (name -> `PopulationProfile`): the shares and the [Fe/H] peaks of
    the place itself, the exact ones where a cluster or a satellite reaches it.
    """
    total = sum(pop.share for pop in populations.values())
    if total <= 0.0:
        return 0.0
    return sum(pop.share * _formed(name) * _own_share(kind, structure, float(_failed_mean_of_peaks(structure, pop.fe_h)), True)
               for name, pop in populations.items() if pop.share > 0.0) / total


def component_factor(shares, structure=None, kind=BLACK_HOLE):
    """`place_factor` for a component whose stars have the population `shares` (`{population: share}`): the cluster or
    satellite's own mix, with the generic peaks (the middle draw's structure without one)."""
    from planetgen.galaxy import structure as galaxyStructure
    factors = population_factors(galaxyStructure.REFERENCE_STRUCTURE if structure is None else structure, 0.0, kind)
    return sum(share * float(factors[name]) for name, share in shares.items())


def bound_factor(structure, kind=BLACK_HOLE):
    """An upper bound on any population's `C_p F_p` anywhere: the metal-free share (the failed share never rises with the
    metals, and the metal-free share is not below the solar one) over the generator's stars for the black holes; for the
    neutron stars what the solar share leaves, less the unbound share. What the object-first scatter bounds a weight by."""
    _check_kind(kind)
    if kind == BLACK_HOLE:
        return structure.black_hole_metal_free_share * _per_generator_star()
    return (1.0 - structure.black_hole_solar_share) * (1.0 - structure.neutron_star_unbound_share) * _per_generator_star()


def standalone_factor(structure=None, kind=BLACK_HOLE):
    """The factor of a sector with no galaxy around it, taken at solar metallicity with every remnant formed: the
    failed share at [Fe/H] = 0 (the galaxy's solar share; the middle draw without a structure), or what it leaves."""
    from planetgen.galaxy import structure as galaxyStructure
    _check_kind(kind)
    structure = galaxyStructure.REFERENCE_STRUCTURE if structure is None else structure
    failed = float(black_hole_formation.failed_share(0.0, structure.black_hole_solar_share,
                                                     structure.black_hole_metal_free_share))
    return _own_share(kind, structure, failed, True) * _per_generator_star()


def rate_per_star(kind_factors=None, kind=BLACK_HOLE):
    """The `kind` per star a sector with these `kind_factors` (the profile's, none for a standalone sector) rolls."""
    factor = standalone_factor(kind=kind) if not kind_factors or kind not in kind_factors else kind_factors[kind]
    return base_rate_per_star(kind) * factor


@functools.lru_cache(maxsize=128)
def mean_factor(shape, kind=BLACK_HOLE, bound=True):
    """
    The galaxy's mean `sum w_p C_p F_p` per star (the neutron stars' with `bound` False: the count formed):
    `place_factor` averaged over the stars, on a grid of log-spaced radii
    (from 1 pc to a thousand disk scale lengths), equal-area polar angles and a few azimuths, weighted by the density and
    the volume. The halo and the bulge are in it; the clusters' and satellites' stars (a percent) are not.
    """
    radial = tuning.BLACK_HOLE_GALAXY_MEAN_RADIAL_STEPS
    polar = tuning.BLACK_HOLE_GALAXY_MEAN_POLAR_STEPS
    azimuthal = tuning.BLACK_HOLE_GALAXY_MEAN_AZIMUTH_STEPS
    top = 1000.0 * shape.disk_scale_length_pc
    log_r = (np.arange(radial) + 0.5) / radial * math.log(top)
    r = np.exp(log_r)
    cos_theta = -1.0 + 2.0 * (np.arange(polar) + 0.5) / polar
    phi = 2.0 * math.pi * (np.arange(azimuthal) + 0.5) / azimuthal
    rr, cc, pp = np.meshgrid(r, cos_theta, phi, indexing="ij")
    sin = np.sqrt(1.0 - cc * cc)
    x, y, z = rr * sin * np.cos(pp), rr * sin * np.sin(pp), rr * cc
    volume = rr ** 3                    # r^3 d(ln r) d(cos theta) d(phi): the common steps drop out of the ratio
    mix = galaxyDensity.population_densities_array(x, y, z, shape)
    factors = population_factors(shape.structure, pc_to_ly(rr), kind, bound)
    total = sum(np.maximum(density, 0.0) for density in mix.values())
    weighted = sum(np.maximum(mix[name], 0.0) * factors[name] for name in mix)
    return float((volume * weighted).sum() / (volume * total).sum())


def formed_per_star(shape, kind=BLACK_HOLE):
    """The `kind` the galaxy has formed per star, the bound ones and the ejected (`ejecta_field`)."""
    _check_kind(kind)
    if kind == BLACK_HOLE:
        return base_rate_per_star() / retained_share() * mean_factor(shape)
    return base_rate_per_star(kind) * mean_factor(shape, kind, False)
