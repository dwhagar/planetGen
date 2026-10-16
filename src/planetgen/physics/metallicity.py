# planetgen/physics/metallicity.py

"""
Metallicity (GEN.158): one [Fe/H] per star, in dex against the Sun.

`sample_fe_h` draws it from the star's stellar population and its distance
from the galactic centre: the disc populations follow a mild radial gradient
that is zero at the Sun's distance, the bulge and the halo have their own mean
and scatter, and a cluster carries its own value (the caller passes it). The
first use is the planet-occurrence factor: a protoplanetary disk's solids scale
as `solids_factor(fe_h)` and giant-planet occurrence as
`giant_planet_factor(fe_h)`. The constants are in `tuning` ("Metallicity").

The planet cull (GEN.236) is the occurrence of planets against the disk's: `planet_occurrence_factor` for planets of
any kind and `giant_occurrence_factor` for giants, both 1 at and above `floor_fe_h()`;
`expected_factor` averages one over a population's [Fe/H] distribution for a place's profile.
"""

import math

from planetgen import tuning
from planetgen.physics import constants
from planetgen.util import draw

_KPC_IN_LY = constants.KILOPARSEC_M / constants.LY_TO_M
_SUN_DISTANCE_KPC = constants.GALACTIC_CENTER_DISTANCE_LY / _KPC_IN_LY


def disc_gradient_fe_h(galactic_center_dist_ly):
    """The mean [Fe/H] of the disc at this distance from the centre: zero at the Sun's distance, rising
    inward and falling outward by `METALLICITY_DISC_FALL_DEX_PER_KPC`."""
    return -tuning.METALLICITY_DISC_FALL_DEX_PER_KPC * (galactic_center_dist_ly / _KPC_IN_LY - _SUN_DISTANCE_KPC) + 0.0


def halo_mean_fe_h(galactic_center_dist_ly):
    """The mean [Fe/H] of the halo at this distance from the centre (-1.6 inside 25 kpc, -2.2 outside)."""
    kpc = galactic_center_dist_ly / _KPC_IN_LY
    halo = tuning.METALLICITY_HALO
    return halo["inner"] if kpc < halo["break_kpc"] else halo["outer"]


def mean_fe_h(population, galactic_center_dist_ly):
    """The mean and scatter `(mean, sigma)` of [Fe/H] for a population at this distance from the centre."""
    if population == "halo":
        return halo_mean_fe_h(galactic_center_dist_ly), tuning.METALLICITY_HALO["sigma"]
    if population == "cluster":
        raise ValueError("a cluster carries its own [Fe/H]; pass it instead of sampling one")
    try:
        row = tuning.METALLICITY_BY_POPULATION[population]
    except KeyError:
        raise ValueError(f"unknown stellar population {population!r}; expected one of "
                         f"{sorted(tuning.METALLICITY_BY_POPULATION) + ['halo']}") from None
    gradient = disc_gradient_fe_h(galactic_center_dist_ly) if row["gradient"] else 0.0
    return row["offset"] + gradient, row["sigma"]


def sample_fe_h(population, galactic_center_dist_ly, rng=draw):
    """A star's [Fe/H] (dex): normal about the population's mean at this distance, clamped to
    `METALLICITY_RANGE_FE_H`. `population` None means the whole disc (`intermediate`)."""
    mean, sigma = mean_fe_h("intermediate" if population is None else population, galactic_center_dist_ly)
    low, high = tuning.METALLICITY_RANGE_FE_H
    return min(max(rng.gauss(mean, sigma), low), high)


def solids_factor(fe_h):
    """How much a disk's solid surface density scales against the Sun's: `10^[Fe/H]`."""
    return 10.0 ** fe_h


def giant_planet_factor(fe_h):
    """How much giant-planet occurrence scales against the Sun's: `10^(2 [Fe/H])`, capped at 1."""
    return min(1.0, 10.0 ** (tuning.GIANT_PLANET_METALLICITY_EXPONENT * fe_h))


def floor_fe_h():
    """The [Fe/H] at and above which planet occurrence is the disk's: `-tuning.PLANET_OCCURRENCE_FLOOR_DEX_BELOW_SOLAR`."""
    return -tuning.PLANET_OCCURRENCE_FLOOR_DEX_BELOW_SOLAR


def planet_occurrence_factor(fe_h):
    """Planet occurrence at `fe_h` against the disk's, which the generator's counts are those of: 1 at and above
    `floor_fe_h()`, `10^(exponent (fe_h - floor))` below it."""
    below = fe_h - floor_fe_h()
    return 1.0 if below >= 0.0 else 10.0 ** (tuning.PLANET_OCCURRENCE_METALLICITY_EXPONENT * below)


def giant_occurrence_factor(fe_h):
    """Giant-planet occurrence at `fe_h` against the disk's: `giant_planet_factor` taken from the same floor, so 1 at
    and above it."""
    return min(1.0, giant_planet_factor(fe_h) / giant_planet_factor(floor_fe_h()))


_HERMITE_NODES = ((0.0, 8.0 / 15.0), (1.3556261799742659, 0.2220759220056126), (-1.3556261799742659, 0.2220759220056126),
                  (2.8569700138728056, 0.011257411327720691), (-2.8569700138728056, 0.011257411327720691))
"""tuple: `(node, weight)` of the five-point Gauss-Hermite rule for a standard normal (weights add up to 1)."""


def expected_factor(function, peaks):
    """
    The mean of `function([Fe/H])` over a stellar population's [Fe/H] distribution `peaks`, `((weight, mean, sigma),
    ...)` as the sector profile holds it, each draw clamped to `tuning.METALLICITY_RANGE_FE_H` as `sample_fe_h` clamps
    it. A closed form (five-point Gauss-Hermite on each normal), so the profile stays a pure function.
    """
    low, high = tuning.METALLICITY_RANGE_FE_H
    total = sum(weight for weight, _mean, _sigma in peaks)
    if total <= 0.0:
        return function(0.0)
    return sum(weight * node_weight * function(min(max(mean + sigma * node, low), high))
               for weight, mean, sigma in peaks for node, node_weight in _HERMITE_NODES) / total
