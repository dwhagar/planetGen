# planetgen/galaxy/system_mass.py

"""
The mean mass of a galaxy's systems, and of the systems at a place (GEN.545).

`generation.system_mass` works out the mean mass of a system whose stars formed at one age from the generator's own
mass function, remnants and multiplicity. This module averages it over a population's ages (the disk's by its star
formation history, `galaxy.age_law`; the others uniformly over their windows) and then over the populations of a galaxy
or a cell. The average is by systems, because the density counts systems (`galaxy.density`): a population that holds
more of them counts more, whatever it weighs.
"""

import functools

import numpy as np

from planetgen import tuning
from planetgen.galaxy import disk_history
from planetgen.generation import system_mass


def _law_mean_solar(law):
    """The mean system mass of a population whose ages follow `law` (`galaxy.age_law.AgeLaw`)."""
    low, high = law.window
    edges = low + (high - low) * np.arange(system_mass.COHORT_AGES + 1) / system_mass.COHORT_AGES
    weights = [law.weight(float(a), float(b)) for a, b in zip(edges[:-1], edges[1:])]
    middles = 0.5 * (edges[:-1] + edges[1:])
    return system_mass.cohort_mean_system_mass_solar(middles, weights)


@functools.lru_cache(maxsize=64)
def population_means_solar(shape):
    """
    The mean system mass of each population of `shape`'s galaxy, solar masses: `tuning.STELLAR_POPULATION_AGE_RANGES_GY`'s
    names to the mean over the window, with the age law of the disk's three populations (`galaxy.disk_history`).
    """
    laws = disk_history.population_age_laws(shape)
    means = {}
    for name, window in tuning.STELLAR_POPULATION_AGE_RANGES_GY.items():
        means[name] = _law_mean_solar(laws[name]) if name in laws else system_mass.window_mean_system_mass_solar(*window)
    return means


def mean_over_solar(weights, shape):
    """
    The mean system mass of a mix of populations, solar masses: the populations' means averaged by their numbers of
    systems.

    Args:
        weights (dict): Population name -> its density, or any number proportional to its systems.
        shape (galaxy.density.GalaxyShape): The galaxy, for its disk's star formation history.

    Returns:
        float: The mean; the old disk's when the mix holds nothing.
    """
    means = population_means_solar(shape)
    total = sum(weights.values())
    if total <= 0.0:
        return means["old"]
    return sum(value * means[name] for name, value in weights.items()) / total
