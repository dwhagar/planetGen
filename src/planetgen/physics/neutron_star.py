# planetgen/physics/neutron_star.py

"""
The radius of a neutron star from its mass (GEN.463, GEN.969).

A neutron star's radius is fixed by the equation of state of dense matter. `physics.neutron_star_structure` holds the radius of a
star of each mass for 21 unified equations of state (GEN.968). Which of them is real is not known, only constrained: NICER measured
PSR J0030+0451 at 1.34 solar masses, 12.71 +1.14 -1.19 km [R3766], and PSR J0740+6620 at 2.08 solar
masses, 13.7 +2.6 -1.5 km, and with the gravitational-wave and nuclear constraints [R3075] give R(1.4) = 12.45 +- 0.65 km and R(2.08) = 12.35 +- 0.75 km.

So each equation of state carries a weight from how well it fits those constraints (`eos_weight`): a normal weight on its R(1.4)
and on its R(2.08) about the Miller et al. values, and for one too soft to hold PSR J0740+6620 (2.08 +- 0.07 solar masses, Fonseca et al.
2021) a normal weight on the mass it falls short by. A star's radius is a weighted quantile of the radii, at its mass, of the equations
of state that hold such a star: the standard-normal `deviation` the caller draws becomes the probability that picks the quantile
(0 gives the weighted median, as the band's centre did before), interpolated between the radii so that it is continuous.
[Stand-in NS4, tier T2, tweak L1: the weights are a likelihood on the three numbers above, not a fit.]
"""

import math
from functools import lru_cache

from planetgen.datasets import neutron_star_matter
from planetgen.physics import neutron_star_structure as structure

MILLER_R14 = (12.45, 0.65)
MILLER_R208 = (12.35, 0.75)
J0740_MASS = (2.08, 0.07)
"""tuples: `(centre, sigma)` of the constraints the weights use: R(1.4) and R(2.08), km [R3075], and the mass of PSR
J0740+6620, solar masses [R1404]."""


def _normal_weight(value, centre, sigma):
    return math.exp(-0.5 * ((value - centre) / sigma) ** 2)


def eos_weight(eos):
    """The weight, relative, of equation of state `eos` (1 to 21) from the constraints in the module note."""
    weight = _normal_weight(structure.radius_km(eos, 1.4), *MILLER_R14)
    radius_heavy = structure.radius_km(eos, J0740_MASS[0])
    if radius_heavy is None:
        return weight * _normal_weight(structure.maximum_mass_solar(eos), *J0740_MASS)
    return weight * _normal_weight(radius_heavy, *MILLER_R208)


@lru_cache(maxsize=None)
def _weights():
    return {eos: eos_weight(eos) for eos in range(1, neutron_star_matter.EOS_COUNT + 1)}


def _radii_and_weights(mass_solar):
    """The radii (ascending) and the matching weights, summing to 1, of the equations of state that hold a star of this mass."""
    pairs = sorted((structure.radius_km(eos, mass_solar), _weights()[eos]) for eos in structure.supporting(mass_solar))
    if not pairs:
        raise ValueError(f"no equation of state holds a neutron star of {mass_solar} solar masses")
    total = sum(weight for _, weight in pairs)
    return [radius for radius, _ in pairs], [weight / total for _, weight in pairs]


def radius_km(mass_solar, deviation=0.0):
    """
    The radius, km, of a neutron star of this mass: the weighted quantile of the radii of the equations of state that hold it, at
    the probability a standard normal `deviation` stands for (0 for the weighted median). Every `deviation` lies between the
    lowest and highest radius (`radius_limits_km`).
    """
    radii, weights = _radii_and_weights(mass_solar)
    probability = 0.5 * (1.0 + math.erf(deviation / math.sqrt(2.0)))
    # Each radius sits at the middle of its share of the probability, and the radius is straight between neighbours.
    middle, cumulative = [], 0.0
    for weight in weights:
        middle.append(cumulative + 0.5 * weight)
        cumulative += weight
    if probability <= middle[0]:
        return radii[0]
    if probability >= middle[-1]:
        return radii[-1]
    for index in range(1, len(radii)):
        if probability <= middle[index]:
            fraction = (probability - middle[index - 1]) / (middle[index] - middle[index - 1])
            return radii[index - 1] + (radii[index] - radii[index - 1]) * fraction
    return radii[-1]


def radius_limits_km(mass_solar):
    """The lowest and highest radius, km, of the equations of state that hold a star of this mass."""
    radii, _ = _radii_and_weights(mass_solar)
    return radii[0], radii[-1]
