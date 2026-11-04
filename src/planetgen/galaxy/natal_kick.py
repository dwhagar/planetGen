# planetgen/galaxy/natal_kick.py

"""
The share of a galaxy's neutron stars its natal kicks send out of it (GEN.554; docs/design/systems/natal-kick.md)
=========================================================================================================

A neutron star is born with the speed of the star that made it, the galaxy's circular speed `Vc` along the disc, plus a
kick in a random direction. The kick's speed follows a Maxwellian of one-dimensional dispersion `sigma`
(`tuning.NEUTRON_STAR_KICK_SIGMA_KMS`), so the vector sum is a three-dimensional Gaussian about the progenitor's velocity,
and the star leaves the galaxy when its speed is above the escape speed at its birth radius `R`:

    unbound(R) = P(|V + k| > v_esc(R)),    k ~ N(0, sigma^2) in each of three axes, |V| = Vc

which has the closed form of the noncentral chi distribution of three degrees of freedom (`survival`). The massive stars
that leave the neutron stars are born in the thin disc, whose stars at radius `R` are `R exp(-R / L)` per unit radius for
the scale length `L`, so the galaxy's share is the mean over that law (`tuning.NEUTRON_STAR_KICK_BIRTH_RADII` equal-share
radii).

The escape speed is the galaxy's: the halo's Navarro, Frenk and White potential from the mass chain (`mass_chain`: the halo
mass, the virial radius and the concentration) and the stars' mass as a Plummer sphere of the disc's scale length,

    v_esc^2 = 2 (V200^2 c / m(c) ln(1 + x) / x  +  G M_* / sqrt(R^2 + L^2)),   x = R / r_s,  m(c) = ln(1 + c) - c / (1 + c)

with `V200^2 = G M200 / r200` and `r_s = r200 / c`. A halo and disc of the Milky Way's masses and size give an escape speed
of 548 km/s at 8.2 kpc, and 26% of the neutron stars unbound; the shares of the galaxies the generator draws run from none to
three quarters, the small ones losing most.

The neutron stars it ejects leave with the speed they have left once they are out of the well, `sqrt(v^2 - v_esc^2)`: the
kick speeds above the escape speed at the birth radius, weighted by the same law (`ejecta_speeds_kms`, GEN.554).

Every function is pure.
"""

import functools
import math

from planetgen import tuning
from planetgen.physics import constants, unit_registry as units

_msun_to_kg = units.converter("Msun", "kg")

_KM2_PER_S2_PER_SOLAR_MASS_PER_KPC = _msun_to_kg(constants.G) / constants.KILOPARSEC_M / 1.0e6
"""float: G in kpc (km/s)^2 per solar mass."""
_GAMMA2_BISECTIONS = 80


def survival(speed_over_sigma, threshold_over_sigma):
    """
    `P(|a + k| > b)` for a vector `a` of length `speed_over_sigma` and a kick `k` with a standard normal component on each
    of three axes: the upper tail of the noncentral chi distribution of three degrees of freedom,
    `Q(b - a) + Q(b + a) + (phi(b - a) - phi(b + a)) / a` (`Q` the normal upper tail, `phi` its density), and the central
    chi's `erfc(b / sqrt 2) + sqrt(2 / pi) b exp(-b^2 / 2)` for `a` that is 0.
    """
    a, b = speed_over_sigma, threshold_over_sigma
    if b <= 0.0:
        return 1.0
    if a < 1.0e-6:
        return math.erfc(b / math.sqrt(2.0)) + math.sqrt(2.0 / math.pi) * b * math.exp(-0.5 * b * b)
    tail = lambda x: 0.5 * math.erfc(x / math.sqrt(2.0))
    density = lambda x: math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)
    return min(1.0, max(0.0, tail(b - a) + tail(b + a) + (density(b - a) - density(b + a)) / a))


@functools.lru_cache(maxsize=1)
def birth_radii_over_scale_length():
    """The radii, in scale lengths, that divide the disc's stars (`R exp(-R / L)` per unit radius) into
    `tuning.NEUTRON_STAR_KICK_BIRTH_RADII` equal shares, each radius the middle of its share."""
    count = tuning.NEUTRON_STAR_KICK_BIRTH_RADII
    radii = []
    for i in range(count):
        share = (i + 0.5) / count
        low, high = 0.0, 100.0
        for _ in range(_GAMMA2_BISECTIONS):
            middle = 0.5 * (low + high)
            if 1.0 - (1.0 + middle) * math.exp(-middle) < share:
                low = middle
            else:
                high = middle
        radii.append(0.5 * (low + high))
    return tuple(radii)


def escape_speed_kms(chain, disc_scale_length_kpc, radius_kpc):
    """The speed, km/s, a star needs to leave the galaxy from `radius_kpc` (the halo's and the stars' potentials above)."""
    if chain.halo_mass_solar is None:
        halo = 0.0      # a galaxy under the halo threshold holds only its stars
    else:
        concentration = chain.concentration
        halo_scale_kpc = chain.virial_radius_kpc / concentration
        enclosed = math.log1p(concentration) - concentration / (1.0 + concentration)
        x = radius_kpc / halo_scale_kpc
        halo = (_KM2_PER_S2_PER_SOLAR_MASS_PER_KPC * chain.halo_mass_solar / chain.virial_radius_kpc
                * concentration / enclosed * math.log1p(x) / x)
    stars = (_KM2_PER_S2_PER_SOLAR_MASS_PER_KPC * chain.stellar_mass_solar
             / math.sqrt(radius_kpc ** 2 + disc_scale_length_kpc ** 2))
    return math.sqrt(2.0 * (halo + stars))


def neutron_star_unbound_share(chain, circular_speed_kms, disc_scale_length_kpc, sigma_kms=None):
    """
    The share of the neutron stars a galaxy forms that its kicks send out of it: `survival` of the kick over the escape
    speed, averaged over the disc's birth radii (`birth_radii_over_scale_length`).

    Args:
        chain (mass_chain.MassChain): The galaxy's stellar and halo masses, virial radius and concentration.
        circular_speed_kms (float): The progenitor's speed along the disc, the galaxy's `Vc`.
        disc_scale_length_kpc (float): The thin disc's scale length.
        sigma_kms (float, optional): The kick's one-dimensional dispersion; `tuning.NEUTRON_STAR_KICK_SIGMA_KMS` when omitted.
    """
    sigma = tuning.NEUTRON_STAR_KICK_SIGMA_KMS if sigma_kms is None else sigma_kms
    total = 0.0
    for radius in birth_radii_over_scale_length():
        speed = escape_speed_kms(chain, disc_scale_length_kpc, radius * disc_scale_length_kpc)
        total += survival(circular_speed_kms / sigma, speed / sigma)
    return total / len(birth_radii_over_scale_length())


_EJECTA_SPEED_STEP = 0.025
"""float: The step, in kick dispersions, of the speeds integrated above the escape speed."""
_EJECTA_SPEED_SPAN = 12.0
"""float: How far above the escape speed, in kick dispersions, the integral runs (the tail is under 1e-30 by then)."""


def _speed_density(speed_over_sigma, progenitor_over_sigma):
    """The density of `|a + k|` (unit dispersion) at `s`: `s / a (phi(s - a) - phi(s + a))`, and the chi density of three
    degrees of freedom for an `a` of 0."""
    a, s = progenitor_over_sigma, speed_over_sigma
    root = math.sqrt(2.0 * math.pi)
    if a < 1.0e-6:
        return math.sqrt(2.0 / math.pi) * s * s * math.exp(-0.5 * s * s)
    return s / a * (math.exp(-0.5 * (s - a) ** 2) - math.exp(-0.5 * (s + a) ** 2)) / root


def ejecta_speeds_kms(chain, circular_speed_kms, disc_scale_length_kpc, sigma_kms=None):
    """
    The speeds, km/s, the ejected neutron stars of a galaxy have far out, as `tuning.EJECTA_SPEED_QUANTILES` equal-share
    quantiles: each kick speed `s` above the escape speed `v` at its birth radius leaves with `sqrt(s^2 - v^2)`, weighted by
    the speed law (`survival`'s density) and the birth radii (`birth_radii_over_scale_length`). `()` for a galaxy that
    ejects none.
    """
    sigma = tuning.NEUTRON_STAR_KICK_SIGMA_KMS if sigma_kms is None else sigma_kms
    a = circular_speed_kms / sigma
    pairs = []
    for radius in birth_radii_over_scale_length():
        escape = escape_speed_kms(chain, disc_scale_length_kpc, radius * disc_scale_length_kpc) / sigma
        steps = int(_EJECTA_SPEED_SPAN / _EJECTA_SPEED_STEP)
        for i in range(steps):
            s = escape + (i + 0.5) * _EJECTA_SPEED_STEP
            weight = _speed_density(s, a) * _EJECTA_SPEED_STEP
            if weight > 0.0:
                pairs.append((sigma * math.sqrt(max(s * s - escape * escape, 0.0)), weight))
    total = sum(weight for _speed, weight in pairs)
    if total <= 1.0e-12 * len(birth_radii_over_scale_length()):
        return ()
    pairs.sort()
    count = tuning.EJECTA_SPEED_QUANTILES
    speeds, cumulative, index = [], 0.0, 0
    for i in range(count):
        target = (i + 0.5) / count * total
        while index < len(pairs) - 1 and cumulative + pairs[index][1] < target:
            cumulative += pairs[index][1]
            index += 1
        speeds.append(pairs[index][0])
    return tuple(speeds)
