# planetgen/physics/belt_forcing.py

"""
A companion star's secular forcing of a belt's members (GEN.1492, docs/design/systems/belt-secular-excitation-by-companions.md)
===============================================================================================================================

A belt that circles a star with another star on a bound orbit does not keep the eccentricities it was made with. Over many orbits the companion's pull averages into a slow turning of
each member's periapsis and a **forced** eccentricity along the companion's periapsis: the member's eccentricity vector is a free part (the drawn one) circling a fixed forced part. This module
is the first-order secular theory of that forcing for a body of no mass (Murray and Dermott, *Solar System Dynamics*, sections 7.2 and 7.3, [R3220]), the same theory
`secular_chain` gives the planets, with a star as the perturber.

**The theory.** For a member at `a` round a host of mass `M` with mean motion `n`, and a companion of mass `m'` on an orbit of semi-major axis `a'` and eccentricity `e'` (`alpha = a / a'`):

    A_c = (n / 4) (m' / M) alpha^2 b_{3/2}^{(1)}(alpha) (1 - e'^2)^(-3/2)       how fast the companion turns the member's periapsis
    B_c = (n / 4) (m' / M) alpha^2 b_{3/2}^{(2)}(alpha) e' (1 - e'^2)^(-5/2)     the forcing along the companion's periapsis
    e_forced = | sum_c B_c exp(i varpi_c) | / (A_planets + A_relativity + sum_c A_c)

with the Laplace coefficients `resonant_chains._laplace` and, for a small `alpha`, `e_forced ~ (5/4) alpha e' / (1 - e'^2)`. The rest of the apsidal rate (`A_planets`, with relativity) comes from the
star's planets as a body of no mass among them (`secular_chain.matrices`); a heavy giant beside the belt turns it faster than the companion can force it, which **dilutes** the forcing.
The forced parts of several companions add as vectors (linear theory). The member's orbit at the star's age is the vector sum of its drawn free eccentricity (`e_free`, its drawn
periapsis `varpi_free`) and the forced vector: `combined`.

**What it leaves out (stated; stand-in rows in the gap plan, section 6.256).** Inclination (the orbits of stars are coplanar here, so the forced inclination is zero); Kozai-Lidov pumping; resonances with the
companion; the orbit the companion had while the disk lasted; planets are not forced by the companion here. The theory is held to `alpha` of `ALPHA_LIMIT` and an `e_forced` under
`secular_chain.LINEAR_ECCENTRICITY_LIMIT`; past them a member keeps its drawn orbit (`forced_vector` gives zero).
"""

import math

import numpy as np

from planetgen.physics import constants, resonant_chains, secular_chain
from planetgen.util import draw

ALPHA_LIMIT = 0.3
"""float: The largest `a / a'` the first-order theory is used at (the belt is cut inside the companion's tide and the Holman and Wiegert window, `disk_truncation.star_disk_limit_au`, which keeps a belt under
about this; beyond it the series in `alpha` and the 1:1 terms of the disturbing function no longer converge fast). [Stand-in FE4, tier Marginal educated guess from physics, tweak L1.]"""


def _mean_motion_rad_yr(a_au, host_mass_kg):
    return math.sqrt(constants.G * host_mass_kg / (a_au * constants.AU_TO_M) ** 3) * constants.SECONDS_PER_YEAR


def companion_rates(a_au, host_mass_kg, companion_mass_kg, axis_au, eccentricity):
    """`(A_c, B_c)`, radians per year, of one companion for a member at `a_au` (zero, zero outside the theory: `alpha` beyond `ALPHA_LIMIT`, or an orbit that is not bound)."""
    if not (a_au > 0.0 and axis_au > 0.0 and host_mass_kg > 0.0 and companion_mass_kg > 0.0 and 0.0 <= eccentricity < 1.0):
        return 0.0, 0.0
    alpha = a_au / axis_au
    if alpha > ALPHA_LIMIT:
        return 0.0, 0.0
    weight = 0.25 * _mean_motion_rad_yr(a_au, host_mass_kg) * companion_mass_kg / host_mass_kg * alpha * alpha
    b1 = resonant_chains._laplace(1, alpha, s=1.5)
    b2 = resonant_chains._laplace(2, alpha, s=1.5)
    return weight * b1 * (1.0 - eccentricity ** 2) ** -1.5, weight * b2 * eccentricity * (1.0 - eccentricity ** 2) ** -2.5


def planet_rate(a_au, planet_axes_au, planet_masses_kg, host_mass_kg):
    """The apsidal rate the star's planets and relativity give a member of no mass at `a_au`, radians per year (`secular_chain.matrices` with the member as a body of no mass)."""
    if not planet_axes_au:
        return float(secular_chain.relativistic_rates_rad_yr([a_au], host_mass_kg, [0.0])[0])
    A, _B = secular_chain.matrices([a_au, *planet_axes_au], [0.0, *planet_masses_kg], host_mass_kg)
    return float(A[0, 0])


def periapsis_rad(host_mass_kg, member_no, axis_au):
    """The direction of a companion's periapsis, radians (one angle for a system, drawn from a stream of its own: the generator stores no orientation of a star's orbit). [Stand-in FE2, tier Marginal educated guess from physics, tweak L1.]"""
    return draw.Stream(f"belt-forcing:{host_mass_kg!r}:{member_no}:{axis_au!r}").uniform(0.0, 2.0 * math.pi)


def forced_vector(a_au, host_mass_kg, companions, other_rate_rad_yr):
    """
    The forced eccentricity of a member at `a_au` as a complex number `e_forced * exp(i varpi)`, from `companions` (each `(mass_kg, axis_au, eccentricity, periapsis_rad)`) and the apsidal
    rate `other_rate_rad_yr` the planets and relativity give it. Zero for a member outside the theory (`ALPHA_LIMIT`, or a forced eccentricity past `secular_chain.LINEAR_ECCENTRICITY_LIMIT`).
    """
    total = max(other_rate_rad_yr, 0.0)
    forcing = 0j
    for mass, axis, eccentricity, periapsis in companions:
        rate, push = companion_rates(a_au, host_mass_kg, mass, axis, eccentricity)
        total += rate
        forcing += push * complex(math.cos(periapsis), math.sin(periapsis))
    if total <= 0.0 or forcing == 0j:
        return 0j
    vector = forcing / total
    return 0j if abs(vector) > secular_chain.LINEAR_ECCENTRICITY_LIMIT else vector


def combined(e_free, varpi_free_rad, vector):
    """`(e, varpi)` of the member: the vector sum of the free part and the forced `vector`, `varpi` in [0, 2 pi)."""
    total = e_free * complex(math.cos(varpi_free_rad), math.sin(varpi_free_rad)) + vector
    return abs(total), math.atan2(total.imag, total.real) % (2.0 * math.pi)
