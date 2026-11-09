# planetgen/physics/secular_chain.py

"""
The secular exchange of a star's planets: Laplace-Lagrange theory (GEN.1450)
============================================================================

Every planet of a star pulls on every other, and over many orbits the pulls average into a slow exchange of eccentricity and inclination among the planets (the orbits keep their
semi-major axes). A planet beside a heavy giant does not keep the circular, flat orbit it was made on: the giant's orbit pushes its periapsis round and raises its eccentricity to the
**forced** value, and the planets of a crowded chain share an angular-momentum deficit between them. This module is the first-order theory of that exchange (Laplace 1784 and Lagrange 1782;
Murray and Dermott, *Solar System Dynamics*, sections 7.2 to 7.4, [R3220], recalled and derived), for the planets of one star.

**The theory.** With `z = e exp(i varpi)` and `w = sin(inc) exp(i node)` of each planet (`varpi` the longitude of periapsis), to first order in the eccentricities and inclinations

    dz/dt = i A z,      dw/dt = i B w,

where, for a planet `j` of period frequency `n_j` and a neighbour `k` of mass `m_k`, with `alpha = a_inner / a_outer` and `alphabar = alpha` when `k` is the outer of the two and 1 when it is the inner,

    A_jj = n_j/4 sum_k (m_k / M) alpha alphabar b_{3/2}^{(1)}(alpha),    A_jk = -n_j/4 (m_k / M) alpha alphabar b_{3/2}^{(2)}(alpha),
    B_jj = -A_jj,                                                        B_jk = n_j/4 (m_k / M) alpha alphabar b_{3/2}^{(1)}(alpha),

and `b_{3/2}^{(j)}` are Laplace coefficients (`resonant_chains._laplace`, the same trapezoid rule as for the resonant terms). General relativity adds `3 n GM/(c^2 a)` to the
apsidal rate of each planet, which matters for the inner planets of a heavy chain. The matrices have real eigenvalues `g_k` (apsidal) and `s_k` (nodal; one is zero, the invariable plane),
and the solution is the sum of the eigenmodes, `z(t) = V exp(i g t) V^-1 z(0)`, so the planets' eccentricities and inclinations at the star's age follow from where they were made.
The exchange conserves the angular momentum deficit, so the sum of the planets' `m sqrt(a) e^2` stays: the exchange moves eccentricity between planets and does not make it.

**What the rest of the generator reads.** `precession_rates` gives the apsidal and nodal rates of each planet from its siblings and relativity, which the Kozai screen of a planet with a distant
companion (`stability.kozai_verdict`, GEN.1447) and the spin-orbit secular drift (GEN.1454) compare their own rates against; `evolve` gives the orbits at an age; `give_secular_orbits`
(called with the orbits final, `StarSystem._settle_to_final_age`) gives each planet its orbit at the star's age and keeps the system's spacing.

**What it does not do.** The theory is linear: it is good for eccentricities and inclinations below about 0.2 and for planets well short of crossing, which the `Exchange` reports and the
caller keeps (an exchange that would make two orbits cross is not applied, `give_secular_orbits`). A companion star at inclinations beyond a few degrees is outside the expansion and is the
Kozai screen's (GEN.1447); resonant terms are `resonant_chains`'; the stellar quadrupole and tidal precession of the closest planets are not included. The registered stand-ins are in
`docs/design/research/gap-filling-plan.md`, section 6.244.
"""

import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import constants, resonant_chains

_SOLAR_MASS_KG = 1.98847e30
"""float: Only to take the solar mass out of the relativistic term; the star's own mass enters through `constants.G`."""
MAXIMUM_EIGEN_IMAGINARY_SHARE = 1.0e-6
"""float: The most the eigenvalues of `A` or `B` may have as an imaginary part, as a share of the largest: the matrices are similar to symmetric ones, so a larger share is a numerical fault."""


class SecularChainError(ValueError):
    """An input the theory cannot use."""


@dataclass(frozen=True)
class Exchange:
    """
    The orbits of a star's planets after the secular exchange.

    Attributes:
        eccentricity (tuple): Each planet's eccentricity.
        periapsis_longitude (tuple): The longitude of periapsis `varpi`, radians in [0, 2 pi).
        inclination (tuple): The inclination to the reference plane, radians.
        node_longitude (tuple): The longitude of the ascending node, radians in [0, 2 pi).
        apsidal_rates (tuple): The frequencies `g_k` of the apsidal eigenmodes, radians per year.
        nodal_rates (tuple): The frequencies `s_k` of the nodal eigenmodes, radians per year.
        linear (bool): Whether every planet the theory moved stayed inside the limits of the first-order theory (a planet that began beyond them is held).
    """

    eccentricity: tuple
    periapsis_longitude: tuple
    inclination: tuple
    node_longitude: tuple
    apsidal_rates: tuple
    nodal_rates: tuple
    linear: bool


LINEAR_ECCENTRICITY_LIMIT = 0.3
"""float: The first-order theory is held to eccentricities below this (the expansion in `e` of the disturbing function converges for `e` below 0.6627, the Laplace limit, and
the first-order terms are good to a few percent well short of it). Source tier T6 (Murray and Dermott, recalled), tweak level L0."""
LINEAR_INCLINATION_LIMIT_RAD = math.radians(30.0)
"""float: The same for the inclination (the expansion is in `sin i`). Source tier T6, tweak level L0."""


def _laplace_pair(alpha):
    """`(b_{3/2}^{(1)}(alpha), b_{3/2}^{(2)}(alpha))`."""
    return (resonant_chains._laplace(1, alpha, s=1.5), resonant_chains._laplace(2, alpha, s=1.5))


def _check(a_au, masses_kg, host_mass_kg):
    if len(a_au) != len(masses_kg):
        raise SecularChainError("a mass for each semi-major axis")
    if not host_mass_kg > 0.0 or any(not a > 0.0 for a in a_au) or any(not m >= 0.0 for m in masses_kg):
        raise SecularChainError("semi-major axes and the host's mass must be positive and masses not negative")


def mean_motions_rad_yr(a_au, host_mass_kg, masses_kg):
    """Each planet's mean motion, radians per year (Kepler's third law with the planet's own mass)."""
    return [math.sqrt(constants.G * (host_mass_kg + m) / (a * constants.AU_TO_M) ** 3) * constants.SECONDS_PER_YEAR for a, m in zip(a_au, masses_kg)]


def relativistic_rates_rad_yr(a_au, host_mass_kg, masses_kg):
    """The apsidal precession general relativity gives each planet, radians per year: `3 n G M / (c^2 a)`."""
    motions = mean_motions_rad_yr(a_au, host_mass_kg, masses_kg)
    scale = constants.G * host_mass_kg / constants.SPEED_OF_LIGHT_M_S ** 2
    return [3.0 * n * scale / (a * constants.AU_TO_M) for n, a in zip(motions, a_au)]


def matrices(a_au, masses_kg, host_mass_kg):
    """
    The matrices `(A, B)` of the first-order secular theory, radians per year, for planets of the semi-major axes (AU, any order) and masses (kg) round a host. `A` includes the
    relativistic precession on its diagonal.
    """
    _check(a_au, masses_kg, host_mass_kg)
    count = len(a_au)
    motions = mean_motions_rad_yr(a_au, host_mass_kg, masses_kg)
    A = np.zeros((count, count))
    B = np.zeros((count, count))
    for j in range(count):
        for k in range(count):
            if j == k:
                continue
            inner, outer = (a_au[j], a_au[k]) if a_au[j] < a_au[k] else (a_au[k], a_au[j])
            alpha = inner / outer
            bar = alpha if a_au[j] < a_au[k] else 1.0
            weight = 0.25 * motions[j] * masses_kg[k] / (host_mass_kg + masses_kg[j]) * alpha * bar
            b1, b2 = _laplace_pair(alpha)
            A[j, j] += weight * b1
            A[j, k] = -weight * b2
            B[j, j] -= weight * b1
            B[j, k] = weight * b1
    for j, rate in enumerate(relativistic_rates_rad_yr(a_au, host_mass_kg, masses_kg)):
        A[j, j] += rate
    return A, B


def precession_rates(a_au, masses_kg, host_mass_kg):
    """
    `(apsidal, nodal)`: each planet's own precession rates from its siblings and relativity, radians per year (the diagonals of `A` and `B`: the rate at which a planet's periapsis and
    node turn when the other planets' eccentricities and inclinations are zero). The rate a Kozai cycle of a distant perturber has to beat to run (GEN.1447) and the rate the spin
    axis follows against (GEN.1454).
    """
    A, B = matrices(a_au, masses_kg, host_mass_kg)
    return tuple(float(x) for x in np.diag(A)), tuple(float(x) for x in np.diag(B))


def _modes(matrix):
    """`(frequencies, V, V^-1)` of a matrix whose eigenvalues are real."""
    values, vectors = np.linalg.eig(matrix)
    scale = float(np.max(np.abs(values))) if values.size else 0.0
    if scale > 0.0 and float(np.max(np.abs(values.imag))) > MAXIMUM_EIGEN_IMAGINARY_SHARE * scale:
        raise SecularChainError("the secular matrix has complex eigenvalues")
    return values.real, vectors, np.linalg.inv(vectors)


def _evolve(matrix, start, years):
    frequencies, vectors, inverse = _modes(matrix)
    return vectors @ (np.exp(1j * frequencies * years) * (inverse @ start))


def _drive(matrix, start, years, held):
    """`_evolve` for the planets not `held`, which keep their `start` and push on the rest as they are: the rest go round the fixed point those held ones set (`z_eq = -A_PP^-1 A_PH z_H`)
    as `z_eq + exp(i A_PP t) (z(0) - z_eq)`."""
    free = [i for i, h in enumerate(held) if not h]
    fixed = [i for i, h in enumerate(held) if h]
    if not free:
        return start.copy()
    if not fixed:
        return _evolve(matrix, start, years)
    block = matrix[np.ix_(free, free)]
    rest = np.linalg.solve(block, matrix[np.ix_(free, fixed)] @ start[fixed])
    out = start.copy()
    out[free] = -rest + _evolve(block, start[free] + rest, years)
    return out


def evolve(a_au, masses_kg, host_mass_kg, eccentricity, periapsis_longitude, inclination, node_longitude, years):
    """
    The orbits `years` later (negative for earlier): the `Exchange` of planets of the semi-major axes (AU), masses (kg), eccentricities, longitudes of periapsis, inclinations and
    longitudes of node (radians). A planet already beyond the limits of the theory (`LINEAR_ECCENTRICITY_LIMIT`, `LINEAR_INCLINATION_LIMIT_RAD`) is held where it is: it keeps its
    orbit and pulls on the others as it is (`_drive`). A lone planet is returned as it is, apart from the relativistic turn of its periapsis.
    """
    A, B = matrices(a_au, masses_kg, host_mass_kg)
    e0, i0 = np.array(eccentricity, float), np.array(inclination, float)
    z0 = e0 * np.exp(1j * np.array(periapsis_longitude))
    w0 = np.sin(i0) * np.exp(1j * np.array(node_longitude))
    held_e = e0 > LINEAR_ECCENTRICITY_LIMIT
    held_i = i0 > LINEAR_INCLINATION_LIMIT_RAD
    z = _drive(A, z0, years, held_e)
    w = _drive(B, w0, years, held_i)
    e_new = np.abs(z)
    inc_new = np.arcsin(np.clip(np.abs(w), 0.0, 1.0))
    linear = bool(np.all(e_new[~held_e] <= LINEAR_ECCENTRICITY_LIMIT) and np.all(inc_new[~held_i] <= LINEAR_INCLINATION_LIMIT_RAD))
    return Exchange(
        tuple(float(x) for x in e_new), tuple(float(x) % (2.0 * math.pi) for x in np.angle(z)),
        tuple(float(x) for x in inc_new), tuple(float(x) % (2.0 * math.pi) for x in np.angle(w)),
        tuple(float(x) for x in _modes(A)[0]), tuple(float(x) for x in _modes(B)[0]), linear)
