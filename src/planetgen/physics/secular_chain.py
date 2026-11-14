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

**Belts (GEN.1464).** A belt's mass pulls on the planets like any body's: it is cut into rings (`belt_rings`, eight equal-width rings with the belt's mass spread over its annulus) that add their pull
to each planet's apsidal and nodal rate. The rings keep a flat circular orbit and are not moved by the planets (`evolve`, `fixed_axes_au`; a belt is thousandths to a few Earth masses against planets of tens to
thousands), so a belt changes how fast the planets' orbits turn and what the neighbours' pulls raise, and takes no part in the exchange of eccentricity. A ring within `BELT_ALPHA_LIMIT` of a planet is left
out (`clear_rings`). Registry BR1 to BR3, section 6.244.

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
BELT_RINGS = 8
"""int: How many equal-width rings a belt is cut into when it pulls on the planets (`belt_rings`, GEN.1464). The pull of a ring sums, so the cut only has to resolve the way the Laplace coefficients change across the belt;
eight rings stay within a few percent of sixty-four for a belt a factor of two wide. Registry BR1, tier T9, tweak L1."""
BELT_ALPHA_LIMIT = 0.9
"""float: A ring of a belt closer than this share to a planet (`a_inner / a_outer` above it) is left out of that planet's secular pull: the Laplace coefficients grow without limit as the ratio nears 1 and
the first-order series stops converging, but the planets' spacing keeps a belt several Hill radii from its neighbours, so only the ring nearest a giant is ever dropped. Registry BR2, tier T9, tweak L1."""


def belt_rings(lower_au, upper_au, mass_kg, rings=BELT_RINGS):
    """
    `(axes_au, masses_kg)` of a belt cut into `rings` equal-width rings (GEN.1464): the belt's mass spread over its annulus at one surface density (`belt_mass.surface_density_earth_per_au2`), so a ring's
    share of the mass goes as the area it covers, `hi^2 - lo^2`, and its axis is the centroid of that area in `a`, `(2 / 3) (hi^3 - lo^3) / (hi^2 - lo^2)`. `([], [])` for a belt with no mass or no width.
    """
    if not (mass_kg > 0.0 and upper_au > lower_au > 0.0 and rings > 0):
        return [], []
    edges = [lower_au + (upper_au - lower_au) * i / rings for i in range(rings + 1)]
    total = upper_au ** 2 - lower_au ** 2
    axes, masses = [], []
    for lo, hi in zip(edges, edges[1:]):
        axes.append(2.0 / 3.0 * (hi ** 3 - lo ** 3) / (hi ** 2 - lo ** 2))
        masses.append(mass_kg * (hi ** 2 - lo ** 2) / total)
    return axes, masses


def clear_rings(axes_au, masses_kg, planet_axes_au, limit=BELT_ALPHA_LIMIT):
    """The rings (`belt_rings`) that are not within `limit` (as `a_inner / a_outer`) of any planet at `planet_axes_au`: `(axes, masses)`."""
    kept = [(a, m) for a, m in zip(axes_au, masses_kg) if all(min(a, p) / max(a, p) <= limit for p in planet_axes_au)]
    return [a for a, _ in kept], [m for _, m in kept]


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


def matrices(a_au, masses_kg, host_mass_kg, rows=None):
    """
    The matrices `(A, B)` of the first-order secular theory, radians per year, for planets of the semi-major axes (AU, any order) and masses (kg) round a host. `A` includes the
    relativistic precession on its diagonal. `rows` (default all) is how many of the leading bodies get a row: the bodies after them (the rings of a belt, GEN.1464) pull on the rows but are not
    moved by them, so their own rows stay zero.
    """
    _check(a_au, masses_kg, host_mass_kg)
    count = len(a_au)
    rows = count if rows is None else rows
    motions = mean_motions_rad_yr(a_au, host_mass_kg, masses_kg)
    A = np.zeros((count, count))
    B = np.zeros((count, count))
    for j in range(rows):
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
        if j < rows:
            A[j, j] += rate
    return A, B


def precession_rates(a_au, masses_kg, host_mass_kg, fixed_axes_au=(), fixed_masses_kg=()):
    """
    `(apsidal, nodal)`: each planet's own precession rates from its siblings and relativity, radians per year (the diagonals of `A` and `B`: the rate at which a planet's periapsis and
    node turn when the other planets' eccentricities and inclinations are zero). The rate a Kozai cycle of a distant perturber has to beat to run (GEN.1447) and the rate the spin
    axis follows against (GEN.1454). `fixed_axes_au` and `fixed_masses_kg` are bodies that pull on the planets and are not moved by them (the rings of a belt, `belt_rings`, GEN.1464); the
    rates returned are the planets'.
    """
    count = len(a_au)
    A, B = matrices([*a_au, *fixed_axes_au], [*masses_kg, *fixed_masses_kg], host_mass_kg, rows=count)
    return tuple(float(x) for x in np.diag(A)[:count]), tuple(float(x) for x in np.diag(B)[:count])


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
    as `z_eq + exp(i A_PP t) (z(0) - z_eq)`. Held ones at rest (`start` zero: the flat rings of a belt) set `z_eq = 0`. The block of planets far from every held body can be singular: the rows of `B_PP` sum to zero without them (the invariable plane
    is a mode of the planets alone) and a distant, light belt's pull on it is below roundoff (GEN.1661); with a push on such a block, `z_eq` is the shortest solution."""
    free = [i for i, h in enumerate(held) if not h]
    fixed = [i for i, h in enumerate(held) if h]
    if not free:
        return start.copy()
    if not fixed:
        return _evolve(matrix, start, years)
    block = matrix[np.ix_(free, free)]
    push = matrix[np.ix_(free, fixed)] @ start[fixed]
    if not np.any(push):
        rest = np.zeros(len(free), complex)         # the held ones are at rest (the flat rings of a belt, GEN.1661): nothing sets a fixed point but the origin
    elif np.linalg.matrix_rank(block) == len(free):
        rest = np.linalg.solve(block, push)
    else:
        rest = np.linalg.lstsq(block, push, rcond=None)[0]     # a block with a zero mode (the invariable plane of the free planets) the held ones barely touch: the shortest solution
    out = start.copy()
    out[free] = -rest + _evolve(block, start[free] + rest, years)
    return out


def evolve(a_au, masses_kg, host_mass_kg, eccentricity, periapsis_longitude, inclination, node_longitude, years, fixed_axes_au=(), fixed_masses_kg=()):
    """
    The orbits `years` later (negative for earlier): the `Exchange` of planets of the semi-major axes (AU), masses (kg), eccentricities, longitudes of periapsis, inclinations and
    longitudes of node (radians). A planet already beyond the limits of the theory (`LINEAR_ECCENTRICITY_LIMIT`, `LINEAR_INCLINATION_LIMIT_RAD`) is held where it is: it keeps its
    orbit and pulls on the others as it is (`_drive`). A lone planet is returned as it is, apart from the relativistic turn of its periapsis.

    `fixed_axes_au` and `fixed_masses_kg` (GEN.1464) are bodies that pull on the planets and keep a flat circular orbit whatever the planets do: the rings of a belt (`belt_rings`). They add their
    pull to each planet's apsidal and nodal rate; being circular and flat they force nothing, and being light (a belt is thousandths of an Earth mass) they are not moved visibly by the planets
    they pull on, so they are held and not evolved (registry BR3). The result is the planets' only.
    """
    count = len(a_au)
    A, B = matrices([*a_au, *fixed_axes_au], [*masses_kg, *fixed_masses_kg], host_mass_kg, rows=count)
    pad = len(fixed_axes_au)
    e0, i0 = np.array([*eccentricity, *([0.0] * pad)], float), np.array([*inclination, *([0.0] * pad)], float)
    z0 = e0 * np.exp(1j * np.array([*periapsis_longitude, *([0.0] * pad)]))
    w0 = np.sin(i0) * np.exp(1j * np.array([*node_longitude, *([0.0] * pad)]))
    held_e = e0 > LINEAR_ECCENTRICITY_LIMIT
    held_i = i0 > LINEAR_INCLINATION_LIMIT_RAD
    held_e[count:] = True
    held_i[count:] = True
    z = _drive(A, z0, years, held_e)[:count]
    w = _drive(B, w0, years, held_i)[:count]
    held_e, held_i = held_e[:count], held_i[:count]
    e_new = np.abs(z)
    inc_new = np.arcsin(np.clip(np.abs(w), 0.0, 1.0))
    linear = bool(np.all(e_new[~held_e] <= LINEAR_ECCENTRICITY_LIMIT) and np.all(inc_new[~held_i] <= LINEAR_INCLINATION_LIMIT_RAD))
    return Exchange(
        tuple(float(x) for x in e_new), tuple(float(x) % (2.0 * math.pi) for x in np.angle(z)),
        tuple(float(x) for x in inc_new), tuple(float(x) % (2.0 * math.pi) for x in np.angle(w)),
        tuple(float(x) for x in _modes(A[:count, :count])[0]), tuple(float(x) for x in _modes(B[:count, :count])[0]), linear)
