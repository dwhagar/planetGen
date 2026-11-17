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

**The tilt (GEN.1492, second part).** A companion whose orbit plane is tilted to the belt's turns each member's orbit plane round the companion's own orbit normal at the nodal rate
`A_c` (the same first-order rate as the apsidal one, `companion_rates`), while the star's planets turn it round the planets' plane at their nodal rate `R` (`planet_nodal_rate`). To first order in the angles the two pulls
add as the rate-weighted sum of the two normals, so the member's plane precesses round the **forced axis** `n_f = (R z + sum A_c n_c) / |R z + sum A_c n_c|` (`z` the planets' plane, `n_c` the companion's normal,
`Companion.normal`): a heavy giant holds a belt flat and a bare belt follows the companion. The angle between the member's drawn orbit normal and `n_f` does not change as it precesses, so after a long time the normal
stands anywhere on the cone of that half-angle about `n_f`, at an azimuth the draw gives it (`tilted_orbit`): the member's inclination is no longer a draw about the planets' plane alone. A forced axis farther than
`INCLINATION_LIMIT_RAD` from `z` is outside the first-order theory (and near the Kozai-Lidov limit), and the member keeps its drawn orbit.

**What it leaves out (stated; stand-in rows in the gap plan, section 6.256).** Kozai-Lidov pumping (a tilt past the limit above); resonances with the
companion; the orbit the companion had while the disk lasted; planets are not forced by the companion here. The theory is held to `alpha` of `ALPHA_LIMIT` and an `e_forced` under
`secular_chain.LINEAR_ECCENTRICITY_LIMIT`; past them a member keeps its drawn orbit (`forced_vector` gives zero).
"""

import math

import numpy as np

from planetgen.physics import constants, resonant_chains, secular_chain, spin, vec3
from planetgen.util import draw

ALPHA_LIMIT = 0.3
"""float: The largest `a / a'` the first-order theory is used at (the belt is cut inside the companion's tide and the Holman and Wiegert window, `disk_truncation.star_disk_limit_au`, which keeps a belt under
about this; beyond it the series in `alpha` and the 1:1 terms of the disturbing function no longer converge fast). [Stand-in 6.256 FE4, tier T6, tweak L1.]"""


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


def planet_rate(a_au, planet_axes_au, planet_masses_kg, host_mass_kg, ring_axes_au=(), ring_masses_kg=()):
    """The apsidal rate the star's planets and relativity give a member of no mass at `a_au`, radians per year (`secular_chain.matrices` with the member as a body of no mass). The rings of the
    star's belts (`secular_chain.belt_rings`, GEN.1464) pull the member too, those within `secular_chain.BELT_ALPHA_LIMIT` of it left out."""
    ring_axes, ring_masses = secular_chain.clear_rings(list(ring_axes_au), list(ring_masses_kg), [a_au])
    if not planet_axes_au and not ring_axes:
        return float(secular_chain.relativistic_rates_rad_yr([a_au], host_mass_kg, [0.0])[0])
    A, _B = secular_chain.matrices([a_au, *planet_axes_au, *ring_axes], [0.0, *planet_masses_kg, *ring_masses], host_mass_kg, rows=1)
    return float(A[0, 0])


def periapsis_rad(host_mass_kg, member_no, axis_au):
    """The direction of a companion's periapsis, radians (one angle for a system, drawn from a stream of its own: the generator stores no orientation of a star's orbit). [Stand-in 6.256 FE2, tier T6, tweak L1.]"""
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


INCLINATION_LIMIT_RAD = 0.5
"""float: The farthest the forced axis may stand from the planets' plane (radians, 28.6 degrees) for the first-order theory of the tilt to be used; a companion tilted more than this is near the limit
(39.2 degrees) past which the Kozai-Lidov pumping of eccentricity replaces a steady precession, which this module does not model. A member of such a belt keeps its drawn orbit. [Stand-in 6.256 FE5, tier T6, tweak L1.]"""


def planet_nodal_rate(a_au, planet_axes_au, planet_masses_kg, host_mass_kg, ring_axes_au=(), ring_masses_kg=()):
    """The nodal rate, radians per year (a positive number), the star's planets give a member of no mass at `a_au`: the diagonal of the nodal matrix `B` of `secular_chain.matrices` with the member as a
    body of no mass, the rings of the star's belts counted as in `planet_rate`. Zero without planets or rings; relativity does not turn a node."""
    ring_axes, ring_masses = secular_chain.clear_rings(list(ring_axes_au), list(ring_masses_kg), [a_au])
    if not planet_axes_au and not ring_axes:
        return 0.0
    _A, B = secular_chain.matrices([a_au, *planet_axes_au, *ring_axes], [0.0, *planet_masses_kg, *ring_masses], host_mass_kg, rows=1)
    return abs(float(B[0, 0]))


def forced_normal(a_au, host_mass_kg, companions, planet_nodal_rate_rad_yr):
    """
    The unit normal `(x, y, z)` of the plane a member at `a_au` precesses about (the forced axis), from `companions` (each `(mass_kg, axis_au, eccentricity, normal)`, the normal a unit vector in the
    planets' frame) and the nodal rate the planets give it. `None` for a member the theory does not reach: no companion with a normal, none inside `ALPHA_LIMIT`, or a forced axis farther than
    `INCLINATION_LIMIT_RAD` from the planets' plane (the member keeps its drawn orbit).
    """
    pull = (0.0, 0.0, max(planet_nodal_rate_rad_yr, 0.0))
    pulled = False
    for mass, axis, eccentricity, normal in companions:
        if normal is None:
            continue
        rate, _push = companion_rates(a_au, host_mass_kg, mass, axis, eccentricity)
        if rate > 0.0:
            pull = vec3.add(pull, vec3.scale(normal, rate))
            pulled = True
    if not pulled or not vec3.norm(pull) > 0.0:
        return None
    axis = vec3.unit(pull)
    if axis[2] <= 0.0 or math.acos(min(axis[2], 1.0)) > INCLINATION_LIMIT_RAD:
        return None
    return axis


def tilted_orbit(inclination_deg, node_deg, axis, phase_rad):
    """
    `(inclination_deg, node_deg)` of a member whose drawn orbit (`inclination_deg`, `node_deg`, against the planets' plane) has precessed about the forced `axis` for a long time, at `phase_rad` round it: the
    member's orbit normal on the cone about `axis` that holds its drawn normal. The angle to `axis` is kept, so an `axis` of `(0, 0, 1)` gives back the drawn inclination and a phase-turned node.
    """
    normal = spin.orbit_normal(inclination_deg, node_deg)
    along = vec3.dot(normal, axis)
    across = vec3.add(normal, vec3.scale(axis, -along))
    reach = vec3.norm(across)
    if not reach > 0.0:
        return inclination_deg, node_deg
    first = vec3.scale(across, 1.0 / reach)
    second = vec3.cross(axis, first)
    turned = vec3.add(vec3.scale(axis, along), vec3.add(vec3.scale(first, reach * math.cos(phase_rad)), vec3.scale(second, reach * math.sin(phase_rad))))
    x, y, z = vec3.unit(turned)
    inclination = math.degrees(math.acos(max(-1.0, min(1.0, z))))
    node = math.degrees(math.atan2(x, -y)) % 360.0 if inclination > 0.0 else node_deg
    return inclination, node


def companion_terms(axes_au, host_mass_kg, companions):
    """
    `(apsidal, forcing)` of a companion star on planets at `axes_au` (GEN.1499): for each planet the sum of the companions' `A_c` (radians per year, what the companion adds to the planet's own
    apsidal rate, the diagonal of `secular_chain`'s `A`) and the sum of `B_c exp(i varpi_c)` (the push along the companions' periapses, a complex number, radians per year). `companions` is as for
    `forced_vector`. A planet beyond `ALPHA_LIMIT` of a companion gets nothing from it (`companion_rates`). The pull is that of a body of no mass (`companion_rates`): a planet's own mass changes the companion's pull on it by the planet-to-host mass ratio, left out at first order.
    """
    apsidal = [0.0] * len(axes_au)
    forcing = [0j] * len(axes_au)
    for j, axis in enumerate(axes_au):
        for mass, companion_axis, eccentricity, periapsis in companions:
            rate, push = companion_rates(axis, host_mass_kg, mass, companion_axis, eccentricity)
            apsidal[j] += rate
            forcing[j] += push * complex(math.cos(periapsis), math.sin(periapsis))
    return tuple(apsidal), tuple(forcing)


def forced_planet_vectors(axes_au, masses_kg, host_mass_kg, companions, ring_axes_au=(), ring_masses_kg=()):
    """
    The forced eccentricity vectors `z_eq` (complex, one for each planet) of planets at `axes_au` with masses `masses_kg` under `companions` (GEN.1499): the fixed point `A z = g` of the planets'
    secular exchange (`secular_chain.matrices`, relativity and the rings of the star's belts included) with the companions' `A_c` added to its diagonal and their push `g` (`companion_terms`) on the
    right. One planet and no rings give `forced_vector` exactly. The planets stay circling `z_eq` with their free eccentricity as it was; `None` for a system the theory does not reach: no
    companion pull on any planet, a singular exchange, or a forced eccentricity past `secular_chain.LINEAR_ECCENTRICITY_LIMIT` for any planet.
    """
    apsidal, forcing = companion_terms(axes_au, host_mass_kg, companions)
    if not any(abs(f) > 0.0 for f in forcing):
        return None
    count = len(axes_au)
    ring_axes, ring_masses = secular_chain.clear_rings(list(ring_axes_au), list(ring_masses_kg), list(axes_au))
    A, _B = secular_chain.matrices([*axes_au, *ring_axes], [*masses_kg, *ring_masses], host_mass_kg, rows=count)
    block = A[:count, :count] + np.diag(apsidal)
    push = np.array(forcing, dtype=complex)
    # the rings are held at rest (zero eccentricity, as in `secular_chain._drive`): they pull on the planets and are not moved by them, so they only add to the planets' diagonal
    if np.linalg.matrix_rank(block) < count:
        return None
    vectors = np.linalg.solve(block, push)
    if float(np.max(np.abs(vectors))) > secular_chain.LINEAR_ECCENTRICITY_LIMIT:
        return None
    return tuple(complex(v) for v in vectors)


def stellar_terms(axes_au, host_mass_kg, companions):
    """
    The `secular_chain.StellarTerms` of `companions` for planets at `axes_au` (GEN.1499): the apsidal rate and the push of `companion_terms`, and, for each companion whose plane is known, the same rate
    on the node and the push `A_c w_c` of its plane (`w_c = sin(inc_c) exp(i node_c) = -n_y + i n_x` for its unit normal `n`). `companions` are `(mass_kg, axis_au, eccentricity, periapsis_rad, normal)`, the
    normal a unit vector in the planets' frame or None. A companion with no plane pulls the planets' periapses and no node; the member-of-no-mass theory of `companion_rates` holds for the planets.
    """
    apsidal, push = companion_terms(axes_au, host_mass_kg, [c[:4] for c in companions])
    nodal = [0.0] * len(axes_au)
    tilt = [0j] * len(axes_au)
    for mass, companion_axis, eccentricity, _periapsis, normal in companions:
        if normal is None:
            continue
        w = complex(-normal[1], normal[0])
        for j, axis in enumerate(axes_au):
            rate, _push = companion_rates(axis, host_mass_kg, mass, companion_axis, eccentricity)
            nodal[j] += rate
            tilt[j] += rate * w
    return secular_chain.StellarTerms(apsidal, push, tuple(nodal), tuple(tilt))
