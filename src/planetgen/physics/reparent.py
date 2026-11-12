# planetgen/physics/reparent.py

"""
Reparent: the state-preserving change of a body's parent (GEN.583, docs/design/architecture/unified-system-barycentre.md sections 10.4 and 11.2)
==================================================================================================================================================

A body's truth is its position and velocity **relative to its parent**, in the shared non-rotating axes (`MassObject`, `SpatialPosition3D`'s
"system" frame). When the parent changes (the heaviest member of a system changes, a planet is captured, a planet is thrown out, two bodies merge)
nothing in the sky moves: the body keeps its absolute position and velocity, and only the point it is measured from changes. The relative state about
the new parent is therefore the difference of two states, and the elements follow from it (`state_vectors.elements_from_state`):

    r' = r - R_p        v' = v - V_p        (both in the old parent's frame, R_p and V_p the new parent's state in it)

This module is the pure part, in AU, AU a year and kilograms (the units `MassObject` keeps its motion in): the difference, the new pull, the elements and
whether the body is still bound (`reparent_state`), and the three quantities that must not change in the process, the centre of mass, the total linear
momentum and the orbital angular momentum about the centre of mass (`barycentre`, `momentum`, `angular_momentum`). `generation/reparent.py` applies it to
a `Planet` or moon and makes what depends on the new orbit again (the climate, the class, the spin, the Hill sphere) through the code GEN.1597 already uses.

**The pull** follows the generator's own convention (`MassObject.orbit_mu_au3_per_year2`): the mass of the primary only, `4 pi^2 M / M_sun` in AU^3/yr^2, so
the period `reparent_state` gives is the `calculate_orbital_period_years` of the new semi-major axis and primary and the state read back from the body's
elements is the state it was given.

**Chains.** A body's state in any ancestor's frame is the sum of the hops from the body up to that ancestor, leaf first (`chain_sum`), so a moon's velocity in
its star's frame is its velocity about its planet plus the planet's about the star, each hop small against the sum.
"""

import collections
import math

from planetgen.physics import constants, kepler, state_vectors

Reparented = collections.namedtuple("Reparented", ("position", "velocity", "mu", "orbit", "bound"))
"""The state of a body about its new parent: `position` (AU) and `velocity` (AU/yr) relative to it, `mu` (AU^3/yr^2) of the new primary, the osculating `orbit`
(`state_vectors.elements_from_state`) and `bound` (whether the orbit is elliptical, so the body is still with its new parent)."""


def primary_mu_au3_per_year2(primary_mass_kg):
    """The pull of a primary of `primary_mass_kg` in AU^3/yr^2, by the generator's convention (the primary's mass alone)."""
    return kepler.gravitational_parameter_au3_yr2(primary_mass_kg / constants.SOLAR_MASS_TO_KG)


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scale(a, k):
    return (a[0] * k, a[1] * k, a[2] * k)


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def chain_sum(hops):
    """
    The state in an ancestor's frame of a body whose hops up to it are `hops`, leaf first: the sum of the vectors. Each hop is a position (AU) or a velocity (AU/yr) relative to the parent
    above it, so the sum is the body's position or velocity relative to the last parent. The leaf-first order keeps the small numbers together.
    """
    total = (0.0, 0.0, 0.0)
    for hop in hops:
        total = _add(total, hop)
    return total


def relative_state(position, velocity, parent_position, parent_velocity):
    """The state of a body about another point: the difference of the two states, both given in one frame (AU, AU/yr)."""
    return _sub(position, parent_position), _sub(velocity, parent_velocity)


def reparent_state(position, velocity, new_parent_position, new_parent_velocity, new_primary_mass_kg):
    """
    The state of a body about its new parent, the new parent's state given in the frame the body's own is in (the old parent's), and the orbit it makes.

    Args:
        position (tuple): The body's position about its present parent, AU.
        velocity (tuple): Its velocity about its present parent, AU/yr.
        new_parent_position (tuple): The new parent's position about the present parent, AU.
        new_parent_velocity (tuple): The new parent's velocity about the present parent, AU/yr.
        new_primary_mass_kg (float): The mass of the new primary, kg: the new parent, or for a pair its total.

    Returns:
        Reparented: The relative state, the pull, the osculating orbit and whether the body is bound.

    Raises:
        ValueError: For a body on its new parent or moving straight along the line to it (`state_vectors.elements_from_state`: no orbital plane).
    """
    relative_position, relative_velocity = relative_state(position, velocity, new_parent_position, new_parent_velocity)
    mu = primary_mu_au3_per_year2(new_primary_mass_kg)
    orbit = state_vectors.elements_from_state(relative_position, relative_velocity, mu)
    return Reparented(relative_position, relative_velocity, mu, orbit, orbit["kind"] == "elliptical")


def barycentre(masses, positions):
    """The centre of mass of bodies of `masses` at `positions`, in the frame of the positions: `sum(m r) / sum(m)`."""
    total = math.fsum(masses)
    weighted = (math.fsum(m * p[axis] for m, p in zip(masses, positions)) for axis in range(3))
    return tuple(w / total for w in weighted)


def momentum(masses, velocities):
    """The total linear momentum of bodies of `masses` moving at `velocities`: `sum(m v)`."""
    return tuple(math.fsum(m * v[axis] for m, v in zip(masses, velocities)) for axis in range(3))


def angular_momentum(masses, positions, velocities, about=None):
    """
    The orbital angular momentum of bodies about the point `about` (the centre of mass when not given): `sum(m (r - about) x (v - v_about))`, where `v_about` is the velocity of the centre of mass
    for that point, so the figure is the one of the system's own motion and does not change with the frame the states are given in.
    """
    centre = barycentre(masses, positions) if about is None else about
    centre_velocity = _scale(momentum(masses, velocities), 1.0 / math.fsum(masses))
    total = (0.0, 0.0, 0.0)
    for m, r, v in zip(masses, positions, velocities):
        total = _add(total, _scale(_cross(_sub(r, centre), _sub(v, centre_velocity)), m))
    return total
