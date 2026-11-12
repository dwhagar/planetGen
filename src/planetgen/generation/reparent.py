# generation/reparent.py

"""
Reparent a planet or a moon (GEN.583, docs/design/architecture/unified-system-barycentre.md sections 10.4 and 11.2)
===================================================================================================================

The state-preserving change of a body's parent. The body keeps its absolute position and velocity; only the point they are measured from changes
(`physics/reparent.py`, which has the arithmetic and the quantities that must not change). What is made again is what the new orbit changes, and nothing else:

* The orbit: the semi-major axis, the eccentricity, the plane and the period come from the new relative state, the Hill sphere from the new primary.
* A planet: its zone, air, class, spin, day and night, tides and the moons' own are made again by `validation.reconcile_moved_planet`, the code GEN.1597
  made keep the spin a planet was born with. A moon: its zone, air, class and spin by `validation.reconcile_moved_moon`, the planet's by its new pull.
* What the body is, its mass, radius, composition, birth spin (`_birth_spin`) and birth orbit (`_orbit_birth`), is kept: where a body sits does not change what it is made of.

A body that comes out of the change **not bound** to its new parent (an orbit that is not an ellipse) is raised as `Freed` with the state it would have, before anything
is changed: the caller (an ejection) decides what a free body becomes (`phenomena/starless_host.py`). The three cases of the design note use it: the heaviest member
of a system changes at generation or at an event (a planet passes to the other star of a pair), a captured moon passes to another planet, a planet becomes free.
A planet that becomes a moon, or a moon a planet, changes kind as well as parent and is not made here.

The units are the body's own: AU, km/s and kg.
"""

import math

from planetgen.generation import validation
from planetgen.physics import constants, daylight, ocean_tides, planets as planetPhysics, reparent as kinematics

_AU_YR_TO_KM_S = constants.AU_M / constants.SECONDS_PER_YEAR / 1000.0
"""float: Kilometres a second in one AU a year."""


class Freed(ValueError):
    """A body that would not be bound to its new parent. `state` is the `physics.reparent.Reparented` it would have had."""

    def __init__(self, body, state):
        super().__init__(f"{getattr(body, 'name', 'the body')} would be on a {state.orbit['kind']} orbit about its new parent (eccentricity {state.orbit['eccentricity']:.3f})")
        self.state = state


def _state(body):
    """`(position, velocity)` of `body` about its present parent, AU and AU/yr, in the system's frame."""
    position = (body.position_x, body.position_y, body.position_z)
    velocity = tuple(v / _AU_YR_TO_KM_S for v in (body.velocity_x_kms, body.velocity_y_kms, body.velocity_z_kms))
    return position, velocity


def _apply(body, new_primary_mass_kg, state):
    """Writes the relative state and the orbit it makes onto `body`, with what follows from the orbit alone (period, Hill sphere, speed)."""
    orbit = state.orbit
    body._frame_axis = None            # the elements below are measured in the system's frame, not in an equator's (GEN.1352)
    body.primary_mass_kg = new_primary_mass_kg
    body.distance = orbit["semi_major_axis"]
    body.eccentricity = orbit["eccentricity"]
    body.orbital_inclination_deg = math.degrees(orbit["inclination"])
    body.orbital_ascending_node_deg = math.degrees(orbit["ascending_node"])
    body.period = planetPhysics.calculate_orbital_period_years(body.distance, new_primary_mass_kg)
    body.set_position_au(*state.position)
    body.set_velocity_kms(*(v * _AU_YR_TO_KM_S for v in state.velocity))
    body.orbital_speed_kms = math.sqrt(sum(v * v for v in state.velocity)) * _AU_YR_TO_KM_S
    body.min_update_interval_years = planetPhysics.minimum_update_interval_years(body.period)
    planetPhysics.update_hill_sphere(body)


def reparent_planet(planet, new_parent_position_au, new_parent_velocity_au_yr, new_primary_mass_kg, new_star=None):
    """
    Passes a planet (not a moon) to a new parent: its absolute position and velocity are kept, its orbit and what depends on it follow.

    Args:
        planet (Planet): The planet; its state is about its present star.
        new_parent_position_au (tuple): The new parent's position about the present star, AU.
        new_parent_velocity_au_yr (tuple): The new parent's velocity about the present star, AU/yr.
        new_primary_mass_kg (float): The new primary's mass, kg (the new star, or for a pair the total of the two).
        new_star: The new host `Star`, when it is a different object from the present one (the planet's `star` is set to it).

    Returns:
        physics.reparent.Reparented: The state about the new parent.

    Raises:
        Freed: For a planet that would not be bound to it; the planet is not changed.
        ValueError: For a moon (`reparent_moon`) or a body with no state.
    """
    if planet.is_moon:
        raise ValueError("a moon passes to a new planet with reparent_moon")
    position, velocity = _state(planet)
    state = kinematics.reparent_state(position, velocity, new_parent_position_au, new_parent_velocity_au_yr, new_primary_mass_kg)
    if not state.bound:
        raise Freed(planet, state)
    if new_star is not None:
        planet.star = new_star
    _apply(planet, new_primary_mass_kg, state)
    validation.reconcile_moved_planet(planet)
    return state


def reparent_moon(moon, old_planet, new_planet, new_parent_position_au, new_parent_velocity_au_yr):
    """
    Passes a moon from `old_planet` to `new_planet`: its absolute position and velocity are kept, its orbit about the new planet and what depends on it follow.

    Args:
        moon (Planet): The moon; its state is about `old_planet`, in the system's frame.
        old_planet (Planet): The planet it orbits now.
        new_planet (Planet): The planet it will orbit.
        new_parent_position_au (tuple): `new_planet`'s position about `old_planet`, AU.
        new_parent_velocity_au_yr (tuple): Its velocity about `old_planet`, AU/yr.

    Returns:
        physics.reparent.Reparented: The state about the new planet.

    Raises:
        Freed: For a moon that would not be bound to it; nothing is changed.
        ValueError: For a body that is not one of `old_planet`'s moons.
    """
    if moon not in old_planet.moons:
        raise ValueError("the moon is not one of the old planet's moons")
    position, velocity = _state(moon)
    state = kinematics.reparent_state(position, velocity, new_parent_position_au, new_parent_velocity_au_yr, new_planet.mass)
    if not state.bound:
        raise Freed(moon, state)
    old_planet.moons.remove(moon)
    new_planet.moons.append(moon)
    _apply(moon, new_planet.mass, state)
    validation.reconcile_moved_moon(moon, new_planet)
    for planet in (old_planet, new_planet):
        validation.refresh_planet_reflex(planet, force=True)
    daylight.generate(moon, new_planet)
    ocean_tides.generate(moon, new_planet)
    return state
