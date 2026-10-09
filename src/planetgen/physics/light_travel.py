# planetgen/physics/light_travel.py

"""
Light-travel positions (VIEW.5): where an object appears to be to an
observer. Light takes time to cross the distance, so an observer never
sees where a distant object is now, only where it was when the light left
it: at the retarded time `t - tau`, where `tau` solves
`|position(t - tau) - observer| = c * tau`.

`apparent_position` takes the object's position as a function of time and
finds `tau` by iteration, which converges for any object slower than light
(each step shrinks the error by about `speed / c`). A star 1 ly away at
rest is seen where it is, a year ago; one moving away from the observer is
seen closer than it is.
"""

import math

from planetgen.physics.constants import SPEED_OF_LIGHT_M_S

MAX_ITERATIONS = 200
TOLERANCE = 1e-12
"""float: The iteration stops when `tau` changes by less than this fraction of itself (or of a second)."""


def uniform_motion(position, velocity, epoch_s=0.0):
    """
    A `position_at(t)` for an object moving in a straight line: at
    `position` (metres) at time `epoch_s` (seconds), moving at `velocity`
    (metres a second).
    """
    return lambda t: tuple(p + v * (t - epoch_s) for p, v in zip(position, velocity))


def apparent_position(position_at, observer, now_s=0.0, speed_m_s=SPEED_OF_LIGHT_M_S):
    """
    Where an observer at `observer` (metres) sees an object at time
    `now_s`.

    Args:
        position_at (callable): `position_at(t)` -> `(x, y, z)` metres, the
            object's position at time `t` seconds; must be defined for
            times before `now_s`.
        observer (tuple): The observer's `(x, y, z)`, metres, held fixed.
        now_s (float): The time of the observation, seconds.
        speed_m_s (float): The speed of light, metres a second.

    Returns:
        tuple: `(position, light_time_s)` -- the apparent position
            (metres) and how long ago the light left it.

    Raises:
        ValueError: If the light time does not settle (an object moving at
            or above light speed toward the observer).
    """
    tau = 0.0
    for _ in range(MAX_ITERATIONS):
        position = position_at(now_s - tau)
        new_tau = math.dist(position, observer) / speed_m_s
        if abs(new_tau - tau) <= TOLERANCE * max(new_tau, 1.0):
            return position_at(now_s - new_tau), new_tau
        tau = new_tau
    raise ValueError("the light time did not settle: the object is moving too fast toward the observer")


def apparent_position_of(body, observer, now_s=None, speed_m_s=SPEED_OF_LIGHT_M_S):
    """
    `apparent_position` for a `SpatialPosition3D` (GEN.74) moving in a
    straight line at its galactic velocity: its galactic position holds at
    its epoch (`now_s` when it has none), and the observation is made at
    `now_s` (its epoch by default).

    Args:
        body (SpatialPosition3D): The object.
        observer (tuple): The observer's galactic `(x, y, z)`, metres.
        now_s (float, optional): The time of the observation, unix seconds.

    Returns:
        tuple: `(position, light_time_s)` as `apparent_position`, the
            position in metres in the galactic frame.
    """
    epoch = body.epoch_unix
    if now_s is None:
        now_s = 0.0 if epoch is None else epoch
    if epoch is None:
        epoch = now_s
    unit = body.length_unit_m
    position = tuple(v * unit for v in body.get_coordinates("galactic", "cartesian"))
    velocity = body.get_velocity_vector("galactic")
    return apparent_position(uniform_motion(position, velocity, epoch), observer, now_s, speed_m_s)
