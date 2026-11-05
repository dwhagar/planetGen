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

The sky of a system (VIEW.21, S3 of reports/system-sky-from-the-system-2026-11-04.md) needs this for tens of thousands of stars at once, so the
second half of the module is the same law on arrays (`numpy`, metres and seconds): `retarded_offsets` for the apparent offset of every star from
the observer, `seen_age` for the age a star had when its light left it, and `remnant_seen` for whether the light of an explosion has arrived.
"""

import numpy as np

from planetgen.physics import unit_registry
from planetgen.physics.constants import SPEED_OF_LIGHT_M_S
from planetgen.util import detmath

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
        new_tau = detmath.dist(position, observer) / speed_m_s
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
    unit = unit_registry.unit(body.length_unit).factor
    position = tuple(v * unit for v in body.get_coordinates("galactic", "cartesian"))
    velocity = body.get_velocity_vector("galactic")
    return apparent_position(uniform_motion(position, velocity, epoch), observer, now_s, speed_m_s)


# --- the same law on arrays (VIEW.21) ---------------------------------------------------------------------------------------------------

MIN_PASSES = 2
"""int: Fewest light-time passes on arrays (each shrinks the error by about `speed / c`: 1e-3 for 300 km/s, so two give 1e-6 of the light time)."""

MAX_PASSES = 60
"""int: Most passes before the light time is called unsettled (a star at or above light speed toward the observer)."""

PASS_TOLERANCE = 1e-12
"""float: The passes stop when no star's light time changes by more than this share of itself (or of a second)."""


def _length(vectors):
    """The length of each row of an `(N, 3)` array, summed left to right (no BLAS, `util.detmath`)."""
    return np.sqrt(vectors[:, 0] * vectors[:, 0] + vectors[:, 1] * vectors[:, 1] + vectors[:, 2] * vectors[:, 2])


def retarded_offsets(offset_m, velocity_m_s, observer_velocity_m_s, elapsed_s, speed_m_s=SPEED_OF_LIGHT_M_S):
    """
    Where each star appears, as an offset from the observer, for all stars at once (the law of `apparent_position` on arrays).

    Every star is a straight-line mover from its stored epoch (position, velocity, epoch). With `offset_m` the star's offset from the observer at
    that epoch, `elapsed_s` the time from it to the sky's time `t`, and `tau` the light time:

        apparent offset = offset + (velocity - observer_velocity) * elapsed - velocity * tau       tau = |apparent offset| / c

    The passes start from `tau = 0` and run until none changes (two at least; more for a star faster than 1,000 km/s).

    Args:
        offset_m (array): `(N, 3)` metres, each star's offset from the observer at the star's own epoch.
        velocity_m_s (array): `(N, 3)` metres a second, each star's velocity.
        observer_velocity_m_s (sequence): `(3,)` metres a second, the observer's velocity.
        elapsed_s (float or array): Seconds from each star's epoch to the sky's time (a number for one shared epoch, or `(N,)`).
        speed_m_s (float): The speed of light.

    Returns:
        tuple: `(offsets, light_time_s)`: the `(N, 3)` apparent offsets in metres and the `(N,)` seconds the light took.

    Raises:
        ValueError: If the light time of some star does not settle (it moves at or above light speed toward the observer).
    """
    offset = np.asarray(offset_m, dtype=float).reshape(-1, 3)
    velocity = np.asarray(velocity_m_s, dtype=float).reshape(-1, 3)
    relative = velocity - np.asarray(observer_velocity_m_s, dtype=float).reshape(1, 3)
    elapsed = np.broadcast_to(np.asarray(elapsed_s, dtype=float), (offset.shape[0],))
    advanced = offset + relative * elapsed[:, None]       # where the star is seen from the observer's side at the sky's time, if the light were instant
    tau = np.zeros(offset.shape[0])
    for number in range(MAX_PASSES):
        apparent = advanced - velocity * tau[:, None]
        new_tau = _length(apparent) / speed_m_s
        settled = bool(np.all(np.abs(new_tau - tau) <= PASS_TOLERANCE * np.maximum(new_tau, 1.0)))
        tau = new_tau
        if settled and number + 1 >= MIN_PASSES:
            return advanced - velocity * tau[:, None], tau
    raise ValueError("the light time did not settle: a star is moving too fast toward the observer")


def seen_age(age_s, light_time_s, phase_end_age_s=None):
    """
    The age each star had when its light left it (VIEW.21): `age - tau`.

    Args:
        age_s (array): `(N,)` seconds, each star's age at the sky's time.
        light_time_s (array): `(N,)` seconds, from `retarded_offsets`.
        phase_end_age_s (array, optional): `(N,)` seconds, the age at which each star leaves its present phase (`nan` or `None` for none stored).

    Returns:
        tuple: `(age_seen_s, formed, other_phase)`. `formed` is False for a star too young to have been seen yet (its light left before it formed): leave it out
            of the sky. `other_phase` is True where the star's phase ended within the light time, so it is seen in the phase before it (its luminosity and
            colour at `age_seen_s` come from the stellar track, only for these few stars).
    """
    age = np.asarray(age_s, dtype=float)
    seen = age - np.asarray(light_time_s, dtype=float)
    formed = seen >= 0.0
    if phase_end_age_s is None:
        return seen, formed, np.zeros(age.shape, dtype=bool)
    end = np.asarray(phase_end_age_s, dtype=float)
    other_phase = formed & (end <= age) & (seen < end)       # nan compares False: a star with no stored end is never in another phase
    return seen, formed, other_phase


def remnant_seen(age_s, light_time_s):
    """Whether the light of each explosion (a supernova remnant of `age_s` seconds) has reached the observer: it has when the explosion is older than the light time."""
    return np.asarray(age_s, dtype=float) >= np.asarray(light_time_s, dtype=float)


def curvature_bound_m(speed_m_s, light_time_s, galactic_radius_m):
    """
    The most a galactic orbit leaves a straight line over a light time, `(v T)^2 / (2 R)` metres: why a star's stored straight-line motion is exact
    enough for its apparent position. A star 10 kpc away (32,600 years of light, 220 km/s, orbit radius 8.2 kpc) is off by 0.003 pc, 3e-7 radians, 0.07 arcsecond (the plan's 7e-8 radians was an arithmetic slip; a screen pixel of the system view is 200 arcseconds).
    """
    return (speed_m_s * light_time_s) ** 2 / (2.0 * galactic_radius_m)
