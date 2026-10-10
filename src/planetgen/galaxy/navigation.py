# planetgen/galaxy/navigation.py

"""
Navigation
==========

Pure math for the NAV feature: the course (distance, bearing, mark)
between two absolute 3D positions, and warp and fold travel times along that
distance. This module has no idea where a position came from -- it doesn't
know about sectors, systems, or the database -- callers (the NAV DB query
layer, `queryDb.nav_between`) are responsible for resolving two endpoints
down to a comparable pair of `(x, y, z)` light-year positions in the same
frame, and naming that frame's center, before calling in here.

Course convention (nested reference frames)
--------------------------------------------
Boss's design, `docs/design/navigation-frames.md`. A course reads
"bearing mark mark", both 000-359, where North (bearing 000) points from
the ship toward the current frame's center, flattened onto the frame's
reference plane:

    - Galactic Standard Frame (`FRAME_GALACTIC`): center the galactic core
      (galaxy-frame `(0, 0, 0)`), up galactic +Z. Used between sectors.
    - Sector Local Frame (`FRAME_SECTOR`): center the sector cell's center
      (sector-local `(0, 0, 0)`), up galactic +Z. Used inside one sector.
    - System Local Frame (`FRAME_SYSTEM`): center the central star, up the
      system's angular momentum (ecliptic normal). Used inside a system's
      heliopause. NAV endpoints today are whole systems and phenomena, so
      every course leaves the heliopause and this frame is never chosen;
      `compute_course` supports it for in-system navigation. The edge is
      the system's own heliopause, squeezed by any cloud around it
      (`queryDb.system_detail`'s `heliopause_au`).

With D = target - ship, U the frame's unit up vector, N the ship-to-center
vector with its U part removed (normalized), and E = N x U:

    - Bearing = atan2(D.E, D.N), 0-360 degrees, 000 = North, 090 = East.
    - Elevation = atan2(D.U, sqrt((D.N)^2 + (D.E)^2)), -90 to +90 degrees.
    - Mark = elevation mod 360: 000-090 is up, 270-359 is down (270 is
      straight down); nothing between 091 and 269 appears.

"0 mark 0" therefore points straight at the frame's center. When the ship
sits on the frame's up axis through the center (North undefined), North
falls back to the frame's +X (the galaxy's existing zero meridian, ring
slot 0), or +Y if the up vector is itself +X.

Travel times
------------
`warp_travel_times` and `fold_travel_times` report how long a given
distance takes at several reference factors, from Boss's two speed curves
(`warp_speed_c` and `fold_speed_c`; every coefficient is a named constant
in `program_constants`). At `v` times light-speed, 1 light-year takes
`1 / v` years to cross.

    - Warp `w` (0 < w < 10): speed in c = w^(10/3) + 1 / (1 + e^(-9.3575
      (w - 9.5))) * (198.9 / (10 - w)^0.75 + 1721.7 - w^(10/3)). Plain
      w^(10/3) below about warp 9, then a logistic hand-off to a term that
      climbs without bound toward warp 10.
    - Dimensional fold `F` (0 < F < 10): speed in c = 6 F^4 / (10 - F).
"""

import math
from collections import namedtuple

from planetgen import tuning
from planetgen.util.format import format_period_years

FRAME_GALACTIC = "galactic"
"""str: The Galactic Standard Frame (North = the galactic core)."""

FRAME_SECTOR = "sector"
"""str: The Sector Local Frame (North = the sector's center)."""

FRAME_SYSTEM = "system"
"""str: The System Local Frame (North = the central star)."""

_UP = (0.0, 0.0, 1.0)
_ORIGIN = (0.0, 0.0, 0.0)

Course = namedtuple("Course", ["distance_ly", "bearing_deg", "mark_deg", "elevation_deg", "frame"])
"""
The course from one position to another, in one reference frame.

Attributes:
    distance_ly (float): Straight-line distance, in light-years.
    bearing_deg (float): 0-360 degrees, 0 toward the frame's center
                         (flattened onto its plane), 90 to the East.
    mark_deg (float): `elevation_deg` mod 360 (0-90 up, 270-360 down).
    elevation_deg (float): -90 to +90 degrees above/below the frame's
                           plane.
    frame (str): `FRAME_GALACTIC`, `FRAME_SECTOR` or `FRAME_SYSTEM`.

Bearing, mark and elevation are all `0.0` when the two positions coincide
(direction is undefined at zero distance).
"""

FoldLeg = namedtuple("FoldLeg", ["fold_factor", "velocity_multiple_of_c", "years", "formatted"])
"""
Travel time at one reference dimensional fold factor -- `WarpLeg`'s
counterpart, from `fold_travel_times`.
"""

WarpLeg = namedtuple("WarpLeg", ["warp_factor", "velocity_multiple_of_c", "years", "formatted"])
"""
Travel time at one reference warp factor.

Attributes:
    warp_factor (float): The warp factor (e.g. 1, 3, 6, 9).
    velocity_multiple_of_c (float): Speed at this warp factor, as a
                                    multiple of light-speed.
    years (float): Travel time in years for the distance given to
                   `warp_travel_times`.
    formatted (str): `years`, formatted via `format.format_period_years`.
"""


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _scale(v, k):
    return (v[0] * k, v[1] * k, v[2] * k)


def _norm(v):
    return math.sqrt(_dot(v, v))


def _flatten(v, up):
    """`v` with its component along the unit vector `up` removed."""
    return _sub(v, _scale(up, _dot(v, up)))


def course_between(origin, destination, frame=FRAME_GALACTIC, center=_ORIGIN, up=_UP):
    """
    Computes the course from `origin` to `destination` in one reference
    frame (see the module docstring's "Course convention" section) --
    Boss's `compute_course`.

    All three positions must be in the same coordinate system and unit
    (light-years); this function does no unit conversion or frame
    combination of its own.

    Args:
        origin (tuple): The ship's `(x, y, z)` position.
        destination (tuple): The target's `(x, y, z)` position.
        frame (str): The frame's name, carried into the result. Defaults
            to `FRAME_GALACTIC`.
        center (tuple): The frame's center, which bearing 000 points at.
            Defaults to `(0, 0, 0)`.
        up (tuple): The frame's up vector (any nonzero length). Defaults
            to +Z.

    Returns:
        Course: The distance, bearing, mark and elevation.

    Raises:
        ValueError: If `up` is the zero vector.
    """
    up_length = _norm(up)
    if up_length == 0:
        raise ValueError("a frame's up vector can't be zero")
    u_hat = _scale(up, 1 / up_length)

    displacement = _sub(destination, origin)
    distance = _norm(displacement)
    if distance == 0:
        return Course(distance_ly=0.0, bearing_deg=0.0, mark_deg=0.0, elevation_deg=0.0, frame=frame)

    north = _flatten(_sub(center, origin), u_hat)
    # Relative to the ship-to-center distance, so a ship a few light-years
    # from the axis of a galaxy-sized frame still counts as off-axis.
    if _norm(north) <= 1e-9 * max(_norm(_sub(center, origin)), 1.0):
        fallback = (0.0, 1.0, 0.0) if abs(u_hat[0]) > 0.9 else (1.0, 0.0, 0.0)
        north = _flatten(fallback, u_hat)
    n_hat = _scale(north, 1 / _norm(north))
    e_hat = _cross(n_hat, u_hat)

    d_north = _dot(displacement, n_hat)
    d_east = _dot(displacement, e_hat)
    d_up = _dot(displacement, u_hat)

    bearing = math.degrees(math.atan2(d_east, d_north)) % 360.0
    elevation = math.degrees(math.atan2(d_up, math.hypot(d_north, d_east)))
    return Course(distance_ly=distance, bearing_deg=bearing, mark_deg=elevation % 360.0,
                  elevation_deg=elevation, frame=frame)


def format_course(bearing_deg, mark_deg):
    """
    A course as "000 mark 000": each angle rounded to a whole degree,
    wrapped so 359.6 reads 000, zero-padded to three digits.

    Args:
        bearing_deg (float): 0-360 degrees.
        mark_deg (float): 0-360 degrees (`Course.mark_deg`).

    Returns:
        str: e.g. "045 mark 330".
    """
    return f"{round(bearing_deg) % 360:03d} mark {round(mark_deg) % 360:03d}"


def warp_speed_c(warp_factor):
    """
    Speed at `warp_factor`, as a multiple of light-speed, on Boss's warp
    curve (see the module docstring's "Travel times" section).

    Args:
        warp_factor (float): The warp factor, above 0 and below
            `tuning.WARP_FACTOR_LIMIT` (10).

    Returns:
        float: The speed in multiples of c. Warp 1 is c to within 1e-30.

    Raises:
        ValueError: If `warp_factor` is outside (0, 10).
    """
    if not 0 < warp_factor < tuning.WARP_FACTOR_LIMIT:
        raise ValueError(f"warp factor must be above 0 and below {tuning.WARP_FACTOR_LIMIT}")
    base = warp_factor ** tuning.WARP_VELOCITY_EXPONENT
    blend = 1 / (1 + math.exp(-tuning.WARP_TRANSITION_STEEPNESS
                              * (warp_factor - tuning.WARP_TRANSITION_MIDPOINT)))
    asymptote = (tuning.WARP_ASYMPTOTE_COEFFICIENT
                 / (tuning.WARP_FACTOR_LIMIT - warp_factor) ** tuning.WARP_ASYMPTOTE_EXPONENT)
    return base + blend * (asymptote + tuning.WARP_HIGH_WARP_OFFSET - base)


def fold_speed_c(fold_factor):
    """
    Speed at dimensional fold factor `fold_factor`, as a multiple of
    light-speed: `6 F^4 / (10 - F)`.

    Args:
        fold_factor (float): The fold factor, above 0 and below
            `tuning.FOLD_FACTOR_LIMIT` (10).

    Returns:
        float: The speed in multiples of c.

    Raises:
        ValueError: If `fold_factor` is outside (0, 10).
    """
    if not 0 < fold_factor < tuning.FOLD_FACTOR_LIMIT:
        raise ValueError(f"fold factor must be above 0 and below {tuning.FOLD_FACTOR_LIMIT}")
    return (tuning.FOLD_SPEED_COEFFICIENT * fold_factor ** tuning.FOLD_SPEED_EXPONENT
            / (tuning.FOLD_FACTOR_LIMIT - fold_factor))


def _travel_times(distance_ly, factors, speed_c, make_leg):
    """One leg per factor: `make_leg(factor, speed, years, formatted)`."""
    legs = []
    for factor in factors:
        speed = speed_c(factor)
        years = distance_ly / speed
        legs.append(make_leg(factor, speed, years, format_period_years(years)))
    return legs


MINUTES_PER_YEAR = 365.25 * 24 * 60
"""float: Minutes in a Julian year, for turning a stay at a stop into years."""


def route_travel_times(hop_distances_ly, stay_minutes, factors, speed_c, make_leg):
    """
    NAV.11: one leg per factor for a whole route, the hops timed rest to
    rest at that factor's speed plus `stay_minutes` at every stop between
    the two ends (a route of N hops has N - 1 of them).
    `make_leg(factor, speed, years, formatted)` as `_travel_times`.
    """
    stops = max(len(hop_distances_ly) - 1, 0)
    stay_years = stay_minutes * stops / MINUTES_PER_YEAR
    legs = []
    for factor in factors:
        speed = speed_c(factor)
        years = sum(distance / speed for distance in hop_distances_ly) + stay_years
        legs.append(make_leg(factor, speed, years, format_period_years(years)))
    return legs


def route_warp_times(hop_distances_ly, stay_minutes=0.0, warp_factors=tuning.WARP_FACTORS_FOR_NAV):
    """`route_travel_times` on `warp_speed_c`'s curve: a list of `WarpLeg`."""
    return route_travel_times(hop_distances_ly, stay_minutes, warp_factors, warp_speed_c, WarpLeg)


def route_fold_times(hop_distances_ly, stay_minutes=0.0, fold_factors=tuning.FOLD_FACTORS_FOR_NAV):
    """`route_travel_times` on `fold_speed_c`'s curve: a list of `FoldLeg`."""
    return route_travel_times(hop_distances_ly, stay_minutes, fold_factors, fold_speed_c, FoldLeg)


def warp_travel_times(distance_ly, warp_factors=tuning.WARP_FACTORS_FOR_NAV):
    """
    Computes travel time across `distance_ly` at each of `warp_factors`,
    on `warp_speed_c`'s curve.

    Args:
        distance_ly (float): The distance to travel, in light-years.
        warp_factors (iterable): The warp factors to report travel time
                                 at. Defaults to
                                 `tuning.WARP_FACTORS_FOR_NAV`.

    Returns:
        list: A `WarpLeg` for each of `warp_factors`, in the given order.
    """
    return _travel_times(distance_ly, warp_factors, warp_speed_c, WarpLeg)


def fold_travel_times(distance_ly, fold_factors=tuning.FOLD_FACTORS_FOR_NAV):
    """
    Computes travel time across `distance_ly` at each of `fold_factors`,
    on `fold_speed_c`'s curve.

    Args:
        distance_ly (float): The distance to travel, in light-years.
        fold_factors (iterable): The fold factors to report travel time
                                 at. Defaults to
                                 `tuning.FOLD_FACTORS_FOR_NAV`.

    Returns:
        list: A `FoldLeg` for each of `fold_factors`, in the given order.
    """
    return _travel_times(distance_ly, fold_factors, fold_speed_c, FoldLeg)
