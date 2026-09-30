# stellarObjects/navigation.py

"""
Navigation
==========

Pure math for the NAV feature: the course (distance, azimuth, altitude)
between two absolute 3D positions, and warp and fold travel times along that
distance. This module has no idea where a position came from -- it doesn't
know about sectors, systems, or the database -- callers (the NAV DB query
layer, `stellarObjects.navigation` DB layer to come in a later task) are
responsible for resolving two systems down to a comparable pair of `(x, y,
z)` positions in the same frame before calling in here. That keeps this
module usable for both same-sector course-plotting (sector-local positions,
light-years) and cross-sector/galactic course-plotting (galaxy-frame
positions, also converted to light-years by the caller) without needing to
know which case it's in.

Course convention (galactic-plane-relative)
--------------------------------------------
Azimuth and altitude are both measured relative to the galactic plane (the
shared X-Y plane every position in this package's coordinate system --
sector-local and galaxy-frame alike -- is defined against; see
`docs/design/galaxy-coordinate-system.md`), not relative to whatever
direction a ship happens to be facing:

    - Azimuth: the angle of the destination's direction projected onto the
      X-Y plane, measured counterclockwise from the +X axis (`atan2(dy,
      dx)`), 0-360 degrees.
    - Altitude: the angle of elevation of the destination above (positive)
      or below (negative) the X-Y plane (`asin(dz / distance)`), -90 to +90
      degrees.

This is the same "azimuth, mark altitude" framing as a horizon-relative
bearing on a planet's surface, just anchored to the galaxy's own plane
instead of a local horizon -- appropriate since these coordinates already
share one absolute frame across every sector.

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

from . import program_constants
from .utils import years_to_time_string

Course = namedtuple("Course", ["distance_ly", "azimuth_deg", "altitude_deg"])
"""
The course from one position to another.

Attributes:
    distance_ly (float): Straight-line distance, in light-years.
    azimuth_deg (float): 0-360 degrees, counterclockwise from +X in the
                         galactic (X-Y) plane. Undefined (0.0) when the two
                         positions coincide.
    altitude_deg (float): -90 to +90 degrees, elevation above/below the
                          galactic plane. Undefined (0.0) when the two
                          positions coincide.
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
    formatted (str): `years`, formatted via `utils.years_to_time_string`.
"""


# TODO(nav #33): courses become "bearing mark mark-angle", both 0-359, with
# 0 mark 0 toward the frame's center. Rework on Boss's nested reference
# frames (docs/design/navigation-frames.md has the model and pseudocode):
# - Galactic Standard Frame between sectors: North = the galactic core.
# - Sector Local Frame: North = the sector's barycenter, up = galactic +Z.
# - System Local Frame: North = the central star/barycenter, up = the
#   system's angular momentum (ecliptic normal). Hand-off: star -> sector
#   past the heliopause (~120 AU); sector -> galactic when the course
#   crosses a sector boundary (> 4 pc). compute_course(ship, target, frame,
#   center, up) builds the N/E/U basis (N = center direction flattened onto
#   the plane, E = N x U) and returns bearing = atan2(D.E, D.N) mod 360 and
#   mark = atan2(D.U, horizontal).
def course_between(origin, destination):
    """
    Computes the course from `origin` to `destination`: straight-line
    distance plus galactic-plane-relative azimuth and altitude. See the
    module docstring's "Course convention" section.

    Both positions must already be in the same frame and the same unit
    (whichever the caller is working in -- sector-local or galaxy-frame,
    light-years either way) -- this function does no unit conversion or
    frame combination of its own.

    Args:
        origin (tuple): The starting `(x, y, z)` position.
        destination (tuple): The destination `(x, y, z)` position.

    Returns:
        Course: The distance/azimuth/altitude from `origin` to
               `destination`. Azimuth and altitude are both `0.0` if the
               two positions coincide (direction is undefined at zero
               distance).
    """
    dx = destination[0] - origin[0]
    dy = destination[1] - origin[1]
    dz = destination[2] - origin[2]
    distance = math.sqrt(dx * dx + dy * dy + dz * dz)

    if distance == 0:
        return Course(distance_ly=0.0, azimuth_deg=0.0, altitude_deg=0.0)

    azimuth = math.degrees(math.atan2(dy, dx)) % 360
    altitude = math.degrees(math.asin(dz / distance))
    return Course(distance_ly=distance, azimuth_deg=azimuth, altitude_deg=altitude)


def warp_speed_c(warp_factor):
    """
    Speed at `warp_factor`, as a multiple of light-speed, on Boss's warp
    curve (see the module docstring's "Travel times" section).

    Args:
        warp_factor (float): The warp factor, above 0 and below
            `program_constants.WARP_FACTOR_LIMIT` (10).

    Returns:
        float: The speed in multiples of c. Warp 1 is c to within 1e-30.

    Raises:
        ValueError: If `warp_factor` is outside (0, 10).
    """
    if not 0 < warp_factor < program_constants.WARP_FACTOR_LIMIT:
        raise ValueError(f"warp factor must be above 0 and below {program_constants.WARP_FACTOR_LIMIT}")
    base = warp_factor ** program_constants.WARP_VELOCITY_EXPONENT
    blend = 1 / (1 + math.exp(-program_constants.WARP_TRANSITION_STEEPNESS
                              * (warp_factor - program_constants.WARP_TRANSITION_MIDPOINT)))
    asymptote = (program_constants.WARP_ASYMPTOTE_COEFFICIENT
                 / (program_constants.WARP_FACTOR_LIMIT - warp_factor) ** program_constants.WARP_ASYMPTOTE_EXPONENT)
    return base + blend * (asymptote + program_constants.WARP_HIGH_WARP_OFFSET - base)


def fold_speed_c(fold_factor):
    """
    Speed at dimensional fold factor `fold_factor`, as a multiple of
    light-speed: `6 F^4 / (10 - F)`.

    Args:
        fold_factor (float): The fold factor, above 0 and below
            `program_constants.FOLD_FACTOR_LIMIT` (10).

    Returns:
        float: The speed in multiples of c.

    Raises:
        ValueError: If `fold_factor` is outside (0, 10).
    """
    if not 0 < fold_factor < program_constants.FOLD_FACTOR_LIMIT:
        raise ValueError(f"fold factor must be above 0 and below {program_constants.FOLD_FACTOR_LIMIT}")
    return (program_constants.FOLD_SPEED_COEFFICIENT * fold_factor ** program_constants.FOLD_SPEED_EXPONENT
            / (program_constants.FOLD_FACTOR_LIMIT - fold_factor))


def _travel_times(distance_ly, factors, speed_c, make_leg):
    """One leg per factor: `make_leg(factor, speed, years, formatted)`."""
    legs = []
    for factor in factors:
        speed = speed_c(factor)
        years = distance_ly / speed
        legs.append(make_leg(factor, speed, years, years_to_time_string(years)))
    return legs


def warp_travel_times(distance_ly, warp_factors=program_constants.WARP_FACTORS_FOR_NAV):
    """
    Computes travel time across `distance_ly` at each of `warp_factors`,
    on `warp_speed_c`'s curve.

    Args:
        distance_ly (float): The distance to travel, in light-years.
        warp_factors (iterable): The warp factors to report travel time
                                 at. Defaults to
                                 `program_constants.WARP_FACTORS_FOR_NAV`.

    Returns:
        list: A `WarpLeg` for each of `warp_factors`, in the given order.
    """
    return _travel_times(distance_ly, warp_factors, warp_speed_c, WarpLeg)


def fold_travel_times(distance_ly, fold_factors=program_constants.FOLD_FACTORS_FOR_NAV):
    """
    Computes travel time across `distance_ly` at each of `fold_factors`,
    on `fold_speed_c`'s curve.

    Args:
        distance_ly (float): The distance to travel, in light-years.
        fold_factors (iterable): The fold factors to report travel time
                                 at. Defaults to
                                 `program_constants.FOLD_FACTORS_FOR_NAV`.

    Returns:
        list: A `FoldLeg` for each of `fold_factors`, in the given order.
    """
    return _travel_times(distance_ly, fold_factors, fold_speed_c, FoldLeg)
