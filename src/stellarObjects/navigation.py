# stellarObjects/navigation.py

"""
Navigation
==========

Pure math for the NAV feature: the course (distance, azimuth, altitude)
between two absolute 3D positions, and warp travel times along that
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

Warp travel time
-----------------
`warp_travel_times` reports how long a given distance takes at several
reference warp factors, using `velocity_multiple_of_c = warp_factor **
program_constants.WARP_VELOCITY_EXPONENT` (i.e. warp factor to the
3.33... power) -- so 1 light-year takes `1 / velocity_multiple_of_c`
years to cross.
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


def warp_travel_times(distance_ly, warp_factors=program_constants.WARP_FACTORS_FOR_NAV):
    """
    Computes travel time across `distance_ly` at each of `warp_factors`,
    using `velocity_multiple_of_c = warp_factor **
    program_constants.WARP_VELOCITY_EXPONENT` -- 1 light-year takes
    `1 / velocity_multiple_of_c` years to cross at that speed.

    Args:
        distance_ly (float): The distance to travel, in light-years.
        warp_factors (iterable): The warp factors to report travel time
                                 at. Defaults to
                                 `program_constants.WARP_FACTORS_FOR_NAV`
                                 (1, 3, 6, 9).

    Returns:
        list: A `WarpLeg` for each of `warp_factors`, in the given order.
    """
    legs = []
    for warp_factor in warp_factors:
        velocity_multiple_of_c = warp_factor ** program_constants.WARP_VELOCITY_EXPONENT
        years = distance_ly / velocity_multiple_of_c
        legs.append(WarpLeg(
            warp_factor=warp_factor,
            velocity_multiple_of_c=velocity_multiple_of_c,
            years=years,
            formatted=years_to_time_string(years),
        ))
    return legs
