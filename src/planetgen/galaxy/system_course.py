# planetgen/galaxy/system_course.py

"""
A course inside a star system (NAV.27, docs/design/orbits/course-avoidance.md section 6)
=======================================================================================

Inside a system the planets and moons move, so the course is planned with every body where it is when the ship leaves
(`physics.body_positions.positions_at`), and a destination that moves is aimed at where it will be when the ship gets there.
The bending itself is the planner of NAV.26 (`galaxy.avoid.plan_course`); this module places its spheres and its two ends:

- each body's keep-out sphere (its Hill radius, a star's radiation radius: `db.obstacles.system_keep_outs`) is centered on the body at
  the departure time. The sphere of the destination is centered where the destination will be, and a sphere that holds an end is
  dropped by the planner, so a moon's course inside its planet's Hill sphere is a straight line;
- the ship flies at one speed (`tuning.COURSE_SYSTEM_WARP_FACTOR`), so the arrival time is the course's length over that speed. The end
  is aimed at the destination at that time and the course planned again, `tuning.COURSE_SYSTEM_AIM_STEPS` times at most;
- a slow ship (below about 0.01 c) would need the obstacles at the time the ship passes each, which this does not do: warp and fold
  speeds are 1 c and more, and the planets move 5 to 48 km/s.

Everything is kilometers in the scene's frame (the barycenter of the system's stars at the origin), `frame` "system".
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import avoid
from planetgen.physics import constants
from planetgen.physics.body_positions import positions_at

_SECONDS_PER_YEAR = constants.SECONDS_PER_YEAR


def _position(end, positions):
    """Where an end is: the body's place at the departure time if the scene has it, else the `position` given with the end (a system
    itself is the origin, a belt is a point on its ring, the edge of the heliopause is a point), in kilometers."""
    ref = end.get("ref")
    if ref in positions:
        return np.asarray(positions[ref], dtype=float)
    return np.asarray(end.get("position") or (0.0, 0.0, 0.0), dtype=float)


def _spheres(keep_outs, positions, end_ref, end_at):
    """The planner's spheres: every body with a keep-out radius, centered on where it is at the departure time and, for the destination,
    where it will be on arrival."""
    out = []
    for ref, radius in keep_outs.items():
        if ref not in positions or not radius or radius <= 0.0:
            continue
        center = end_at if ref == end_ref else positions[ref]
        out.append({"ref": ref, "center": list(map(float, center)), "radius": float(radius)})
    return out


def plan_system_course(scene, keep_outs, start, end, speed_kms, departure_years=0.0):
    """
    The `adjusted` course (course-avoidance.md section 7) inside one system, in kilometers.

    Args:
        scene (dict): The system's scene (`web.maps.systemscene.build_scene`).
        keep_outs (dict): `{ref: radius_km}` for the bodies that keep a ship out (`db.obstacles.system_keep_outs`); a body without one is
            not an obstacle.
        start, end (dict): `{ref, position}`: the body's reference (`planet:7`) and, for an end that is not a body of the scene (the
            system itself, a belt, a point on the heliopause), its place in km. A body of the scene is placed by the scene.
        speed_kms (float): The speed the course is flown at, km/s.
        departure_years (float): When the ship leaves, in years after the scene's epoch (0: when the bodies were last stepped).

    Returns:
        dict: `galaxy.avoid.plan_course`'s `adjusted` (frame "system", unit "km") plus `departure_years`, `arrival_years`, `speed_kms`,
            `travel_s` (the time the course takes at that speed) and `aimed_at` (where the destination is aimed at, km).
    """
    if not speed_kms or speed_kms <= 0.0:
        raise ValueError("a course needs a speed above zero")
    positions = positions_at(scene, departure_years)
    start_at = _position(start, positions)
    end_ref = end.get("ref")
    moves = end_ref in positions
    end_at = _position(end, positions)
    def plan(aimed_at):
        return avoid.plan_course(start_at, aimed_at, _spheres(keep_outs, positions, end_ref, aimed_at), frame="system", unit="km")
    result = plan(end_at)
    for _step in range(tuning.COURSE_SYSTEM_AIM_STEPS if moves else 0):
        travel_years = result["length"] / speed_kms / _SECONDS_PER_YEAR
        aimed = np.asarray(positions_at(scene, departure_years + travel_years)[end_ref], dtype=float)
        moved = float(np.linalg.norm(aimed - end_at))
        end_at = aimed
        result = plan(end_at)
        if moved <= 1.0e-6 * max(result["length"], 1.0) or math.isclose(moved, 0.0):
            break
    travel_s = result["length"] / speed_kms
    result.update(departure_years=departure_years, arrival_years=departure_years + travel_s / _SECONDS_PER_YEAR,
                  speed_kms=speed_kms, travel_s=travel_s, aimed_at=end_at.tolist())
    return result
