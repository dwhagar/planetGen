# planetgen/api/nav_system.py

"""
The adjusted course inside a star system (NAV.27), for `GET /api/nav?adjusted=1`.

`nav_course` gives a leg inside one system for two ends in one system (`within`), and for leaving (`out`) or entering (`into`) a system on
a course between systems. Each is a straight line in the System Local Frame. With `adjusted=1` this plans it again around the keep-out
spheres of the system's stars, planets and moons, with every body where it is when the ship leaves and the destination aimed at where it
will be on arrival (`galaxy.system_course`). A leg gains `adjusted`; for two ends in one system the course's own `adjusted` is the same.
The scene is read once for each system.
"""

from planetgen import tuning
from planetgen.db import obstacles
from planetgen.galaxy import navigation
from planetgen.galaxy.system_course import plan_system_course
from planetgen.physics import constants, unit_registry as units
from planetgen.web.maps.systemscene import build_scene

_km_to_ly = units.converter("km", "ly")

_SCENE_KINDS = ("star", "planet", "moon", "comet")


def _end(ref, position_km):
    """An end of a leg for the planner: a body of the scene is placed by the scene, anything else (the system itself, a belt, the edge of
    the heliopause) by the position the leg gave."""
    kind = ref.split(":", 1)[0] if ref else None
    return {"ref": ref if kind in _SCENE_KINDS else None, "position": list(position_km)}


def _in_light_years(course):
    """The lengths also in light-years and the warp and fold times for the adjusted length, as `db.query._adjusted_course` does between systems."""
    for name in ("length", "straight_length", "extra"):
        course[name + "_ly"] = _km_to_ly(course[name])
    course["warp_times"] = [leg._asdict() for leg in navigation.warp_travel_times(course["length_ly"])]
    course["fold_times"] = [leg._asdict() for leg in navigation.fold_travel_times(course["length_ly"])]
    return course


def adjust_legs(conn, result):
    """
    Plans the adjusted course of every leg of `result` (`nav_course`'s answer) that lies inside a system and adds it as the leg's
    `adjusted`; for two ends in one system (`scope` "system") it is also the answer's `adjusted`. A leg that cannot be planned (a system
    with no scene) is left as it was.
    """
    scenes = {}
    speed = navigation.warp_speed_c(tuning.COURSE_SYSTEM_WARP_FACTOR) * constants.SPEED_OF_LIGHT_KMS
    for leg in result["legs"]:
        if "system_id" not in leg:
            continue
        system_id = leg["system_id"]
        if system_id not in scenes:
            try:
                scene = build_scene(conn, system_id)
            except ValueError:
                scene = None
            scenes[system_id] = (scene, obstacles.system_keep_outs(conn, scene) if scene else {})
        scene, keep_outs = scenes[system_id]
        if scene is None:
            continue
        course = plan_system_course(scene, keep_outs, _end(leg["from"], leg["from_km"]), _end(leg["to"], leg["to_km"]), speed)
        course = _in_light_years(course)
        course["moving_end"] = leg["to"].split(":", 1)[0] in _SCENE_KINDS
        names = {body["ref"]: body["name"] for body in _bodies(scene)}
        course["names"] = {ref: names[ref] for ref in course["avoided"] + course["ignored_below_resolution"] if ref in names}
        leg["adjusted"] = course
        if result["scope"] == "system":
            result["adjusted"] = course
    return result


def _bodies(scene):
    yield from scene["stars"]
    for planet in scene["planets"]:
        yield planet
        yield from planet["moons"]
