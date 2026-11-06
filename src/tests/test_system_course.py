# tests/test_system_course.py

"""
NAV.27: the course inside a star system (`galaxy/system_course.plan_system_course`): every body where it is when the ship leaves, the
destination aimed at where it will be on arrival, the bending the planner of NAV.26 does; and `db.obstacles.system_keep_outs`, the radius each
body keeps out. The planner tests build a small scene by hand (`tests/test_body_positions.state_orbit`); the keep-out test needs MySQL.
"""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.db import obstacles, store
from planetgen.galaxy import keepout
from planetgen.galaxy.system_course import plan_system_course
from planetgen.physics.body_positions import positions_at
from planetgen.physics.constants import AU_TO_KM as AU_KM
from planetgen.galaxy.navigation import warp_speed_c
from planetgen.physics import constants
from tests.test_body_positions import SCENE_EPOCH, state_orbit

C_KMS = constants.SPEED_OF_LIGHT_KMS


def scene():
    """One star at the origin and three planets on circular orbits: planet:1 at +1 AU, planet:2 at -1 AU (the star lies between them),
    planet:3 at 0.5 AU on the +y side; planet:2 has a moon."""
    return {
        "epoch_unix": SCENE_EPOCH,
        "tree": {"nodes": [], "members": [{"ref": "star:1", "member_no": 0, "parent_node_id": None, "parent_node_no": None, "slot": 0,
                                           "mass_solar": 1.0}], "origin": None},
        "stars": [{"ref": "star:1", "id": 1, "orbit": None}],
        "planets": [
            {"ref": "planet:1", "id": 1, "orbit": state_orbit("star:1", 1.0, 1.0, 0.0, 0.0, 0.0, SCENE_EPOCH), "moons": []},
            {"ref": "planet:2", "id": 2, "orbit": state_orbit("star:1", 1.0, 1.0, 0.0, 0.0, 180.0, SCENE_EPOCH),
             "moons": [{"ref": "moon:1", "id": 1, "orbit": state_orbit("planet:2", 0.002, 0.01, 0.0, 0.0, 0.0, SCENE_EPOCH)}]},
            {"ref": "planet:3", "id": 3, "orbit": state_orbit("star:1", 0.5, 0.35, 0.0, 0.0, 90.0, SCENE_EPOCH), "moons": []},
        ],
        "comets": [],
    }


KEEP = {"star:1": 0.2 * AU_KM, "planet:1": 1.0e6, "planet:2": 1.0e6, "planet:3": 5.0e5, "moon:1": 1.0e4}
WARP_1 = warp_speed_c(tuning.COURSE_SYSTEM_WARP_FACTOR) * C_KMS


def test_the_star_between_two_planets_bends_the_course_and_it_is_kilometers_in_the_system_frame():
    course = plan_system_course(scene(), KEEP, {"ref": "planet:1"}, {"ref": "planet:2"}, WARP_1)
    assert (course["frame"], course["unit"]) == ("system", "km")
    assert course["avoided"] == ["star:1"]
    assert course["length"] > course["straight_length"] == pytest.approx(2.0 * AU_KM, rel=1e-3)
    assert course["warnings"] == []
    star = np.array([0.0, 0.0, 0.0])
    for point in course["polyline"]:
        assert np.linalg.norm(np.array(point) - star) >= 0.2 * AU_KM * (1.0 - 1e-6)


def test_a_body_the_course_does_not_cross_leaves_it_straight():
    course = plan_system_course(scene(), KEEP, {"ref": "planet:1"}, {"ref": "planet:3"}, WARP_1)
    assert course["avoided"] == [] and course["extra"] == 0.0
    assert course["waypoints"][0]["p"] == pytest.approx(positions_at(scene(), 0.0)["planet:1"])


def test_a_moon_inside_its_planets_sphere_is_a_straight_line():
    # The moon is inside planet:2's keep-out sphere: that sphere holds an end, so it is ignored (and the star is not in the way).
    course = plan_system_course(scene(), KEEP, {"ref": "planet:3"}, {"ref": "moon:1"}, WARP_1)
    assert course["avoided"] == []
    inside = plan_system_course(scene(), KEEP, {"ref": "planet:2"}, {"ref": "moon:1"}, WARP_1)
    assert inside["extra"] == 0.0 and inside["avoided"] == []


def test_the_destination_is_aimed_at_where_it_will_be_on_arrival():
    slow = 300.0    # km/s: ten times a planet's speed, so the destination moves a long way while the ship flies
    near = plan_system_course(scene(), {}, {"ref": "planet:3"}, {"ref": "planet:1"}, slow)
    here = positions_at(scene(), 0.0)["planet:1"]
    assert np.linalg.norm(np.array(near["aimed_at"]) - np.array(here)) > 0.1 * AU_KM
    assert near["arrival_years"] > near["departure_years"] == 0.0
    # The aim is where the planet is at the arrival time the course says (to the steps the planner takes).
    later = positions_at(scene(), near["arrival_years"])["planet:1"]
    assert np.linalg.norm(np.array(near["aimed_at"]) - np.array(later)) < 0.05 * np.linalg.norm(np.array(later) - np.array(here))
    assert near["travel_s"] == pytest.approx(near["length"] / slow)


def test_at_warp_speed_the_destination_has_hardly_moved():
    fast = plan_system_course(scene(), {}, {"ref": "planet:3"}, {"ref": "planet:1"}, WARP_1)
    here = positions_at(scene(), 0.0)["planet:1"]
    assert np.linalg.norm(np.array(fast["aimed_at"]) - np.array(here)) < 1.0e-3 * AU_KM
    assert fast["travel_s"] < 3600.0


def test_departure_time_moves_every_body():
    now = plan_system_course(scene(), KEEP, {"ref": "planet:1"}, {"ref": "planet:2"}, WARP_1)
    later = plan_system_course(scene(), KEEP, {"ref": "planet:1"}, {"ref": "planet:2"}, WARP_1, departure_years=0.25)
    assert later["departure_years"] == 0.25
    assert np.linalg.norm(np.array(later["waypoints"][0]["p"]) - np.array(now["waypoints"][0]["p"])) > 0.5 * AU_KM


def test_an_end_that_is_not_a_body_of_the_scene_takes_the_position_it_is_given():
    edge = [3.0e10, 0.0, 0.0]
    course = plan_system_course(scene(), KEEP, {"ref": "planet:1"}, {"ref": None, "position": edge}, WARP_1)
    assert course["waypoints"][-1]["p"] == pytest.approx(edge)
    system = plan_system_course(scene(), KEEP, {"ref": None}, {"ref": "planet:3"}, WARP_1)   # the system itself is the origin
    assert system["waypoints"][0]["p"] == [0.0, 0.0, 0.0] and system["avoided"] == []


def test_a_course_needs_a_speed():
    with pytest.raises(ValueError):
        plan_system_course(scene(), KEEP, {"ref": "planet:1"}, {"ref": "planet:3"}, 0.0)


def test_every_body_of_a_saved_system_keeps_out_its_own_radius(mysql_config):
    from tests.test_system_scene import _system
    from planetgen.galaxy.sector import SpaceSector
    from planetgen.web.maps.systemscene import build_scene

    sector = SpaceSector("Course Sector", edge_ly=11.5)
    system, cfg = _system(False)
    sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems ORDER BY id").fetchone()["id"]
        built = build_scene(conn, system_id)
        found = obstacles.system_keep_outs(conn, built)
        star = conn.execute("SELECT id, radius_km, luminosity_w FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchone()
        planet = conn.execute("SELECT id, hill_radius_km, radius_km FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchone()
    finally:
        conn.close()
    expected = keepout.star_in_system_keep_out(star["radius_km"], star["luminosity_w"]).radius_km
    assert found[f"star:{star['id']}"] == pytest.approx(expected)
    assert found[f"planet:{planet['id']}"] == pytest.approx(planet["hill_radius_km"] or planet["radius_km"])
    assert all(ref.split(":")[0] in ("star", "planet", "moon") and radius > 0 for ref, radius in found.items())
    # Every body of the scene has a place the planner can use.
    assert set(found) <= set(positions_at(built, 0.0))
