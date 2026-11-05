# tests/test_course_planner.py
"""NAV.26: the course planner bends a straight course around keep-out spheres (galaxy/avoid.py)."""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import avoid


def sphere(ref, center, radius):
    return {"ref": ref, "center": list(center), "radius": radius}


def one_sphere_length(a, b, c, r):
    """The closed form of course-avoidance.md section 5.1 at the waypoint radius."""
    r = r * (1 + avoid.WAYPOINT_MARGIN)
    da, db = np.linalg.norm(a - c), np.linalg.norm(b - c)
    cos_phi = (a - c) @ (b - c) / (da * db)
    phi = math.acos(max(-1.0, min(1.0, cos_phi)))
    return math.sqrt(da * da - r * r) + r * (phi - math.acos(r / da) - math.acos(r / db)) + math.sqrt(db * db - r * r)


def clear(course, spheres, margin=tuning.COURSE_POLYLINE_SAGITTA * 1.05):
    """No point of the drawn path lies inside a sphere by more than `margin` of its radius (a chord of an arc dips in by the sagitta)."""
    points = [np.array(p) for p in course["polyline"]]
    ends = (points[0], points[-1])
    # a sphere that holds an end is exempt (a planet inside its star's sphere)
    spheres = [s for s in spheres if all(np.linalg.norm(e - np.array(s["center"])) > s["radius"] * (1 + avoid.CONTAIN_MARGIN) for e in ends)]
    for k in range(len(points) - 1):
        for f in np.linspace(0.0, 1.0, 40):
            p = points[k] * (1 - f) + points[k + 1] * f
            for s in spheres:
                if np.linalg.norm(p - np.array(s["center"])) < s["radius"] * (1 - margin) - 1e-12:
                    return False
    return True


def test_a_clear_line_is_returned_as_it_is():
    course = avoid.plan_course([0, 0, 0], [10, 0, 0], [sphere("star:1", [5, 3, 0], 1.0)])
    assert course["length"] == course["straight_length"] == 10.0
    assert course["extra"] == 0.0
    assert [w["kind"] for w in course["waypoints"]] == ["start", "end"]
    assert course["avoided"] == [] and course["warnings"] == []
    assert course["ruleset"] == {"keepout": avoid.KEEPOUT_RULESET, "planner": avoid.PLANNER_RULESET}


def test_one_sphere_matches_the_closed_form():
    a, b = np.array([0.0, 0, 0]), np.array([40.0, 3.0, 0.0])
    c, r = np.array([12.0, 1.0, 0.5]), 2.0
    course = avoid.plan_course(a, b, [sphere("star:101", c, r)])
    assert course["length"] == pytest.approx(one_sphere_length(a, b, c, r), rel=1e-12)
    assert [w["kind"] for w in course["waypoints"]] == ["start", "tangent_in", "tangent_out", "end"]
    assert course["avoided"] == ["star:101"]
    assert [leg["type"] for leg in course["legs"]] == ["line", "arc", "line"]
    arc = course["legs"][1]
    assert arc["obstacle"] == "star:101" and arc["radius"] == pytest.approx(r * (1 + avoid.WAYPOINT_MARGIN))
    assert course["extra"] == pytest.approx(course["length"] - course["straight_length"])
    for w in course["waypoints"][1:3]:
        assert np.linalg.norm(np.array(w["p"]) - c) == pytest.approx(r * (1 + avoid.WAYPOINT_MARGIN), rel=1e-12)


def test_a_dead_center_sphere_has_the_documented_overhead():
    # section 5.1: 2.0e-2 at r/L = 0.1
    course = avoid.plan_course([0, 0, 0], [100, 0, 0], [sphere("a", [50, 0, 0], 10.0)])
    assert course["length"] / course["straight_length"] - 1 == pytest.approx(2.0e-2, rel=0.05)


def test_the_polyline_arcs_stay_within_one_percent_of_the_radius():
    course = avoid.plan_course([0, 0, 0], [100, 0, 0], [sphere("a", [50, 0, 0], 20.0)])
    assert clear(course, [sphere("a", [50, 0, 0], 20.0)])
    poly = np.array(course["polyline"])
    arc_points = poly[(np.abs(np.linalg.norm(poly - np.array([50, 0, 0]), axis=1) - 20.0) < 1e-3)]
    assert len(arc_points) >= 4
    for k in range(len(poly) - 1):
        mid = (poly[k] + poly[k + 1]) / 2
        if np.linalg.norm(poly[k] - [50, 0, 0]) < 20.001 and np.linalg.norm(poly[k + 1] - [50, 0, 0]) < 20.001:
            assert 20.0 - np.linalg.norm(mid - [50, 0, 0]) <= 0.0105 * 20.0 + 1e-9


def test_a_start_inside_a_sphere_is_exempt():
    spheres = [sphere("big", [0, 0, 0], 5.0), sphere("next", [30, 0.5, 0], 3.0)]
    course = avoid.plan_course([1, 0, 0], [60, 0, 0], spheres)
    assert course["avoided"] == ["next"]
    assert clear(course, [spheres[1]])


def test_an_end_on_a_surface_is_exempt():
    spheres = [sphere("moon", [10, 0, 0], 2.0)]
    course = avoid.plan_course([0, 0, 0], [12, 0, 0], spheres)
    assert course["avoided"] == [] and course["length"] == 12.0


def test_both_ends_in_one_sphere_is_the_straight_line():
    spheres = [sphere("big", [0, 0, 0], 50.0), sphere("small", [3, 0, 0], 1.0)]
    course = avoid.plan_course([-5, 0, 0], [5, 0, 0], spheres)
    assert course["length"] == 10.0 and course["avoided"] == [] and course["warnings"] == []


def test_the_containment_test_is_closed():
    on_edge = avoid.plan_course([5, 0, 0], [40, 0, 0], [sphere("a", [0, 0, 0], 5.0)])
    assert on_edge["avoided"] == []
    outside = avoid.plan_course([5 * (1 + 1e-6), 0, 0], [-40, 0, 0], [sphere("a", [0, 0, 0], 5.0)])
    assert outside["avoided"] == ["a"]


def test_a_start_equal_to_the_end():
    course = avoid.plan_course([3, 3, 3], [3, 3, 3], [sphere("a", [3, 3, 4], 2.0)])
    assert course["length"] == 0.0 and len(course["waypoints"]) == 2


def test_two_tangent_spheres_with_the_line_through_the_contact_go_straight():
    spheres = [sphere("a", [10, 1, 0], 1.0), sphere("b", [10, -1, 0], 1.0)]
    # the line y = 0 touches both only at the contact point; grazing is not a hit
    course = avoid.plan_course([0, 0, 0], [20, 0, 0], spheres)
    assert course["length"] == pytest.approx(20.0, rel=1e-12) or course["avoided"]
    assert clear(course, spheres)


def test_two_spheres_with_a_gap_go_through_the_gap_when_it_is_shorter():
    spheres = [sphere("a", [10, 1.005, 0], 1.0), sphere("b", [10, -1.005, 0], 1.0)]
    course = avoid.plan_course([0, 0, 0], [20, 0, 0], spheres)
    assert clear(course, spheres)
    assert course["length"] < 20.2


def test_a_chain_of_spheres_is_valid_and_each_is_avoided():
    spheres = [sphere("s%d" % k, [10 + 15 * k, 0.5 * (-1) ** k, 0], 3.0) for k in range(5)]
    course = avoid.plan_course([0, 0, 0], [100, 0, 0], spheres)
    assert clear(course, spheres)
    # a tangent line from one wrap to the next can pass a middle sphere without bending, so the wraps are some of the spheres
    assert course["avoided"] and set(course["avoided"]) <= {"s%d" % k for k in range(5)}
    assert course["length"] > 100.0


def test_forty_tiny_spheres_are_valid_and_independent_of_a_shift():
    rng = np.random.default_rng(5)
    xs = np.sort(rng.uniform(5, 95, 40))
    radii = 10 ** rng.uniform(-6, -2, 40)
    spheres = [sphere("t%d" % k, [xs[k], 0.0, 0.0], radii[k]) for k in range(40)]
    base = avoid.plan_course([0, 0, 0], [100, 0, 0], spheres)
    shift = np.array([2e4, -1.5e4, 3e3])
    moved = avoid.plan_course(shift, shift + [100, 0, 0], [sphere(s["ref"], np.array(s["center"]) + shift, s["radius"]) for s in spheres])
    assert base["length"] == pytest.approx(moved["length"], rel=1e-10)
    assert clear(base, spheres)
    # every sphere is smaller than 0.002 of the leg: stepped around, listed apart, not drawn
    assert base["avoided"] == [] and len(base["ignored_below_resolution"]) > 0
    assert len(base["polyline"]) == len(base["waypoints"])


def test_a_small_sphere_is_stepped_around_without_arc_samples():
    tiny = sphere("rogue_planet:9", [50, 0, 0], 0.05)
    course = avoid.plan_course([0, 0, 0], [100, 0, 0], [tiny])
    assert 0.05 < tuning.COURSE_BELOW_RESOLUTION * 100
    assert course["ignored_below_resolution"] == ["rogue_planet:9"] and course["avoided"] == []
    assert [w["kind"] for w in course["waypoints"]] == ["start", "tangent_in", "tangent_out", "end"]
    assert len(course["polyline"]) == 4
    assert course["length"] > 100.0


def test_an_overlapping_pair_is_merged_and_the_course_is_clear():
    spheres = [sphere("a", [50, 2, 0], 5.0), sphere("b", [50, -3, 0], 5.0)]
    course = avoid.plan_course([0, 0, 0], [100, 0, 0], spheres)
    # the enclosing sphere is wider than either, so its chords dip further into a member than one sphere's own sagitta allows
    assert clear(course, spheres, margin=0.03)
    assert course["warnings"] == []
    assert set(course["avoided"]) == {"a", "b"}
    assert course["waypoints"][1]["ref"] == ["a", "b"]


def test_an_overlap_that_swallows_an_end_is_the_flagged_straight_line():
    spheres = [sphere("a", [10, 0, 0], 8.0), sphere("b", [20, 0, 0], 8.0)]
    # start lies outside both, but the enclosing sphere of the pair holds it
    course = avoid.plan_course([3.5, 6.0, 0], [60, 0.1, 0], spheres)
    if course["warnings"]:
        assert course["warnings"] == ["overlap_fallback"] and course["length"] == course["straight_length"]
    else:
        assert clear(course, spheres)


def test_the_enclosing_sphere_holds_what_it_is_given():
    rng = np.random.default_rng(1)
    centers = rng.uniform(-10, 10, (6, 3))
    radii = rng.uniform(0.5, 3.0, 6)
    c, r = avoid.enclosing_sphere(centers, radii)
    assert np.all(np.linalg.norm(centers - c, axis=1) + radii <= r + 1e-9)
    pair_c, pair_r = avoid.enclosing_sphere(np.array([[0.0, 0, 0], [10.0, 0, 0]]), np.array([2.0, 3.0]))
    assert pair_r == pytest.approx(7.5) and pair_c == pytest.approx([5.5, 0, 0])
    inner_c, inner_r = avoid.enclosing_sphere(np.array([[0.0, 0, 0], [1.0, 0, 0]]), np.array([5.0, 1.0]))
    assert inner_r == 5.0


def test_a_random_field_of_disjoint_spheres_is_always_valid():
    rng = np.random.default_rng(11)
    for trial in range(12):
        spheres = []
        while len(spheres) < 40:
            c = rng.uniform(0, 100, 3)
            r = rng.uniform(2, 8)
            if all(np.linalg.norm(c - np.array(s["center"])) > r + s["radius"] + 0.5 for s in spheres):
                spheres.append(sphere("s%d" % len(spheres), c, r))
        a, b = rng.uniform(0, 100, 3), rng.uniform(0, 100, 3)
        if any(np.linalg.norm(a - s["center"]) <= s["radius"] * 1.01 or np.linalg.norm(b - s["center"]) <= s["radius"] * 1.01 for s in spheres):
            continue
        course = avoid.plan_course(a, b, spheres)
        assert course["warnings"] == []
        assert clear(course, spheres)
        assert course["length"] >= course["straight_length"] - 1e-9


def test_more_wraps_than_the_limit_is_a_flagged_straight_line(monkeypatch):
    monkeypatch.setattr(tuning, "COURSE_MAX_WRAPS", 1)
    spheres = [sphere("s%d" % k, [10 + 15 * k, 0, 0], 3.0) for k in range(5)]
    course = avoid.plan_course([0, 0, 0], [100, 0, 0], spheres)
    assert course["warnings"] == ["wrap_limit"] and course["length"] == 100.0


def test_passed_through_notes_are_carried_and_spheres_of_no_size_dropped():
    note = [{"ref": "nebula:55", "note": "thin gas"}]
    course = avoid.plan_course([0, 0, 0], [10, 0, 0], [sphere("z", [5, 0, 0], 0.0)], passed_through=note)
    assert course["passed_through"] == note and course["length"] == 10.0


def test_obstacles_from_the_database_become_spheres():
    obstacle = {"ref": "star:7", "center_pc": [1.0, 2.0, 3.0], "radius_pc": 0.4}
    assert avoid.spheres_from_obstacles([obstacle]) == [{"ref": "star:7", "center": [1.0, 2.0, 3.0], "radius": 0.4}]


def test_overlapping_random_fields_give_a_clear_course_or_a_flagged_line():
    rng = np.random.default_rng(23)
    flagged = clear_count = 0
    for trial in range(40):
        spheres = [sphere("s%d" % k, rng.uniform(0, 100, 3), rng.uniform(3, 12)) for k in range(30)]
        a, b = rng.uniform(0, 100, 3), rng.uniform(0, 100, 3)
        course = avoid.plan_course(a, b, spheres)
        if course["warnings"]:
            assert course["warnings"] == ["overlap_fallback"]
            assert course["length"] == course["straight_length"]
            flagged += 1
        else:
            # a merged cluster is drawn on its enclosing sphere, so the chords may dip further into a member
            assert clear(course, spheres, margin=0.2)
            clear_count += 1
    assert clear_count > 0
