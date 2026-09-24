"""
html/lib/starmap.py regression tests.

Covers `_rotate_to_galaxy_frame` -- the fix for the Sector Map's star dots
being plotted as if their own sector-local (x, y, z) axes already ran
parallel to the galaxy frame's, when the wedge outline/"Galactic Center"
compass arrow were always correctly expressed in the galaxy frame
directly -- via `stellarObjects.sectorGeometry.cube_orientation`'s already-
established "Cube orientation" convention (radial-outward local +Z,
projected-galactic-north local +X), applied at render time from the
sector's own stored `center_x/y/z_pc` rather than needing any new stored
orientation.

Also covers `render_map_panel`'s scene-data contract -- since the Sector
Map moved from server-rendered CSS `<div>`s to a `<script
type="application/json">` block a WebGL client reads (see starmap.py's own
module docstring), these tests parse that JSON instead of regex-matching
HTML, but check the same underlying behavior the old assertions did
(rotation agreement, phenomenon labels/kinds, radius scaling, click hint
text).

Same `sys.path` setup as `test_navmap.py` (`html/lib` isn't part of the
installed `stellarObjects` package, CGI-only plumbing); no database
needed, since `render_map_panel`/`_rotate_to_galaxy_frame` take plain
dicts/tuples.

Run with: pytest src/tests/test_starmap.py
"""
import json
import math
import os
import re
import sys

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_SRC_DIR, "html", "lib"))

import pytest  # noqa: E402

from starmap import _rotate_to_galaxy_frame, render_map_panel  # noqa: E402
from starmap import _phenomenon_cloud_radius_px  # noqa: E402


def _vec_norm(v):
    return math.sqrt(sum(c * c for c in v))


def _make_system(x=100.0, y=50.0, z=-30.0):
    return {
        "id": 1, "name": "Test System", "quadrant": "Q1", "location": "loc",
        "x": x, "y": y, "z": z,
        "stars": [{
            "star_type": "G2V", "temperature_k": 5772, "radius_km": 696000,
            "luminosity_w": 3.828e26, "temp_display": "5772 K",
        }],
    }


def test_no_placement_leaves_the_local_vector_unchanged():
    # No galaxy placement means no galaxy frame to rotate into at all --
    # the map's own long-standing fallback behavior (plain, unrotated cube).
    local_vec = (10.0, -5.0, 3.0)
    assert _rotate_to_galaxy_frame(None, local_vec) == local_vec
    assert _rotate_to_galaxy_frame((None, None, None), local_vec) == local_vec


@pytest.mark.parametrize("center_pc", [
    (500.0, -200.0, 800.0),
    (500.0, 200.0, -100.0),
    (-1234.5, 67.8, 90.1),
])
def test_rotation_preserves_vector_length(center_pc):
    # A rotation must not stretch or shrink -- only reorient -- the vector.
    local_vec = (10.0, -5.0, 3.0)
    rotated = _rotate_to_galaxy_frame(center_pc, local_vec)
    assert _vec_norm(rotated) == pytest.approx(_vec_norm(local_vec))


def test_on_axis_degeneracy_does_not_crash():
    # center_pc exactly along galactic +Z is the one direction where
    # "project galactic north onto the plane perpendicular to radial"
    # degenerates (see cube_orientation's own fallback) -- must still
    # produce a valid, length-preserving rotation, not raise or return
    # something degenerate.
    local_vec = (10.0, -5.0, 3.0)
    rotated = _rotate_to_galaxy_frame((0.0, 0.0, 999.0), local_vec)
    assert _vec_norm(rotated) == pytest.approx(_vec_norm(local_vec))


def _scene_data(html):
    """Extracts and parses `render_map_panel`'s embedded `#starmap-data`
    JSON payload -- the client-side scene-data contract these tests check
    against instead of the old version's server-rendered HTML `<div>`s."""
    match = re.search(
        r'<script type="application/json" id="starmap-data">(.*?)</script>',
        html, re.DOTALL,
    )
    assert match, f"no #starmap-data script found in:\n{html}"
    return json.loads(match.group(1))


def _star_position(scene, index=0):
    star = scene["stars"][index]
    return (star["x"], star["y"], star["z"])


def test_render_map_panel_rotates_star_dots_only_when_placed():
    """
    End-to-end: the same system's dot must land in a different scene
    position depending on whether the sector has an (off-axis) galaxy
    placement, since a real placement now rotates it into the galaxy
    frame first -- and must land back at the original, unrotated position
    once the placement is removed again (`center_pc=None`), confirming the
    unplaced fallback is exactly the pre-fix behavior, not a regression.
    """
    system = _make_system()

    scene_unplaced = _scene_data(render_map_panel("db", 1000.0, None, None, None, [system]))
    unplaced_pos = _star_position(scene_unplaced)

    scene_placed = _scene_data(render_map_panel("db", 1000.0, 3, 42, (500.0, 200.0, -100.0), [system]))
    placed_pos = _star_position(scene_placed)

    assert placed_pos != unplaced_pos

    scene_unplaced_again = _scene_data(render_map_panel("db", 1000.0, None, None, None, [system]))
    assert _star_position(scene_unplaced_again) == unplaced_pos


def test_render_map_panel_placed_on_axis_matches_unplaced():
    """
    A galaxy placement exactly along galactic +Z is the one direction
    `cube_orientation` treats as its own local +Z too (see
    `sectorGeometry.cube_orientation`'s degeneracy fallback), so rotating
    into the galaxy frame there is the identity transform -- confirms the
    rotation is doing real work in the general (off-axis) case above, not
    coincidentally matching for an unrelated reason.
    """
    system = _make_system()
    scene_unplaced = _scene_data(render_map_panel("db", 1000.0, None, None, None, [system]))
    scene_placed_on_axis = _scene_data(render_map_panel("db", 1000.0, 3, 42, (0.0, 0.0, 999.0), [system]))
    # Only the star's own position should match; the placed scene also
    # carries a wedge outline/compass arrow the unplaced one doesn't.
    assert _star_position(scene_placed_on_axis) == _star_position(scene_unplaced)


def test_render_map_panel_outline_is_a_wedge_when_placed_and_a_cube_otherwise():
    system = _make_system()
    scene_unplaced = _scene_data(render_map_panel("db", 1000.0, None, None, None, [system]))
    scene_placed = _scene_data(render_map_panel("db", 1000.0, 3, 42, (500.0, 200.0, -100.0), [system]))
    assert scene_unplaced["outline"]["kind"] == "cube"
    assert scene_placed["outline"]["kind"] == "wedge"
    # 12 edges, each a pair of 3D points, for either shape.
    for scene in (scene_unplaced, scene_placed):
        assert len(scene["outline"]["edges"]) == 12
        for edge in scene["outline"]["edges"]:
            assert len(edge) == 2
            assert all(len(point) == 3 for point in edge)


def test_render_map_panel_compass_present_only_when_placed():
    system = _make_system()
    scene_unplaced = _scene_data(render_map_panel("db", 1000.0, None, None, None, [system]))
    scene_placed = _scene_data(render_map_panel("db", 1000.0, 3, 42, (500.0, 200.0, -100.0), [system]))
    assert scene_unplaced["compass"] is None
    assert scene_placed["compass"] is not None
    assert scene_placed["compass"]["label"] == "N"


def test_render_map_panel_binary_system_gets_two_star_entries():
    system = _make_system()
    system["stars"].append({
        "star_type": "M4V", "temperature_k": 3200, "radius_km": 200000,
        "luminosity_w": 1.0e24, "temp_display": "3200 K",
    })
    scene = _scene_data(render_map_panel("db", 1000.0, None, None, None, [system]))
    assert len(scene["stars"]) == 2
    assert scene["stars"][0]["name"] == "Test System A"
    assert scene["stars"][1]["name"] == "Test System B"
    # The secondary is offset from (not stacked exactly on) the primary.
    assert (scene["stars"][1]["x"], scene["stars"][1]["y"]) != (scene["stars"][0]["x"], scene["stars"][0]["y"])


# --- nebula/asteroid-field clouds (schema v18) ------------------------------

def _phenomenon(type_="nebula", descriptor="emission", radius_ly=10.0, offset=(1.0, 2.0, -0.5), distance_ly=2.3):
    return {
        "id": 1, "type": type_, "name": "Test Cloud", "descriptor": descriptor, "radius_ly": radius_ly,
        "offset_x_ly": offset[0], "offset_y_ly": offset[1], "offset_z_ly": offset[2], "distance_ly": distance_ly,
    }


def test_phenomenon_cloud_radius_grows_with_radius_ly():
    small = _phenomenon_cloud_radius_px(1.0, half_edge=5000.0)
    large = _phenomenon_cloud_radius_px(100.0, half_edge=5000.0)
    assert 0 < small < large


def test_phenomenon_cloud_radius_is_capped_for_a_nebula_far_larger_than_the_sector():
    # A 200 ly emission nebula next to a small sector must not blow out
    # into an unbounded sprite size -- see `_MAX_CLOUD_RADIUS_PX`.
    from starmap import _MAX_CLOUD_RADIUS_PX
    huge = _phenomenon_cloud_radius_px(200.0, half_edge=100.0)
    assert huge <= _MAX_CLOUD_RADIUS_PX


def test_render_map_panel_draws_a_nebula_cloud_with_its_own_kind_and_colors():
    system = _make_system()
    phenomenon = _phenomenon(type_="nebula", descriptor="reflection")
    html = render_map_panel("db", 1000.0, None, None, None, [system], phenomena=[phenomenon])
    scene = _scene_data(html)
    assert len(scene["clouds"]) == 1
    cloud = scene["clouds"][0]
    assert cloud["kind"] == "nebula"
    assert cloud["typeLabel"] == "Reflection Nebula"
    assert cloud["coreColor"].startswith("#6fa8ff")
    assert "Click a star system or cloud for details." in html


def test_render_map_panel_draws_an_asteroid_field_cloud():
    system = _make_system()
    phenomenon = _phenomenon(type_="asteroid_field", descriptor="dense")
    scene = _scene_data(render_map_panel("db", 1000.0, None, None, None, [system], phenomena=[phenomenon]))
    cloud = scene["clouds"][0]
    assert cloud["kind"] == "asteroidField"
    assert cloud["typeLabel"] == "Asteroid Field (Dense)"
    assert "coreColor" not in cloud


def test_render_map_panel_without_phenomena_matches_omitting_the_argument():
    system = _make_system()
    html_default = render_map_panel("db", 1000.0, None, None, None, [system])
    html_explicit_empty = render_map_panel("db", 1000.0, None, None, None, [system], phenomena=[])
    html_none = render_map_panel("db", 1000.0, None, None, None, [system], phenomena=None)
    assert _scene_data(html_default)["clouds"] == []
    assert html_default == html_explicit_empty == html_none


def test_render_map_panel_places_a_phenomenon_directly_in_the_galaxy_frame_unrotated():
    # A phenomenon's offset_*_ly is already galaxy-frame (see
    # queryDb.phenomena_near_sector) -- unlike a star system's sector-local
    # x/y/z, it must NOT be re-rotated by `_rotate_to_galaxy_frame` even
    # when the sector itself has a galaxy placement.
    phenomenon = _phenomenon(offset=(3.0, 0.0, 0.0))
    scene_unplaced = _scene_data(render_map_panel("db", 1000.0, None, None, None, [], phenomena=[phenomenon]))
    scene_placed = _scene_data(render_map_panel("db", 1000.0, 3, 42, (500.0, 200.0, -100.0), [], phenomena=[phenomenon]))

    def _cloud_position(scene):
        cloud = scene["clouds"][0]
        return (cloud["x"], cloud["y"], cloud["z"])

    assert _cloud_position(scene_unplaced) == _cloud_position(scene_placed)


def test_phenomenon_cloud_radius_floors_at_minimum_for_a_point_like_object():
    # black_hole/neutron_star always query radius_ly as a literal 0 (see
    # queryDb._PHENOMENON_TABLES) -- must still floor at a visible minimum,
    # not collapse to an invisible 0px marker.
    from starmap import _MAX_CLOUD_RADIUS_PX
    radius = _phenomenon_cloud_radius_px(0.0, half_edge=5000.0)
    assert 0 < radius <= _MAX_CLOUD_RADIUS_PX


def test_render_map_panel_draws_an_accreting_black_hole_point():
    system = _make_system()
    phenomenon = _phenomenon(type_="black_hole", descriptor="accreting", radius_ly=0)
    scene = _scene_data(render_map_panel("db", 1000.0, None, None, None, [system], phenomena=[phenomenon]))
    cloud = scene["clouds"][0]
    assert cloud["kind"] == "blackHoleAccreting"
    assert cloud["typeLabel"] == "Black Hole (Accreting)"
    assert cloud["radiusText"] == "0.00 ly"


def test_render_map_panel_draws_a_quiescent_black_hole_point():
    system = _make_system()
    phenomenon = _phenomenon(type_="black_hole", descriptor="quiescent", radius_ly=0)
    scene = _scene_data(render_map_panel("db", 1000.0, None, None, None, [system], phenomena=[phenomenon]))
    cloud = scene["clouds"][0]
    assert cloud["kind"] == "blackHoleQuiescent"
    assert cloud["typeLabel"] == "Black Hole (Quiescent)"


def test_render_map_panel_draws_a_neutron_star_point():
    system = _make_system()
    phenomenon = _phenomenon(type_="neutron_star", descriptor="millisecond", radius_ly=0)
    scene = _scene_data(render_map_panel("db", 1000.0, None, None, None, [system], phenomena=[phenomenon]))
    cloud = scene["clouds"][0]
    assert cloud["kind"] == "neutronStar"
    assert cloud["typeLabel"] == "Neutron Star (Millisecond)"


def test_render_map_panel_neutron_star_with_no_descriptor_falls_back_to_plain_label():
    system = _make_system()
    phenomenon = _phenomenon(type_="neutron_star", descriptor=None, radius_ly=0)
    scene = _scene_data(render_map_panel("db", 1000.0, None, None, None, [system], phenomena=[phenomenon]))
    assert scene["clouds"][0]["typeLabel"] == "Neutron Star"


def test_json_script_escapes_script_close_tag_in_a_name():
    # A system name is arbitrary user-supplied text (see `--name`) -- one
    # containing "</script>" must not be able to break out of the
    # embedded JSON block (see starmap.py's own `_json_script`).
    system = _make_system()
    system["name"] = "Evil</script><script>alert(1)</script>"
    html = render_map_panel("db", 1000.0, None, None, None, [system])
    assert "</script><script>alert" not in html
    scene = _scene_data(html)
    assert scene["stars"][0]["name"].startswith("Evil</script>")
