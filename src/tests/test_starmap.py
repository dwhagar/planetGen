"""
html/lib/starmap.py regression tests.

Covers `_rotate_to_galaxy_frame` -- the fix for the Sector Map's star dots
being plotted as if their own sector-local (x, y, z) axes already ran
parallel to the galaxy frame's, when the wedge outline/"Galactic Center"
compass arrow were always correctly expressed in the galaxy frame
directly (see docs/TODO.md's now-resolved "Sector Map's on-shell wedge
shape and 'Galactic Center' compass arrow ... assume a sector's own local
(x, y, z) axes run parallel to the galaxy frame's axes" entry) -- via
`stellarObjects.sectorGeometry.cube_orientation`'s already-established
"Cube orientation" convention (radial-outward local +Z, projected-
galactic-north local +X), applied at render time from the sector's own
stored `center_x/y/z_pc` rather than needing any new stored orientation.

Same `sys.path` setup as `test_navmap.py` (`html/lib` isn't part of the
installed `stellarObjects` package, CGI-only plumbing); no database
needed, since `render_map_panel`/`_rotate_to_galaxy_frame` take plain
dicts/tuples.

Run with: pytest src/tests/test_starmap.py
"""
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


def _dot_anchor_style(html):
    match = re.search(
        r'star-dot-anchor" style="left:([\-0-9.]+)px; top:([\-0-9.]+)px.*?translateZ\(([\-0-9.]+)px\)',
        html,
    )
    assert match, f"no star-dot-anchor found in:\n{html}"
    return tuple(float(g) for g in match.groups())


def test_render_map_panel_rotates_star_dots_only_when_placed():
    """
    End-to-end: the same system's dot must land in a different on-screen
    position depending on whether the sector has an (off-axis) galaxy
    placement, since a real placement now rotates it into the galaxy
    frame first -- and must land back at the original, unrotated position
    once the placement is removed again (`center_pc=None`), confirming the
    unplaced fallback is exactly the pre-fix behavior, not a regression.
    """
    system = _make_system()

    html_unplaced = render_map_panel("db", 1000.0, None, None, None, [system])
    unplaced_px = _dot_anchor_style(html_unplaced)

    html_placed = render_map_panel("db", 1000.0, 3, 42, (500.0, 200.0, -100.0), [system])
    placed_px = _dot_anchor_style(html_placed)

    assert placed_px != unplaced_px

    html_unplaced_again = render_map_panel("db", 1000.0, None, None, None, [system])
    assert _dot_anchor_style(html_unplaced_again) == unplaced_px


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
    html_unplaced = render_map_panel("db", 1000.0, None, None, None, [system])
    html_placed_on_axis = render_map_panel("db", 1000.0, 3, 42, (0.0, 0.0, 999.0), [system])
    # Only the dot position should match; the placed scene also draws a
    # wedge outline/compass arrow the unplaced one doesn't, so compare
    # just the anchors, not the whole HTML.
    assert _dot_anchor_style(html_placed_on_axis) == _dot_anchor_style(html_unplaced)


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
    # into an unbounded DOM element size -- see `_MAX_CLOUD_RADIUS_PX`.
    from starmap import _MAX_CLOUD_RADIUS_PX
    huge = _phenomenon_cloud_radius_px(200.0, half_edge=100.0)
    assert huge <= _MAX_CLOUD_RADIUS_PX


def test_render_map_panel_draws_a_nebula_cloud_with_its_own_data_attrs():
    system = _make_system()
    phenomenon = _phenomenon(type_="nebula", descriptor="reflection")
    html = render_map_panel("db", 1000.0, None, None, None, [system], phenomena=[phenomenon])
    assert 'class="phenomenon-cloud nebula-cloud billboard"' in html
    assert 'data-kind="phenomenon"' in html
    assert 'data-phenomenon-type="Reflection Nebula"' in html
    assert "Click a star system or cloud for details." in html


def test_render_map_panel_draws_an_asteroid_field_cloud():
    system = _make_system()
    phenomenon = _phenomenon(type_="asteroid_field", descriptor="dense")
    html = render_map_panel("db", 1000.0, None, None, None, [system], phenomena=[phenomenon])
    assert 'class="phenomenon-cloud asteroid-cloud billboard"' in html
    assert 'data-phenomenon-type="Asteroid Field (Dense)"' in html


def test_render_map_panel_without_phenomena_matches_omitting_the_argument():
    system = _make_system()
    html_default = render_map_panel("db", 1000.0, None, None, None, [system])
    html_explicit_empty = render_map_panel("db", 1000.0, None, None, None, [system], phenomena=[])
    html_none = render_map_panel("db", 1000.0, None, None, None, [system], phenomena=None)
    assert "phenomenon-cloud" not in html_default
    assert html_default == html_explicit_empty == html_none


def test_render_map_panel_places_a_phenomenon_directly_in_the_galaxy_frame_unrotated():
    # A phenomenon's offset_*_ly is already galaxy-frame (see
    # queryDb.phenomena_near_sector) -- unlike a star system's sector-local
    # x/y/z, it must NOT be re-rotated by `_rotate_to_galaxy_frame` even
    # when the sector itself has a galaxy placement.
    phenomenon = _phenomenon(offset=(3.0, 0.0, 0.0))
    html_unplaced = render_map_panel("db", 1000.0, None, None, None, [], phenomena=[phenomenon])
    html_placed = render_map_panel("db", 1000.0, 3, 42, (500.0, 200.0, -100.0), [], phenomena=[phenomenon])

    def _cloud_anchor_style(html):
        match = re.search(
            r'cloud-anchor" style="left:([\-0-9.]+)px; top:([\-0-9.]+)px.*?translateZ\(([\-0-9.]+)px\)',
            html,
        )
        assert match, f"no cloud-anchor found in:\n{html}"
        return tuple(float(g) for g in match.groups())

    assert _cloud_anchor_style(html_unplaced) == _cloud_anchor_style(html_placed)
