"""
html/lib/galaxymap.py regression tests -- covers `_phenomenon_elements`/
`render_galaxy_map_panel`'s new `phenomena` parameter (schema v18's
galaxy-placed nebulae/asteroid fields, plotted as small fixed-size dots
alongside sector "star" dots -- see `starmap.py`'s translucent clouds for
where a phenomenon's own real physical size is depicted instead). Same
`sys.path` setup as `test_starmap.py`/`test_navmap.py`; no database
needed, since `render_galaxy_map_panel` takes plain dicts, the same shape
`queryDb.galaxy_placed_sectors`/`galaxy_placed_phenomena` return.

Run with: pytest src/tests/test_galaxymap.py
"""
import os
import sys

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_SRC_DIR, "html", "lib"))

from galaxymap import render_galaxy_map_panel  # noqa: E402


def _sector(id_=1, x=10.0, y=5.0, name="Test Sector"):
    return {
        "id": id_, "name": name, "x": x, "y": y, "z": 0.0,
        "galactic_radius_pc": (x ** 2 + y ** 2) ** 0.5, "shell_index": 1, "system_count": 3,
    }


def _phenomenon(id_=1, type_="nebula", descriptor="dark", x=20.0, y=-15.0, radius_ly=8.0):
    return {
        "id": id_, "type": type_, "name": "Test Phenomenon", "descriptor": descriptor,
        "radius_ly": radius_ly, "x": x, "y": y, "z": 0.0,
        "galactic_radius_pc": (x ** 2 + y ** 2) ** 0.5,
    }


def test_render_galaxy_map_panel_draws_a_nebula_dot():
    html = render_galaxy_map_panel("db", [_sector()], phenomena=[_phenomenon(type_="nebula")])
    assert 'class="galaxymap-phenomenon"' in html
    assert "Test Phenomenon" in html
    assert "Nebula" in html


def test_render_galaxy_map_panel_draws_an_asteroid_field_dot():
    html = render_galaxy_map_panel("db", [_sector()], phenomena=[_phenomenon(type_="asteroid_field")])
    assert 'class="galaxymap-phenomenon"' in html
    assert "Asteroid Field" in html


def test_render_galaxy_map_panel_without_phenomena_matches_omitting_the_argument():
    sectors = [_sector()]
    html_default = render_galaxy_map_panel("db", sectors)
    html_none = render_galaxy_map_panel("db", sectors, phenomena=None)
    html_empty = render_galaxy_map_panel("db", sectors, phenomena=[])
    assert "galaxymap-phenomenon" not in html_default
    assert html_default == html_none == html_empty


def test_render_galaxy_map_panel_handles_no_sectors_and_a_phenomenon():
    # Must not crash even when there's nothing else placed in the galaxy
    # yet -- a phenomenon can exist with no generated sectors around it.
    html = render_galaxy_map_panel("db", [], phenomena=[_phenomenon()])
    assert 'class="galaxymap-phenomenon"' in html
    assert "No sectors have been placed" in html


def test_phenomenon_dots_use_a_distinct_color_per_type():
    from galaxymap import _PHENOMENON_COLORS
    assert _PHENOMENON_COLORS["nebula"] != _PHENOMENON_COLORS["asteroid_field"]


def test_render_galaxy_map_panel_draws_a_black_hole_dot():
    html = render_galaxy_map_panel(
        "db", [_sector()], phenomena=[_phenomenon(type_="black_hole", descriptor="accreting", radius_ly=0)]
    )
    assert 'class="galaxymap-phenomenon"' in html
    assert "Accreting Black Hole" in html


def test_render_galaxy_map_panel_draws_a_neutron_star_dot():
    html = render_galaxy_map_panel(
        "db", [_sector()], phenomena=[_phenomenon(type_="neutron_star", descriptor="young", radius_ly=0)]
    )
    assert 'class="galaxymap-phenomenon"' in html
    assert "Young Neutron Star" in html


def test_black_hole_and_neutron_star_dots_omit_the_misleading_zero_radius_across_line():
    # A point-like object's real "size" (radius_ly=0) isn't worth showing
    # -- the tooltip should read "~N ly from core" without a "~0.0 ly
    # across" prefix a nebula/asteroid field's own tooltip does show.
    nebula_html = render_galaxy_map_panel("db", [], phenomena=[_phenomenon(type_="nebula", radius_ly=8.0)])
    bh_html = render_galaxy_map_panel(
        "db", [], phenomena=[_phenomenon(type_="black_hole", descriptor="quiescent", radius_ly=0)]
    )
    assert "ly across" in nebula_html
    assert "ly across" not in bh_html


def test_all_four_phenomenon_types_have_distinct_colors_and_labels():
    from galaxymap import _PHENOMENON_COLORS, _PHENOMENON_TYPE_LABELS
    types = ("nebula", "asteroid_field", "black_hole", "neutron_star")
    colors = {_PHENOMENON_COLORS[t] for t in types}
    labels = {_PHENOMENON_TYPE_LABELS[t] for t in types}
    assert len(colors) == 4
    assert len(labels) == 4


def _galaxy_shape():
    """A real, small `queryDb.galaxy_density_shape`-shaped dict -- same
    field set `html/galaxy.py` passes through from `GET /api/galaxy/shape`,
    built via the actual `galaxyDensity.build_galaxy_shape` rather than
    hand-picked numbers, so these tests exercise the real model."""
    from stellarObjects.galaxyDensity import build_galaxy_shape

    shape = build_galaxy_shape(
        disk_scale_length_pc=2800.0, disk_scale_height_pc=350.0,
        bulge_scale_radius_pc=200.0, bulge_amplitude=1.0,
        arm_count=2, pitch_angle_rad=0.2618, arm_amplitude=0.4,
    )
    shape_dict = dict(shape._asdict())
    shape_dict.update({"edge_pc": 3.5, "outer_shell_index": 4000, "expected_system_count_at_density_1": 4.32})
    return shape_dict


def test_render_galaxy_map_panel_shades_real_density_when_shape_given():
    html = render_galaxy_map_panel("db", [_sector()], galaxy_shape=_galaxy_shape())
    assert "galaxyDensityClip" in html
    assert "<rect" in html
    assert "galaxyCloud" not in html
    assert "real predicted stellar density" in html
    assert "density skeleton hasn't been built yet" not in html


def test_render_galaxy_map_panel_falls_back_to_illustrative_gradient_without_shape():
    html = render_galaxy_map_panel("db", [_sector()], galaxy_shape=None)
    assert "galaxyCloud" in html
    assert "galaxyDensityClip" not in html
    assert "illustrative expected density, not real data" in html
    assert "generate.py plan" in html


def test_render_galaxy_map_panel_zooms_out_to_show_the_spiral_when_little_is_placed():
    # A single sector close to the core would otherwise keep the default
    # view pinned to just a few hundred ly (_MIN_RINGS_SHOWN) -- deep
    # inside this shape's own "trivially solid" bulge core, nowhere near
    # its spiral structure. With a real galaxy_shape, the default view
    # should reach out toward the disk scale length instead.
    shape = _galaxy_shape()
    with_shape = render_galaxy_map_panel("db", [_sector()], galaxy_shape=shape)
    without_shape = render_galaxy_map_panel("db", [_sector()], galaxy_shape=None)

    def px_per_ly(html):
        marker = 'data-px-per-ly="'
        start = html.index(marker) + len(marker)
        return float(html[start:html.index('"', start)])

    assert px_per_ly(with_shape) < px_per_ly(without_shape)


def test_render_galaxy_map_panel_handles_a_real_shape_with_no_sectors():
    html = render_galaxy_map_panel("db", [], galaxy_shape=_galaxy_shape())
    assert "<rect" in html
    assert "No sectors have been placed" in html
