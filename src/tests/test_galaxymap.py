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
