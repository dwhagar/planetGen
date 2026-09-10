"""
html/lib/navmap.py regression tests.

Covers the pure geometry (`_scale`/`_project_all`) and the rendered SVG
panel's shape -- same `sys.path` setup as `test_mdconvert.py` (`html/lib`
isn't part of the installed `stellarObjects` package, CGI-only plumbing),
and no database needed: `render_nav_map_panel` takes plain waypoint dicts,
the same shape `html/nav.py` builds from `queryDb.nav_between`'s result.

Run with: pytest src/tests/test_navmap.py
"""
import os
import sys

# This file lives at src/tests/ (src layout) -- one dirname() call reaches
# src/, then down into html/lib (src/html/lib, not a top-level html/).
_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_SRC_DIR, "html", "lib"))

import pytest  # noqa: E402

from navmap import _project_all, _scale, render_nav_map_panel  # noqa: E402


def _waypoint(id_, name, position, role):
    return {"id": id_, "name": name, "position": position, "role": role}


def test_scale_uses_the_larger_axis_for_a_uniform_ratio():
    # A wide, short bounding box: x spans 10, y spans 0 -- the ratio must
    # come from the x span so a uniform scale doesn't stretch y.
    px_per_ly, center_x, center_y = _scale([(0.0, 0.0), (10.0, 0.0)])
    assert center_x == pytest.approx(5.0)
    assert center_y == pytest.approx(0.0)
    assert px_per_ly > 0


def test_scale_floors_span_for_coincident_or_near_points():
    # Two effectively-identical points must not divide by (near) zero.
    px_per_ly, _cx, _cy = _scale([(1.0, 1.0), (1.0, 1.0)])
    assert px_per_ly > 0


def test_project_all_flips_y_for_screen_space():
    # +Y is "up" in this project's galactic-plane convention; SVG grows
    # downward, so a positive-y point must land above center on screen
    # (a smaller svg_y).
    px_per_ly, center_x, center_y = _scale([(0.0, 0.0), (0.0, 10.0)])
    projected = _project_all([(0.0, 0.0), (0.0, 10.0)], px_per_ly, center_x, center_y)
    (_x0, y0), (_x1, y1) = projected
    assert y1 < y0


def test_render_nav_map_panel_direct_only_has_no_route_line():
    waypoints = [
        _waypoint(1, "Origin", (0.0, 0.0, 0.0), "origin"),
        _waypoint(2, "Destination", (4.0, 0.0, 0.0), "destination"),
    ]
    html = render_nav_map_panel("test.db", waypoints, has_route=False)

    assert "navmap-direct-line" in html
    assert "navmap-route-line" not in html
    assert "Origin" in html
    assert "Destination" in html
    assert "no route via adjacent systems was found" in html


def test_render_nav_map_panel_with_route_draws_route_line_and_hops():
    waypoints = [
        _waypoint(1, "Origin", (0.0, 0.0, 0.0), "origin"),
        _waypoint(2, "Waystation", (2.0, 0.0, 0.0), "hop"),
        _waypoint(3, "Destination", (4.0, 0.0, 0.0), "destination"),
    ]
    html = render_nav_map_panel("test.db", waypoints, has_route=True)

    assert "navmap-direct-line" in html
    assert "navmap-route-line" in html
    assert "navmap-hop" in html
    assert "Waystation" in html
    assert 'href="system.py?db=test.db&amp;id=1"' in html
    assert 'href="system.py?db=test.db&amp;id=3"' in html


def test_render_nav_map_panel_escapes_names():
    waypoints = [
        _waypoint(1, "<script>alert(1)</script>", (0.0, 0.0, 0.0), "origin"),
        _waypoint(2, "Destination", (1.0, 0.0, 0.0), "destination"),
    ]
    html = render_nav_map_panel("test.db", waypoints, has_route=False)

    assert "<script>" not in html
    assert "&lt;script&gt;" in html
