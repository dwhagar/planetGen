"""
planetgen/web/maps/navmap.py regression tests.

Covers the pure geometry (`_scale`/`_project_all`) and the rendered SVG
panel's shape, with no database needed: `render_nav_map_panel` takes plain waypoint dicts,
the same shape `html/nav.py` builds from `queryDb.nav_between`'s result.

Run with: pytest src/tests/test_navmap.py
"""

import pytest  # noqa: E402

from planetgen.web.maps.navmap import _project_all, _scale, render_nav_map_panel  # noqa: E402


def _link(name, **params):
    """Stands in for `web.helpers.page_url`: `/<name>?k=v&...`."""
    return f"/{name}?" + "&".join(f"{key}={value}" for key, value in sorted(params.items()))


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
    html = render_nav_map_panel(_link, waypoints, has_route=False)

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
    html = render_nav_map_panel(_link, waypoints, has_route=True)

    assert "navmap-direct-line" in html
    assert "navmap-route-line" in html
    assert "navmap-hop" in html
    assert "Waystation" in html
    assert '<a class="navmap-point navmap-origin" href="/system?system_id=1">' in html
    assert '<a class="navmap-point navmap-hop" href="/system?system_id=2">' in html
    assert '<a class="navmap-point navmap-destination" href="/system?system_id=3">' in html
    assert "data-nav" not in html


def test_render_nav_map_panel_draws_a_jump_through_unknown_space_red_dashed_and_glowing():
    """NAV.36: a hop flagged as a jump through unknown space is its own path (red, dashed, glowing, a title in words), the known hops stay
    one solid path, and the panel and its help say what the red hop is; with no such hop none of it is drawn."""
    waypoints = [
        _waypoint(1, "Origin", (0.0, 0.0, 0.0), "origin"),
        {**_waypoint(2, "Waystation", (2.0, 0.0, 0.0), "hop"), "hop_unknown": True},
        _waypoint(3, "Far", (4.0, 1.0, 0.0), "hop"),
        _waypoint(4, "Destination", (6.0, 1.0, 0.0), "destination"),
    ]
    html = render_nav_map_panel(_link, waypoints, has_route=True)
    assert html.count("navmap-route-unknown") == 2 + 0          # the path, and the legend's sample line
    assert html.count('class="navmap-route-line navmap-route-unknown"') == 1
    assert 'filter="url(#navmap-glow)"' in html and '<filter id="navmap-glow"' in html
    assert "A jump through unknown space: no charted system lies on this hop" in html
    assert "A jump through unknown space: red, dashed and glowing." in html and "A red, dashed, glowing hop" in html
    assert html.count('class="navmap-route-line"') == 1          # Far to Destination is one known run
    known = render_nav_map_panel(_link, [{**w, "hop_unknown": False} for w in waypoints], has_route=True)
    assert "navmap-route-unknown" not in known and "navmap-glow" not in known and "navmap-key" not in known
    assert known.count('class="navmap-route-line"') == 1


def test_render_nav_map_panel_escapes_names():
    waypoints = [
        _waypoint(1, "<script>alert(1)</script>", (0.0, 0.0, 0.0), "origin"),
        _waypoint(2, "Destination", (1.0, 0.0, 0.0), "destination"),
    ]
    html = render_nav_map_panel(_link, waypoints, has_route=False)

    assert "<script>" not in html
    assert "&lt;script&gt;" in html


def test_render_nav_map_panel_links_a_phenomenon_endpoint():
    waypoints = [
        {"id": 4, "kind": "phenomenon", "type": "black_hole", "name": "Maw",
         "position": (0.0, 0.0, 0.0), "role": "origin"},
        _waypoint(2, "Destination", (1.0, 0.0, 0.0), "destination"),
    ]
    html = render_nav_map_panel(_link, waypoints, has_route=False)
    assert 'href="/phenomenon?phenomenon_id=4&amp;phenomenon_type=black_hole"' in html


def _compass_tip(html):
    import re
    match = re.search(r'class="navmap-compass" x1="([\d.]+)" y1="([\d.]+)" x2="([\d.]+)" y2="([\d.]+)"', html)
    x1, y1, x2, y2 = map(float, match.groups())
    return x2 - x1, y2 - y1


def test_compass_points_from_the_origin_toward_the_frame_center():
    # Origin out on +Y: bearing 000 is toward the center, -Y, which is
    # down on screen (SVG y grows downward).
    waypoints = [
        _waypoint(1, "Origin", (0.0, 10.0, 0.0), "origin"),
        _waypoint(2, "Destination", (4.0, 10.0, 0.0), "destination"),
    ]
    html = render_nav_map_panel(_link, waypoints, has_route=False)
    dx, dy = _compass_tip(html)
    assert dx == pytest.approx(0.0, abs=0.1)
    assert dy > 0
    assert "Bearing 000" in html


def test_compass_uses_the_given_frame_center_and_falls_back_to_plus_x():
    waypoints = [
        _waypoint(1, "Origin", (0.0, 0.0, 0.0), "origin"),
        _waypoint(2, "Destination", (4.0, 0.0, 0.0), "destination"),
    ]
    toward_minus_x = render_nav_map_panel(_link, waypoints, has_route=False, frame_center=(-5.0, 0.0, 0.0))
    dx, dy = _compass_tip(toward_minus_x)
    assert dx < 0 and dy == pytest.approx(0.0, abs=0.1)

    # Origin on the center itself: +X, as navigation.course_between does.
    dx, dy = _compass_tip(render_nav_map_panel(_link, waypoints, has_route=False))
    assert dx > 0 and dy == pytest.approx(0.0, abs=0.1)


def test_render_nav_map_panel_tooltips_name_the_course_to_the_next_stop():
    """NAV.42: a stop's tooltip carries its course to the next stop; the last stop has none."""
    first = _waypoint(1, "Alpha", (0.0, 0.0, 0.0), "origin")
    first["course"] = "045 mark 012"
    last = _waypoint(2, "Omega", (3.0, 4.0, 0.0), "destination")
    html = render_nav_map_panel(_link, [first, last], has_route=True)
    assert "<title>Alpha (next stop: 045 mark 012)</title>" in html
    assert "<title>Omega</title>" in html


def test_a_course_bent_round_a_sphere_is_drawn_over_the_route_and_the_key_says_so():  # NAV.28
    waypoints = [_waypoint(1, "Alpha", (0.0, 0.0, 0.0), "origin"), _waypoint(2, "Beta", (10.0, 0.0, 0.0), "destination")]
    bent = [[(0.0, 0.0, 0.0), (4.0, 1.5, 0.0), (6.0, 1.5, 0.0), (10.0, 0.0, 0.0)]]
    plain = render_nav_map_panel(_link, waypoints, has_route=True)
    html = render_nav_map_panel(_link, waypoints, has_route=True, bent=bent)
    assert "navmap-bent-line" not in plain and "bent round" not in plain
    assert html.count('<path class="navmap-bent-line"') == 2 and "The course bent round a keep-out sphere." in html
    assert "amber line" in html
    # one in the drawing, one in the key
    assert html.index("navmap-route-line") < html.index('<path class="navmap-bent-line"')       # over the route


def test_the_map_is_scaled_to_hold_a_bend_that_leaves_the_straight_line():
    waypoints = [_waypoint(1, "Alpha", (0.0, 0.0, 0.0), "origin"), _waypoint(2, "Beta", (10.0, 0.0, 0.0), "destination")]
    flat = render_nav_map_panel(_link, waypoints, has_route=False, bent=[[(0.0, 0.0, 0.0), (10.0, 0.0, 0.0)]])
    bulged = render_nav_map_panel(_link, waypoints, has_route=False, bent=[[(0.0, 0.0, 0.0), (5.0, 40.0, 0.0), (10.0, 0.0, 0.0)]])
    assert flat != bulged and "navmap-bent-line" in bulged
