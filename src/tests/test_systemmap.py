"""
html/lib/systemmap.py regression tests.

Covers the true-position layout engine that replaced the old schematic
"every planet due east of its star" System Map: the shared log-radial
scale (`_radial_scale_bounds`/`_radial_px`), real-angle placement
(`_polar_to_px`), the mass-weighted binary barycenter split
(`_binary_star_positions_km`), marker de-overlap (`_relax_markers`), 2D
label collision avoidance (`_label_sides_2d`), and the full rendered panel
for single-star/'close'-binary/'wide'-binary/empty systems. Same `sys.path`
setup as `test_starmap.py`/`test_navmap.py` (`html/lib` isn't part of the
installed `stellarObjects` package, CGI-only plumbing); no database
needed, since `render_system_map_panel` takes plain dicts, the same shape
`queryDb.system_detail` returns.

Run with: pytest src/tests/test_systemmap.py
"""
import math
import os
import sys

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_SRC_DIR, "html", "lib"))

import pytest  # noqa: E402

import systemmap as sm  # noqa: E402

AU_KM = 1.496e8


def _star(id_, mass_kg=1.989e30, radius_km=696_000, star_type="G2V", temperature_k=5778.0, luminosity_w=3.828e26):
    return {
        "id": id_, "mass_kg": mass_kg, "radius_km": radius_km, "star_type": star_type,
        "temperature_k": temperature_k, "luminosity_w": luminosity_w,
    }


def _planet(id_, star_id, x_km, y_km, distance_km=None, radius_km=6371.0, planet_class="G", body_type="t",
            moons=None, life_chemical=None, zone="e", period_years=1.0, gravity_g=1.0,
            atmosphere=None, composition=None, surface_temperature_k=None):
    return {
        "id": id_, "star_id": star_id, "name": f"Planet {id_}",
        "position_x_km": x_km, "position_y_km": y_km, "position_z_km": 0.0,
        "distance_km": distance_km if distance_km is not None else math.hypot(x_km, y_km),
        "radius_km": radius_km, "planet_class": planet_class, "body_type": body_type,
        "zone": zone, "period_years": period_years, "gravity_g": gravity_g,
        "life_chemical": life_chemical, "moons": moons or [],
        "atmosphere": atmosphere, "composition": composition, "surface_temperature_k": surface_temperature_k,
    }


def _moon(id_, x_km, y_km, radius_km=1737.0):
    return {
        "id": id_, "name": f"Moon {id_}", "position_x_km": x_km, "position_y_km": y_km, "position_z_km": 0.0,
        "distance_km": math.hypot(x_km, y_km), "radius_km": radius_km, "planet_class": "D", "body_type": "t",
        "zone": "e", "period_years": 0.05, "gravity_g": 0.1, "life_chemical": None,
    }


def _belt(id_, star_id, distance_km, lower_km, upper_km, density="typical"):
    return {
        "id": id_, "star_id": star_id, "distance_km": distance_km,
        "lower_limit_km": lower_km, "upper_limit_km": upper_km,
        "density": density, "composition_summary": "rock and ice",
    }


# --- _radial_scale_bounds / _radial_px ------------------------------------

def test_radial_scale_bounds_spans_from_below_the_smallest_to_the_largest():
    lo, hi = sm._radial_scale_bounds([1.0 * AU_KM, 10.0 * AU_KM, 100.0 * AU_KM])
    assert lo < 1.0 * AU_KM
    assert hi == pytest.approx(100.0 * AU_KM)


def test_radial_scale_bounds_handles_a_single_distance():
    lo, hi = sm._radial_scale_bounds([5.0 * AU_KM])
    assert 0 < lo < 5.0 * AU_KM <= hi


def test_radial_scale_bounds_handles_no_distances_without_crashing():
    lo, hi = sm._radial_scale_bounds([])
    assert lo == hi  # nothing to spread across


def test_radial_px_is_monotonic_in_distance():
    lo, hi = sm._radial_scale_bounds([1.0 * AU_KM, 50.0 * AU_KM])
    r1 = sm._radial_px(1.0 * AU_KM, lo, hi)
    r2 = sm._radial_px(10.0 * AU_KM, lo, hi)
    r3 = sm._radial_px(50.0 * AU_KM, lo, hi)
    assert r1 < r2 < r3


def test_radial_px_is_zero_for_a_body_on_its_own_anchor():
    lo, hi = sm._radial_scale_bounds([1.0 * AU_KM, 50.0 * AU_KM])
    assert sm._radial_px(0.0, lo, hi) == 0.0
    assert sm._radial_px(None, lo, hi) == 0.0


def test_radial_px_stays_within_the_configured_pixel_budget():
    lo, hi = sm._radial_scale_bounds([0.05 * AU_KM, 500.0 * AU_KM])
    for distance_au in (0.05, 1.0, 30.0, 500.0):
        r_px = sm._radial_px(distance_au * AU_KM, lo, hi)
        assert 0.0 <= r_px <= sm._MIN_RADIUS_PX + sm._RADIUS_SPREAD_PX + 1e-6


# --- _polar_to_px ----------------------------------------------------------

def test_polar_to_px_places_a_body_at_its_real_angle():
    # Due "east" (+x, y=0) -> straight right of the anchor, same height.
    x, y = sm._polar_to_px(100.0, 100.0, 50.0, 1.0, 0.0)
    assert x == pytest.approx(150.0)
    assert y == pytest.approx(100.0)

    # Due "north" (+y) -> up the screen, i.e. a *smaller* pixel y (SVG's
    # own y axis grows downward -- the one sign flip this module documents).
    x, y = sm._polar_to_px(100.0, 100.0, 50.0, 0.0, 1.0)
    assert x == pytest.approx(100.0)
    assert y == pytest.approx(50.0)


def test_polar_to_px_returns_the_anchor_itself_for_zero_radius():
    assert sm._polar_to_px(42.0, 7.0, 0.0, 1.0, 1.0) == (42.0, 7.0)


# --- _binary_star_positions_km ---------------------------------------------

def test_binary_star_positions_are_a_true_mass_weighted_barycenter():
    primary = _star(1, mass_kg=2.0e30)
    secondary = _star(2, mass_kg=1.0e30)
    system = {
        "binary_mutual_position_x_km": 3.0e5, "binary_mutual_position_y_km": 4.0e5,
        "binary_mutual_position_z_km": 0.0,
    }
    positions = sm._binary_star_positions_km(system, [primary, secondary])

    px, py = positions[1]
    sx, sy = positions[2]
    # secondary - primary must reproduce the stored separation exactly.
    assert (sx - px) == pytest.approx(3.0e5)
    assert (sy - py) == pytest.approx(4.0e5)
    # The heavier star sits closer to the shared barycenter (weighted by
    # the *other* star's mass): m1*p1 + m2*p2 == 0 around it.
    assert 2.0e30 * px + 1.0e30 * sx == pytest.approx(0.0, abs=1e-3)
    assert 2.0e30 * py + 1.0e30 * sy == pytest.approx(0.0, abs=1e-3)
    assert math.hypot(px, py) < math.hypot(sx, sy)


def test_binary_star_positions_falls_back_to_an_even_split_without_mass():
    system = {"binary_mutual_position_x_km": 10.0, "binary_mutual_position_y_km": 0.0,
              "binary_mutual_position_z_km": 0.0}
    positions = sm._binary_star_positions_km(system, [_star(1, mass_kg=0.0), _star(2, mass_kg=0.0)])
    assert positions[1] == pytest.approx((-5.0, 0.0))
    assert positions[2] == pytest.approx((5.0, 0.0))


# --- _relax_markers ---------------------------------------------------------

def test_relax_markers_separates_two_overlapping_markers():
    markers = [
        {"cx": 100.0, "cy": 100.0, "r": 20.0},
        {"cx": 105.0, "cy": 100.0, "r": 20.0},
    ]
    sm._relax_markers(markers, min_gap=5.0)
    dist = math.hypot(markers[1]["cx"] - markers[0]["cx"], markers[1]["cy"] - markers[0]["cy"])
    assert dist >= 20.0 + 20.0 + 5.0 - 1e-6


def test_relax_markers_leaves_already_separated_markers_untouched():
    markers = [{"cx": 0.0, "cy": 0.0, "r": 5.0}, {"cx": 100.0, "cy": 0.0, "r": 5.0}]
    original = [dict(m) for m in markers]
    sm._relax_markers(markers, min_gap=5.0)
    assert markers == original


def test_relax_markers_never_moves_a_fixed_marker():
    markers = [
        {"cx": 100.0, "cy": 100.0, "r": 30.0, "fixed": True},
        {"cx": 105.0, "cy": 100.0, "r": 10.0},
    ]
    sm._relax_markers(markers, min_gap=2.0)
    assert markers[0]["cx"] == 100.0 and markers[0]["cy"] == 100.0
    dist = math.hypot(markers[1]["cx"] - 100.0, markers[1]["cy"] - 100.0)
    assert dist >= 30.0 + 10.0 + 2.0 - 1e-6


def test_relax_markers_handles_identical_positions_without_crashing():
    markers = [{"cx": 50.0, "cy": 50.0, "r": 10.0}, {"cx": 50.0, "cy": 50.0, "r": 10.0}]
    sm._relax_markers(markers, min_gap=4.0)
    dist = math.hypot(markers[1]["cx"] - markers[0]["cx"], markers[1]["cy"] - markers[0]["cy"])
    assert dist >= 10.0 + 10.0 + 4.0 - 1e-6


# --- _label_sides_2d ---------------------------------------------------------

def test_label_sides_2d_gives_isolated_labels_the_below_band():
    placements = sm._label_sides_2d([(0.0, 0.0, 10.0, "A"), (500.0, 500.0, 10.0, "B")])
    assert [p["direction"] for p in placements] == ["below", "below"]
    assert all(not p["pushed"] for p in placements)


def test_label_sides_2d_cycles_through_all_four_cardinal_directions():
    # Four markers stacked at (nearly) the same point -- claims below,
    # above, right, left in that order before needing to push any of them
    # into a diagonal spot.
    placements = sm._label_sides_2d([
        (100.0, 100.0, 10.0, "AAAAAAAAAA"),
        (100.0, 100.0, 10.0, "BBBBBBBBBB"),
        (100.0, 100.0, 10.0, "CCCCCCCCCC"),
        (100.0, 100.0, 10.0, "DDDDDDDDDD"),
    ])
    assert [p["direction"] for p in placements] == ["below", "above", "right", "left"]
    assert all(not p["pushed"] for p in placements)


def test_label_sides_2d_pushes_a_fifth_label_further_out_instead_of_dropping_it():
    # A fifth marker at the same point has no room left in any of the four
    # cardinal directions at the normal gap -- it should still get a real
    # placement (a further-out "below"/"above" tier) rather than being
    # dropped outright. Uses long (10-char) names specifically because a
    # width-dependent push (e.g. a diagonal offset) would NOT reliably
    # clear a wide label -- the "below"/"above" push here is
    # width-independent by construction (see _LABEL_PUSH_DIRECTIONS'
    # own comment), so this must still succeed regardless of name length.
    placements = sm._label_sides_2d([
        (100.0, 100.0, 10.0, "AAAAAAAAAA"),
        (100.0, 100.0, 10.0, "BBBBBBBBBB"),
        (100.0, 100.0, 10.0, "CCCCCCCCCC"),
        (100.0, 100.0, 10.0, "DDDDDDDDDD"),
        (100.0, 100.0, 10.0, "EEEEEEEEEE"),
    ])
    fifth = placements[4]
    assert fifth is not None
    assert fifth["pushed"] is True
    assert fifth["direction"] in sm._LABEL_PUSH_DIRECTIONS


def test_label_sides_2d_drops_a_label_only_once_every_direction_is_exhausted():
    # Seven markers at the same point exhausts all four cardinal directions
    # at the normal gap plus the "below"/"above" pushed tier (6 real
    # placements) -- the seventh has nowhere left, and only then is it
    # actually dropped.
    entries = [(100.0, 100.0, 10.0, f"NAME{i}") for i in range(7)]
    placements = sm._label_sides_2d(entries)
    assert all(p is not None for p in placements[:6])
    assert placements[6] is None


# --- _body_marker_svg label leader lines ------------------------------------

def test_leader_line_omitted_for_default_below_and_above():
    # "below"/"above" at the normal gap are immediately adjacent to the
    # marker -- the obvious default reading, same as this map's original
    # (pre-collision-avoidance) behavior -- so neither draws a leader line.
    for direction in ("below", "above"):
        svg = sm._body_marker_svg(
            100.0, 100.0, 10.0, "M", "t", "Test Body", "sysmap-planet", {"id": 1},
            label_direction=direction, label_pushed=False,
        )
        assert "sysmap-label-leader" not in svg


def test_leader_line_drawn_for_right_left_and_pushed_placements():
    # "right"/"left" (a less self-evidently-connected direction than
    # below/above) and anything in the further-out pushed tier both get a
    # leader line connecting them back to the marker.
    for direction, pushed in (("right", False), ("left", False), ("below", True), ("above", True)):
        svg = sm._body_marker_svg(
            100.0, 100.0, 10.0, "M", "t", "Test Body", "sysmap-planet", {"id": 1},
            label_direction=direction, label_pushed=pushed,
        )
        assert "sysmap-label-leader" in svg


# --- render_system_map_panel (end-to-end HTML) ------------------------------

def test_single_star_system_renders_star_and_planet_markers():
    system = {"name": "Test System", "binary_configuration": None}
    star = _star(10)
    planet = _planet(1, 10, AU_KM, 0.0)
    html = sm.render_system_map_panel(system, [star], [planet], [])
    assert '<svg class="sysmap-svg" data-scene="system"' in html
    assert 'data-kind="star"' in html
    assert 'data-kind="planet"' in html
    assert "sysmap-hidden" not in html.split("</svg>")[0]  # the system scene itself is never hidden


def test_planet_with_moons_gets_its_own_hidden_drilldown_scene():
    moon = _moon(101, 3.8e5, 0.0)
    planet = _planet(1, 10, AU_KM, 0.0, moons=[moon])
    system = {"name": "Moon Test", "binary_configuration": None}
    html = sm.render_system_map_panel(system, [_star(10)], [planet], [])
    assert 'data-scene="planet-1"' in html
    assert '<svg class="sysmap-svg sysmap-hidden" data-scene="planet-1"' in html
    assert 'data-kind="moon"' in html


def test_belt_renders_as_a_ring_not_a_directional_band():
    belt = _belt(1, 10, AU_KM * 2.7, AU_KM * 2.1, AU_KM * 3.3)
    system = {"name": "Belt Test", "binary_configuration": None}
    html = sm.render_system_map_panel(system, [_star(10)], [], [belt])
    assert '<circle class="sysmap-belt"' in html
    assert 'data-kind="belt"' in html


def test_close_binary_places_both_stars_and_shares_planets_at_the_barycenter():
    primary, secondary = _star(20, mass_kg=2.5e30), _star(21, mass_kg=1.0e30)
    system = {
        "name": "Close Binary", "binary_configuration": "close",
        "binary_mutual_position_x_km": 0.2 * AU_KM, "binary_mutual_position_y_km": 0.0,
        "binary_mutual_position_z_km": 0.0,
    }
    # A 'close' pair's planets orbit the merged proxy -- star_id is NULL.
    planet = _planet(1, None, AU_KM * 3.0, AU_KM * 0.5)
    html = sm.render_system_map_panel(system, [primary, secondary], [planet], [])
    assert html.count('data-kind="star"') == 2
    assert 'data-kind="planet"' in html
    assert "Primary" in html and "Secondary" in html


def test_wide_binary_draws_both_stars_own_independent_planets():
    primary, secondary = _star(30, mass_kg=1.8e30), _star(31, mass_kg=0.6e30)
    system = {
        "name": "Wide Binary", "binary_configuration": "wide",
        "binary_mutual_position_x_km": 50 * AU_KM, "binary_mutual_position_y_km": 20 * AU_KM,
        "binary_mutual_position_z_km": 0.0,
    }
    primary_planet = _planet(1, 30, AU_KM * 1.2, AU_KM * 0.1)
    secondary_planet = _planet(2, 31, AU_KM * 0.4, AU_KM * 0.05)
    html = sm.render_system_map_panel(system, [primary, secondary], [primary_planet, secondary_planet], [])
    # Both planets must actually be drawn -- the old map only ever showed
    # the primary's own planets for a 'wide' pair; this must not regress.
    assert html.count('data-id="1"') >= 1
    assert html.count('data-id="2"') >= 1


def test_empty_system_renders_without_crashing():
    html = sm.render_system_map_panel({"name": "Empty", "binary_configuration": None}, [], [], [])
    assert "<svg" in html
    assert "Nothing to show yet" in html


def test_life_chemical_planet_gets_a_life_badge():
    planet = _planet(1, 10, AU_KM, 0.0, life_chemical="carbon")
    system = {"name": "Life Test", "binary_configuration": None}
    html = sm.render_system_map_panel(system, [_star(10)], [planet], [])
    assert "sysmap-life-badge" in html


def test_many_clustered_planets_do_not_crash_and_still_place_every_body():
    star = _star(50)
    planets = []
    for i in range(15):
        angle = 0.01 * i  # nearly identical angle -> forces real overlap
        distance = AU_KM * (2.0 + 0.001 * i)  # nearly identical distance
        planets.append(_planet(100 + i, 50, distance * math.cos(angle), distance * math.sin(angle)))
    system = {"name": "Cluster Test", "binary_configuration": None}
    html = sm.render_system_map_panel(system, [star], planets, [])
    assert html.count('class="sysmap-body sysmap-planet"') == 15


# --- body preview data (sysmap-preview / static/systemmap.js's own input) --

def test_planet_attrs_carries_a_resolved_class_color():
    planet = _planet(1, 10, AU_KM, 0.0, planet_class="J")
    attrs = sm._planet_attrs(planet)
    assert attrs["color"] == sm._class_color("J")
    assert attrs["color"].startswith("#")


def test_planet_attrs_flags_a_real_atmosphere_but_not_a_missing_one():
    with_atmosphere = sm._planet_attrs(_planet(1, 10, AU_KM, 0.0, atmosphere="Nitrogen-Oxygen"))
    assert with_atmosphere["hasatmosphere"] == "true"
    assert with_atmosphere["atmosphere"] == "Nitrogen-Oxygen"

    literal_none = sm._planet_attrs(_planet(2, 10, AU_KM, 0.0, atmosphere="None"))
    assert literal_none["hasatmosphere"] is None
    assert literal_none["atmosphere"] == "None (airless)"

    unset = sm._planet_attrs(_planet(3, 10, AU_KM, 0.0, atmosphere=None))
    assert unset["hasatmosphere"] is None
    assert unset["atmosphere"] == "None (airless)"


def test_planet_attrs_formats_surface_temperature():
    attrs = sm._planet_attrs(_planet(1, 10, AU_KM, 0.0, surface_temperature_k=287.6))
    assert attrs["surfacetemp"] == "288 K"
    assert sm._planet_attrs(_planet(2, 10, AU_KM, 0.0, surface_temperature_k=None))["surfacetemp"] is None


def test_render_system_map_panel_includes_the_preview_container_only_with_planets():
    system = {"name": "Preview Test", "binary_configuration": None}
    with_planet = sm.render_system_map_panel(system, [_star(10)], [_planet(1, 10, AU_KM, 0.0)], [])
    assert 'id="sysmap-preview"' in with_planet
    assert 'id="sysmap-preview-canvas"' in with_planet

    without_planets = sm.render_system_map_panel(system, [_star(10)], [], [])
    assert 'id="sysmap-preview"' not in without_planets


def test_3d_preview_lives_inside_the_map_viewport_not_a_separate_box():
    # Regression guard: the rotating 3D sphere preview used to sit in the
    # separate `.starmap-side` column next to the map, not in the map
    # itself. It should now be a sibling of the SVG scenes inside
    # `.sysmap-viewport` -- i.e. the viewport's own closing `</div>` comes
    # AFTER `#sysmap-preview`, not before it.
    system = {"name": "Preview Location Test", "binary_configuration": None}
    html = sm.render_system_map_panel(system, [_star(10)], [_planet(1, 10, AU_KM, 0.0)], [])

    viewport_start = html.index('class="starmap-viewport sysmap-viewport"')
    preview_start = html.index('id="sysmap-preview"')
    side_start = html.index('class="starmap-side"')
    assert viewport_start < preview_start < side_start, (
        "#sysmap-preview should be inside .sysmap-viewport, before .starmap-side opens"
    )


def test_render_system_map_panel_planet_marker_carries_preview_data_attrs():
    planet = _planet(1, 10, AU_KM, 0.0, planet_class="J", body_type="g", atmosphere="Hydrogen-Helium",
                      surface_temperature_k=165.0, composition="hydrogen and helium")
    system = {"name": "Preview Data Test", "binary_configuration": None}
    html = sm.render_system_map_panel(system, [_star(10)], [planet], [])
    assert f'data-color="{sm._class_color("J")}"' in html
    assert 'data-hasatmosphere="true"' in html
    assert 'data-atmosphere="Hydrogen-Helium"' in html
    assert 'data-surfacetemp="165 K"' in html
    assert 'data-composition="hydrogen and helium"' in html
