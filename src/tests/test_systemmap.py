"""
planetgen/web/maps/systemmap.py regression tests.

Covers the true-position layout engine that replaced the old schematic
"every planet due east of its star" System Map: the shared log-radial
scale (`_radial_scale_bounds`/`_radial_px`), real-angle placement
(`_polar_to_px`), the mass-weighted binary barycenter split
(`_binary_star_positions_km`), marker de-overlap (`_relax_markers`), 2D
label collision avoidance (`_label_sides_2d`), and the full rendered panel
for single-star/'close'-binary/'wide'-binary/empty systems. no database
needed, since `render_system_map_panel` takes plain dicts, the same shape
`queryDb.system_detail` returns.

Run with: pytest src/tests/test_systemmap.py
"""
import html as html_lib
import math
import random
import re

import pytest  # noqa: E402

from planetgen.web.maps import systemmap as sm  # noqa: E402

AU_KM = 1.496e8


def _star(id_, mass_kg=1.989e30, radius_km=696_000, star_type="G2V", temperature_k=5778.0, luminosity_w=3.828e26):
    return {
        "id": id_, "mass_kg": mass_kg, "radius_km": radius_km, "star_type": star_type,
        "temperature_k": temperature_k, "luminosity_w": luminosity_w,
    }


def _planet(id_, star_id, x_km, y_km, distance_km=None, radius_km=6371.0, planet_class="G", body_type="t",
            moons=None, life_chemical=None, zone="e", period_years=1.0, gravity_g=1.0,
            atmosphere=None, composition=None, surface_temperature_k=None, atmospheric_pressure_pa=None):
    return {
        "id": id_, "star_id": star_id, "name": f"Planet {id_}",
        "position_x_km": x_km, "position_y_km": y_km, "position_z_km": 0.0,
        "distance_km": distance_km if distance_km is not None else math.hypot(x_km, y_km),
        "radius_km": radius_km, "planet_class": planet_class, "body_type": body_type,
        "zone": zone, "period_years": period_years, "gravity_g": gravity_g,
        "life_chemical": life_chemical, "moons": moons or [],
        "atmosphere": atmosphere, "composition": composition, "surface_temperature_k": surface_temperature_k,
        "atmospheric_pressure_pa": atmospheric_pressure_pa,
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
        (350.0, 350.0, 10.0, "AAAAAAAAAA"),
        (350.0, 350.0, 10.0, "BBBBBBBBBB"),
        (350.0, 350.0, 10.0, "CCCCCCCCCC"),
        (350.0, 350.0, 10.0, "DDDDDDDDDD"),
    ])
    assert [p["direction"] for p in placements] == ["below", "above", "right", "left"]
    assert all(not p["pushed"] for p in placements)


@pytest.mark.parametrize("cx, cy", [
    (8.0, 350.0), (692.0, 350.0), (350.0, 8.0), (350.0, 692.0),
    (8.0, 8.0), (692.0, 8.0), (8.0, 692.0), (692.0, 692.0),
])
def test_label_sides_2d_keeps_every_label_inside_the_view(cx, cy):
    # MAP.50: names near any edge or corner of the map stay on screen,
    # even when neighbors take the first choices (a label with no room
    # left is dropped, never drawn off the edge).
    placements = sm._label_sides_2d([(cx, cy, 10.0, "A Very Long Planet Name %d" % n) for n in range(4)])
    assert placements[0] is not None and placements[1] is not None
    for placement in filter(None, placements):
        x1, y1, x2, y2 = placement["rect"]
        assert 0 <= x1 and x2 <= sm._VIEW_SIZE_PX and 0 <= y1 and y2 <= sm._VIEW_SIZE_PX, placement


def test_rendered_label_is_slid_inside_the_view():
    # An outer planet due west with a long name: its label's drawn text
    # (middle-anchored below/above, or start/end-anchored beside) stays
    # within the scene's frame.
    name = "Very Long Outer Planet Name"
    planets = [_planet(1, 10, -30.0 * AU_KM, 0.0), _planet(2, 10, 1.0 * AU_KM, 0.0)]
    planets[0]["name"] = name
    html = sm.render_system_map_panel({"name": "Edge", "binary_configuration": None}, [_star(10)], planets, [])
    match = re.search(r'<text class="sysmap-label" x="([\d.-]+)" y="[\d.-]+" text-anchor="(\w+)">'
                      + re.escape(name), html)
    assert match
    x, anchor = float(match.group(1)), match.group(2)
    width = 2 * sm._label_half_width_px(name)
    left = {"start": x, "middle": x - width / 2, "end": x - width}[anchor]
    lo, _, size, _ = map(float, re.search(r'viewBox="([^"]+)"', html).group(1).split())
    assert left >= lo and left + width <= lo + size


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
        (350.0, 350.0, 10.0, "AAAAAAAAAA"),
        (350.0, 350.0, 10.0, "BBBBBBBBBB"),
        (350.0, 350.0, 10.0, "CCCCCCCCCC"),
        (350.0, 350.0, 10.0, "DDDDDDDDDD"),
        (350.0, 350.0, 10.0, "EEEEEEEEEE"),
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


def _drawn_orbits_and_belts(html):
    orbits = [float(r) for r in re.findall(r'class="sysmap-orbit" cx="[\d.]+" cy="[\d.]+" r="([\d.]+)"', html)]
    belts = [(float(r), float(w)) for r, w in
             re.findall(r'class="sysmap-belt"[^>]*? r="([\d.]+)" stroke-width="([\d.]+)"', html)]
    return orbits, belts


def test_belt_ring_runs_from_its_inner_to_its_outer_edge_on_the_orbit_scale():
    # A belt's distance_km is its inner edge (as generation stores it).
    belt = _belt(1, 10, AU_KM * 2.0, AU_KM * 2.0, AU_KM * 3.6)
    lo, hi = sm._radial_scale_bounds([AU_KM, 2.0 * AU_KM, 3.6 * AU_KM, 6.0 * AU_KM])
    center, width = sm._belt_band(belt, lo, hi)
    assert center - width / 2 == pytest.approx(sm._radial_px(2.0 * AU_KM, lo, hi))
    assert center + width / 2 == pytest.approx(sm._radial_px(3.6 * AU_KM, lo, hi))


def test_thin_belt_is_widened_but_never_over_a_neighboring_orbit():
    belt = _belt(1, 10, AU_KM * 2.0, AU_KM * 2.0, AU_KM * 2.001)
    lo, hi = sm._radial_scale_bounds([AU_KM, 2.0 * AU_KM, 6.0 * AU_KM])
    inner = sm._radial_px(2.0 * AU_KM, lo, hi)
    center, width = sm._belt_band(belt, lo, hi)
    assert width == pytest.approx(sm._BELT_MIN_BAND_PX)
    center, width = sm._belt_band(belt, lo, hi, [inner - 3.0, inner + 5.0])
    assert center - width / 2 >= inner - 3.0 + sm._BELT_ORBIT_GAP_PX - 1e-9
    assert center + width / 2 <= inner + 5.0 - sm._BELT_ORBIT_GAP_PX + 1e-9


@pytest.mark.parametrize("seed", range(40))
def test_no_planet_orbit_is_drawn_inside_a_belt(seed):
    # MAP.49: planets that really orbit clear of a belt (generation keeps
    # them clear) are never drawn inside its ring -- including planets
    # on inclined orbits just past the belt's outer edge, and belts close
    # to their inner neighbor.
    rng = random.Random(seed)
    star_id = 10
    planets, belts, r_au, next_id = [], [], 0.3, 1
    for _ in range(rng.randint(3, 9)):
        r_au *= rng.uniform(1.3, 2.2)
        if rng.random() < 0.35:
            upper = r_au * rng.uniform(1.001, 2.0)
            belts.append(_belt(next_id, star_id, r_au * AU_KM, r_au * AU_KM, upper * AU_KM))
            r_au = upper
        else:
            angle, tilt = rng.uniform(0, 2 * math.pi), math.radians(rng.uniform(0, 8))
            d = r_au * AU_KM
            planet = _planet(next_id, star_id, d * math.cos(tilt) * math.cos(angle),
                             d * math.cos(tilt) * math.sin(angle), distance_km=d)
            planet["position_z_km"] = d * math.sin(tilt)
            planets.append(planet)
        r_au *= 1.02
        next_id += 1
    system = {"name": "Belt Test", "binary_configuration": None}
    html = sm.render_system_map_panel(system, [_star(star_id)], planets, belts)
    orbits, drawn_belts = _drawn_orbits_and_belts(html)
    assert len(drawn_belts) == len(belts)
    for center, width in drawn_belts:
        for orbit in orbits:
            assert not (center - width / 2 < orbit < center + width / 2), (orbit, center, width)


# --- MAP.88: the whole scene fits the frame ---------------------------------

_SCENE_SVG = re.compile(r'<svg class="sysmap-svg[^"]*"[^>]*?viewBox="([^"]+)"[^>]*>(.*?)</svg>', re.S)
_DRAWN = re.compile(r'<(circle|ellipse|polygon|line|text)\b([^>]*)>([^<]*)')
_ATTR = re.compile(r'([\w-]+)="([^"]*)"')


def drawn_outside_view(html, slack=0.15):
    """Every element of every System Map scene that reaches past its own
    scene's viewBox: circles (with a belt's band), ellipses (a gas giant's
    ring, by its long radius), polygons, lines, and names (by the same
    estimated width the layout reserves for them)."""
    outside = []
    for view, inner in _SCENE_SVG.findall(html):
        x0, y0, width, height = map(float, view.split())
        for tag, raw, text in _DRAWN.findall(inner):
            a = dict(_ATTR.findall(raw))
            if tag == "circle":
                r = float(a["r"]) + float(a.get("stroke-width", 0)) / 2
                box = (float(a["cx"]) - r, float(a["cy"]) - r, float(a["cx"]) + r, float(a["cy"]) + r)
            elif tag == "ellipse":
                r = float(a["rx"])
                box = (float(a["cx"]) - r, float(a["cy"]) - r, float(a["cx"]) + r, float(a["cy"]) + r)
            elif tag == "polygon":
                points = [tuple(map(float, p.split(","))) for p in a["points"].split()]
                xs, ys = [p[0] for p in points], [p[1] for p in points]
                box = (min(xs), min(ys), max(xs), max(ys))
            elif tag == "line":
                xs, ys = (float(a["x1"]), float(a["x2"])), (float(a["y1"]), float(a["y2"]))
                box = (min(xs), min(ys), max(xs), max(ys))
            elif "sysmap-label" in a.get("class", ""):
                x, y, w = float(a["x"]), float(a["y"]), 2 * sm._label_half_width_px(html_lib.unescape(text))
                left = {"start": x, "middle": x - w / 2, "end": x - w}[a["text-anchor"]]
                box = (left, y - 12, left + w, y + 4)
            else:
                continue
            if (box[0] < x0 - slack or box[1] < y0 - slack or box[2] > x0 + width + slack
                    or box[3] > y0 + height + slack):
                outside.append((tag, a.get("class"), box, view))
    return outside


def _random_system(rng):
    """A crowded random system: up to 14 planets (gas giants with rings,
    moons) out to hundreds of AU, belts, facilities, and a single star, a
    close pair or a wide pair."""
    kind = rng.choice([None, "close", "wide"])
    stars = [_star(1, radius_km=696_000 * rng.uniform(0.2, 40))]
    system = {"name": "Fit %d" % rng.randint(0, 999), "binary_configuration": kind}
    if kind:
        stars.append(_star(2, radius_km=696_000 * rng.uniform(0.2, 10)))
        sep = (0.1 if kind == "close" else 300) * AU_KM * rng.uniform(0.5, 3)
        angle = rng.uniform(0, 2 * math.pi)
        system.update(binary_mutual_position_x_km=sep * math.cos(angle),
                      binary_mutual_position_y_km=sep * math.sin(angle), binary_mutual_position_z_km=0.0)
    planets, belts, facilities, next_id = [], [], [], 1
    for star_id in ([1, 2] if kind == "wide" else [None if kind else 1]):
        r_au = rng.uniform(0.05, 1.0)
        for _ in range(rng.randint(1, 14)):
            r_au *= rng.uniform(1.2, 2.5)
            if rng.random() < 0.2:
                upper = r_au * rng.uniform(1.05, 2.0)
                belts.append(_belt(next_id, star_id, r_au * AU_KM, r_au * AU_KM, upper * AU_KM))
                r_au = upper
            else:
                angle = rng.uniform(0, 2 * math.pi)
                giant = rng.random() < 0.4
                moons = [_moon(100 * next_id + m, *(lambda d, t: (d * math.cos(t), d * math.sin(t)))(
                    rng.uniform(2e5, 5e7), rng.uniform(0, 2 * math.pi)), radius_km=rng.uniform(5, 5000))
                    for m in range(rng.randint(0, 6))]
                planets.append(_planet(next_id, star_id, r_au * AU_KM * math.cos(angle),
                                       r_au * AU_KM * math.sin(angle),
                                       radius_km=rng.uniform(70_000, 250_000) if giant else rng.uniform(500, 20_000),
                                       planet_class="J" if giant else "G", body_type="g" if giant else "t",
                                       moons=moons))
                planets[-1]["name"] = rng.choice(["Ib", "Very Long Outermost Planet Name", "Kepler Prime"])
                if rng.random() < 0.3:
                    facilities.append({"id": next_id, "name": "Dock", "kind": "station", "host_type": "planet",
                                       "host_id": next_id, "host_name": None, "placement": "orbital",
                                       "orbit_distance_km": 4e4, "orbit_period_years": 0.01,
                                       "orbital_speed_kms": 3.0, "orbit_phase_deg": rng.uniform(0, 360)})
            next_id += 1
    facilities.append({"id": 999, "name": "Far Gate", "kind": "outpost", "host_type": "star", "host_id": 1,
                       "host_name": None, "placement": "orbital", "orbit_distance_km": 2000 * AU_KM,
                       "orbit_period_years": 9e4, "orbital_speed_kms": 0.5, "orbit_phase_deg": 45.0})
    return system, stars, planets, belts, facilities


@pytest.mark.parametrize("seed", range(60))
def test_every_drawn_element_sits_inside_the_view(seed):
    # MAP.88: outer planets, their rings and moons, belts, facilities and
    # names never run past the edge of any scene.
    html = sm.render_system_map_panel(*_random_system(random.Random(seed)))
    assert drawn_outside_view(html) == []


def test_a_scene_that_fits_keeps_the_fixed_frame():
    # The outermost orbit is 335 px out; a small world there fits.
    planet = _planet(1, 10, AU_KM, 0.0, radius_km=100.0)
    planet["name"] = "B"
    html = sm.render_system_map_panel({"name": "Small", "binary_configuration": None}, [_star(10)], [planet], [])
    assert set(re.findall(r'viewBox="([^"]+)"', html)) == {"0 0 700 700"}


def test_a_scene_that_overflows_widens_its_view_evenly_around_the_center():
    # An outermost gas giant's ring reaches past the 700 px frame, so the
    # view grows (the scene shrinks on screen) and stays centered on the star.
    giant = _planet(2, 10, 30 * AU_KM, 0.0, radius_km=200_000, planet_class="J", body_type="g")
    planets = [_planet(1, 10, AU_KM, 0.0), giant]
    html = sm.render_system_map_panel({"name": "Wide", "binary_configuration": None}, [_star(10)], planets, [])
    x0, y0, width, height = map(float, re.search(r'viewBox="([^"]+)"', html).group(1).split())
    assert x0 < 0 and x0 == y0 and width == height
    assert x0 + width / 2 == pytest.approx(sm._CENTER_PX, abs=0.1)
    assert drawn_outside_view(html) == []


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


# --- body preview data (sysmap-spheres-canvas / static/systemmap.js's own input) --

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
    assert attrs["surfacetemp"] == "288 K (14 °C, 58 °F)"
    assert sm._planet_attrs(_planet(2, 10, AU_KM, 0.0, surface_temperature_k=None))["surfacetemp"] is None


def test_planet_attrs_formats_surface_pressure():
    # MAP.117: under the temperature in the side panel.
    attrs = sm._planet_attrs(_planet(1, 10, AU_KM, 0.0, atmosphere="Nitrogen-Oxygen",
                                     atmospheric_pressure_pa=101_325.0))
    assert attrs["surfacepressure"] == "101 kPa (1 atm, 14.7 psi)"
    assert sm._planet_attrs(_planet(2, 10, AU_KM, 0.0, atmosphere="None"))["surfacepressure"] == "None"
    no_value = sm._planet_attrs(_planet(3, 10, AU_KM, 0.0, atmosphere="Nitrogen", atmospheric_pressure_pa=None))
    assert no_value["surfacepressure"] is None


def test_render_system_map_panel_includes_the_spheres_canvas_whenever_there_are_stars():
    system = {"name": "Preview Test", "binary_configuration": None}
    with_planet = sm.render_system_map_panel(system, [_star(10)], [_planet(1, 10, AU_KM, 0.0)], [])
    assert 'id="sysmap-spheres-canvas"' in with_planet

    # A star alone (no planets/belts) still gets a marker of its own in the
    # "system" scene -- the shared canvas is guarded on `stars`, not
    # `planets`, since every star's marker also gets a live sphere now, not
    # just planets/moons.
    star_only = sm.render_system_map_panel(system, [_star(10)], [], [])
    assert 'id="sysmap-spheres-canvas"' in star_only

    no_stars = sm.render_system_map_panel(system, [], [], [])
    assert 'id="sysmap-spheres-canvas"' not in no_stars


def test_spheres_canvas_lives_inside_the_map_viewport_not_a_separate_box():
    # Regression guard: the rotating 3D sphere layer must be a sibling of
    # the SVG scenes inside `.sysmap-viewport` (so `style.css` can stack it
    # directly behind them by z-index), not off in the separate
    # `.starmap-side` info column -- i.e. the viewport's own closing
    # `</div>` comes AFTER `#sysmap-spheres-canvas`, not before it.
    system = {"name": "Preview Location Test", "binary_configuration": None}
    html = sm.render_system_map_panel(system, [_star(10)], [_planet(1, 10, AU_KM, 0.0)], [])

    viewport_start = html.index('class="starmap-viewport sysmap-viewport"')
    canvas_start = html.index('id="sysmap-spheres-canvas"')
    side_start = html.index('class="starmap-side"')
    assert viewport_start < canvas_start < side_start, (
        "#sysmap-spheres-canvas should be inside .sysmap-viewport, before .starmap-side opens"
    )


def test_render_system_map_panel_planet_marker_carries_preview_data_attrs():
    planet = _planet(1, 10, AU_KM, 0.0, planet_class="J", body_type="g", atmosphere="Hydrogen-Helium",
                      surface_temperature_k=165.0, composition="hydrogen and helium")
    system = {"name": "Preview Data Test", "binary_configuration": None}
    html = sm.render_system_map_panel(system, [_star(10)], [planet], [])
    assert f'data-color="{sm._class_color("J")}"' in html
    assert 'data-hasatmosphere="true"' in html
    assert 'data-atmosphere="Hydrogen-Helium"' in html
    assert 'data-surfacetemp="165 K (-108 °C, -163 °F)"' in html
    assert 'data-composition="hydrogen and helium"' in html


# --- Measure distance (real km data attrs + toggle button) -----------------

def test_planet_marker_carries_real_km_position_not_the_log_scaled_pixel_one():
    # static/systemmap.js's own "measure distance" feature needs the real,
    # un-log-scaled position -- confirms _planet_attrs hands over the exact
    # stored position_x/y_km, not something derived from the drawn marker's
    # own (log-scaled) pixel radius.
    planet = _planet(1, 10, 3.0 * AU_KM, 4.0 * AU_KM)
    attrs = sm._planet_attrs(planet)
    assert attrs["xkm"] == pytest.approx(3.0 * AU_KM)
    assert attrs["ykm"] == pytest.approx(4.0 * AU_KM)
    assert attrs["radiuskm"] == planet["radius_km"]


def test_moon_marker_also_carries_real_km_position():
    moon = _moon(1, 1.0e5, -2.0e5)
    attrs = sm._planet_attrs(moon, kind="moon", parent_name="Planet 1")
    assert attrs["xkm"] == pytest.approx(1.0e5)
    assert attrs["ykm"] == pytest.approx(-2.0e5)
    assert attrs["radiuskm"] == moon["radius_km"]


def test_render_system_map_panel_star_marker_carries_radiuskm():
    system = {"name": "Radius Test", "binary_configuration": None}
    html = sm.render_system_map_panel(system, [_star(10, radius_km=696_000)], [_planet(1, 10, AU_KM, 0.0)], [])
    assert 'data-radiuskm="696000"' in html


def test_render_system_map_panel_wide_binary_star_marker_carries_radiuskm():
    system = {
        "name": "Wide Binary Radius Test", "binary_configuration": "wide",
        "binary_mutual_position_x_km": 5000.0 * AU_KM, "binary_mutual_position_y_km": 0.0,
        "binary_mutual_position_z_km": 0.0,
    }
    stars = [_star(10, radius_km=696_000), _star(11, radius_km=350_000)]
    html = sm.render_system_map_panel(system, stars, [_planet(1, 10, AU_KM, 0.0), _planet(2, 11, AU_KM, 0.0)], [])
    assert 'data-radiuskm="696000"' in html
    assert 'data-radiuskm="350000"' in html


def test_render_system_map_panel_includes_measure_distance_button_whenever_there_are_stars():
    system = {"name": "Measure Button Test", "binary_configuration": None}
    with_star = sm.render_system_map_panel(system, [_star(10)], [_planet(1, 10, AU_KM, 0.0)], [])
    assert 'id="sysmap-measure-btn"' in with_star

    no_stars = sm.render_system_map_panel(system, [], [], [])
    assert 'id="sysmap-measure-btn"' not in no_stars


def test_close_binary_star_markers_carry_their_own_real_km_offset_not_the_barycenter_origin():
    # Regression guard: static/systemmap.js's own obstacle-avoidance math
    # assumes the scene's "self" .sysmap-star marker's own data-xkm/
    # data-ykm IS the obstacle circle's real center -- for a close binary
    # that's each star's own small real offset from the shared
    # barycenter (see _binary_star_positions_km), never (0, 0) for both
    # (only a single, non-binary star sits exactly at the scene's own
    # local origin).
    primary, secondary = _star(20, mass_kg=2.5e30), _star(21, mass_kg=1.0e30)
    system = {
        "name": "Close Binary", "binary_configuration": "close",
        "binary_mutual_position_x_km": 0.2 * AU_KM, "binary_mutual_position_y_km": 0.0,
        "binary_mutual_position_z_km": 0.0,
    }
    expected = sm._binary_star_positions_km(system, [primary, secondary])
    html = sm.render_system_map_panel(system, [primary, secondary], [], [])

    for star_id in (20, 21):
        ex, ey = expected[star_id]
        assert f'data-xkm="{ex}"' in html
        assert f'data-ykm="{ey}"' in html


def test_wide_binary_companion_marker_carries_the_real_separation_not_its_drawn_position():
    # The companion star's marker is drawn at a fixed, merely
    # representative pixel position (see _render_wide_binary_scenes's own
    # docstring -- the real separation is routinely 10-1000x a planet's
    # own orbit and would blow out the frame) -- but a "measure distance"
    # click needs the REAL separation, not that drawn position, so
    # data-xkm/data-ykm must carry binary_mutual_position_x/y_km exactly,
    # decoupled from wherever the marker itself was actually drawn.
    primary, secondary = _star(30, mass_kg=1.8e30), _star(31, mass_kg=0.6e30)
    system = {
        "name": "Wide Binary", "binary_configuration": "wide",
        "binary_mutual_position_x_km": 50 * AU_KM, "binary_mutual_position_y_km": 20 * AU_KM,
        "binary_mutual_position_z_km": 0.0,
    }
    html = sm.render_system_map_panel(system, [primary, secondary], [], [])
    assert f'data-xkm="{50 * AU_KM}"' in html
    assert f'data-ykm="{20 * AU_KM}"' in html
    # The primary, in its own "system" scene, is always the local origin.
    assert 'data-xkm="0.0"' in html
    assert 'data-ykm="0.0"' in html


# --- MAP.92: radius and mass in the side panel --------------------------------

def _marker_attrs(html, kind, body_id):
    tag = re.search(r'<g class="sysmap-body[^"]*"[^>]*data-kind="%s" data-id="%d"[^>]*>' % (kind, body_id), html)
    assert tag, (kind, body_id)
    return dict(_ATTR.findall(tag.group(0)))


def test_planet_and_moon_markers_carry_radius_and_mass():
    moon = dict(_moon(7, 4e5, 0.0, radius_km=1737.0), mass_kg=7.35e22)
    earth = dict(_planet(1, 10, AU_KM, 0.0, moons=[moon]), mass_kg=5.972e24)
    giant = dict(_planet(2, 10, 5 * AU_KM, 0.0, radius_km=69_911.0, planet_class="J", body_type="g"),
                 mass_kg=1.898e27)
    html = sm.render_system_map_panel({"name": "Sol", "binary_configuration": None}, [_star(10)], [earth, giant], [])
    planet = _marker_attrs(html, "planet", 1)
    assert planet["data-radius"] == "6.37 × 10³ km (1.00 Earth radii)"
    assert planet["data-mass"] == "5.97 × 10²⁴ kg (1.00 Earth masses)"
    assert _marker_attrs(html, "planet", 2)["data-mass"] == "1.90 × 10²⁷ kg (1.00 Jupiter masses)"
    moon_attrs = _marker_attrs(html, "moon", 7)
    assert moon_attrs["data-radius"] == "1.74 × 10³ km (0.273 Earth radii)"
    assert moon_attrs["data-mass"] == "7.35 × 10²² kg (0.0123 Earth masses)"
    # The planet drilled into, at the center of its moon scene, too.
    center = re.search(r'<g class="sysmap-body sysmap-planet" [^>]*data-self="true"[^>]*>', html).group(0)
    assert 'data-mass="5.97 × 10²⁴ kg (1.00 Earth masses)"' in center


def test_unknown_radius_and_mass_show_a_dash():
    planet = dict(_planet(1, 10, AU_KM, 0.0, radius_km=None), mass_kg=None)
    html = sm.render_system_map_panel({"name": "Sol", "binary_configuration": None}, [_star(10)], [planet], [])
    attrs = _marker_attrs(html, "planet", 1)
    assert attrs["data-radius"] == attrs["data-mass"] == "–"
