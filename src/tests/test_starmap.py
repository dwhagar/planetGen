"""
planetgen/web/maps/starmap.py regression tests.

Covers `_rotate_to_galaxy_frame` -- the fix for the Sector Map's star dots
being plotted as if their own sector-local (x, y, z) axes already ran
parallel to the galaxy frame's, when the wedge outline/"Galactic Center"
compass arrow were always correctly expressed in the galaxy frame
directly -- via `planetgen.galaxy.geometry.sector_orientation`'s
convention (local +X radially outward from the galactic axis, +Y along the
ring, +Z galactic north), applied at render time from the sector's own
stored `center_x/y/z_pc` rather than needing any new stored orientation.

Also covers `render_map_panel`'s scene-data contract -- since the Sector
Map moved from server-rendered CSS `<div>`s to a `<script
type="application/json">` block a WebGL client reads (see starmap.py's own
module docstring), these tests parse that JSON instead of regex-matching
HTML, but check the same underlying behavior the old assertions did
(rotation agreement, phenomenon labels/kinds, radius scaling, click hint
text).

no database
needed, since `render_map_panel`/`_rotate_to_galaxy_frame` take plain
dicts/tuples.

Run with: pytest src/tests/test_starmap.py
"""
import json
import math
import os
import re

import pytest  # noqa: E402

from planetgen.web.maps.starmap import _rotate_to_galaxy_frame, render_map_panel  # noqa: E402


def _link(name, **params):
    """Stands in for `web.helpers.page_url`: `/<name>?k=v&...`."""
    return f"/{name}?" + "&".join(f"{key}={value}" for key, value in sorted(params.items()))
from planetgen.web.maps.starmap import _phenomenon_cloud_radius_px  # noqa: E402


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
    # center_pc exactly on the galactic axis is the one place "radially
    # outward" is undefined (see sector_orientation's own fallback) -- must still
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

    scene_unplaced = _scene_data(render_map_panel(_link, 1000.0, None, None, [system]))
    unplaced_pos = _star_position(scene_unplaced)

    scene_placed = _scene_data(render_map_panel(_link, 1000.0, (3, 0, 7), (500.0, 200.0, -100.0), [system]))
    placed_pos = _star_position(scene_placed)

    assert placed_pos != unplaced_pos

    scene_unplaced_again = _scene_data(render_map_panel(_link, 1000.0, None, None, [system]))
    assert _star_position(scene_unplaced_again) == unplaced_pos


def test_scene_stars_carry_a_point_of_light_sized_by_luminosity():
    """MAP.15: each star in the scene data carries the point of light
    `sectormap.js` draws (`light`, in screen pixels): a supergiant's halo
    is far wider and stronger than a red dwarf's, its core bigger, and
    its color follows its temperature."""
    dwarf = _make_system()
    giant = _make_system(x=-100.0)
    giant["id"] = 2
    giant["stars"] = [{
        "star_type": "M2IA", "temperature_k": 3600, "radius_km": 696000 * 800,
        "luminosity_w": 3.828e26 * 2e5, "temp_display": "3600 K",
    }]
    scene = _scene_data(render_map_panel(_link, 1000.0, None, None, [dwarf, giant]))
    small, big = scene["stars"]
    for star in (small, big):
        assert set(star["light"]) == {"color", "corePx", "sizePx", "glow", "bright"}
        assert star["r"] <= 8
        # A point a few pixels across, never a ball.
        assert star["light"]["corePx"] <= 5
        assert star["light"]["sizePx"] <= 40
        assert star["light"]["sizePx"] >= 2 * star["light"]["corePx"] + 2
    # (A Sun's halo is a little wider since MAP.87 drew the faint end
    # brighter.)
    assert big["light"]["sizePx"] > 1.6 * small["light"]["sizePx"]
    assert big["light"]["corePx"] > small["light"]["corePx"]
    # The halo's light (strength over its area).
    assert big["light"]["glow"] * big["light"]["sizePx"] ** 2 > 2 * small["light"]["glow"] * small["light"]["sizePx"] ** 2
    assert big["light"]["bright"] > small["light"]["bright"]
    # A 3600 K supergiant is orange-red: more red than blue.
    color = big["light"]["color"]
    assert int(color[1:3], 16) > int(color[5:7], 16)


@pytest.mark.parametrize("type_, descriptor, lit", [
    ("quasar", "radio-loud", True),
    ("neutron_star", "pulsar", True),
    ("black_hole", "accreting", True),
    ("black_hole", "quiescent", False),
    ("rogue_planet", "terrestrial", True),
    ("interstellar_comet", "icy", False),
    ("nebula", "emission", False),
    ("supernova_remnant", "shell", False),
    ("asteroid_field", "dense", False),
])
def test_only_light_giving_phenomena_are_points_of_light(type_, descriptor, lit):
    """MAP.15: quasars, neutron stars and accreting black holes are drawn
    as points of light, and so is a rogue planet, faintly (MAP.82); a
    quiescent black hole and an interstellar comet keep their spheres,
    and clouds stay clouds."""
    phenomenon = _phenomenon(type_=type_, descriptor=descriptor)
    scene = _scene_data(render_map_panel(_link, 1000.0, None, None, [_make_system()], phenomena=[phenomenon]))
    cloud = scene["clouds"][0]
    assert ("light" in cloud) is lit
    if lit:
        assert set(cloud["light"]) >= {"color", "corePx", "sizePx", "glow", "bright"}


def test_render_map_panel_placed_on_axis_matches_unplaced():
    """
    A galaxy placement exactly on the galactic axis is where
    `sector_orientation` falls back to the galaxy's own axes, so rotating
    into the galaxy frame there is the identity transform -- confirms the
    rotation is doing real work in the general (off-axis) case above, not
    coincidentally matching for an unrelated reason.
    """
    system = _make_system()
    scene_unplaced = _scene_data(render_map_panel(_link, 1000.0, None, None, [system]))
    scene_placed_on_axis = _scene_data(render_map_panel(_link, 1000.0, (3, 0, 7), (0.0, 0.0, 999.0), [system]))
    # Only the star's own position should match; the placed scene also
    # carries a cell outline/compass arrow the unplaced one doesn't.
    assert _star_position(scene_placed_on_axis) == _star_position(scene_unplaced)


def test_render_map_panel_outline_is_a_cell_when_placed_and_a_cube_otherwise():
    system = _make_system()
    scene_unplaced = _scene_data(render_map_panel(_link, 1000.0, None, None, [system]))
    scene_placed = _scene_data(render_map_panel(_link, 1000.0, (3, 0, 7), (500.0, 200.0, -100.0), [system]))
    assert scene_unplaced["outline"]["kind"] == "cube"
    assert scene_placed["outline"]["kind"] == "cell"
    # 12 edges, each a line of 3D points, for either shape: the cube's are
    # straight, and the cell's 4 slot-angle edges are sampled arcs.
    for scene in (scene_unplaced, scene_placed):
        assert len(scene["outline"]["edges"]) == 12
        for edge in scene["outline"]["edges"]:
            assert len(edge) >= 2
            assert all(len(point) == 3 for point in edge)
    assert all(len(edge) == 2 for edge in scene_unplaced["outline"]["edges"])
    lengths = sorted(len(edge) for edge in scene_placed["outline"]["edges"])
    assert lengths[:8] == [2] * 8
    assert all(n > 2 for n in lengths[8:])


def test_cell_outline_arcs_follow_the_ring_radius():
    from planetgen.web.maps import starmap
    from planetgen.galaxy.geometry import ring_bounds_pc, sector_position_pc
    from stellarObjects.utils import mpc_to_pc

    address, edge_mpc, half_edge = (3, 0, 7), 4000.0, 2000.0
    edges = starmap._cell_arc_edges_px(address, edge_mpc, half_edge)
    edge_pc = mpc_to_pc(edge_mpc)
    center = sector_position_pc(*address, edge_pc)
    radii = ring_bounds_pc(3, edge_pc)
    scale = starmap._SCENE_HALF_PX / half_edge
    for i, (a, b) in enumerate(starmap._EDGE_PAIRS):
        if a ^ b != 1:
            continue
        # Every sample of an arc edge lies on the inner or outer ring radius.
        for x, y, _z in edges[i]:
            gx = center[0] + mpc_to_pc(x / scale)
            gy = center[1] + mpc_to_pc(-y / scale)
            assert min(abs(math.hypot(gx, gy) - r) for r in radii) < 1e-6


def test_render_map_panel_compass_present_only_when_placed():
    system = _make_system()
    scene_unplaced = _scene_data(render_map_panel(_link, 1000.0, None, None, [system]))
    scene_placed = _scene_data(render_map_panel(_link, 1000.0, (3, 0, 7), (500.0, 200.0, -100.0), [system]))
    assert scene_unplaced["compass"] is None
    assert scene_placed["compass"] is not None
    assert scene_placed["compass"]["label"] == "N"


def test_render_map_panel_binary_system_gets_two_star_entries():
    system = _make_system()
    system["stars"][0]["name"] = "Test System Kelmoor"
    system["stars"].append({
        "name": "Test System Ostra", "star_type": "M4V", "temperature_k": 3200, "radius_km": 200000,
        "luminosity_w": 1.0e24, "temp_display": "3200 K",
    })
    scene = _scene_data(render_map_panel(_link, 1000.0, None, None, [system]))
    assert len(scene["stars"]) == 2
    # Each star shows its own name, with no A/B letters (bodyNames.py).
    assert scene["stars"][0]["name"] == "Test System Kelmoor"
    assert scene["stars"][1]["name"] == "Test System Ostra"
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
    from planetgen.web.maps.starmap import _MAX_CLOUD_RADIUS_PX
    huge = _phenomenon_cloud_radius_px(200.0, half_edge=100.0)
    assert huge <= _MAX_CLOUD_RADIUS_PX


def test_render_map_panel_draws_a_nebula_cloud_with_its_own_kind_and_colors():
    system = _make_system()
    phenomenon = _phenomenon(type_="nebula", descriptor="reflection")
    html = render_map_panel(_link, 1000.0, None, None, [system], phenomena=[phenomenon])
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
    scene = _scene_data(render_map_panel(_link, 1000.0, None, None, [system], phenomena=[phenomenon]))
    cloud = scene["clouds"][0]
    assert cloud["kind"] == "asteroidField"
    assert cloud["typeLabel"] == "Asteroid Field (Dense)"
    assert "coreColor" not in cloud


def test_render_map_panel_without_phenomena_matches_omitting_the_argument():
    system = _make_system()
    html_default = render_map_panel(_link, 1000.0, None, None, [system])
    html_explicit_empty = render_map_panel(_link, 1000.0, None, None, [system], phenomena=[])
    html_none = render_map_panel(_link, 1000.0, None, None, [system], phenomena=None)
    assert _scene_data(html_default)["clouds"] == []
    assert html_default == html_explicit_empty == html_none


def test_render_map_panel_places_a_phenomenon_directly_in_the_galaxy_frame_unrotated():
    # A phenomenon's offset_*_ly is already galaxy-frame (see
    # queryDb.phenomena_near_sector) -- unlike a star system's sector-local
    # x/y/z, it must NOT be re-rotated by `_rotate_to_galaxy_frame` even
    # when the sector itself has a galaxy placement.
    phenomenon = _phenomenon(offset=(3.0, 0.0, 0.0))
    scene_unplaced = _scene_data(render_map_panel(_link, 1000.0, None, None, [], phenomena=[phenomenon]))
    scene_placed = _scene_data(render_map_panel(_link, 1000.0, (3, 0, 7), (500.0, 200.0, -100.0), [], phenomena=[phenomenon]))

    def _cloud_position(scene):
        cloud = scene["clouds"][0]
        return (cloud["x"], cloud["y"], cloud["z"])

    assert _cloud_position(scene_unplaced) == _cloud_position(scene_placed)


def test_phenomenon_cloud_radius_floors_at_minimum_for_a_point_like_object():
    # black_hole/neutron_star always query radius_ly as a literal 0 (see
    # queryDb._PHENOMENON_TABLES) -- must still floor at a visible minimum,
    # not collapse to an invisible 0px marker.
    from planetgen.web.maps.starmap import _MAX_CLOUD_RADIUS_PX
    radius = _phenomenon_cloud_radius_px(0.0, half_edge=5000.0)
    assert 0 < radius <= _MAX_CLOUD_RADIUS_PX


def test_render_map_panel_draws_an_accreting_black_hole_point():
    system = _make_system()
    phenomenon = _phenomenon(type_="black_hole", descriptor="accreting", radius_ly=0)
    scene = _scene_data(render_map_panel(_link, 1000.0, None, None, [system], phenomena=[phenomenon]))
    cloud = scene["clouds"][0]
    assert cloud["kind"] == "blackHoleAccreting"
    assert cloud["typeLabel"] == "Black Hole (Accreting)"
    assert cloud["radiusText"] is None  # a point has no radius row


def test_render_map_panel_draws_a_quasar_point():
    system = _make_system()
    phenomenon = _phenomenon(type_="quasar", descriptor="radio-loud", radius_ly=0)
    scene = _scene_data(render_map_panel(_link, 1000.0, None, None, [system], phenomena=[phenomenon]))
    cloud = scene["clouds"][0]
    assert cloud["kind"] == "quasar"
    assert cloud["typeLabel"] == "Quasar (Radio-loud)"


def test_render_map_panel_draws_a_quiescent_black_hole_point():
    system = _make_system()
    phenomenon = _phenomenon(type_="black_hole", descriptor="quiescent", radius_ly=0)
    scene = _scene_data(render_map_panel(_link, 1000.0, None, None, [system], phenomena=[phenomenon]))
    cloud = scene["clouds"][0]
    assert cloud["kind"] == "blackHoleQuiescent"
    assert cloud["typeLabel"] == "Black Hole (Quiescent)"


def test_render_map_panel_draws_a_neutron_star_point():
    system = _make_system()
    phenomenon = _phenomenon(type_="neutron_star", descriptor="millisecond", radius_ly=0)
    scene = _scene_data(render_map_panel(_link, 1000.0, None, None, [system], phenomena=[phenomenon]))
    cloud = scene["clouds"][0]
    assert cloud["kind"] == "neutronStar"
    assert cloud["typeLabel"] == "Neutron Star (Millisecond)"


def test_render_map_panel_neutron_star_with_no_descriptor_falls_back_to_plain_label():
    system = _make_system()
    phenomenon = _phenomenon(type_="neutron_star", descriptor=None, radius_ly=0)
    scene = _scene_data(render_map_panel(_link, 1000.0, None, None, [system], phenomena=[phenomenon]))
    assert scene["clouds"][0]["typeLabel"] == "Neutron Star"


def _neighbor(address=(6, -2, 9), direction_pc=(1.0, 0.0, 0.0), exists=False, sector_id=None,
              sector_name=None, designation="600000063"):
    return {
        "ring_index": address[0], "layer_index": address[1], "ring_slot_index": address[2],
        "direction_pc": list(direction_pc),
        "exists": exists, "sector_id": sector_id, "sector_name": sector_name, "designation": designation,
    }


def test_render_map_panel_neighbors_empty_by_default():
    system = _make_system()
    scene = _scene_data(render_map_panel(_link, 1000.0, (5, 0, 17), (500.0, 200.0, -100.0), [system]))
    assert scene["neighbors"] == []


def test_render_map_panel_existing_neighbor_carries_its_href():
    system = _make_system()
    neighbor = _neighbor(exists=True, sector_id=42, sector_name="Neighboring Sector")
    scene = _scene_data(render_map_panel(
        _link, 1000.0, (5, 0, 17), (500.0, 200.0, -100.0), [system], neighbors=[neighbor],
    ))
    entry = scene["neighbors"][0]
    assert entry["exists"] is True
    assert entry["name"] == "Neighboring Sector"
    assert entry["href"] == "/sector?sector_id=42"
    assert "navTarget" not in entry and "navParams" not in entry
    assert (entry["ringIndex"], entry["layerIndex"], entry["ringSlotIndex"]) == (6, -2, 9)


def test_render_map_panel_missing_neighbor_carries_no_href():
    system = _make_system()
    neighbor = _neighbor(exists=False, designation="ABCDEF")
    scene = _scene_data(render_map_panel(
        _link, 1000.0, (5, 0, 17), (500.0, 200.0, -100.0), [system], neighbors=[neighbor],
    ))
    entry = scene["neighbors"][0]
    assert entry["exists"] is False
    assert "href" not in entry
    assert entry["designation"] == "ABCDEF"


def test_render_map_panel_neighbor_indicator_sits_along_its_own_direction():
    # Placed at a fixed reach past the scene's own edge, in the given
    # direction -- along +x here, so y/z should stay at (approximately)
    # zero and x should be positive and clearly past the scene's own
    # half-edge (see planetgen/web/maps/starmap.py's own `_NEIGHBOR_INDICATOR_REACH`).
    system = _make_system()
    neighbor = _neighbor(direction_pc=(1.0, 0.0, 0.0))
    scene = _scene_data(render_map_panel(
        _link, 1000.0, (5, 0, 17), (500.0, 200.0, -100.0), [system], neighbors=[neighbor],
    ))
    entry = scene["neighbors"][0]
    assert entry["x"] > scene["sceneHalfPx"]
    assert entry["y"] == pytest.approx(0.0, abs=1e-9)
    assert entry["z"] == pytest.approx(0.0, abs=1e-9)


def test_json_script_escapes_script_close_tag_in_a_name():
    # A system name is arbitrary user-supplied text (see `--name`) -- one
    # containing "</script>" must not be able to break out of the
    # embedded JSON block (see starmap.py's own `_json_script`).
    system = _make_system()
    system["name"] = "Evil</script><script>alert(1)</script>"
    html = render_map_panel(_link, 1000.0, None, None, [system])
    assert "</script><script>alert" not in html
    scene = _scene_data(html)
    assert scene["stars"][0]["name"].startswith("Evil</script>")


def test_render_map_panel_links_are_plain_hrefs():
    system = _make_system()
    phenomenon = {"id": 7, "type": "nebula", "name": "Veil & <Co>", "descriptor": "emission",
                  "radius_ly": 2.0, "distance_ly": 1.0,
                  "offset_x_ly": 0.0, "offset_y_ly": 0.0, "offset_z_ly": 0.0}
    html = render_map_panel(_link, 1000.0, None, None, [system], phenomena=[phenomenon])
    scene = _scene_data(html)
    assert scene["stars"][0]["href"] == "/system?system_id=1"
    assert scene["clouds"][0]["href"] == "/phenomenon?phenomenon_id=7&phenomenon_type=nebula"
    # The no-JavaScript fallback list: escaped <a href> links, no forms.
    assert '<li><a href="/system?system_id=1">Test System</a></li>' in html
    assert '<a href="/phenomenon?phenomenon_id=7&amp;phenomenon_type=nebula">Veil &amp; &lt;Co&gt;</a>' in html
    assert "<form" not in html and "data-nav" not in html


# --- MAP.87: the faint end drawn brighter ---------------------------------------------

def _light_for(luminosity_solar, radius_solar=1.0, temperature_k=5772.0):
    from planetgen.web.maps import starmap

    return starmap._star_light(luminosity_solar * starmap.SOLAR_LUMINOSITY, radius_solar * 696000.0, temperature_k)


def _unboosted(luminosity_solar, radius_solar=1.0):
    """`_star_light`'s sizes before MAP.87 (no boost)."""
    from planetgen.web.maps import starmap

    share = starmap._log_share(luminosity_solar, starmap._LIGHT_LOG_LUMINOSITY)
    core = starmap._lerp(starmap._LIGHT_CORE_PX, starmap._log_share(radius_solar, starmap._LIGHT_LOG_RADIUS))
    return {
        "sizePx": round(max(starmap._lerp(starmap._LIGHT_SIZE_PX, share * share), 2 * core + 2), 2),
        "glow": round(starmap._lerp(starmap._LIGHT_GLOW, share), 3),
        "bright": round(starmap._lerp(starmap._LIGHT_BRIGHT, share), 3),
    }


def test_light_boost_is_four_at_the_dim_end_and_none_from_1000_suns():
    from planetgen.web.maps.starmap import star_light_boost

    assert star_light_boost(1e-4) == pytest.approx(4.0)
    assert star_light_boost(1e-6) == pytest.approx(4.0)
    assert star_light_boost(1000.0) == pytest.approx(1.0)
    assert star_light_boost(1e6) == pytest.approx(1.0)
    assert star_light_boost(None) == 1.0 and star_light_boost(0.0) == 1.0


def test_light_boost_tapers_without_a_jump():
    from planetgen.web.maps.starmap import star_light_boost

    luminosities = [10 ** (e / 20) for e in range(-100, 81)]
    boosts = [star_light_boost(lum) for lum in luminosities]
    assert all(a >= b for a, b in zip(boosts, boosts[1:])), "brighter stars never get more boost"
    assert max(a - b for a, b in zip(boosts, boosts[1:])) < 0.1, "no jump anywhere"
    assert star_light_boost(1.0) == pytest.approx(2.2, abs=0.05)


def test_a_brighter_star_is_never_drawn_fainter():
    """The halo's light (strength x area) never falls as luminosity
    rises, on the Sector Map and with the Galaxy Map's own ranges
    (galaxymap3d.js's STAR_MIN_PX..STAR_MAX_PX and STAR_GLOW)."""
    from planetgen.web.maps import starmap

    def galaxy_light(lum):
        share = starmap._log_share(lum, (-4.0, 6.0))
        size, glow = starmap._boost_light(6 + 24 * share * share, 0.2 + 0.4 * share, starmap.star_light_boost(lum))
        return glow * size * size

    luminosities = [10 ** (e / 50) for e in range(-250, 301)]
    sector = [_light_for(lum, radius_solar=0.1) for lum in luminosities]
    sector_light = [light["glow"] * light["sizePx"] ** 2 for light in sector]
    assert all(b >= a * 0.995 for a, b in zip(sector_light, sector_light[1:]))
    galaxy = [galaxy_light(lum) for lum in luminosities]
    assert all(b >= a for a, b in zip(galaxy, galaxy[1:]))


@pytest.mark.parametrize("luminosity", [1000.0, 5e4, 1e6])
def test_stars_of_1000_suns_and_up_look_as_before(luminosity):
    light = _light_for(luminosity, radius_solar=50.0)
    expected = _unboosted(luminosity, radius_solar=50.0)
    for key in ("sizePx", "glow", "bright"):
        assert light[key] == pytest.approx(expected[key], abs=0.011), key


def test_the_dimmest_red_dwarf_gives_four_times_the_light():
    """Light ~ halo strength x halo area: twice the strength over twice
    the area; the core is the star's own."""
    light = _light_for(1e-4, radius_solar=0.1, temperature_k=2600.0)
    before = _unboosted(1e-4, radius_solar=0.1)
    assert light["glow"] == pytest.approx(2 * before["glow"], rel=0.01)
    assert (light["sizePx"] / before["sizePx"]) ** 2 == pytest.approx(2.0, rel=0.02)
    assert light["glow"] * light["sizePx"] ** 2 == pytest.approx(4 * before["glow"] * before["sizePx"] ** 2, rel=0.03)
    assert light["bright"] == before["bright"]


def test_a_sun_is_boosted_part_way():
    light = _light_for(1.0)
    before = _unboosted(1.0)
    assert light["glow"] / before["glow"] == pytest.approx(2.2 ** 0.5, rel=0.02)


# --- MAP.82 to MAP.84: rogue planets ---------------------------------------------------

def _rogue_scene():
    phenomenon = _phenomenon(type_="rogue_planet", descriptor="terrestrial")
    html = render_map_panel(_link, 1000.0, None, None, [_make_system()], phenomena=[phenomenon])
    return html, _scene_data(html)["clouds"][0]


def test_an_unmarked_rogue_planet_is_a_faint_speck_with_no_glow():
    _html, rogue = _rogue_scene()
    light = rogue["light"]
    assert light["glow"] == 0
    assert light["bright"] <= 0.25
    assert light["corePx"] <= 2 and light["sizePx"] <= 4
    star = _scene_data(render_map_panel(_link, 1000.0, None, None, [_make_system()]))["stars"][0]["light"]
    assert light["sizePx"] < star["corePx"] * 2, "smaller than a star"


def test_a_marked_rogue_planet_is_bigger_brighter_and_glows():
    _html, rogue = _rogue_scene()
    unmarked, marked = rogue["light"], rogue["markedLight"]
    assert marked["corePx"] >= 3 * unmarked["corePx"]
    assert marked["bright"] == 1.0 and marked["glow"] > 0.3
    assert marked["sizePx"] > 4 * unmarked["sizePx"]


def test_the_mark_rogue_planets_button_starts_off_and_can_show_it():
    html, _rogue = _rogue_scene()
    button = re.search(r'<button[^>]*data-action="toggle-rogue-markers"[^>]*>', html).group(0)
    assert 'aria-pressed="false"' in button
    assert "starmap-toggle" in button


def test_the_toggle_highlight_follows_aria_pressed_in_the_stylesheet():
    css = open(os.path.join(os.path.dirname(__file__), "..", "html", "static", "style.css"), encoding="utf-8").read()
    rule = re.search(r'\.starmap-toggle\[aria-pressed="true"\]\s*\{([^}]*)\}', css)
    assert rule and "var(--accent)" in rule.group(1)
