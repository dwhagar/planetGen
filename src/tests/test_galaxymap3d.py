"""
html/lib/galaxymap3d.py regression tests -- the interactive 3D Galaxy
Map's server-side panel builder (the `/galaxy` view, `web/galaxy_views.py`, calls `view_radius_bounds`
to pick its starting radius, `initial_tile_request` for the first
frame's tiles, then `render_galaxy_map3d_panel` to build the page). No
database needed -- every function takes plain data, the same shape
`apiclient.get_galaxy_shape`/`tilecache.fetch_tiles` return.

Run with: pytest src/tests/test_galaxymap3d.py
"""
import json
import math
import os
import sys

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_SRC_DIR, "html", "lib"))
sys.path.insert(0, _SRC_DIR)

import pytest  # noqa: E402

from galaxymap3d import (  # noqa: E402
    CAMERA_FOV_DEG,
    CLICK_ZOOM_FACTOR_MAX,
    CLICK_ZOOM_FACTOR_MIN,
    DENSITY_SHAPE_FIELDS,
    FETCH_RADIUS_FACTOR,
    MIN_VIEW_RADIUS_FLOOR_PC,
    galaxy_extent_pc,
    initial_tile_request,
    render_galaxy_map3d_panel,
    view_radius_bounds,
)

from stellarObjects.galaxyViewport import (  # noqa: E402
    parse_tile_key,
    tile_level_for_view_radius,
    tiles_intersecting_sphere,
)
from stellarObjects.program_constants import GALAXY_RADIUS_PC  # noqa: E402

EDGE_PC = 4.0


def _empty_view(edge_pc=EDGE_PC, has_shape=False):
    return {"stamp": "0123456789abcdef", "tiles": {}, "edge_pc": edge_pc, "has_shape": has_shape}


def _json_payload(html):
    marker = '<script type="application/json" id="galaxymap3d-data">'
    start = html.index(marker) + len(marker)
    end = html.index("</script>", start)
    return json.loads(html[start:end])


# --- view_radius_bounds -----------------------------------------------------

def test_view_radius_bounds_without_a_shape_uses_the_default_galaxy_radius():
    min_radius, max_radius = view_radius_bounds(EDGE_PC, None)
    assert max_radius == pytest.approx(GALAXY_RADIUS_PC * 1.05 / math.tan(math.radians(CAMERA_FOV_DEG / 2)))
    assert min_radius == pytest.approx(EDGE_PC * 1.5)


def test_view_radius_bounds_with_a_shape_uses_its_own_outer_ring_index():
    min_radius, max_radius = view_radius_bounds(EDGE_PC, {"outer_ring_index": 99})
    assert max_radius == pytest.approx((99 + 1) * EDGE_PC * 1.05 / math.tan(math.radians(CAMERA_FOV_DEG / 2)))


def test_view_radius_bounds_max_fits_the_whole_galaxy_in_view():
    # At the zoomed-all-the-way-out radius, half the field of view must
    # span at least the galaxy's own padded outer edge.
    shape = {"outer_ring_index": 99}
    _min_radius, max_radius = view_radius_bounds(EDGE_PC, shape)
    visible_half_height = max_radius * math.tan(math.radians(CAMERA_FOV_DEG / 2))
    assert visible_half_height >= galaxy_extent_pc(EDGE_PC, shape) - 1e-9


def test_view_radius_bounds_min_has_an_absolute_floor():
    # A pathologically tiny edge_pc must not compute a floor smaller than
    # MIN_VIEW_RADIUS_FLOOR_PC.
    min_radius, _max_radius = view_radius_bounds(0.001, None)
    assert min_radius == pytest.approx(MIN_VIEW_RADIUS_FLOOR_PC)


def test_view_radius_bounds_max_is_always_well_past_min():
    min_radius, max_radius = view_radius_bounds(EDGE_PC, {"outer_ring_index": 0})
    assert max_radius >= min_radius * 10


# --- render_galaxy_map3d_panel -----------------------------------------------

def test_panel_includes_the_canvas_and_controls():
    html = render_galaxy_map3d_panel("mydb", None, EDGE_PC, _empty_view())
    assert 'id="galaxymap3d-canvas"' in html
    # Seen from above only (MAP.17): Back, Forward, Up and Whole galaxy,
    # no zoom buttons and no free camera.
    for action in ("back", "forward", "up", "reset", "reset-view"):
        assert f'data-action="{action}"' in html
    for gone in ("zoom-in", "zoom-out", "free-look", "slice"):
        assert f'data-action="{gone}"' not in html
    assert 'id="galaxymap3d-info"' in html
    # The drill-down's breadcrumb, slab strip, tooltip and toggles.
    for element in ("galaxymap3d-crumbs", "galaxymap3d-slabs", "galaxymap3d-tooltip", "galaxymap3d-notice"):
        assert f'id="{element}"' in html
    assert 'data-action="territories"' in html
    assert 'id="galaxymap3d-territories"' in html
    assert 'data-action="generated-only"' in html


def test_panel_leaves_territories_out_when_there_are_none():
    """No polities yet: no Territories button, no legend, no endpoint."""
    html = render_galaxy_map3d_panel("mydb", None, EDGE_PC, _empty_view(), territory_path=None)
    assert 'data-action="territories"' not in html
    assert 'id="galaxymap3d-territories"' not in html
    assert _json_payload(html)["territoryPath"] is None


def test_panel_json_payload_has_every_field_the_client_reads():
    view = _empty_view(has_shape=True)
    html = render_galaxy_map3d_panel("mydb", {"outer_ring_index": 50}, EDGE_PC, view,
                                     fetch_path="/galaxy/tiles", sector_url="/sector/{id}")
    data = _json_payload(html)

    assert data["storageKey"] == "mydb"
    assert "db" not in data
    assert data["fetchPath"] == "/galaxy/tiles"
    assert data["stagePath"] == "/galaxy/stage"
    assert data["territoryPath"] == "/galaxy/territories"
    assert data["sectorUrl"] == "/sector/{id}"
    assert data["hasShape"] is True
    for field in ("tileRootEdgePc", "tileMaxLevel", "fetchRadiusFactor", "maxTilesPerRequest",
                  "blockMinPx", "blockBudget"):
        assert field in data
    assert data["edgePc"] == pytest.approx(EDGE_PC)
    assert data["edgeLy"] > 0
    assert data["clickZoomFactorMin"] == CLICK_ZOOM_FACTOR_MIN
    assert data["clickZoomFactorMax"] == CLICK_ZOOM_FACTOR_MAX
    assert data["initialCenter"] == [0.0, 0.0, 0.0]
    assert data["initialRadiusPc"] == data["maxViewRadiusPc"]
    assert data["fovDeg"] == CAMERA_FOV_DEG
    assert data["galaxyRadiusPc"] == pytest.approx(galaxy_extent_pc(EDGE_PC, {"outer_ring_index": 50}))
    # The wedge lines stop at the outermost ring's outside, not the padded view radius (MAP.43).
    assert data["galaxyEdgePc"] == pytest.approx(51 * EDGE_PC)
    assert data["galaxyEdgePc"] < data["galaxyRadiusPc"]
    assert data["initial"] == view
    assert data["generate"] is None
    assert data["phenomenonUrl"] is None
    assert data["systemUrl"] is None


def test_panel_passes_the_admin_generate_target_through():
    target = {"url": "/admin/generate", "csrfField": "csrf_token", "csrfToken": "abc"}
    data = _json_payload(render_galaxy_map3d_panel("mydb", None, EDGE_PC, _empty_view(), generate=target))
    assert data["generate"] == target


def test_panel_shows_a_hint_when_no_shape_has_been_built():
    html = render_galaxy_map3d_panel("mydb", None, EDGE_PC, _empty_view(has_shape=False))
    assert "density skeleton hasn&#x27;t been built yet" in html or "density skeleton hasn't been built yet" in html
    assert "generate.py plan" in html


def test_panel_omits_the_hint_when_a_shape_exists():
    html = render_galaxy_map3d_panel("mydb", {"outer_ring_index": 50}, EDGE_PC, _empty_view(has_shape=True))
    assert "density skeleton hasn't been built yet" not in html


def test_panel_embeds_the_density_shape_for_the_prisms():
    shape = {field: float(i + 1) for i, field in enumerate(DENSITY_SHAPE_FIELDS)}
    shape.update({"outer_ring_index": 50, "edge_pc": EDGE_PC, "expected_system_count_at_density_1": 9.0})
    data = _json_payload(render_galaxy_map3d_panel("mydb", shape, EDGE_PC, _empty_view(has_shape=True)))
    assert data["densityShape"] == dict({field: shape[field] for field in DENSITY_SHAPE_FIELDS},
                                        sector_min_density=pytest.approx(1 / 9.0))


def test_panel_has_no_density_shape_without_a_skeleton():
    data = _json_payload(render_galaxy_map3d_panel("mydb", None, EDGE_PC, _empty_view()))
    assert data["densityShape"] is None


def test_panel_json_is_safely_escaped_against_script_breakout():
    # A database name is arbitrary operator-supplied text -- confirms the
    # same </script>-breakout mitigation lib/starmap.py's own
    # _json_script uses is applied here too.
    html = render_galaxy_map3d_panel("weird</script><script>alert(1)</script>db", None, EDGE_PC, _empty_view())
    assert "</script><script>alert" not in html
    data = _json_payload(html)
    assert data["storageKey"] == "weird</script><script>alert(1)</script>db"


def test_panel_includes_the_initial_views_tiles():
    view = _empty_view(has_shape=True)
    view["tiles"]["1/1/1/1"] = {
        "placed": [{"id": 1, "name": "Real Sector", "x": 1.0, "y": 2.0, "z": 3.0,
                     "galactic_radius_pc": 3.7, "ring_index": 1, "layer_index": 0, "ring_slot_index": 0,
                     "designation": "ABC", "system_count": 4}],
        "planned": [],
        "filled": {"g": 1, "cells": [[1, 0, 0, 1, 4, "Real Sector"]]},
    }
    html = render_galaxy_map3d_panel("mydb", {"outer_ring_index": 10}, EDGE_PC, view)
    data = _json_payload(html)
    assert data["initial"]["tiles"]["1/1/1/1"]["placed"][0]["name"] == "Real Sector"
    assert data["initial"]["tiles"]["1/1/1/1"]["filled"]["cells"][0][5] == "Real Sector"
    assert data["initial"]["stamp"] == "0123456789abcdef"


# --- initial_tile_request ----------------------------------------------------

@pytest.mark.parametrize("orbit_radius", [16000.0, 900.0, 120.0, 30.0, 5.3])
def test_initial_tile_request_matches_the_shared_tile_math(orbit_radius):
    center = (321.0, -45.0, 6.0)
    keys = initial_tile_request(orbit_radius, center)
    view_radius = orbit_radius * FETCH_RADIUS_FACTOR
    level = tile_level_for_view_radius(view_radius)
    view_keys = tiles_intersecting_sphere(level, center, view_radius)
    # One level only: the map no longer fetches the finest (planned) tiles.
    assert keys == view_keys
    for key in keys:
        parse_tile_key(key)
    assert len(keys) <= 128

