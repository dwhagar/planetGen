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
    FETCH_RADIUS_FACTOR,
    MIN_VIEW_RADIUS_FLOOR_PC,
    PLANNED_MAX_VIEW_RADIUS_PC,
    galaxy_extent_pc,
    initial_tile_request,
    render_galaxy_map3d_panel,
    view_radius_bounds,
)

from stellarObjects.galaxyViewport import (  # noqa: E402
    TILE_MAX_LEVEL,
    parse_tile_key,
    tile_level_for_view_radius,
    tiles_intersecting_sphere,
)
from stellarObjects.program_constants import GALAXY_RADIUS_PC  # noqa: E402

EDGE_PC = 3.526


def _empty_view(edge_pc=EDGE_PC, has_shape=False):
    return {"stamp": "0123456789abcdef", "tiles": {}, "density": None, "edge_pc": edge_pc, "has_shape": has_shape}


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


def test_view_radius_bounds_with_a_shape_uses_its_own_outer_shell_index():
    min_radius, max_radius = view_radius_bounds(EDGE_PC, {"outer_shell_index": 99})
    assert max_radius == pytest.approx((99 + 1) * EDGE_PC * 1.05 / math.tan(math.radians(CAMERA_FOV_DEG / 2)))


def test_view_radius_bounds_max_fits_the_whole_galaxy_in_view():
    # At the zoomed-all-the-way-out radius, half the field of view must
    # span at least the galaxy's own padded outer edge.
    shape = {"outer_shell_index": 99}
    _min_radius, max_radius = view_radius_bounds(EDGE_PC, shape)
    visible_half_height = max_radius * math.tan(math.radians(CAMERA_FOV_DEG / 2))
    assert visible_half_height >= galaxy_extent_pc(EDGE_PC, shape) - 1e-9


def test_view_radius_bounds_min_has_an_absolute_floor():
    # A pathologically tiny edge_pc must not compute a floor smaller than
    # MIN_VIEW_RADIUS_FLOOR_PC.
    min_radius, _max_radius = view_radius_bounds(0.001, None)
    assert min_radius == pytest.approx(MIN_VIEW_RADIUS_FLOOR_PC)


def test_view_radius_bounds_max_is_always_well_past_min():
    min_radius, max_radius = view_radius_bounds(EDGE_PC, {"outer_shell_index": 0})
    assert max_radius >= min_radius * 10


# --- render_galaxy_map3d_panel -----------------------------------------------

def test_panel_includes_the_canvas_and_controls():
    html = render_galaxy_map3d_panel("mydb", None, EDGE_PC, _empty_view())
    assert 'id="galaxymap3d-canvas"' in html
    assert 'data-action="zoom-in"' in html
    assert 'data-action="zoom-out"' in html
    assert 'data-action="reset"' in html
    assert 'id="galaxymap3d-info"' in html


def test_panel_json_payload_has_every_field_the_client_reads():
    view = _empty_view(has_shape=True)
    html = render_galaxy_map3d_panel("mydb", {"outer_shell_index": 50}, EDGE_PC, view,
                                     fetch_path="/galaxy/tiles", sector_url="/sector/{id}")
    data = _json_payload(html)

    assert data["storageKey"] == "mydb"
    assert "db" not in data
    assert data["fetchPath"] == "/galaxy/tiles"
    assert data["sectorUrl"] == "/sector/{id}"
    assert data["hasShape"] is True
    for field in ("tileRootEdgePc", "tileMaxLevel", "plannedTileMaxEdgePc", "plannedMaxViewRadiusPc",
                  "fetchRadiusFactor", "maxTilesPerRequest"):
        assert field in data
    assert data["edgePc"] == pytest.approx(EDGE_PC)
    assert data["edgeLy"] > 0
    assert data["clickZoomFactorMin"] == CLICK_ZOOM_FACTOR_MIN
    assert data["clickZoomFactorMax"] == CLICK_ZOOM_FACTOR_MAX
    assert data["initialCenter"] == [0.0, 0.0, 0.0]
    assert data["initialRadiusPc"] == data["maxViewRadiusPc"]
    assert data["fovDeg"] == CAMERA_FOV_DEG
    assert data["galaxyRadiusPc"] == pytest.approx(galaxy_extent_pc(EDGE_PC, {"outer_shell_index": 50}))
    assert data["initial"] == view


def test_panel_shows_a_hint_when_no_shape_has_been_built():
    html = render_galaxy_map3d_panel("mydb", None, EDGE_PC, _empty_view(has_shape=False))
    assert "density skeleton hasn&#x27;t been built yet" in html or "density skeleton hasn't been built yet" in html
    assert "generate.py plan" in html


def test_panel_omits_the_hint_when_a_shape_exists():
    html = render_galaxy_map3d_panel("mydb", {"outer_shell_index": 50}, EDGE_PC, _empty_view(has_shape=True))
    assert "density skeleton hasn't been built yet" not in html


def test_panel_json_is_safely_escaped_against_script_breakout():
    # A database name is arbitrary operator-supplied text -- confirms the
    # same </script>-breakout mitigation lib/starmap.py's own
    # _json_script uses is applied here too.
    html = render_galaxy_map3d_panel("weird</script><script>alert(1)</script>db", None, EDGE_PC, _empty_view())
    assert "</script><script>alert" not in html
    data = _json_payload(html)
    assert data["storageKey"] == "weird</script><script>alert(1)</script>db"


def test_panel_includes_real_placed_and_planned_data_from_the_initial_view():
    view = _empty_view(has_shape=True)
    view["tiles"]["1/1/1/1"] = {
        "placed": [{"id": 1, "name": "Real Sector", "x": 1.0, "y": 2.0, "z": 3.0,
                     "galactic_radius_pc": 3.7, "shell_index": 1, "shell_slot_index": 0,
                     "designation": "ABC", "system_count": 4}],
        "planned": [],
    }
    html = render_galaxy_map3d_panel("mydb", {"outer_shell_index": 10}, EDGE_PC, view)
    data = _json_payload(html)
    assert data["initial"]["tiles"]["1/1/1/1"]["placed"][0]["name"] == "Real Sector"
    assert data["initial"]["stamp"] == "0123456789abcdef"


# --- initial_tile_request ----------------------------------------------------

@pytest.mark.parametrize("orbit_radius", [16000.0, 900.0, 120.0, 30.0, 5.3])
def test_initial_tile_request_matches_the_shared_tile_math(orbit_radius):
    center = (321.0, -45.0, 6.0)
    keys, density_key = initial_tile_request(orbit_radius, True, center)
    view_radius = orbit_radius * FETCH_RADIUS_FACTOR
    level = tile_level_for_view_radius(view_radius)
    view_keys = tiles_intersecting_sphere(level, center, view_radius)
    assert keys[:len(view_keys)] == view_keys
    for key in keys:
        parse_tile_key(key)
    if view_radius <= PLANNED_MAX_VIEW_RADIUS_PC:
        assert any(key.startswith(f"{TILE_MAX_LEVEL}/") for key in keys)
        assert density_key is None
    else:
        assert all(key.startswith(f"{level}/") for key in keys)
        assert density_key is not None and density_key.startswith(f"{level}/")
    assert len(keys) <= 128


def test_initial_tile_request_has_no_density_without_a_shape():
    _keys, density_key = initial_tile_request(16000.0, False)
    assert density_key is None


def test_planned_tiles_cover_the_whole_view_not_a_smaller_ball():
    # A zoomed-in view fetches planned tiles out to the full view radius,
    # so the dots fill the screen instead of clustering into a ball around
    # the target.
    center = (8000.0, 3.0, -2.0)
    orbit_radius = PLANNED_MAX_VIEW_RADIUS_PC / FETCH_RADIUS_FACTOR
    keys, density_key = initial_tile_request(orbit_radius, True, center)
    view_radius = orbit_radius * FETCH_RADIUS_FACTOR
    planned_keys = tiles_intersecting_sphere(TILE_MAX_LEVEL, center, view_radius)
    assert set(planned_keys) <= set(keys)
    assert density_key is None
