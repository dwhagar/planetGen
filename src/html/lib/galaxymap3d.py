# html/lib/galaxymap3d.py

"""
Interactive 3D Galaxy Map panel: a real perspective-camera WebGL scene
(three.js, `static/galaxymap3d.js`) a visitor can fly freely through --
this replaced an earlier flat, face-on SVG projection (a fixed camera,
whose whole-galaxy dataset was baked into one page load, and whose
fixed-radius markers grew relative to the view as you zoomed in with no
camera to shrink them the opposite way). This map's camera can travel
anywhere in the galaxy instead, so most of what it draws is fetched live
as the camera moves, one fixed cube of space ("tile") at a time
(`html/galaxy_tiles.py`, this page's own client-side JS `fetch()` target
-- see that script's own docstring, and `stellarObjects.galaxyViewport`'s
"Cube tiles" section) rather than server-rendered once. Tiles are cached
on the server's disk (`lib/tilecache.py`) and in the visitor's browser.

This module's job mirrors `lib/starmap.py`'s division of labor:
`galaxy.py` (the page) fetches the first frame's tiles
(`initial_tile_request` says which); this module only ever turns
already-fetched plain data into the panel's HTML and its one starting
JSON payload -- every *later* payload (`static/galaxymap3d.js`'s own live
tile fetches as the camera moves) never passes through this module at
all.

Three content tiers, matching each tile's `placed`/`planned` lists and
the separate density cloud (see `stellarObjects.galaxyViewport`'s module
docstring for what each means):

- **Placed**: real, already-generated sectors -- bright, clickable, sized
  by `system_count`, colored by real stellar density (`system_count /
  edge_ly ** 3`, relative to `physical_constants.LOCAL_STELLAR_DENSITY_LY3`
  -- see `render_galaxy_map3d_panel`'s own `referenceDensityPerLy3`),
  navigates to `sector.py`.
- **Planned**: real, not-yet-generated sector addresses this galaxy's own
  density model (when built) predicts would qualify -- small, dim,
  clickable, shows its designation/address (copyable straight into
  `generate.py galaxy --shell K --slot N`) rather than navigating
  anywhere (there's nothing to navigate to yet).
- **Density**: a coarse illustrative point cloud, for whatever part of
  the current view is too wide to enumerate individual planned
  addresses -- not clickable, carries no info of its own.

Unlike `lib/starmap.py` (every dot's size/color/position is computed
once, server-side, and the client only ever draws exactly what it's
handed), per-dot *styling* (radius/opacity/color
formulas) for these three tiers lives entirely in `static/galaxymap3d.js`
instead: the overwhelming majority of what gets drawn arrives through
this page's own live `fetch()` re-queries as the camera moves, which
never pass through this Python module again after the first paint, so a
styling formula written here would only ever apply to one initial frame
and silently diverge from whatever the client re-derives for every frame
after it. Keeping the one formula client-side (applied identically to
this module's own embedded starting payload and to every live re-fetch)
avoids that split -- this module still computes every *position*
(`sectors.center_x/y/z_pc`, already real galaxy-frame parsecs, passed
through as-is) and every *zoom-range number* (below), since those don't
change meaning between the first paint and a later live fetch.

Zoom range and the click-to-zoom step size are computed here, not
hardcoded client-side, matching the rest of this project's convention
that a client script visualizes numbers a `lib/` module already worked
out:

- `min_view_radius_pc`/`max_view_radius_pc` bound how far the camera can
  dolly -- the floor sits just past a couple of sector-widths (so
  approaching one specific sector never has to overshoot past clicking
  distance), the ceiling is far enough back that this galaxy's own real
  outer edge (its stored skeleton's `outer_shell_index`, or a real
  Milky-Way-scale radius as a starting-point default when no skeleton has
  been built yet) fits inside the camera's field of view.
- Left/right-click zoom is **logarithmic**, not a flat step: the closer
  the camera already is, the smaller each click's own multiplicative jump
  gets. `static/galaxymap3d.js`'s own `clickZoomFactor` interpolates
  between `click_zoom_factor_min` (near `min_view_radius_pc`) and
  `click_zoom_factor_max` (near `max_view_radius_pc`) by the camera's
  *current* distance in log space, recomputed fresh on every click -- so
  a handful of clicks from the full-galaxy view still closes most of the
  distance (big steps while far out), while a click near one sector
  nudges in gently instead of blowing past it (small steps while close
  in), rather than a fixed "2x every click" factor that's either too slow
  to cross the galaxy or too coarse to land on one sector.
"""

import json

import math

try:
    from stellarObjects.galaxyViewport import (
        PLANNED_TILE_MAX_EDGE_PC,
        TILE_MAX_LEVEL,
        TILE_ROOT_EDGE_PC,
    )
    from stellarObjects.physical_constants import LOCAL_STELLAR_DENSITY_LY3
    from stellarObjects.program_constants import GALAXY_RADIUS_PC
    from stellarObjects.utils import ly_to_pc, pc_to_ly
except ImportError:
    TILE_ROOT_EDGE_PC = 65536.0
    TILE_MAX_LEVEL = 12
    PLANNED_TILE_MAX_EDGE_PC = 16.0
    # The planetGen package isn't on the import path in this deployment --
    # duplicated fallback, matching every other lib/ module's identical
    # pattern (see e.g. galaxymap.py's own top-of-file try/except).
    GALAXY_RADIUS_PC = 15000.0
    LOCAL_STELLAR_DENSITY_LY3 = 0.00284
    def ly_to_pc(ly):
        return ly / 3.2616
    def pc_to_ly(pc):
        return pc * 3.2616

DEFAULT_SECTOR_EDGE_LY_FALLBACK = 11.5
"""float: Used only if `stellarObjects.program_constants` itself isn't
importable (see the top-of-file fallback above) -- matches that module's
own `DEFAULT_SECTOR_EDGE_LY`."""

CLICK_ZOOM_FACTOR_MIN = 1.15
"""float: The click-to-zoom multiplicative step size once the camera is
already at (or near) `min_view_radius_pc` -- fine enough to approach one
specific sector without overshooting past it. See the module docstring's
own "logarithmic zoom" explanation."""

CLICK_ZOOM_FACTOR_MAX = 4.0
"""float: The click-to-zoom step size at (or near) `max_view_radius_pc`
(the full-galaxy view) -- large enough that a handful of clicks closes
most of a real galaxy's own huge dynamic range, instead of the many more
clicks a flat, non-adaptive factor would need from that far out."""

MIN_VIEW_RADIUS_FLOOR_PC = 2.0
"""float: Absolute floor under `view_radius_bounds`'s own
edge-length-scaled minimum, in case a pathologically small `edge_pc`
would otherwise compute a smaller one."""

MAX_VIEW_RADIUS_MARGIN = 1.05
"""float: `view_radius_bounds` pads the galaxy's own real outer edge by
this factor for the zoomed-all-the-way-out ceiling -- a hair of headroom
so the outermost real content isn't sitting exactly on the view's own
edge."""

CAMERA_FOV_DEG = 50.0
"""float: The map camera's vertical field of view, degrees (passed to
`static/galaxymap3d.js` as `fovDeg`). `view_radius_bounds` needs it to
back the camera off far enough that the whole galaxy fits in the square
viewport at the zoomed-all-the-way-out view."""


def view_radius_bounds(edge_pc, galaxy_shape):
    """
    `(min_view_radius_pc, max_view_radius_pc)` -- see the module
    docstring's own explanation of what these bound. A pure function of
    already-fetched data (no I/O), so `galaxy.py` (the page) can call it
    directly to pick the radius its own first `get_galaxy_view` call uses,
    before this module's own panel-rendering function ever runs.

    Args:
        edge_pc (float): The sector edge length, parsecs.
        galaxy_shape (dict or None): `apiclient.get_galaxy_shape`'s own
            return shape, or `None` if `generate.py plan` has never been
            run against this database.

    Returns:
        tuple[float, float]: `(min_view_radius_pc, max_view_radius_pc)`.
    """
    min_radius = max(MIN_VIEW_RADIUS_FLOOR_PC, edge_pc * 1.5)
    max_radius = galaxy_extent_pc(edge_pc, galaxy_shape) / math.tan(math.radians(CAMERA_FOV_DEG / 2))
    return min_radius, max(max_radius, min_radius * 10)


def galaxy_extent_pc(edge_pc, galaxy_shape):
    """
    The galaxy's own outer edge, parsecs, padded by
    `MAX_VIEW_RADIUS_MARGIN`: its stored skeleton's `outer_shell_index`,
    or `GALAXY_RADIUS_PC` when no skeleton has been built yet. The camera
    target is kept inside this radius (`galaxyRadiusPc`), and
    `view_radius_bounds` backs the camera off far enough to fit it.

    Args:
        edge_pc (float): The sector edge length, parsecs.
        galaxy_shape (dict or None): `apiclient.get_galaxy_shape`'s own
            return shape, or `None`.

    Returns:
        float: The padded outer radius, parsecs.
    """
    if galaxy_shape and galaxy_shape.get("outer_shell_index") is not None:
        return (galaxy_shape["outer_shell_index"] + 1) * edge_pc * MAX_VIEW_RADIUS_MARGIN
    return GALAXY_RADIUS_PC * MAX_VIEW_RADIUS_MARGIN


FETCH_RADIUS_FACTOR = 1.6
"""float: The map fetches content out to this multiple of the camera's
orbit radius, so what's just off-screen is already there when the camera
turns."""

PLANNED_VIEW_RADIUS_PC = 20.0
"""float: Planned (not-yet-generated) slots are shown out to this far from
the camera target -- a sphere of ~400 slots at the default sector size,
from up to ~64 of the smallest (`PLANNED_TILE_MAX_EDGE_PC`) tiles."""

PLANNED_MAX_VIEW_RADIUS_PC = 200.0
"""float: Planned slots are only fetched while the view radius is at most
this; the density cloud takes over for wider views (the same 200 pc
switch-over the map has always had)."""

MAX_TILES_PER_REQUEST = 128
"""int: Mirrors `queryDb.MAX_TILES_PER_REQUEST`."""


def _tile_level_for_view_radius(radius_pc):
    """Same as `galaxyViewport.tile_level_for_view_radius` (duplicated
    so this module keeps working under its ImportError fallback)."""
    if radius_pc <= 0:
        return TILE_MAX_LEVEL
    return max(0, min(TILE_MAX_LEVEL, math.floor(math.log2(TILE_ROOT_EDGE_PC / radius_pc))))


def _tiles_intersecting_sphere(level, center_pc, radius_pc):
    """Same as `galaxyViewport.tiles_intersecting_sphere`, duplicated for
    the same reason as `_tile_level_for_view_radius`."""
    edge = TILE_ROOT_EDGE_PC / (2 ** level)
    origin = -TILE_ROOT_EDGE_PC / 2.0
    span = 2 ** level
    ranges = []
    for axis in range(3):
        first = math.floor((center_pc[axis] - radius_pc - origin) / edge)
        last = math.floor((center_pc[axis] + radius_pc - origin) / edge)
        ranges.append(range(max(0, first), min(span - 1, last) + 1))
    found = []
    for ix in ranges[0]:
        for iy in ranges[1]:
            for iz in ranges[2]:
                index = (ix, iy, iz)
                distance_sq = 0.0
                for axis in range(3):
                    lo = origin + index[axis] * edge
                    nearest = min(max(center_pc[axis], lo), lo + edge)
                    distance_sq += (center_pc[axis] - nearest) ** 2
                if distance_sq <= radius_pc * radius_pc:
                    found.append((distance_sq, f"{level}/{ix}/{iy}/{iz}"))
    found.sort()
    return [key for _distance, key in found]


def _tile_containing(level, point_pc):
    edge = TILE_ROOT_EDGE_PC / (2 ** level)
    origin = -TILE_ROOT_EDGE_PC / 2.0
    span = 2 ** level
    index = [max(0, min(span - 1, math.floor((point_pc[axis] - origin) / edge))) for axis in range(3)]
    return f"{level}/{index[0]}/{index[1]}/{index[2]}"


def initial_tile_request(orbit_radius_pc, has_shape, center_pc=(0.0, 0.0, 0.0)):
    """
    The tiles (and density anchor) the map's first frame needs, computed
    exactly the way `static/galaxymap3d.js`'s own `neededTiles` does for
    every later camera position, so the browser's first live fetch finds
    the first frame's tiles already cached.

    Args:
        orbit_radius_pc (float): The starting camera orbit radius
            (`view_radius_bounds`' max).
        has_shape (bool): Whether the galaxy has a density skeleton (no
            density cloud without one).
        center_pc (tuple): The starting camera target.

    Returns:
        tuple: `(tile_keys, density_key)` -- `density_key` is `None` when
            no density cloud is wanted.
    """
    view_radius = orbit_radius_pc * FETCH_RADIUS_FACTOR
    level = _tile_level_for_view_radius(view_radius)
    keys = _tiles_intersecting_sphere(level, center_pc, view_radius)
    if view_radius <= PLANNED_MAX_VIEW_RADIUS_PC:
        planned_level = _tile_level_for_view_radius(PLANNED_TILE_MAX_EDGE_PC)
        planned_radius = min(view_radius, PLANNED_VIEW_RADIUS_PC)
        keys += [key for key in _tiles_intersecting_sphere(planned_level, center_pc, planned_radius) if key not in keys]
    density_key = None
    if has_shape and view_radius > PLANNED_MAX_VIEW_RADIUS_PC:
        density_key = _tile_containing(level, center_pc)
    return keys, density_key

def _json_script(data):
    """Same `<script type="application/json">`-safe escaping
    `lib/starmap.py`'s own `_json_script` uses -- see that function's
    docstring for why (a database name/label can contain `</script>`)."""
    return (
        json.dumps(data)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
    )


def render_galaxy_map3d_panel(db_name, galaxy_shape, edge_pc, initial_view):
    """
    Builds the "Galaxy Map (3D)" panel: a `<canvas>` `static/
    galaxymap3d.js` renders an interactive WebGL scene into (drag to
    rotate, scroll or the zoom buttons to zoom, click to center/select,
    double-click to center/select AND zoom in -- no right-click action),
    plus a `<script type="application/json">` block carrying the
    zoom-range numbers (`view_radius_bounds`), the tile settings the
    client needs to pick tiles the same way `initial_tile_request` does,
    and `initial_view`'s own payload for the first frame -- everything
    after that first frame comes from the client's own live `fetch()`
    calls to `galaxy_tiles.py`.

    Args:
        db_name (str): The current `?db=` value -- carried in the JSON
                       payload so the client's own fetch calls (and any
                       navigation to a clicked placed sector's `sector.py`)
                       can build their URLs/params without needing it
                       threaded through separately.
        galaxy_shape (dict or None): `apiclient.get_galaxy_shape`'s own
            return shape, or `None` if `generate.py plan` has never been
            run -- when `None`, the panel shows a hint that planned-sector
            qualification/density shading isn't real yet (see
            `queryDb.galaxy_view`'s own `has_shape` field, which
            `initial_view` already carries through).
        edge_pc (float): The sector edge length, parsecs
                         (`initial_view["edge_pc"]`, passed separately
                         since `view_radius_bounds` needs it before
                         `initial_view` itself is fetched).
        initial_view (dict): `tilecache.fetch_tiles`' own return shape
            (`stamp`/`tiles`/`density`/`edge_pc`/`has_shape`), fetched by
            `galaxy.py` for `initial_tile_request`'s tiles -- the
            zoomed-all-the-way-out starting view. The client seeds its
            tile cache with it.

    Returns:
        str: A complete `<section class="panel">` block.
    """
    min_radius, max_radius = view_radius_bounds(edge_pc, galaxy_shape)
    edge_ly = pc_to_ly(edge_pc)

    scene_data = {
        "db": db_name,
        "fetchPath": "galaxy_tiles.py",
        "tileRootEdgePc": TILE_ROOT_EDGE_PC,
        "tileMaxLevel": TILE_MAX_LEVEL,
        "plannedTileMaxEdgePc": PLANNED_TILE_MAX_EDGE_PC,
        "plannedViewRadiusPc": PLANNED_VIEW_RADIUS_PC,
        "plannedMaxViewRadiusPc": PLANNED_MAX_VIEW_RADIUS_PC,
        "fetchRadiusFactor": FETCH_RADIUS_FACTOR,
        "maxTilesPerRequest": MAX_TILES_PER_REQUEST,
        "hasShape": bool(initial_view.get("has_shape")),
        "edgePc": edge_pc,
        "edgeLy": edge_ly,
        "fovDeg": CAMERA_FOV_DEG,
        "galaxyRadiusPc": galaxy_extent_pc(edge_pc, galaxy_shape),
        "minViewRadiusPc": min_radius,
        "maxViewRadiusPc": max_radius,
        "clickZoomFactorMin": CLICK_ZOOM_FACTOR_MIN,
        "clickZoomFactorMax": CLICK_ZOOM_FACTOR_MAX,
        "initialCenter": [0.0, 0.0, 0.0],
        "initialRadiusPc": max_radius,
        "initial": initial_view,
        # The real, sampled-in-the-solar-neighborhood average this
        # project's own generation already calibrates against (see
        # physical_constants.LOCAL_STELLAR_DENSITY_LY3's own citations) --
        # static/galaxymap3d.js colors each placed sector by its own real
        # density (system_count / edge_ly ** 3) RELATIVE to this same
        # reference, so "denser than real average" / "sparser than real
        # average" means the same thing on this map that it does in the
        # generator itself, not an arbitrary client-side scale.
        "referenceDensityPerLy3": LOCAL_STELLAR_DENSITY_LY3,
    }

    shape_hint = (
        ""
        if galaxy_shape
        else (
            '<p class="hint">The galaxy\'s density skeleton hasn\'t been built yet '
            "(<code>generate.py plan</code>) -- every enumerable sector address is shown as "
            "\"planned\" regardless of predicted density, and no illustrative density cloud is "
            "shown for the wider view.</p>"
        )
    )

    return f"""
<section class="panel">
<div class="panel-header">
  <h2>Galaxy Map (3D)</h2>
  <span class="hint">Drag to rotate &middot; scroll or the +/&minus; buttons to zoom &middot; click a dot or
  empty space to center the view there and select it &middot; double-click to do the same AND zoom in (bigger
  steps while zoomed out, finer near a single sector) &middot; generated-sector dots colored by their own real
  stellar density (dim &rarr; bright, relative to the real local average) &middot; small dim dots &asymp; real,
  not-yet-generated sector addresses &middot; faint spheres &asymp; illustrative predicted density</span>
</div>
{shape_hint}
<div class="starmap-layout">
<div class="starmap-viewport">
<canvas id="galaxymap3d-canvas" class="starmap-canvas" tabindex="0" role="application"
     aria-label="Interactive 3D Galaxy Map. Drag or use arrow keys to rotate, scroll or the zoom buttons to
     zoom, click a dot or empty space to center the view there and select it, double-click to do the same and
     zoom in."></canvas>
<div class="starmap-scale" id="galaxymap3d-scale"></div>
</div>
<div class="starmap-side">
<div class="starmap-controls" id="galaxymap3d-controls">
  <button type="button" class="starmap-btn" data-action="zoom-out" aria-label="Zoom out">&minus;</button>
  <button type="button" class="starmap-btn" data-action="zoom-in" aria-label="Zoom in">+</button>
  <button type="button" class="starmap-btn" data-action="reset">Reset view</button>
</div>
<aside class="starmap-info" id="galaxymap3d-info">
<p class="hint">Click a sector dot for details, or double-click a dot or empty space to zoom in there.</p>
</aside>
</div>
</div>
<script type="application/json" id="galaxymap3d-data">{_json_script(scene_data)}</script>
</section>
"""
