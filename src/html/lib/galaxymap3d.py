# html/lib/galaxymap3d.py

"""
Interactive 3D Galaxy Map panel: a real perspective-camera WebGL scene
(three.js, `static/galaxymap3d.js`) a visitor can fly freely through --
this replaced an earlier flat, face-on SVG projection (a fixed camera,
whose whole-galaxy dataset was baked into one page load, and whose
fixed-radius markers grew relative to the view as you zoomed in with no
camera to shrink them the opposite way). This map's camera can travel
anywhere in the galaxy instead, so most of what it draws is fetched live
as the camera moves (`html/galaxy_view.py`, this page's own client-side
JS `fetch()` target -- see that script's own docstring) rather than
server-rendered once.

This module's job mirrors `lib/starmap.py`'s division of labor:
`galaxy.py` (the page) makes every `apiclient` call (`get_galaxy_shape`/
`get_galaxy_view`); this module only ever turns already-fetched plain
data into the panel's HTML and its one starting JSON payload -- every
*later* payload (`static/galaxymap3d.js`'s own live re-fetches as the
camera moves) never passes through this module at all.

Three content tiers, matching `queryDb.galaxy_view`'s own three lists
(`placed`/`planned`/`density` -- see `stellarObjects.galaxyViewport`'s
module docstring for what each means):

- **Placed**: real, already-generated sectors -- bright, clickable,
  sized/colored by `system_count`, navigates to `sector.py`.
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
  distance), the ceiling is this galaxy's own real outer edge (its
  stored skeleton's `outer_shell_index`, or a real Milky-Way-scale radius
  as a starting-point default when no skeleton has been built yet).
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

try:
    from stellarObjects.program_constants import GALAXY_RADIUS_PC
    from stellarObjects.utils import ly_to_pc, pc_to_ly
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # duplicated fallback, matching every other lib/ module's identical
    # pattern (see e.g. galaxymap.py's own top-of-file try/except).
    GALAXY_RADIUS_PC = 15000.0
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
    if galaxy_shape and galaxy_shape.get("outer_shell_index") is not None:
        max_radius = (galaxy_shape["outer_shell_index"] + 1) * edge_pc * MAX_VIEW_RADIUS_MARGIN
    else:
        max_radius = GALAXY_RADIUS_PC * MAX_VIEW_RADIUS_MARGIN
    return min_radius, max(max_radius, min_radius * 10)


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
    zoom-range numbers (`view_radius_bounds`) and `initial_view`'s own
    payload for the first frame -- everything after that first frame comes
    from the client's own live `fetch()` calls to `galaxy_view.py`.

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
        initial_view (dict): `apiclient.get_galaxy_view`'s own return
            shape (`placed`/`planned`/`density`/`edge_pc`/`has_shape`),
            fetched by `galaxy.py` for the galactic origin at
            `view_radius_bounds`'s own `max_view_radius_pc` -- the
            zoomed-all-the-way-out starting view.

    Returns:
        str: A complete `<section class="panel">` block.
    """
    min_radius, max_radius = view_radius_bounds(edge_pc, galaxy_shape)
    edge_ly = pc_to_ly(edge_pc)

    scene_data = {
        "db": db_name,
        "fetchPath": "galaxy_view.py",
        "edgePc": edge_pc,
        "edgeLy": edge_ly,
        "minViewRadiusPc": min_radius,
        "maxViewRadiusPc": max_radius,
        "clickZoomFactorMin": CLICK_ZOOM_FACTOR_MIN,
        "clickZoomFactorMax": CLICK_ZOOM_FACTOR_MAX,
        "initialCenter": [0.0, 0.0, 0.0],
        "initialRadiusPc": max_radius,
        "initial": initial_view,
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
  steps while zoomed out, finer near a single sector) &middot; bright dots &asymp; generated sectors &middot;
  small dim dots &asymp; real, not-yet-generated sector addresses &middot; faint spheres &asymp; illustrative
  predicted density</span>
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
