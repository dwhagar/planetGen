# planetgen/web/maps/galaxymap3d.py

"""
Interactive 3D Galaxy Map panel: a real perspective-camera WebGL scene
(three.js, `static/galaxymap3d.js`) a visitor can fly freely through --
this replaced an earlier flat, face-on SVG projection (a fixed camera,
whose whole-galaxy dataset was baked into one page load, and whose
fixed-radius markers grew relative to the view as you zoomed in with no
camera to shrink them the opposite way). This map's camera can travel
anywhere in the galaxy instead, so most of what it draws is fetched live
as the camera moves, one fixed cube of space ("tile") at a time
(`/galaxy/tiles`, `planetgen/web/galaxy_views.py`, this page's own client-side
JS `fetch()` target -- see that view's own docstring, and `planetgen.galaxy.viewport`'s
"Cube tiles" section) rather than server-rendered once. Tiles are cached
on the server's disk (`planetgen/web/lib/tilecache.py`) and in the visitor's browser.

This module's job mirrors `planetgen/web/maps/starmap.py`'s division of labor:
`/galaxy` (`planetgen/web/galaxy_views.py`) fetches the first frame's tiles
(`initial_tile_request` says which); this module only ever turns
already-fetched plain data into the panel's HTML and its one starting
JSON payload -- every *later* payload (`static/galaxymap3d.js`'s own live
tile fetches as the camera moves) never passes through this module at
all.

One solid of blocks shows everything (`static/galaxyprisms.js`): each
block is a power-of-3 cube of whole sectors, sized to the view and shaded
by the galaxy's predicted density, which the browser computes itself from
the shape parameters the panel embeds as `densityShape`. Unfilled space is
translucent; a block grows more solid with its share of generated
("filled") sectors, counted from each tile's `filled` summary
(`queryDb.galaxy_filled_in_box`), and is fully solid once every sector in
it is generated. At one sector per block a filled sector takes its real
stellar density's color (`system_count / edge_ly ** 3`, relative to
`physical_constants.LOCAL_STELLAR_DENSITY_LY3`, see
`render_galaxy_map3d_panel`'s `referenceDensityPerLy3`) and links to its
sector page; an unfilled one shows its designation and address, plus, for
a logged-in admin, the same Generate buttons the Sector Map gives a
neighbor that isn't generated yet (`static/generatebuttons.js`).

Unlike `planetgen/web/maps/starmap.py` (every dot's size/color/position is computed
once, server-side, and the client only ever draws exactly what it's
handed), per-block *styling* (opacity/color
formulas) lives entirely in `static/galaxymap3d.js`
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
  outer edge (its stored skeleton's `outer_ring_index`, or a real
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
    from planetgen.galaxy.viewport import (
        TILE_MAX_LEVEL,
        TILE_ROOT_EDGE_PC,
    )
    from planetgen.physics.constants import LOCAL_STELLAR_DENSITY_LY3
    from planetgen.tuning import GALAXY_RADIUS_PC
    from planetgen.physics.units import ly_to_pc, pc_to_ly
except ImportError:
    TILE_ROOT_EDGE_PC = 65536.0
    TILE_MAX_LEVEL = 12
    # The planetGen package isn't on the import path in this deployment --
    # duplicated fallback, matching every other lib/ module's identical
    # pattern (see e.g. galaxymap.py's own top-of-file try/except).
    GALAXY_RADIUS_PC = 15000.0
    LOCAL_STELLAR_DENSITY_LY3 = 0.00284
    def ly_to_pc(ly):
        return ly / 3.2616
    def pc_to_ly(pc):
        return pc * 3.2616

DEFAULT_SECTOR_EDGE_LY_FALLBACK = 13.046  # 4 pc
"""float: Used only if `planetgen.tuning` itself isn't
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

BLOCK_MIN_PX = 4
"""int: The smallest a mega-block is allowed to be on screen, CSS pixels
across at the view's focus (passed to `static/galaxymap3d.js` as
`blockMinPx`). `galaxyprisms.blockSizeForScale` picks the smallest power
of 3 sectors a side that is at least this wide, so one block always covers
at least a pixel's worth of sectors."""

BLOCK_BUDGET = 60000
"""int: The most blocks one view draws (`blockBudget`). A view whose
surface blocks would overflow it uses blocks three times bigger."""

CAMERA_FOV_DEG = 50.0
"""float: The map camera's vertical field of view, degrees (passed to
`static/galaxymap3d.js` as `fovDeg`). `view_radius_bounds` needs it to
back the camera off far enough that the whole galaxy fits in the square
viewport at the zoomed-all-the-way-out view."""


def view_radius_bounds(edge_pc, galaxy_shape):
    """
    `(min_view_radius_pc, max_view_radius_pc)` -- see the module
    docstring's own explanation of what these bound. A pure function of
    already-fetched data (no I/O), so the `/galaxy` page view can call it
    directly to pick the radius its own first tile request uses,
    before this module's own panel-rendering function ever runs.

    Args:
        edge_pc (float): The sector edge length, parsecs.
        galaxy_shape (dict or None): `apiclient.get_galaxy_shape`'s own
            return shape, or `None` if `planetgen plan` has never been
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
    `MAX_VIEW_RADIUS_MARGIN`: its stored skeleton's `outer_ring_index`,
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
    return galaxy_edge_pc(edge_pc, galaxy_shape) * MAX_VIEW_RADIUS_MARGIN


def galaxy_edge_pc(edge_pc, galaxy_shape):
    """
    The galaxy's own outer edge, parsecs, unpadded: the outside of its
    stored skeleton's outermost ring (`outer_ring_index`), or
    `GALAXY_RADIUS_PC` when no skeleton has been built yet. Every ring is
    a full circle, so this is the edge at every bearing; the wedge lines
    stop here (`galaxyEdgePc`, MAP.43).

    Args:
        edge_pc (float): The sector edge length, parsecs.
        galaxy_shape (dict or None): `apiclient.get_galaxy_shape`'s own
            return shape, or `None`.

    Returns:
        float: The outer radius, parsecs.
    """
    if galaxy_shape and galaxy_shape.get("outer_ring_index") is not None:
        return (galaxy_shape["outer_ring_index"] + 1) * edge_pc
    return GALAXY_RADIUS_PC


FETCH_RADIUS_FACTOR = 1.6
"""float: The map fetches content out to this multiple of the camera's
orbit radius, so what's just off-screen is already there when the camera
turns."""

MAX_TILES_PER_REQUEST = 128
"""int: Mirrors `queryDb.MAX_TILES_PER_REQUEST`."""


def _tile_level_for_view_radius(radius_pc):
    """Same as `galaxyViewport.tile_level_for_view_radius` (duplicated
    so this module keeps working under its ImportError fallback)."""
    if not radius_pc > 0:  # also NaN
        return TILE_MAX_LEVEL
    if math.isinf(radius_pc):
        return 0
    ratio = TILE_ROOT_EDGE_PC / radius_pc
    if math.isinf(ratio):
        return TILE_MAX_LEVEL
    return max(0, min(TILE_MAX_LEVEL, math.floor(math.log2(ratio))))


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


def initial_tile_request(orbit_radius_pc, center_pc=(0.0, 0.0, 0.0)):
    """
    The tiles the map's first frame needs, computed exactly the way
    `static/galaxymap3d.js`'s own `neededTiles` does for every later
    camera position, so the browser's first live fetch finds the first
    frame's tiles already cached.

    Args:
        orbit_radius_pc (float): The starting camera orbit radius
            (`view_radius_bounds`' max).
        center_pc (tuple): The starting camera target.

    Returns:
        list: Tile keys, nearest first.
    """
    view_radius = orbit_radius_pc * FETCH_RADIUS_FACTOR
    level = _tile_level_for_view_radius(view_radius)
    return _tiles_intersecting_sphere(level, center_pc, view_radius)


DENSITY_SHAPE_FIELDS = (
    "disk_scale_length_pc", "disk_scale_height_pc",
    "bulge_scale_radius_pc", "bulge_amplitude",
    "arm_count", "pitch_angle_rad", "arm_amplitude",
    "spiral_reference_radius_pc", "spiral_reference_angle_rad",
    "k_norm",
)
"""The `GalaxyShape` fields `static/galaxyprisms.js`'s `relativeDensity`
reads."""

MODEL_TERM_FIELDS = (
    "thick_disk_amplitude", "thick_disk_scale_length_pc", "thick_disk_scale_height_pc",
    "bulge_scale_y_pc", "bulge_scale_z_pc", "bar_cos", "bar_sin",
)
"""The `galaxyDensity.model_terms` it reads too (the thick disk and the
bar), which `queryDb.galaxy_density_shape` serves alongside the shape."""


def _density_shape(galaxy_shape):
    """Just the density model's own fields and terms from
    `apiclient.get_galaxy_shape`'s dict, or `None` without a shape (or with one missing a field), plus
    `sector_min_density`: the relative density a sector needs to expect
    one star (`1 / expected_system_count_at_density_1`), the same
    threshold the generator's skeleton uses, so the prisms outline exactly
    the galaxy's layers. Left out if the shape doesn't carry it."""
    if not galaxy_shape or any(galaxy_shape.get(field) is None for field in DENSITY_SHAPE_FIELDS + MODEL_TERM_FIELDS):
        return None
    shape = {field: galaxy_shape[field] for field in DENSITY_SHAPE_FIELDS + MODEL_TERM_FIELDS}
    expected = galaxy_shape.get("expected_system_count_at_density_1")
    if expected:
        shape["sector_min_density"] = 1.0 / expected
    return shape


def _escape(text):
    """Text as HTML, for the course hint's names and URLs (a sector or
    system name is generated, but never trusted as markup)."""
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _json_script(data):
    """Same `<script type="application/json">`-safe escaping
    `planetgen/web/maps/starmap.py`'s own `_json_script` uses -- see that function's
    docstring for why (a database name/label can contain `</script>`)."""
    return (
        json.dumps(data)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
    )


GALAXY_MAP_HELP = """<h3>Moving</h3>
  <ul>
    <li>Drag to turn the view, right-drag (or Shift-drag) to move it, scroll or pinch to zoom.</li>
    <li>Hover over the disk to see its arcs (each about 40&deg; of bearing by a third of the radius, top to
    bottom of the disk) and click one to zoom into it. Then pick a slab (a layer of the arc) on the map or with
    the buttons beside it, then a block of that slab, and so on down to single sectors.</li>
    <li>Back and Forward retrace your steps; the round Steps button lists them all, with Jump to newest. Up (or
    Esc) goes one step out, Whole galaxy (or Home) starts over, and Re-center puts the camera back to the
    current step's own view.</li>
    <li>Arrow keys and Enter pick too. The keys 1 to 9 open a bookmark while the map has focus.</li>
  </ul>
  <h3>Reading the map</h3>
  <ul>
    <li>Blocks are colored by predicted density (brighter is denser). Unfilled space is see-through; a block
    with generated sectors is more solid the more of them are generated.</li>
    <li>Glowing points are stars, sized by the star, colored by its temperature and brighter the more luminous:
    the brightest (1,000 L&#9737; and up on a new galaxy) everywhere, placed before their sectors are generated,
    and generated systems' stars fainter and fainter as you zoom in. Click one (inside a block) to see it.</li>
  </ul>
  <h3>Bookmarks</h3>
  <ul>
    <li>&#9734; on the breadcrumb bookmarks the view or the selected sector; Bookmarks opens one.</li>
  </ul>"""
"""The Galaxy Map's help dialog (UX.50): the gestures, keys and legend that used to sit under the map."""

SECTOR_MAP_HELP = """<h3>Moving</h3>
  <ul>
    <li>Drag to turn the view, right-drag (or Shift-drag) to move it, scroll or pinch to zoom, and Re-center (in
    Menu) to come back.</li>
    <li>Click a neighboring sector to open that sector's page.</li>
  </ul>
  <h3>Reading the map</h3>
  <ul>
    <li>A point of light is a star: its halo size shows brightness and its color temperature.</li>
    <li>Bright points are quasars, neutron stars and accreting black holes; small spheres are quiet black holes
    and interstellar comets.</li>
    <li>Translucent clouds are nebulae, asteroid fields and supernova remnants; faint clouds reach in from a
    neighboring sector.</li>
    <li>Faint points are rogue planets; Mark rogue planets rings them.</li>
  </ul>"""
"""The sector page's map (MAP.68): the help dialog (UX.50)."""

ADDRESS_BLOCK = """<form class="galaxy-address" id="galaxymap3d-address" role="search" hidden>
  <label for="galaxymap3d-address-input">Go to</label>
  <input type="text" id="galaxymap3d-address-input" name="address" autocomplete="off" spellcheck="false"
         placeholder="Designation, ring/layer/slot, x, y, z pc, or a name">
  <button type="submit" class="starmap-btn">Go</button>
</form>
<div class="galaxy-address-matches" id="galaxymap3d-matches" hidden></div>
"""
"""The Go-to form (the stage view's address bar)."""

CRUMBS_BLOCK = """<nav class="crumbs galaxy-crumbs" id="galaxymap3d-crumbs" aria-label="Map position" hidden></nav>
"""
"""The breadcrumb the stage view fills."""

SLABS_BLOCK = """<div class="galaxy-slabs" id="galaxymap3d-slabs" role="group" aria-labelledby="galaxymap3d-slabs-heading"></div>
"""
"""The slab buttons beside the map."""

ZOOM_BUTTONS = """  <button type="button" class="starmap-btn" data-action="zoom-out" aria-label="Zoom out">&minus;</button>
  <button type="button" class="starmap-btn" data-action="zoom-in" aria-label="Zoom in">+</button>
"""
"""A sector page's map: the zoom buttons."""

HISTORY_BUTTONS = """  <button type="button" class="starmap-btn" data-action="back" data-icon="back" disabled>Back</button>
  <details class="galaxy-steps" id="galaxymap3d-steps">
    <summary class="starmap-btn galaxy-steps-button" data-icon="steps" aria-label="Steps to here"
             title="Steps to here: go back to any of them"><span aria-hidden="true">&#9679;</span><span class="galaxy-steps-label">Steps</span></summary>
    <div class="galaxy-steps-panel" data-steps-panel></div>
  </details>
  <button type="button" class="starmap-btn" data-action="forward" data-icon="forward" disabled>Forward</button>
  <button type="button" class="starmap-btn" data-action="up" data-icon="up" disabled
          title="One step back out (Esc)">Up</button>
  <button type="button" class="starmap-btn" data-action="reset" data-icon="reset"
          title="Back to the whole galaxy (Home)">Whole galaxy</button>
"""
"""Back, the steps menu (with Jump to newest), Forward, Up and Whole galaxy."""

def render_galaxy_map3d_panel(db_name, galaxy_shape, edge_pc, initial_view, fetch_path="/galaxy/tiles",
                              sector_url=None, generate=None, phenomenon_url=None,
                              system_url=None, stage_path="/galaxy/stage",
                              locate_path="/galaxy/locate", course=None,
                              territory_path="/galaxy/territories", pick=None, nav_url=None, pinned=None,
                              nebula_shape_path="/galaxy/nebula/{id}/shape"):
    """
    Builds the "Galaxy Map (3D)" panel: a `<canvas>` `static/
    galaxymap3d.js` renders an interactive WebGL scene into (always
    seen from above and driven by the drill-down's picks: quarter,
    layer, arc, ..., sector -- no free camera, MAP.17),
    plus a `<script type="application/json">` block carrying the
    zoom-range numbers (`view_radius_bounds`), the tile settings the
    client needs to pick tiles the same way `initial_tile_request` does,
    and `initial_view`'s own payload for the first frame -- everything
    after that first frame comes from the client's own live `fetch()`
    calls to `fetch_path` (the Flask `/galaxy/tiles`).

    Args:
        db_name (str): The database the page shows (from config, never
                       the URL). Only names the browser's own
                       `localStorage` tile cache (`storageKey`) and
                       bookmark list (the Bookmarks menu's
                       `data-bookmark-db`, `static/bookmarks.js`), so each
                       database's are kept apart; it is never sent back
                       to the server.
        galaxy_shape (dict or None): `apiclient.get_galaxy_shape`'s own
            return shape, or `None` if `planetgen plan` has never been
            run -- when `None`, the panel shows a hint that density
            shading isn't real yet (see
            `queryDb.galaxy_tiles`'s own `has_shape` field, which
            `initial_view` already carries through).
        edge_pc (float): The sector edge length, parsecs
                         (`initial_view["edge_pc"]`, passed separately
                         since `view_radius_bounds` needs it before
                         `initial_view` itself is fetched).
        initial_view (dict): `tilecache.fetch_tiles`' own return shape
            (`stamp`/`tiles`/`density`/`edge_pc`/`has_shape`), fetched by
            the `/galaxy` view for `initial_tile_request`'s tiles -- the
            zoomed-all-the-way-out starting view. The client seeds its
            tile cache with it.
        fetch_path (str): The tile endpoint's URL (`/galaxy/tiles`).
        sector_url (str, optional): A sector page URL with `{id}` where
            the id goes, for the info panel's "View sector" link (a real
            `<a href>`); without it the panel shows no link.
        generate (dict or None): For a logged-in admin only:
            `{"url", "csrfField", "csrfToken"}`
            (`web.helpers.generate_target`), where an unfilled sector's
            Generate buttons post. `None` (every visitor) shows no buttons.
        phenomenon_url (str, optional): A phenomenon page URL with
            `{type}` and `{id}` where they go, for a cloud's "View
            phenomenon" link; without it the panel shows no link.
        system_url (str, optional): A star system page URL with `{id}`
            where the id goes, for a filled bright star's "View system"
            link; without it the panel shows no link.
        stage_path (str): The drill-down's stage endpoint
            (`/galaxy/stage`), which the map opens on
            (`static/galaxystageview.js`).
        locate_path (str): The address bar's name lookup
            (`/galaxy/locate`).
        territory_path (str or None): The territory overlay's endpoint
            (`/galaxy/territories`), fetched when the Territories button
            is pressed. `None` (no polities generated yet) leaves the
            button and its legend out.
        nebula_shape_path (str): The endpoint of a nebula's mesh (`/galaxy/
            nebula/{id}/shape`, `{id}` where the nebula's id goes), which the
            map draws a nebula from once it is big enough on screen
            (`static/nebulamesh.js`, MAP.103).
        pick (dict or None): The NAV page's pick mode
            (`web/galaxy_views._pick_from_args`): `pick` ("from" or
            "to"), `banner` and `cancel` (the NAV page with the other
            endpoint kept), `other` (that endpoint) and `keep_name`,
            `keep_value` and `nav_url` (for the Bookmarks menu, which
            keeps the pick, NAV.40). It shows a banner, and only a choice
            holding something generated can be taken. The page's script carries the pick
            on (`static/navpick.js`): it adds the pick to the map's own
            URLs and to a sector click, and moves it on when "Start
            Here" or "End Here" is pressed.
        nav_url (str, optional): The NAV page's URL (`/nav`), for a
            the "Start Here" and "End Here" buttons (NAV.29), which with
            one end chosen open it; without it the panel shows none.
        pinned (dict or None): A sector page's map (MAP.68): the one
            sector it is locked to, `{"ring", "layer", "slot", "center_pc"}`
            (`center_pc`: its centre in the galaxy frame). The panel then
            opens that sector in place, with no steps, breadcrumb, address
            or history of its own (the page has its own title, bookmark
            and banners), and a click on another sector opens that
            sector's page.
        course (dict or None): A NAV course to draw over the map
            (`web/nav_page.galaxy_course`): `scope`, `points` (galaxy-
            frame parsecs), `sector` and `navUrl`. `None` draws none.

    Returns:
        str: A complete `<section class="panel">` block.
    """
    min_radius, max_radius = view_radius_bounds(edge_pc, galaxy_shape)
    edge_ly = pc_to_ly(edge_pc)

    scene_data = {
        "storageKey": db_name,
        "fetchPath": fetch_path,
        "stagePath": stage_path,
        "locatePath": locate_path,
        "territoryPath": territory_path,
        "nebulaShapePath": nebula_shape_path,
        "course": course,
        "pick": pick["pick"] if pick else None,
        # The NAV pick the page opened with (static/navpick.js carries on
        # from it): the other end already chosen and where Cancel goes.
        "pickOther": pick["other"] if pick else None,
        "pickCancel": pick["cancel"] if pick else None,
        "navUrl": nav_url,
        "sectorUrl": sector_url,
        "generate": generate,
        "phenomenonUrl": phenomenon_url,
        "systemUrl": system_url,
        "tileRootEdgePc": TILE_ROOT_EDGE_PC,
        "tileMaxLevel": TILE_MAX_LEVEL,
        "fetchRadiusFactor": FETCH_RADIUS_FACTOR,
        "maxTilesPerRequest": MAX_TILES_PER_REQUEST,
        "hasShape": bool(initial_view.get("has_shape")),
        # The galaxy's density model parameters -- static/galaxyprisms.js
        # evaluates planetgen.galaxy.density.relative_density from
        # these itself to shade the density prisms.
        "densityShape": _density_shape(galaxy_shape),
        "blockMinPx": BLOCK_MIN_PX,
        "blockBudget": BLOCK_BUDGET,
        "edgePc": edge_pc,
        "edgeLy": edge_ly,
        "fovDeg": CAMERA_FOV_DEG,
        "galaxyRadiusPc": galaxy_extent_pc(edge_pc, galaxy_shape),
        "galaxyEdgePc": galaxy_edge_pc(edge_pc, galaxy_shape),
        "minViewRadiusPc": min_radius,
        "maxViewRadiusPc": max_radius,
        "clickZoomFactorMin": CLICK_ZOOM_FACTOR_MIN,
        "clickZoomFactorMax": CLICK_ZOOM_FACTOR_MAX,
        "initialCenter": list(pinned["center_pc"]) if pinned else [0.0, 0.0, 0.0],
        "initialRadiusPc": min(max_radius, max(min_radius, 6 * edge_pc)) if pinned else max_radius,
        "pinned": {k: pinned[k] for k in ("ring", "layer", "slot")} if pinned else None,
        # The first frame's tiles, without the tile cache's own count of how
        # many it had cached: the page is cached and compared as text, and
        # that count differs between two loads of the same page.
        "initial": {key: value for key, value in initial_view.items() if key != "cached"},
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

    pick_banner = ""
    bookmark_pick = ""
    if pick:
        # The Bookmarks menu keeps the pick (NAV.40, static/bookmarks.js).
        bookmark_pick = "".join(
            f' data-{name}="{_escape(pick[key])}"'
            for name, key in (("pick", "pick"), ("keep-name", "keep_name"), ("keep-value", "keep_value"),
                              ("nav-url", "nav_url"))
        )
        pick_banner = (
            '<p class="pick-banner" role="status"><strong>' + _escape(pick["banner"]) + ":</strong> drill down"
            " to a generated sector and click it, then pick a system or phenomenon there &middot; "
            '<a href="' + _escape(pick["cancel"]) + '">Cancel</a></p>\n'
        )
    course_hint = ""
    if course and course.get("points"):
        ends = [course["points"][0]["name"], course["points"][-1]["name"]]
        stops = len(course["points"]) - 2
        course_hint = (
            '<p class="hint">Showing the course from <strong>' + _escape(ends[0]) + "</strong> to <strong>"
            + _escape(ends[1]) + "</strong>"
            + (f" ({stops} stop{'s' if stops != 1 else ''} on the way)" if stops else "")
            + ' &middot; <a href="' + _escape(course["navUrl"]) + '">back to NAV</a></p>'
        )
    elif course and course.get("sector"):
        course_hint = (
            '<p class="hint">That whole course sits inside <strong>' + _escape(course["sector"]["name"])
            + '</strong>, so the map shows that sector &middot; <a href="' + _escape(course["navUrl"])
            + '">back to NAV</a></p>'
        )

    shape_hint = (
        ""
        if galaxy_shape
        else (
            '<p class="hint">The galaxy\'s density skeleton hasn\'t been built yet '
            "(<code>planetgen plan</code>) -- only generated sectors are shown, as blocks, "
            "and no density shading is shown.</p>"
        )
    )

    rogue_button = (
        '    <button type="button" class="starmap-btn starmap-toggle" data-action="toggle-rogue-markers" data-icon="rogue-markers"'
        ' aria-pressed="false"\n'
        '            title="Ring each rogue planet, so it is easy to find among the stars">Mark rogue planets</button>\n'
        if pinned else ""
    )
    # MAP.131: what the blocks' fill shows, and the legend of the choice.
    color_block = legend_block = ""
    if not pinned:
        color_block = (
            '    <label class="galaxy-color-by">Color by <select id="galaxymap3d-color-by" data-color-by\n'
            '            title="What the blocks and sectors are colored by">\n'
            '      <option value="default">Age, density and luminosity</option>\n'
            '      <option value="density">Density</option>\n'
            '      <option value="age">Mean age</option>\n'
            '      <option value="luminosity">Luminosity</option>\n'
            '      <option value="stars">Star count</option>\n'
            '    </select></label>\n')
        legend_block = (
            '<div class="galaxy-legend" id="galaxymap3d-legend" role="img" hidden>'
            '<span class="galaxy-legend-title"></span>'
            '<span class="galaxy-legend-bar"></span>'
            '<span class="galaxy-legend-low"></span><span class="galaxy-legend-high"></span></div>\n')
    territory_button = territory_box = ""
    if territory_path:
        territory_button = (
            '    <button type="button" class="starmap-btn" data-action="territories" data-icon="territories"'
            ' aria-pressed="false"\n'
            '            title="Show which polity holds what: each one\'s reach and the systems it owns">'
            'Territories</button>\n')
        territory_box = '<div class="galaxy-territories" id="galaxymap3d-territories" hidden></div>\n'

    title = "Sector Map" if pinned else "Galaxy Map"
    info_hint = (
        "Click a star, cloud or body for details." if pinned
        else "Click an arc of the galaxy (a piece of the disk, top to bottom) to look at it more closely."
    )
    help_html = SECTOR_MAP_HELP if pinned else GALAXY_MAP_HELP
    if pinned:
        # A sector page's map has no steps, breadcrumb, address or slab
        # buttons; the page carries the pick banner and the bookmark.
        pick_banner = bookmark_pick = ""
        zoom_buttons = ZOOM_BUTTONS
        address_block = crumbs_block = slabs_block = history_buttons = bookmarks_block = ""
    else:
        zoom_buttons = ""
        address_block, crumbs_block, slabs_block, history_buttons = ADDRESS_BLOCK, CRUMBS_BLOCK, SLABS_BLOCK, HISTORY_BUTTONS
        bookmarks_block = f"""  <details class="bookmarks-menu" data-bookmarks-menu data-bookmarks-keys="map" data-bookmark-db="{_escape(db_name)}"{bookmark_pick}>
    <summary class="starmap-btn" data-icon="bookmarks"
             title="Places saved with the breadcrumb's &#9734; (1 to 9 open the first nine while the map has focus)">Bookmarks</summary>
    <div class="bookmarks-panel" data-bookmarks-panel></div>
  </details>
"""

    return f"""
<section class="panel galaxymap3d-panel" id="map">
<div class="panel-header">
  <h2 class="sr-only">{title}</h2>
</div>
{pick_banner}{shape_hint}{course_hint}
{address_block}<p class="hint galaxy-stage-notice" id="galaxymap3d-notice" role="status" hidden></p>
{crumbs_block}<div class="starmap-layout">
<div class="galaxy-map-row">
<div class="galaxy-map-main">
<div class="starmap-viewport">
<canvas id="galaxymap3d-canvas" class="starmap-canvas" tabindex="0" role="application"
     aria-label="Interactive Galaxy Map. Arrow keys move among the parts you can pick and Enter takes
     one; slabs are also picked with the buttons beside the map; Escape or Backspace goes one step back out
     and Home returns to the whole galaxy. Dragging, or Shift and the arrow keys, turns the view and the wheel zooms."></canvas>
<div class="starmap-scale" id="galaxymap3d-scale" aria-live="polite"></div>
{legend_block}
<div class="map-tooltip" id="galaxymap3d-tooltip" hidden></div>
</div>
{slabs_block}</div>
<aside class="starmap-info" id="galaxymap3d-info">
<p class="hint">{info_hint}</p>
</aside>
</div>
<div class="starmap-side">
<div class="starmap-controls" id="galaxymap3d-controls">
{zoom_buttons}{history_buttons}{bookmarks_block}  <details class="galaxy-menu" id="galaxymap3d-menu">
    <summary class="starmap-btn" data-icon="menu" title="More map controls">Menu</summary>
    <div class="galaxy-menu-panel" role="group" aria-label="More map controls">
    <div class="galaxy-kinds" id="galaxymap3d-kinds" role="group" aria-label="Show on the map" hidden></div>
    <button type="button" class="starmap-btn" data-action="reset-view" data-icon="reset-view"
            title="Back to this step's own view after turning, moving or zooming it">Re-center</button>
    <button type="button" class="starmap-btn" data-action="charted-only" data-icon="charted-only" aria-pressed="false"
            title="Dim the stars and blocks outside charted sectors and outline the charted ones">Charted only</button>
{color_block}{rogue_button}{territory_button}    <button type="button" class="starmap-btn" data-action="map-help" data-icon="help"
            title="How to move around the map and read it">Map help</button>
    </div>
  </details>
</div>
{territory_box}</div>
</div>
<sl-dialog id="galaxymap3d-help" class="map-help-dialog" label="{title} help">
  {help_html}
  <sl-button slot="footer" data-dialog-close>Close</sl-button>
</sl-dialog>
<script type="application/json" id="galaxymap3d-data">{_json_script(scene_data)}</script>
</section>
"""
