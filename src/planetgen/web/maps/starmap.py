# planetgen/web/maps/starmap.py

"""
The Sector Map's scene data: every placed star system in a sector (two,
overlapping, for a binary), every nearby standalone phenomenon
(nebula/asteroid field/supernova remnant/black hole/neutron star/rogue
planet/interstellar comet/quasar), and a small clickable indicator toward
each immediately surrounding sector (`_neighbor_indicator_data`, see
`queryDb.sector_neighbors`), as the JSON `map_scene_data` returns. The
sector page answers it at `/sector/<id>/scene`, and the Galaxy Map draws it
(`static/galaxysector.js`, `static/sectorscene.js`, via three.js, vendored
at `static/vendor/three.module.min.js`, see that directory's
`THIRD_PARTY_NOTICES.txt`) in place, the sector page's map being the Galaxy
Map locked to the sector (MAP.68). This module owns every astrophysical and
layout decision: the client only positions, colors and labels what it is
handed.

Every star's (x, y, z) -- rotated once, server-side, from its own
sector-local axes into the galaxy frame when the sector has a galaxy
placement (`_rotate_to_galaxy_frame`, so it agrees with the cell outline/
compass arrow below, both already galaxy-frame quantities) -- is expressed
in fixed pixel-scale world units (`_SCENE_HALF_PX` per half the sector's
own edge); the client scales them by `halfEdgePc` to put the sector where
it sits in the galaxy.

`_outline_data` computes the scene's own bounding shape -- one of two,
depending on whether the sector has a galaxy placement: a sector with a
grid address gets its real cylindrical cell (see
`galaxyGeometry.sector_cell_vertices_pc` -- one ring wide, one layer tall,
one slot's angle across, its curved edges sampled as arcs); one without a
placement falls back to a plain axis-aligned cube of edge `edge_mpc`.

Links: `map_scene_data` takes a `link_url(name, **params)` callable (the
Flask pages pass `web.helpers.page_url`) and calls it as
`link_url("system", system_id=...)`, `link_url("phenomenon",
phenomenon_type=..., phenomenon_id=...)` and `link_url("sector",
sector_id=...)`; every scene entry carries the result as `href`.
"""

import colorsys
import math

from planetgen.web.lib.fmt import format_distance_ly, format_number

try:
    from planetgen.physics.constants import SPECTRAL_CLASS_COLORS, TEMP_RANGES, SOLAR_LUMINOSITY, SOLAR_RADIUS_M
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # duplicated literally rather than left unimportable, since these
    # drive the map's color math directly.
    SPECTRAL_CLASS_COLORS = {'O': 'Blue', 'B': 'Blue-White', 'A': 'White', 'F': 'Yellow-White', 'G': 'Yellow', 'K': 'Orange', 'M': 'Red'}
    TEMP_RANGES = {
        'O': (30000, 60000), 'B': (10000, 30000), 'A': (7500, 10000), 'F': (6000, 7500),
        'G': (5200, 6000), 'K': (3700, 5200), 'M': (2400, 3700),
    }
    SOLAR_LUMINOSITY = 3.82e26
    SOLAR_RADIUS_M = 6.957e8

try:
    from planetgen.galaxy.geometry import (
        layer_bounds_pc,
        ring_bounds_pc,
        sector_cell_vertices_pc,
        sector_orientation,
        sector_position_pc,
        slot_angle_bounds,
    )
    from planetgen.physics.units import ly_to_milliparsecs, milliparsecs_to_ly, mpc_to_pc, pc_to_mpc
except ImportError:
    # Same deployment gap as above -- without these, a sector with no
    # galaxy placement (or one this deployment can't reach the geometry
    # module for) just keeps the plain axis-aligned cube outline and no
    # scale bar; see `_cell_edges_px`/`_ly_per_px_at_zoom_1` below. Without
    # `sector_orientation` specifically, star dots fall back to being
    # plotted in their own unrotated local axes (see `_rotate_to_galaxy_frame`).
    # Without `ly_to_milliparsecs`, nebula/asteroid-field clouds can't be
    # scaled/positioned at all and are omitted.
    sector_position_pc = None
    sector_cell_vertices_pc = None
    sector_orientation = None
    ring_bounds_pc = None
    layer_bounds_pc = None
    slot_angle_bounds = None
    mpc_to_pc = None
    pc_to_mpc = None
    milliparsecs_to_ly = None
    ly_to_milliparsecs = None

# The 3D scene's world-unit scale -- every position/radius below is in
# these units (a sector's own half-edge maps to this many of them), handed
# to `sectorscene.js` as-is; a real perspective camera doesn't otherwise care
# what unit "1" means, unlike the old CSS version where this was also a
# literal pixel count.
_SCENE_HALF_PX = 160.0

# A plain axis-aligned fallback cube's own corners, at distance
# `_SCENE_HALF_PX` out along all three axes at once -- used by
# `_default_zoom` as the fallback shape's extent when there's no cell
# wireframe to measure instead.
_CUBE_CORNER_RADIUS_PX = _SCENE_HALF_PX * math.sqrt(3)

# Never start a sector further zoomed out than this, however large its
# cell/cube/plotted content gets (an extreme far-out placement could
# otherwise compute an unusably tiny initial view) -- `sectorscene.js`'s own
# zoom-out control/scroll remains available past this floor regardless.
_MIN_DEFAULT_ZOOM = 0.2


def _default_zoom(extent_radii_px):
    """
    Picks the zoom level the sector map should *open* at, so its real
    content -- the cell/cube outline, every plotted star/cloud -- actually
    fits in frame on first paint, instead of always starting at a flat
    `zoom = 1` regardless of how much bigger the real shape is.

    A galaxy-placed sector's cell outline is allowed to extend past
    `_SCENE_HALF_PX` (see `_cell_edges_px`'s docstring -- a curved cell's
    corners don't coincide with a cube's),
    and even the plain-cube fallback's own corners sit `_SCENE_HALF_PX *
    sqrt(3)` out -- both already past a "fits at zoom 1" frame.
    `sectorscene.js`'s camera distance at zoom 1 is calibrated so a sphere of
    radius `_SCENE_HALF_PX` exactly fills the frame (see its own
    `_referenceDistance`), so this fraction is what that client-side
    calibration is relative to.

    Args:
        extent_radii_px (list[float]): Distance from the scene's own
                                       center, for every point that matters
                                       (cell/cube vertices, star/cloud
                                       positions) -- the empty-scene case
                                       (no systems, no clouds, cube
                                       fallback) still always has at least
                                       `_CUBE_CORNER_RADIUS_PX` in this list.

    Returns:
        float: A zoom factor in `(0, 1]` -- `1.0` when everything already
              fits at the scene's native size, smaller the more the real
              content overflows it, floored at `_MIN_DEFAULT_ZOOM`.
    """
    max_radius = max(extent_radii_px, default=_SCENE_HALF_PX)
    if max_radius <= _SCENE_HALF_PX:
        return 1.0
    return max(_MIN_DEFAULT_ZOOM, _SCENE_HALF_PX / max_radius)

_SUN_RADIUS_KM = SOLAR_RADIUS_M / 1000.0
# A star's core radius in scene units: no longer drawn (the map draws a
# point of light a fixed number of pixels across, `_star_light`), but it
# still spaces a binary's pair apart and sizes the highlight ring.
_MIN_DOT_R = 1.5
_MAX_DOT_R = 6.0

# Secondary-star offset (binary systems), as a fraction of the primary's
# own dot radius -- "down and to the right", overlapping the primary
# rather than sitting fully clear of it. This is a fixed delta in the
# scene's own local coordinate space (not screen space), so the pair stays
# rigidly together as one unit regardless of camera angle instead of
# needing to be recomputed as the view turns.
_BINARY_OFFSET_FRACTION = 0.85

# The secondary's dot never exceeds this fraction of the primary's own
# radius, even when the secondary star is physically the larger of the
# two (a real possibility -- role is "which formed/is named first", not
# "which is bigger") -- the secondary is meant to read as a small partner
# badge on the main dot, not compete with it for visual weight.
_SECONDARY_MAX_RATIO = 0.65

# Anchor RGB for each of `SPECTRAL_CLASS_COLORS`' values -- hand-picked to
# read unambiguously as their name (a real O-star's blackbody tint is a
# much subtler blue than this) since the whole point is "a White Giant
# looks white, a Blue Giant looks blue" at a glance, not colorimetric
# accuracy.
_COLOR_NAME_RGB = {
    "Blue": (94, 140, 255),
    "Blue-White": (176, 202, 255),
    "White": (255, 255, 255),
    "Yellow-White": (255, 244, 214),
    "Yellow": (255, 224, 102),
    "Orange": (255, 154, 77),
    "Red": (255, 90, 77),
}
_NEUTRAL_RGB = (200, 200, 200)

# Spectral letter -> anchor RGB, derived from SPECTRAL_CLASS_COLORS so a
# change to that mapping (or TEMP_RANGES) stays in sync here automatically.
_SPECTRAL_LETTER_RGB = {
    letter: _COLOR_NAME_RGB.get(color_name, _NEUTRAL_RGB)
    for letter, color_name in SPECTRAL_CLASS_COLORS.items()
}

# Which pairs of 8-point vertex-list indices form a cube (or cell)'s 12
# edges -- indices are `4*a + 2*b + c` for whichever three binary choices
# the 8 corners vary over (see `sector_cell_vertices_pc`'s docstring for
# the cell's own r/z/theta bits; `_cube_corners_px` uses the same
# convention for its x/y/z sign bits below), so two vertices share an edge
# exactly when their indices differ in a single bit -- this is just a
# cube's own edge topology, reused for the cell's 8 corners regardless of
# how curved its faces really are, and for the plain fallback cube's own
# corners too.
_EDGE_PAIRS = (
    (0, 1), (0, 2), (0, 4),
    (1, 3), (1, 5),
    (2, 3), (2, 6),
    (3, 7),
    (4, 5), (4, 6),
    (5, 7),
    (6, 7),
)


def _cube_corners_px(half_edge_px):
    """
    The plain axis-aligned fallback cube's own 8 corners, in the same
    index convention `_EDGE_PAIRS` expects (`4*x_bit + 2*y_bit + z_bit`) --
    the wireframe-cube replacement for the old CSS version's 6 filled,
    bordered `.cube-face` divs (see this module's own docstring for why a
    wireframe reads just as clearly and lets both this and the cell case
    share one edges-list shape).
    """
    return [
        (
            half_edge_px if x_bit else -half_edge_px,
            half_edge_px if y_bit else -half_edge_px,
            half_edge_px if z_bit else -half_edge_px,
        )
        for x_bit in (0, 1)
        for y_bit in (0, 1)
        for z_bit in (0, 1)
    ]


def _rotate_to_galaxy_frame(center_pc, local_vec):
    """
    Re-expresses `local_vec` -- a star system's position relative to its
    sector's own center (`star_systems.position_x/y/z_mpc`, in whatever
    axes the sector's own local frame uses) -- along the galaxy frame's
    axes instead, via `galaxyGeometry.sector_orientation`'s fixed
    convention (local `+X` radially outward from the galactic axis, `+Y`
    along the ring, `+Z` galactic north), computed from the stored
    `center_pc` alone.

    Without this, a star dot's local (x, y, z) was plotted as if it were
    already expressed in the galaxy frame -- consistent with itself, but
    not with the cell outline or the "Galactic Center" compass arrow
    (`_cell_edges_px`/`_compass_data`), which read the sector's *actual*
    galaxy-frame placement directly and always were correct. Rotating the
    star dots into that same frame is what makes all three agree.

    Args:
        center_pc (tuple or None): `(center_x_pc, center_y_pc,
                                    center_z_pc)` -- `None` for a sector
                                    with no galaxy placement, in which
                                    case there's no galaxy frame to
                                    rotate into at all and `local_vec` is
                                    returned unchanged (the map's fallback
                                    plain-cube behavior, same as before
                                    this rotation existed).
        local_vec (tuple): `(x, y, z)`, in the sector's own local axes --
                           any consistent unit (this only rotates
                           direction, never rescales).

    Returns:
        tuple: `(x, y, z)`, re-expressed along the galaxy frame's axes,
              in `local_vec`'s original units.
    """
    if center_pc is None or any(c is None for c in center_pc) or sector_orientation is None:
        return local_vec
    if math.hypot(center_pc[0], center_pc[1]) < 1e-9:
        return local_vec

    axis_x, axis_y, axis_z = sector_orientation(center_pc)
    lx, ly, lz = local_vec
    return (
        lx * axis_x[0] + ly * axis_y[0] + lz * axis_z[0],
        lx * axis_x[1] + ly * axis_y[1] + lz * axis_z[1],
        lx * axis_x[2] + ly * axis_y[2] + lz * axis_z[2],
    )


def _cell_edges_px(address, edge_mpc, half_edge):
    """
    The sector's real cylindrical cell (`sector_cell_vertices_pc`) as 8
    `(x, y, z)` scene-space points, in the same convention
    `map_scene_data` uses for star dots (sector-center-relative parsecs
    -> milliparsecs -> normalized by `half_edge` -> scene units, y flipped
    since +y is "up" on screen) -- so the outline and the stars share one
    frame. Not clamped: a core-ring cell is a pie wedge noticeably bigger
    than the cube.

    Args:
        address (tuple or None): `(ring_index, layer_index,
            ring_slot_index)`; any `None` means no grid address.

    Returns:
        list[tuple] or None: 8 points, or `None` for a sector with no grid
            address, or when the geometry helpers aren't importable --
            callers fall back to the plain cube.
    """
    if (
        address is None
        or any(value is None for value in address)
        or not edge_mpc
        or sector_cell_vertices_pc is None
        or sector_position_pc is None
        or mpc_to_pc is None
        or pc_to_mpc is None
    ):
        return None

    edge_pc = mpc_to_pc(edge_mpc)
    center_pc = sector_position_pc(*address, edge_pc)
    vertices_pc = sector_cell_vertices_pc(*address, edge_pc)

    points = []
    for vx, vy, vz in vertices_pc:
        nx = pc_to_mpc(vx - center_pc[0]) / half_edge
        ny = pc_to_mpc(vy - center_pc[1]) / half_edge
        nz = pc_to_mpc(vz - center_pc[2]) / half_edge
        points.append((nx * _SCENE_HALF_PX, -ny * _SCENE_HALF_PX, nz * _SCENE_HALF_PX))
    return points


_ARC_STEP_RAD = math.pi / 90
"""Angle between samples along a cell's curved edges (2 degrees): a slot's
inner and outer faces are arcs of the ring's radius, so each arc edge is
drawn as a polyline of this many-degree steps rather than one chord."""


def _cell_arc_edges_px(address, edge_mpc, half_edge):
    """
    The sector's cell outline as 12 polylines in the same scene frame as
    `_cell_edges_px`: the 8 edges between corners that differ in radius or
    height are straight, and the 4 that differ in slot angle (corners `i`
    and `i ^ 1`) follow the ring's arc, sampled every `_ARC_STEP_RAD`.

    Returns:
        list[list[tuple]] or None: One point list per `_EDGE_PAIRS` edge
            (2 points for a straight edge, more for an arc), or `None`
            when `_cell_edges_px` would be.
    """
    corners = _cell_edges_px(address, edge_mpc, half_edge)
    if corners is None or ring_bounds_pc is None:
        return None

    ring_index, layer_index, slot_index = address
    edge_pc = mpc_to_pc(edge_mpc)
    center_pc = sector_position_pc(*address, edge_pc)
    r_bounds = ring_bounds_pc(ring_index, edge_pc)
    z_bounds = layer_bounds_pc(layer_index, edge_pc)
    t_start, t_end = slot_angle_bounds(ring_index, slot_index)
    steps = max(2, int(math.ceil(abs(t_end - t_start) / _ARC_STEP_RAD)))

    def to_px(r, z, theta):
        nx = pc_to_mpc(r * math.cos(theta) - center_pc[0]) / half_edge
        ny = pc_to_mpc(r * math.sin(theta) - center_pc[1]) / half_edge
        nz = pc_to_mpc(z - center_pc[2]) / half_edge
        return (nx * _SCENE_HALF_PX, -ny * _SCENE_HALF_PX, nz * _SCENE_HALF_PX)

    edges = []
    for i, j in _EDGE_PAIRS:
        if i ^ j == 1:
            r, z = r_bounds[(i >> 2) & 1], z_bounds[(i >> 1) & 1]
            edges.append([to_px(r, z, t_start + (t_end - t_start) * k / steps) for k in range(steps + 1)])
        else:
            edges.append([corners[i], corners[j]])
    return edges


def _outline_data(address, edge_mpc, half_edge):
    """
    Builds the scene's outline as `{"kind": "cell"|"cube", "edges": [...]}`
    -- 12 edges (see `_EDGE_PAIRS`), each a list of points to draw as one
    line. A sector with a grid address gets its real cylindrical cell
    (`_cell_arc_edges_px`), whose inner and outer faces' edges are sampled
    arcs; anything else gets the plain axis-aligned fallback cube's
    straight edges (`_cube_corners_px`). `sectorscene.js` draws it as a
    faint wireframe.

    Returns:
        tuple[dict, list[float]]: The outline data, and the extent (from
                                  the scene's own center) of every one of
                                  its vertices -- for `_default_zoom`'s
                                  empty-scene fallback.
    """
    edges = _cell_arc_edges_px(address, edge_mpc, half_edge)
    kind = "cell"
    if edges is None:
        vertices, kind = _cube_corners_px(_SCENE_HALF_PX), "cube"
        edges = [[vertices[i], vertices[j]] for i, j in _EDGE_PAIRS]

    points = [point for edge in edges for point in edge]
    extent_radii_px = [math.sqrt(vx * vx + vy * vy + vz * vz) for vx, vy, vz in points]
    edges = [[[round(c, 2) for c in point] for point in edge] for edge in edges]
    return {"kind": kind, "edges": edges}, extent_radii_px


_COMPASS_ARROW_REACH = 1.15
"""How far the compass arrow reaches past `_SCENE_HALF_PX` -- past 1.0 so
its tip clears a full-size cube/cell instead of ending right at (or
inside) its own boundary."""


def _compass_data(center_pc):
    """
    Points from the sector's own local origin toward the galactic center --
    the sector-map's own "north" arrow (labeled plain "N", the same
    convention a real map's compass rose uses), except the direction it
    points is computed exactly from this sector's own stored
    `sectors.center_x/y/z_pc` (the negative of the sector's own outward
    radial direction, `-normalize(center_pc)`) rather than fixed to a
    constant screen direction.

    This arrow and `_outline_data`'s cell outline are both computed
    directly from galaxy-frame quantities (`sectors.center_x/y/z_pc`,
    `sector_cell_vertices_pc`), so they were always correct on their own
    terms. Star dots are rotated into that same frame at render time
    instead (`_rotate_to_galaxy_frame`), so this arrow, the cell outline,
    and the star dots it surrounds all agree on one frame.

    Args:
        center_pc (tuple or None): `(center_x_pc, center_y_pc,
                                    center_z_pc)`, or `None` if this
                                    sector has no galaxy placement.

    Returns:
        dict or None: `{"tip": [x, y, z], "label": "N"}`, or `None` if
                      `center_pc` is `None` or (within floating-point
                      tolerance) the galactic center itself, which has no
                      meaningful direction to point.
    """
    if center_pc is None or any(c is None for c in center_pc):
        return None

    cx, cy, cz = center_pc
    norm = math.sqrt(cx * cx + cy * cy + cz * cz)
    if norm < 1e-9:
        return None

    ux, uy, uz = -cx / norm, -cy / norm, -cz / norm
    reach = _SCENE_HALF_PX * _COMPASS_ARROW_REACH
    return {"tip": [ux * reach, -uy * reach, uz * reach], "label": "N"}


_NEIGHBOR_INDICATOR_REACH = 1.35
"""How far a neighboring-sector indicator sits past `_SCENE_HALF_PX`, as a
multiple of it -- past `_COMPASS_ARROW_REACH` so these never sit right on
top of the compass arrow's own tip, and (like that arrow) placed by
direction alone rather than at the neighbor's real, wildly varying
distance -- see `_neighbor_indicator_data`'s own docstring."""

_NEIGHBOR_INDICATOR_RADIUS_PX = 7.0
"""A neighboring-sector indicator's own drawn size -- fixed, unlike a star
dot's radius (`_star_dot_radius`), since there's no real "size" a sector
address has; about a bright star's halo so it reads as a marker of
similar visual weight."""


def _neighbor_indicator_data(link_url, neighbor):
    """
    Builds one neighboring-sector indicator's plain-dict scene entry --
    `sectorscene.js` draws it as a small clickable marker just outside this
    sector's own cell/cube, in the real direction (already galaxy-frame,
    same as the cell outline and compass arrow -- see this module's own
    docstring) of that neighbor's actual center, but placed at a fixed
    `_NEIGHBOR_INDICATOR_REACH` rather than that real (and highly
    variable -- a lateral neighbor sits about one sector-edge away, a
    radial one likewise, but neither exactly) distance, the same "point
    by direction, not by true distance" convention `_compass_data`'s own
    arrow already uses.

    Args:
        link_url (callable): Builds an existing neighbor's own `href`
                       (see the module docstring).
        neighbor (dict): One entry from `queryDb.sector_neighbors`.

    Returns:
        dict or None: The indicator's scene entry (`x`/`y`/`z`,
            `ringIndex`/`layerIndex`/`ringSlotIndex`, `designation`, `exists`, and --
            only when `exists` -- `name`/`href`), or
            `None` if this neighbor's own direction is (within
            floating-point tolerance) undefined -- not expected in
            practice (two distinct sector centers are never that close),
            but the same defensive-`None` convention `_compass_data` uses
            for its own zero-length case.
    """
    dx, dy, dz = neighbor["direction_pc"]
    norm = math.sqrt(dx * dx + dy * dy + dz * dz)
    if norm < 1e-9:
        return None

    ux, uy, uz = dx / norm, dy / norm, dz / norm
    reach = _SCENE_HALF_PX * _NEIGHBOR_INDICATOR_REACH
    data = {
        "x": ux * reach, "y": -uy * reach, "z": uz * reach,
        "r": _NEIGHBOR_INDICATOR_RADIUS_PX,
        "isNeighbor": True,
        "ringIndex": neighbor["ring_index"], "layerIndex": neighbor["layer_index"],
        "ringSlotIndex": neighbor["ring_slot_index"],
        "designation": neighbor["designation"],
        "exists": neighbor["exists"],
    }
    if neighbor["exists"]:
        data["name"] = neighbor["sector_name"]
        data["href"] = link_url("sector", sector_id=neighbor["sector_id"])
        data["sectorId"] = neighbor["sector_id"]
    bright = neighbor.get("bright_stars")
    if bright:
        data["brightStarCount"] = len(bright)
        data["brightStars"] = [
            f'{star["star_type"]}, {format_number(star["luminosity_sol"])} L\u2609' for star in bright[:_NEIGHBOR_BRIGHT_STARS_SHOWN]
        ]
    return data


_NEIGHBOR_BRIGHT_STARS_SHOWN = 5
"""int: How many of an unfilled neighbor's waiting bright stars its info
panel lists (brightest first); the rest are counted."""


def _ly_per_px_at_zoom_1(half_edge):
    """
    Light-years per world unit at zoom factor 1 -- an exact ratio derived
    straight from `half_edge` (half of `edge_mpc`, the sector's real,
    stored size) mapping to `_SCENE_HALF_PX` units, the same normalizing
    divisor every star dot and cell vertex on this map is already placed
    by. `sectorscene.js`'s scale-bar legend divides this by the live zoom
    factor and picks a round bar length from it, so the bar always
    reflects the sector's actual physical scale rather than an arbitrary
    fixed guess.

    Returns:
        float or None: ly per unit, or `None` if `milliparsecs_to_ly`
                       isn't importable in this deployment (the scale bar
                       is then omitted entirely rather than shown in raw
                       milliparsecs).
    """
    if milliparsecs_to_ly is None:
        return None
    mpc_per_unit = half_edge / _SCENE_HALF_PX
    return milliparsecs_to_ly(mpc_per_unit)


def _kelvin_to_hex(temp_k):
    """
    Approximates a blackbody color for `temp_k` as a `#rrggbb` string --
    the standard Tanner Helland fit (clamped to its 1000-40000 K valid
    range). Used only as a fallback for `star_color` when `star_type`
    doesn't start with a recognized spectral letter (malformed/unexpected
    data) -- every normal star gets its color from the named spectral
    color instead, not this.
    """
    temp = max(1000.0, min(40000.0, temp_k)) / 100.0

    if temp <= 66:
        red = 255.0
    else:
        red = 329.698727446 * ((temp - 60) ** -0.1332047592)

    if temp <= 66:
        green = 99.4708025861 * math.log(temp) - 161.1195681661
    else:
        green = 288.1221695283 * ((temp - 60) ** -0.0755148492)

    if temp >= 66:
        blue = 255.0
    elif temp <= 19:
        blue = 0.0
    else:
        blue = 138.5177312231 * math.log(temp - 10) - 305.0447927307

    def _clamp(value):
        return max(0, min(255, round(value)))

    return "#{:02x}{:02x}{:02x}".format(_clamp(red), _clamp(green), _clamp(blue))


def _rgb_hex(rgb_float):
    r, g, b = rgb_float
    return "#{:02x}{:02x}{:02x}".format(
        max(0, min(255, round(r * 255))),
        max(0, min(255, round(g * 255))),
        max(0, min(255, round(b * 255))),
    )


def star_color(star_type, temperature_k, luminosity_w):
    """
    Picks a dot's fill and stroke color for one star, factoring in all
    four of color/temperature/brightness/size this map is meant to
    convey (size is handled separately by `_star_dot_radius`):

    - Color: the named spectral color baked into `star_type`'s leading
      letter (O/B/A/F/G/K/M -> Blue/Blue-White/White/Yellow-White/
      Yellow/Orange/Red, see `SPECTRAL_CLASS_COLORS`) -- so "White
      Giant" and "White Dwarf" both render white, "Blue Giant" blue,
      etc., regardless of luminosity class, matching the generator's own
      convention (color is a function of spectral letter only, never of
      giant/dwarf/supergiant).
    - Temperature: nudges lightness within that color's band by where
      `temperature_k` actually falls between the spectral class's own
      `TEMP_RANGES` bounds (hotter end of the band = slightly brighter),
      instead of only ever using the same shade for an entire class.
    - Brightness: `luminosity_w` (relative to `SOLAR_LUMINOSITY`, log
      scaled) pushes lightness/saturation up for a genuinely luminous
      star (giants/supergiants read as more vivid/glowing) and down for
      a dim one (red dwarfs, subdwarfs read as duller/more muted).

    Returns:
        tuple[str, str]: `(fill_hex, stroke_hex)` -- the stroke is a
                         darkened version of the same hue, so every dot
                         has a visible edge regardless of the page theme
                         or how light its own fill is (a pure white dot
                         would otherwise vanish against a light-mode
                         panel), and so two overlapping binary dots read
                         as distinct circles rather than one blob.
    """
    letter = (star_type or "")[:1].upper()
    anchor = _SPECTRAL_LETTER_RGB.get(letter)
    if anchor is None:
        # Unrecognized/malformed star_type -- fall back to a pure
        # temperature blackbody guess rather than a hardcoded neutral.
        fill = _kelvin_to_hex(temperature_k or 5778.0)
        return fill, fill

    hue, lightness, saturation = colorsys.rgb_to_hls(anchor[0] / 255, anchor[1] / 255, anchor[2] / 255)

    luminosity_solar = (luminosity_w / SOLAR_LUMINOSITY) if luminosity_w else 1.0
    log_luminosity = math.log10(max(luminosity_solar, 1e-6))
    brightness = max(-1.0, min(1.0, log_luminosity / 5.0))

    temp_lo, temp_hi = TEMP_RANGES.get(letter, (3000.0, 10000.0))
    if temperature_k and temp_hi > temp_lo:
        temp_fraction = max(0.0, min(1.0, (temperature_k - temp_lo) / (temp_hi - temp_lo)))
    else:
        temp_fraction = 0.5

    # Proportional to remaining headroom rather than a flat add -- every
    # anchor here is already fully saturated (S=1), so a flat lightness
    # add on a very luminous star (brightness -> 1) washes any hue out to
    # near-white well before it gets there (an "orange" supergiant
    # reading as pale peach). Scaling by (1 - lightness)/lightness keeps
    # a color recognizably itself while still visibly brightening/
    # dimming, and leaves White's l=1 anchor untouched either way (no
    # headroom to move in, up or down).
    if brightness >= 0:
        lightness = lightness + (1.0 - lightness) * brightness * 0.35
    else:
        lightness = lightness + lightness * brightness * 0.5
    lightness = max(0.08, min(0.97, lightness + (temp_fraction - 0.5) * 0.08))
    # Only nudge saturation up/down for an anchor that already has some --
    # "White" is deliberately achromatic (saturation 0), and the 0.2 floor
    # below would otherwise inject an arbitrary hue tint into it (colorsys
    # still returns *some* hue for a fully desaturated color, even though
    # it's meaningless) the moment any brightness adjustment applied.
    if saturation > 0:
        saturation = max(0.2, min(1.0, saturation + brightness * 0.15))

    fill = _rgb_hex(colorsys.hls_to_rgb(hue, lightness, saturation))
    stroke = _rgb_hex(colorsys.hls_to_rgb(hue, max(0.05, lightness - 0.28), saturation))
    return fill, stroke


_SUN_DOT_R = 3.0
"""The Sun's own dot radius -- the pivot of `_star_dot_radius`'s log scale."""
_DOT_R_PER_DECADE_BELOW_SUN = 0.75
_DOT_R_PER_DECADE_ABOVE_SUN = 1.0
"""How much a dot grows per factor of ten in radius: gently below the Sun
(0.01 solar radii, a white dwarf, lands on `_MIN_DOT_R`), a little faster
above it so 10, 100 and 1000 solar radii (giants, bright giants,
supergiants) each read a little larger, 1000 reaching `_MAX_DOT_R`."""


def _star_dot_radius(radius_km):
    """Maps a star's physical radius to its core's radius in scene units on
    a log scale pivoting at the Sun (`_SUN_DOT_R`), so a white dwarf, a red
    dwarf, the Sun, a giant and a supergiant are still different sizes;
    clamped to the narrow `_MIN_DOT_R`..`_MAX_DOT_R` so even a supergiant
    is a small point (its brightness shows in its halo, `_star_light`)."""
    if not radius_km or radius_km <= 0:
        return _MIN_DOT_R
    decades = math.log10(radius_km / _SUN_RADIUS_KM)
    per_decade = _DOT_R_PER_DECADE_ABOVE_SUN if decades > 0 else _DOT_R_PER_DECADE_BELOW_SUN
    dot_r = _SUN_DOT_R + per_decade * decades
    return max(_MIN_DOT_R, min(_MAX_DOT_R, dot_r))


# MAP.15: every star, and every phenomenon that gives off light, is drawn
# as a point of light the way the Galaxy Map draws its bright stars
# (`static/galaxymap3d.js`, "Bright stars"): a tiny bright core in a soft
# halo, both a fixed number of screen pixels across (`sectorscene.js` grows
# them a little, up to `POINT_CLOSE_GROWTH`, as the camera closes in),
# never a ball sized in scene units. `_star_light` works out each star's
# from its radius, luminosity and temperature with the Galaxy Map's own
# ranges, its faint end lifted a little so a red dwarf still shows on the
# smaller Sector Map.
_LIGHT_LOG_LUMINOSITY = (-4.0, 6.0)
"""log10(L / L_sun) mapped onto 0..1 for the halo's size and strength."""
_LIGHT_LOG_RADIUS = (-1.0, 3.0)
"""log10(R / R_sun) mapped onto 0..1 for the core's size."""
_LIGHT_CORE_PX = (1.8, 5.0)
"""The core's diameter in pixels, from a red dwarf to a supergiant."""
_LIGHT_SIZE_PX = (13.0, 40.0)
"""The halo's full diameter in pixels, from the faintest star to the most
luminous (grown by the square of the luminosity share, as on the Galaxy
Map, so only a really bright star gets a big halo)."""
_LIGHT_GLOW = (0.55, 0.85)
"""The halo's strength at its middle, faint star to bright star."""
_LIGHT_BRIGHT = (0.9, 1.0)
"""The core's opacity and how far it whitens, faint star to bright star."""


# MAP.87: the dim end of the luminosity range drawn brighter, on the
# Sector Map here and on the Galaxy Map (`static/starlight.js`, this
# curve's JavaScript twin; tests/js/starlight.test.mjs checks they agree).
# Display only: nothing stored changes.
_BOOST_AT_DIM = 4.0
"""How many times as bright the dimmest stars are drawn."""
_BOOST_LOG_LUMINOSITY = (-4.0, 3.0)
"""log10(L / L_sun) where the boost is `_BOOST_AT_DIM` (1e-4 L_sun and
fainter) and where it has tapered to none (1000 L_sun and brighter)."""


_BOOST_TAPER_POWER = 1.5
"""The boost's exponent falls as `share ** 1.5` over the log range: flat
at the dim end, and slow enough that a brighter star never ends up drawn
fainter than a dimmer one on either map (the halo's own light grows only
about 5.6 times from 1e-4 to 1000 L_sun on the Sector Map, so a faster
taper, a smoothstep say, would dip in the middle)."""


def star_light_boost(luminosity_solar):
    """
    How many times as bright a star's point of light is drawn (MAP.87):
    `_BOOST_AT_DIM` at the dim end of `_BOOST_LOG_LUMINOSITY`, 1 from its
    bright end up, tapering in between (`_BOOST_TAPER_POWER`); a Sun gets
    about 2.2. 1 for a star with no luminosity on record.
    """
    if not luminosity_solar or luminosity_solar <= 0:
        return 1.0
    lo, hi = _BOOST_LOG_LUMINOSITY
    share = max(0.0, min(1.0, (math.log10(luminosity_solar) - lo) / (hi - lo)))
    return _BOOST_AT_DIM ** (1.0 - share ** _BOOST_TAPER_POWER)


def _boost_light(size_px, glow, boost):
    """A point of light drawn `boost` times as bright: its halo's strength
    and area each grow by sqrt(boost) (its width by boost ** 0.25); the
    core, sized and lit by the star itself, is left alone. Returns
    `(size_px, glow)`."""
    root = math.sqrt(boost)
    return size_px * math.sqrt(root), glow * root


def _log_share(value, value_range):
    """Where `value`'s log10 falls in `value_range`, clamped to 0..1."""
    lo, hi = value_range
    share = (math.log10(max(value, 1e-12)) - lo) / (hi - lo)
    return max(0.0, min(1.0, share))


def _lerp(value_range, share):
    return value_range[0] + (value_range[1] - value_range[0]) * share


def _star_light(luminosity_w, radius_km, temperature_k):
    """
    A star's point of light (MAP.15), as `sectorscene.js` draws it: a core
    sized by the star's radius, a halo whose width and strength grow with
    its luminosity, the core a little dimmer for a faint star, the faint
    end drawn brighter (`star_light_boost`, MAP.87), all in the star's
    blackbody color (`_kelvin_to_hex`, the same fit the Galaxy
    Map's `starColor` uses).

    Returns:
        dict: `color` (`#rrggbb`), `corePx` (the core's diameter),
            `sizePx` (the halo's full diameter, never less than the core
            plus a pixel's rim each side), `glow` (the halo's strength)
            and `bright` (the core's opacity and whitening), every size
            in screen pixels.
    """
    luminosity_solar = (luminosity_w / SOLAR_LUMINOSITY) if luminosity_w else 1.0
    if radius_km and radius_km > 0:
        radius_solar = radius_km / _SUN_RADIUS_KM
    else:
        # Without a stored radius, guess one from the luminosity.
        radius_solar = luminosity_solar ** 0.35
    lum_share = _log_share(luminosity_solar, _LIGHT_LOG_LUMINOSITY)
    core_px = _lerp(_LIGHT_CORE_PX, _log_share(radius_solar, _LIGHT_LOG_RADIUS))
    size_px, glow = _boost_light(
        _lerp(_LIGHT_SIZE_PX, lum_share * lum_share), _lerp(_LIGHT_GLOW, lum_share),
        star_light_boost(luminosity_solar),
    )
    return {
        "color": _kelvin_to_hex(temperature_k or 5778.0),
        "corePx": round(core_px, 2),
        "sizePx": round(max(size_px, 2 * core_px + 2), 2),
        "glow": round(glow, 3),
        "bright": round(_lerp(_LIGHT_BRIGHT, lum_share), 3),
    }


# The phenomena that give off light, each a fixed point of light (their
# rows carry no luminosity): a quasar outshines everything on the map, an
# accreting black hole is a hot orange point (its disc, not the hole,
# shines, so its core is not whitened), a neutron star a small blue-white
# one. A quiescent black hole and an interstellar comet give off no light
# of their own and keep their spheres; a rogue planet is a faint point of
# its own (`_ROGUE_LIGHT`); nebulae, supernova remnants and asteroid
# fields stay clouds.
_PHENOMENON_LIGHTS = {
    "quasar": {"color": "#dfe6ff", "corePx": 5.0, "sizePx": 48.0, "glow": 0.9, "bright": 1.0},
    "blackHoleAccreting": {
        "color": "#ff9d4d", "corePx": 3.0, "sizePx": 28.0, "glow": 0.75, "bright": 0.95, "whiten": 0.0,
    },
    "neutronStar": {"color": "#a9d4ff", "corePx": 2.2, "sizePx": 22.0, "glow": 0.7, "bright": 1.0},
}

# A rogue planet gives off no light of its own, so unmarked it is a dim,
# tiny point with no glow, barely noticeable (MAP.82): Boss, "Rogue planet
# detail is for them to be dim barely noticeable". The "Mark rogue
# planets" button (off by default, MAP.83) swaps in `_ROGUE_MARKED_LIGHT`:
# bigger, fully lit and glowing, ringed, and easy to pick (MAP.84; the
# pick reach is sectorscene.js's ROGUE_PICK_PX).
_ROGUE_LIGHT = {"color": "#a993f0", "corePx": 1.5, "sizePx": 3.5, "glow": 0.0, "bright": 0.2, "whiten": 0.0}
_ROGUE_MARKED_LIGHT = {
    "color": "#b59cff", "corePx": 5.0, "sizePx": 18.0, "glow": 0.6, "bright": 1.0, "whiten": 0.3,
}


def _star_data(link_url, system, star, x_px, y_px, z_px, max_r=None):
    """
    Builds one star's plain-dict scene entry -- `sectorscene.js` draws it as
    a point of light (`light`, see `_star_light`); `r` is its core radius
    in scene units, still what spaces a binary's pair apart and sizes its
    highlight ring's fallback. Labeled with the star's own name (a
    binary's stars are `<system> <word>` -- see `bodyNames.py`), falling
    back to the system's.
    """
    dot_r = _star_dot_radius(star["radius_km"])
    if max_r is not None:
        dot_r = min(dot_r, max_r)
    fill, stroke = star_color(star["star_type"], star["temperature_k"], star["luminosity_w"])
    return {
        "x": x_px, "y": y_px, "z": z_px, "r": dot_r,
        "light": _star_light(star["luminosity_w"], star["radius_km"], star["temperature_k"]),
        "fill": fill, "stroke": stroke,
        "name": star.get("name") or system["name"],
        "starType": star["star_type"],
        # The star's luminosity in suns, for the map's luminosity floor (MAP.123).
        "luminositySol": (star["luminosity_w"] or 0) / SOLAR_LUMINOSITY,
        "temp": star["temp_display"],
        "quadrant": system["quadrant"],
        "location": system["location"],
        # `sectorscene.js`'s info panel links here with a plain `<a href>`.
        "href": link_url("system", system_id=system["id"]),
        # The system's NAV endpoint (`nav_page.endpoint`), what its ☆
        # Bookmark saves (`static/mappick.js`).
        "endpoint": f'system:{int(system["id"])}',
    }


# A cloud's own radius (`radius_ly`) is frequently far larger than a
# single sector -- a big emission nebula can dwarf the whole scene -- so
# unlike a star dot's radius (always tiny next to the scene), this needs
# its own generous cap: large enough to visibly engulf/overflow the frame
# without an unbounded sphere size for a pathological radius_ly value.
_MAX_CLOUD_RADIUS_PX = 6 * (2 * _SCENE_HALF_PX)

# Fill color (and base opacity, baked into the alpha channel below) per
# `nebulae.nebula_type` -- not spectral-accurate the way `star_color` is
# (a nebula's visible color really does vary this much by type: emission
# nebulae genuinely glow reddish-pink from ionized hydrogen's H-alpha
# line, reflection nebulae blue from scattered starlight, planetary
# nebulae teal/cyan from doubly-ionized oxygen, and a dark nebula is by
# definition an opaque silhouette, not a glow at all -- hence its own
# near-black, higher-opacity fill instead of a lighter translucent one).
NEBULA_SURROUNDINGS_PATH = "/galaxy/nebula/{id}/surroundings"
"""str: The site's endpoint for the stars round a nebula
(`web/galaxy_views.galaxy_nebula_surroundings`), for its page's 3D view."""

NEBULA_SHAPE_PATH = "/galaxy/nebula/{id}/shape"
"""str: The site's nebula mesh endpoint (`web/galaxy_views.galaxy_nebula_shape`),
`{id}` where the nebula's id goes."""

_NEBULA_TYPE_COLORS = {
    "diffuse": "#e3a6c8",
    "emission": "#ff6f91",
    "reflection": "#6fa8ff",
    "planetary": "#5be8c9",
    "dark": "#1c1c24",
}
_NEBULA_TYPE_ALPHA = {
    "diffuse": (0x78, 0x24),
    "emission": (0xB0, 0x40),
    "reflection": (0xA0, 0x38),
    "planetary": (0xA8, 0x3c),
    "dark": (0xE8, 0x90),
}
_DEFAULT_NEBULA_COLOR = "#c9a8e0"
_DEFAULT_NEBULA_ALPHA = (0x90, 0x38)

# Every other phenomenon kind's own look is a fixed recipe (never a
# function of per-instance data beyond which kind/descriptor it is), so
# unlike a nebula's per-instance core/edge color, `sectorscene.js` owns those
# recipes directly (see its own `CLOUD_KIND_RECIPES`) -- this module only
# ever needs to say *which* fixed kind a given phenomenon is.


def _phenomenon_cloud_radius_px(radius_ly, half_edge):
    if ly_to_milliparsecs is None or not half_edge:
        return 0.0
    radius_mpc = ly_to_milliparsecs(radius_ly)
    radius_px = (radius_mpc / half_edge) * _SCENE_HALF_PX
    return max(4.0, min(_MAX_CLOUD_RADIUS_PX, radius_px))


def nebula_look(descriptor):
    """`(color "#rrggbb", core opacity, edge opacity)` -- opacities 0 to 1 --
    of a nebula type, the look every map gives it."""
    color = _NEBULA_TYPE_COLORS.get(descriptor, _DEFAULT_NEBULA_COLOR)
    core_alpha, edge_alpha = _NEBULA_TYPE_ALPHA.get(descriptor, _DEFAULT_NEBULA_ALPHA)
    return color, core_alpha / 255, edge_alpha / 255


def _cloud_data(link_url, phenomenon, x_px, y_px, z_px, radius_px):
    """
    Builds one phenomenon's plain-dict scene entry. A nebula gets its own
    per-instance `coreColor`/`edgeColor` (two gradient stops, `#rrggbbaa`)
    computed from `_NEBULA_TYPE_COLORS`/`_NEBULA_TYPE_ALPHA`; every other
    kind carries no color data at all -- `sectorscene.js`'s own fixed
    recipes (mirroring this module's old `_ASTEROID_FIELD_BACKGROUND`/
    `_BLACK_HOLE_*_BACKGROUND`/`_NEUTRON_STAR_BACKGROUND` constants, now
    retired from here) draw those from `kind` alone. Carries `href`, the
    phenomenon's own detail page, the same as `_star_data` does for a star
    system -- `sectorscene.js`'s info panel links there.

    Args:
        link_url (callable): Builds `href` (see the module docstring).
        phenomenon (dict): One entry from `queryDb.phenomena_near_sector`
                           (`id`, `type`, `name`, `descriptor`, `radius_ly`,
                           `distance_ly`, and `home`: `False` marks the
                           entry `neighbor`).
        x_px, y_px, z_px (float): Already-normalized scene-space position.
        radius_px (float): This cloud's drawn radius (`_phenomenon_cloud_radius_px`).

    Returns:
        dict: The cloud's scene entry.
    """
    phenomenon_type = phenomenon["type"]
    descriptor = phenomenon["descriptor"] or ""

    data = {
        "x": x_px, "y": y_px, "z": z_px, "r": radius_px,
        "name": phenomenon["name"],
        "radiusText": format_distance_ly(phenomenon["radius_ly"]) if phenomenon["radius_ly"] else None,
        "distanceText": f'{format_distance_ly(phenomenon["distance_ly"])} from sector center',
        "href": link_url("phenomenon", phenomenon_type=phenomenon["type"], phenomenon_id=phenomenon["id"]),
        # What the sector page's Contents "Show on map" buttons name.
        "key": f'{phenomenon["type"]}:{phenomenon["id"]}',
    }
    if not phenomenon.get("home", True):
        # A neighboring sector's cloud reaching into this one (MAP.45):
        # sectorscene.js draws it fainter, and the info panel says so.
        data["neighbor"] = True
        data["distanceText"] += " (from a neighboring sector)"

    if phenomenon_type == "nebula":
        color = _NEBULA_TYPE_COLORS.get(descriptor, _DEFAULT_NEBULA_COLOR)
        core_alpha, edge_alpha = _NEBULA_TYPE_ALPHA.get(descriptor, _DEFAULT_NEBULA_ALPHA)
        data["kind"] = "nebula"
        # sectorscene.js swaps the sphere for the nebula's shape (MAP.103).
        data["nebulaId"] = phenomenon["id"]
        data["coreColor"] = f"{color}{core_alpha:02x}"
        data["edgeColor"] = f"{color}{edge_alpha:02x}"
        data["typeLabel"] = f'{descriptor.capitalize()} Nebula'
    elif phenomenon_type == "asteroid_field":
        data["kind"] = "asteroidField"
        data["typeLabel"] = f'Asteroid Field ({descriptor.capitalize()})'
    elif phenomenon_type == "black_hole":
        data["kind"] = "blackHoleAccreting" if descriptor == "accreting" else "blackHoleQuiescent"
        data["typeLabel"] = f'Black Hole ({descriptor.capitalize()})' if descriptor else "Black Hole"
    elif phenomenon_type == "neutron_star":
        data["kind"] = "neutronStar"
        data["typeLabel"] = f'Neutron Star ({descriptor.replace("-", " ").capitalize()})' if descriptor else "Neutron Star"
    elif phenomenon_type == "supernova_remnant":
        data["kind"] = "supernovaRemnant"
        data["typeLabel"] = f'Supernova Remnant ({descriptor.capitalize()})' if descriptor else "Supernova Remnant"
    elif phenomenon_type == "rogue_planet":
        data["kind"] = "roguePlanet"
        # GEN.8: its planet class first, e.g. "Rogue Planet (Class C, terrestrial)".
        bits = [f'Class {phenomenon["class"]}' if phenomenon.get("class") else None,
                descriptor.capitalize() if not phenomenon.get("class") else descriptor]
        bits = [bit for bit in bits if bit]
        data["typeLabel"] = f'Rogue Planet ({", ".join(bits)})' if bits else "Rogue Planet"
    elif phenomenon_type == "interstellar_comet":
        data["kind"] = "interstellarComet"
        data["typeLabel"] = f'Interstellar Comet ({descriptor.capitalize()})' if descriptor else "Interstellar Comet"
    elif phenomenon_type == "quasar":
        data["kind"] = "quasar"
        data["typeLabel"] = f'Quasar ({descriptor.capitalize()})' if descriptor else "Quasar"
    else:
        # Defensive fallback for a future phenomenon type this function
        # doesn't know about yet -- drawn the same as a default-colored
        # nebula rather than a crash or a silently wrong label.
        data["kind"] = "nebula"
        data["coreColor"] = f"{_DEFAULT_NEBULA_COLOR}90"
        data["edgeColor"] = f"{_DEFAULT_NEBULA_COLOR}30"
        data["typeLabel"] = phenomenon_type.replace("_", " ").title()

    light = _PHENOMENON_LIGHTS.get(data["kind"])
    if light is not None:
        # Drawn as a point of light rather than a sphere (MAP.15).
        data["light"] = dict(light)
    elif data["kind"] == "roguePlanet":
        data["light"] = dict(_ROGUE_LIGHT)
        data["markedLight"] = dict(_ROGUE_MARKED_LIGHT)
    return data


def map_scene_data(
    link_url, edge_mpc, address, center_pc, systems, phenomena=None, neighbors=None, generate=None,
):
    """
    The Sector Map's scene JSON (the `#starmap-data` block `map_scene_data`
    embeds, and what `/sector/<id>/scene` answers for the Galaxy Map's
    sector stage, MAP.66). The arguments are `map_scene_data`'s; see it.

    Returns:
        dict: `sceneHalfPx`, `defaultZoom`, `lyPerPxAtZoom1`, `outline`,
            `compass`, `stars`, `clouds`, `neighbors`, `generate`,
            `edgeLy`, `centerPc` and `halfEdgePc`.
    """
    half_edge = (edge_mpc / 2) if edge_mpc else 1.0

    stars_data = []
    # Every plotted star/cloud's distance from the scene's own center --
    # kept separate from the outline's own extent (see `_outline_data`)
    # rather than one combined list, so `_default_zoom` can fit *this*,
    # real content on its own whenever there is any. A core-ring sector's
    # pie-wedge cell can be wider than `_SCENE_HALF_PX` (see
    # `_cell_edges_px`), and a combined extent list would let that shape
    # force the whole default zoom down to its own floor
    # (`_MIN_DEFAULT_ZOOM`) even when every star fits. The cell/cube
    # shape only ever decides the default zoom when there's no real
    # content to fit instead (see `default_zoom` below).
    extent_radii_px = []
    for system in systems:
        # +y is "up" on screen; the stored convention increases downward,
        # matching the old CSS scene's own layout axes, so the sign flips
        # here once, at the one place normalized position becomes a scene
        # coordinate -- everything downstream (including the binary
        # offset below) works in already-screen-oriented units. No depth
        # sort is needed here -- the client's real depth buffer handles
        # occlusion.
        galaxy_x, galaxy_y, galaxy_z = _rotate_to_galaxy_frame(
            center_pc, (system["x"] or 0, system["y"] or 0, system["z"] or 0)
        )
        nx = max(-1.05, min(1.05, galaxy_x / half_edge))
        ny = max(-1.05, min(1.05, galaxy_y / half_edge))
        nz = max(-1.05, min(1.05, galaxy_z / half_edge))
        x_px = nx * _SCENE_HALF_PX
        y_px = -ny * _SCENE_HALF_PX
        z_px = nz * _SCENE_HALF_PX
        extent_radii_px.append(math.sqrt(x_px * x_px + y_px * y_px + z_px * z_px))

        stars = system["stars"]
        is_binary = len(stars) > 1
        primary_index = len(stars_data)
        primary_entry = _star_data(link_url, system, stars[0], x_px, y_px, z_px)
        if is_binary:
            primary_entry["name"] = system["name"]  # the system is what gets picked
        stars_data.append(primary_entry)

        if is_binary:
            primary_r = _star_dot_radius(stars[0]["radius_km"])
            offset = primary_r * _BINARY_OFFSET_FRACTION
            companion = _star_data(
                link_url, system, stars[1], x_px + offset, y_px + offset, z_px,
                max_r=primary_r * _SECONDARY_MAX_RATIO,
            )
            # A binary is one pickable system (MAP.136): the map resolves a
            # pick of the companion to this entry.
            companion["companionOf"] = primary_index
            stars_data.append(companion)

    clouds_data = []
    for phenomenon in (phenomena or []) if ly_to_milliparsecs is not None else ():
        # Already galaxy-frame (see this function's own `phenomena`
        # docstring) -- no `_rotate_to_galaxy_frame` step, unlike a
        # system's sector-local x/y/z above. Position is intentionally
        # NOT clamped the way a star system's normalized position is
        # (+-1.05): a neighbor's cloud is allowed to sit mostly outside
        # this sector's own cell (that's the whole point of
        # `phenomena_near_sector`'s bounding-sphere overlap test) and/or
        # be far larger than it.
        nx = ly_to_milliparsecs(phenomenon["offset_x_ly"]) / half_edge
        ny = ly_to_milliparsecs(phenomenon["offset_y_ly"]) / half_edge
        nz = ly_to_milliparsecs(phenomenon["offset_z_ly"]) / half_edge
        x_px = nx * _SCENE_HALF_PX
        y_px = -ny * _SCENE_HALF_PX
        z_px = nz * _SCENE_HALF_PX
        radius_px = _phenomenon_cloud_radius_px(phenomenon["radius_ly"], half_edge)
        if radius_px:
            clouds_data.append(_cloud_data(link_url, phenomenon, x_px, y_px, z_px, radius_px))
            # The cloud's own edge, not just its center -- a large nebula
            # can dwarf the scene (see `_MAX_CLOUD_RADIUS_PX`), and its
            # center alone would understate how far out it actually reaches.
            extent_radii_px.append(
                math.sqrt(x_px * x_px + y_px * y_px + z_px * z_px) + radius_px
            )

    outline, shape_extent_radii_px = _outline_data(address, edge_mpc, half_edge)

    compass = _compass_data(center_pc)
    # Fit the real content (star systems/clouds) whenever there is any --
    # see `extent_radii_px`'s own comment above for why the sector shape's
    # own, potentially much larger extent must NOT also be allowed to drag
    # that zoom down with it. Only an entirely empty sector (nothing
    # plotted at all) falls back to fitting that shape's own extent
    # instead, so the camera still opens at a sensible distance rather than
    # zoomed to fit nothing at all.
    default_zoom = _default_zoom(extent_radii_px or shape_extent_radii_px)

    ly_per_px = _ly_per_px_at_zoom_1(half_edge)

    neighbors_data = [
        entry for entry in (_neighbor_indicator_data(link_url, neighbor) for neighbor in (neighbors or []))
        if entry is not None
    ]

    scene_data = {
        "sceneHalfPx": _SCENE_HALF_PX,
        "defaultZoom": default_zoom,
        "lyPerPxAtZoom1": ly_per_px,
        "outline": outline,
        "compass": compass,
        "stars": stars_data,
        "clouds": clouds_data,
        # Where a nebula's mesh comes from (`{id}` is its id; MAP.103).
        "nebulaShapePath": NEBULA_SHAPE_PATH,
        "neighbors": neighbors_data,
        "generate": generate,
        # The sector's own edge in light years, for the Generate buttons'
        # "up to about N sectors" estimate (static/generatebuttons.js).
        "edgeLy": milliparsecs_to_ly(edge_mpc) if edge_mpc and milliparsecs_to_ly else None,
        # Where the sector is, for drawing its scene inside the Galaxy Map's
        # (`static/sectorscene.js`): its center in galaxy-frame parsecs
        # and half its edge in parsecs, which `sceneHalfPx` scene units span.
        "centerPc": list(center_pc) if center_pc else None,
        "halfEdgePc": (edge_mpc / 2) / 1000 if edge_mpc else None,
    }
    return scene_data


UNCHARTED_LABEL = "Uncharted"
"""str: The mark an opened sector nothing was generated in carries (MAP.162):
in its title, its info panel and the Galaxy Map's frame."""

# A scattered black hole is a purple point and a neutron star a dark blue
# one, both lit so they show before their sector is made (MAP.164); a
# quasar, a supernova remnant and a planetary nebula keep the looks every
# map gives them.
_UNCHARTED_LIGHTS = {
    "black_hole": {"color": "#a06bff", "corePx": 4.0, "sizePx": 26.0, "glow": 0.7, "bright": 1.0, "whiten": 0.0},
    "neutron_star": {"color": "#3f63e0", "corePx": 3.0, "sizePx": 22.0, "glow": 0.7, "bright": 1.0, "whiten": 0.0},
}

# What a scatter row has no size for: its drawn radius, in light years.
_UNCHARTED_RADIUS_LY = {"nebula": 1.5, "supernova_remnant": 1.0}
_UNCHARTED_COMPACT_RADIUS_LY = 0.02

# A hypervelocity star is a bright blue-white star ejected from the core.
_HYPERVELOCITY_STAR = {"star_type": "B", "temperature_k": 15000.0, "luminosity_sol": 2000.0, "radius_sol": 5.0}

_UNCHARTED_LABELS = {
    "black_hole": "Black hole", "neutron_star": "Neutron star", "nebula": "Planetary nebula",
    "supernova_remnant": "Supernova remnant", "quasar": "Quasar", "hypervelocity_star": "Hypervelocity star",
}


def _uncharted_position_px(position_pc, center_pc, half_edge_pc):
    """A galaxy-frame position (parsecs) as scene units about the sector's center."""
    nx, ny, nz = ((position_pc[i] - center_pc[i]) / half_edge_pc for i in range(3))
    return nx * _SCENE_HALF_PX, -ny * _SCENE_HALF_PX, nz * _SCENE_HALF_PX


def _uncharted_star_entry(star, x_px, y_px, z_px, name):
    """One waiting star (a `bright_stars` row, or `_HYPERVELOCITY_STAR`) as a scene star
    entry; it has no system page, so no `href`."""
    luminosity_w = star["luminosity_sol"] * SOLAR_LUMINOSITY
    radius_km = star["radius_sol"] * _SUN_RADIUS_KM if star.get("radius_sol") else None
    temperature_k = star["temperature_k"] or 5778.0
    fill, stroke = star_color(star["star_type"], temperature_k, luminosity_w)
    return {
        "x": x_px, "y": y_px, "z": z_px, "r": _star_dot_radius(radius_km),
        "light": _star_light(luminosity_w, radius_km, temperature_k),
        "fill": fill, "stroke": stroke, "name": name, "starType": star["star_type"],
        "luminositySol": star["luminosity_sol"], "temp": f"{int(temperature_k)} K", "uncharted": True,
    }


def uncharted_scene_data(address, designation, center_pc, edge_pc, contents, generate=None):
    """
    The scene JSON of a sector nothing was generated in (MAP.162): what the
    scatters left in the cell, drawn like a generated sector's contents so
    the Galaxy Map can open it in place, with `uncharted` (True) and
    `designation` for the mark it carries. Same keys as `map_scene_data`
    (no outline or neighbors: the Galaxy Map draws the cell itself).

    Args:
        address (tuple): The cell's `(ring, layer, slot)`.
        designation (str): Its galaxy designation.
        center_pc (sequence): Its center, galaxy-frame parsecs.
        edge_pc (float): The grid's sector edge.
        contents (dict): `queryDb.uncharted_sector_contents`: `stars` and `scattered`.
        generate (dict, optional): `generate_target`, for the Generate buttons.
    """
    half_edge_pc = edge_pc / 2.0
    stars, clouds = [], []
    for star in contents["stars"]:
        x_px, y_px, z_px = _uncharted_position_px((star["x"], star["y"], star["z"]), center_pc, half_edge_pc)
        stars.append(_uncharted_star_entry(
            star, x_px, y_px, z_px, f'{star["star_type"]} star ({UNCHARTED_LABEL.lower()})'))
    for row in contents["scattered"]:
        x_px, y_px, z_px = _uncharted_position_px((row["x"], row["y"], row["z"]), center_pc, half_edge_pc)
        if row["type"] == "hypervelocity_star":
            stars.append(_uncharted_star_entry(
                _HYPERVELOCITY_STAR, x_px, y_px, z_px, f"Hypervelocity star ({UNCHARTED_LABEL.lower()})"))
            continue
        kind_label = _UNCHARTED_LABELS.get(row["type"], row["type"].replace("_", " ").capitalize())
        radius_ly = _UNCHARTED_RADIUS_LY.get(row["type"], _UNCHARTED_COMPACT_RADIUS_LY)
        distance_pc = math.sqrt(sum((row[axis] - center_pc[i]) ** 2 for i, axis in enumerate("xyz")))
        phenomenon = {
            "type": row["type"], "id": f's{row["id"]}', "name": f"{kind_label} ({UNCHARTED_LABEL.lower()})",
            "descriptor": "planetary" if row["type"] == "nebula" else (row["subtype"] or ""),
            "radius_ly": radius_ly, "distance_ly": milliparsecs_to_ly(distance_pc * 1000.0),
        }
        cloud = _cloud_data(lambda *args, **kwargs: None, phenomenon, x_px, y_px, z_px,
                            _phenomenon_cloud_radius_px(radius_ly, half_edge_pc * 1000.0))
        cloud.pop("href", None)
        cloud.pop("nebulaId", None)        # a scatter row has no mesh to fetch
        cloud["uncharted"] = True
        light = _UNCHARTED_LIGHTS.get(row["type"])
        if light is not None:
            cloud["light"] = dict(light)
        clouds.append(cloud)
    return {
        "sceneHalfPx": _SCENE_HALF_PX, "defaultZoom": 1.0, "lyPerPxAtZoom1": None, "outline": None, "compass": None,
        "stars": stars, "clouds": clouds, "nebulaShapePath": NEBULA_SHAPE_PATH, "neighbors": [],
        "generate": generate, "edgeLy": milliparsecs_to_ly(edge_pc * 1000.0), "centerPc": list(center_pc), "halfEdgePc": half_edge_pc,
        "uncharted": True, "designation": designation, "address": list(address),
    }
