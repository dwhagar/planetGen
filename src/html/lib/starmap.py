# html/lib/starmap.py

"""
Interactive 3D sector starmap: every placed star system in a sector (two,
overlapping, for a binary), plus every nearby nebula/asteroid field/black
hole/neutron star, rendered as a real WebGL scene (`static/sectormap.js`,
via three.js -- vendored at `static/vendor/three.module.min.js`, see that
directory's `THIRD_PARTY_NOTICES.txt`) instead of the CSS
`transform-style: preserve-3d` scene this module used to build directly as
HTML `<div>`s. This module's job is now just the data: every position,
size, and color is still computed here exactly as before (this file owns
all of it -- the client only ever positions/colors/labels what it's handed,
it makes no astrophysical or layout decisions of its own), serialized as
JSON into one `<script type="application/json">` block `sectormap.js`
reads on page load, alongside a `<canvas>` for it to render into and a
`<noscript>` fallback list (plain links, no map) for a browser that can't
run it.

Every star's (x, y, z) -- rotated once, server-side, from its own
sector-local axes into the galaxy frame when the sector has a galaxy
placement (`_rotate_to_galaxy_frame`, so it agrees with the wedge outline/
compass arrow below, both already galaxy-frame quantities) -- is expressed
in the same fixed pixel-scale world units as before (`_SCENE_HALF_PX` per
half the sector's own edge), so a change here needs no matching change to
how the client interprets a position: it's still "this many units from the
scene's own center," just handed to a real perspective camera instead of a
flat CSS transform now.

`_outline_data` still computes the scene's own bounding shape -- one of
two, depending on whether the sector has a galaxy placement: a sector with
`shell_index`/`shell_slot_index` set gets its real, approximate on-shell
wedge (see `sector_wedge_vertices_pc` -- bounded by the shell's own radial
thickness and roughly how much angular "real estate" this slot's
Fibonacci placement owns among its neighbors); one without a placement
falls back to a plain axis-aligned cube of edge `edge_mpc`. `sectormap.js`
no longer draws this shape as a wireframe (dropped as visual clutter), but
its extent still drives `_default_zoom`'s fallback framing for a sector
with nothing else plotted in it (see that function's own docstring).

Clicking a star/cloud (or activating one of the `<noscript>`/accessible-
fallback list's own buttons) doesn't navigate straight to `system.py`/
`phenomenon.py` -- it populates the info side panel first, so a click
shows details before the panel's own link is what navigates away.
"""

import colorsys
import json
import math

from fmt import esc, post_link

try:
    from stellarObjects.physical_constants import SPECTRAL_CLASS_COLORS, TEMP_RANGES, SOLAR_LUMINOSITY, SOLAR_RADIUS_M
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
    from stellarObjects.galaxyGeometry import sector_position_pc, sector_wedge_vertices_pc
    from stellarObjects.sectorGeometry import cube_orientation
    from stellarObjects.utils import ly_to_milliparsecs, milliparsecs_to_ly, mpc_to_pc, pc_to_mpc
except ImportError:
    # Same deployment gap as above -- without these, a sector with no
    # galaxy placement (or one this deployment can't reach the geometry
    # module for) just keeps the plain axis-aligned cube outline and no
    # scale bar; see `_wedge_edges_px`/`_ly_per_px_at_zoom_1` below. Without
    # `cube_orientation` specifically, star dots fall back to being
    # plotted in their own unrotated local axes (see `_rotate_to_galaxy_frame`).
    # Without `ly_to_milliparsecs`, nebula/asteroid-field clouds can't be
    # scaled/positioned at all and are omitted.
    sector_position_pc = None
    sector_wedge_vertices_pc = None
    cube_orientation = None
    mpc_to_pc = None
    pc_to_mpc = None
    milliparsecs_to_ly = None
    ly_to_milliparsecs = None

# The 3D scene's world-unit scale -- every position/radius below is in
# these units (a sector's own half-edge maps to this many of them), handed
# to `sectormap.js` as-is; a real perspective camera doesn't otherwise care
# what unit "1" means, unlike the old CSS version where this was also a
# literal pixel count.
_SCENE_HALF_PX = 160.0

# A plain axis-aligned fallback cube's own corners, at distance
# `_SCENE_HALF_PX` out along all three axes at once -- used by
# `_default_zoom` as the fallback shape's extent when there's no wedge
# wireframe to measure instead.
_CUBE_CORNER_RADIUS_PX = _SCENE_HALF_PX * math.sqrt(3)

# Never start a sector further zoomed out than this, however large its
# wedge/cube/plotted content gets (an extreme far-out placement could
# otherwise compute an unusably tiny initial view) -- `sectormap.js`'s own
# zoom-out control/scroll remains available past this floor regardless.
_MIN_DEFAULT_ZOOM = 0.2


def _default_zoom(extent_radii_px):
    """
    Picks the zoom level the sector map should *open* at, so its real
    content -- the wedge/cube outline, every plotted star/cloud -- actually
    fits in frame on first paint, instead of always starting at a flat
    `zoom = 1` regardless of how much bigger the real shape is.

    A galaxy-placed sector's wedge wireframe is deliberately allowed to
    extend well past `_SCENE_HALF_PX` (see `_wedge_edges_px`'s docstring --
    the wedge's angular patch doesn't coincide with a cube's flat sides),
    and even the plain-cube fallback's own corners sit `_SCENE_HALF_PX *
    sqrt(3)` out -- both already past a "fits at zoom 1" frame.
    `sectormap.js`'s camera distance at zoom 1 is calibrated so a sphere of
    radius `_SCENE_HALF_PX` exactly fills the frame (see its own
    `_referenceDistance`), so this fraction is what that client-side
    calibration is relative to.

    Args:
        extent_radii_px (list[float]): Distance from the scene's own
                                       center, for every point that matters
                                       (wedge/cube vertices, star/cloud
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
_MIN_DOT_R = 3.0
_MAX_DOT_R = 14.0

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

# Which pairs of 8-point vertex-list indices form a cube (or wedge)'s 12
# edges -- indices are `4*a + 2*b + c` for whichever three binary choices
# the 8 corners vary over (see `sector_wedge_vertices_pc`'s docstring for
# the wedge's own r/phi/theta bits; `_cube_corners_px` uses the same
# convention for its x/y/z sign bits below), so two vertices share an edge
# exactly when their indices differ in a single bit -- this is just a
# cube's own edge topology, reused for the wedge's 8 corners regardless of
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
    wireframe reads just as clearly and lets both this and the wedge case
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
    axes instead, via `sectorGeometry.cube_orientation`'s fixed convention
    (radial-outward local `+Z`, projected-galactic-north local `+X`; see
    `docs/design/galaxy-coordinate-system.md`'s "Cube orientation"
    section). This is the same convention `sectorGeometry.py` already
    uses as the tangent-plane basis for this sector's own wedge vertices
    -- computed once per sector there at generation time (baked into the
    stored vertex positions), and once per render here, from the stored
    `center_pc` alone, with no extra state of its own.

    Without this, a star dot's local (x, y, z) was plotted as if it were
    already expressed in the galaxy frame -- consistent with itself, but
    not with the wedge outline or the "Galactic Center" compass arrow
    (`_wedge_edges_px`/`_compass_data`), which read the sector's *actual*
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
    if center_pc is None or any(c is None for c in center_pc) or cube_orientation is None:
        return local_vec

    axis_x, axis_y, axis_z = cube_orientation(center_pc)
    lx, ly, lz = local_vec
    return (
        lx * axis_x[0] + ly * axis_y[0] + lz * axis_z[0],
        lx * axis_x[1] + ly * axis_y[1] + lz * axis_z[1],
        lx * axis_x[2] + ly * axis_y[2] + lz * axis_z[2],
    )


def _wedge_edges_px(shell_index, shell_slot_index, edge_mpc, half_edge):
    """
    Computes the sector's approximate on-shell wedge (see
    `sector_wedge_vertices_pc`) as 8 `(x, y, z)` scene-space unit points,
    in the same coordinate convention `render_map_panel` uses for star
    dots (sector-center-relative parsecs -> milliparsecs -> normalized by
    `half_edge` -> scene units, with the y-axis flipped since +y is "up"
    on screen but down in database/layout convention) -- so the wireframe
    and the star dots it surrounds always share one consistent frame.

    Unlike a star dot's normalized position, these are deliberately *not*
    clamped to the +-1.05 the cube fallback uses: the whole point of this
    shape is that it doesn't stay inside the sector's own edge_mpc cube,
    since the shell's angular patch and the cube's flat sides don't
    coincide (see `sector_wedge_vertices_pc`'s docstring).

    Returns:
        list[tuple] or None: 8 `(x, y, z)` points, or `None` if this
                             sector has no galaxy placement (no
                             `shell_index`/`shell_slot_index`) or the
                             geometry helpers aren't importable in this
                             deployment -- either way, callers fall back
                             to the plain cube.
    """
    if (
        shell_index is None
        or shell_slot_index is None
        or not edge_mpc
        or sector_wedge_vertices_pc is None
        or sector_position_pc is None
        or mpc_to_pc is None
        or pc_to_mpc is None
    ):
        return None

    edge_pc = mpc_to_pc(edge_mpc)
    center_pc = sector_position_pc(shell_index, shell_slot_index, edge_pc)
    vertices_pc = sector_wedge_vertices_pc(shell_index, shell_slot_index, edge_pc)

    points = []
    for vx, vy, vz in vertices_pc:
        nx = pc_to_mpc(vx - center_pc[0]) / half_edge
        ny = pc_to_mpc(vy - center_pc[1]) / half_edge
        nz = pc_to_mpc(vz - center_pc[2]) / half_edge
        points.append((nx * _SCENE_HALF_PX, -ny * _SCENE_HALF_PX, nz * _SCENE_HALF_PX))
    return points


def _outline_data(shell_index, shell_slot_index, edge_mpc, half_edge):
    """
    Builds the scene's outline as `{"kind": "wedge"|"cube", "edges": [...]}`
    -- 12 `[point, point]` pairs either way (see `_EDGE_PAIRS`), from the
    sector's real on-shell wedge vertices (`_wedge_edges_px`) when it has a
    galaxy placement, or the plain axis-aligned fallback cube's own corners
    (`_cube_corners_px`) otherwise.

    Returns:
        tuple[dict, list[float]]: The outline data, and the extent (from
                                  the scene's own center) of every one of
                                  its vertices -- for `_default_zoom`'s
                                  empty-scene fallback.
    """
    wedge_vertices = _wedge_edges_px(shell_index, shell_slot_index, edge_mpc, half_edge)
    if wedge_vertices is not None:
        vertices, kind = wedge_vertices, "wedge"
    else:
        vertices, kind = _cube_corners_px(_SCENE_HALF_PX), "cube"

    edges = [[list(vertices[i]), list(vertices[j])] for i, j in _EDGE_PAIRS]
    extent_radii_px = [math.sqrt(vx * vx + vy * vy + vz * vz) for vx, vy, vz in vertices]
    return {"kind": kind, "edges": edges}, extent_radii_px


_COMPASS_ARROW_REACH = 1.15
"""How far the compass arrow reaches past `_SCENE_HALF_PX` -- past 1.0 so
its tip clears a full-size cube/wedge instead of ending right at (or
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

    This arrow and `_outline_data`'s wedge outline are both computed
    directly from galaxy-frame quantities (`sectors.center_x/y/z_pc`,
    `sector_wedge_vertices_pc`), so they were always correct on their own
    terms. Star dots are rotated into that same frame at render time
    instead (`_rotate_to_galaxy_frame`), so this arrow, the wedge outline,
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


def _ly_per_px_at_zoom_1(half_edge):
    """
    Light-years per world unit at zoom factor 1 -- an exact ratio derived
    straight from `half_edge` (half of `edge_mpc`, the sector's real,
    stored size) mapping to `_SCENE_HALF_PX` units, the same normalizing
    divisor every star dot and wedge vertex on this map is already placed
    by. `sectormap.js`'s scale-bar legend divides this by the live zoom
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
    range). Used only as a fallback for `_star_color` when `star_type`
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


def _star_color(star_type, temperature_k, luminosity_w):
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


def _star_dot_radius(radius_km):
    """Maps a star's physical radius to a dot radius in scene units --
    square-root scaled against the Sun's radius (linear scaling would make
    red dwarfs invisible next to giants, which differ by 2+ orders of
    magnitude in radius_km) and clamped so the map stays legible at either
    extreme."""
    if not radius_km or radius_km <= 0:
        return _MIN_DOT_R
    ratio = radius_km / _SUN_RADIUS_KM
    dot_r = _MIN_DOT_R + 4.0 * math.sqrt(ratio)
    return max(_MIN_DOT_R, min(_MAX_DOT_R, dot_r))


def _star_data(db_name, system, star, x_px, y_px, z_px, label_suffix, max_r=None):
    """
    Builds one star's plain-dict scene entry -- `sectormap.js` draws it as
    a real, textured, glowing 3D sphere (`static/bodyRendering.js`'s
    granulation texture, tinted to `fill` below), unlike the old CSS
    version's `_dot_html`, which had to counter-rotate a flat disc by hand
    every frame to fake a billboard.
    """
    dot_r = _star_dot_radius(star["radius_km"])
    if max_r is not None:
        dot_r = min(dot_r, max_r)
    fill, stroke = _star_color(star["star_type"], star["temperature_k"], star["luminosity_w"])
    return {
        "x": x_px, "y": y_px, "z": z_px, "r": dot_r,
        "fill": fill, "stroke": stroke,
        "name": f'{system["name"]}{label_suffix}',
        "starType": star["star_type"],
        "temp": star["temp_display"],
        "quadrant": system["quadrant"],
        "location": system["location"],
        # Read by `sectormap.js`'s info panel and passed straight to
        # `window.planetgenSubmitNav` (`static/navform.js`) -- a plain
        # object, not `fmt.data_nav_params`'s escaped-JSON-*string* form,
        # since this is embedded directly into `_json_script`'s own JSON
        # rather than an HTML attribute (see that convention's other use,
        # `lib/galaxymap.py`/`lib/navmap.py`'s `<svg>` markers, which --
        # unlike this scene -- really is HTML `data-nav-params`).
        "navTarget": "system.py",
        "navParams": {"db": db_name, "id": system["id"]},
    }


# A cloud's own radius (`radius_ly`) is frequently far larger than a
# single sector -- a big emission nebula can dwarf the whole scene -- so
# unlike a star dot's radius (always tiny next to the scene), this needs
# its own generous cap: large enough to visibly engulf/overflow the frame
# without an unbounded sphere size for a pathological radius_ly value.
_MAX_CLOUD_RADIUS_PX = 6 * (2 * _SCENE_HALF_PX)

# Fill color (and base opacity, baked into the alpha channel below) per
# `nebulae.nebula_type` -- not spectral-accurate the way `_star_color` is
# (a nebula's visible color really does vary this much by type: emission
# nebulae genuinely glow reddish-pink from ionized hydrogen's H-alpha
# line, reflection nebulae blue from scattered starlight, planetary
# nebulae teal/cyan from doubly-ionized oxygen, and a dark nebula is by
# definition an opaque silhouette, not a glow at all -- hence its own
# near-black, higher-opacity fill instead of a lighter translucent one).
_NEBULA_TYPE_COLORS = {
    "emission": "#ff6f91",
    "reflection": "#6fa8ff",
    "planetary": "#5be8c9",
    "dark": "#1c1c24",
}
_NEBULA_TYPE_ALPHA = {
    "emission": (0xB0, 0x40),
    "reflection": (0xA0, 0x38),
    "planetary": (0xA8, 0x3c),
    "dark": (0xE8, 0x90),
}
_DEFAULT_NEBULA_COLOR = "#c9a8e0"
_DEFAULT_NEBULA_ALPHA = (0x90, 0x38)

# Every other phenomenon kind's own look is a fixed recipe (never a
# function of per-instance data beyond which kind/descriptor it is), so
# unlike a nebula's per-instance core/edge color, `sectormap.js` owns those
# recipes directly (see its own `CLOUD_KIND_RECIPES`) -- this module only
# ever needs to say *which* fixed kind a given phenomenon is.


def _phenomenon_cloud_radius_px(radius_ly, half_edge):
    if ly_to_milliparsecs is None or not half_edge:
        return 0.0
    radius_mpc = ly_to_milliparsecs(radius_ly)
    radius_px = (radius_mpc / half_edge) * _SCENE_HALF_PX
    return max(4.0, min(_MAX_CLOUD_RADIUS_PX, radius_px))


def _cloud_data(db_name, phenomenon, x_px, y_px, z_px, radius_px):
    """
    Builds one phenomenon's plain-dict scene entry. A nebula gets its own
    per-instance `coreColor`/`edgeColor` (two gradient stops, `#rrggbbaa`)
    computed from `_NEBULA_TYPE_COLORS`/`_NEBULA_TYPE_ALPHA`; every other
    kind carries no color data at all -- `sectormap.js`'s own fixed
    recipes (mirroring this module's old `_ASTEROID_FIELD_BACKGROUND`/
    `_BLACK_HOLE_*_BACKGROUND`/`_NEUTRON_STAR_BACKGROUND` constants, now
    retired from here) draw those from `kind` alone. Carries `navTarget`/
    `navParams` to `phenomenon.py` (this project's detail page for a
    standalone phenomenon), the same plain-object convention `_star_data`
    uses for a star system -- `sectormap.js`'s info panel passes either
    straight to `window.planetgenSubmitNav` (`static/navform.js`) so
    following it never puts `phenomenon.py?...` in the browser's own
    address bar.

    Args:
        db_name (str): The current `?db=` value, for `navParams`.
        phenomenon (dict): One entry from `queryDb.phenomena_near_sector`
                           (`id`, `type`, `name`, `descriptor`, `radius_ly`,
                           `distance_ly`).
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
        "radiusText": f'{phenomenon["radius_ly"]:,.2f} ly',
        "distanceText": f'{phenomenon["distance_ly"]:,.1f} ly from sector center',
        "navTarget": "phenomenon.py",
        "navParams": {"db": db_name, "type": phenomenon["type"], "id": phenomenon["id"]},
    }

    if phenomenon_type == "nebula":
        color = _NEBULA_TYPE_COLORS.get(descriptor, _DEFAULT_NEBULA_COLOR)
        core_alpha, edge_alpha = _NEBULA_TYPE_ALPHA.get(descriptor, _DEFAULT_NEBULA_ALPHA)
        data["kind"] = "nebula"
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
    else:
        # Defensive fallback for a future phenomenon type this function
        # doesn't know about yet -- drawn the same as a default-colored
        # nebula rather than a crash or a silently wrong label.
        data["kind"] = "nebula"
        data["coreColor"] = f"{_DEFAULT_NEBULA_COLOR}90"
        data["edgeColor"] = f"{_DEFAULT_NEBULA_COLOR}30"
        data["typeLabel"] = phenomenon_type.replace("_", " ").title()

    return data


def _json_script(data):
    """
    Serializes `data` for safe embedding inside a `<script
    type="application/json">` block: escapes `<`, `>`, and `&` as Unicode
    escapes (the standard "JSON in an HTML script tag" mitigation, e.g.
    Django's `json_script`) so a database value containing `</script>` (a
    system/phenomenon/sector name is arbitrary user-supplied text -- see
    `--name`) can't break out of the tag, despite this module no longer
    running every such value through `fmt.esc` the way its HTML-attribute-
    building predecessor did.
    """
    return (
        json.dumps(data)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
    )


def _noscript_list_html(db_name, systems, phenomena):
    """
    A plain, always-present (no JS required) list of links -- the
    `<noscript>` fallback for a browser that can't run the WebGL scene
    `sectormap.js` builds, so the sector's own systems/phenomena are still
    reachable rather than the panel being entirely blank without
    JavaScript. Not a substitute for the map itself (no position/size/
    color -- just names and links), same spirit as any other progressive-
    enhancement fallback list. Built from `fmt.post_link` -- a real
    `<form>` submit button, not a plain `<a href>` -- same as every other
    in-app link now, so `db`/an id doesn't show up in the address bar even
    here, and it needs no JavaScript of its own to work either.
    """
    items = []
    for system in systems:
        items.append(f'<li>{post_link("system.py", {"db": db_name, "id": system["id"]}, esc(system["name"]))}</li>')
    for phenomenon in (phenomena or []):
        link = post_link(
            "phenomenon.py", {"db": db_name, "type": phenomenon["type"], "id": phenomenon["id"]},
            esc(phenomenon["name"]),
        )
        items.append(f'<li>{link}</li>')
    if not items:
        return ""
    return f'<noscript><ul class="starmap-noscript-list">{"".join(items)}</ul></noscript>'


def render_map_panel(db_name, edge_mpc, shell_index, shell_slot_index, center_pc, systems, phenomena=None):
    """
    Builds the "Sector Map" panel: a `<canvas>` `sectormap.js` renders an
    interactive WebGL scene into (drag to rotate, scroll/button to zoom,
    click for info), plus a `<script type="application/json">` block
    carrying every position/size/color/label that scene needs -- one entry
    per placed star system (two, overlapping, for a binary -- the primary
    at the system's actual position, the secondary offset down-and-right
    from it) and one per nearby nebula/asteroid field/black hole/neutron
    star -- and an info side panel the same script fills in when something
    is clicked.

    The scene's own bounding shape -- the sector's approximate on-shell
    wedge (see `_outline_data`/`sector_wedge_vertices_pc`) when this sector
    has a galaxy placement (`shell_index`/`shell_slot_index` both set) and
    the geometry helpers are importable, or the plain axis-aligned cube
    otherwise -- is no longer drawn as a wireframe (`sectormap.js` dropped
    it as visual clutter); it's still computed here purely to drive
    `_default_zoom`'s fallback framing for a sector with nothing else
    plotted in it.

    Args:
        db_name (str): The current `?db=` value, used to build each
                       entry's `navParams` (`static/navform.js` posts
                       these to `system.py`/`phenomenon.py` on click, same
                       as every other in-app navigation -- see that
                       file's own docstring).
        edge_mpc (float): The sector's cube edge (`sectors.edge_mpc`) --
                          every system's `position_*_mpc` is relative to
                          the sector's cubic center (see schema.sql's
                          `star_systems` comment), so half of this is the
                          normalizing divisor for each axis.
        shell_index (int or None): The sector's `sectors.shell_index`.
        shell_slot_index (int or None): The sector's
                                        `sectors.shell_slot_index`.
        center_pc (tuple or None): `(center_x_pc, center_y_pc,
                                   center_z_pc)` -- drives the "Galactic
                                   Center" compass arrow (`_compass_data`)
                                   and rotates each system's local
                                   position into the galaxy frame
                                   (`_rotate_to_galaxy_frame`) before it's
                                   plotted, so star dots agree with the
                                   arrow and the wedge outline on one
                                   frame; `None` for a sector with no
                                   galaxy placement, which omits the
                                   arrow and leaves positions unrotated.
        systems (list[dict]): One entry per placed system (position not
                              NULL), each with `id`, `name`, `quadrant`,
                              `location`, `x`/`y`/`z` (the raw
                              `position_*_mpc` columns), and `stars`: a
                              list of 1 dict (single star) or 2 (primary,
                              then secondary), each with `star_type`,
                              `temperature_k`, `radius_km`,
                              `luminosity_w`, `temp_display`.
        phenomena (list[dict] or None): `queryDb.phenomena_near_sector`'s
                              return shape -- every nebula/asteroid field/
                              black hole/neutron star whose sphere could
                              plausibly reach into this sector's cube. Its
                              `offset_x/y/z_ly` are already galaxy-frame
                              (computed directly from two galaxy-frame
                              centers -- see `schema.sql`'s "v18" note), so
                              -- unlike `systems`' sector-local `x`/`y`/`z`
                              -- these are placed directly with no
                              `_rotate_to_galaxy_frame` step. `None`/empty
                              draws no clouds at all.

    Returns:
        str: A complete `<section class="panel">` block.
    """
    half_edge = (edge_mpc / 2) if edge_mpc else 1.0

    stars_data = []
    # Every plotted star/cloud's distance from the scene's own center --
    # kept separate from the outline's own extent (see `_outline_data`)
    # rather than one combined list, so `_default_zoom` can fit *this*,
    # real content on its own whenever there is any. A galaxy-placed
    # sector far out on its shell can have a wedge wireframe many times
    # wider than `_SCENE_HALF_PX` (deliberately allowed to extend well
    # past it -- see `_wedge_edges_px`'s own docstring), and a combined
    # extent list used to let that single huge shape force the *whole*
    # default zoom down to its own floor (`_MIN_DEFAULT_ZOOM`) even when
    # every actual star sat well within frame on its own. The wedge/cube
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
        # "A"/"B", not "Primary"/"Secondary" -- matches how the generator
        # already names the stars themselves (the secondary's own stored
        # name is "<system name> B"; see systemData.StarSystem.__init__),
        # and reads as a real star name rather than an internal role label.
        primary_suffix = " A" if is_binary else ""
        stars_data.append(_star_data(db_name, system, stars[0], x_px, y_px, z_px, primary_suffix))

        if is_binary:
            primary_r = _star_dot_radius(stars[0]["radius_km"])
            offset = primary_r * _BINARY_OFFSET_FRACTION
            stars_data.append(_star_data(
                db_name, system, stars[1], x_px + offset, y_px + offset, z_px, " B",
                max_r=primary_r * _SECONDARY_MAX_RATIO,
            ))

    clouds_data = []
    for phenomenon in (phenomena or []) if ly_to_milliparsecs is not None else ():
        # Already galaxy-frame (see this function's own `phenomena`
        # docstring) -- no `_rotate_to_galaxy_frame` step, unlike a
        # system's sector-local x/y/z above. Position is intentionally
        # NOT clamped the way a star system's normalized position is
        # (+-1.05): a cloud is allowed to sit mostly outside this sector's
        # own cube (that's the whole point of `phenomena_near_sector`'s
        # bounding-sphere overlap test) and/or be far larger than it.
        nx = ly_to_milliparsecs(phenomenon["offset_x_ly"]) / half_edge
        ny = ly_to_milliparsecs(phenomenon["offset_y_ly"]) / half_edge
        nz = ly_to_milliparsecs(phenomenon["offset_z_ly"]) / half_edge
        x_px = nx * _SCENE_HALF_PX
        y_px = -ny * _SCENE_HALF_PX
        z_px = nz * _SCENE_HALF_PX
        radius_px = _phenomenon_cloud_radius_px(phenomenon["radius_ly"], half_edge)
        if radius_px:
            clouds_data.append(_cloud_data(db_name, phenomenon, x_px, y_px, z_px, radius_px))
            # The cloud's own edge, not just its center -- a large nebula
            # can dwarf the scene (see `_MAX_CLOUD_RADIUS_PX`), and its
            # center alone would understate how far out it actually reaches.
            extent_radii_px.append(
                math.sqrt(x_px * x_px + y_px * y_px + z_px * z_px) + radius_px
            )

    outline, shape_extent_radii_px = _outline_data(shell_index, shell_slot_index, edge_mpc, half_edge)

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

    scene_data = {
        "sceneHalfPx": _SCENE_HALF_PX,
        "defaultZoom": default_zoom,
        "lyPerPxAtZoom1": ly_per_px,
        "outline": outline,
        "compass": compass,
        "stars": stars_data,
        "clouds": clouds_data,
    }

    if systems or clouds_data:
        click_hint = "Click a star system or cloud for details." if clouds_data else "Click a star system for details."
        info_panel = (
            '<aside class="starmap-info" id="starmap-info">'
            f'<p class="hint">{click_hint}</p></aside>'
        )
    else:
        info_panel = '<aside class="starmap-info" id="starmap-info"><p class="hint">No systems placed in this sector.</p></aside>'

    scale_bar_html = (
        '<div class="starmap-scale" id="starmap-scale">'
        '<span class="starmap-scale-bar" id="starmap-scale-bar"></span>'
        '<span class="starmap-scale-label" id="starmap-scale-label"></span>'
        "</div>"
        if ly_per_px else ""
    )

    noscript_html = _noscript_list_html(db_name, systems, phenomena)

    return f"""
<section class="panel">
<div class="panel-header">
  <h2>Sector Map</h2>
  <span class="hint">Drag to rotate &middot; scroll to zoom &middot; dot size &asymp; star radius &middot; color &asymp; spectral type &amp; brightness &middot; translucent clouds &asymp; nebulae/asteroid fields, glowing points &asymp; black holes/neutron stars, near this sector</span>
</div>
<div class="starmap-layout">
<div class="starmap-viewport">
<canvas id="starmap-canvas" class="starmap-canvas" tabindex="0" role="application"
     aria-label="Interactive 3D sector map. Drag or use arrow keys to rotate, scroll or the zoom buttons to zoom."></canvas>
{scale_bar_html}
{noscript_html}
</div>
<div class="starmap-side">
<div class="starmap-controls" id="starmap-controls">
  <button type="button" class="starmap-btn" data-action="zoom-out" aria-label="Zoom out">&minus;</button>
  <button type="button" class="starmap-btn" data-action="zoom-in" aria-label="Zoom in">+</button>
  <button type="button" class="starmap-btn" data-action="reset">Reset view</button>
</div>
{info_panel}
</div>
</div>
<script type="application/json" id="starmap-data">{_json_script(scene_data)}</script>
</section>
"""
