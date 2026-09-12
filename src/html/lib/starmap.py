# html/lib/starmap.py

"""
Interactive 3D sector starmap: every placed star system in a sector
rendered as one plain `<div>` per star (two, overlapping, for a binary),
positioned in a real CSS 3D scene (`transform-style: preserve-3d`) built
from the system's (x, y, z) position within the sector -- sized by that
star's physical radius, and colored by its spectral color (from
`star_type`, e.g. "White" for an A-class star, "Blue" for an O-class star
-- see `physical_constants.SPECTRAL_CLASS_COLORS`), shaded by its
luminosity (brighter = more vivid/lighter, dimmer = more muted/darker) and
nudged by where its exact temperature falls within its spectral class's
range.

The scene's own outline is drawn around those dots too, and is one of two
shapes depending on whether the sector has a galaxy placement: a sector
with `shell_index`/`shell_slot_index` set gets a 12-edge wireframe of its
real, approximate on-shell wedge (see `sector_wedge_vertices_pc` --
bounded by the shell's own radial thickness and roughly how much angular
"real estate" this slot's Fibonacci placement owns among its neighbors,
so the shape reflects where the sector actually sits and points on its
shell, not a generic box); one without a placement falls back to the
sector's plain axis-aligned `edge_mpc` cube (6 bordered `.cube-face`
`<div>`s), same as before this wedge shape existed.

Unlike a hand-rolled JS rotation-matrix/projection routine, this hands
the actual 3D math to the browser: each star's (x, y, z) -- rotated once,
server-side, from its own sector-local axes into the galaxy frame when
the sector has a galaxy placement (`_rotate_to_galaxy_frame`, so it
agrees with the wedge outline/compass arrow below, both already
galaxy-frame quantities) -- is placed via plain layout position
(`left`/`top`) plus `transform: translateZ()` for depth (the wedge
wireframe's edges instead use a single
combined `translate3d()` + two rotations each, since an edge's endpoints
are two arbitrary points rather than one point plus a flat XY-plane
circle -- see `_wedge_wireframe_html`) -- rotating the whole
`.starmap-scene` element (`static/sectormap.js`, via drag) and letting
`preserve-3d` composite every descendant (star dots and the outline's own
`<div>`s alike) in true 3D, occlusion included, is what the browser's own
compositor is already built to do. Zoom is a separate, plain 2D `scale()`
on an *outer* wrapper (`.starmap-zoom`, kept outside the `perspective`
element rather than sandwiched between it and the rotating scene, so it
never disturbs the perspective math) -- since this is vector/DOM content
rather than a raster image, scaling it is lossless.

Clicking a dot doesn't navigate straight to `system.py` -- it populates
the info side panel via the `data-*` attributes read from the clicked
element (see `static/sectormap.js`), so a click shows details first and
the panel's own link is what navigates away.
"""

import colorsys
import math

from fmt import esc

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
    from stellarObjects.utils import milliparsecs_to_ly, mpc_to_pc, pc_to_mpc
except ImportError:
    # Same deployment gap as above -- without these, a sector with no
    # galaxy placement (or one this deployment can't reach the geometry
    # module for) just keeps the plain axis-aligned cube outline and no
    # scale bar; see `_wedge_edges_px`/`_scale_bar_attrs` below. Without
    # `cube_orientation` specifically, star dots fall back to being
    # plotted in their own unrotated local axes (see `_rotate_to_galaxy_frame`).
    sector_position_pc = None
    sector_wedge_vertices_pc = None
    cube_orientation = None
    mpc_to_pc = None
    pc_to_mpc = None
    milliparsecs_to_ly = None

# The 3D scene's on-screen footprint, in pixels -- a fixed size (unlike
# the old responsive SVG viewBox) since a CSS 3D scene needs an explicit
# width/height for its children's `left`/`top`/`translateZ` coordinates
# to mean anything; `sectormap.js`'s zoom control is what makes this not
# a hard ceiling for the viewer.
_SCENE_SIZE_PX = 320
_SCENE_HALF_PX = _SCENE_SIZE_PX / 2

_SUN_RADIUS_KM = SOLAR_RADIUS_M / 1000.0
_MIN_DOT_R = 3.0
_MAX_DOT_R = 14.0

# Secondary-star offset (binary systems), as a fraction of the primary's
# own dot radius -- "down and to the right", overlapping the primary
# rather than sitting fully clear of it. This is a fixed delta in the
# scene's own local coordinate space (not screen space), so the pair
# rotates rigidly together as one unit under `.starmap-scene`'s rotation
# instead of needing to be recomputed as the view turns.
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

# The standard "CSS 3D cube" recipe: 6 identically-sized, border-only
# faces, each pushed out along its own now-rotated local Z axis by half
# the cube's edge length. Every face is axis-aligned by construction (a
# plain, unskewed cube), so this needs no per-edge alignment math the way
# an arbitrary wireframe would -- just these 6 fixed transforms.
_CUBE_FACE_TRANSFORMS = {
    "front": f"translateZ({_SCENE_HALF_PX}px)",
    "back": f"rotateY(180deg) translateZ({_SCENE_HALF_PX}px)",
    "right": f"rotateY(90deg) translateZ({_SCENE_HALF_PX}px)",
    "left": f"rotateY(-90deg) translateZ({_SCENE_HALF_PX}px)",
    "top": f"rotateX(90deg) translateZ({_SCENE_HALF_PX}px)",
    "bottom": f"rotateX(-90deg) translateZ({_SCENE_HALF_PX}px)",
}


def _cube_faces_html():
    return "".join(
        f'<div class="cube-face" style="transform:{transform}"></div>'
        for transform in _CUBE_FACE_TRANSFORMS.values()
    )


# Which pairs of `sector_wedge_vertices_pc` vertex indices form the wedge's
# 12 edges -- indices are `4*r_bit + 2*phi_bit + theta_bit` (see that
# function's docstring), so two vertices share an edge exactly when their
# indices differ in a single bit (this is just a cube's own edge topology,
# reused for the wedge's 8 corners regardless of how curved its faces
# really are).
_WEDGE_EDGE_PAIRS = (
    (0, 1), (0, 2), (0, 4),
    (1, 3), (1, 5),
    (2, 3), (2, 6),
    (3, 7),
    (4, 5), (4, 6),
    (5, 7),
    (6, 7),
)


def _line_html(p1_px, p2_px, css_class):
    """
    Draws one straight `<div>` line between two arbitrary 3D scene-space
    pixel points -- shared by the wedge wireframe's 12 edges (below) and
    the "toward galactic center" compass arrow (`_compass_html`), since
    both are just "connect point A to point B" with no flatness
    assumption needed (unlike the plain cube's 6 filled `.cube-face`
    rectangles, which only work because a cube's faces genuinely are
    flat, axis-aligned rectangles).

    The line is one `<div>` stretched to its own length and pointed from
    `p1_px` to `p2_px` via two rotations (`rotateY` then `rotateZ`, the
    standard two-angle solution for aiming a local +X-axis segment at an
    arbitrary 3D point -- `rotateZ` picks the elevation off the XZ-plane,
    `rotateY` then picks the azimuth within it), the same "point A at B
    with plain trigonometry" idea as `_dot_html`'s billboarding, just for
    a line instead of a circle.
    """
    x1, y1, z1 = p1_px
    x2, y2, z2 = p2_px
    dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
    length = math.sqrt(dx * dx + dy * dy + dz * dz)
    if length < 1e-6:
        return ""
    yaw = -math.degrees(math.atan2(dz, dx))
    elevation = math.degrees(math.atan2(dy, math.hypot(dx, dz)))
    return (
        f'<div class="{css_class}" style="'
        f'width:{length:.2f}px; '
        "transform:"
        f"translate3d({_SCENE_HALF_PX + x1:.1f}px, {_SCENE_HALF_PX + y1:.1f}px, {z1:.1f}px) "
        f"rotateY({yaw:.3f}deg) rotateZ({elevation:.3f}deg);"
        '"></div>'
    )


def _wedge_wireframe_html(vertices_px):
    """
    Draws the sector's approximate on-shell cell as a 12-edge wireframe
    connecting `vertices_px` (8 `(x, y, z)` scene-space pixel points, see
    `_wedge_edges_px`) -- unlike the plain cube's 6 filled `.cube-face`
    rectangles, most of this shape's faces aren't flat (a wedge bounded by
    two shell radii and an angular patch has faces that are pieces of a
    sphere or a cone, not a plane -- see `sector_wedge_vertices_pc`), so
    there's no flat quad `transform` that would draw them correctly.
    Edges, in contrast, are just straight lines between two exact points
    and need no flatness assumption at all, so the shape is rendered as a
    wireframe instead of a solid.
    """
    return "".join(
        _line_html(vertices_px[i], vertices_px[j], "wedge-edge")
        for i, j in _WEDGE_EDGE_PAIRS
    )


_COMPASS_LABEL_W = 130.0
_COMPASS_LABEL_H = 20.0
# How far the compass arrow reaches past the scene's own half-size --
# past 1.0 so its tip clears a full-size cube/wedge instead of ending
# right at (or inside) its own boundary, but not by much more than that:
# the viewport clips anything that lands outside it (`overflow: hidden`
# on `.starmap-viewport`), and dragging can point this arrow in any
# direction, so a larger reach makes the label likelier to get clipped
# at the default view angle.
_COMPASS_ARROW_REACH = 1.15


def _compass_html(center_pc):
    """
    Draws an arrow from the sector's own local origin toward the galactic
    center, plus a billboarded "Galactic Center" label at its tip -- the
    sector-map equivalent of a map's north arrow, except the direction it
    points is computed exactly from this sector's own stored
    `sectors.center_x/y/z_pc` (the negative of the sector's own outward
    radial direction, `-normalize(center_pc)`) rather than fixed to a
    constant screen direction.

    This arrow and `_wedge_edges_px`'s wedge outline are both computed
    directly from galaxy-frame quantities (`sectors.center_x/y/z_pc`,
    `sector_wedge_vertices_pc`), so they were always correct on their own
    terms and need no rotation of their own. What used to be wrong is that
    star dots (`render_map_panel`, `star_systems.position_x/y/z_mpc`) were
    plotted as if their own local (x, y, z) axes already ran parallel to
    the galaxy frame's -- not a design convention this project actually
    enforces at generation time (`galaxyGen.py` never rotates a sector's
    local star positions to align with its galaxy-frame placement).
    `render_map_panel` now rotates star positions into the galaxy frame
    at render time instead (`_rotate_to_galaxy_frame`, the "Cube
    orientation" convention docs/design/galaxy-coordinate-system.md's
    section 3 already proposes and `sectorGeometry.py` already applies to
    this sector's own wedge vertices), so this arrow, the wedge outline,
    and the star dots it surrounds all agree on one frame.

    Args:
        center_pc (tuple or None): `(center_x_pc, center_y_pc,
                                    center_z_pc)`, or `None` if this
                                    sector has no galaxy placement.

    Returns:
        str: The arrow's and label's HTML, or `""` if `center_pc` is
             `None` or (within floating-point tolerance) the galactic
             center itself, which has no meaningful direction to point.
    """
    if center_pc is None or any(c is None for c in center_pc):
        return ""

    cx, cy, cz = center_pc
    norm = math.sqrt(cx * cx + cy * cy + cz * cz)
    if norm < 1e-9:
        return ""

    ux, uy, uz = -cx / norm, -cy / norm, -cz / norm
    reach = _SCENE_HALF_PX * _COMPASS_ARROW_REACH
    tip_px = (ux * reach, -uy * reach, uz * reach)

    origin_px = (0.0, 0.0, 0.0)
    arrow = _line_html(origin_px, tip_px, "compass-arrow")
    if not arrow:
        return ""

    tip_x, tip_y, tip_z = tip_px
    left = _SCENE_HALF_PX + tip_x - _COMPASS_LABEL_W / 2
    top = _SCENE_HALF_PX + tip_y - _COMPASS_LABEL_H / 2
    label = (
        '<div class="compass-label-anchor" '
        f'style="left:{left:.1f}px; top:{top:.1f}px; '
        f'width:{_COMPASS_LABEL_W:.0f}px; height:{_COMPASS_LABEL_H:.0f}px; '
        f'transform:translateZ({tip_z:.1f}px);">'
        '<div class="compass-label billboard">Galactic Center &rarr;</div>'
        "</div>"
    )
    return arrow + label


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
    (`_wedge_edges_px`/`_compass_html`), which read the sector's *actual*
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
    `sector_wedge_vertices_pc`) as 8 `(x, y, z)` scene-space pixel points,
    in the same coordinate convention `render_map_panel` uses for star
    dots (sector-center-relative parsecs -> milliparsecs -> normalized by
    `half_edge` -> scene pixels, with the y-axis flipped since +y is "up"
    on screen but down in CSS layout) -- so the wireframe and the star
    dots it surrounds always share one consistent frame.

    Unlike a star dot's normalized position, these are deliberately *not*
    clamped to the +-1.05 the cube fallback uses: the whole point of this
    shape is that it doesn't stay inside the sector's own edge_mpc cube,
    since the shell's angular patch and the cube's flat sides don't
    coincide (see `sector_wedge_vertices_pc`'s docstring).

    Returns:
        list[tuple] or None: 8 `(x_px, y_px, z_px)` points, or `None` if
                             this sector has no galaxy placement (no
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


def _ly_per_px_at_zoom_1(half_edge):
    """
    Light-years per on-screen pixel at zoom factor 1 -- an exact ratio
    derived straight from `half_edge` (half of `edge_mpc`, the sector's
    real, stored size) mapping to `_SCENE_HALF_PX` pixels, the same
    normalizing divisor every star dot and wedge vertex on this map is
    already placed by. `sectormap.js`'s scale-bar legend divides this by
    the live zoom factor and picks a round bar length from it, so the
    bar always reflects the sector's actual physical scale rather than
    an arbitrary fixed pixel-per-lightyear guess.

    Returns:
        float or None: ly per pixel, or `None` if `milliparsecs_to_ly`
                       isn't importable in this deployment (the scale bar
                       is then omitted entirely rather than shown in raw
                       milliparsecs).
    """
    if milliparsecs_to_ly is None:
        return None
    mpc_per_px = half_edge / _SCENE_HALF_PX
    return milliparsecs_to_ly(mpc_per_px)


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
    """Maps a star's physical radius to a dot radius in pixels -- square-
    root scaled against the Sun's radius (linear scaling would make red
    dwarfs invisible next to giants, which differ by 2+ orders of
    magnitude in radius_km) and clamped so the map stays legible at
    either extreme."""
    if not radius_km or radius_km <= 0:
        return _MIN_DOT_R
    ratio = radius_km / _SUN_RADIUS_KM
    dot_r = _MIN_DOT_R + 4.0 * math.sqrt(ratio)
    return max(_MIN_DOT_R, min(_MAX_DOT_R, dot_r))


def _dot_html(db_name, system, star, x_px, y_px, z_px, label_suffix, max_r=None):
    """
    Builds one star as two nested `<div>`s -- an outer "anchor" that just
    positions it (plain layout `left`/`top`, recentered on the scene's
    middle, for x/y; `transform: translateZ()` for z) and an inner
    `.star-dot` that draws the actual circle. The split exists for
    billboarding: a lone dot positioned this way is a flat disc lying in
    the scene's local plane, so once the scene rotates far enough it's
    seen edge-on and all but disappears -- exactly wrong for a map whose
    entire point is looking at these from any angle. `sectormap.js`
    counter-rotates the *inner* div by the scene's current rotation on
    every drag frame -- `rotateY(-rotateY) rotateX(-rotateX)`, the plain
    algebraic inverse (reverse function order, negated angles) of
    `.starmap-scene`'s own `rotateX(rotateX) rotateY(rotateY)` -- so the
    circle always faces the camera; the outer anchor is what still
    carries the correct, unrotated (x, y, z) position through the scene's
    rotation. This inverse is only correct because `.starmap-stage` has no
    `perspective`: a vanishing-point projection would make the real
    composition genuinely projective rather than pure rotation, and a
    plain inverse stops cancelling it correctly (confirmed directly --
    with `perspective` still set, this same formula only worked at the
    one rotation angle it happened to be tested at, and broke, in a
    different way, at another).
    """
    dot_r = _star_dot_radius(star["radius_km"])
    if max_r is not None:
        dot_r = min(dot_r, max_r)
    fill, stroke = _star_color(star["star_type"], star["temperature_k"], star["luminosity_w"])
    name = f'{system["name"]}{label_suffix}'
    left = _SCENE_HALF_PX + x_px - dot_r
    top = _SCENE_HALF_PX + y_px - dot_r
    return (
        '<div class="star-dot-anchor" '
        f'style="left:{left:.1f}px; top:{top:.1f}px; width:{dot_r * 2:.1f}px; height:{dot_r * 2:.1f}px; '
        f'transform:translateZ({z_px:.1f}px);">'
        '<div class="star-dot billboard" tabindex="0" role="button" '
        f'style="background:{fill}; border-color:{stroke};" '
        f'data-name="{esc(name)}" '
        f'data-type="{esc(star["star_type"])}" '
        f'data-temp="{esc(star["temp_display"])}" '
        f'data-quadrant="{esc(system["quadrant"])}" '
        f'data-location="{esc(system["location"])}" '
        f'data-href="system.py?db={esc(db_name)}&amp;id={system["id"]}" '
        f'aria-label="{esc(name)}" title="{esc(name)}"></div>'
        '</div>'
    )


def render_map_panel(db_name, edge_mpc, shell_index, shell_slot_index, center_pc, systems):
    """
    Builds the "Sector Map" panel: a draggable/zoomable 3D scene (see
    `static/sectormap.js` for the rotate/zoom/click wiring) with one dot
    per placed star system (two, overlapping, for a binary -- the primary
    at the system's actual position, the secondary offset down-and-right
    from it), plus an info side panel that the same script fills in when
    a dot is clicked.

    The scene's own outline is the sector's approximate on-shell wedge
    (see `_wedge_edges_px`/`sector_wedge_vertices_pc`) when this sector
    has a galaxy placement (`shell_index`/`shell_slot_index` both set) and
    the geometry helpers are importable -- a 12-edge wireframe reflecting
    where the sector actually sits on its shell, not a generic cube. Any
    sector without a placement (standalone `sectorGen.py` tooling, or one
    migrated from a pre-v4 database) falls back to the plain axis-aligned
    `.cube-face` cube instead, same as before this shape existed.

    Args:
        db_name (str): The current `?db=` value, used to build each dot's
                       `data-href` (`system.py?db=...&id=...`).
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
                                   Center" compass arrow (`_compass_html`)
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

    Returns:
        str: A complete `<section class="panel">` block.
    """
    half_edge = (edge_mpc / 2) if edge_mpc else 1.0

    dots = []
    for system in systems:
        # +y is "up" on screen; CSS's own y axis increases downward, so
        # the sign flips here once, at the one place normalized position
        # becomes a pixel coordinate -- everything downstream (including
        # the binary offset below) works in already-screen-oriented
        # pixels. No depth sort is needed here (unlike the old fixed
        # isometric SVG) -- `preserve-3d` composites every dot and cube
        # face by its real depth as the scene rotates, live, in the
        # browser itself.
        galaxy_x, galaxy_y, galaxy_z = _rotate_to_galaxy_frame(
            center_pc, (system["x"] or 0, system["y"] or 0, system["z"] or 0)
        )
        nx = max(-1.05, min(1.05, galaxy_x / half_edge))
        ny = max(-1.05, min(1.05, galaxy_y / half_edge))
        nz = max(-1.05, min(1.05, galaxy_z / half_edge))
        x_px = nx * _SCENE_HALF_PX
        y_px = -ny * _SCENE_HALF_PX
        z_px = nz * _SCENE_HALF_PX

        stars = system["stars"]
        is_binary = len(stars) > 1
        primary_suffix = " -- Primary" if is_binary else ""
        dots.append(_dot_html(db_name, system, stars[0], x_px, y_px, z_px, primary_suffix))

        if is_binary:
            primary_r = _star_dot_radius(stars[0]["radius_km"])
            offset = primary_r * _BINARY_OFFSET_FRACTION
            dots.append(_dot_html(
                db_name, system, stars[1], x_px + offset, y_px + offset, z_px, " -- Secondary",
                max_r=primary_r * _SECONDARY_MAX_RATIO,
            ))

    wedge_vertices = _wedge_edges_px(shell_index, shell_slot_index, edge_mpc, half_edge)
    if wedge_vertices is not None:
        outline_html = _wedge_wireframe_html(wedge_vertices)
        shape_hint = "outline &asymp; sector's real position/orientation on its shell"
    else:
        outline_html = _cube_faces_html()
        shape_hint = "cube edge &asymp; sector size (not placed in a galaxy)"

    compass_html = _compass_html(center_pc)

    # No role/aria-label here -- `role="img"` on an ancestor would flatten
    # every descendant (each star dot's own `role="button"`/`tabindex`)
    # into a single opaque image for assistive tech, breaking keyboard
    # access to the dots. The accessible description lives on
    # `.starmap-stage` below instead, one level up.
    scene = (
        '<div class="starmap-scene" id="starmap-scene" '
        f'style="width:{_SCENE_SIZE_PX}px; height:{_SCENE_SIZE_PX}px;">'
        f"{outline_html}{compass_html}{''.join(dots)}</div>"
    )

    if systems:
        info_panel = (
            '<aside class="starmap-info" id="starmap-info">'
            '<p class="hint">Click a star system for details.</p></aside>'
        )
    else:
        info_panel = '<aside class="starmap-info" id="starmap-info"><p class="hint">No systems placed in this sector.</p></aside>'

    # The scale bar's own live pixel width is computed and kept up to date
    # by sectormap.js (it changes with zoom) -- this only hands over the
    # one fixed, exact ratio (ly per pixel *at zoom 1*) it needs to do
    # that; see `_ly_per_px_at_zoom_1`. Omitted (and the bar left blank)
    # when that ratio isn't computable in this deployment.
    ly_per_px = _ly_per_px_at_zoom_1(half_edge)
    scale_attr = f' data-ly-per-px="{ly_per_px:.10g}"' if ly_per_px else ""
    scale_bar_html = (
        f'<div class="starmap-scale" id="starmap-scale"{scale_attr}>'
        '<span class="starmap-scale-bar" id="starmap-scale-bar"></span>'
        '<span class="starmap-scale-label" id="starmap-scale-label"></span>'
        "</div>"
        if ly_per_px else ""
    )

    return f"""
<section class="panel">
<div class="panel-header">
  <h2>Sector Map</h2>
  <span class="hint">Drag to rotate &middot; scroll to zoom &middot; dot size &asymp; star radius &middot; color &asymp; spectral type &amp; brightness &middot; {shape_hint}</span>
</div>
<div class="starmap-layout">
<div class="starmap-viewport">
<div class="starmap-zoom" id="starmap-zoom">
<div class="starmap-stage" id="starmap-stage" tabindex="0" role="application"
     aria-label="Interactive 3D sector map. Drag or use arrow keys to rotate, scroll or the zoom buttons to zoom.">
{scene}
</div>
</div>
{scale_bar_html}
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
</section>
"""
