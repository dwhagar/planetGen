# html/lib/starmap.py

"""
Interactive 3D sector-cube starmap: every placed star system in a sector
rendered as one plain `<div>` per star (two, overlapping, for a binary),
positioned in a real CSS 3D scene (`perspective` + `transform-style:
preserve-3d`) built from the system's (x, y, z) position within the
sector's cube -- sized by that star's physical radius, and colored by its
spectral color (from `star_type`, e.g. "White" for an A-class star,
"Blue" for an O-class star -- see `physical_constants.SPECTRAL_CLASS_COLORS`),
shaded by its luminosity (brighter = more vivid/lighter, dimmer = more
muted/darker) and nudged by where its exact temperature falls within its
spectral class's range.

Unlike a hand-rolled JS rotation-matrix/projection routine, this hands
the actual 3D math to the browser: each star's (x, y, z) is placed once,
as-is, via plain layout position (`left`/`top`) plus `transform:
translateZ()` for depth -- rotating the whole `.starmap-scene` element
(`static/sectormap.js`, via drag) and letting `preserve-3d` composite
every descendant (star dots and the 6 cube-face `<div>`s alike) in true
3D, occlusion included, is what the browser's own compositor is already
built to do. Zoom is a separate, plain 2D `scale()` on an *outer*
wrapper (`.starmap-zoom`, kept outside the `perspective` element rather
than sandwiched between it and the rotating scene, so it never disturbs
the perspective math) -- since this is vector/DOM content rather than a
raster image, scaling it is lossless.

Clicking a dot doesn't navigate straight to `system.py` -- it populates
the info side panel via the `data-*` attributes read from the clicked
element (see `static/sectormap.js`), so a click shows details first and
the panel's own link is what navigates away.
"""

import colorsys
import math

from dbutil import esc

try:
    from stellarObjects.physical_constants import SPECTRAL_CLASS_COLORS, TEMP_RANGES, SOLAR_LUMINOSITY, SOLAR_RADIUS_M
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # duplicated literally (matching dbutil.py's own fallback constants)
    # rather than left unimportable, since these drive the map's color
    # math directly.
    SPECTRAL_CLASS_COLORS = {'O': 'Blue', 'B': 'Blue-White', 'A': 'White', 'F': 'Yellow-White', 'G': 'Yellow', 'K': 'Orange', 'M': 'Red'}
    TEMP_RANGES = {
        'O': (30000, 60000), 'B': (10000, 30000), 'A': (7500, 10000), 'F': (6000, 7500),
        'G': (5200, 6000), 'K': (3700, 5200), 'M': (2400, 3700),
    }
    SOLAR_LUMINOSITY = 3.82e26
    SOLAR_RADIUS_M = 6.957e8

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
        '<div class="star-dot" tabindex="0" role="button" '
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


def render_map_panel(db_name, edge_mpc, systems):
    """
    Builds the "Sector Map" panel: a draggable/zoomable 3D cube (see
    `static/sectormap.js` for the rotate/zoom/click wiring) with one dot
    per placed star system (two, overlapping, for a binary -- the primary
    at the system's actual position, the secondary offset down-and-right
    from it), plus an info side panel that the same script fills in when
    a dot is clicked.

    Args:
        db_name (str): The current `?db=` value, used to build each dot's
                       `data-href` (`system.py?db=...&id=...`).
        edge_mpc (float): The sector's cube edge (`sectors.edge_mpc`) --
                          every system's `position_*_mpc` is relative to
                          the sector's cubic center (see schema.sql's
                          `star_systems` comment), so half of this is the
                          normalizing divisor for each axis.
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
        nx = max(-1.05, min(1.05, (system["x"] or 0) / half_edge))
        ny = max(-1.05, min(1.05, (system["y"] or 0) / half_edge))
        nz = max(-1.05, min(1.05, (system["z"] or 0) / half_edge))
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

    # No role/aria-label here -- `role="img"` on an ancestor would flatten
    # every descendant (each star dot's own `role="button"`/`tabindex`)
    # into a single opaque image for assistive tech, breaking keyboard
    # access to the dots. The accessible description lives on
    # `.starmap-stage` below instead, one level up.
    scene = (
        '<div class="starmap-scene" id="starmap-scene" '
        f'style="width:{_SCENE_SIZE_PX}px; height:{_SCENE_SIZE_PX}px;">'
        f"{_cube_faces_html()}{''.join(dots)}</div>"
    )

    if systems:
        info_panel = (
            '<aside class="starmap-info" id="starmap-info">'
            '<p class="hint">Click a star system for details.</p></aside>'
        )
    else:
        info_panel = '<aside class="starmap-info" id="starmap-info"><p class="hint">No systems placed in this sector.</p></aside>'

    return f"""
<section class="panel">
<div class="panel-header">
  <h2>Sector Map</h2>
  <span class="hint">Drag to rotate &middot; scroll to zoom &middot; dot size &asymp; star radius &middot; color &asymp; spectral type &amp; brightness</span>
</div>
<div class="starmap-layout">
<div class="starmap-zoom" id="starmap-zoom">
<div class="starmap-stage" id="starmap-stage" tabindex="0" role="application"
     aria-label="Interactive 3D sector map. Drag or use arrow keys to rotate, scroll or the zoom buttons to zoom.">
{scene}
</div>
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
