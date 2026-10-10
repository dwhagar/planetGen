# planetgen/web/maps/systemmap.py

"""
Interactive 2D "System Map" for `system.py`: a true top-down plot of every
body's *real* position, not a schematic. Each star, planet, moon, and
asteroid belt is placed at its actual angle (from `position_x/y_km`,
already stored relative to whatever it orbits -- see
`planetPhysics.update_orbital_position`'s docstring) and a distance from
its own anchor that's log-scaled into a shared pixel budget the whole
scene uses (`_radial_scale_bounds`/`_radial_px`) -- a log scale is what
lets one diagram span both a close binary's ~0.05-0.3 AU separation and an
outer planet's tens-of-AU orbit without either collapsing to a point or
blowing out the frame, the same reasoning `docs/design/
galaxy-coordinate-system.md` gives for using parsecs (not km/ly) one scale
level up. A once-only pass of `_relax_markers` then nudges apart any two
markers real angle happened to place too close together (the exact,
minimal-displacement fix, not a wholesale relayout) -- true position
first, decluttering only where it's actually needed.

A single-star or 'close' (P-type) binary system's planets all orbit a
shared origin (the barycenter -- see `doubleStar.BinaryStarProxy`'s
docstring: a 'close' pair's own two stars are placed around that same
origin too, split by their real mass ratio from
`star_systems.binary_mutual_position_x/y/z_km`, the secondary's position
relative to the primary). A 'wide' (S-type) pair's two stars each anchor
their own independent set of planets/belts (`planets.star_id`/
`asteroid_belts.star_id` says which) -- unlike the old schematic map, this
one draws both stars' full systems in the same true-position scene (no
need for the old "primary only, see the tables below for the secondary"
compromise: once real absolute positions exist, a wide pair's second star
and its planets are just more points in the same frame).

Each body is colored by planet class and sized (log-scaled, so a Class D
moonlet and a Class J gas giant both stay visible on the same diagram) by
`radius_km` -- unrelated to the *radial position* log scale above, which
only ever concerns where a body sits, never how big its own marker is.

Unlike `starmap.py`'s draggable 3D cube (positions come straight from
stored (x, y, z) coordinates -- there's real 3D data to project), this
map only ever needs a flat top-down projection (x, y; z is dropped, the
same way a real solar-system diagram ignores the small inclination every
planet's orbit actually has) -- so instead of a rotatable scene, this
renders one flat `<svg>` per "scene": the whole system, plus one more for
every planet that has moons (moons orbit their own parent planet, a
completely different distance scale, so they get their own scene and
their own local log-radial fit rather than sharing the system scene's).
`static/systemmap.js` just toggles which scene is visible. Clicking a
planet with moons swaps the view so that planet sits at the scene's own
origin with its own moons arranged around it (mirroring a real "zoom into
this planet's moon system" diagram); clicking any other object or a belt
instead fills the info side panel from its `data-*` attributes, exactly
like `starmap.py`/`sectorscene.js`'s own click-for-info pattern.

Every star/planet/moon marker's flat SVG circle is also live-rendered as a
small rotating shaded sphere (three.js, same vendored build
`sectorscene.js` uses -- `static/systemmap.js`'s `#sysmap-spheres-canvas`),
colored by the body's own class (`class_color`, handed over pre-resolved
as `data-color` so this module stays the one place that mapping lives --
a star instead gets its own real spectral-type color from `star_color`,
same as `starmap.py`), banded with a tilted ring for a gas giant
(`data-bodytype`), and wrapped in a soft fresnel-glow atmosphere shell,
tinted by `data-surfacetemp`, when `data-hasatmosphere` is set. One shared
WebGL context draws every visible marker's own sphere each frame via
scissored sub-viewports (not one `<canvas>`/context per body -- browsers
cap how many WebGL contexts can exist at once, easily blown through by a
crowded system), positioned and sized to exactly cover that marker's own
`<circle>` (whose own fill becomes transparent once its sphere is live,
via `sysmap-sphere-active` -- see `static/systemmap.js`), so the sphere
reads as *replacing* the flat marker rather than a separate preview
floating beside it. This is an appearance layer only (no real position
data goes into it) -- the true-position diagram (marker placement, labels,
click-for-info, scene switching) remains this map's actual subject and is
untouched by whether a sphere successfully renders on top of it.
"""

import colorsys
import math
import statistics

from planetgen.web.lib.fmt import esc, format_distance_km, format_pressure_pa, format_speed_kms, format_temperature_k
from planetgen.web.maps.starmap import star_color, _SUN_RADIUS_KM
from planetgen.web.lib.tabledisplay import (
    format_body_distance, format_body_mass, format_body_radius, format_period, format_star_luminosity, format_star_mass, format_star_radius,
    to_plain_text,
)

try:
    from planetgen.tuning import PLANET_CLASSES
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # the map still works, just without each class's one-line flavor
    # description in the info panel.
    PLANET_CLASSES = {}

# The whole diagram's fixed coordinate space -- square (unlike the old
# short/wide letterboxed strip), since a true top-down position plot needs
# room in every direction, not just "away from the star". Never scrolled
# or zoomed, so every system's worth of bodies has to fit inside it
# regardless of how many there are (`_radial_px`'s shared log scale is
# what keeps that bounded).
_VIEW_SIZE_PX = 700.0
_CENTER_PX = _VIEW_SIZE_PX / 2

# A body's *ideal* radial distance from its own anchor is log-scaled
# between these two bounds (see `_radial_scale_bounds`/`_radial_px`) --
# `_MIN_RADIUS_PX` leaves room near the center for the star marker(s)
# themselves, `_RADIUS_SPREAD_PX` is how much room is left after that
# before running into the fixed frame's own edge (with a small margin for
# labels).
_MIN_RADIUS_PX = 55.0
_RADIUS_SPREAD_PX = 280.0

# Minimum clearance (beyond simple non-overlap) enforced between any two
# markers by `_relax_markers` -- stars and planets/moons alike share one
# list, so a giant/supergiant star's own marker gets exactly the same
# "push the other one out of the way" treatment as two crowded planets,
# rather than the separate, star-specific clearance hack the old
# fixed-layout map needed.
_MARKER_GAP_PX = 6.0
_RELAX_ITERATIONS = 30

# Planet/moon marker radius is log-scaled against these two bounds --
# Class D's own floor and Class J's own ceiling (see
# `program_constants.PLANET_CLASSES`) -- rather than each system's own
# min/max, so a given radius_km always maps to the same dot size across
# every system's map, not just relative to its own system-mates.
_BODY_RADIUS_LO_KM = 50.0
_BODY_RADIUS_HI_KM = 250_000.0
_PLANET_MIN_R = 6.0
_PLANET_MAX_R = 24.0
_MOON_MIN_R = 4.5
_MOON_MAX_R = 19.0

_STAR_MIN_R = 22.0
_STAR_MAX_R = 48.0

# A drilled-into planet acting as its own moon-system's "star" is drawn at
# one fixed size regardless of its actual radius_km -- it's the diagram's
# new focal point, not a body being measured against its own moons on the
# same log scale.
_CENTER_PLANET_R = 36.0

_ZONE_LABELS = {"h": "Hot zone", "e": "Ecosphere", "c": "Cold zone"}

# An asteroid belt is a *range* of orbit radii -- drawn as a full ring
# around its own anchor, the natural true-position shape for "every angle
# at roughly this distance", rather than the old schematic map's
# directional shaded band. The ring runs from the belt's inner edge to its
# outer edge, both mapped through the same log scale as every orbit
# (`_belt_band`), so a planet orbit outside the belt is drawn outside the
# ring. A belt too thin to see is widened to `_BELT_MIN_BAND_PX`, but never
# past a neighboring orbit (less `_BELT_ORBIT_GAP_PX`).
_BELT_MIN_BAND_PX = 9.0
_BELT_ORBIT_GAP_PX = 2.0

# Anchor color per planet class, hand-picked (not derived from any
# physical model, unlike `starmap.py`'s spectral colors -- there's no
# single physical quantity a planet class maps to the way a star's
# spectral letter maps to a blackbody hue) to read as distinct from its
# neighbors while still gesturing at the class's own flavor: A/B/E/K's
# reds-and-rust for volcanic/molten worlds, C/D's muted greys for dead/icy
# small bodies, F/G/H/L's earthy greens-and-tans for the barren-to-
# vegetated progression, M/O/P's blues for the habitable "water world"
# spread (M dark, O a brighter ocean blue, P pale glacier blue), I/J/T's
# gas-giant palette, N/Q/V's outliers (Venus-hot mustard, an eccentric-
# orbit violet, a high-gravity super-Earth magenta), and S's dark slate for
# the barren rocky super-Earth.
_CLASS_COLORS = {
    "A": "#d9483f",
    "B": "#e8823a",
    "C": "#8a7a66",
    "D": "#c9d6e0",
    "E": "#8b3a2b",
    "F": "#5f8a72",
    "G": "#a89a5c",
    "H": "#d9b568",
    "I": "#7fd6d0",
    "J": "#d9a441",
    "K": "#b5651d",
    "L": "#4f9e5f",
    "M": "#1c3f7a",
    "N": "#cbb03a",
    "O": "#2f8fc9",
    "P": "#a9d8f0",
    "Q": "#7a5cc9",
    "S": "#5b6878",
    "T": "#7a86d9",
    "V": "#9c3f63",
}
_DEFAULT_CLASS_COLOR = "#8a8f9c"


def class_color(planet_class):
    return _CLASS_COLORS.get((planet_class or "").upper(), _DEFAULT_CLASS_COLOR)


def _class_description(planet_class):
    entry = PLANET_CLASSES.get((planet_class or "").upper())
    return entry["description"] if entry else ""


def _hex_to_rgb01(hex_color):
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i + 2], 16) / 255 for i in (0, 2, 4))


def _darken_hex(hex_color, amount):
    r, g, b = _hex_to_rgb01(hex_color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    l = max(0.0, l - amount)
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return "#{:02x}{:02x}{:02x}".format(round(r * 255), round(g * 255), round(b * 255))


def _text_color_for(hex_color):
    """Picks black or white text for a class letter drawn on top of
    `hex_color`, via the standard relative-luminance threshold -- so a
    pale class (D, P) gets dark text and a dark class (M) gets light text
    instead of one hardcoded color washing out on half the palette."""
    r, g, b = _hex_to_rgb01(hex_color)
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "#1b1e28" if luminance > 0.55 else "#ffffff"


def _log_scale(value_km, lo_km, hi_km, min_px, max_px):
    if not value_km or value_km <= 0:
        return min_px
    value_km = max(lo_km, min(hi_km, value_km))
    frac = (math.log10(value_km) - math.log10(lo_km)) / (math.log10(hi_km) - math.log10(lo_km))
    return min_px + (max_px - min_px) * frac


def _planet_radius_px(radius_km):
    return _log_scale(radius_km, _BODY_RADIUS_LO_KM, _BODY_RADIUS_HI_KM, _PLANET_MIN_R, _PLANET_MAX_R)


def _moon_radius_px(radius_km):
    return _log_scale(radius_km, _BODY_RADIUS_LO_KM, _BODY_RADIUS_HI_KM, _MOON_MIN_R, _MOON_MAX_R)


def _star_radius_px(radius_km):
    if not radius_km or radius_km <= 0:
        return _STAR_MIN_R
    ratio = radius_km / _SUN_RADIUS_KM
    dot_r = _STAR_MIN_R + 9.0 * math.sqrt(ratio)
    return max(_STAR_MIN_R, min(_STAR_MAX_R, dot_r))


def _radial_scale_bounds(distances_km):
    """
    Picks the `(lo, hi)` log bounds every body's own orbital distance is
    scaled between (`_radial_px`) -- shared across the *whole* scene
    (every star's own separation from the barycenter, every planet's own
    distance from its anchor, every belt's own nominal distance), which is
    exactly what makes this one continuous log scale rather than several
    incompatible ones: a close binary's ~0.05 AU separation and an outer
    planet's 30+ AU orbit both land somewhere sane on it instead of either
    collapsing together or blowing out the frame.

    Mirrors the old fixed-layout map's own dynamic-range derivation (a
    typical log-step below the smallest real distance as the floor, rather
    than anchoring the floor to the smallest distance itself, which would
    force it to `frac=0` by construction and collapse the star-to-first-
    orbit gap to nothing) -- just no longer tied to any particular ordering
    of bodies, since real angle (not rank) now decides direction.

    Args:
        distances_km (list[float]): Every real, meaningful distance in the
            scene (0 or negative entries are ignored -- a body sitting
            exactly on its own anchor has nothing to scale).

    Returns:
        tuple[float, float]: `(lo_km, hi_km)`.
    """
    valid = sorted(d for d in distances_km if d and d > 0)
    if len(valid) >= 2:
        log_vals = [math.log10(d) for d in valid]
        typical_step = statistics.median(b - a for a, b in zip(log_vals, log_vals[1:]))
        if typical_step <= 0:
            typical_step = 0.3
        lo = 10 ** (log_vals[0] - typical_step)
        hi = 10 ** log_vals[-1]
    elif valid:
        lo, hi = valid[0] / 3.0, valid[0]
    else:
        lo = hi = 1.0
    return lo, hi


def _radial_px(distance_km, lo_km, hi_km):
    """Maps `distance_km` to a pixel radius via the shared log scale
    (`_radial_scale_bounds`) -- `0.0` for a body sitting on its own anchor
    (nothing to place it away from), `_MIN_RADIUS_PX + _RADIUS_SPREAD_PX`
    (the scale's own far edge) when every real distance in the scene
    happens to be identical (`hi_km <= lo_km`, nothing to spread across)."""
    if not distance_km or distance_km <= 0:
        return 0.0
    if hi_km <= lo_km:
        return _MIN_RADIUS_PX + _RADIUS_SPREAD_PX
    distance_km = max(lo_km, min(hi_km, distance_km))
    frac = (math.log10(distance_km) - math.log10(lo_km)) / (math.log10(hi_km) - math.log10(lo_km))
    return _MIN_RADIUS_PX + frac * _RADIUS_SPREAD_PX


def _orbit_radius_km(planet):
    """A planet's distance from its anchor in 3D, which is its orbit's
    radius. Not the top-down `hypot(x, y)`: on an inclined orbit that is
    shorter, and drew a planet just outside a belt (whose edges are real
    radii) with its orbit inside the belt's ring."""
    return math.hypot(planet.get("position_x_km") or 0.0, planet.get("position_y_km") or 0.0,
                      planet.get("position_z_km") or 0.0) or planet.get("distance_km") or 0.0


def _polar_to_px(anchor_x_px, anchor_y_px, radius_px, x_km, y_km):
    """
    Places a body at `radius_px` from `(anchor_x_px, anchor_y_px)`, in its
    *real* direction (`atan2(y_km, x_km)`) -- the one place this module
    turns a real position into a screen point. `+y` is "up" in the
    generator's own convention; SVG's own `y` grows downward, so the sign
    flips here once, same as `starmap.py`'s identical convention one scale
    level up.
    """
    if radius_px <= 0 or (not x_km and not y_km):
        return anchor_x_px, anchor_y_px
    theta = math.atan2(y_km, x_km)
    return anchor_x_px + radius_px * math.cos(theta), anchor_y_px - radius_px * math.sin(theta)


def _relax_markers(markers, min_gap):
    """
    Nudges apart any two markers (each a dict with `cx`/`cy`/`r`, mutated
    in place) that ended up closer than `r_a + r_b + min_gap` after real
    angle/distance placement -- the standard pairwise-repulsion micro-
    adjustment (move each of a colliding pair half the overlap apart,
    along the line between them), repeated a fixed number of times so a
    tightly-packed run settles into a stable, minimally-displaced, non-
    overlapping layout instead of one collision cascading arbitrarily far.
    Real position is the primary layout signal; this only ever moves a
    marker as far as legibility actually requires, unlike the old fixed
    schematic map's rank-based spacing (which never had real angles to
    preserve in the first place).

    A marker with `fixed: True` (a star, already relaxed against its own
    companion in an earlier pass -- see `_render_system_scene`) still
    repels others but never moves itself: planets/moons are anchored on a
    star's own *final* drawn position, so letting a later pass keep
    nudging that star around here would silently pull its own orbit rings
    out from under it.

    Args:
        markers (list[dict]): Each with `cx`, `cy`, `r` (pixels), and
            optionally `fixed` (bool, default False) -- updated in place.
        min_gap (float): Minimum clearance, beyond simple non-overlap,
            between any two markers' own edges.
    """
    for _ in range(_RELAX_ITERATIONS):
        moved = False
        for i in range(len(markers)):
            for j in range(i + 1, len(markers)):
                a, b = markers[i], markers[j]
                dx = b["cx"] - a["cx"]
                dy = b["cy"] - a["cy"]
                dist = math.hypot(dx, dy)
                min_dist = a["r"] + b["r"] + min_gap
                if dist >= min_dist:
                    continue
                if dist < 1e-6:
                    dx, dy, dist = 1.0, 0.0, 1.0
                ux, uy = dx / dist, dy / dist
                a_fixed, b_fixed = a.get("fixed", False), b.get("fixed", False)
                if a_fixed and b_fixed:
                    continue  # both pinned -- nothing this pass can do
                moved = True
                overlap = min_dist - dist
                if a_fixed:
                    b["cx"] += ux * overlap
                    b["cy"] += uy * overlap
                elif b_fixed:
                    a["cx"] -= ux * overlap
                    a["cy"] -= uy * overlap
                else:
                    a["cx"] -= ux * overlap / 2.0
                    a["cy"] -= uy * overlap / 2.0
                    b["cx"] += ux * overlap / 2.0
                    b["cy"] += uy * overlap / 2.0
        if not moved:
            break


_LABEL_CHAR_WIDTH_PX = 8.5
_LABEL_MIN_HALFWIDTH_PX = 20.0
_LABEL_HALF_HEIGHT_PX = 9.0
_LABEL_GAP_PX = 4.0

# Tried in this order for every body -- "below"/"above" first (the
# original two, cheapest to read since they sit on the marker's own
# vertical axis), then the two horizontal sides, before finally trying
# "below"/"above" again at _LABEL_PUSH_GAP_PX clearance (a second, further-
# out tier in the same two directions, rather than new directions or a
# multiplied gap in all four): a "below" rect pushed out by at least one
# full label-height clears any earlier "below" rect regardless of either
# label's own text width (their x-ranges never have to be compared at all
# once their y-ranges are made disjoint by construction), which a
# multiplied-gap push in "right"/"left" can't promise the same way (there,
# the two rects' clearance instead depends on the OTHER label's own width,
# which isn't known when placing this one). A body whose label still
# collides in all 6 tries (an extremely dense cluster) is the only
# remaining case that drops its label entirely -- rare enough now that
# it's a last resort, not the common outcome the old 2-direction version
# made it.
_LABEL_DIRECTIONS = ("below", "above", "right", "left")
_LABEL_PUSH_DIRECTIONS = ("below", "above")
_LABEL_PUSH_GAP_PX = _LABEL_GAP_PX + _LABEL_HALF_HEIGHT_PX * 2 + 2.0

# A label never runs off the fixed viewBox (MAP.50): a "below"/"above"
# label crossing a side edge slides sideways back inside, a "right"/"left"
# one crossing the top or bottom slides up or down, and a candidate that
# would cross an edge it can't slide along (a "right" label at the right
# edge) isn't used.
_LABEL_EDGE_MARGIN_PX = 3.0

# A scene's drawn area, `(lo, hi)` on both axes: the fixed 700 px frame
# unless MAP.88's fit (`_fit_bounds`) widened the viewBox around the center
# so that everything drawn stays inside it.
_DEFAULT_BOUNDS = (0.0, _VIEW_SIZE_PX)


def _slide_into_view(lo, hi, bounds=_DEFAULT_BOUNDS):
    """How far to move the span `[lo, hi]` to sit inside `bounds` (less
    `_LABEL_EDGE_MARGIN_PX`), 0.0 when it already does."""
    if lo < bounds[0] + _LABEL_EDGE_MARGIN_PX:
        return bounds[0] + _LABEL_EDGE_MARGIN_PX - lo
    if hi > bounds[1] - _LABEL_EDGE_MARGIN_PX:
        return bounds[1] - _LABEL_EDGE_MARGIN_PX - hi
    return 0.0


def _edge_anchor(dx, bounds=_DEFAULT_BOUNDS):
    """`(x, text_anchor)` pinning a label slid by `dx` against the side
    edge it was slid away from."""
    if dx > 0:
        return bounds[0] + _LABEL_EDGE_MARGIN_PX, "start"
    return bounds[1] - _LABEL_EDGE_MARGIN_PX, "end"


def _fit_label_rect(rect, direction, bounds=_DEFAULT_BOUNDS):
    """`(rect, dx, dy)` with `rect` slid back inside the viewBox along the
    label's own line, or `None` when it crosses an edge across that line."""
    x1, y1, x2, y2 = rect
    if direction in ("below", "above"):
        dx, dy = _slide_into_view(x1, x2, bounds), 0.0
        if _slide_into_view(y1, y2, bounds):
            return None
    else:
        dx, dy = 0.0, _slide_into_view(y1, y2, bounds)
        if _slide_into_view(x1, x2, bounds):
            return None
    return (x1 + dx, y1 + dy, x2 + dx, y2 + dy), dx, dy


# MAP.88: nothing drawn runs past the frame. Once a scene's markers,
# orbits, belts and facilities are placed, `_fit_bounds` widens the square
# viewBox evenly around the center until every one of them (and each
# star's name below it) sits inside with `_FIT_MARGIN_PX` to spare, which
# shrinks the whole scene on screen; planet and moon names are then placed
# inside those bounds. A scene that already fits keeps the 700 px frame.
_FIT_MARGIN_PX = 6.0
_STAR_LABEL_DROP_PX = 20.0  # the star name's baseline (r + 15) plus its descenders
_HALO_EXTRA_PX = 9.0
_RING_REACH = 1.7  # a gas giant's ring ellipse, `rx` = 1.7 r


def _circle_box(cx, cy, r):
    return (cx - r, cy - r, cx + r, cy + r)


def _body_box(cx, cy, r_px, body_type=None, is_self=False):
    """Everything `_body_marker_svg` draws for a planet or moon but its
    name: the circle, a gas giant's ring, a moon scene's halo."""
    reach = r_px * _RING_REACH if body_type == "g" else r_px
    if is_self:
        reach = max(reach, r_px + _HALO_EXTRA_PX)
    return _circle_box(cx, cy, reach)


def _star_box(cx, cy, r_px):
    """A star marker and the name `_star_marker_svg` writes under it
    (sideways, that name slides back inside like any label)."""
    return (cx - r_px, cy - r_px, cx + r_px, cy + r_px + _STAR_LABEL_DROP_PX)


def _fit_bounds(boxes):
    """The square `(lo, hi)` around `_CENTER_PX` that holds every
    `(x1, y1, x2, y2)` in `boxes` with `_FIT_MARGIN_PX` to spare, and
    never less than the fixed frame."""
    half = _CENTER_PX
    for x1, y1, x2, y2 in boxes:
        half = max(half, _CENTER_PX - x1 + _FIT_MARGIN_PX, x2 - _CENTER_PX + _FIT_MARGIN_PX,
                   _CENTER_PX - y1 + _FIT_MARGIN_PX, y2 - _CENTER_PX + _FIT_MARGIN_PX)
    return (_CENTER_PX - half, _CENTER_PX + half)


def _label_half_width_px(text):
    """A cheap stand-in for actually measuring `text` at `.sysmap-label`'s
    font size (not available server-side, since this is a static SVG) --
    just enough to decide whether two labels would collide, not to lay out
    precisely."""
    return max(_LABEL_MIN_HALFWIDTH_PX, len(text or "") * _LABEL_CHAR_WIDTH_PX / 2)


def _rects_overlap(a, b):
    ax1, ay1, ax2, ay2 = a
    bx1, by1, bx2, by2 = b
    return ax1 < bx2 and ax2 > bx1 and ay1 < by2 and ay2 > by1


def _star_label_rect(cx, cy, r_px, name, bounds=_DEFAULT_BOUNDS):
    """The bounding box `_star_marker_svg` always draws a star's own name
    label in -- unconditionally below the marker, never collision-checked
    the way a planet/moon's own label is (`_label_sides_2d`) since a scene
    only ever has one or two of these. Exposed so a *planet's* label can
    still be kept out of it via `_label_sides_2d`'s `seed_rects` (see
    `_render_wide_binary_scenes`'s companion marker, which sits close
    enough to some of the primary's own far-out planets to otherwise
    collide with one)."""
    half_w = _label_half_width_px(name)
    top = cy + r_px + _LABEL_GAP_PX
    dx = _slide_into_view(cx - half_w, cx + half_w, bounds)
    return (cx - half_w + dx, top, cx + half_w + dx, top + _LABEL_HALF_HEIGHT_PX * 2)


def _label_candidate_rect(cx, cy, marker_r, half_w, direction, gap):
    """One candidate placement's bounding box for `direction` (one of
    `_LABEL_DIRECTIONS`), at `gap` clearance from the marker's own edge."""
    if direction == "below":
        top = cy + marker_r + gap
        return (cx - half_w, top, cx + half_w, top + _LABEL_HALF_HEIGHT_PX * 2)
    if direction == "above":
        bottom = cy - marker_r - gap
        return (cx - half_w, bottom - _LABEL_HALF_HEIGHT_PX * 2, cx + half_w, bottom)
    if direction == "right":
        left = cx + marker_r + gap
        return (left, cy - _LABEL_HALF_HEIGHT_PX, left + half_w * 2, cy + _LABEL_HALF_HEIGHT_PX)
    # "left"
    right = cx - marker_r - gap
    return (right - half_w * 2, cy - _LABEL_HALF_HEIGHT_PX, right, cy + _LABEL_HALF_HEIGHT_PX)


# The server's placement is the no-script layout; widths are estimated
# (_label_half_width_px), so static/systemmap.js's layoutLabels measures
# the real text once a scene is shown and nudges or hides any label that
# still overlaps.
def _label_sides_2d(entries, seed_rects=None, bounds=_DEFAULT_BOUNDS):
    """
    Given `[(cx, cy, marker_r, name), ...]`, returns one
    `{"direction", "rect", "pushed"}` dict (or `None` to skip drawing
    that label entirely) per entry: where that body's name label should
    be drawn. Unlike the old fixed schematic map's `_label_sides` (which
    only ever had to de-collide labels along one shared horizontal band,
    since every marker sat on the same line), markers here can be
    anywhere in the plane, so this checks each candidate label's own
    bounding box against every *other* label already placed, in true 2D,
    rather than assuming any shared axis.

    Tries all four cardinal directions (`_LABEL_DIRECTIONS`) at the normal
    gap first, then "below"/"above" again at `_LABEL_PUSH_GAP_PX` clearance
    before giving up -- a "pushed" placement is far enough from its own
    marker that `_body_marker_svg` draws a thin leader line connecting
    them, so a repositioned label still visibly belongs to its body. Only
    a body whose label collides in all 6 tries (an extremely dense
    cluster) drops its label entirely -- that body's name, and everything
    else about it, is still one click away in the info panel.

    Args:
        entries (list[tuple]): `(cx, cy, marker_r, name)`, in the order
            markers should be given placement priority (earlier entries
            never yield to a later one).
        seed_rects (list[tuple], optional): Extra already-occupied label
            rects (e.g. `_star_label_rect`'s) no candidate here may
            collide with, even though they belong to no entry in this
            call.
        bounds (tuple, optional): The scene's `(lo, hi)` drawn area
            (`_fit_bounds`) every label must stay inside.

    Returns:
        list[dict or None]: One entry per input, in the same order. Each
            dict has `direction` (one of `_LABEL_DIRECTIONS`), `rect` (the
            accepted bounding box), `pushed` (bool, whether this used
            the further-out "below"/"above" tier and so needs a leader
            line), and `shift` (`(dx, dy)`, how far the label slid to stay
            inside the view, see `_fit_label_rect`).
    """
    accepted_rects = list(seed_rects or [])
    placements = []
    for cx, cy, marker_r, name in entries:
        half_w = _label_half_width_px(name)
        chosen = None
        for pushed, gap, directions in (
            (False, _LABEL_GAP_PX, _LABEL_DIRECTIONS),
            (True, _LABEL_PUSH_GAP_PX, _LABEL_PUSH_DIRECTIONS),
        ):
            for direction in directions:
                fitted = _fit_label_rect(_label_candidate_rect(cx, cy, marker_r, half_w, direction, gap), direction, bounds)
                if fitted is None:
                    continue
                rect, dx, dy = fitted
                if not any(_rects_overlap(rect, other) for other in accepted_rects):
                    chosen = {"direction": direction, "rect": rect, "pushed": pushed, "shift": (dx, dy)}
                    break
            if chosen:
                break
        if chosen:
            accepted_rects.append(chosen["rect"])
        placements.append(chosen)
    return placements


def _data_attrs(attrs):
    return "".join(f' data-{key}="{esc(value)}"' for key, value in attrs.items() if value is not None)


_LIFE_BADGE_FILL = "#3ecf6e"
_LIFE_BADGE_STROKE = "#1c7a3e"


def _label_position(cx, cy, r_px, direction, gap):
    """The `(x, y, text_anchor)` a `.sysmap-label` `<text>` should render
    at for `direction`/`gap` (see `_label_candidate_rect`, which this
    mirrors so the rendered text actually lands inside the rect that was
    reserved for it during collision checking). "below"/"above" keep this
    module's original hand-tuned baseline offsets (10/-4 past the rect's
    own near edge) at the default gap; "right"/"left" are new."""
    if direction == "below":
        return cx, cy + r_px + gap + 10, "middle"
    if direction == "above":
        return cx, cy - r_px - gap - 4, "middle"
    if direction == "right":
        return cx + r_px + gap, cy + 4, "start"
    # "left"
    return cx - r_px - gap, cy + 4, "end"


def _leader_line_svg(cx, cy, r_px, direction, gap):
    """A short stub `<line>` from the marker's own edge to just short of
    where its label begins -- drawn whenever a label isn't in its default,
    immediately-adjacent "below" spot (see `_body_marker_svg`), so a label
    `_label_sides_2d` had to reposition to avoid a collision still reads as
    visually connected to its own body rather than looking like a stray,
    unrelated label floating nearby."""
    if direction == "below":
        x1, y1, x2, y2 = cx, cy + r_px, cx, cy + r_px + gap
    elif direction == "above":
        x1, y1, x2, y2 = cx, cy - r_px, cx, cy - r_px - gap
    elif direction == "right":
        x1, y1, x2, y2 = cx + r_px, cy, cx + r_px + gap, cy
    else:  # "left"
        x1, y1, x2, y2 = cx - r_px, cy, cx - r_px - gap, cy
    return f'<line class="sysmap-label-leader" x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}"></line>'


def _body_marker_svg(cx, cy, r_px, planet_class, body_type, label_text, extra_class, attrs, is_self=False,
                      label_direction="below", label_pushed=False, show_label=True, has_life=False,
                      label_shift=(0.0, 0.0), bounds=_DEFAULT_BOUNDS):
    """
    Builds one clickable `<g>` for a planet or moon: a filled/stroked
    circle colored by `planet_class` (see `_CLASS_COLORS`), the class
    letter centered inside it once the circle is big enough to hold text
    legibly, a small tilted ring behind gas giants (`body_type == "g"`)
    for an at-a-glance silhouette cue beyond just color, a small green
    "supports life" badge at the marker's own edge when `has_life` is set
    (`planets`/`moons`.`life_chemical` is non-NULL), and (when
    `show_label` is set) a name label at `label_direction` (one of
    `_LABEL_DIRECTIONS` -- "below", this function's own default, needs no
    further reasoning; any other value, or `label_pushed` being set,
    means `_label_sides_2d` had to move this label to avoid colliding with
    a neighbor, so a short leader line is drawn connecting it back to this
    marker). `show_label` alone going False (also from `_label_sides_2d`,
    when no direction at any distance could clear a tight enough cluster)
    only drops the visible text, not the body's name from `aria-label` --
    it's still reachable, just via a click on the marker rather than a
    glance at the diagram. `is_self` marks the one body a moon-scene is
    *about* (the planet drilled into, now standing in for the scene's own
    origin) -- drawn with a soft halo and picked out by
    `static/systemmap.js` as the default info-panel content when that
    scene opens, via its `data-self="true"` marker.
    """
    fill = class_color(planet_class)
    stroke = _darken_hex(fill, 0.22)
    text_color = _text_color_for(fill)

    parts = []
    if is_self:
        parts.append(f'<circle class="sysmap-halo" cx="{cx:.1f}" cy="{cy:.1f}" r="{r_px + 9:.1f}" fill="{fill}"></circle>')
    if body_type == "g":
        parts.append(
            f'<ellipse class="sysmap-ring" cx="{cx:.1f}" cy="{cy:.1f}" rx="{r_px * 1.7:.1f}" ry="{r_px * 0.5:.1f}" '
            f'transform="rotate(-20 {cx:.1f} {cy:.1f})"></ellipse>'
        )
    parts.append(
        f'<circle class="sysmap-body-fill" cx="{cx:.1f}" cy="{cy:.1f}" r="{r_px:.1f}" '
        f'fill="{fill}" stroke="{stroke}"></circle>'
    )
    if planet_class and r_px >= 8:
        parts.append(
            f'<text class="sysmap-classletter" x="{cx:.1f}" y="{cy:.1f}" fill="{text_color}" '
            f'text-anchor="middle" dominant-baseline="central">{esc(planet_class.upper())}</text>'
        )
    if has_life:
        badge_r = max(2.5, r_px * 0.32)
        badge_x = cx + r_px * 0.72
        badge_y = cy - r_px * 0.72
        parts.append(
            f'<circle class="sysmap-life-badge" cx="{badge_x:.1f}" cy="{badge_y:.1f}" r="{badge_r:.1f}" '
            f'fill="{_LIFE_BADGE_FILL}" stroke="{_LIFE_BADGE_STROKE}"><title>Supports life</title></circle>'
        )
    if label_text and show_label:
        gap = _LABEL_PUSH_GAP_PX if label_pushed else _LABEL_GAP_PX
        # "below"/"above" at the normal gap need no leader -- immediately
        # adjacent, above/below is the obvious default reading, exactly
        # like this map's original (pre-collision-avoidance) behavior.
        # "right"/"left" get one even at the normal gap, since a label
        # floating to a marker's side isn't as self-evidently "this
        # marker's own name" without one; the further-out pushed tier
        # (below/above only -- see _LABEL_PUSH_DIRECTIONS) always does too.
        needs_leader = label_pushed or label_direction in ("right", "left")
        if needs_leader:
            parts.append(_leader_line_svg(cx, cy, r_px, label_direction, gap))
        label_x, label_y, anchor = _label_position(cx, cy, r_px, label_direction, gap)
        label_y += label_shift[1]
        if label_shift[0]:
            # Slid off a side edge: pin the text to that edge, so the real
            # (unestimated) width still sits right against it.
            label_x, anchor = _edge_anchor(label_shift[0], bounds)
        parts.append(
            f'<text class="sysmap-label" x="{label_x:.1f}" y="{label_y:.1f}" text-anchor="{anchor}">{esc(label_text)}</text>'
        )

    self_attr = ' data-self="true"' if is_self else ""
    return (
        f'<g class="sysmap-body {extra_class}" tabindex="0" role="button"{self_attr}{_data_attrs(attrs)} '
        f'aria-label="{esc(label_text or "body")}">{"".join(parts)}</g>'
    )


def _star_marker_svg(cx, cy, r_px, star, attrs, bounds=_DEFAULT_BOUNDS):
    fill, stroke = star_color(star["star_type"], star["temperature_k"], star["luminosity_w"])
    name = attrs.get("name", "star")
    # Handed to the client as `data-color` too (like a planet/moon's own
    # "color" attr -- see `_planet_attrs`) so `static/systemmap.js`'s 3D
    # sphere renderer can color a star's own live-rendered sphere without
    # duplicating `star_color`'s spectral-type logic in JS.
    star_attrs = dict(attrs)
    star_attrs["color"] = fill
    half_w = _label_half_width_px(name)
    label_x, anchor = cx, "middle"
    if _slide_into_view(cx - half_w, cx + half_w, bounds):
        label_x, anchor = _edge_anchor(_slide_into_view(cx - half_w, cx + half_w, bounds), bounds)
    return (
        f'<g class="sysmap-body sysmap-star" tabindex="0" role="button"{_data_attrs(star_attrs)} '
        f'aria-label="{esc(name)}">'
        f'<circle class="sysmap-body-fill" cx="{cx:.1f}" cy="{cy:.1f}" r="{r_px:.1f}" '
        f'fill="{fill}" stroke="{stroke}"></circle>'
        f'<text class="sysmap-label sysmap-star-label" x="{label_x:.1f}" y="{cy + r_px + 15:.1f}" '
        f'text-anchor="{anchor}">{esc(name)}</text>'
        "</g>"
    )


def _atmosphere_text(planet):
    atmosphere = planet.get("atmosphere")
    return atmosphere if atmosphere and atmosphere != "None" else "None (airless)"


def _pressure_text(planet, has_atmosphere):
    """MAP.117: the side panel's surface pressure, "None" for an airless
    body, left out when the row has no value."""
    if not has_atmosphere:
        return "None"
    pressure_pa = planet.get("atmospheric_pressure_pa")
    return None if pressure_pa is None else to_plain_text(format_pressure_pa(pressure_pa))


def _planet_attrs(planet, kind="planet", parent_name=None, scene_target=None):
    atmosphere = planet.get("atmosphere")
    has_atmosphere = bool(atmosphere) and atmosphere != "None"
    surface_temp_k = planet.get("surface_temperature_k")
    attrs = {
        "kind": kind,
        "id": planet["id"],
        "name": planet["name"],
        "class": (planet["planet_class"] or "").upper(),
        "classdesc": _class_description(planet["planet_class"]),
        # The same class->color lookup the marker's own fill uses
        # (`class_color`) -- handed to the client as a resolved hex string
        # (like `starmap.py`'s star colors) rather than duplicating
        # `_CLASS_COLORS` in JS, so this module stays the one place a
        # planet class's color is decided. Drives this marker's own live 3D
        # sphere in `static/systemmap.js` (`#sysmap-spheres-canvas`).
        "color": class_color(planet["planet_class"]),
        "bodytype": "Gas Giant" if planet["body_type"] == "g" else "Terrestrial",
        # MAP.92: the side panel's radius and mass, formatted here like
        # every other field.
        "radius": to_plain_text(format_body_radius(planet.get("radius_km"))),
        "mass": to_plain_text(format_body_mass(planet.get("mass_kg"), planet["body_type"] == "g")),
        "zone": _ZONE_LABELS.get(planet["zone"], ""),
        "distance": to_plain_text(format_body_distance(planet["distance_km"], planet.get("_is_moon", False))),
        "period": format_period(planet["period_years"]),
        "gravity": f'{round(planet["gravity_g"], 3) if planet["gravity_g"] is not None else ""} g',
        "life": planet.get("life_chemical"),
        "atmosphere": _atmosphere_text(planet),
        "composition": planet.get("composition"),
        "surfacetemp": format_temperature_k(surface_temp_k) if surface_temp_k is not None else None,
        "surfacepressure": _pressure_text(planet, has_atmosphere),
        # Presence alone (not the description text) is what the 3D preview
        # needs to decide whether to draw an atmosphere glow shell at all --
        # `_data_attrs` omits a `None` value entirely, so this attribute's
        # mere presence in the DOM is the boolean.
        "hasatmosphere": "true" if has_atmosphere else None,
        # Raw, real (not log-scaled, not pixel) coordinates -- relative to
        # this body's own anchor (the star it orbits, or the planet it
        # moons), same units `_star_scene_svg`/`_render_system_scene`
        # themselves compute `cx`/`cy` from. `static/systemmap.js`'s own
        # "measure distance" feature reads these directly rather than
        # trying to reconstruct a real km distance from two markers'
        # drawn pixel positions, which the shared log radial scale
        # (`_radial_px`) makes lossy for that purpose -- real angle is
        # preserved, but real *distance* isn't recoverable from pixels
        # alone. `radiuskm` (this body's own physical size) is included
        # too so the same feature can treat it as a "can't route through
        # this body" obstacle when it's the scene's own center (a moon
        # scene's `data-self="true"` planet) -- unused, harmless, for
        # every other body.
        "xkm": planet.get("position_x_km") or 0.0,
        "ykm": planet.get("position_y_km") or 0.0,
        "radiuskm": planet.get("radius_km"),
        "note": planet.get("_position_note"),
    }
    if parent_name is not None:
        attrs["parent"] = parent_name
    if scene_target is not None:
        attrs["scene"] = scene_target
        attrs["moons"] = len(planet.get("moons") or [])
    return attrs


def _belt_band(belt, lo_km, hi_km, orbit_radii_px=()):
    """
    The drawn ring for `belt`: `(center_px, width_px)`.

    Its inner and outer edges are the belt's own `lower_limit_km` and
    `upper_limit_km` through `_radial_px`, the scale every orbit uses (a
    belt's `distance_km` is its inner edge, not its middle, so centering
    the ring on it drew half the ring inside the belt). A ring thinner than
    `_BELT_MIN_BAND_PX` is widened about its middle, but each edge stops
    `_BELT_ORBIT_GAP_PX` short of the nearest orbit in `orbit_radii_px`
    (the drawn orbit radii around the same anchor) that lies outside it.
    """
    inner = _radial_px(belt["lower_limit_km"] or belt["distance_km"], lo_km, hi_km)
    outer = _radial_px(belt["upper_limit_km"] or belt["distance_km"], lo_km, hi_km)
    inner, outer = min(inner, outer), max(inner, outer)
    if outer - inner < _BELT_MIN_BAND_PX:
        middle = (inner + outer) / 2
        inner_limit = max((r + _BELT_ORBIT_GAP_PX for r in orbit_radii_px if r < inner), default=0.0)
        outer_limit = min((r - _BELT_ORBIT_GAP_PX for r in orbit_radii_px if r > outer), default=math.inf)
        inner = max(min(inner, middle - _BELT_MIN_BAND_PX / 2), min(inner, inner_limit))
        outer = min(max(outer, middle + _BELT_MIN_BAND_PX / 2), max(outer, outer_limit))
    return (inner + outer) / 2, outer - inner


def _belt_ring_svg(cx, cy, radius_px, band_px, belt):
    """A belt is a *range* of orbit radii, drawn as a full ring (a plain
    stroked circle, `stroke-width` = its own radial spread) around its own
    anchor -- every angle at roughly this distance, the natural true-
    position shape, unlike the old schematic map's one-directional shaded
    band (which only ever made sense when every body was forced due-east
    of its star to begin with)."""
    attrs = {
        "kind": "belt",
        "id": belt["id"],
        "name": "Asteroid Belt",
        "density": belt["density"].capitalize(),
        "distance": f'{format_distance_km(belt["lower_limit_km"])} - {format_distance_km(belt["upper_limit_km"])}',
        "composition": belt["composition_summary"],
    }
    label = f'Asteroid belt ({belt["density"]})'
    return (
        f'<circle class="sysmap-belt" tabindex="0" role="button"{_data_attrs(attrs)} '
        f'aria-label="{esc(label)}" cx="{cx:.1f}" cy="{cy:.1f}" r="{radius_px:.1f}" '
        f'stroke-width="{band_px:.1f}"></circle>'
    )


# --- Facilities (schema v42) -------------------------------------------------
#
# A facility is a small diamond at its host: an orbital one on a dashed
# orbit around its star (on the scene's own log scale) or just outside its
# planet's or moon's marker (a facility's orbit is far too small to show on
# the system's scale, so it hugs its host the way the host's real orbit
# would look from that far out), a terrestrial one on its world's rim, and
# a belt one on the belt's ring. Its angle is its stored `orbit_phase_deg`;
# a facility with none gets a fixed angle from its id. Clicking one fills
# the info panel from its `data-*` attributes like any other marker; it is
# not a `.sysmap-body`, so it gets no 3D sphere and isn't measurable.

_FACILITY_HALF_PX = 4.5
_FACILITY_HOST_GAP_PX = 6.0


def _facility_angle_deg(facility):
    phase = facility.get("orbit_phase_deg")
    return phase if phase is not None else (int(str(facility["id"]).replace("-", ""), 16) % 4096 * 137.5) % 360.0


def _facility_point(anchor_x_px, anchor_y_px, radius_px, facility):
    theta = math.radians(_facility_angle_deg(facility))
    return anchor_x_px + radius_px * math.cos(theta), anchor_y_px - radius_px * math.sin(theta)


def _facility_marker_svg(cx, cy, facility):
    """One facility's clickable diamond."""
    kind_label = facility["kind"].replace("-", " ").capitalize()
    orbital = facility.get("orbit_distance_km") is not None
    attrs = {
        "kind": "facility",
        "name": facility["name"],
        "facilitykind": kind_label,
        "host": facility.get("host_name") or facility["host_type"].replace("_", " ").capitalize(),
        "placement": {"terrestrial": "On the surface", "orbital": "In orbit", "asteroid": "Among the asteroids"}
        .get(facility["placement"], facility["placement"]),
        "distance": format_distance_km(facility["orbit_distance_km"]) if orbital else None,
        "period": format_period(facility["orbit_period_years"])
        if orbital and facility.get("orbit_period_years") is not None else None,
        "speed": format_speed_kms(facility["orbital_speed_kms"])
        if orbital and facility.get("orbital_speed_kms") is not None else None,
    }
    h = _FACILITY_HALF_PX
    points = f"{cx:.1f},{cy - h:.1f} {cx + h:.1f},{cy:.1f} {cx:.1f},{cy + h:.1f} {cx - h:.1f},{cy:.1f}"
    label = f'{facility["name"]} ({kind_label})'
    return (
        f'<g class="sysmap-facility sysmap-facility-{esc(facility["kind"])}" tabindex="0" role="button"'
        f'{_data_attrs(attrs)} aria-label="{esc(label)}">'
        f'<polygon points="{points}"></polygon><title>{esc(label)}</title></g>'
    )


def _facility_orbit_svg(cx, cy, radius_px):
    return (f'<circle class="sysmap-facility-orbit" cx="{cx:.1f}" cy="{cy:.1f}" '
            f'r="{radius_px:.1f}"></circle>')


def _facilities_svg(facilities, hosts, scale):
    """
    Every facility whose host is drawn in this scene.

    Args:
        facilities (list[dict]): `GET /api/systems/<id>/facilities` items.
        hosts (dict): `{(host_type, host_id): (cx, cy, r_px, centered)}`
            for each marker drawn in this scene, `r_px` its drawn radius
            (a belt's ring radius). `centered` marks the host the scene is
            built around (a star, or a moon scene's planet): its orbital
            facilities go on the scene's own log scale.
        scale (tuple): The scene's `(lo_km, hi_km)`.

    Returns:
        tuple[str, str, list]: `(orbit lines, markers, boxes)`, `boxes`
            the drawn extent of each (`_fit_bounds`).
    """
    orbits, markers, boxes = [], [], []
    for facility in facilities or ():
        host = hosts.get((facility["host_type"], facility["host_id"]))
        if host is None:
            continue
        cx, cy, r_px, centered = host
        if facility["host_type"] == "asteroid_belt":
            radius_px = r_px
        elif facility["placement"] != "orbital":
            radius_px = r_px
        elif centered:
            radius_px = max(_radial_px(facility["orbit_distance_km"], *scale), r_px + _FACILITY_HOST_GAP_PX)
            orbits.append(_facility_orbit_svg(cx, cy, radius_px))
            boxes.append(_circle_box(cx, cy, radius_px))
        else:
            radius_px = r_px + _FACILITY_HOST_GAP_PX
            orbits.append(_facility_orbit_svg(cx, cy, radius_px))
            boxes.append(_circle_box(cx, cy, radius_px))
        fx, fy = _facility_point(cx, cy, radius_px, facility)
        markers.append(_facility_marker_svg(fx, fy, facility))
        boxes.append(_circle_box(fx, fy, _FACILITY_HALF_PX))
    return "".join(orbits), "".join(markers), boxes


def _binary_star_positions_km(system, stars):
    """
    Splits `star_systems.binary_mutual_position_x/y/z_km` (the secondary's
    position relative to the primary) into each star's own position
    relative to their shared mass-weighted barycenter -- `primary * m2 +
    secondary * m1 = 0` around that point, by definition of a barycenter --
    rather than treating the primary as a fixed, unmoving anchor. Used both
    to actually place each star's own marker and, for a 'close' (P-type)
    pair, as the shared origin every planet/belt in the system also orbits
    (`doubleStar.BinaryStarProxy` merges the pair into one *effective*
    star sitting at that same barycenter -- see its own docstring).

    Args:
        system (dict): `queryDb.system_detail`'s return shape -- reads
            `binary_mutual_position_x/y/z_km`.
        stars (list[dict]): Exactly 2 entries (primary, then secondary),
            each with `id`/`mass_kg`.

    Returns:
        dict: `{star_id: (x_km, y_km)}` for both stars.
    """
    primary, secondary = stars[0], stars[1]
    m1 = primary.get("mass_kg") or 0.0
    m2 = secondary.get("mass_kg") or 0.0
    total = m1 + m2
    bx = system.get("binary_mutual_position_x_km") or 0.0
    by = system.get("binary_mutual_position_y_km") or 0.0
    primary_frac = (m2 / total) if total else 0.5
    secondary_frac = (m1 / total) if total else 0.5
    return {
        primary["id"]: (-bx * primary_frac, -by * primary_frac),
        secondary["id"]: (bx * secondary_frac, by * secondary_frac),
    }


def _scene_svg_pair(scene_id, aria_label, hidden, orbits_inner, bodies_inner, scale=None, bounds=_DEFAULT_BOUNDS):
    """
    Wraps a scene's own orbit-line markup and body-marker markup as TWO
    sibling `<svg data-scene="...">` elements sharing the same
    `data-scene`/`viewBox`/hidden state, instead of the one combined
    `<svg>` every scene-builder function used to return.

    Why: `static/systemmap.js`'s per-marker sphere canvas
    (`#sysmap-spheres-canvas`) sits at a fixed z-index relative to the
    SVG content -- but within a single `<svg>`, a body's own opaque
    sphere (drawn on an entirely separate `<canvas>` element, never as
    SVG content of its own) has no way to occlude an orbit line that's a
    sibling SVG element painted in that very same stacking context, which
    is why an orbit line used to visibly cut across every rendered sphere
    it passed under/through. Splitting each scene into an orbits-only
    layer and a body-marker layer lets `style.css` stack them either side
    of that canvas (`.sysmap-orbits-layer`'s own z-index rule) instead.

    `static/systemmap.js`'s `showScene` already toggles every
    `.sysmap-svg` sharing one `data-scene` value together (see that
    function's own `orbitLayers` handling) -- both `<svg>`s below carry
    the same `class="sysmap-svg"` (the orbits one with the additional
    `sysmap-orbits-layer` modifier `style.css` positions/z-indexes by),
    so no further JS wiring is needed to keep them in lockstep.

    Args:
        scene_id (str): Shared `data-scene` value for both `<svg>`s.
        aria_label (str): The body-marker layer's own `aria-label` (raw,
            not pre-escaped -- this function escapes it). The orbits
            layer is `aria-hidden` instead -- it carries no information
            of its own beyond what the body markers already expose.
        hidden (bool): Whether both `<svg>`s start `sysmap-hidden`.
        orbits_inner (str): This scene's own decorative `<circle
            class="sysmap-orbit">` markup (already-built SVG fragments) --
            NOT a belt ring, which (unlike a plain orbit line) is itself a
            real focusable/clickable marker (`tabindex="0" role="button"`,
            see `_belt_ring_svg`) and belongs in `bodies_inner` instead,
            alongside the other body markers -- this layer is
            `aria-hidden`, which must never contain focusable content.
        bodies_inner (str): This scene's own star/planet/moon/belt marker
            markup (already-built SVG fragments).
        scale (tuple, optional): `(lo_km, hi_km)`, the scene's own
            `_radial_scale_bounds`. Written onto the body-marker layer as
            `data-lokm`/`data-hikm` (with `data-cpx`, the scene's center)
            so `static/systemmap.js` can draw a measured route with the
            same log scale the markers were placed with.
        bounds (tuple, optional): The scene's `(lo, hi)` drawn area
            (`_fit_bounds`), which becomes both `<svg>`s' viewBox.

    Returns:
        str: Two concatenated sibling `<svg>` elements.
    """
    scale_attrs = ""
    if scale is not None:
        scale_attrs = (f' data-lokm="{scale[0]:.6g}" data-hikm="{scale[1]:.6g}" data-cpx="{_CENTER_PX:.1f}"'
                       f' data-minpx="{_MIN_RADIUS_PX:.1f}" data-spreadpx="{_RADIUS_SPREAD_PX:.1f}"')
    hidden_class = " sysmap-hidden" if hidden else ""
    lo, size = round(bounds[0], 1), round(bounds[1] - bounds[0], 1)
    view_box = f'viewBox="{lo:g} {lo:g} {size:g} {size:g}"'
    orbits_svg = (
        f'<svg class="sysmap-svg sysmap-orbits-layer{hidden_class}" data-scene="{esc(scene_id)}" '
        f'{view_box} aria-hidden="true">{orbits_inner}</svg>'
    )
    bodies_svg = (
        f'<svg class="sysmap-svg{hidden_class}" data-scene="{esc(scene_id)}" '
        f'{view_box}{scale_attrs} role="group" aria-label="{esc(aria_label)}">{bodies_inner}</svg>'
    )
    return orbits_svg + bodies_svg


def _star_scene_svg(scene_id, aria_label, hidden, star, planets, belts, star_attrs, extra_svg="", extra_obstacle=None,
                    facilities=None):
    """
    Builds one star-centered scene: the given `star` fixed at this scene's
    own origin with its own `planets`/`belts` arranged around it on their
    own dedicated log scale (`_radial_scale_bounds`, computed from nothing
    but this star's own bodies) -- shared by the "system" scene a
    single-star or 'wide' (S-type) binary's primary gets
    (`_render_system_scene`) and the secondary's own drilled-into scene
    (`_render_wide_secondary_scene`), so each star's planetary system
    always gets this diagram's full radial pixel budget to itself, never
    sharing it with a companion star's own (potentially vastly larger,
    tens-to-thousands-of-AU) orbital separation the way one shared scale
    across the whole pair used to force.

    Args:
        scene_id (str): This `<svg>`'s own `data-scene` value.
        aria_label (str): This `<svg>`'s own `aria-label`.
        hidden (bool): Whether this scene starts hidden (`sysmap-hidden`) --
            `False` only for the one scene `static/systemmap.js` shows by
            default (`data-scene="system"`).
        star (dict): The star anchoring this scene.
        planets (list[dict]): Only this star's own planets.
        belts (list[dict]): Only this star's own asteroid belts.
        star_attrs (dict): `_data_attrs`-ready attributes for the star's
            own marker (name/role/kind/etc, plus `scene`/`self` when this
            star should itself be clickable into a deeper scene, or is the
            scene's own "you are here" body -- see `_data_attrs`).
        extra_svg (str): Extra raw SVG appended after every body in this
            scene (e.g. a wide pair's companion-star marker, in the
            primary's own "system" scene only).
        extra_obstacle (dict, optional): `{"cx", "cy", "r", "label_rect"}`
            for a marker drawn separately (via `extra_svg`) that this
            scene's own planets/belts must still be kept clear of -- the
            companion-star marker `extra_svg` draws isn't one of this
            scene's own relaxed `markers`, so without this it would be
            invisible to `_relax_markers`/`_label_sides_2d` and a
            far-enough-out planet could still be placed right on top of
            it (or its label) instead of being pushed aside like any two
            of this scene's own bodies already are from each other.
            `label_rect` (a `_star_label_rect`) is optional within this
            dict -- omitted when the obstacle has no label of its own to
            avoid. `box` (also optional) is the drawn extent the scene's
            fit (`_fit_bounds`) must keep inside the frame.
        facilities (list[dict], optional): The system's facilities; those
            on this star, its planets and its belts are drawn
            (`_facilities_svg`).

    Returns:
        str: A complete `<svg class="sysmap-svg">` scene.
    """
    local_r_list = []
    for planet in planets:
        local_r_list.append(_orbit_radius_km(planet))
    for belt in belts:
        local_r_list.extend((belt["lower_limit_km"], belt["upper_limit_km"]))
    lo, hi = _radial_scale_bounds(local_r_list)

    star_r = _star_radius_px(star["radius_km"])
    orbit_paths = []
    # Kept OUT of orbit_paths (which _scene_svg_pair puts in the
    # aria-hidden, non-interactive orbits layer) -- unlike a plain
    # sysmap-orbit circle, a belt ring is itself a real clickable/
    # focusable marker (tabindex="0" role="button", see
    # _belt_ring_svg), so it belongs with the other body markers in the
    # scene's normal, accessible layer instead.
    belt_svgs = []
    markers = [{"type": "star", "cx": _CENTER_PX, "cy": _CENTER_PX, "r": star_r, "fixed": True}]
    if extra_obstacle is not None:
        markers.append({
            "type": "obstacle", "cx": extra_obstacle["cx"], "cy": extra_obstacle["cy"],
            "r": extra_obstacle["r"], "fixed": True,
        })
    orbit_radii = []
    boxes = [_star_box(_CENTER_PX, _CENTER_PX, star_r)]
    if extra_obstacle is not None and extra_obstacle.get("box"):
        boxes.append(extra_obstacle["box"])
    for planet in planets:
        lx_km, ly_km = planet.get("position_x_km") or 0.0, planet.get("position_y_km") or 0.0
        r_px = _radial_px(_orbit_radius_km(planet), lo, hi)
        orbit_radii.append(r_px)
        boxes.append(_circle_box(_CENTER_PX, _CENTER_PX, r_px))
        cx, cy = _polar_to_px(_CENTER_PX, _CENTER_PX, r_px, lx_km, ly_km)
        orbit_paths.append(f'<circle class="sysmap-orbit" cx="{_CENTER_PX:.1f}" cy="{_CENTER_PX:.1f}" r="{r_px:.1f}"></circle>')
        markers.append({"type": "planet", "cx": cx, "cy": cy, "r": _planet_radius_px(planet["radius_km"]), "row": planet})
    facility_hosts = {("star", star["id"]): (_CENTER_PX, _CENTER_PX, star_r, True)}
    for belt in belts:
        r_px, band_px = _belt_band(belt, lo, hi, orbit_radii)
        belt_svgs.append(_belt_ring_svg(_CENTER_PX, _CENTER_PX, r_px, band_px, belt))
        boxes.append(_circle_box(_CENTER_PX, _CENTER_PX, r_px + band_px / 2))
        facility_hosts[("asteroid_belt", belt["id"])] = (_CENTER_PX, _CENTER_PX, r_px, False)

    _relax_markers(markers, _MARKER_GAP_PX)
    planet_markers = [m for m in markers if m["type"] == "planet"]
    for marker in planet_markers:
        facility_hosts[("planet", marker["row"]["id"])] = (marker["cx"], marker["cy"], marker["r"], False)
        boxes.append(_body_box(marker["cx"], marker["cy"], marker["r"], marker["row"]["body_type"]))
    facility_orbits, facility_markers, facility_boxes = _facilities_svg(facilities, facility_hosts, (lo, hi))
    bounds = _fit_bounds(boxes + facility_boxes)

    star_svg = _star_marker_svg(_CENTER_PX, _CENTER_PX, star_r, star, star_attrs, bounds)

    seed_rects = [extra_obstacle["label_rect"]] if extra_obstacle and extra_obstacle.get("label_rect") else None
    placements = _label_sides_2d(
        [(m["cx"], m["cy"], m["r"], m["row"]["name"]) for m in planet_markers], seed_rects=seed_rects, bounds=bounds,
    )
    body_svgs = []
    for marker, placement in zip(planet_markers, placements):
        row = marker["row"]
        scene_target = f'planet-{row["id"]}' if row.get("moons") else None
        attrs = _planet_attrs(row, kind="planet", scene_target=scene_target)
        body_svgs.append(_body_marker_svg(
            marker["cx"], marker["cy"], marker["r"], row["planet_class"], row["body_type"], row["name"],
            "sysmap-planet", attrs,
            label_direction=(placement["direction"] if placement else "below"),
            label_pushed=bool(placement and placement["pushed"]),
            label_shift=(placement["shift"] if placement else (0.0, 0.0)),
            show_label=(placement is not None),
            has_life=bool(row.get("life_chemical")),
            bounds=bounds,
        ))

    return _scene_svg_pair(
        scene_id, aria_label, hidden,
        "".join(orbit_paths) + facility_orbits,
        # belt_svgs drawn first, same relative "underneath" stacking its
        # old spot (inside orbit_paths, before star_svg/body_svgs) had.
        f'{"".join(belt_svgs)}{star_svg}{"".join(body_svgs)}{extra_svg}{facility_markers}',
        scale=(lo, hi), bounds=bounds,
    )


def _star_temp_text(star):
    temperature_k = star.get("temperature_k")
    return f"{int(temperature_k)} K" if temperature_k is not None else "\u2013"


def _star_label(system, star):
    """A star's shown name: its own (a binary's stars are `<system> <word>`
    -- see `bodyNames.py`), falling back to the system's."""
    return star.get("name") or system["name"]


def _wide_binary_star_attrs(system, star, is_primary, scene_target=None):
    """Shared `star_attrs` dict for `_star_scene_svg`'s own star marker, in
    either of a wide (S-type) pair's two scenes."""
    attrs = {
        "kind": "star", "id": star["id"], "name": _star_label(system, star),
        "role": "Primary" if is_primary else "Secondary",
        "type": star["star_type"], "temp": _star_temp_text(star),
        "mass": to_plain_text(format_star_mass(star["mass_kg"])),
        "radius": to_plain_text(format_star_radius(star["radius_km"])),
        "lum": to_plain_text(format_star_luminosity(star["luminosity_w"])),
        # Raw km, alongside the formatted display string above -- see
        # `_planet_attrs`'s identical `radiuskm` for why (the "measure
        # distance" feature's own around-the-star obstacle radius).
        "radiuskm": star["radius_km"],
        # This star, viewed in its own scene, always sits at that scene's
        # own local origin -- exactly like a moon scene's own drilled-
        # into planet. `_render_wide_binary_scenes` overrides both keys
        # for the one case that isn't this star's own scene: the
        # companion marker it draws *inside the other star's* scene,
        # where the real separation (not 0, 0) is what a "measure
        # distance" click needs.
        "xkm": 0.0, "ykm": 0.0,
    }
    if scene_target is not None:
        attrs["scene"] = scene_target
    return attrs


def _render_wide_binary_scenes(system, stars, planets, belts, facilities=None):
    """
    Builds a wide (S-type) binary's two scenes: the primary's own "system"
    scene (default-visible) and the secondary's own scene, reached by
    clicking its small companion marker in the primary's scene -- mirroring
    `_render_moon_scene`'s "drill into it" pattern one level up (a star,
    not a planet), per the module docstring's note on why a wide pair's
    true star-to-star separation (tens to thousands of AU -- see
    `wideBinary.py`'s own module docstring) can't share one radial pixel
    budget with either star's own, much smaller planetary system: sharing
    one budget either crushed both stars' planets down near the frame's
    center to make room for the real separation, or -- since the two
    stars' own drawn positions themselves then sat close to the frame's
    outer edge -- pushed their planets (and label text) straight off the
    visible canvas. Each star's own planets get this diagram's full
    radial budget in its own scene instead; only the companion *marker*
    (no planets of its own drawn in this scene) uses a fixed, merely
    representative distance along the real direction to it.

    Args:
        system (dict): As `render_system_map_panel` receives it.
        stars (list[dict]): Exactly 2 entries, primary then secondary.
        planets (list[dict]): Every planet in the system (both stars'
            own -- distinguished by `star_id`).
        belts (list[dict]): Ditto, for asteroid belts.

    Returns:
        list[str]: `[primary_scene_svg, secondary_scene_svg]`.
    """
    primary, secondary = stars[0], stars[1]
    primary_planets = [p for p in planets if p.get("star_id") == primary["id"]]
    primary_belts = [b for b in belts if b.get("star_id") == primary["id"]]
    secondary_planets = [p for p in planets if p.get("star_id") == secondary["id"]]
    secondary_belts = [b for b in belts if b.get("star_id") == secondary["id"]]

    secondary_scene_id = f'star-{secondary["id"]}'

    # The companion marker's real direction (from `binary_mutual_position_
    # x/y_km`, "the secondary's position relative to the primary" -- see
    # `queryDb.system_detail`'s docstring), but a fixed, merely
    # representative distance -- this scene's own outer edge -- rather
    # than the real separation itself, which is routinely 10-1000x any
    # planet's own distance from its star (see `wideBinary.py`'s sampled
    # range) and would otherwise place the marker far outside this
    # diagram's fixed frame.
    bx = system.get("binary_mutual_position_x_km") or 0.0
    by = system.get("binary_mutual_position_y_km") or 0.0
    companion_r = _star_radius_px(secondary["radius_km"])
    companion_radius_px = _MIN_RADIUS_PX + _RADIUS_SPREAD_PX
    companion_cx, companion_cy = _polar_to_px(_CENTER_PX, _CENTER_PX, companion_radius_px, bx, by)
    companion_attrs = _wide_binary_star_attrs(system, secondary, is_primary=False, scene_target=secondary_scene_id)
    # Override the "own scene" (0, 0) default `_wide_binary_star_attrs`
    # sets -- this marker represents the secondary at its real position
    # *relative to the primary*, not its own scene's origin, and that
    # real separation (bx, by) is exactly what a "measure distance" click
    # needs, decoupled from `companion_cx`/`companion_cy`'s own fixed,
    # merely-representative drawn position (see this function's own
    # docstring on why those two can't be the same value here).
    companion_attrs["xkm"], companion_attrs["ykm"] = bx, by
    companion_marker_svg = _star_marker_svg(companion_cx, companion_cy, companion_r, secondary, companion_attrs)
    companion_obstacle = {
        "cx": companion_cx, "cy": companion_cy,
        # Inflated past the marker's own drawn radius to also cover its
        # label's footprint (always directly below it -- see
        # `_star_marker_svg`) -- `_relax_markers` only ever repels by
        # circle, so without this a planet's own circle could clear the
        # companion's circle only to still land its *label* on top of
        # "Esarer B"'s own name underneath it.
        "r": companion_r + _LABEL_GAP_PX + _LABEL_HALF_HEIGHT_PX * 2,
        "label_rect": _star_label_rect(companion_cx, companion_cy, companion_r, companion_attrs["name"]),
        "box": _star_box(companion_cx, companion_cy, companion_r),
    }

    primary_scene = _star_scene_svg(
        "system", f'System map for {system["name"]}', False,
        primary, primary_planets, primary_belts,
        _wide_binary_star_attrs(system, primary, is_primary=True),
        extra_svg=companion_marker_svg, extra_obstacle=companion_obstacle, facilities=facilities,
    )
    secondary_attrs = _wide_binary_star_attrs(system, secondary, is_primary=False)
    secondary_attrs["self"] = "true"
    secondary_scene = _star_scene_svg(
        secondary_scene_id, f'Planets of {_star_label(system, secondary)}', True,
        secondary, secondary_planets, secondary_belts, secondary_attrs, facilities=facilities,
    )
    return [primary_scene, secondary_scene]


def _render_system_scene(system, stars, planets, belts, facilities=None):
    """Builds the "whole system" scene -- see the module docstring for the
    barycenter/anchor model this uses for a close binary pair (a wide
    pair instead gets two separate scenes -- see
    `_render_wide_binary_scenes`)."""
    is_binary = len(stars) > 1
    # A 'close' (P-type) pair's planets orbit the merged effective proxy,
    # which sits at the shared barycenter -- `planets.star_id`/
    # `asteroid_belts.star_id` are both NULL for these (see schema.sql's
    # "v15" note), so `anchor_px` below never needs to distinguish them by
    # id at all, only by whether this system is a 'close' pair.
    close_binary = is_binary and system.get("binary_configuration") != "wide"

    if is_binary:
        star_pos_km = _binary_star_positions_km(system, stars)
    elif stars:
        star_pos_km = {stars[0]["id"]: (0.0, 0.0)}
    else:
        star_pos_km = {}

    local_r_list = [math.hypot(*star_pos_km[star["id"]]) for star in stars] if is_binary else []
    for planet in planets:
        local_r_list.append(_orbit_radius_km(planet))
    for belt in belts:
        local_r_list.extend((belt["lower_limit_km"], belt["upper_limit_km"]))

    lo, hi = _radial_scale_bounds(local_r_list)

    # Phase 1: place and relax the star(s) alone. A close binary's tiny
    # real separation (as little as ~0.05 AU) routinely log-compresses to
    # near `_MIN_RADIUS_PX` for both stars, so their drawn markers
    # frequently *do* overlap even though the real stars themselves never
    # do -- this settles that on its own, before anything else anchors to
    # the result.
    star_markers = []
    for index, star in enumerate(stars):
        sx, sy = star_pos_km[star["id"]]
        r_px = _radial_px(math.hypot(sx, sy), lo, hi)
        cx, cy = _polar_to_px(_CENTER_PX, _CENTER_PX, r_px, sx, sy)
        is_primary = index == 0
        star_markers.append({
            "type": "star", "cx": cx, "cy": cy, "r": _star_radius_px(star["radius_km"]),
            "star": star, "is_primary": is_primary,
        })
    _relax_markers(star_markers, _MARKER_GAP_PX)
    for marker in star_markers:
        marker["fixed"] = True

    star_anchor_px = {m["star"]["id"]: (m["cx"], m["cy"]) for m in star_markers}

    def anchor_px(star_id):
        if close_binary or not is_binary:
            return _CENTER_PX, _CENTER_PX
        return star_anchor_px.get(star_id, (_CENTER_PX, _CENTER_PX))

    # Phase 2: place every planet/belt against those now-final star
    # positions, then relax planets/moons against each other AND against
    # the (already-fixed) stars together -- a large star marker pushes a
    # too-close planet out of the way, but the star itself never moves
    # again here (see `_relax_markers`'s own `fixed` handling).
    orbit_paths = []
    # See _star_scene_svg's identical belt_svgs comment -- kept out of
    # orbit_paths (the aria-hidden orbits layer) since a belt ring is a
    # real focusable/clickable marker, not decorative.
    belt_svgs = []
    markers = list(star_markers)
    orbit_radii = {}  # anchor -> drawn orbit radii around it
    boxes = [_star_box(m["cx"], m["cy"], m["r"]) for m in star_markers]
    for planet in planets:
        ax_px, ay_px = anchor_px(planet.get("star_id"))
        lx_km, ly_km = planet.get("position_x_km") or 0.0, planet.get("position_y_km") or 0.0
        r_px = _radial_px(_orbit_radius_km(planet), lo, hi)
        orbit_radii.setdefault((ax_px, ay_px), []).append(r_px)
        boxes.append(_circle_box(ax_px, ay_px, r_px))
        cx, cy = _polar_to_px(ax_px, ay_px, r_px, lx_km, ly_km)
        orbit_paths.append(f'<circle class="sysmap-orbit" cx="{ax_px:.1f}" cy="{ay_px:.1f}" r="{r_px:.1f}"></circle>')
        markers.append({"type": "planet", "cx": cx, "cy": cy, "r": _planet_radius_px(planet["radius_km"]), "row": planet})

    facility_hosts = {}
    for belt in belts:
        ax_px, ay_px = anchor_px(belt.get("star_id"))
        r_px, band_px = _belt_band(belt, lo, hi, orbit_radii.get((ax_px, ay_px), ()))
        belt_svgs.append(_belt_ring_svg(ax_px, ay_px, r_px, band_px, belt))
        boxes.append(_circle_box(ax_px, ay_px, r_px + band_px / 2))
        facility_hosts[("asteroid_belt", belt["id"])] = (ax_px, ay_px, r_px, False)

    _relax_markers(markers, _MARKER_GAP_PX)
    # A close pair is orbited as one (`_db._facility_host`): a facility on
    # either star circles the barycenter, clear of both markers.
    pair_r = max((math.hypot(m["cx"] - _CENTER_PX, m["cy"] - _CENTER_PX) + m["r"] for m in star_markers),
                 default=0.0)
    for marker in markers:
        if marker["type"] == "star":
            host = (_CENTER_PX, _CENTER_PX, pair_r, True) if close_binary else (
                marker["cx"], marker["cy"], marker["r"], True)
            facility_hosts[("star", marker["star"]["id"])] = host
        else:
            facility_hosts[("planet", marker["row"]["id"])] = (marker["cx"], marker["cy"], marker["r"], False)
            boxes.append(_body_box(marker["cx"], marker["cy"], marker["r"], marker["row"]["body_type"]))
    facility_orbits, facility_markers, facility_boxes = _facilities_svg(facilities, facility_hosts, (lo, hi))
    bounds = _fit_bounds(boxes + facility_boxes)

    star_svgs = []
    planet_markers = []
    for marker in markers:
        if marker["type"] != "star":
            planet_markers.append(marker)
            continue
        star = marker["star"]
        is_primary = marker["is_primary"]
        sx, sy = star_pos_km[star["id"]]
        star_svgs.append(_star_marker_svg(marker["cx"], marker["cy"], marker["r"], star, {
            "kind": "star",
            "id": star["id"],
            "name": _star_label(system, star),
            "role": ("Primary" if is_primary else "Secondary") if is_binary else "Single",
            "type": star["star_type"], "temp": _star_temp_text(star),
            "mass": to_plain_text(format_star_mass(star["mass_kg"])),
            "radius": to_plain_text(format_star_radius(star["radius_km"])),
            "lum": to_plain_text(format_star_luminosity(star["luminosity_w"])),
            "radiuskm": star["radius_km"],
            # Real km (a close binary's own true, small separation from
            # the shared barycenter -- (0, 0) for a single star, per
            # `star_pos_km`'s own construction above), not this marker's
            # drawn `cx`/`cy` -- same reasoning as `_planet_attrs`'s
            # identical `xkm`/`ykm`.
            "xkm": sx, "ykm": sy,
        }, bounds))

    placements = _label_sides_2d([(m["cx"], m["cy"], m["r"], m["row"]["name"]) for m in planet_markers],
                                 bounds=bounds)
    body_svgs = []
    for marker, placement in zip(planet_markers, placements):
        row = marker["row"]
        scene_target = f'planet-{row["id"]}' if row.get("moons") else None
        attrs = _planet_attrs(row, kind="planet", scene_target=scene_target)
        body_svgs.append(_body_marker_svg(
            marker["cx"], marker["cy"], marker["r"], row["planet_class"], row["body_type"], row["name"],
            "sysmap-planet", attrs,
            label_direction=(placement["direction"] if placement else "below"),
            label_pushed=bool(placement and placement["pushed"]),
            label_shift=(placement["shift"] if placement else (0.0, 0.0)),
            show_label=(placement is not None),
            has_life=bool(row.get("life_chemical")),
            bounds=bounds,
        ))

    return _scene_svg_pair(
        "system", f'System map for {system["name"]}', False,
        "".join(orbit_paths) + facility_orbits,
        # belt_svgs drawn first, same relative "underneath" stacking its
        # old spot (inside orbit_paths, before star_svgs/body_svgs) had.
        f'{"".join(belt_svgs)}{"".join(star_svgs)}{"".join(body_svgs)}{facility_markers}',
        scale=(lo, hi), bounds=bounds,
    )


def _render_moon_scene(planet, facilities=None):
    moons = planet.get("moons") or []
    for moon in moons:
        moon["_is_moon"] = True
        moon["_parent_name"] = planet["name"]

    local_r_list = [
        math.hypot(moon.get("position_x_km") or 0.0, moon.get("position_y_km") or 0.0) or moon.get("distance_km") or 0.0
        for moon in moons
    ]
    lo, hi = _radial_scale_bounds(local_r_list)

    orbit_paths = []
    orbit_radii = []
    # The drilled-into planet itself acts as a fixed obstacle here (same
    # `fixed` treatment `_render_system_scene` gives its own star(s)) --
    # `_CENTER_PLANET_R` (36px) is large enough to otherwise swallow a
    # close-orbiting first moon outright once real angle/distance placement
    # puts it nearby.
    markers = [{"type": "center", "cx": _CENTER_PX, "cy": _CENTER_PX, "r": _CENTER_PLANET_R, "fixed": True}]
    for moon in moons:
        lx, ly = moon.get("position_x_km") or 0.0, moon.get("position_y_km") or 0.0
        r_px = _radial_px(math.hypot(lx, ly), lo, hi)
        orbit_radii.append(r_px)
        cx, cy = _polar_to_px(_CENTER_PX, _CENTER_PX, r_px, lx, ly)
        orbit_paths.append(f'<circle class="sysmap-orbit" cx="{_CENTER_PX:.1f}" cy="{_CENTER_PX:.1f}" r="{r_px:.1f}"></circle>')
        markers.append({"type": "moon", "cx": cx, "cy": cy, "r": _moon_radius_px(moon["radius_km"]), "row": moon})

    _relax_markers(markers, _MARKER_GAP_PX)
    markers = [m for m in markers if m["type"] == "moon"]
    facility_hosts = {("planet", planet["id"]): (_CENTER_PX, _CENTER_PX, _CENTER_PLANET_R, True)}
    boxes = [_body_box(_CENTER_PX, _CENTER_PX, _CENTER_PLANET_R, planet["body_type"], is_self=True)]
    boxes.extend(_circle_box(_CENTER_PX, _CENTER_PX, r_px) for r_px in orbit_radii)
    for marker in markers:
        facility_hosts[("moon", marker["row"]["id"])] = (marker["cx"], marker["cy"], marker["r"], False)
        boxes.append(_body_box(marker["cx"], marker["cy"], marker["r"], marker["row"]["body_type"]))
    facility_orbits, facility_markers, facility_boxes = _facilities_svg(facilities, facility_hosts, (lo, hi))
    bounds = _fit_bounds(boxes + facility_boxes)

    center_svg = _body_marker_svg(
        _CENTER_PX, _CENTER_PX, _CENTER_PLANET_R, planet["planet_class"], planet["body_type"], planet["name"],
        # The drilled-into planet is this scene's origin: its moons'
        # xkm/ykm are relative to it, so a measurement needs (0, 0) here,
        # not its position around the star.
        "sysmap-planet", dict(_planet_attrs(planet), xkm=0.0, ykm=0.0), is_self=True,
        has_life=bool(planet.get("life_chemical")), bounds=bounds,
    )

    placements = _label_sides_2d([(m["cx"], m["cy"], m["r"], m["row"]["name"]) for m in markers], bounds=bounds)
    body_svgs = []
    for marker, placement in zip(markers, placements):
        row = marker["row"]
        attrs = _planet_attrs(row, kind="moon", parent_name=row.get("_parent_name"))
        body_svgs.append(_body_marker_svg(
            marker["cx"], marker["cy"], marker["r"], row["planet_class"], row["body_type"], row["name"],
            "sysmap-moon", attrs,
            label_direction=(placement["direction"] if placement else "below"),
            label_pushed=bool(placement and placement["pushed"]),
            label_shift=(placement["shift"] if placement else (0.0, 0.0)),
            show_label=(placement is not None),
            has_life=bool(row.get("life_chemical")),
            bounds=bounds,
        ))

    # A moon scene starts hidden -- `static/systemmap.js` reveals it on
    # demand. This is a CSS class (`.sysmap-hidden`, toggled by
    # `classList`), not the HTML `hidden` attribute/`.hidden` IDL property:
    # `SVGElement` doesn't reflect that property to the attribute the way
    # `HTMLElement` does, so setting `.hidden` in JS silently only sets a
    # plain, attribute-less expando property on an `<svg>` -- confirmed
    # directly (`hasAttribute('hidden')` stayed false/true opposite of
    # what `.hidden` itself reported) -- leaving both this attribute and
    # any `[hidden]` CSS rule permanently out of sync with it.
    return _scene_svg_pair(
        f'planet-{planet["id"]}', f'Moons of {planet["name"]}', True,
        "".join(orbit_paths) + facility_orbits,
        f'{center_svg}{"".join(body_svgs)}{facility_markers}',
        scale=(lo, hi), bounds=bounds,
    )


# MAP.57: a NaN or infinite number from the database never reaches the
# SVG. Every row is copied with such numbers turned into `None` (which
# every formatter and placement step already treats as "not recorded");
# a planet or moon whose position is lost that way but whose orbit
# distance survives is drawn due east at that distance with a note in its
# info panel, and one with neither is left out.
_POSITION_KEYS = ("position_x_km", "position_y_km", "position_z_km")
_POSITION_NOTE = "Position not recorded: drawn at its orbit distance, due east of what it orbits."


def _is_bad_number(value):
    return isinstance(value, float) and not math.isfinite(value)


def _finite_row(row):
    """A copy of `row` with every NaN or infinite float set to `None`."""
    return {key: None if _is_bad_number(value) else value for key, value in row.items()}


def _finite_body(row):
    """`_finite_row` for a planet or moon (and its moons), or `None` when
    it has no usable position or distance left to draw it at."""
    lost = any(_is_bad_number(row.get(key)) for key in _POSITION_KEYS)
    body = _finite_row(row)
    if lost:
        distance = body.get("distance_km")
        if not distance or distance <= 0:
            return None
        body.update(position_x_km=distance, position_y_km=0.0, position_z_km=0.0, _position_note=_POSITION_NOTE)
    if body.get("moons"):
        body["moons"] = [moon for moon in map(_finite_body, body["moons"]) if moon is not None]
    return body


_NEBULA_WASH_TYPES = ("diffuse", "emission", "reflection", "planetary", "dark")
"""tuple: The nebula types `style.css` gives a wash color
(`.sysmap-nebula-<type>`); any other gets `.sysmap-nebula-other`."""


def render_system_map_panel(system, stars, planets, belts, facilities=None, scene_url=None, nav=None):
    """
    Builds the "System Map" panel embedded in `system.py`: TWO sibling
    `<svg>`s per scene (an orbits-only layer plus a body-marker layer --
    see `_scene_svg_pair`'s own docstring for why two, not one), for the
    whole system plus one more pair for every planet with moons (see the
    module docstring for why this doesn't need to recurse any deeper than
    that), and an info side panel that `static/systemmap.js` fills in on
    click and swaps between scenes.

    Args:
        system (dict): The `star_systems` row (and, when it sits inside a
                      nebula, `inside`: `{"type": "nebula", "name",
                      "descriptor"}` -- the map is tinted its color),
                      including
                      `binary_configuration` (`'close'`, `'wide'`, or
                      `None` -- see `schema.sql`'s "v15" note) and
                      `binary_mutual_position_x/y/z_km` (the secondary's
                      position relative to the primary -- see this
                      module's own docstring for how a 'close' vs. 'wide'
                      pair each use it).
        stars (list[dict]): 1 entry (single star) or 2 (primary, then
                            secondary), each with `id`, `mass_kg`,
                            `star_type`, `temperature_k`, `radius_km`,
                            `luminosity_w`.
        planets (list[dict]): Each a `planets` row (including its own
                              `star_id` -- NULL for a 'close' pair, the
                              specific owning star for a single star or a
                              'wide' pair -- and `position_x/y/z_km`) plus
                              a `moons` key (list of `moons` rows, possibly
                              empty).
        belts (list[dict]): `asteroid_belts` rows, including `star_id`.
        facilities (list[dict], optional): `GET /api/systems/<id>/
            facilities` items, drawn at their hosts (`_facilities_svg`).
        nav (dict, optional): The NAV buttons of a body's info panel
            (NAV.50, `system_pages.body_nav`): `take` and `label` while a
            NAV start or destination is being picked, else `start` and
            `end`; each a URL with `{ref}` where the body's object
            reference goes. Handed to `static/systemnav.js` as `data-nav-*`.
        scene_url (str, optional): Where the 3D view fetches the system's
            scene JSON (`GET /api/systems/<id>/scene`, relayed by the
            site). Without it the panel is the diagram only, with no 3D
            switch (MAP.74).

    Returns:
        str: A complete `<section class="panel">` block.
    """
    system = _finite_row(system)
    stars = [_finite_row(star) for star in stars]
    planets = [planet for planet in map(_finite_body, planets) if planet is not None]
    belts = [_finite_row(belt) for belt in belts]
    facilities = [_finite_row(facility) for facility in facilities] if facilities else facilities
    is_wide_binary = len(stars) > 1 and system.get("binary_configuration") == "wide"
    if is_wide_binary:
        scenes = _render_wide_binary_scenes(system, stars, planets, belts, facilities)
    else:
        scenes = [_render_system_scene(system, stars, planets, belts, facilities)]
    scenes.extend(
        _render_moon_scene(planet, facilities) for planet in planets if planet.get("moons")
    )

    if planets or belts or stars:
        info_panel = (
            '<aside class="starmap-info" id="sysmap-info">'
            '<p class="hint">Click a star, planet, moon, asteroid belt or facility for details.</p></aside>'
        )
    else:
        info_panel = '<aside class="starmap-info" id="sysmap-info"><p class="hint">Nothing to show yet.</p></aside>'

    # One shared WebGL canvas, sized to cover the whole `.sysmap-viewport`
    # and drawn BETWEEN each scene's two sibling `<svg>`s -- below the
    # body-marker layer, above the orbits-only layer (see `style.css`'s
    # `.sysmap-spheres-canvas`/`.sysmap-svg`/`.sysmap-orbits-layer`
    # z-index rules, and `_scene_svg_pair`'s own docstring for why each
    # scene is two `<svg>`s rather than one) -- that `static/systemmap.js`
    # uses to render every visible star/planet/moon marker's own live 3D
    # sphere in place, each scissored to exactly that marker's own
    # `<circle>`, so a sphere now actually occludes the orbit line drawn
    # under it instead of the line always painting on top regardless of
    # the sphere's own footprint. A single `<canvas>`, not one per marker
    # (or an `<svg>`, the one part of the System Map that's genuinely 3D).
    # `stars` (never empty whenever this panel renders at all -- see
    # `render_system_map_panel`'s caller) is this block's presence guard,
    # same reasoning the old single-body preview used `planets` for.
    spheres_html = (
        '<canvas id="sysmap-spheres-canvas" class="sysmap-spheres-canvas" aria-hidden="true"></canvas>'
        if stars else ""
    )

    # "Measure distance" toggle -- lets a visitor click any two bodies in
    # the currently visible scene (two planets, two moons, or a binary's
    # two stars) for the real straight-line distance between them, plus a
    # route-around-the-obstacle distance when that straight line would
    # otherwise pass through the scene's own center body (see
    # `static/systemmap.js`'s own `computeMeasurement`). Gated on `stars`
    # the same way `spheres_html` above is -- a system with no stars at
    # all (shouldn't happen -- see this function's own docstring) has
    # nothing to measure between regardless.
    # UX.61, UX.50: under the map, as on the other maps, beside Map help.
    controls_html = (
        '<div class="starmap-controls" id="sysmap-controls">'
        '<button type="button" class="starmap-btn" id="sysmap-measure-btn" aria-pressed="false">'
        "Measure distance</button>"
        '<button type="button" class="starmap-btn" data-dialog-open="sysmap-help" '
        'title="How to read the map">Map help</button></div>'
        if stars else ""
    )

    view_switch = ""
    view3d_html = ""
    if scene_url and stars:
        view_switch = (
            '<div class="sysview-switch" role="group" aria-label="Map view">'
            '<button type="button" class="starmap-btn" id="sysmap-view-diagram" aria-pressed="true">Diagram</button>'
            '<button type="button" class="starmap-btn" id="sysmap-view-3d" aria-pressed="false">3D</button></div>'
        )
        view3d_html = (
            f'<div class="sysview3d" id="sysview3d" data-scene-url="{esc(scene_url)}" hidden>'
            '<div class="starmap-viewport sysview-viewport">'
            '<canvas id="sysview3d-canvas" class="sysview-canvas" aria-label="3D system view" role="img"></canvas>'
            '<div class="sysview-labels" id="sysview3d-labels" aria-hidden="true"></div>'
            '<div class="map-tooltip" id="sysview3d-tip" hidden></div>'
            '</div>'
            '<div class="starmap-controls sysview-controls" id="sysview3d-controls">'
            '<button type="button" class="starmap-btn" id="sysview3d-play" aria-pressed="true">Pause</button>'
            '<button type="button" class="starmap-btn" id="sysview3d-slower">Slower</button>'
            '<button type="button" class="starmap-btn" id="sysview3d-faster">Faster</button>'
            '<button type="button" class="starmap-btn" id="sysview3d-now">Now</button>'
            '<output class="sysview-rate" id="sysview3d-rate" aria-live="polite"></output>'
            '<label class="sysview-scale">Scale <select id="sysview3d-scale">'
            '<option value="compressed">Compressed</option><option value="true">True scale</option></select></label>'
            '<button type="button" class="starmap-btn" id="sysview3d-recenter" '
            'title="Center the view on the selected body, or the middle (C)">Center</button>'
            '<button type="button" class="starmap-btn" id="sysview3d-reset">Reset view</button>'
            '<button type="button" class="starmap-btn" data-dialog-open="sysview3d-help" '
            'title="How to use the 3D view">Map help</button></div>'
            '<p class="hint sysview-scale-note" id="sysview3d-scale-note"></p>'
            '<details class="sysview-list"><summary>Bodies in this view</summary>'
            '<ul id="sysview3d-list"></ul></details></div>'
        )
        view3d_help = (
            '<sl-dialog id="sysview3d-help" class="map-help-dialog" label="3D system view help"><ul>'
            "<li>Drag to turn, shift-drag or right-drag to pan, the wheel or a pinch to zoom.</li>"
            "<li>W, A, S, D and the arrow keys fly; C or Center puts the center the view turns about on the selected body; Home or Reset view goes back.</li>"
            "<li>Click a body for details; double-click flies to it. Follow rides along with it.</li>"
            "<li>Pause, Slower and Faster move time; Now returns to the present.</li>"
            "<li>Compressed scale shows everything; True scale shows real sizes and distances.</li></ul>"
            '<sl-button slot="footer" data-dialog-close>Close</sl-button></sl-dialog>'
        )
    else:
        view3d_help = ""

    facility_legend = (
        '<span class="sysmap-legend-facility" aria-hidden="true"></span> marks a facility;' if facilities else ""
    )

    # A system inside a nebula sits in a faint wash of its color (MAP.103):
    # at this scale (AU) the cloud, light-years across, is only the sky.
    inside = system.get("inside")
    in_nebula = bool(inside) and inside.get("type") == "nebula"
    wash_class = ""
    wash_hint = ""
    if in_nebula:
        descriptor = inside.get("descriptor")
        wash_class = " sysmap-in-nebula sysmap-nebula-" + (descriptor if descriptor in _NEBULA_WASH_TYPES else "other")
        wash_hint = f". The tint is the gas of {esc(inside.get('name') or 'a nebula')} around the system"

    nav_attrs = "".join(
        f' data-nav-{key}="{esc(value)}"' for key, value in (nav or {}).items() if value)

    return f"""
<section class="panel">
<div class="panel-header">
  <h2 class="sr-only">System Map</h2>
</div>
<div class="starmap-layout sysmap-layout" id="sysmap-root"{nav_attrs}>
<div class="sysmap-map-col">
{view_switch}
<div class="sysmap-diagram" id="sysmap-diagram">
<div class="starmap-viewport sysmap-viewport{wash_class}">
{spheres_html}
{''.join(scenes)}
</div>
{controls_html}
</div>
{view3d_html}
</div>
<div class="starmap-side">
<div class="sysmap-crumb" id="sysmap-crumb"></div>
{info_panel}
</div>
</div>
<sl-dialog id="sysmap-help" class="map-help-dialog" label="System Map help">
  <ul>
    <li>True top-down positions: the real angle, with distance on a log scale.</li>
    <li>Click a planet with moons to view its moon system.</li>
    <li>Circle size is the body's radius (log scale); color is the planet class;
    <span class="sysmap-legend-life-badge" aria-hidden="true"></span> marks a world that supports life;
    {facility_legend}</li>
    <li>Measure distance, then click two bodies for the real distance between them{wash_hint}.</li>
  </ul>
  <sl-button slot="footer" data-dialog-close>Close</sl-button>
</sl-dialog>
{view3d_help}
</section>
"""
