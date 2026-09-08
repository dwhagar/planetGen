# html/lib/systemmap.py

"""
Interactive 2D "System Map" for `system.py`: a letterboxed concentric-
orbit diagram. Every orbit (a planet's or a belt's distance from the
star) is a real circle centered on the star -- exactly the classical top-
down solar-system picture -- but the star sits at the far-left edge of a
short, wide, fixed-size window, so only the sliver of each circle nearest
the star is ever visible (see `_orbit_arc_points`): a close-in orbit's
tight circle still shows real curvature within that sliver, while a far
orbit's much larger circle reads as almost a flat vertical line, the same
way a small coin's edge looks more curved up close than a car tyre's.
Each body sits at the apex of its own orbit's visible arc (the single
point closest to the star, directly "ahead" of it) -- since every arc's
apex sits on the same horizontal line through the star, and successively
larger orbits place that apex further right, the net effect is exactly
"planets lined up left to right in distance order" with a nested-bracket
orbit indicator behind each one. The whole diagram is a fixed size (no
scrolling or zooming): each orbit's distance is log-scaled into a fixed
pixel budget shared by every body in the scene, uniformly compressed
further if that scaling alone would still overflow it (see
`_orbit_radii_px`) -- a crowded system reads as tightly packed rather
than as a diagram that needs to scroll.

Each body is colored by planet class and sized (log-scaled, so a Class D
moonlet and a Class J gas giant both stay visible on the same diagram) by
`radius_km`.

Unlike `starmap.py`'s draggable 3D cube (positions come straight from
stored (x, y, z) coordinates -- there's real 3D data to project), a
system has no meaningful 2D layout data of its own: only a distance from
the star and an orbital index. So instead of a rotatable scene, this
renders one flat `<svg>` per "scene" -- the whole system, plus one more
for every planet that has moons -- and `static/systemmap.js` just toggles
which scene is visible. Clicking a planet with moons swaps the view so
that planet sits where the star was, with its own moons lined up the same
way after it (mirroring a real "zoom into this planet's moon system"
diagram); clicking any other object or a belt instead fills the info side
panel from its `data-*` attributes, exactly like `starmap.py`/
`sectormap.js`'s own click-for-info pattern.
"""

import colorsys
import math
import statistics

from dbutil import esc
from starmap import _star_color, _SUN_RADIUS_KM
from tabledisplay import (
    format_body_distance, format_period, format_star_luminosity, format_star_mass, format_star_radius,
)

try:
    from stellarObjects.program_constants import PLANET_CLASSES
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # the map still works, just without each class's one-line flavor
    # description in the info panel.
    PLANET_CLASSES = {}

# The whole diagram's fixed coordinate space -- short and wide (see the
# module docstring for why), and never scrolled or zoomed, so every
# system's worth of orbits has to fit inside it regardless of how many
# there are (`_orbit_radii_px` is what guarantees that).
_VIEW_W_PX = 1000.0
_VIEW_H_PX = 300.0
_STAR_X_PX = 60.0
_STAR_Y_PX = _VIEW_H_PX / 2

# How far right each orbit's apex (see the module docstring) is allowed
# to sit, in the same rank+log-blended units `_orbit_radii_px` computes
# before rescaling -- `_ORBIT_SPREAD_PX` is the *target* spread when
# there's room for it; `_MAX_ORBIT_R_PX` is the hard budget the whole
# scene gets rescaled down to fit within when there isn't (many bodies,
# e.g. `+max_planets +asteroid_belt` forcing 30+ orbital slots).
_MIN_ORBIT_GAP_PX = 55.0
_ORBIT_SPREAD_PX = 620.0
_MAX_ORBIT_R_PX = 880.0

# Vertical half-extent of the sliver of each orbit circle that's actually
# drawn -- see `_orbit_arc_points`. An orbit whose own radius is smaller
# than this is small enough that its *entire* circle fits in that
# vertical band, so the visible "sliver" is really the whole circle.
_ARC_HALF_HEIGHT_PX = 62.0

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

# Minimum clearance (beyond simple non-overlap) between the center body's
# own edge and the innermost orbiting body's edge -- `_MIN_ORBIT_GAP_PX`
# alone assumes both are small, but the center can be a giant/supergiant
# star (`_STAR_MAX_R` = 48) or the fixed-size drilled-into planet
# (`_CENTER_PLANET_R` = 36), either of which can otherwise swallow a
# close, large-radius first body outright -- confirmed directly: a
# supergiant star's own circle fully overlapped a nearby Class A planet's
# marker before `_render_row_svg` started enforcing this.
_CENTER_CLEARANCE_PADDING_PX = 8.0

# Secondary-star offset/size-cap for a binary system's own two dots,
# mirroring `starmap.py`'s `_BINARY_OFFSET_FRACTION`/`_SECONDARY_MAX_RATIO`
# (see that module's own comments for the reasoning) -- duplicated rather
# than imported since these are laid out in flat 2D pixels here, not
# `starmap.py`'s local 3D scene-coordinate space.
_BINARY_OFFSET_FRACTION = 0.9
_SECONDARY_MAX_RATIO = 0.65

_ZONE_LABELS = {"h": "Hot zone", "e": "Ecosphere", "c": "Cold zone"}

# An asteroid belt is a *range* of orbit radii rather than one -- drawn as
# a shaded band between its own lower/upper limits (see `_belt_band_svg`)
# instead of one more single-radius arc, clamped so a belt with almost no
# radial spread is still visibly a band and one with an enormous spread
# doesn't swallow its neighbors.
_BELT_MIN_BAND_PX = 9.0
_BELT_MAX_BAND_PX = 40.0

# Anchor color per planet class, hand-picked (not derived from any
# physical model, unlike `starmap.py`'s spectral colors -- there's no
# single physical quantity a planet class maps to the way a star's
# spectral letter maps to a blackbody hue) to read as distinct from its
# neighbors while still gesturing at the class's own flavor: A/B/E/K's
# reds-and-rust for volcanic/molten worlds, C/D's muted greys for dead/icy
# small bodies, F/G/H/L's earthy greens-and-tans for the barren-to-
# vegetated progression, M/O/P's blues for the habitable "water world"
# spread (M dark, O a brighter ocean blue, P pale glacier blue), I/J/T's
# gas-giant palette, and N/Q/V/W's outliers (Venus-hot mustard,
# eccentric/tidally-locked violets, a high-gravity super-Earth magenta).
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
    "T": "#7a86d9",
    "V": "#9c3f63",
    "W": "#c06bb0",
}
_DEFAULT_CLASS_COLOR = "#8a8f9c"


def _class_color(planet_class):
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


def _orbit_radii_px(distances_km):
    """
    Returns one orbit radius (pixels, distance from the star) per entry
    in `distances_km`, in the same order, assuming the list is already
    sorted ascending by distance.

    Each body's *ideal* radius is log-scaled between a floor and
    `_ORBIT_SPREAD_PX` (so relative spacing tracks relative distance -- a
    log scale compresses the huge inner/outer ratio a real system spans
    far less harshly than a sqrt scale would, while still keeping every
    ratio finite). That floor is *not* the innermost body's own distance --
    doing that would anchor it to frac=0 by construction, which is what
    used to put the first orbit right against the star regardless of how
    far out it actually was, and (with few bodies) blow the entire pixel
    budget on whatever ratio happened to separate the two extremes even
    when neighboring bodies were genuinely close together. Instead the
    floor sits one *typical* log-step below the innermost body -- the
    median gap between this same system's own consecutive bodies -- so the
    star-to-first-orbit gap reads on the same scale as every other gap
    instead of always collapsing to zero (confirmed directly: a two-planet
    system 2.8x apart in real distance used to place them at the two
    extreme edges of the frame, indistinguishable from a pair 1000x apart).
    Those ideal radii would still let a tight cluster of close-in bodies
    collapse into a few pixels of each other, so
    `_resolve_min_gap` nudges only the bodies that actually collide apart
    -- symmetrically, around their own shared average, not by ratcheting
    every later body rightward off of one early collision -- rather than
    applying a blanket minimum gap that would flatten real log-scale
    separation into uniform rank spacing almost everywhere (confirmed
    directly: with a flat per-step clamp instead of this, most consecutive
    gaps in a real generated system came out exactly equal to the clamp,
    despite genuinely different distances). Finally -- since this diagram
    never scrolls or zooms -- the whole resolved layout is uniformly
    rescaled down if it still doesn't fit within `_MAX_ORBIT_R_PX`, so an
    unusually crowded system (many bodies forced via `+max_planets`)
    always fits the fixed frame, just more tightly packed rather than
    needing more room.
    """
    valid = [d for d in distances_km if d and d > 0]
    if len(valid) >= 2:
        log_vals = sorted(math.log10(d) for d in valid)
        typical_step = statistics.median(b - a for a, b in zip(log_vals, log_vals[1:]))
        lo = 10 ** (log_vals[0] - typical_step)
        hi = 10 ** log_vals[-1]
    elif valid:
        lo = hi = valid[0]
    else:
        lo, hi = 0.0, 0.0
    ideal = []
    for index, distance_km in enumerate(distances_km):
        if hi > lo and distance_km and distance_km > 0:
            distance_km = max(lo, min(hi, distance_km))
            frac = (math.log10(distance_km) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
        elif len(distances_km) > 1:
            frac = index / (len(distances_km) - 1)
        else:
            frac = 0.0
        ideal.append(_MIN_ORBIT_GAP_PX + frac * _ORBIT_SPREAD_PX)

    radii = _resolve_min_gap(ideal, _MIN_ORBIT_GAP_PX)

    # A cluster's resolved center can pull its own leftmost member closer
    # to the star than `_MIN_ORBIT_GAP_PX` (e.g. two bodies at nearly
    # identical distances, both belonging to the same cluster) -- shift
    # everything uniformly right to restore that clearance rather than
    # letting the innermost marker crowd the star itself, same as every
    # other gap in this list, this preserves every relative spacing.
    if radii and radii[0] < _MIN_ORBIT_GAP_PX:
        shift = _MIN_ORBIT_GAP_PX - radii[0]
        radii = [r + shift for r in radii]

    if radii and radii[-1] > _MAX_ORBIT_R_PX:
        shrink = _MAX_ORBIT_R_PX / radii[-1]
        radii = [r * shrink for r in radii]
    return radii


def _resolve_min_gap(ideal_positions, min_gap):
    """
    Given `ideal_positions` (already sorted ascending), returns positions
    in the same order that (a) never sit closer than `min_gap` to their
    neighbor and (b) stay as close as possible to their own ideal
    position, in the least-squares sense -- the standard "pool adjacent
    violators" approach to isotonic-with-minimum-spacing regression, also
    used for de-overlapping a column of sorted chart labels.

    Each maximal run of positions that collide once spaced `min_gap` apart
    (a "cluster") is re-centered on the *average* of its own members' own
    ideal positions, then laid out evenly `min_gap` apart around that
    center -- so within a cluster, members below the average shift right
    and members above it shift left, both by as little as the constraint
    allows, rather than one early collision permanently displacing every
    later position (what a simple left-to-right "at least min_gap past
    the previous point" clamp does instead).
    """
    # Each cluster is (sum of its members' own ideal positions, member
    # count) -- from which its center-of-mass anchor (sum / count) and
    # its evenly-`min_gap`-spaced span around that anchor are derived on
    # demand, both while merging below and when expanding back out at the
    # end. Processed as a stack: only ever the top one or two clusters
    # can be in violation of the gap constraint at any point, since
    # everything below the top was already fully resolved against its own
    # neighbors on a previous iteration.
    clusters = []
    for position in ideal_positions:
        clusters.append((position, 1))
        while len(clusters) > 1:
            (sum1, n1), (sum2, n2) = clusters[-2], clusters[-1]
            right_edge_1 = sum1 / n1 + (n1 - 1) * min_gap / 2
            left_edge_2 = sum2 / n2 - (n2 - 1) * min_gap / 2
            if left_edge_2 - right_edge_1 >= min_gap - 1e-9:
                break
            clusters[-2:] = [(sum1 + sum2, n1 + n2)]

    resolved = []
    for total, count in clusters:
        center = total / count
        start = center - (count - 1) * min_gap / 2
        resolved.extend(start + i * min_gap for i in range(count))
    return resolved


def _theta_max_for_radius(radius_px):
    """The half-angle (radians, either side of due-"east") of the visible
    sliver of an orbit circle of `radius_px` -- see the module docstring.
    A circle smaller than `_ARC_HALF_HEIGHT_PX` fits entirely within the
    vertical band regardless of angle, so its own full right half
    (+/-90 degrees) is what's visible."""
    return math.asin(min(1.0, _ARC_HALF_HEIGHT_PX / radius_px)) if radius_px > 0 else math.pi / 2


def _orbit_arc_points(star_x, star_y, radius_px, samples=22):
    theta_max = _theta_max_for_radius(radius_px)
    return [
        (
            star_x + radius_px * math.cos(-theta_max + 2 * theta_max * i / samples),
            star_y + radius_px * math.sin(-theta_max + 2 * theta_max * i / samples),
        )
        for i in range(samples + 1)
    ]


def _orbit_arc_svg(star_x, star_y, radius_px):
    points = _orbit_arc_points(star_x, star_y, radius_px)
    d = f"M {points[0][0]:.1f} {points[0][1]:.1f} " + " ".join(f"L {x:.1f} {y:.1f}" for x, y in points[1:])
    return f'<path class="sysmap-orbit" d="{d}"></path>'


def _belt_band_svg(star_x, star_y, r_lower_px, r_upper_px, attrs, label_text):
    """A belt is a *range* of orbit radii -- drawn as a shaded band
    between the visible slivers of its lower and upper bounds (both
    slivers built at the outer radius's own angular range, so the two
    edges stay a consistent height apart rather than the inner edge
    narrowing on its own tighter angle) rather than one more single-point
    marker."""
    samples = 18
    theta_max = _theta_max_for_radius(r_upper_px)
    outer = [
        (
            star_x + r_upper_px * math.cos(-theta_max + 2 * theta_max * i / samples),
            star_y + r_upper_px * math.sin(-theta_max + 2 * theta_max * i / samples),
        )
        for i in range(samples + 1)
    ]
    inner = [
        (
            star_x + r_lower_px * math.cos(-theta_max + 2 * theta_max * i / samples),
            star_y + r_lower_px * math.sin(-theta_max + 2 * theta_max * i / samples),
        )
        for i in range(samples + 1)
    ]
    d = (
        f"M {outer[0][0]:.1f} {outer[0][1]:.1f} "
        + " ".join(f"L {x:.1f} {y:.1f}" for x, y in outer[1:])
        + " " + " ".join(f"L {x:.1f} {y:.1f}" for x, y in reversed(inner))
        + " Z"
    )
    return (
        f'<g class="sysmap-belt" tabindex="0" role="button"{_data_attrs(attrs)} aria-label="{esc(label_text)}">'
        f'<path d="{d}"></path>'
        "</g>"
    )


def _data_attrs(attrs):
    return "".join(f' data-{key}="{esc(value)}"' for key, value in attrs.items() if value is not None)


def _body_marker_svg(cx, cy, r_px, planet_class, body_type, label_text, extra_class, attrs, is_self=False,
                      label_above=False, show_label=True):
    """
    Builds one clickable `<g>` for a planet or moon: a filled/stroked
    circle colored by `planet_class` (see `_CLASS_COLORS`), the class
    letter centered inside it once the circle is big enough to hold text
    legibly, a small tilted ring behind gas giants (`body_type == "g"`)
    for an at-a-glance silhouette cue beyond just color, and (when
    `show_label` is set) a name label underneath -- or, when `label_above`
    is also set (see `_label_sides`), above, to keep it clear of a
    tightly-packed neighbor's own label. `show_label` alone going False
    (also from `_label_sides`, when even alternating sides couldn't clear
    a tight enough cluster) only drops the visible text, not the body's
    name from `aria-label` -- it's still reachable, just via a click on
    the marker rather than a glance at the diagram. `is_self` marks the
    one body a moon-scene is *about* (the planet drilled into, now
    standing in for the scene's star) -- drawn with a soft halo and picked
    out by `static/systemmap.js` as the default info-panel content when
    that scene opens, via its `data-self="true"` marker.
    """
    fill = _class_color(planet_class)
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
    parts.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r_px:.1f}" fill="{fill}" stroke="{stroke}"></circle>')
    if planet_class and r_px >= 8:
        parts.append(
            f'<text class="sysmap-classletter" x="{cx:.1f}" y="{cy:.1f}" fill="{text_color}" '
            f'text-anchor="middle" dominant-baseline="central">{esc(planet_class.upper())}</text>'
        )
    if label_text and show_label:
        label_y = cy - r_px - 8 if label_above else cy + r_px + 14
        parts.append(
            f'<text class="sysmap-label" x="{cx:.1f}" y="{label_y:.1f}" text-anchor="middle">{esc(label_text)}</text>'
        )

    self_attr = ' data-self="true"' if is_self else ""
    return (
        f'<g class="sysmap-body {extra_class}" tabindex="0" role="button"{self_attr}{_data_attrs(attrs)} '
        f'aria-label="{esc(label_text or "body")}">{"".join(parts)}</g>'
    )


def _star_marker_svg(cx, cy, r_px, star, attrs):
    fill, stroke = _star_color(star["star_type"], star["temperature_k"], star["luminosity_w"])
    name = attrs.get("name", "star")
    return (
        f'<g class="sysmap-body sysmap-star" tabindex="0" role="button"{_data_attrs(attrs)} '
        f'aria-label="{esc(name)}">'
        f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r_px:.1f}" fill="{fill}" stroke="{stroke}"></circle>'
        f'<text class="sysmap-label sysmap-star-label" x="{cx:.1f}" y="{cy + r_px + 15:.1f}" '
        f'text-anchor="middle">{esc(name)}</text>'
        "</g>"
    )


def _stars_svg(cx, cy, system, stars):
    """Returns `(svg, primary_radius_px)` -- the caller needs the
    primary's own drawn radius separately, to size the clearance gap to
    the first orbiting body (see `_render_row_svg`)."""
    if not stars:
        return "", 0.0
    primary = stars[0]
    is_binary = len(stars) > 1
    primary_r = _star_radius_px(primary["radius_km"])
    parts = [_star_marker_svg(cx, cy, primary_r, primary, {
        "kind": "star",
        "name": f'{system["name"]}{" -- Primary" if is_binary else ""}',
        "role": "Primary" if is_binary else "Single",
        "type": primary["star_type"],
        "temp": f'{int(primary["temperature_k"])} K',
        "mass": format_star_mass(primary["mass_kg"]),
        "radius": format_star_radius(primary["radius_km"]),
        "lum": format_star_luminosity(primary["luminosity_w"]),
    })]
    if is_binary:
        secondary = stars[1]
        offset = primary_r * _BINARY_OFFSET_FRACTION
        secondary_r = min(_star_radius_px(secondary["radius_km"]), primary_r * _SECONDARY_MAX_RATIO)
        parts.append(_star_marker_svg(cx - offset * 0.3, cy + offset, secondary_r, secondary, {
            "kind": "star",
            "name": f'{system["name"]} -- Secondary',
            "role": "Secondary",
            "type": secondary["star_type"],
            "temp": f'{int(secondary["temperature_k"])} K',
            "mass": format_star_mass(secondary["mass_kg"]),
            "radius": format_star_radius(secondary["radius_km"]),
            "lum": format_star_luminosity(secondary["luminosity_w"]),
        }))
    return "".join(parts), primary_r


def _planet_attrs(planet, kind="planet", parent_name=None, scene_target=None):
    attrs = {
        "kind": kind,
        "id": planet["id"],
        "name": planet["name"],
        "class": (planet["planet_class"] or "").upper(),
        "classdesc": _class_description(planet["planet_class"]),
        "bodytype": "Gas Giant" if planet["body_type"] == "g" else "Terrestrial",
        "zone": _ZONE_LABELS.get(planet["zone"], ""),
        "distance": format_body_distance(planet["distance_km"], planet.get("_is_moon", False)),
        "period": format_period(planet["period_years"]),
        "gravity": f'{round(planet["gravity_g"], 3) if planet["gravity_g"] is not None else ""} g',
    }
    if parent_name is not None:
        attrs["parent"] = parent_name
    if scene_target is not None:
        attrs["scene"] = scene_target
        attrs["moons"] = len(planet.get("moons") or [])
    return attrs


def _belt_band_px(belt, radius_px):
    spread_fraction = (belt["upper_limit_km"] - belt["lower_limit_km"]) / max(belt["distance_km"], 1.0)
    return max(_BELT_MIN_BAND_PX, min(_BELT_MAX_BAND_PX, radius_px * spread_fraction))


def _first_body_extent_px(kind, row, radius_px):
    """The first orbiting body's own "radius" for clearance purposes --
    the drawn marker radius for a planet/moon, or half the drawn band
    width for a belt (its near edge sits that far inside its own nominal
    orbit radius). Needed by `_render_row_svg` to keep the center body
    from overlapping whatever sits on the innermost orbit."""
    if kind == "belt":
        return _belt_band_px(row, radius_px) / 2
    radius_fn = _moon_radius_px if row.get("_is_moon") else _planet_radius_px
    return radius_fn(row["radius_km"])


_LABEL_CHAR_WIDTH_PX = 8.5
_LABEL_MIN_HALFWIDTH_PX = 20.0
_LABEL_GAP_PAD_PX = 6.0


def _label_half_width_px(text):
    """A cheap stand-in for actually measuring `text` at `.sysmap-label`'s
    font size (not available server-side, since this is a static SVG) --
    just enough to decide whether two neighboring names would collide, not
    to lay out precisely."""
    return max(_LABEL_MIN_HALFWIDTH_PX, len(text or "") * _LABEL_CHAR_WIDTH_PX / 2)


def _label_sides(cx_and_names):
    """
    Given `[(cx, name), ...]` in ascending-cx order, returns one
    `"below"`/`"above"`/`None` per entry: which band that body's name
    label should be drawn in, or `None` to skip drawing it at all.
    `_orbit_radii_px`'s own collision resolution only keeps *markers* from
    overlapping (`_MIN_ORBIT_GAP_PX`, far narrower than most names render
    at) -- a run of tightly-packed bodies would otherwise stack their
    labels into an unreadable smear along the one horizontal band every
    default "below" label shares. Each label is checked against the
    nearest earlier label still in the same band (a "below" band freed up
    by the previous label going "above" is fair game again, and vice
    versa), so a moderately tight run alternates below/above/below/... --
    the standard fix for crowded labels along one axis. A run tight enough
    to still collide in *both* bands (three or more names packed within
    barely `_MIN_ORBIT_GAP_PX` of each other) drops the label rather than
    drawing overlapping text -- that body's name, and everything else
    about it, is still one click away in the info panel.
    """
    sides = []
    prev_edge = {"below": None, "above": None}
    for cx, name in cx_and_names:
        half = _label_half_width_px(name)
        fits = {
            band: prev_edge[band] is None or cx - half - prev_edge[band] >= _LABEL_GAP_PAD_PX
            for band in ("below", "above")
        }
        side = "below" if fits["below"] else "above" if fits["above"] else None
        if side is not None:
            prev_edge[side] = cx + half
        sides.append(side)
    return sides


def _render_row_svg(scene_id, aria_label, star_x, star_y, star_svg, orbit_entries, center_radius_px):
    """
    Shared layout engine for both scene kinds (the whole system, and one
    planet's moons): places `orbit_entries` -- `(kind, row)` pairs,
    "planet" or "belt", already sorted ascending by distance -- each on
    its own concentric orbit circle centered at `(star_x, star_y)`, and
    returns the complete, fixed-size `<svg>`. `star_svg` is pre-built
    markup for whatever sits at the center (the system's own star(s), or
    a drilled-into planet standing in for one) -- this function only
    positions what orbits it. `center_radius_px` is that center body's
    own drawn radius (the primary star's, or `_CENTER_PLANET_R`), used
    below to keep it clear of the innermost orbiting body.
    """
    radii = _orbit_radii_px([row["distance_km"] for _kind, row in orbit_entries])

    # `_orbit_radii_px` only ever reasons in flat pixel gaps between
    # bodies of unknown size -- it has no way to know the center it's
    # laying orbits around might be a giant/supergiant star (up to
    # `_STAR_MAX_R`) or the fixed-size drilled-into planet
    # (`_CENTER_PLANET_R`), either of which can otherwise overlap a
    # close, large first body outright. Shift the whole layout right by
    # any shortfall -- preserving every relative gap already computed --
    # rather than resizing anything.
    if radii and orbit_entries:
        first_kind, first_row = orbit_entries[0]
        needed = center_radius_px + _first_body_extent_px(first_kind, first_row, radii[0]) + _CENTER_CLEARANCE_PADDING_PX
        if radii[0] < needed:
            shift = needed - radii[0]
            radii = [r + shift for r in radii]

    label_sides = iter(_label_sides([
        (star_x + radius_px, row["name"])
        for (kind, row), radius_px in zip(orbit_entries, radii) if kind != "belt"
    ]))

    orbits = []
    bodies = []
    for (kind, row), radius_px in zip(orbit_entries, radii):
        cx, cy = star_x + radius_px, star_y
        if kind == "belt":
            half_band = _belt_band_px(row, radius_px) / 2
            orbits.append(_belt_band_svg(star_x, star_y, radius_px - half_band, radius_px + half_band, {
                "kind": "belt",
                "id": row["id"],
                "name": "Asteroid Belt",
                "density": row["density"].capitalize(),
                "distance": f'{row["lower_limit_km"]:,.0f} - {row["upper_limit_km"]:,.0f} km',
                "composition": row["composition_summary"],
            }, f'Asteroid belt ({row["density"]})'))
            continue

        orbits.append(_orbit_arc_svg(star_x, star_y, radius_px))
        moons = row.get("moons") or []
        scene_target = f'planet-{row["id"]}' if moons else None
        radius_fn = _moon_radius_px if row.get("_is_moon") else _planet_radius_px
        side = next(label_sides)
        bodies.append(_body_marker_svg(
            cx, cy, radius_fn(row["radius_km"]), row["planet_class"], row["body_type"],
            row["name"], "sysmap-moon" if row.get("_is_moon") else "sysmap-planet",
            _planet_attrs(row, kind="moon" if row.get("_is_moon") else "planet",
                          parent_name=row.get("_parent_name"), scene_target=scene_target),
            label_above=(side == "above"), show_label=(side is not None),
        ))

    return (
        f'<svg class="sysmap-svg" data-scene="{scene_id}" '
        f'viewBox="0 0 {_VIEW_W_PX:.0f} {_VIEW_H_PX:.0f}" role="group" aria-label="{esc(aria_label)}">'
        f'{"".join(orbits)}{star_svg}{"".join(bodies)}'
        "</svg>"
    )


def _render_system_scene(system, stars, planets, belts):
    orbit_entries = sorted(
        [("planet", planet) for planet in planets] + [("belt", belt) for belt in belts],
        key=lambda entry: entry[1]["distance_km"],
    )
    star_svg, primary_r = _stars_svg(_STAR_X_PX, _STAR_Y_PX, system, stars)
    return _render_row_svg(
        "system", f'System map for {system["name"]}', _STAR_X_PX, _STAR_Y_PX, star_svg, orbit_entries, primary_r,
    )


def _render_moon_scene(planet):
    moons = sorted(planet.get("moons") or [], key=lambda moon: moon["distance_km"])
    for moon in moons:
        moon["_is_moon"] = True
        moon["_parent_name"] = planet["name"]
    orbit_entries = [("planet", moon) for moon in moons]
    center_svg = _body_marker_svg(
        _STAR_X_PX, _STAR_Y_PX, _CENTER_PLANET_R, planet["planet_class"], planet["body_type"], planet["name"],
        "sysmap-planet", _planet_attrs(planet), is_self=True,
    )
    svg = _render_row_svg(
        f'planet-{planet["id"]}', f'Moons of {planet["name"]}', _STAR_X_PX, _STAR_Y_PX, center_svg, orbit_entries,
        _CENTER_PLANET_R,
    )
    # A moon scene starts hidden -- `static/systemmap.js` reveals it on
    # demand. This is a CSS class (`.sysmap-hidden`, toggled by
    # `classList`), not the HTML `hidden` attribute/`.hidden` IDL property:
    # `SVGElement` doesn't reflect that property to the attribute the way
    # `HTMLElement` does, so setting `.hidden` in JS silently only sets a
    # plain, attribute-less expando property on an `<svg>` -- confirmed
    # directly (`hasAttribute('hidden')` stayed false/true opposite of
    # what `.hidden` itself reported) -- leaving both this attribute and
    # any `[hidden]` CSS rule permanently out of sync with it.
    return svg.replace('class="sysmap-svg"', 'class="sysmap-svg sysmap-hidden"', 1)


def render_system_map_panel(system, stars, planets, belts):
    """
    Builds the "System Map" panel embedded in `system.py`: an `<svg>` per
    scene (the whole system, plus one more for every planet with moons --
    see the module docstring for why this doesn't need to recurse any
    deeper than that) and an info side panel that `static/systemmap.js`
    fills in on click and swaps between scenes.

    Args:
        system (sqlite3.Row): The `star_systems` row.
        stars (list[dict]): 1 entry (single star) or 2 (primary, then
                            secondary), each with `star_type`,
                            `temperature_k`, `radius_km`, `luminosity_w`,
                            and the `table_*` display strings.
        planets (list[dict]): Each a `planets` row plus a `moons` key
                              (list of `moons` rows, possibly empty).
        belts (list[dict]): `asteroid_belts` rows.

    Returns:
        str: A complete `<section class="panel">` block.
    """
    scenes = [_render_system_scene(system, stars, planets, belts)]
    scenes.extend(
        _render_moon_scene(planet) for planet in planets if planet.get("moons")
    )

    if planets or belts or stars:
        info_panel = (
            '<aside class="starmap-info" id="sysmap-info">'
            '<p class="hint">Click a star, planet, moon, or asteroid belt for details.</p></aside>'
        )
    else:
        info_panel = '<aside class="starmap-info" id="sysmap-info"><p class="hint">Nothing to show yet.</p></aside>'

    return f"""
<section class="panel">
<div class="panel-header">
  <h2>System Map</h2>
  <span class="hint">Click a planet with moons to view its moon system &middot; circle size &asymp; body radius (log scale) &middot; color &asymp; planet class</span>
</div>
<div class="starmap-layout" id="sysmap-root">
<div class="starmap-viewport sysmap-viewport">
{''.join(scenes)}
</div>
<div class="starmap-side">
<div class="sysmap-crumb" id="sysmap-crumb"></div>
{info_panel}
</div>
</div>
</section>
"""
