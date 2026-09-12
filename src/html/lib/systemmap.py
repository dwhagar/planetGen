# html/lib/systemmap.py

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
like `starmap.py`/`sectormap.js`'s own click-for-info pattern.
"""

import colorsys
import math
import statistics

from fmt import esc
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
# (stroke width = its own radial spread) around its own anchor, the
# natural true-position shape for "every angle at roughly this distance",
# rather than the old schematic map's directional shaded band. Clamped so
# a belt with almost no radial spread is still visibly a ring and one with
# an enormous spread doesn't swallow its neighbors.
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
# gas-giant palette, and N/Q/V's outliers (Venus-hot mustard, an eccentric-
# orbit violet, a high-gravity super-Earth magenta).
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


def _label_sides_2d(entries):
    """
    Given `[(cx, cy, marker_r, name), ...]`, returns one `"below"`/
    `"above"`/`None` per entry: which band that body's name label should
    be drawn in, or `None` to skip drawing it. Unlike the old fixed
    schematic map's `_label_sides` (which only ever had to de-collide
    labels along one shared horizontal band, since every marker sat on the
    same line), markers here can be anywhere in the plane, so this checks
    each candidate label's own bounding box against every *other* label
    already placed, in true 2D, rather than assuming any shared axis.
    A body whose label collides in both bands drops its label entirely --
    that body's name, and everything else about it, is still one click
    away in the info panel.

    Args:
        entries (list[tuple]): `(cx, cy, marker_r, name)`, in the order
            markers should be given placement priority (earlier entries
            never yield to a later one).

    Returns:
        list[str or None]: One entry per input, in the same order.
    """
    accepted_rects = []
    sides = []
    for cx, cy, marker_r, name in entries:
        half_w = _label_half_width_px(name)
        candidates = {
            "below": (
                cx - half_w, cy + marker_r + _LABEL_GAP_PX,
                cx + half_w, cy + marker_r + _LABEL_GAP_PX + _LABEL_HALF_HEIGHT_PX * 2,
            ),
            "above": (
                cx - half_w, cy - marker_r - _LABEL_GAP_PX - _LABEL_HALF_HEIGHT_PX * 2,
                cx + half_w, cy - marker_r - _LABEL_GAP_PX,
            ),
        }
        chosen = None
        for side in ("below", "above"):
            rect = candidates[side]
            if not any(_rects_overlap(rect, other) for other in accepted_rects):
                chosen = side
                accepted_rects.append(rect)
                break
        sides.append(chosen)
    return sides


def _data_attrs(attrs):
    return "".join(f' data-{key}="{esc(value)}"' for key, value in attrs.items() if value is not None)


_LIFE_BADGE_FILL = "#3ecf6e"
_LIFE_BADGE_STROKE = "#1c7a3e"


def _body_marker_svg(cx, cy, r_px, planet_class, body_type, label_text, extra_class, attrs, is_self=False,
                      label_above=False, show_label=True, has_life=False):
    """
    Builds one clickable `<g>` for a planet or moon: a filled/stroked
    circle colored by `planet_class` (see `_CLASS_COLORS`), the class
    letter centered inside it once the circle is big enough to hold text
    legibly, a small tilted ring behind gas giants (`body_type == "g"`)
    for an at-a-glance silhouette cue beyond just color, a small green
    "supports life" badge at the marker's own edge when `has_life` is set
    (`planets`/`moons`.`life_chemical` is non-NULL), and (when
    `show_label` is set) a name label underneath -- or, when `label_above`
    is also set (see `_label_sides_2d`), above, to keep it clear of a
    tightly-packed neighbor's own label. `show_label` alone going False
    (also from `_label_sides_2d`, when neither band could clear a tight
    enough cluster) only drops the visible text, not the body's name from
    `aria-label` -- it's still reachable, just via a click on the marker
    rather than a glance at the diagram. `is_self` marks the one body a
    moon-scene is *about* (the planet drilled into, now standing in for
    the scene's own origin) -- drawn with a soft halo and picked out by
    `static/systemmap.js` as the default info-panel content when that
    scene opens, via its `data-self="true"` marker.
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
    if has_life:
        badge_r = max(2.5, r_px * 0.32)
        badge_x = cx + r_px * 0.72
        badge_y = cy - r_px * 0.72
        parts.append(
            f'<circle class="sysmap-life-badge" cx="{badge_x:.1f}" cy="{badge_y:.1f}" r="{badge_r:.1f}" '
            f'fill="{_LIFE_BADGE_FILL}" stroke="{_LIFE_BADGE_STROKE}"><title>Supports life</title></circle>'
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
        "life": planet.get("life_chemical"),
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
        "distance": f'{belt["lower_limit_km"]:,.0f} - {belt["upper_limit_km"]:,.0f} km',
        "composition": belt["composition_summary"],
    }
    label = f'Asteroid belt ({belt["density"]})'
    return (
        f'<circle class="sysmap-belt" tabindex="0" role="button"{_data_attrs(attrs)} '
        f'aria-label="{esc(label)}" cx="{cx:.1f}" cy="{cy:.1f}" r="{radius_px:.1f}" '
        f'stroke-width="{band_px:.1f}"></circle>'
    )


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


def _render_system_scene(system, stars, planets, belts):
    """Builds the "whole system" scene -- see the module docstring for the
    barycenter/anchor model this uses for a binary pair."""
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
        lx, ly = planet.get("position_x_km") or 0.0, planet.get("position_y_km") or 0.0
        local_r_list.append(math.hypot(lx, ly) or planet.get("distance_km") or 0.0)
    for belt in belts:
        local_r_list.append(belt["distance_km"])

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
        star_markers.append({
            "type": "star", "cx": cx, "cy": cy, "r": _star_radius_px(star["radius_km"]),
            "star": star, "suffix": (" -- Primary" if index == 0 else " -- Secondary") if is_binary else "",
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
    markers = list(star_markers)
    for planet in planets:
        ax_px, ay_px = anchor_px(planet.get("star_id"))
        lx_km, ly_km = planet.get("position_x_km") or 0.0, planet.get("position_y_km") or 0.0
        r_px = _radial_px(math.hypot(lx_km, ly_km), lo, hi)
        cx, cy = _polar_to_px(ax_px, ay_px, r_px, lx_km, ly_km)
        orbit_paths.append(f'<circle class="sysmap-orbit" cx="{ax_px:.1f}" cy="{ay_px:.1f}" r="{r_px:.1f}"></circle>')
        markers.append({"type": "planet", "cx": cx, "cy": cy, "r": _planet_radius_px(planet["radius_km"]), "row": planet})

    for belt in belts:
        ax_px, ay_px = anchor_px(belt.get("star_id"))
        r_px = _radial_px(belt["distance_km"], lo, hi)
        band_px = _belt_band_px(belt, r_px)
        orbit_paths.append(_belt_ring_svg(ax_px, ay_px, r_px, band_px, belt))

    _relax_markers(markers, _MARKER_GAP_PX)

    star_svgs = []
    planet_markers = []
    for marker in markers:
        if marker["type"] != "star":
            planet_markers.append(marker)
            continue
        star = marker["star"]
        is_primary = not marker["suffix"] or marker["suffix"].endswith("Primary")
        star_svgs.append(_star_marker_svg(marker["cx"], marker["cy"], marker["r"], star, {
            "kind": "star",
            "name": f'{system["name"]}{marker["suffix"]}',
            "role": ("Primary" if is_primary else "Secondary") if is_binary else "Single",
            "type": star["star_type"], "temp": f'{int(star["temperature_k"])} K',
            "mass": format_star_mass(star["mass_kg"]), "radius": format_star_radius(star["radius_km"]),
            "lum": format_star_luminosity(star["luminosity_w"]),
        }))

    sides = _label_sides_2d([(m["cx"], m["cy"], m["r"], m["row"]["name"]) for m in planet_markers])
    body_svgs = []
    for marker, side in zip(planet_markers, sides):
        row = marker["row"]
        scene_target = f'planet-{row["id"]}' if row.get("moons") else None
        attrs = _planet_attrs(row, kind="planet", scene_target=scene_target)
        body_svgs.append(_body_marker_svg(
            marker["cx"], marker["cy"], marker["r"], row["planet_class"], row["body_type"], row["name"],
            "sysmap-planet", attrs, label_above=(side == "above"), show_label=(side is not None),
            has_life=bool(row.get("life_chemical")),
        ))

    return (
        f'<svg class="sysmap-svg" data-scene="system" '
        f'viewBox="0 0 {_VIEW_SIZE_PX:.0f} {_VIEW_SIZE_PX:.0f}" role="group" '
        f'aria-label="System map for {esc(system["name"])}">'
        f'{"".join(orbit_paths)}{"".join(star_svgs)}{"".join(body_svgs)}'
        "</svg>"
    )


def _render_moon_scene(planet):
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
    # The drilled-into planet itself acts as a fixed obstacle here (same
    # `fixed` treatment `_render_system_scene` gives its own star(s)) --
    # `_CENTER_PLANET_R` (36px) is large enough to otherwise swallow a
    # close-orbiting first moon outright once real angle/distance placement
    # puts it nearby.
    markers = [{"type": "center", "cx": _CENTER_PX, "cy": _CENTER_PX, "r": _CENTER_PLANET_R, "fixed": True}]
    for moon in moons:
        lx, ly = moon.get("position_x_km") or 0.0, moon.get("position_y_km") or 0.0
        r_px = _radial_px(math.hypot(lx, ly), lo, hi)
        cx, cy = _polar_to_px(_CENTER_PX, _CENTER_PX, r_px, lx, ly)
        orbit_paths.append(f'<circle class="sysmap-orbit" cx="{_CENTER_PX:.1f}" cy="{_CENTER_PX:.1f}" r="{r_px:.1f}"></circle>')
        markers.append({"type": "moon", "cx": cx, "cy": cy, "r": _moon_radius_px(moon["radius_km"]), "row": moon})

    _relax_markers(markers, _MARKER_GAP_PX)
    markers = [m for m in markers if m["type"] == "moon"]

    center_svg = _body_marker_svg(
        _CENTER_PX, _CENTER_PX, _CENTER_PLANET_R, planet["planet_class"], planet["body_type"], planet["name"],
        "sysmap-planet", _planet_attrs(planet), is_self=True,
        has_life=bool(planet.get("life_chemical")),
    )

    sides = _label_sides_2d([(m["cx"], m["cy"], m["r"], m["row"]["name"]) for m in markers])
    body_svgs = []
    for marker, side in zip(markers, sides):
        row = marker["row"]
        attrs = _planet_attrs(row, kind="moon", parent_name=row.get("_parent_name"))
        body_svgs.append(_body_marker_svg(
            marker["cx"], marker["cy"], marker["r"], row["planet_class"], row["body_type"], row["name"],
            "sysmap-moon", attrs, label_above=(side == "above"), show_label=(side is not None),
            has_life=bool(row.get("life_chemical")),
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
    return (
        f'<svg class="sysmap-svg sysmap-hidden" data-scene="planet-{planet["id"]}" '
        f'viewBox="0 0 {_VIEW_SIZE_PX:.0f} {_VIEW_SIZE_PX:.0f}" role="group" '
        f'aria-label="Moons of {esc(planet["name"])}">'
        f'{"".join(orbit_paths)}{center_svg}{"".join(body_svgs)}'
        "</svg>"
    )


def render_system_map_panel(system, stars, planets, belts):
    """
    Builds the "System Map" panel embedded in `system.py`: an `<svg>` per
    scene (the whole system, plus one more for every planet with moons --
    see the module docstring for why this doesn't need to recurse any
    deeper than that) and an info side panel that `static/systemmap.js`
    fills in on click and swaps between scenes.

    Args:
        system (dict): The `star_systems` row, including
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
  <span class="hint">True top-down positions (real angle, log-scaled distance) &middot; click a planet with moons to view its moon system &middot; circle size &asymp; body radius (log scale) &middot; color &asymp; planet class &middot; <span class="sysmap-legend-life-badge" aria-hidden="true"></span> supports life</span>
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
