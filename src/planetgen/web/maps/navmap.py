# planetgen/web/maps/navmap.py

"""
NAV map: a flat, top-down SVG plot of the galactic X-Y plane showing a
NAV request's origin, destination, and (when one was found) the optimal
route's intermediate hops -- the "actual plotted image/diagram of the two
points" `docs/api.md`/`docs/TODO.md` flagged as still open once
`queryDb.nav_between` started returning course/route data as numbers and
links only (the NAV page, `/nav`).

Deliberately modeled on `galaxymap.py`'s flat 2D SVG rather than
`starmap.py`'s rotatable 3D CSS scene: like a galaxy Quadrant, this map is
blind to height by design. Both NAV frames (`planetgen.galaxy.navigation`'s
Galactic and Sector frames) keep galactic +Z as up, so a course's bearing
lives in this X-Y plane, and the mark the NAV page's course panel reports
covers the third axis -- a route's
hops are typically nowhere near coplanar with the direct line in practice,
so this map is honest about showing a galactic-plane projection rather
than faking a 3D perspective a flat, static SVG (no drag-to-rotate script,
unlike starmap.py) can't render usefully anyway. A small arrow marks the
origin's bearing 000 (toward the frame's center: the galactic core or the
sector's center) so the plotted picture and the course panel's bearing
read as the same thing.

Auto-scaled to whatever positions it's given -- there's no fixed "sector
size" to normalize against the way starmap.py's per-sector `edge_mpc`
gives it -- so the map's own extent is the bounding box of every point
plotted, padded, with one uniform light-years-per-pixel ratio on both axes
(never stretched independently, which would visually distort bearing
angles away from what they actually are).

NAV.41: the panel is wide (`_SVG_WIDTH` by `_SVG_HEIGHT`) and spans the
page's width, and its text -- stop names, the bearing arrow's label and
the scale bar's length -- is HTML laid over the SVG at percentage
positions, so it stays at the page's body text size however wide the
panel is drawn, where SVG text would shrink with the drawing. A hop label
that would overlap one already placed is dropped (`_place_labels`); the
route list below the map names every stop anyway.
"""

import math

from planetgen.web.lib.fmt import esc, format_distance_ly

_SVG_WIDTH = 600.0
_SVG_HEIGHT = 360.0
_CENTER_X = _SVG_WIDTH / 2
_CENTER_Y = _SVG_HEIGHT / 2
_MARGIN_X = 60.0  # room for labels beside the outermost points
_MARGIN_Y = 44.0  # room for labels above and the scale bar below
_HALF_W = _SVG_WIDTH / 2 - _MARGIN_X
_HALF_H = _SVG_HEIGHT / 2 - _MARGIN_Y

# Label size, in SVG units, at the narrowest panel a label must not
# overlap at (a phone's panel, about 320 px across): 1rem text with an average
# glyph taken as 0.62em wide (a little over the real average, so estimates err wide). A hop label overlapping another at this width
# but not at `_WIDE_PANEL_PX` only shows on wider screens (CSS hides
# `.navmap-label-wide` on phones); one that overlaps even there is left
# out.
_NARROW_PANEL_PX = 320.0
_WIDE_PANEL_PX = 960.0
_LABEL_FONT_PX = 16.0
_LABEL_CHAR_EM = 0.62

_MIN_SPAN_LY = 1.0  # keeps a map of two very close (or identical) points sane
_PADDING_FRACTION = 0.08  # extra room around the tightest bounding box

_ORIGIN_R = 7.0
_DEST_R = 7.0
_HOP_R = 4.0

# The bearing arrow sits in the top-left corner, under its label, out of
# the route's way.
_COMPASS_X = 30.0
_COMPASS_Y = 64.0
_COMPASS_REACH = 22.0
_SCALE_BAR_FRACTION = 0.4  # of _HALF_H


def _project_all(points, px_per_ly, center_x, center_y):
    """
    Projects every `(x, y)` in `points` (galactic-plane light-year
    coordinates, z ignored -- see the module docstring) into SVG pixel
    space, centered on `(center_x, center_y)` and scaled by `px_per_ly`.

    Args:
        points (list[tuple]): `(x, y)` light-year positions.
        px_per_ly (float): Pixels per light-year, uniform on both axes.
        center_x (float): Light-year x of the map's own center.
        center_y (float): Light-year y of the map's own center.

    Returns:
        list[tuple]: `(svg_x, svg_y)` for each input point, in order.
    """
    projected = []
    for x, y in points:
        # +Y is "up" in this project's galactic-plane convention (matching
        # navigation.py's frames and galaxymap.py's own projection);
        # SVG's y axis grows downward, so the sign flips here once, at the
        # one place a light-year position becomes a screen coordinate.
        svg_x = _CENTER_X + (x - center_x) * px_per_ly
        svg_y = _CENTER_Y - (y - center_y) * px_per_ly
        projected.append((svg_x, svg_y))
    return projected


def _scale(points):
    """
    Picks one uniform pixels-per-light-year ratio and a light-year center
    point that fits every `(x, y)` in `points` inside the drawing area
    (`_HALF_W` by `_HALF_H` either side of center) with `_PADDING_FRACTION`
    of padding.

    Args:
        points (list[tuple]): `(x, y)` light-year positions to fit.

    Returns:
        tuple: `(px_per_ly, center_x, center_y)`.
    """
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    # One ratio for both axes, from whichever axis is the tighter fit, keeps
    # the scale uniform -- a stretched map would distort bearing angles.
    padding = 1.0 + 2 * _PADDING_FRACTION
    span_x = max(max_x - min_x, _MIN_SPAN_LY) * padding
    span_y = max(max_y - min_y, _MIN_SPAN_LY) * padding
    px_per_ly = min(2 * _HALF_W / span_x, 2 * _HALF_H / span_y)

    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    return px_per_ly, center_x, center_y


def _compass_html(origin_xy, frame_center_xy):
    """An arrow in the panel's corner along the origin's bearing 000 --
    from the origin toward the frame's center, flattened onto the galactic
    plane (`navigation.course_between`) -- so the plot's orientation reads
    the same as the course panel's bearing. Falls back to +X, as
    `course_between` does, when the origin sits right over the center."""
    dx = frame_center_xy[0] - origin_xy[0]
    dy = frame_center_xy[1] - origin_xy[1]
    length = math.hypot(dx, dy)
    if length <= 1e-9 * max(math.hypot(*frame_center_xy), math.hypot(*origin_xy), 1.0):
        dx, dy, length = 1.0, 0.0, 1.0
    # SVG's y axis grows downward, hence the sign flip on dy.
    tip_x = _COMPASS_X + _COMPASS_REACH * dx / length
    tip_y = _COMPASS_Y - _COMPASS_REACH * dy / length
    line = (
        f'<circle class="navmap-compass-ring" cx="{_COMPASS_X:.1f}" cy="{_COMPASS_Y:.1f}" r="{_COMPASS_REACH:.1f}"/>'
        f'<line class="navmap-compass" x1="{_COMPASS_X:.1f}" y1="{_COMPASS_Y:.1f}" '
        f'x2="{tip_x:.1f}" y2="{tip_y:.1f}"/>'
    )
    label = {"text": "Bearing 000", "x": 8.0, "y": _COMPASS_Y - _COMPASS_REACH - 4, "anchor": "start"}
    return line, label


def _scale_bar_html(px_per_ly):
    """A fixed-length legend bar (see the module docstring on why this
    map's extent isn't a fixed, round physical size the way a sector's
    `edge_mpc` is) labeled with its own real light-year length, so the
    plotted picture still carries a sense of physical scale."""
    bar_px = _HALF_H * _SCALE_BAR_FRACTION
    bar_ly = bar_px / px_per_ly
    x0 = _SVG_WIDTH - _MARGIN_X - bar_px
    x1 = _SVG_WIDTH - _MARGIN_X
    y = _SVG_HEIGHT - 14
    lines = (
        f'<line class="navmap-scale-bar" x1="{x0:.1f}" y1="{y:.1f}" x2="{x1:.1f}" y2="{y:.1f}"/>'
        f'<line class="navmap-scale-tick" x1="{x0:.1f}" y1="{y - 4:.1f}" x2="{x0:.1f}" y2="{y + 4:.1f}"/>'
        f'<line class="navmap-scale-tick" x1="{x1:.1f}" y1="{y - 4:.1f}" x2="{x1:.1f}" y2="{y + 4:.1f}"/>'
    )
    label = {"text": format_distance_ly(bar_ly), "x": x1, "y": y - 6, "anchor": "end"}
    return lines, label


def _label_html(css_class, svg_x, svg_y, anchor, text, href=None):
    """An HTML label over the SVG whose bottom edge sits at `(svg_x,
    svg_y)`, aligned to it by `anchor` ("start", "middle" or "end", as
    SVG's `text-anchor`); a link when `href` is given."""
    tag = "a" if href else "span"
    href_attr = f' href="{esc(href)}"' if href else ""
    left = 100 * svg_x / _SVG_WIDTH
    top = 100 * svg_y / _SVG_HEIGHT
    return (
        f'<{tag} class="navmap-label navmap-anchor-{anchor} {css_class}"{href_attr} '
        f'style="left: {left:.2f}%; top: {top:.2f}%">{esc(text)}</{tag}>'
    )


def _label_box(text, svg_x, svg_y, anchor, panel_px):
    """The `(x0, y0, x1, y1)` box, in SVG units, a label takes up when the
    panel is drawn `panel_px` wide."""
    units_per_px = _SVG_WIDTH / panel_px
    width = len(text) * _LABEL_FONT_PX * _LABEL_CHAR_EM * units_per_px
    height = _LABEL_FONT_PX * 1.25 * units_per_px
    x0 = {"start": svg_x, "middle": svg_x - width / 2, "end": svg_x - width}[anchor]
    return x0, svg_y - height, x0 + width, svg_y


def _overlaps(box, boxes):
    return any(box[0] < b[2] and b[0] < box[2] and box[1] < b[3] and b[1] < box[3] for b in boxes)


def _anchor_for(svg_x):
    """Labels near the panel's sides grow inward so they stay on it."""
    if svg_x < _MARGIN_X + 40:
        return "start"
    if svg_x > _SVG_WIDTH - _MARGIN_X - 40:
        return "end"
    return "middle"


def _place_labels(labels, fixed=()):
    """
    Picks which point labels to show (NAV.41): the origin and destination
    always; each hop, in route order, only where it doesn't overlap a label
    already placed. Returns `(label, width_class)` pairs, `width_class`
    being "" (every width), "navmap-label-wide" (fits only on wider
    screens) -- labels that fit at neither width are left out.

    Args:
        labels (list[dict]): Each with `text`, `x`, `y`, `anchor` and
            `always` (true for the origin and destination).
        fixed (iterable[dict]): Labels always shown (the bearing and
            scale labels) that hops must not cover, same keys.
    """
    narrow = [_label_box(f["text"], f["x"], f["y"], f["anchor"], _NARROW_PANEL_PX) for f in fixed]
    wide = [_label_box(f["text"], f["x"], f["y"], f["anchor"], _WIDE_PANEL_PX) for f in fixed]
    placed = []
    for label in sorted(labels, key=lambda item: not item["always"]):
        args = (label["text"], label["x"], label["y"], label["anchor"])
        narrow_box = _label_box(*args, _NARROW_PANEL_PX)
        wide_box = _label_box(*args, _WIDE_PANEL_PX)
        if label["always"] or not (_overlaps(narrow_box, narrow) or _overlaps(wide_box, wide)):
            narrow.append(narrow_box)
            wide.append(wide_box)
            placed.append((label, ""))
        elif not _overlaps(wide_box, wide):
            wide.append(wide_box)
            placed.append((label, "navmap-label-wide"))
    return placed


def _point_url(link_url, waypoint):
    # A route's own intermediate hops are always real systems (see
    # queryDb.nav_between's docstring), but the origin/destination
    # themselves can each be a standalone phenomenon instead -- linked to
    # its phenomenon page rather than a system page.
    if waypoint.get("kind") == "phenomenon":
        return link_url("phenomenon", phenomenon_type=waypoint["type"], phenomenon_id=waypoint["id"])
    return link_url("system", system_id=waypoint["id"])


def _course_note(waypoint):
    """NAV.42: a stop's tooltip line for the course to the next stop, or nothing (the last stop, or no route)."""
    return f" (next stop: {waypoint['course']})" if waypoint.get("course") else ""


def _point_html(url, waypoint, svg_x, svg_y, css_class, radius):
    return (
        f'<a class="navmap-point {css_class}" href="{esc(url)}">'
        f'<title>{esc(waypoint["name"])}{esc(_course_note(waypoint))}</title>'
        f'<circle cx="{svg_x:.1f}" cy="{svg_y:.1f}" r="{radius:.1f}"/></a>'
    )


def render_nav_map_panel(link_url, waypoints, has_route, frame_center=(0.0, 0.0, 0.0)):
    """
    Builds the "NAV Map" panel: a flat, top-down SVG plot of an origin, a
    destination, and (when `has_route` is true) the optimal route's
    intermediate hops, all in the galactic X-Y plane (see the module
    docstring on why altitude/z is left out).

    Args:
        link_url (callable): `link_url(name, **params)` -> URL (the
                       Flask NAV page passes `web.helpers.page_url`);
                       called as `link_url("system", system_id=...)` or,
                       for a phenomenon origin/destination,
                       `link_url("phenomenon", phenomenon_type=...,
                       phenomenon_id=...)`, for each point's plain
                       `<a href>` link.
        waypoints (list[dict]): Ordered origin-to-destination, each with
                                `id`, `name`, `position` (an `(x, y, z)`
                                light-year tuple in `nav_between`'s scope
                                frame -- `origin_position`/
                                `destination_position`/a route's
                                `positions` entries, see `queryDb.
                                nav_between`), `role` (`"origin"`,
                                `"destination"`, or `"hop"`), and `kind`
                                (`"system"`, the default read when the key
                                is absent, or `"phenomenon"` -- a `"hop"`
                                is always `"system"`, only the origin/
                                destination can be a phenomenon; when
                                `kind == "phenomenon"`, `type` must also be
                                set to its `queryDb._PHENOMENON_TYPE_TO_TABLE`
                                key). Always at least the origin and destination; any
                                `"hop"` entries in between are the
                                optimal route's intermediate systems, in
                                path order.
        has_route (bool): Whether an optimal route (possibly with no
                          intermediate hops, i.e. adjacent) was found --
                          draws a solid polyline through every waypoint in
                          order when true, in addition to the dashed
                          direct line between origin and destination
                          (always drawn, even with no route: same-system
                          NAV and an unreachable pair both still have a
                          direct course).
        frame_center (tuple): The course frame's center, in the same
                          frame as the waypoints' positions -- the compass
                          arrow points from the origin toward it. Defaults
                          to `(0, 0, 0)`, the center of both NAV frames
                          (sector-local and galaxy-frame).

    Returns:
        str: A complete `<section class="panel">` block.
    """
    points_xy = [(w["position"][0], w["position"][1]) for w in waypoints]
    px_per_ly, center_x, center_y = _scale(points_xy)
    projected = _project_all(points_xy, px_per_ly, center_x, center_y)

    origin_xy, destination_xy = projected[0], projected[-1]
    direct_line = (
        f'<line class="navmap-direct-line" x1="{origin_xy[0]:.1f}" y1="{origin_xy[1]:.1f}" '
        f'x2="{destination_xy[0]:.1f}" y2="{destination_xy[1]:.1f}"/>'
    )

    route_line = ""
    if has_route and len(waypoints) > 1:
        path_d = "M " + " L ".join(f"{x:.1f} {y:.1f}" for x, y in projected)
        route_line = f'<path class="navmap-route-line" d="{path_d}" fill="none"/>'

    points_html = []
    labels = []
    for waypoint, (svg_x, svg_y) in zip(waypoints, projected):
        css_class, radius = {
            "origin": ("navmap-origin", _ORIGIN_R),
            "destination": ("navmap-destination", _DEST_R),
        }.get(waypoint["role"], ("navmap-hop", _HOP_R))
        url = _point_url(link_url, waypoint)
        points_html.append(_point_html(url, waypoint, svg_x, svg_y, css_class, radius))
        labels.append({
            "text": waypoint["name"], "x": svg_x, "y": svg_y - radius - 3, "anchor": _anchor_for(svg_x),
            "always": waypoint["role"] != "hop", "url": url, "class": css_class,
        })
    compass_line, compass_label = _compass_html(points_xy[0], frame_center[:2])
    scale_lines, scale_label = _scale_bar_html(px_per_ly)
    labels_html = [
        _label_html(f'{label["class"]}-label {width_class}'.strip(), label["x"], label["y"], label["anchor"],
                    label["text"], href=label["url"])
        for label, width_class in _place_labels(labels, fixed=(compass_label, scale_label))
    ]
    labels_html += [
        _label_html(css_class, label["x"], label["y"], label["anchor"], label["text"])
        for css_class, label in (("navmap-compass-label", compass_label), ("navmap-scale-label", scale_label))
    ]

    body = f"{direct_line}{route_line}{''.join(points_html)}{compass_line}{scale_lines}"

    svg = (
        f'<svg class="navmap-svg" viewBox="0 0 {_SVG_WIDTH:.0f} {_SVG_HEIGHT:.0f}" role="img" '
        'aria-label="Top-down plot of the origin, destination, and optimal route in the galactic plane. '
        'Click a point for its system.">'
        f"{body}</svg>"
    )

    route_hint = (
        "the solid line is the optimal route via adjacent systems"
        if has_route
        else "no route via adjacent systems was found, so only the direct line is shown"
    )

    return f"""
<section class="panel">
<div class="panel-header">
  <h2 class="sr-only">NAV Map</h2>
</div>
<div class="navmap-viewport">
{svg}<div class="navmap-labels">{''.join(labels_html)}</div>
</div>
<div class="starmap-controls">
  <button type="button" class="starmap-btn" data-dialog-open="navmap-help" title="How to read the map">Map help</button>
</div>
<sl-dialog id="navmap-help" class="map-help-dialog" label="NAV Map help">
  <ul>
    <li>Top-down view (the galactic X-Y plane; height is not shown, see the course's mark above).</li>
    <li>The dashed line is the direct course; {route_hint}.</li>
    <li>Click a point for its system.</li>
  </ul>
  <sl-button slot="footer" data-dialog-close>Close</sl-button>
</sl-dialog>
</section>
"""
