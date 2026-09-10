# html/lib/navmap.py

"""
NAV map: a flat, top-down SVG plot of the galactic X-Y plane showing a
NAV request's origin, destination, and (when one was found) the optimal
route's intermediate hops -- the "actual plotted image/diagram of the two
points" `docs/api.md`/`docs/TODO.md` flagged as still open once
`queryDb.nav_between` started returning course/route data as numbers and
links only (`html/nav.py`).

Deliberately modeled on `galaxymap.py`'s flat 2D SVG rather than
`starmap.py`'s rotatable 3D CSS scene: like a galaxy Quadrant, this map is
blind to height (altitude) by design. `stellarObjects.navigation`'s own
azimuth convention (see that module's docstring) is already defined purely
within the galactic X-Y plane, and the numeric Altitude figure `html/
nav.py`'s course panel already reports covers the third axis -- a route's
hops are typically nowhere near coplanar with the direct line in practice,
so this map is honest about showing an azimuth-plane projection rather
than faking a 3D perspective a flat, static SVG (no drag-to-rotate script,
unlike starmap.py) can't render usefully anyway. A small arrow marks the
scene's own +X direction (azimuth 0 degrees) so the plotted picture and
the course panel's azimuth figure read as the same thing.

Auto-scaled to whatever positions it's given -- there's no fixed "sector
size" to normalize against the way starmap.py's per-sector `edge_mpc`
gives it -- so the map's own extent is the bounding box of every point
plotted, padded, with one uniform light-years-per-pixel ratio on both axes
(never stretched independently, which would visually distort azimuth
angles away from what they actually are).
"""

import math

from dbutil import esc

_SVG_SIZE = 360.0
_CENTER = _SVG_SIZE / 2
_MAP_RADIUS_PX = _SVG_SIZE / 2 - 36  # leaves margin for labels/compass/scale bar

_MIN_SPAN_LY = 1.0  # keeps a map of two very close (or identical) points sane
_PADDING_FRACTION = 0.22  # extra room around the tightest bounding box

_ORIGIN_R = 7.0
_DEST_R = 7.0
_HOP_R = 4.0

_COMPASS_REACH_FRACTION = 0.9
_SCALE_BAR_FRACTION = 0.4  # of _MAP_RADIUS_PX


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
        # navigation.py's azimuth math and galaxymap.py's own projection);
        # SVG's y axis grows downward, so the sign flips here once, at the
        # one place a light-year position becomes a screen coordinate.
        svg_x = _CENTER + (x - center_x) * px_per_ly
        svg_y = _CENTER - (y - center_y) * px_per_ly
        projected.append((svg_x, svg_y))
    return projected


def _scale(points):
    """
    Picks one uniform pixels-per-light-year ratio and a light-year center
    point that fits every `(x, y)` in `points` inside `_MAP_RADIUS_PX` with
    `_PADDING_FRACTION` of padding.

    Args:
        points (list[tuple]): `(x, y)` light-year positions to fit.

    Returns:
        tuple: `(px_per_ly, center_x, center_y)`.
    """
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    # A single span (the larger of the two axes) keeps the scale uniform --
    # a tall, narrow bounding box still gets a square map rather than a
    # stretched one that would distort azimuth angles.
    span = max(max_x - min_x, max_y - min_y, _MIN_SPAN_LY)
    padded_span = span * (1.0 + 2 * _PADDING_FRACTION)
    px_per_ly = (_MAP_RADIUS_PX * 2) / padded_span

    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    return px_per_ly, center_x, center_y


def _compass_html():
    """An arrow from the map's center toward +X, labeled with the azimuth
    it marks -- ties this map's orientation to the course panel's own
    azimuth figure (`navigation.course_between`'s "counterclockwise from
    +X" convention) rather than leaving the plot's rotation ambiguous."""
    reach = _MAP_RADIUS_PX * _COMPASS_REACH_FRACTION
    tip_x = _CENTER + reach
    tip_y = _CENTER
    return (
        f'<line class="navmap-compass" x1="{_CENTER:.1f}" y1="{_CENTER:.1f}" '
        f'x2="{tip_x:.1f}" y2="{tip_y:.1f}"/>'
        f'<text class="navmap-compass-label" x="{tip_x - 4:.1f}" y="{tip_y - 6:.1f}" '
        'text-anchor="end">+X (azimuth 0&deg;)</text>'
    )


def _scale_bar_html(px_per_ly):
    """A fixed-length legend bar (see the module docstring on why this
    map's extent isn't a fixed, round physical size the way a sector's
    `edge_mpc` is) labeled with its own real light-year length, so the
    plotted picture still carries a sense of physical scale."""
    bar_px = _MAP_RADIUS_PX * _SCALE_BAR_FRACTION
    bar_ly = bar_px / px_per_ly
    x0 = _CENTER - bar_px / 2
    x1 = _CENTER + bar_px / 2
    y = _SVG_SIZE - 14
    return (
        f'<line class="navmap-scale-bar" x1="{x0:.1f}" y1="{y:.1f}" x2="{x1:.1f}" y2="{y:.1f}"/>'
        f'<line class="navmap-scale-tick" x1="{x0:.1f}" y1="{y - 4:.1f}" x2="{x0:.1f}" y2="{y + 4:.1f}"/>'
        f'<line class="navmap-scale-tick" x1="{x1:.1f}" y1="{y - 4:.1f}" x2="{x1:.1f}" y2="{y + 4:.1f}"/>'
        f'<text class="navmap-scale-label" x="{_CENTER:.1f}" y="{y - 8:.1f}" text-anchor="middle">'
        f"{bar_ly:,.2f} ly</text>"
    )


def _point_html(db_name, waypoint, svg_x, svg_y, css_class, radius):
    href = f'system.py?db={esc(db_name)}&amp;id={waypoint["id"]}'
    return (
        f'<a class="navmap-point {css_class}" href="{href}"><title>{esc(waypoint["name"])}</title>'
        f'<circle cx="{svg_x:.1f}" cy="{svg_y:.1f}" r="{radius:.1f}"/>'
        f'<text x="{svg_x:.1f}" y="{svg_y - radius - 5:.1f}" text-anchor="middle">'
        f'{esc(waypoint["name"])}</text></a>'
    )


def render_nav_map_panel(db_name, waypoints, has_route):
    """
    Builds the "NAV Map" panel: a flat, top-down SVG plot of an origin, a
    destination, and (when `has_route` is true) the optimal route's
    intermediate hops, all in the galactic X-Y plane (see the module
    docstring on why altitude/z is left out).

    Args:
        db_name (str): The current `?db=` value, used to build each
                       point's link back to `system.py`.
        waypoints (list[dict]): Ordered origin-to-destination, each with
                                `id`, `name`, `position` (an `(x, y, z)`
                                light-year tuple in `nav_between`'s scope
                                frame -- `origin_position`/
                                `destination_position`/a route's
                                `positions` entries, see `queryDb.
                                nav_between`), and `role` (`"origin"`,
                                `"destination"`, or `"hop"`). Always at
                                least the origin and destination; any
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
    for waypoint, (svg_x, svg_y) in zip(waypoints, projected):
        if waypoint["role"] == "origin":
            points_html.append(_point_html(db_name, waypoint, svg_x, svg_y, "navmap-origin", _ORIGIN_R))
        elif waypoint["role"] == "destination":
            points_html.append(_point_html(db_name, waypoint, svg_x, svg_y, "navmap-destination", _DEST_R))
        else:
            points_html.append(_point_html(db_name, waypoint, svg_x, svg_y, "navmap-hop", _HOP_R))

    body = (
        f"{direct_line}{route_line}{''.join(points_html)}"
        f"{_compass_html()}{_scale_bar_html(px_per_ly)}"
    )

    svg = (
        f'<svg class="navmap-svg" viewBox="0 0 {_SVG_SIZE:.0f} {_SVG_SIZE:.0f}" role="img" '
        'aria-label="Top-down plot of the origin, destination, and optimal route in the galactic plane. '
        'Click a point for its system.">'
        f"{body}</svg>"
    )

    route_hint = (
        "solid line &asymp; optimal route via adjacent systems"
        if has_route
        else "no route via adjacent systems was found -- only the direct line is shown"
    )

    return f"""
<section class="panel">
<div class="panel-header">
  <h2>NAV Map</h2>
  <span class="hint">Top-down (galactic X-Y plane, altitude not shown -- see Altitude above) &middot; dashed line &asymp; direct course &middot; {route_hint}</span>
</div>
<div class="navmap-viewport">
{svg}
</div>
</section>
"""
