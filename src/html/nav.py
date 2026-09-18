#!/usr/bin/env python3
# html/nav.py

"""
NAV page: course, distance, and optimal route between two systems.

Calls `GET /api/nav` (the same `queryDb.nav_between` the API itself
calls) rather than re-implementing the availability rules or the
course/route math -- this page is purely a form + rendering shell
around it.

Reaching this page: `system.py` links here (`?from=<id>`) whenever that
system has a `sector_id` -- the same first gate `nav_between` itself
checks (see its docstring's availability rules); the sidenav's own "Nav"
link reaches it with no `from=` at all instead. Without a `from=` yet,
this renders a sector-then-system picker to choose one:

    - Step 1 (`?from_sector=<id>` not yet given): a `<select>` of every
      sector in the database (`GET /api/sectors`, capped at that route's
      own max page size -- see `docs/api.md`'s "Pagination").
    - Step 2 (`?from_sector=<id>` given): a `<select>` of every system
      placed in that one sector (`GET /api/sectors/<id>`'s own `systems`
      list) -- submitting sets `from=` and this function's normal
      destination-picking flow takes over from there, same as arriving
      via `system.py`'s link.

Without a `to=` yet (origin already known), this renders a destination
picker instead of a result:

    - A `<select>` of every other system in the origin's own sector --
      always offered when NAV is available at all, and always a small,
      bounded list (a sector "will typically hold only a handful of
      systems", see `spaceSector.py`'s module docstring), so a dropdown
      scales fine here in a way it wouldn't across an entire galaxy.
    - When the origin's sector itself has a galaxy placement, the same
      two-step sector-then-system picker the origin itself uses above
      (`?to_sector=<id>` then `to=<id>`) -- there is no bounded,
      dropdown-friendly way to offer "any system in any galaxy-placed
      sector" in one list (that set can be arbitrarily large), so this
      narrows to one sector first, exactly like choosing the origin does.

Once `to=` is present, `GET /api/nav` decides the rest: same-sector vs.
cross-sector scope, or a 400 if the two systems turn out not to support
NAV together (rendered here as an inline message on this same page, a
200 response -- the visitor picked a real, existing system that just
doesn't have a valid route, not a broken URL -- while a `to=` that isn't
a real system id at all is the one case treated as a 404, matching every
other page's convention for a bad id).
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import ApiError, NotFoundError, get_nav, get_sector, get_sectors, get_system
from fmt import esc
from navmap import render_nav_map_panel
from page import query_params, run

_SECTOR_PICKER_LIMIT = 500
"""int: Caps how many sectors the picker's `<select>` offers -- matches
`GET /api/sectors`'s own max page size (`docs/api.md`'s "Pagination"), so
this never silently asks the API for more than it would ever return."""


def _picker_form_html(db_name, options_html, field_name, label, hidden, heading, back_href=None):
    """
    One `<select>` + submit form, shared by every step of the sector-then-
    system picker (the origin picker, reached with no `from=` at all, and
    the cross-sector destination picker) -- see the module docstring.

    Args:
        db_name (str): The current `?db=` value.
        options_html (str): Pre-built `<option>` tags.
        field_name (str): The `<select>`'s own `name` (`from_sector`,
            `from`, `to_sector`, or `to`).
        label (str): The `<select>`'s visible label.
        hidden (dict): Extra `name: value` pairs carried forward as hidden
            fields (e.g. an already-chosen `from_sector`/`from`).
        heading (str): This panel's `<h2>` text.
        back_href (str, optional): A "start over" link shown above the
            form (omitted for the very first step, which has nothing to
            go back to).

    Returns:
        str: A complete `<section class="panel">` block.
    """
    hidden_html = "".join(
        f'<input type="hidden" name="{esc(str(name))}" value="{esc(str(value))}">'
        for name, value in hidden.items()
    )
    back_html = f'<p class="hint"><a href="{esc(back_href)}">&larr; Start over</a></p>' if back_href else ""
    return f"""
{back_html}
<section class="panel">
<h2>{esc(heading)}</h2>
<form method="get" action="nav.py" class="search-form">
  <input type="hidden" name="db" value="{esc(db_name)}">
  {hidden_html}
  <div class="search-fields">
    <label class="search-field">{esc(label)}
      <select name="{esc(field_name)}">{options_html}</select>
    </label>
  </div>
  <div class="search-actions">
    <button type="submit" class="btn">Continue</button>
  </div>
</form>
</section>
"""


def _sector_options_html(sectors):
    return "".join(f'<option value="{row["id"]}">{esc(row["name"])}</option>' for row in sectors)


def _system_options_html(systems):
    return "".join(f'<option value="{row["id"]}">{esc(row["name"])}</option>' for row in systems)


def _origin_picker_html(db_name, from_sector_raw):
    """
    The two-step "choose a starting sector, then a starting system"
    picker shown when `nav.py` is reached with no `from=` at all (the
    sidenav's own "Nav" link) -- see the module docstring.
    """
    if not from_sector_raw:
        sectors = get_sectors(db_name, limit=_SECTOR_PICKER_LIMIT)["items"]
        if not sectors:
            return '<section class="panel"><p class="hint">No sectors have been generated in this database yet.</p></section>'
        return _picker_form_html(
            db_name, _sector_options_html(sectors), "from_sector", "Sector", {},
            "Choose a starting sector",
        )

    try:
        from_sector_id = int(from_sector_raw)
    except ValueError:
        raise NotFoundError(f"No such sector: {from_sector_raw!r}")
    sector = get_sector(db_name, from_sector_id)
    systems = sector["systems"]
    if not systems:
        return (
            f'<p class="hint"><a href="nav.py?db={esc(db_name)}">&larr; Start over</a></p>'
            '<section class="panel"><p class="hint">No systems are placed in this sector yet.</p></section>'
        )
    return _picker_form_html(
        db_name, _system_options_html(systems), "from", "Starting system",
        {"from_sector": from_sector_id}, f"Choose a starting system in {sector['name']}",
        back_href=f"nav.py?db={esc(db_name)}",
    )


def _cross_sector_destination_html(db_name, from_id, origin_sector_id, to_sector_raw):
    """
    The cross-sector half of the destination picker: the same two-step
    sector-then-system cascade `_origin_picker_html` uses, scoped to
    `?to_sector=`/`to=` and excluding the origin's own sector (already
    covered by the same-sector `<select>` alongside this one).
    """
    if not to_sector_raw:
        sectors = [
            row for row in get_sectors(db_name, limit=_SECTOR_PICKER_LIMIT)["items"]
            if row["id"] != origin_sector_id
        ]
        if not sectors:
            return '<p class="hint">No other galaxy-placed sectors exist in this database yet.</p>'
        return _picker_form_html(
            db_name, _sector_options_html(sectors), "to_sector", "Sector", {"from": from_id},
            "Choose a destination sector",
        )

    try:
        to_sector_id = int(to_sector_raw)
    except ValueError:
        raise NotFoundError(f"No such sector: {to_sector_raw!r}")
    sector = get_sector(db_name, to_sector_id)
    systems = sector["systems"]
    back_href = f"nav.py?db={esc(db_name)}&from={from_id}"
    if not systems:
        return (
            f'<p class="hint"><a href="{back_href}">&larr; Start over</a></p>'
            '<p class="hint">No systems are placed in this sector yet.</p>'
        )
    return _picker_form_html(
        db_name, _system_options_html(systems), "to", "Destination system",
        {"from": from_id, "to_sector": to_sector_id}, f"Choose a destination system in {sector['name']}",
        back_href=back_href,
    )


def _destination_form_html(db_name, from_id, origin_sector_id, sector_systems, cross_sector_available, to_sector_raw):
    """
    Builds the destination-picker form shown before a `to=` is chosen --
    see the module docstring for the same-sector `<select>` vs. the
    cross-sector sector-then-system cascade.
    """
    options = _system_options_html(sector_systems)
    same_sector_html = ""
    if options:
        same_sector_html = f"""
<form method="get" action="nav.py" class="search-form">
  <input type="hidden" name="db" value="{esc(db_name)}">
  <input type="hidden" name="from" value="{from_id}">
  <div class="search-fields">
    <label class="search-field">Destination (same sector)
      <select name="to">{options}</select>
    </label>
  </div>
  <div class="search-actions">
    <button type="submit" class="btn">Plot course</button>
  </div>
</form>
"""
    else:
        same_sector_html = '<p class="hint">No other systems are placed in this sector yet.</p>'

    cross_sector_html = (
        _cross_sector_destination_html(db_name, from_id, origin_sector_id, to_sector_raw)
        if cross_sector_available else ""
    )

    return f"""
<section class="panel">
<h2>Same-sector destination</h2>
{same_sector_html}
</section>
<section class="panel">
<h2>Cross-sector destination</h2>
{cross_sector_html or '<p class="hint">This sector has no galaxy placement, so NAV is only available to other systems in this same sector.</p>'}
</section>
"""


def _course_panel_html(direct, warp_times, scope):
    rows = "".join(
        f"<tr><td>Warp {leg['warp_factor']}</td><td>{leg['velocity_multiple_of_c']:,.2f}&times; c</td>"
        f"<td>{esc(leg['formatted'])}</td></tr>"
        for leg in warp_times
    )
    scope_label = "Same sector" if scope == "sector" else "Cross-sector (galaxy)"
    return f"""
<section class="panel">
<h2>Direct Course</h2>
<p class="badges"><span class="badge">{esc(scope_label)}</span></p>
<div class="table-scroll"><table>
  <tbody>
    <tr><td>Distance</td><td>{direct['distance_ly']:,.2f} ly</td></tr>
    <tr><td>Azimuth</td><td>{direct['azimuth_deg']:.2f}&deg;</td></tr>
    <tr><td>Altitude</td><td>{direct['altitude_deg']:+.2f}&deg;</td></tr>
  </tbody>
</table></div>
<h3>Travel Time</h3>
<div class="table-scroll"><table>
  <thead><tr><th>Warp Factor</th><th>Speed</th><th>Travel Time</th></tr></thead>
  <tbody>{rows}</tbody>
</table></div>
</section>
"""


def _route_names(db_name, route, known_names):
    """
    Looks up every route hop's display name -- shared by
    `_route_panel_html` (the hop-by-hop list) and `render_nav_map_panel`
    (each waypoint's label/tooltip). `known_names` (`{id: name}`) is
    consulted first, since the origin/destination are always already
    known by the time this runs -- only genuinely new intermediate hop
    ids need their own `GET /api/systems/<id>` lookup (a route's hop
    count is small, bounded by the k-nearest-neighbor adjacency graph
    `queryDb.nav_between` builds it from, so this stays a handful of
    requests at most).

    Args:
        db_name (str): The current `?db=` value.
        route (dict or None): `GET /api/nav`'s `route`.
        known_names (dict): `{system_id: name}` already available (the
            origin and destination), consulted before falling back to
            an API lookup for an intermediate hop.

    Returns:
        dict: `{star_systems.id: name}`, empty if `route` is `None`.
    """
    if route is None or not route["path"]:
        return {}
    names = dict(known_names)
    for system_id in route["path"]:
        if system_id not in names:
            names[system_id] = get_system(db_name, system_id)["name"]
    return names


def _route_panel_html(db_name, route, names):
    if route is None:
        return """
<section class="panel">
<h2>Optimal Route</h2>
<p class="hint">No route via adjacent systems could be found between these two systems.</p>
</section>
"""

    path = route["path"]
    hops = "".join(
        f'<li><a href="system.py?db={esc(db_name)}&id={system_id}">{esc(names.get(system_id, system_id))}</a></li>'
        for system_id in path
    )
    return f"""
<section class="panel">
<h2>Optimal Route</h2>
<p>{len(path)} system(s), {route['distance_ly']:,.2f} ly total.</p>
<ol class="nav-route">{hops}</ol>
</section>
"""


def _nav_map_waypoints(from_id, to_id, origin_name, destination_name, origin_position, destination_position, route, names):
    """
    Builds the ordered origin-to-destination waypoint list
    `render_nav_map_panel` plots: the origin, then any route hops (in path
    order) between the two endpoints, then the destination.
    """
    waypoints = [{"id": from_id, "name": origin_name, "position": origin_position, "role": "origin"}]
    if route is not None:
        for system_id in route["path"][1:-1]:
            waypoints.append({
                "id": system_id,
                "name": names.get(system_id, str(system_id)),
                "position": route["positions"][str(system_id)],
                "role": "hop",
            })
    waypoints.append({"id": to_id, "name": destination_name, "position": destination_position, "role": "destination"})
    return waypoints


def handler():
    params = query_params()
    db_name = params.get("db", "")
    raw_from_id = params.get("from", "")
    raw_to_id = params.get("to")

    if not raw_from_id:
        body = '<p class="breadcrumb">Nav</p>' + _origin_picker_html(db_name, params.get("from_sector", ""))
        return "Nav", body

    origin = get_system(db_name, raw_from_id)
    from_id = origin["id"]

    back_html = (
        f'<p class="breadcrumb"><a href="system.py?db={esc(db_name)}&id={from_id}">{esc(origin["name"])}</a>'
        f" &rarr; Nav</p>"
    )

    if origin["sector_id"] is None:
        body = back_html + (
            '<section class="panel"><p class="error">'
            "NAV is not available: this system isn't assigned to a sector.</p></section>"
        )
        return f"Nav: {origin['name']}", body

    if not raw_to_id:
        sector = get_sector(db_name, origin["sector_id"])
        sector_systems = [row for row in sector["systems"] if row["id"] != from_id]
        cross_sector_available = sector["placed"]

        form_html = _destination_form_html(
            db_name, from_id, origin["sector_id"], sector_systems, cross_sector_available,
            params.get("to_sector", ""),
        )
        return f"Nav: {origin['name']}", back_html + form_html

    try:
        to_id = int(raw_to_id)
    except ValueError:
        raise NotFoundError(f"No such system: {raw_to_id!r}")

    try:
        result = get_nav(db_name, from_id, to_id)
    except ApiError as exc:
        # GET /api/nav responds 400 (wrapped as ApiError, not
        # NotFoundError) when the two systems exist but don't support
        # NAV together -- e.g. different, non-galaxy-placed sectors. A
        # genuinely unknown `to_id` is a 404, raised as NotFoundError
        # instead, and lets `run()`'s own handling take over.
        body = back_html + f'<section class="panel"><p class="error">{esc(str(exc))}</p></section>'
        return f"Nav: {origin['name']}", body

    destination = get_system(db_name, to_id)
    title_html = (
        f'<p class="badges"><span class="badge">To: '
        f'<a href="system.py?db={esc(db_name)}&id={to_id}">{esc(destination["name"])}</a></span></p>'
    )

    course_html = _course_panel_html(result["direct"], result["warp_times"], result["scope"])

    known_names = {from_id: origin["name"], to_id: destination["name"]}
    route_names = _route_names(db_name, result["route"], known_names)
    route_html = _route_panel_html(db_name, result["route"], route_names)

    waypoints = _nav_map_waypoints(
        from_id, to_id, origin["name"], destination["name"],
        result["origin_position"], result["destination_position"],
        result["route"], route_names,
    )
    map_html = render_nav_map_panel(db_name, waypoints, has_route=result["route"] is not None)

    body = back_html + title_html + course_html + map_html + route_html
    return f"Nav: {origin['name']}", body


run(handler)
