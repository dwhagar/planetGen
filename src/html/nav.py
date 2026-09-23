#!/usr/bin/env python3
# html/nav.py

"""
NAV page: course, distance, and optimal route between two endpoints --
each either a star system or a standalone phenomenon (nebula, asteroid
field, black hole, or neutron star).

Calls `GET /api/nav` (the same `queryDb.nav_between` the API itself
calls) rather than re-implementing the availability rules or the
course/route math -- this page is purely a form + rendering shell
around it.

Reaching this page: `system.py` links here with both a "Navigate from
here" (posting a `from` id) and a "Navigate to here" (posting a `to` id)
button whenever that system has a `sector_id`; `phenomenon.py` offers the
same pair of buttons unconditionally, posting `from`/`to` plus
`from_kind`/`to_kind="phenomenon"` and `from_type`/`to_type` (see
`lib/page.py`'s `nav_params`/`lib/fmt.py`'s `post_link`, used for every
navigational link in `html/` now, this page's own picker forms included).
The sidenav's own "Nav" link reaches it with neither a `from` nor a `to`
at all instead.

Without a `from` yet (whether or not a `to` is already known -- see
"Navigate to here" above), this renders a sector-then-system picker to
choose one, carrying the already-known `to`/`to_kind`/`to_type` forward
as hidden fields the whole way so picking an origin lands directly on the
course instead of re-prompting for a destination:

    - Step 1 (`from_sector` not yet given): a `<select>` of every sector
      in the database (`GET /api/sectors`, capped at that route's own max
      page size -- see `docs/api.md`'s "Pagination").
    - Step 2 (`from_sector` given): a `<select>` of every system placed
      in that one sector (`GET /api/sectors/<id>`'s own `systems` list)
      -- submitting sets `from` and this function's normal destination-
      picking flow takes over from there, same as arriving via
      `system.py`'s/`phenomenon.py`'s own "Navigate from here" link. A
      phenomenon is never itself offered in this picker (there is no
      bounded, dropdown-friendly way to browse "any phenomenon in the
      database") -- reachable as an origin/destination only via a direct
      link from its own detail page.

Without a `to` yet (origin already known), this renders a destination
picker instead of a result:

    - When the origin is a system: a `<select>` of every other system in
      its own sector -- always offered when NAV is available at all, and
      always a small, bounded list (a sector "will typically hold only a
      handful of systems", see `spaceSector.py`'s module docstring), so a
      dropdown scales fine here in a way it wouldn't across an entire
      galaxy. Plus, when that sector itself has a galaxy placement, the
      same two-step sector-then-system picker the origin itself uses
      above (`to_sector` then `to`) for a cross-sector destination.
    - When the origin is a phenomenon: only the cross-sector picker above
      (unscoped -- a phenomenon has no sector-local peers to offer as a
      same-sector dropdown at all; see `queryDb.
      _load_nav_phenomenon_endpoint`'s own docstring on why).

Once `to` is present, `GET /api/nav` decides the rest: same-sector vs.
cross-sector scope, or a 400 if the two endpoints turn out not to support
NAV together (rendered here as an inline message on this same page, a
200 response -- the visitor picked a real, existing system/phenomenon
that just doesn't have a valid route, not a broken URL -- while a `to`
that isn't a real id at all is the one case treated as a 404, matching
every other page's convention for a bad id).
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import ApiError, NotFoundError, get_nav, get_phenomenon, get_sector, get_sectors, get_system
from fmt import esc, post_link
from navmap import render_nav_map_panel
from page import nav_params, run

_SECTOR_PICKER_LIMIT = 500
"""int: Caps how many sectors the picker's `<select>` offers -- matches
`GET /api/sectors`'s own max page size (`docs/api.md`'s "Pagination"), so
this never silently asks the API for more than it would ever return."""


def _phenomenon_node_key(phenomenon_type, phenomenon_id):
    """
    Mirrors `queryDb._phenomenon_nav_key`'s exact string format -- the id
    a phenomenon endpoint uses in `GET /api/nav`'s `route.path`/
    `route.positions`, now that it's round-tripped through JSON as a
    plain string rather than the private tuple/string this module has no
    direct access to (`html/` only ever talks to the database through the
    API -- see `docs/html-interface.md`). Duplicated here as the one
    place this module needs to build/recognize that same format itself,
    rather than a shared import across the API boundary.
    """
    return f"phenomenon:{phenomenon_type}:{phenomenon_id}"


def _is_phenomenon_node(node_id):
    return isinstance(node_id, str) and node_id.startswith("phenomenon:")


def _parse_phenomenon_node(node_id):
    _, phenomenon_type, phenomenon_id = node_id.split(":", 2)
    return phenomenon_type, int(phenomenon_id)


def _resolve_endpoint(db_name, entity_id, kind, phenomenon_type):
    """
    Resolves one NAV endpoint's display info -- a system (the default) or
    a standalone phenomenon -- into one common shape the rest of this
    module works from without branching on kind everywhere.

    Args:
        db_name (str): The current `?db=` value.
        entity_id: The raw `star_systems.id` or phenomenon-table id, as
            received from a form/query param (may still be a string).
        kind (str): `"system"` or `"phenomenon"`.
        phenomenon_type (str or None): One of `queryDb.
            _PHENOMENON_TYPE_TO_TABLE`'s keys when `kind ==
            "phenomenon"`, else ignored.

    Returns:
        dict: `id` (int, the resolved row id), `kind`, `type` (the
            phenomenon type, or `None` for a system), `name`, `link_html`
            (a ready link back to `system.py`/`phenomenon.py`),
            `node_key` (the id `nav_between`'s own route/positions dicts
            use for this endpoint -- the raw int for a system,
            `_phenomenon_node_key`'s string for a phenomenon),
            `sector_id` (a system's own, or `None` for a phenomenon --
            it has no sector-local frame to offer here), and
            `has_galaxy_position` (phenomenon-only; `None` for a system,
            where sector placement already covers this).

    Raises:
        NotFoundError: If no such system/phenomenon exists (via
            `apiclient`'s own 404 handling).
    """
    if kind == "phenomenon":
        detail = get_phenomenon(db_name, phenomenon_type, entity_id)
        return {
            "id": detail["id"], "kind": "phenomenon", "type": phenomenon_type,
            "name": detail["name"],
            "link_html": post_link(
                "phenomenon.py", {"db": db_name, "type": phenomenon_type, "id": detail["id"]}, esc(detail["name"]),
            ),
            "node_key": _phenomenon_node_key(phenomenon_type, detail["id"]),
            "sector_id": None,
            "has_galaxy_position": detail.get("galactic_radius_pc") is not None,
        }
    system = get_system(db_name, entity_id)
    return {
        "id": system["id"], "kind": "system", "type": None,
        "name": system["name"],
        "link_html": post_link("system.py", {"db": db_name, "id": system["id"]}, esc(system["name"])),
        "node_key": system["id"],
        "sector_id": system["sector_id"],
        "has_galaxy_position": None,
    }


def _endpoint_hidden_fields(prefix, entity_id, kind, phenomenon_type):
    """
    Hidden-field dict (for `_picker_form_html`'s own `hidden` param, or
    direct interpolation elsewhere) carrying an already-known NAV
    endpoint forward through a picker form -- just `{prefix: id}` for a
    system (unchanged from before phenomenon support existed), plus
    `{prefix}_kind`/`{prefix}_type` when it's a phenomenon.
    """
    fields = {prefix: entity_id}
    if kind == "phenomenon":
        fields[f"{prefix}_kind"] = "phenomenon"
        fields[f"{prefix}_type"] = phenomenon_type
    return fields


def _picker_form_html(db_name, options_html, field_name, label, hidden, heading, back_params=None):
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
            fields (e.g. an already-chosen `from_sector`/`from`, or a
            pre-set `to`/`to_kind`/`to_type` from a "Navigate to here"
            link).
        heading (str): This panel's `<h2>` text.
        back_params (dict, optional): Params for a "start over" link shown
            above the form (omitted for the very first step, which has
            nothing to go back to).

    Returns:
        str: A complete `<section class="panel">` block.
    """
    hidden_html = "".join(
        f'<input type="hidden" name="{esc(str(name))}" value="{esc(str(value))}">'
        for name, value in hidden.items()
    )
    back_html = (
        f'<p class="hint">{post_link("nav.py", back_params, "&larr; Start over")}</p>' if back_params else ""
    )
    return f"""
{back_html}
<section class="panel">
<h2>{esc(heading)}</h2>
<form method="post" action="nav.py" class="search-form">
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


def _origin_picker_html(db_name, from_sector_raw, extra_hidden=None):
    """
    The two-step "choose a starting sector, then a starting system"
    picker shown whenever `nav.py` is reached with no `from=` at all
    (the sidenav's own "Nav" link, or a "Navigate to here" link that only
    set `to=`) -- see the module docstring. `extra_hidden` (typically an
    already-known `to`/`to_kind`/`to_type`, from `_endpoint_hidden_fields`)
    is threaded through every step and every "start over" link so it's
    never lost while choosing an origin.
    """
    extra_hidden = extra_hidden or {}
    if not from_sector_raw:
        sectors = get_sectors(db_name, limit=_SECTOR_PICKER_LIMIT)["items"]
        if not sectors:
            return '<section class="panel"><p class="hint">No sectors have been generated in this database yet.</p></section>'
        return _picker_form_html(
            db_name, _sector_options_html(sectors), "from_sector", "Sector", dict(extra_hidden),
            "Choose a starting sector",
        )

    try:
        from_sector_id = int(from_sector_raw)
    except ValueError:
        raise NotFoundError(f"No such sector: {from_sector_raw!r}")
    sector = get_sector(db_name, from_sector_id)
    systems = sector["systems"]
    start_over_params = {"db": db_name, **extra_hidden}
    if not systems:
        start_over = post_link("nav.py", start_over_params, "&larr; Start over")
        return (
            f'<p class="hint">{start_over}</p>'
            '<section class="panel"><p class="hint">No systems are placed in this sector yet.</p></section>'
        )
    return _picker_form_html(
        db_name, _system_options_html(systems), "from", "Starting system",
        {"from_sector": from_sector_id, **extra_hidden}, f"Choose a starting system in {sector['name']}",
        back_params=start_over_params,
    )


def _cross_sector_destination_html(db_name, from_hidden, origin_sector_id, to_sector_raw):
    """
    The cross-sector half of the destination picker: the same two-step
    sector-then-system cascade `_origin_picker_html` uses, scoped to
    `?to_sector=`/`to=`. Excludes the origin's own sector when it has one
    (already covered by the same-sector `<select>` alongside this one) --
    a phenomenon origin (`origin_sector_id` `None`) excludes nothing,
    since it has no sector of its own to exclude.
    """
    if not to_sector_raw:
        sectors = get_sectors(db_name, limit=_SECTOR_PICKER_LIMIT)["items"]
        if origin_sector_id is not None:
            sectors = [row for row in sectors if row["id"] != origin_sector_id]
        if not sectors:
            return '<p class="hint">No other galaxy-placed sectors exist in this database yet.</p>'
        return _picker_form_html(
            db_name, _sector_options_html(sectors), "to_sector", "Sector", dict(from_hidden),
            "Choose a destination sector",
        )

    try:
        to_sector_id = int(to_sector_raw)
    except ValueError:
        raise NotFoundError(f"No such sector: {to_sector_raw!r}")
    sector = get_sector(db_name, to_sector_id)
    systems = sector["systems"]
    back_params = {"db": db_name, **from_hidden}
    if not systems:
        start_over = post_link("nav.py", back_params, "&larr; Start over")
        return (
            f'<p class="hint">{start_over}</p>'
            '<p class="hint">No systems are placed in this sector yet.</p>'
        )
    return _picker_form_html(
        db_name, _system_options_html(systems), "to", "Destination system",
        {**from_hidden, "to_sector": to_sector_id}, f"Choose a destination system in {sector['name']}",
        back_params=back_params,
    )


def _destination_form_html(db_name, from_hidden, origin_sector_id, sector_systems, cross_sector_available, to_sector_raw):
    """
    Builds the destination-picker form shown before a `to=` is chosen --
    see the module docstring for the same-sector `<select>` (system
    origins only) vs. the cross-sector sector-then-system cascade (every
    origin).
    """
    same_sector_html = ""
    if origin_sector_id is not None:
        options = _system_options_html(sector_systems)
        if options:
            hidden_html = "".join(
                f'<input type="hidden" name="{esc(str(name))}" value="{esc(str(value))}">'
                for name, value in from_hidden.items()
            )
            same_sector_html = f"""
<form method="post" action="nav.py" class="search-form">
  <input type="hidden" name="db" value="{esc(db_name)}">
  {hidden_html}
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
        _cross_sector_destination_html(db_name, from_hidden, origin_sector_id, to_sector_raw)
        if cross_sector_available else ""
    )

    sections = []
    if origin_sector_id is not None:
        sections.append(f"""
<section class="panel">
<h2>Same-sector destination</h2>
{same_sector_html}
</section>
""")
        cross_heading = "Cross-sector destination"
        cross_fallback = '<p class="hint">This sector has no galaxy placement, so NAV is only available to other systems in this same sector.</p>'
    else:
        cross_heading = "Destination"
        cross_fallback = '<p class="hint">NAV is not available from here.</p>'

    sections.append(f"""
<section class="panel">
<h2>{cross_heading}</h2>
{cross_sector_html or cross_fallback}
</section>
""")
    return "".join(sections)


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
    `_route_panel_html` (the hop-by-hop list) and `_nav_map_waypoints`
    (each waypoint's label/tooltip). `known_names` (`{node_key: name}`)
    is consulted first, since the origin/destination are always already
    known by the time this runs -- only genuinely new intermediate hop
    ids need their own `GET /api/systems/<id>` lookup (a route's hop
    count is small, bounded by the k-nearest-neighbor adjacency graph
    `queryDb.nav_between` builds it from, so this stays a handful of
    requests at most).

    Args:
        db_name (str): The current `?db=` value.
        route (dict or None): `GET /api/nav`'s `route`.
        known_names (dict): `{node_key: name}` already available (the
            origin and destination), consulted before falling back to
            an API lookup for an intermediate hop.

    Returns:
        dict: `{node_key: name}`, empty if `route` is `None`.
    """
    if route is None or not route["path"]:
        return {}
    names = dict(known_names)
    for node_id in route["path"]:
        if node_id in names or _is_phenomenon_node(node_id):
            # A phenomenon node is only ever path[0]/path[-1] (see
            # queryDb.nav_between's own docstring), already present in
            # known_names -- unreachable in practice, guarded rather than
            # assumed so get_system() below never misinterprets a
            # phenomenon's own composite string key as a system id.
            continue
        names[node_id] = get_system(db_name, node_id)["name"]
    return names


def _route_panel_html(db_name, route, names):
    if route is None:
        return """
<section class="panel">
<h2>Optimal Route</h2>
<p class="hint">No route via adjacent systems could be found between these two endpoints.</p>
</section>
"""

    path = route["path"]

    def hop_link(node_id):
        label = esc(names.get(node_id, str(node_id)))
        if _is_phenomenon_node(node_id):
            phenomenon_type, phenomenon_id = _parse_phenomenon_node(node_id)
            return post_link("phenomenon.py", {"db": db_name, "type": phenomenon_type, "id": phenomenon_id}, label)
        return post_link("system.py", {"db": db_name, "id": node_id}, label)

    hops = "".join(f"<li>{hop_link(node_id)}</li>" for node_id in path)
    return f"""
<section class="panel">
<h2>Optimal Route</h2>
<p>{len(path)} stop(s), {route['distance_ly']:,.2f} ly total.</p>
<ol class="nav-route">{hops}</ol>
</section>
"""


def _nav_map_waypoints(origin, destination, origin_position, destination_position, route, names):
    """
    Builds the ordered origin-to-destination waypoint list
    `render_nav_map_panel` plots: the origin, then any route hops (in
    path order, always systems -- see `queryDb.nav_between`'s docstring)
    between the two endpoints, then the destination.
    """
    waypoints = [{
        "id": origin["id"], "kind": origin["kind"], "type": origin["type"],
        "name": origin["name"], "position": origin_position, "role": "origin",
    }]
    if route is not None:
        for node_id in route["path"][1:-1]:
            waypoints.append({
                "id": node_id, "kind": "system", "type": None,
                "name": names.get(node_id, str(node_id)),
                "position": route["positions"][str(node_id)],
                "role": "hop",
            })
    waypoints.append({
        "id": destination["id"], "kind": destination["kind"], "type": destination["type"],
        "name": destination["name"], "position": destination_position, "role": "destination",
    })
    return waypoints


def handler():
    params = nav_params()
    db_name = params.get("db", "")

    raw_from_id = params.get("from", "")
    from_kind = params.get("from_kind") or "system"
    from_type = params.get("from_type") or None

    raw_to_id = params.get("to", "")
    to_kind = params.get("to_kind") or "system"
    to_type = params.get("to_type") or None

    if not raw_from_id:
        # No origin yet -- show the origin picker. A "Navigate to here"
        # link (from system.py/phenomenon.py) may have already set `to`;
        # carry it forward as hidden fields so choosing an origin lands
        # directly on the course instead of re-prompting for a destination.
        to_hidden = _endpoint_hidden_fields("to", raw_to_id, to_kind, to_type) if raw_to_id else {}
        body = '<p class="breadcrumb">Nav</p>' + _origin_picker_html(db_name, params.get("from_sector", ""), to_hidden)
        return "Nav", body

    origin = _resolve_endpoint(db_name, raw_from_id, from_kind, from_type)
    from_id = origin["id"]

    back_html = f'<p class="breadcrumb">{origin["link_html"]} &rarr; Nav</p>'

    if origin["kind"] == "phenomenon" and not origin["has_galaxy_position"]:
        body = back_html + (
            '<section class="panel"><p class="error">'
            "NAV is not available: this phenomenon has not been placed in the galaxy.</p></section>"
        )
        return f"Nav: {origin['name']}", body

    if origin["kind"] == "system" and origin["sector_id"] is None:
        body = back_html + (
            '<section class="panel"><p class="error">'
            "NAV is not available: this system isn't assigned to a sector.</p></section>"
        )
        return f"Nav: {origin['name']}", body

    if not raw_to_id:
        from_hidden = _endpoint_hidden_fields("from", from_id, origin["kind"], origin["type"])
        if origin["kind"] == "system":
            sector = get_sector(db_name, origin["sector_id"])
            sector_systems = [row for row in sector["systems"] if row["id"] != from_id]
            cross_sector_available = sector["placed"]
            origin_sector_id = origin["sector_id"]
        else:
            sector_systems = []
            cross_sector_available = True
            origin_sector_id = None

        form_html = _destination_form_html(
            db_name, from_hidden, origin_sector_id, sector_systems, cross_sector_available,
            params.get("to_sector", ""),
        )
        return f"Nav: {origin['name']}", back_html + form_html

    try:
        to_id = int(raw_to_id)
    except ValueError:
        raise NotFoundError(f"No such {'phenomenon' if to_kind == 'phenomenon' else 'system'}: {raw_to_id!r}")

    try:
        result = get_nav(
            db_name, from_id, to_id,
            from_kind=origin["kind"], to_kind=to_kind, from_type=origin["type"], to_type=to_type,
        )
    except ApiError as exc:
        # GET /api/nav responds 400 (wrapped as ApiError, not
        # NotFoundError) when the two endpoints exist but don't support
        # NAV together -- e.g. different, non-galaxy-placed sectors. A
        # genuinely unknown `to_id` is a 404, raised as NotFoundError
        # instead, and lets `run()`'s own handling take over.
        body = back_html + f'<section class="panel"><p class="error">{esc(str(exc))}</p></section>'
        return f"Nav: {origin['name']}", body

    destination = _resolve_endpoint(db_name, to_id, to_kind, to_type)
    title_html = f'<p class="badges"><span class="badge">To: {destination["link_html"]}</span></p>'

    course_html = _course_panel_html(result["direct"], result["warp_times"], result["scope"])

    known_names = {origin["node_key"]: origin["name"], destination["node_key"]: destination["name"]}
    route_names = _route_names(db_name, result["route"], known_names)
    route_html = _route_panel_html(db_name, result["route"], route_names)

    waypoints = _nav_map_waypoints(
        origin, destination, result["origin_position"], result["destination_position"],
        result["route"], route_names,
    )
    map_html = render_nav_map_panel(db_name, waypoints, has_route=result["route"] is not None)

    body = back_html + title_html + course_html + map_html + route_html
    return f"Nav: {origin['name']}", body


run(handler)
