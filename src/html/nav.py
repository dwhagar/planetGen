#!/usr/bin/env python3
# html/nav.py

"""
NAV page: course, distance, and optimal route between two systems.

Reuses `queryDb.nav_between` (the same function `/api/nav` calls) rather
than re-implementing the availability rules or the course/route math a
third time -- this page is purely a form + rendering shell around it.

Reaching this page: `system.py` links here (`?from=<id>`) whenever that
system has a `sector_id` -- the same first gate `nav_between` itself
checks (see its docstring's availability rules). Without a `to=` yet,
this renders a destination picker instead of a result:

    - A `<select>` of every other system in the origin's own sector --
      always offered when NAV is available at all, and always a small,
      bounded list (a sector "will typically hold only a handful of
      systems", see `spaceSector.py`'s module docstring), so a dropdown
      scales fine here in a way it wouldn't across an entire galaxy.
    - When the origin's sector itself has a galaxy placement, an
      additional plain numeric field for a cross-sector destination's
      system id -- there is no bounded, dropdown-friendly way to offer
      "any system in any galaxy-placed sector" (that set can be
      arbitrarily large), so this asks for an id directly rather than
      pretending a full picker would scale. `system.py`'s own page (and
      search results) is where that id would normally come from.

Once `to=` is present, `nav_between` decides the rest: same-sector vs.
cross-sector scope, or `NavUnavailable` if the two systems turn out not
to support NAV together (rendered here as an inline message on this same
page, a 200 response -- the visitor picked a real, existing system that
just doesn't have a valid route, not a broken URL -- while a `to=` that
isn't a real system id at all is the one case treated as a 404, via
`NotFoundError`, matching every other page's convention for a bad id).
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))
# Falls back to src/ (stellarObjects lives at src/stellarObjects/, src
# layout) so `stellarObjects`/`queryDb` are importable even when this
# package hasn't been `pip install`-ed system-wide -- true for the default
# deployment layout (`html/` and `src/` as siblings under /var/lib/planetGen).
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(_HTML_DIR)), "src"))

from dbutil import NotFoundError, esc, fetch_all, fetch_one, open_readonly, resolve_db_name
from page import query_params, run

from queryDb import NavUnavailable, nav_between


def _destination_form_html(db_name, from_id, sector_id, sector_systems, cross_sector_available):
    """
    Builds the destination-picker form shown before a `to=` is chosen --
    see the module docstring for why same-sector destinations get a
    `<select>` but a cross-sector destination is a typed id instead.
    """
    options = "".join(
        f'<option value="{row["id"]}">{esc(row["name"])}</option>'
        for row in sector_systems
    )
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

    cross_sector_html = ""
    if cross_sector_available:
        cross_sector_html = f"""
<form method="get" action="nav.py" class="search-form">
  <input type="hidden" name="db" value="{esc(db_name)}">
  <input type="hidden" name="from" value="{from_id}">
  <div class="search-fields">
    <label class="search-field">Destination system id (cross-sector)
      <input type="number" name="to" min="1" step="1" placeholder="e.g. 42">
    </label>
  </div>
  <div class="search-actions">
    <button type="submit" class="btn">Plot course</button>
  </div>
</form>
"""

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


def _course_panel_html(db_name, direct, warp_times, scope):
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


def _route_panel_html(conn, db_name, route):
    if route is None:
        return """
<section class="panel">
<h2>Optimal Route</h2>
<p class="hint">No route via adjacent systems could be found between these two systems.</p>
</section>
"""

    path = route["path"]
    names = {}
    if path:
        placeholders = ", ".join("?" for _ in path)
        rows = fetch_all(conn, f"SELECT id, name FROM star_systems WHERE id IN ({placeholders})", tuple(path))
        names = {row["id"]: row["name"] for row in rows}

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


def handler():
    params = query_params()
    db_name = params.get("db", "")
    raw_from_id = params.get("from", "")
    raw_to_id = params.get("to")

    config = resolve_db_name(db_name)
    conn = open_readonly(config)
    try:
        origin = fetch_one(conn, "SELECT * FROM star_systems WHERE id = ?", (raw_from_id,))
        if origin is None:
            raise NotFoundError(f"No such system: {raw_from_id!r}")
        from_id = origin["id"]

        back_html = (
            f'<p class="breadcrumb"><a href="index.py">Databases</a>'
            f' &rarr; <a href="system.py?db={esc(db_name)}&id={from_id}">{esc(origin["name"])}</a>'
            f" &rarr; Nav</p>"
        )

        if origin["sector_id"] is None:
            body = back_html + (
                '<section class="panel"><p class="error">'
                "NAV is not available: this system isn't assigned to a sector.</p></section>"
            )
            return f"Nav: {origin['name']}", body

        if not raw_to_id:
            sector_systems = fetch_all(
                conn, "SELECT id, name FROM star_systems WHERE sector_id = ? AND id != ? ORDER BY name",
                (origin["sector_id"], from_id),
            )
            sector_row = fetch_one(conn, "SELECT center_x_pc FROM sectors WHERE id = ?", (origin["sector_id"],))
            cross_sector_available = sector_row is not None and sector_row["center_x_pc"] is not None

            form_html = _destination_form_html(
                db_name, from_id, origin["sector_id"], sector_systems, cross_sector_available,
            )
            return f"Nav: {origin['name']}", back_html + form_html

        try:
            to_id = int(raw_to_id)
        except ValueError:
            raise NotFoundError(f"No such system: {raw_to_id!r}")

        try:
            result = nav_between(conn, from_id, to_id)
        except ValueError as exc:
            raise NotFoundError(str(exc))
        except NavUnavailable as exc:
            body = back_html + f'<section class="panel"><p class="error">{esc(str(exc))}</p></section>'
            return f"Nav: {origin['name']}", body

        destination = fetch_one(conn, "SELECT name FROM star_systems WHERE id = ?", (to_id,))
        title_html = (
            f'<p class="badges"><span class="badge">To: '
            f'<a href="system.py?db={esc(db_name)}&id={to_id}">{esc(destination["name"])}</a></span></p>'
        )

        course_html = _course_panel_html(db_name, result["direct"]._asdict(), [
            leg._asdict() for leg in result["warp_times"]
        ], result["scope"])
        route_html = _route_panel_html(conn, db_name, result["route"])
    finally:
        conn.close()

    body = back_html + title_html + course_html + route_html
    return f"Nav: {origin['name']}", body


run(handler)
