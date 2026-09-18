#!/usr/bin/env python3
# html/sector.py

"""
Sector detail page: the sector's name/size and every system placed in it,
with quadrant and star-type info, linking to `system.py` for each -- plus
an interactive 3D "Sector Map" (see `lib/starmap.py`) of the same systems
plotted by position within the sector, outlined by the sector's real
on-shell wedge shape when it has a galaxy placement (else a plain cube),
plus a translucent cloud for every nebula/asteroid field (and a point
marker for every black hole/neutron star) whose real galaxy-frame sphere
reaches into this sector's own cube (`queryDb.phenomena_near_sector`, via
`GET /api/sectors/<id>`'s `phenomena` key -- see `schema.sql`'s
"v18"/"v21" header notes). That same phenomena list also gets its own
table below the systems one (mirroring `phenomena.py`'s flat listing,
scoped to just this sector's own neighborhood), each row linking to
`phenomenon.py` -- omitted entirely when nothing nearby qualifies.

Once this sector has a wiki page (`sectors.wiki_url` -- either uploaded
from here or set directly via `html/admin.py`'s manual-link admin
section, see `schema.sql`'s "v22" header note), a link to it is shown
(opening in a new tab). An admin session additionally gets an
"Upload to Wiki" form (shown only while `wiki_url` is unset) offering
whichever backend(s) are configured deployment-wide (`GET
/api/wiki-config`) -- unlike a system, a sector has no persisted
generated page of its own, so the page content is built fresh from this
sector's own detail at upload time (see `routes.py`'s
`_sector_wiki_content`).
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import ApiError, auth_me, get_sector, get_wiki_config, upload_sector_to_wiki
from fmt import esc, linkify_location
from galaxymap import sector_quadrant
from page import form_params, incoming_cookie_header, query_params, run
from starmap import render_map_panel

_PHENOMENON_TYPE_LABELS = {
    "nebula": "Nebula", "asteroid_field": "Asteroid Field",
    "black_hole": "Black Hole", "neutron_star": "Neutron Star",
}
"""dict: Same display labels `phenomena.py`'s own flat listing uses, for
`queryDb.phenomena_near_sector`'s `type` values."""


def _wiki_section_html(db_name, sector_id, sector, wiki_config, wiki_message, wiki_error):
    """
    Builds the sector's wiki section: a link to `sector["wiki_url"]`
    (opening in a new tab) once it's set -- from an upload here or a
    manual admin edit alike, this page doesn't distinguish which -- or,
    while unset, an "Upload to Wiki" form offering whichever backend(s)
    `wiki_config` reports configured. Returns just the message/error (no
    link, no form) if `wiki_url` is unset and no backend is configured.
    """
    message_html = f'<p class="hint">{esc(wiki_message)}</p>' if wiki_message else ""
    error_html = f'<p class="error">{esc(wiki_error)}</p>' if wiki_error else ""

    if sector["wiki_url"]:
        return f"""
{message_html}{error_html}
<p class="wiki-link"><a href="{esc(sector['wiki_url'])}" target="_blank" rel="noopener noreferrer">View on Wiki</a></p>
"""

    options = [name for name, configured in (("wikijs", wiki_config.get("wikijs")),
                                              ("mediawiki", wiki_config.get("mediawiki"))) if configured]
    if not options:
        return f"{message_html}{error_html}"

    labels = {"wikijs": "Wiki.js", "mediawiki": "MediaWiki"}
    radios = " ".join(
        f'<label><input type="radio" name="backend" value="{value}"{" checked" if i == 0 else ""}> '
        f'{labels[value]}</label>'
        for i, value in enumerate(options)
    )
    base_url = f"sector.py?db={esc(db_name)}&id={sector_id}"
    return f"""
{message_html}{error_html}
<section class="panel">
<h2>Upload to Wiki</h2>
<form method="post" action="{base_url}" class="search-form">
  <input type="hidden" name="action" value="upload_wiki">
  <div class="search-fields">
    <div class="search-field">{radios}</div>
    <label class="search-field">Path (Wiki.js only -- MediaWiki uses this sector's name)
      <input type="text" name="path" placeholder="e.g. sectors/{esc(sector['name'])}">
    </label>
  </div>
  <div class="search-actions">
    <button type="submit" class="btn">Upload</button>
  </div>
</form>
</section>
"""


def handler():
    params = query_params()
    db_name = params.get("db", "")
    sector_id = params.get("id", "")

    sector = get_sector(db_name, sector_id)
    systems = sector["systems"]
    name_to_id = {row["name"]: row["id"] for row in systems}

    rows = []
    map_systems = []
    for row in systems:
        # binary_type only ever describes a 'close' (P-type) pair's merged
        # effective star -- NULL for a single star and for a 'wide' (S-type)
        # pair (no merged star exists there; see schema.sql's "v15" note),
        # so both fall back to joining each of row["stars"]'s own types
        # instead of showing a blank cell.
        if row["is_binary"] and row.get("binary_type"):
            star_type = row["binary_type"]
        else:
            star_type = " / ".join(star["star_type"] for star in row["stars"]) if row["stars"] else ""
        rows.append(
            "<tr>"
            f'<td><a href="system.py?db={esc(db_name)}&id={row["id"]}">{esc(row["name"])}</a></td>'
            f'<td>{esc(row["quadrant"])}</td>'
            f'<td>{"Yes" if row["is_binary"] else "No"}</td>'
            f'<td>{esc(star_type or "")}</td>'
            f'<td>{linkify_location(db_name, row["location"], name_to_id)}</td>'
            "</tr>"
        )
        if row["position_x_mpc"] is not None and row["stars"]:
            # One entry for a single star, two (primary, then secondary)
            # for a binary -- each with its own star_type/temperature/
            # radius/luminosity, so the map can draw and color each
            # component of a binary independently.
            map_systems.append({
                "id": row["id"],
                "name": row["name"],
                "quadrant": row["quadrant"],
                "location": row["location"],
                "x": row["position_x_mpc"],
                "y": row["position_y_mpc"],
                "z": row["position_z_mpc"],
                "stars": [
                    {
                        "star_type": star["star_type"],
                        "temperature_k": star["temperature_k"],
                        "radius_km": star["radius_km"],
                        "luminosity_w": star["luminosity_w"],
                        "temp_display": f'{int(star["temperature_k"])} K',
                    }
                    for star in row["stars"]
                ],
            })
    rows_html = "".join(rows) or '<tr><td colspan="5"><em>None</em></td></tr>'

    def _phenomenon_row_html(row):
        radius_text = f"{row['radius_ly']:,.2f} ly" if row["radius_ly"] else "&ndash;"
        return (
            "<tr>"
            f'<td><a href="phenomenon.py?db={esc(db_name)}&amp;type={esc(row["type"])}&amp;id={row["id"]}">{esc(row["name"])}</a></td>'
            f'<td>{esc(_PHENOMENON_TYPE_LABELS.get(row["type"], row["type"]))}</td>'
            f'<td>{esc((row["descriptor"] or "").replace("_", " ").capitalize())}</td>'
            f'<td>{radius_text}</td>'
            f'<td>{row["distance_ly"]:,.1f} ly</td>'
            "</tr>"
        )

    phenomena_rows_html = "".join(_phenomenon_row_html(row) for row in (sector.get("phenomena") or []))

    center_pc = (
        (sector["center_x_pc"], sector["center_y_pc"], sector["center_z_pc"])
        if sector["placed"]
        else None
    )
    map_html = render_map_panel(
        db_name, sector["edge_mpc"], sector["shell_index"], sector["shell_slot_index"], center_pc, map_systems,
        phenomena=sector.get("phenomena"),
    )

    edge_text = f"{sector['edge_ly']:,.2f} ly"

    system_count = len(systems)
    phenomenon_count = len(sector.get("phenomena") or [])
    badge_bits = [
        f"Cube edge {esc(edge_text)}",
        f"{system_count} system{'s' if system_count != 1 else ''}",
    ]
    if phenomenon_count:
        phenomenon_word = "phenomenon" if phenomenon_count == 1 else "phenomena"
        badge_bits.append(f"{phenomenon_count} nearby exotic {phenomenon_word}")
    if sector["placed"]:
        quadrant = sector_quadrant(sector["center_x_pc"], sector["center_y_pc"])
        badge_bits.append(
            f'<a href="galaxy.py?db={esc(db_name)}&amp;quadrant={quadrant}">View on Galaxy Map (Quadrant {quadrant})</a>'
        )
    badges_html = "<p class=\"badges\">" + "".join(
        f'<span class="badge">{bit}</span>' for bit in badge_bits
    ) + "</p>"

    phenomena_section_html = ""
    if phenomenon_count:
        phenomena_section_html = f"""
<section class="panel">
<h2>Nearby Exotic Phenomena</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Type</th><th>Descriptor</th><th>Radius</th><th>Distance</th></tr></thead>
  <tbody>{phenomena_rows_html}</tbody>
</table></div>
</section>
"""

    wiki_html = ""
    if identity is not None or sector["wiki_url"]:
        wiki_config = get_wiki_config() if identity is not None else {"wikijs": False, "mediawiki": False}
        wiki_html = _wiki_section_html(db_name, sector_id, sector, wiki_config, wiki_message, wiki_error)

    body = f"""
<p class="breadcrumb"><a href="browse.py?db={esc(db_name)}">{esc(db_name)}</a> &rarr; {esc(sector['name'])}</p>
{badges_html}
{wiki_html}
{map_html}
<section class="panel">
<h2>Systems</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Octant</th><th>Binary</th><th>Star type</th><th>Location</th></tr></thead>
  <tbody>{rows_html}</tbody>
</table></div>
</section>
{phenomena_section_html}
<script src="static/sectormap.js" defer></script>
"""
    return f"Sector: {sector['name']}", body


cookie_header = incoming_cookie_header()
try:
    identity = auth_me(cookie_header)
except ApiError:
    # Fails quiet, same as html/system.py's own admin-session check.
    identity = None

wiki_message = None
wiki_error = None
if identity is not None and os.environ.get("REQUEST_METHOD", "GET").upper() == "POST":
    fields = form_params()
    if fields.get("action") == "upload_wiki":
        params = query_params()
        db_name = params.get("db", "")
        sector_id = params.get("id", "")
        backend = fields.get("backend", "")
        path = fields.get("path", "").strip() or None
        try:
            page = upload_sector_to_wiki(cookie_header, db_name, sector_id, backend, path)
            wiki_message = f"Uploaded to the wiki: {page['url']}"
        except ApiError as exc:
            wiki_error = "A page already exists at that location." if exc.status_code == 409 else str(exc)

run(handler)
