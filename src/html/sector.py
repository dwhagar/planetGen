#!/usr/bin/env python3
# html/sector.py

"""
Sector detail page: the sector's name/size and every system placed in it,
with quadrant and star-type info, linking to `system.py` for each -- plus
an interactive 3D "Sector Map" (see `lib/starmap.py`) of the same systems
plotted by position within the sector, plus a translucent cloud for every
nebula/asteroid field/supernova remnant (and a point marker for every
black hole/neutron star/rogue planet/interstellar comet) whose real
galaxy-frame sphere reaches into this sector's own cube, or that was
generated as part of it (`queryDb.phenomena_near_sector`, via
`GET /api/sectors/<id>`'s `phenomena` key -- see `schema.sql`'s
"v18"/"v21"/"v28" header notes). The "Contents" table lists those
phenomena alongside the sector's systems, nearest the sector's center
first, each row linking to `system.py` or `phenomenon.py`.

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

from apiclient import (
    ApiError,
    NotFoundError,
    auth_me,
    generate_sector_neighborhood,
    get_sector,
    get_wiki_config,
    upload_sector_to_wiki,
)
from fmt import esc, format_distance_ly, linkify_location, post_link
from galaxymap import sector_quadrant
from page import form_params, incoming_cookie_header, nav_params, run
from pagination import page_slice, parse_page, render_pagination
from starmap import render_map_panel

_PHENOMENON_TYPE_LABELS = {
    "nebula": "Nebula", "asteroid_field": "Asteroid Field",
    "black_hole": "Black Hole", "neutron_star": "Neutron Star",
    "supernova_remnant": "Supernova Remnant",
    "rogue_planet": "Rogue Planet", "interstellar_comet": "Interstellar Comet",
    "quasar": "Quasar",
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
    return f"""
{message_html}{error_html}
<section class="panel">
<h2>Upload to Wiki</h2>
<form method="post" action="sector.py" class="search-form">
  <input type="hidden" name="action" value="upload_wiki">
  <input type="hidden" name="db" value="{esc(db_name)}">
  <input type="hidden" name="id" value="{esc(sector_id)}">
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
    params = nav_params()
    db_name = params.get("db", "")
    sector_id = params.get("id", "")

    # Admin-only "generate more sectors around this one" action -- a
    # failed identity check just falls back to the anonymous view rather
    # than failing the whole (otherwise public, read-only) page; see
    # admin.py for the same auth_me/must_change_credentials gate used
    # elsewhere.
    cookie_header = incoming_cookie_header()
    try:
        identity = auth_me(cookie_header)
    except ApiError:
        identity = None
    is_admin = identity is not None and not identity["must_change_credentials"]

    generation_result = None
    generation_error = None
    if is_admin and os.environ.get("REQUEST_METHOD", "GET").upper() == "POST":
        fields = form_params()
        if fields.get("action") == "generate_neighborhood":
            try:
                generation_result = generate_sector_neighborhood(cookie_header, int(sector_id))
            except (ApiError, NotFoundError) as exc:
                # NotFoundError specifically -- not just ApiError -- since
                # apiclient._auth_request raises that instead for a 404
                # (e.g. "this sector was never placed in a galaxy"); left
                # uncaught, it would propagate past this page entirely
                # and render as a full page.run()-level error page
                # instead of this inline message on the sector's own
                # (otherwise perfectly loadable) page.
                generation_error = str(exc)
            except ValueError:
                generation_error = "Invalid sector id."

    sector = get_sector(db_name, sector_id)
    systems = sector["systems"]
    name_to_id = {row["name"]: row["id"] for row in systems}

    # (distance from the sector's center in ly, or None if unplaced; row HTML)
    # for every system and phenomenon, merged into one Contents table below.
    content_rows = []
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
        distance_ly = row.get("center_distance_ly")
        content_rows.append((distance_ly, (
            "<tr>"
            f'<td>{post_link("system.py", {"db": db_name, "id": row["id"]}, esc(row["name"]))}</td>'
            f'<td>{"Binary Star System" if row["is_binary"] else "Star System"}</td>'
            f'<td>{esc(star_type or "")}</td>'
            f'<td>{esc(row["quadrant"])}</td>'
            f'<td>{linkify_location(db_name, row["location"], name_to_id)}</td>'
            f'<td>{format_distance_ly(distance_ly)}</td>'
            "</tr>"
        )))
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

    for row in (sector.get("phenomena") or []):
        details = [(row["descriptor"] or "").replace("_", " ").capitalize()]
        if row["radius_ly"]:
            details.append(f"{row['radius_ly']:,.2f} ly radius")
        details = ", ".join(bit for bit in details if bit)
        phenomenon_link = post_link(
            "phenomenon.py", {"db": db_name, "type": row["type"], "id": row["id"]}, esc(row["name"])
        )
        content_rows.append((row["distance_ly"], (
            "<tr>"
            f'<td>{phenomenon_link}</td>'
            f'<td>{esc(_PHENOMENON_TYPE_LABELS.get(row["type"], row["type"]))}</td>'
            f'<td>{esc(details)}</td>'
            "<td>&ndash;</td>"
            "<td>&ndash;</td>"
            f'<td>{format_distance_ly(row["distance_ly"])}</td>'
            "</tr>"
        )))

    # Nearest the sector's center first; anything with no position last.
    # The map above plots everything, so the full list is already here;
    # the table shows one page of it (see lib/pagination.py).
    content_rows.sort(key=lambda entry: (entry[0] is None, entry[0] or 0.0))
    page_rows, contents_page = page_slice(content_rows, parse_page(params.get("contents_page")))
    contents_html = "".join(html for _distance, html in page_rows) or (
        '<tr><td colspan="6"><em>None</em></td></tr>'
    )
    page_state = {"db": db_name, "id": sector_id, "contents_page": contents_page}

    center_pc = (
        (sector["center_x_pc"], sector["center_y_pc"], sector["center_z_pc"])
        if sector["placed"]
        else None
    )
    map_html = render_map_panel(
        db_name, sector["edge_mpc"], sector["shell_index"], sector["shell_slot_index"], center_pc, map_systems,
        phenomena=sector.get("phenomena"), neighbors=sector.get("neighbors"),
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
        badge_bits.append(f"{phenomenon_count} {phenomenon_word}")
    if sector["placed"]:
        quadrant = sector_quadrant(sector["center_x_pc"], sector["center_y_pc"])
        badge_bits.append(
            post_link("galaxy.py", {"db": db_name, "quadrant": quadrant}, f"View on Galaxy Map (Quadrant {quadrant})")
        )
    badges_html = "<p class=\"badges\">" + "".join(
        f'<span class="badge">{bit}</span>' for bit in badge_bits
    ) + "</p>"

    wiki_html = ""
    if identity is not None or sector["wiki_url"]:
        wiki_config = get_wiki_config() if identity is not None else {"wikijs": False, "mediawiki": False}
        wiki_html = _wiki_section_html(db_name, sector_id, sector, wiki_config, wiki_message, wiki_error)

    admin_panel_html = ""
    if is_admin:
        message_html = ""
        if generation_error:
            message_html = f'<p class="error">{esc(generation_error)}</p>'
        elif generation_result is not None:
            message_html = (
                f'<p class="hint">Generated {generation_result["generated"]} new sector(s) '
                f'({generation_result["already_existed"]} already existed, '
                f'{generation_result["candidates"]} candidate slot(s) within radius).</p>'
            )

        if sector["placed"]:
            action_html = f"""
<p class="hint">Fills in every not-yet-generated sector within a 100 ly sphere
around this one (already-generated sectors are skipped). That sphere can hold
thousands of candidate sectors, so this can take anywhere from a few minutes
to a few hours to finish -- the page will not respond until it completes.</p>
<form method="post" action="sector.py" class="table-form">
  <input type="hidden" name="action" value="generate_neighborhood">
  <input type="hidden" name="db" value="{esc(db_name)}">
  <input type="hidden" name="id" value="{esc(sector_id)}">
  <button type="submit" class="btn">Generate more sectors around this one</button>
</form>
"""
        else:
            action_html = (
                '<p class="hint">This sector has never been placed in a galaxy -- generating a '
                "neighborhood around it requires a galaxy-placed sector (one generated via "
                "'generate.py galaxy', not 'generate.py sector').</p>"
            )

        admin_panel_html = f"""
<section class="panel">
<h2>Admin</h2>
{message_html}
{action_html}
</section>
"""

    breadcrumb_html = f'<p class="breadcrumb">{post_link("browse.py", {"db": db_name}, esc(db_name))} &rarr; {esc(sector["name"])}</p>'
    body = f"""
<div class="page-subhead">{breadcrumb_html}{badges_html}</div>
{wiki_html}
{admin_panel_html}
{map_html}
<section class="panel" id="sector-contents">
<h2>Contents</h2>
<div class="table-scroll"><table>
  <thead><tr><th>Name</th><th>Type</th><th>Details</th><th>Octant</th><th>Location</th><th>From center</th></tr></thead>
  <tbody>{contents_html}</tbody>
</table></div>
{render_pagination("sector.py", page_state, "contents_page", contents_page, len(content_rows),
                   anchor="sector-contents", label="Contents pages")}
</section>
<script type="module" src="static/sectormap.js"></script>
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
        db_name = fields.get("db", "")
        sector_id = fields.get("id", "")
        backend = fields.get("backend", "")
        path = fields.get("path", "").strip() or None
        try:
            page = upload_sector_to_wiki(cookie_header, db_name, sector_id, backend, path)
            wiki_message = f"Uploaded to the wiki: {page['url']}"
        except ApiError as exc:
            wiki_error = "A page already exists at that location." if exc.status_code == 409 else str(exc)

run(handler)
