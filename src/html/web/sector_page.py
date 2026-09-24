# html/web/sector_page.py

"""
The sector page, `/sector/<id>` (was `sector.py`): the sector's size and
badges, its interactive 3D Sector Map (`lib/starmap.py` data, drawn by
`static/sectormap.js`), and one "Contents" table of its systems and the
phenomena near it, nearest the sector's center first, paged with
`?contents_page=N`.

Admin actions are POST forms to the same URL, each carrying
`csrf_field()` and an `action` field:

- `upload_wiki` (any logged-in admin, while the sector has no wiki page
  yet): `POST /api/sectors/<id>/wiki` with the chosen `backend` and
  optional `path`.
- `generate_neighborhood` (an admin whose credentials are current, on a
  galaxy-placed sector): `POST /api/sectors/<id>/generate-neighborhood`.

A successful action flashes its message and redirects (303) back to the
GET page, so reloading never repeats it. A failed one re-renders the page
with the error next to its form. A POST without an admin session, or with
an unknown `action`, just redirects back to the page.
"""

import re

from flask import flash, get_flashed_messages, redirect, request, url_for

import apiclient
from fmt import format_distance_ly, linkify_location
from galaxymap import sector_quadrant
from pagination import page_slice, parse_page
from starmap import render_map_panel

from . import bp
from .helpers import crumb, current_admin, db_name, page_url, pager, render_page, trusted_html

PHENOMENON_TYPE_LABELS = {
    "nebula": "Nebula", "asteroid_field": "Asteroid Field",
    "black_hole": "Black Hole", "neutron_star": "Neutron Star",
    "supernova_remnant": "Supernova Remnant",
    "rogue_planet": "Rogue Planet", "interstellar_comet": "Interstellar Comet",
}
"""dict: Display labels for `queryDb.phenomena_near_sector`'s `type`
values (the same ones the phenomena list uses)."""

_WIKI_BACKENDS = (("wikijs", "Wiki.js"), ("mediawiki", "MediaWiki"))

_FLASH_CATEGORY = "sector"


def _api_message(exc):
    """An `ApiError`'s message without `apiclient`'s "planetGen API error
    (NNN): " prefix."""
    return re.sub(r"^planetGen API error \(\d+\): ", "", str(exc))


def _system_star_type(row):
    # binary_type only describes a 'close' (P-type) pair's merged star --
    # NULL for a single star and for a 'wide' (S-type) pair, which fall
    # back to joining each star's own type (see schema.sql's "v15" note).
    if row["is_binary"] and row.get("binary_type"):
        return row["binary_type"]
    return " / ".join(star["star_type"] for star in row["stars"]) if row["stars"] else ""


def _map_system(row):
    """One placed system in `starmap.render_map_panel`'s input shape."""
    return {
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
    }


def _contents(sector):
    """
    Every system and phenomenon as Contents rows, nearest the sector's
    center first (anything without a position last), plus the placed
    systems for the map.
    """
    systems = sector["systems"]
    name_to_id = {row["name"]: row["id"] for row in systems}

    def system_url(system_id):
        return page_url("system", system_id=system_id)

    rows = []
    map_systems = []
    for row in systems:
        rows.append({
            "distance_ly": row.get("center_distance_ly"),
            "name": row["name"],
            "url": system_url(row["id"]),
            "type": "Binary Star System" if row["is_binary"] else "Star System",
            "details": _system_star_type(row),
            "octant": row["quadrant"],
            "location": trusted_html(linkify_location(None, row["location"], name_to_id, system_url=system_url)),
        })
        if row["position_x_mpc"] is not None and row["stars"]:
            map_systems.append(_map_system(row))

    for row in (sector.get("phenomena") or []):
        details = [(row["descriptor"] or "").replace("_", " ").capitalize()]
        if row["radius_ly"]:
            details.append(f"{row['radius_ly']:,.2f} ly radius")
        rows.append({
            "distance_ly": row["distance_ly"],
            "name": row["name"],
            "url": page_url("phenomenon", phenomenon_type=row["type"], phenomenon_id=row["id"]),
            "type": PHENOMENON_TYPE_LABELS.get(row["type"], row["type"]),
            "details": ", ".join(bit for bit in details if bit),
            "octant": None,
            "location": None,
        })

    rows.sort(key=lambda entry: (entry["distance_ly"] is None, entry["distance_ly"] or 0.0))
    for entry in rows:
        entry["distance"] = trusted_html(format_distance_ly(entry["distance_ly"]))
    return rows, map_systems


def _handle_post(sector_id, admin):
    """
    Runs one admin action. Returns a redirect response on success (or for
    a POST that is not an admin action at all), else `(form, error)` for
    the page to show next to that form.
    """
    page_again = redirect(url_for("web.sector", sector_id=sector_id), code=303)
    if admin is None:
        return page_again
    action = request.form.get("action", "")
    cookie_header = request.headers.get("Cookie")

    if action == "upload_wiki":
        backend = request.form.get("backend", "")
        path = request.form.get("path", "").strip() or None
        try:
            result = apiclient.upload_sector_to_wiki(cookie_header, db_name(), sector_id, backend, path)
        except apiclient.NotFoundError as exc:
            return "wiki", str(exc)
        except apiclient.ApiError as exc:
            if exc.status_code == 409:
                return "wiki", "A page already exists at that location."
            return "wiki", _api_message(exc)
        flash(f"Uploaded to the wiki: {result['url']}", _FLASH_CATEGORY)
        return page_again

    if action == "generate_neighborhood" and not admin["must_change_credentials"]:
        try:
            result = apiclient.generate_sector_neighborhood(cookie_header, sector_id)
        except apiclient.NotFoundError as exc:
            # e.g. "this sector was never placed in a galaxy": shown on
            # the page rather than as a 404 for an existing sector.
            return "neighborhood", str(exc)
        except apiclient.ApiError as exc:
            return "neighborhood", _api_message(exc)
        flash(
            f'Generated {result["generated"]} new sector(s) ({result["already_existed"]} already '
            f'existed, {result["candidates"]} candidate slot(s) within radius).',
            _FLASH_CATEGORY,
        )
        return page_again

    return page_again


@bp.route("/sector/<int:sector_id>", methods=["GET", "POST"])
def sector(sector_id):
    """One sector: badges, Sector Map, Contents (`?contents_page=N`), and
    the admin forms for a logged-in admin."""
    admin = current_admin()
    errors = {}
    if request.method == "POST":
        outcome = _handle_post(sector_id, admin)
        if not isinstance(outcome, tuple):
            return outcome
        form, message = outcome
        errors[form] = message

    detail = apiclient.get_sector(db_name(), sector_id)
    rows, map_systems = _contents(detail)
    page_rows, contents_page = page_slice(rows, parse_page(request.args.get("contents_page")))

    center_pc = (
        (detail["center_x_pc"], detail["center_y_pc"], detail["center_z_pc"]) if detail["placed"] else None
    )
    map_html = render_map_panel(
        page_url, detail["edge_mpc"], detail["shell_index"], detail["shell_slot_index"], center_pc,
        map_systems, phenomena=detail.get("phenomena"), neighbors=detail.get("neighbors"),
    )

    quadrant = sector_quadrant(detail["center_x_pc"], detail["center_y_pc"]) if detail["placed"] else None
    system_count = len(detail["systems"])
    phenomenon_count = len(detail.get("phenomena") or [])

    wiki_backends = []
    if admin is not None and not detail["wiki_url"]:
        wiki_config = apiclient.get_wiki_config()
        wiki_backends = [(name, label) for name, label in _WIKI_BACKENDS if wiki_config.get(name)]

    return render_page(
        "sector.html",
        title=detail["name"],
        section="sectors",
        breadcrumbs=[crumb("Sectors", "sectors"), crumb(detail["name"])],
        description=f"Sector {detail['name']}: {system_count} star system"
                    f"{'' if system_count == 1 else 's'}, its 3D sector map and nearby phenomena.",
        sector=detail,
        sector_id=sector_id,
        edge_text=f"{detail['edge_ly']:,.2f} ly",
        system_count=system_count,
        phenomenon_count=phenomenon_count,
        quadrant=quadrant,
        map_html=trusted_html(map_html),
        rows=page_rows,
        total_rows=len(rows),
        pager=pager("contents_page", contents_page, len(rows), anchor="sector-contents",
                    label="Contents pages"),
        admin=admin,
        can_generate=admin is not None and not admin["must_change_credentials"],
        wiki_backends=wiki_backends,
        errors=errors,
        messages=get_flashed_messages(category_filter=[_FLASH_CATEGORY]),
    )
