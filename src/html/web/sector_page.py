# html/web/sector_page.py

"""
The sector page, `/sector/<id>` (was `sector.py`): the sector's size and
badges, its interactive 3D Sector Map (`lib/starmap.py` data, drawn by
`static/sectormap.js`), and one "Contents" table of its systems and the
phenomena near it, nearest the sector's center first, paged with
`?contents_page=N`. The Contents table also lists the sector's
facilities outside its systems (schema v42): stand-alone ones parked in
open space and those on its asteroid fields.

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

import math
import re

from flask import flash, get_flashed_messages, redirect, request, url_for

import apiclient
from fmt import (
    esc, format_distance_ly, inside_text, linkify_location, nearest_neighbors_location, nearest_systems_html,
    runaway_text,
)
from galaxymap import sector_quadrant
from pagination import page_slice, parse_page
from starmap import render_map_panel
from systempage import facility_kind_label

from api.common import is_http_url
from stellarObjects.galaxyGeometry import provisional_sector_designation
from stellarObjects.utils import pc_to_ly

from . import bp
from .helpers import crumb, current_admin, db_name, generate_target, page_url, pager, render_page, trusted_html

PHENOMENON_TYPE_LABELS = {
    "nebula": "Nebula", "asteroid_field": "Asteroid Field",
    "black_hole": "Black Hole", "neutron_star": "Neutron Star",
    "supernova_remnant": "Supernova Remnant",
    "rogue_planet": "Rogue Planet", "interstellar_comet": "Interstellar Comet",
    "quasar": "Quasar",
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


def debris_html(count):
    """The sector's estimated interstellar comets and planetesimals
    (`queryDb.interstellar_debris_count`, a computed figure, not rows) as
    "About 7&times;10<sup>13</sup> ...", or `None` for an empty sector."""
    if not count or count < 1:
        return None
    exponent = int(math.floor(math.log10(count)))
    mantissa = round(count / 10 ** exponent)
    if mantissa == 10:
        mantissa, exponent = 1, exponent + 1
    figure = f"{round(count):,}" if exponent < 4 else f"{mantissa}&times;10<sup>{exponent}</sup>"
    return trusted_html(f"About {figure} interstellar comets and planetesimals (estimated)")


def _rogue_group_row(rogues):
    """One Contents row holding every rogue planet near the sector (a
    `<details>` list in the Name cell), placed by the nearest one."""
    rogues = sorted(rogues, key=lambda r: (r["distance_ly"] is None, r["distance_ly"] or 0.0, r["name"]))
    members = [
        {
            "name": row["name"],
            "url": page_url("phenomenon", phenomenon_type="rogue_planet", phenomenon_id=row["id"]),
            "details": (row["descriptor"] or "").capitalize(),
            "distance": trusted_html(format_distance_ly(row["distance_ly"])),
        }
        for row in rogues
    ]
    return {
        "distance_ly": rogues[0]["distance_ly"],
        "name": f"{len(rogues)} rogue planets",
        "url": None,
        "members": members,
        "type": PHENOMENON_TYPE_LABELS["rogue_planet"],
        "details": "Unbound planets drifting between the stars",
        "octant": None,
        "location": None,
    }


def _nearest_html(nearest, system_url):
    """A phenomenon's "Nearest: ..." Location cell, or `None` when it has
    no stored neighbors."""
    return trusted_html("Nearest: " + nearest_systems_html(nearest, system_url)) if nearest else None


def _facility_rows(sector, facilities):
    """
    The sector's own facilities (`GET /api/sectors/<id>/facilities`) as
    Contents rows: a stand-alone one at its distance from the sector's
    center, one on an asteroid field at the field's, linked from the
    Location cell. A facility has no page of its own, so its name is plain
    text.
    """
    fields = {row["id"]: row for row in sector.get("phenomena") or [] if row["type"] == "asteroid_field"}
    center = (sector.get("center_x_pc"), sector.get("center_y_pc"), sector.get("center_z_pc"))
    rows = []
    for facility in facilities:
        kind = facility_kind_label(facility["kind"])
        distance_ly, location = None, None
        if facility["host_type"] == "space":
            details = f"Stand-alone {kind.lower()}, parked in open space"
            point = (facility.get("center_x_pc"), facility.get("center_y_pc"), facility.get("center_z_pc"))
            if None not in center and None not in point:
                distance_ly = pc_to_ly(math.dist(center, point))
        else:
            details = f"{kind} on an asteroid field"
            field = fields.get(facility["host_id"])
            if field is not None:
                distance_ly = field["distance_ly"]
            url = page_url("phenomenon", phenomenon_type="asteroid_field", phenomenon_id=facility["host_id"])
            location = trusted_html(f'On <a href="{esc(url)}">{esc(facility["host_name"] or "an asteroid field")}</a>')
        if facility.get("description"):
            details += f": {facility['description']}"
        rows.append({
            "distance_ly": distance_ly,
            "name": facility["name"],
            "url": None,
            "type": f"Facility ({kind})",
            "details": details,
            "octant": None,
            "location": location,
        })
    return rows


def _contents(sector, facilities=()):
    """
    Every system, phenomenon and facility outside a system as Contents
    rows, nearest the sector's center first (anything without a position
    last), plus the placed systems for the map.
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
            "details": ", ".join(
                bit for bit in (_system_star_type(row), runaway_text(row), inside_text(row)) if bit
            ),
            "octant": row["quadrant"],
            "location": trusted_html(
                nearest_neighbors_location(row["location"], row["nearest"], system_url) if row.get("nearest")
                else linkify_location(row["location"], name_to_id, system_url)
            ),
        })
        if row["position_x_mpc"] is not None and row["stars"]:
            map_systems.append(_map_system(row))

    phenomena = sector.get("phenomena") or []
    rogues = [row for row in phenomena if row["type"] == "rogue_planet"]
    if len(rogues) > 1:
        # Boss's choice: many rogue planets read as one folded row, not a
        # screenful of near-identical ones.
        phenomena = [row for row in phenomena if row["type"] != "rogue_planet"]
        rows.append(_rogue_group_row(rogues))

    for row in phenomena:
        details = [(row["descriptor"] or "").replace("_", " ").capitalize()]
        if row.get("class"):
            details.insert(0, f"Class {row['class']}")
        if row["radius_ly"]:
            details.append(f"{format_distance_ly(row['radius_ly'])} radius")
        rows.append({
            "distance_ly": row["distance_ly"],
            "name": row["name"],
            "url": page_url("phenomenon", phenomenon_type=row["type"], phenomenon_id=row["id"]),
            "type": PHENOMENON_TYPE_LABELS.get(row["type"], row["type"]),
            "details": ", ".join(bit for bit in details if bit),
            "octant": row.get("octant"),
            "location": _nearest_html(row.get("nearest"), system_url),
        })

    rows.extend(_facility_rows(sector, facilities))

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


def _with_bright_stars(neighbors):
    """
    The neighbors with each not-yet-generated one's waiting bright stars
    (`apiclient.get_bright_stars_in_cell`, brightest first) added as
    `bright_stars`, so its Sector Map panel can list them. Fails open: a
    cell that can't be read just lists none.
    """
    result = []
    for neighbor in neighbors or ():
        if not neighbor.get("exists"):
            neighbor = dict(neighbor)
            try:
                neighbor["bright_stars"] = apiclient.get_bright_stars_in_cell(
                    db_name(), neighbor["ring_index"], neighbor["layer_index"], neighbor["ring_slot_index"])
            except apiclient.ApiError:
                neighbor["bright_stars"] = []
        result.append(neighbor)
    return result


_PICK_LABELS = {"from": ("Choosing a start", "Use as start"),
                "to": ("Choosing a destination", "Use as destination")}


def _pick_mode(args):
    """
    The NAV pick mode from `?pick=from|to&from=...|to=...` (design doc
    section 9.1), or `None`.

    Returns:
        dict or None: `{"pick", "other", "banner", "label", "cancel"}`:
            which endpoint is being chosen, the other endpoint already
            chosen (an `endpoint()` string, or `None`), the banner text,
            the button label, and the Cancel URL (`/nav` with the other
            endpoint kept). A bad `pick` or other endpoint drops pick
            mode rather than failing the page.
    """
    from .nav_page import endpoint, nav_url, parse_endpoint  # nav_page imports this module

    pick = args.get("pick")
    if pick not in _PICK_LABELS:
        return None
    other_field = "to" if pick == "from" else "from"
    other = args.get(other_field) or None
    if other:
        try:
            other = endpoint(*parse_endpoint(other))
        except apiclient.NotFoundError:
            return None
    banner, label = _PICK_LABELS[pick]
    cancel = nav_url(destination=other) if pick == "from" else nav_url(origin=other)
    return {"pick": pick, "other": other, "banner": banner, "label": label, "cancel": cancel}


def _nav_for(pick):
    """The Sector Map's `nav(kind, id)` hook (`starmap.render_map_panel`)."""
    from .nav_page import endpoint, nav_url  # nav_page imports this module

    def nav(kind, entity_id):
        here = endpoint(kind, entity_id)
        links = {"from": nav_url(origin=here), "to": nav_url(destination=here), "pick": None, "pickLabel": None}
        if pick is not None:
            if pick["pick"] == "from":
                links["pick"] = nav_url(origin=here, destination=pick["other"])
            else:
                links["pick"] = nav_url(origin=pick["other"], destination=here)
            links["pickLabel"] = pick["label"]
        return links
    return nav


def galaxy_map_url(sector):
    """
    "Show on Galaxy Map" for a sector (MAP.2.7): `/galaxy?sector=
    <designation>`, which opens the map's stage 8 holding that sector with
    it selected (design doc section 8.1). `None` for a sector with no
    galaxy address.
    """
    address = (sector.get("ring_index"), sector.get("layer_index"), sector.get("ring_slot_index"))
    if None in address:
        return None
    try:
        designation = provisional_sector_designation(*address)
    except ValueError:
        return None
    return page_url("galaxy", sector=designation)


@bp.route("/sector/<int:sector_id>/galaxy")
def sector_on_galaxy_map(sector_id):
    """Redirects to the sector on the Galaxy Map (for links that know
    only the sector's id, like search results), or to the plain map when
    the sector has no galaxy address."""
    url = galaxy_map_url(apiclient.get_sector(db_name(), sector_id))
    return redirect(url or page_url("galaxy"), code=302)


@bp.route("/system/<int:system_id>/galaxy")
def system_on_galaxy_map(system_id):
    """Redirects to a system's sector on the Galaxy Map, or to the plain
    map for a system outside any sector."""
    system = apiclient.get_system(db_name(), system_id)
    if system["sector_id"] is None:
        return redirect(page_url("galaxy"), code=302)
    return sector_on_galaxy_map(system["sector_id"])


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
    pick = _pick_mode(request.args)
    if detail.get("wiki_url") and not is_http_url(detail["wiki_url"]):
        # Saved before the API checked it: never link to a javascript:/
        # data: URL.
        detail["wiki_url"] = None
    rows, map_systems = _contents(detail, apiclient.get_sector_facilities(db_name(), sector_id))
    page_rows, contents_page = page_slice(rows, parse_page(request.args.get("contents_page")))

    center_pc = (
        (detail["center_x_pc"], detail["center_y_pc"], detail["center_z_pc"]) if detail["placed"] else None
    )
    neighbors = _with_bright_stars(detail.get("neighbors"))
    map_html = render_map_panel(
        page_url, detail["edge_mpc"],
        (detail.get("ring_index"), detail.get("layer_index"), detail.get("ring_slot_index")),
        center_pc, map_systems, phenomena=detail.get("phenomena"), neighbors=neighbors,
        generate=generate_target(admin),
        nav=_nav_for(pick),
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
        edge_text=format_distance_ly(detail['edge_ly']),
        system_count=system_count,
        star_count=detail.get("star_count"),
        debris=debris_html(detail.get("interstellar_debris_count")),
        phenomenon_count=phenomenon_count,
        quadrant=quadrant,
        galaxy_url=galaxy_map_url(detail),
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
        pick=pick,
    )
