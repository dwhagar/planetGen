# planetgen/web/sector_page.py

"""
The sector page, `/sector/<id>` (was `sector.py`): the sector's size and
badges, its interactive 3D Sector Map (`planetgen/web/maps/starmap.py` data, drawn by
`static/sectorscene.js`), and one "Contents" table of its systems and the
phenomena near it, nearest the sector's center first: the shared data
table (UX.41, `lib/datatable.py`; `contents_*` query parameters). The table also lists the sector's
facilities outside its systems (schema v42): stand-alone ones parked in
open space and those on its asteroid fields.

Admin actions are POST forms to the same URL, each carrying
`csrf_field()` and an `action` field:

- `upload_wiki` (any logged-in admin, while the sector has no wiki page
  yet): `POST /api/sectors/<id>/wiki` with the chosen `backend` and
  optional `path`.
- `set_wiki_url` (any logged-in admin): sets the sector's wiki link, or
  clears it when the URL is blank (UX.73; for a hand-written page from
  outside this app, or to fix a link an upload set).
- `generate_neighborhood` (an admin whose credentials are current, on a
  galaxy-placed sector): the size and time first (`POST
  /api/sectors/<id>/generate-neighborhood` with `estimate_only`), then,
  confirmed, a Generate page job (`planetgen galaxy --center-sector`,
  ADM.11) that keeps running when the page is closed.
- The Delete and Regenerate buttons (ADM.8) post an `edit_action`
  instead, handled by `web/edit_actions.py`: delete the sector with
  everything in it, or delete it and generate its slot again.

A successful action flashes its message and redirects (303) back to the
GET page, so reloading never repeats it. A failed one re-renders the page
with the error next to its form. A POST without an admin session, or with
an unknown `action`, just redirects back to the page.
"""

import math
import re
from urllib.parse import urlencode

from flask import abort, current_app, flash, g, get_flashed_messages, jsonify, redirect, request, url_for

from planetgen.web.lib import apiclient
from planetgen.web.lib.fmt import (
    esc, format_number, format_distance_ly, inside_text, linkify_nearest, nearest_systems_html,
    runaway_text,
)
from planetgen.web.maps.galaxymap import sector_quadrant
from planetgen.web.lib.datatable import Column, Facet, Table, in_memory, parts_of, plain
from planetgen.web.lib.tilecache import fetch_tiles
from planetgen.web.maps.galaxymap3d import initial_tile_request, render_galaxy_map3d_panel, view_radius_bounds
from planetgen.web.maps.starmap import map_scene_data
from planetgen.web.lib.systempage import facility_kind_label

from planetgen.api.common import is_http_url
from planetgen.api.limiter import page_limit
from planetgen.admin import activity_log
from planetgen import tuning
from planetgen.util import log
from planetgen.galaxy.geometry import provisional_sector_designation
from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY
from planetgen.physics.units import ly_to_pc, pc_to_ly

from . import bp, edit_actions, generate_page, jobs, tables
from .helpers import bookmark, crumb, current_admin, db_name, generate_target, page_url, render_page, trusted_html

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
    "About 7.00 × 10¹³ ...", or `None` for an empty sector."""
    if not count or count < 1:
        return None
    return trusted_html(f"About {esc(format_number(round(count)))} interstellar comets and planetesimals (estimated)")


def _rogue_group_row(rogues, system_url):
    """One Contents row holding every rogue planet near the sector (a
    `<details>` row the table's full width, UX.24)."""
    rogues = sorted(rogues, key=lambda r: (r["distance_ly"] is None, r["distance_ly"] or 0.0, r["name"]))
    members = [
        {
            "name": row["name"],
            "url": page_url("phenomenon", phenomenon_type="rogue_planet", phenomenon_id=row["id"]),
            "details": ", ".join(bit for bit in (
                f"Class {row['class']}" if row.get("class") else None,
                (row["descriptor"] or "").capitalize(),
            ) if bit),
            "distance": trusted_html(format_distance_ly(row["distance_ly"])),
            "octant": row.get("octant"),
            "location": _nearest_html(row.get("nearest"), system_url),
            "map_target": f"rogue_planet:{row['id']}",
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
        "group": ROGUE_GROUP,
    }


def _system_nearest_html(row, name_to_id, system_url):
    """A system's Nearest cell: its three nearest neighbours linked (UX.52), or
    `None` when it has none."""
    if row.get("nearest"):
        return trusted_html(nearest_systems_html(row["nearest"], system_url))
    text = linkify_nearest(row["location"], name_to_id, system_url)
    return trusted_html(text) if text else None


def _nearest_html(nearest, system_url):
    """A phenomenon's Nearest cell, or `None` when it has no stored
    neighbors."""
    return trusted_html(nearest_systems_html(nearest, system_url)) if nearest else None


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


SYSTEM_GROUP, PHENOMENON_GROUP, ROGUE_GROUP = 0, 1, 2
"""int: The Contents table's groups in the order they're listed (UX.24):
star systems, then other phenomena and facilities, then rogue planets."""


def _contents(sector, facilities=(), fold=True):
    """
    Every system, phenomenon and facility outside a system as Contents
    rows: star systems first, then the other phenomena and facilities,
    then rogue planets (UX.24), each group nearest the sector's center
    first (anything without a position last), plus the placed systems for
    the map. With `fold` false, every rogue planet is a row of its own, not
    one folded row.
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
            "group": SYSTEM_GROUP,
            "type": "Binary Star System" if row["is_binary"] else "Star System",
            "details": ", ".join(
                bit for bit in (_system_star_type(row), runaway_text(row), inside_text(row)) if bit
            ),
            "octant": row["quadrant"],
            "location": _system_nearest_html(row, name_to_id, system_url),
        })
        if row["position_x_mpc"] is not None and row["stars"]:
            map_systems.append(_map_system(row))

    phenomena = sector.get("phenomena") or []
    rogues = [row for row in phenomena if row["type"] == "rogue_planet"]
    if fold and len(rogues) > 1:
        # Boss's choice: many rogue planets read as one folded row, not a
        # screenful of near-identical ones.
        phenomena = [row for row in phenomena if row["type"] != "rogue_planet"]
        rows.append(_rogue_group_row(rogues, system_url))

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
            # MAP.46: a rogue planet's row can point at it on the Sector Map.
            "map_target": f"rogue_planet:{row['id']}" if row["type"] == "rogue_planet" else None,
            "group": ROGUE_GROUP if row["type"] == "rogue_planet" else PHENOMENON_GROUP,
        })

    rows.extend(dict(row, group=PHENOMENON_GROUP) for row in _facility_rows(sector, facilities))

    rows.sort(key=lambda entry: (entry["group"], entry["distance_ly"] is None, entry["distance_ly"] or 0.0))
    for entry in rows:
        entry["distance"] = trusted_html(format_distance_ly(entry["distance_ly"]))
    return rows, map_systems


def _sector_data(sector_id):
    """The sector and its facilities, fetched once per request."""
    cache = g.setdefault("sector_data", {})
    if sector_id not in cache:
        cache[sector_id] = (apiclient.get_sector(db_name(), sector_id),
                            apiclient.get_sector_facilities(db_name(), sector_id))
    return cache[sector_id]


_ROGUE_LABEL = PHENOMENON_TYPE_LABELS["rogue_planet"]
_NO_VALUE = "–"


def _contents_cells(row, sector_id):
    """A Contents row as the data table's cells (UX.41). The folded rogue-planet row links to
    the table filtered to Rogue Planet, which lists them one by one."""
    href = row["url"]
    if row.get("members"):
        href = f"{page_url('sector', sector_id=sector_id)}?{urlencode({'contents_type': _ROGUE_LABEL})}#sector-contents"
    name = {"text": row["name"]}
    if href:
        name["href"] = href
    if row.get("map_target"):
        name["map_target"] = row["map_target"]
    location = {"text": _NO_VALUE}
    if row["location"] is not None:
        location = {"text": plain(row["location"]), "parts": parts_of(row["location"])}
    return [name, {"text": row["type"]}, {"text": row["details"]}, {"text": row["octant"] or _NO_VALUE},
            location, {"text": plain(row["distance"])}]


def _text_key(field):
    def key(row, descending):
        value = row.get(field)
        return str(value).casefold() if value else None
    return key


def _distance_key(row, descending):
    """The listed order: star systems, other phenomena, rogue planets (UX.24), each by distance."""
    distance = row["distance_ly"]
    if distance is None:
        return None
    return (-row["group"] if descending else row["group"], distance)


_CONTENTS_SORTS = {
    "name": _text_key("name"), "type": _text_key("type"), "details": _text_key("details"),
    "octant": _text_key("octant"),
    "location": lambda row, descending: plain(row["location"]).casefold() if row["location"] is not None else None,
    "distance": _distance_key,
}


def _contents_load(state, limit, offset, want_facets):
    """The sector's Contents (its id is the page's, or the table route's `sector`), filtered,
    sorted and sliced."""
    sector_id = request.view_args.get("sector_id") or request.args.get("sector", type=int)
    try:
        detail, facilities = _sector_data(sector_id)
    except apiclient.NotFoundError:
        abort(404)
    rogue_wanted = _ROGUE_LABEL in state.filters["contents_type"]
    rows, _ = _contents(detail, facilities, fold=not rogue_wanted)
    result = in_memory(
        rows, state, limit, offset, want_facets, _CONTENTS_SORTS,
        {"contents_type": lambda row: row["type"], "contents_octant": lambda row: row["octant"] or None},
        lambda row: _contents_cells(row, sector_id))
    if want_facets and not rogue_wanted:
        rogues = sum(1 for row in detail.get("phenomena") or [] if row["type"] == "rogue_planet")
        for option in result.facets["contents_type"]:
            if option["value"] == _ROGUE_LABEL and rogues > 1:
                option["count"] = rogues
    return result


CONTENTS_TABLE = tables.register(Table(
    "sector-contents", "Contents",
    [Column("name", "Name"), Column("type", "Type"), Column("details", "Details"), Column("octant", "Octant"),
     Column("location", "Nearest"), Column("distance", "From center")],
    _contents_load,
    facets=[Facet("contents_type", "Type"), Facet("contents_octant", "Octant")],
    prefix="contents_", noun=("item", "items"), default_sort="distance",
))


def _handle_post(sector_id, admin):
    """
    Runs one admin action. Returns a redirect response on success (or for
    a POST that is not an admin action at all), else `(form, error)` for
    the page to show next to that form.
    """
    page_again = redirect(url_for("web.sector", sector_id=sector_id), code=303)
    if "edit_action" in request.form:
        return edit_actions.handle_post({("sector", sector_id)}, page_url("sector", sector_id=sector_id),
                                        after_delete={"sector": page_url("sectors")})
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
        if result.get("url"):
            flash(f"Uploaded to the wiki: {result['url']}", _FLASH_CATEGORY)
        else:
            flash(f"The wiki upload is queued (job {result['job_id']}); reload in a minute to see the link.",
                  _FLASH_CATEGORY)
        return page_again

    if action == "set_wiki_url":
        wiki_url = request.form.get("wiki_url", "").strip() or None
        try:
            apiclient.admin_set_sector_wiki_url(cookie_header, db_name(), sector_id, wiki_url)
        except apiclient.NotFoundError as exc:
            return "wiki_link", str(exc) or "Not found."
        except apiclient.ApiError as exc:
            if exc.status_code is None or exc.status_code >= 500:
                raise
            return "wiki_link", _api_message(exc)
        flash("Wiki link cleared." if wiki_url is None else f"Wiki link set to {wiki_url}.", _FLASH_CATEGORY)
        return page_again

    if action == "generate_neighborhood" and not admin["must_change_credentials"]:
        if not request.form.get("estimate_ok"):
            # PERF.3: the size and time first, confirmed on the page.
            try:
                result = apiclient.generate_sector_neighborhood(cookie_header, sector_id, estimate_only=True)
            except apiclient.NotFoundError as exc:
                return "neighborhood", str(exc)
            except apiclient.ApiError as exc:
                return "neighborhood", _api_message(exc)
            return "neighborhood_estimate", result["estimate"]
        # ADM.11: a background job, like the Generate page's, so it runs
        # to the end whether or not this page stays open.
        error = _start_neighborhood_job(sector_id, admin,
                                        generate_anyway=bool(request.form.get(generate_page.GENERATE_ANYWAY_FIELD)))
        if error:
            return "neighborhood", error
        flash("Started generating the sectors around this one. It keeps running if you close this page; "
              "follow it on the Generate page.", _FLASH_CATEGORY)
        return page_again

    return page_again


def _start_neighborhood_job(sector_id, admin, generate_anyway=False):
    """Starts `planetgen galaxy --center-sector <id>` over the default
    radius as a Generate page job (`web/jobs.py`). `generate_anyway`: the
    admin chose "Generate anyway" over a no-room refusal (ADM.33), which
    goes in the activity log. Returns an error message, or `None` once it
    has started."""
    database = db_name()
    form = {"mode": "center", "center_sector": str(sector_id),
            "center_radius_pc": str(tuning.DEFAULT_GENERATE_RADIUS_PC)}
    try:
        kind, title, steps = generate_page.build_job("galaxy", form, database)
        env = jobs.mysql_env(current_app.config["MYSQL_CONFIG"], database)
        job_id = jobs.start_job(kind, title, steps, env=env, admin=admin.get("username"), database=database)
    except generate_page.FormError as exc:
        return str(exc)
    except jobs.JobBusy as exc:
        running = exc.job["title"] if exc.job else "Another job"
        return f"{running} is still running. Wait for it to finish, or cancel it on the Generate page."
    except OSError as exc:
        log.error(f"Could not start the neighborhood job for sector {sector_id}: {exc}")
        return f"The job could not be started: {exc}"
    activity_log.event("GEN", "job.start", user=admin.get("username"), job=job_id, kind=kind, db=database,
                      title=title)
    if generate_anyway:
        activity_log.event("GEN", "job.generate_anyway", user=admin.get("username"), job=job_id, db=database,
                           title=title)
    return None


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


_PICK_LABELS = {"from": ("Choosing a start", "Start Here"),
                "to": ("Choosing a destination", "End Here")}


def _pick_mode(args):
    """
    The NAV pick mode from `?pick=from|to&from=...|to=...` (design doc
    section 9.1), or `None`.

    Returns:
        dict or None: `{"pick", "other", "banner", "label", "cancel",
            "query", "keep_name", "keep_value", "nav_url"}`: which
            endpoint is being chosen, the other endpoint already chosen
            (an `endpoint()` string, or `None`), the banner text, the
            button label, the Cancel URL (`/nav` with the other endpoint
            kept), the `?pick=...` query a link carries so the pick
            continues on the page it opens, and what the page's
            bookmarks need to keep the pick (NAV.40, `static/bookmarks.js`
            reads them as `data-pick`, `data-keep-name`,
            `data-keep-value` and `data-nav-url`). A bad `pick` or other
            endpoint drops pick mode rather than failing the page.
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
    params = {"pick": pick}
    if other:
        params[other_field] = other
    return {"pick": pick, "other": other, "banner": banner, "label": label, "cancel": cancel,
            "query": "?" + urlencode(params, safe=":"), "keep_name": other_field, "keep_value": other or "",
            "nav_url": nav_url()}


def _sector_map_html(detail, pick, admin):
    """
    The sector page's map (MAP.68): the Galaxy Map's own engine
    (`render_galaxy_map3d_panel`) locked to this sector, which it opens in
    place (`static/galaxysector.js`). A sector with no place in the galaxy
    has no map; the page says so.
    """
    from .galaxy_views import phenomenon_url_template, sector_url_template, system_url_template  # galaxy_views imports this module

    address = (detail.get("ring_index"), detail.get("layer_index"), detail.get("ring_slot_index"))
    if not detail["placed"] or None in address:
        return ('<section class="panel" id="map"><h2>Sector Map</h2>'
                '<p class="hint">This sector has no place in the galaxy, so it has no map.</p></section>')
    db = db_name()
    galaxy_shape = apiclient.get_galaxy_shape(db)
    edge_pc = galaxy_shape["edge_pc"] if galaxy_shape else ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    min_radius, max_radius = view_radius_bounds(edge_pc, galaxy_shape)
    center = (detail["center_x_pc"], detail["center_y_pc"], detail["center_z_pc"])
    initial_view = fetch_tiles(db, initial_tile_request(min(max_radius, max(min_radius, 6 * edge_pc)), center))
    return render_galaxy_map3d_panel(
        db, galaxy_shape, edge_pc, initial_view,
        fetch_path=url_for("web.galaxy_tiles"),
        stage_path=url_for("web.galaxy_stage"),
        locate_path=url_for("web.galaxy_locate"),
        territory_path=None,
        sector_url=sector_url_template(),
        generate=generate_target(admin),
        pick=pick,
        phenomenon_url=phenomenon_url_template(),
        system_url=system_url_template(),
        nav_url=page_url("nav"),
        pinned={"ring": address[0], "layer": address[1], "slot": address[2], "center_pc": center},
    )


def sector_designation(sector):
    """A sector's galaxy designation (`provisional_sector_designation`),
    or `None` for a sector with no galaxy address."""
    address = (sector.get("ring_index"), sector.get("layer_index"), sector.get("ring_slot_index"))
    if None in address:
        return None
    try:
        return provisional_sector_designation(*address)
    except ValueError:
        return None


def sector_bookmark(sector, sector_id):
    """
    The sector page's ☆ Bookmark entry (MAP.23): kind `"sector"`, keyed
    by its designation so the Galaxy Map's ☆ on the same sector finds it
    (its page URL when it has no galaxy address), with its id for the
    NAV page's system picker.
    """
    url = page_url("sector", sector_id=sector_id)
    return bookmark("sector", sector_designation(sector) or url, sector["name"], url, sector_id=sector_id)


def galaxy_map_url(sector):
    """
    "Show on Galaxy Map" for a sector (MAP.25): `/galaxy?sector=
    <designation>`, which opens the map's stage 8 holding that sector with
    it selected (design doc section 8.1). `None` for a sector with no
    galaxy address.
    """
    designation = sector_designation(sector)
    return page_url("galaxy", sector=designation) if designation else None


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


@bp.route("/sector/<int:sector_id>/scene")
@page_limit("galaxy_tiles")
def sector_scene(sector_id):
    """
    The Sector Map's scene JSON for one sector (`starmap.map_scene_data`,
    the same block the sector page embeds), for the Galaxy Map, which
    opens the sector as the drill-down's last stage on its own page
    (MAP.66). `?pick=...` keeps NAV's pick mode in the entries' links, as
    on the sector page.
    """
    detail = apiclient.get_sector(db_name(), sector_id)
    _rows, map_systems = _contents(detail, apiclient.get_sector_facilities(db_name(), sector_id))
    center_pc = (
        (detail["center_x_pc"], detail["center_y_pc"], detail["center_z_pc"]) if detail["placed"] else None
    )
    scene = map_scene_data(
        page_url, detail["edge_mpc"],
        (detail.get("ring_index"), detail.get("layer_index"), detail.get("ring_slot_index")),
        center_pc, map_systems, phenomena=detail.get("phenomena"),
        neighbors=_with_bright_stars(detail.get("neighbors")),
        generate=generate_target(current_admin()),
    )
    response = jsonify(scene)
    response.headers["Cache-Control"] = "no-store"
    return response


sector_scene.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/sector/<int:sector_id>", methods=["GET", "POST"])
def sector(sector_id):
    """One sector: badges, Sector Map, Contents (`contents_*` parameters), and
    the admin forms for a logged-in admin."""
    from .nearby_page import nearby_url  # nearby_page imports nav_page, which imports this module

    admin = current_admin()
    errors = {}
    estimate = None
    if request.method == "POST":
        outcome = _handle_post(sector_id, admin)
        if not isinstance(outcome, tuple):
            return outcome
        form, message = outcome
        if form == "neighborhood_estimate":
            estimate = message
        else:
            errors[form] = message

    g.setdefault("sector_data", {}).pop(sector_id, None)  # a POST above may have changed it
    detail, facilities = _sector_data(sector_id)
    pick = _pick_mode(request.args)
    if detail.get("wiki_url") and not is_http_url(detail["wiki_url"]):
        # Saved before the API checked it: never link to a javascript:/
        # data: URL.
        detail["wiki_url"] = None
    _rows, map_systems = _contents(detail, facilities)
    contents = tables.render(CONTENTS_TABLE, request.path, anchor="sector-contents", sector=sector_id)

    map_html = _sector_map_html(detail, pick, admin)

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
        nearby_url=nearby_url(f"sector:{sector_id}"),
        edge_text=format_distance_ly(detail['edge_ly']),
        system_count=system_count,
        star_count=detail.get("star_count"),
        debris=debris_html(detail.get("interstellar_debris_count")),
        phenomenon_count=phenomenon_count,
        quadrant=quadrant,
        galaxy_url=galaxy_map_url(detail),
        bookmark=sector_bookmark(detail, sector_id),
        map_html=trusted_html(map_html),
        contents=contents,
        admin=admin,
        can_generate=admin is not None and not admin["must_change_credentials"],
        can_edit=edit_actions.can_edit(admin),
        wiki_backends=wiki_backends,
        errors=errors,
        estimate=estimate,
        messages=get_flashed_messages(category_filter=[_FLASH_CATEGORY]),
        pick=pick,
    )
