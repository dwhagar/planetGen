# planetgen/web/item_pages.py

"""
The dedicated page of an item (UX.143, UX.170): `/object/<ref>` for any object reference
(`planet:FE81000A2B-0000005-001`, a bare system ID). Every kind a reference can name has a page, from one registry:

- a sector, a system and each phenomenon type already have a page of their own (`/sector`, `/system`, `/phenomenon`),
  so their reference goes there;
- a star, planet, moon, belt or comet gets the page here, top to bottom: its name, kind, Object ID, class and PHI-6
  symbols with the bookmark, Navigate, Nearby and (for an admin) the Admin menu; where it is (the Location chain, every
  name a link: galaxy, sector, system, star, planet, and the previous and next body at its own level); the property
  sheet (`describe/sheet.py`, every stored property with its unit); what the system's description says about it (the
  prose, facts and ratings, `describe/`); and the summary tree of everything that belongs to it
  (`web/lib/objecttree.py`, `GET /api/objects/<ref>/children`).

The prose is the description's section for the body (`GET /api/systems/<id>/description`), so the page and the system's
description never say different things. A reference that names nothing is a plain not-found page. A POST with an
`edit_action` is the admin Regenerate or Delete of this body (`web/edit_actions.py`).
"""

from flask import abort, jsonify, redirect, request

from planetgen.describe import sheet as sheet_registry
from planetgen.describe.html import blocks_html, inline_html
from planetgen.describe.model import sections
from planetgen.galaxy import objectref
from planetgen.web.lib import apiclient
from planetgen.web.lib import bodyview, facility_display, life_air_chart, planet_layout, planet_summary
from planetgen.web.lib.habitability_chips import chips_html
from planetgen.web.lib.fmt import esc
from planetgen.web.lib.systempage import facilities_by_host, facilities_html
from planetgen.web.lib.objecttree import CHILD_NOUNS, TreeContext, has_poi_mark, nodes_html, tree_html

from . import edit_actions, static_url
from .class_pages import class_url
from .helpers import bookmark, crumb, current_admin, db_name, page_url, population_status, render_page, trusted_html
from .nav_page import endpoint, nav_url
from .nearby_page import nearby_url

OWN_PAGE_KINDS = ("sector", "system") + objectref.PHENOMENON_KINDS
"""tuple: The kinds whose reference opens a page the site already has."""

ITEM_PAGE_KINDS = objectref.BODY_KINDS + objectref.FACILITY_KINDS + objectref.CONTAINER_KINDS
"""tuple: The kinds this module renders."""

KIND_LABELS = {"star": "Star", "planet": "Planet", "moon": "Moon", "belt": "Asteroid belt", "comet": "Comet",
               "facility": "Facility", "globular_cluster": "Globular cluster", "superbubble": "Superbubble"}


def page_kinds():
    """Every kind a reference can name that has a page: the registry `tests/test_item_pages.py` holds equal to
    `objectref.KINDS`."""
    return OWN_PAGE_KINDS + ITEM_PAGE_KINDS


def own_page_url(kind, printed_id):
    """The page of a sector, system or phenomenon, for its printed Object ID."""
    if kind == "sector":
        return page_url("sector", sector_id=printed_id)
    if kind == "system":
        return page_url("system", system_id=printed_id)
    return page_url("phenomenon", phenomenon_type=kind, phenomenon_id=printed_id)


def item_url(ref):
    """Where a reference (`kind:ID`) opens: its own page, or the item page."""
    kind, printed_id = objectref.parse_public(ref)
    return own_page_url(kind, printed_id) if kind in OWN_PAGE_KINDS else page_url("object_page", ref=f"{kind}:{printed_id}")


def _crumbs(found, name):
    """Home > Galaxy > sector > System > star > ... > this item, from the API's parent chain. The system link reads
    "System" and the star's the star's own name, so a system and a star that share a name do not show twice."""
    crumbs = []
    for parent in found["parents"]:
        if parent["kind"] == "galaxy":
            crumbs.append(crumb("Galaxy", "galaxy"))
        else:
            label = "System" if parent["kind"] == "system" else parent["name"]
            crumbs.append({"label": label, "url": item_url(parent["ref"])})
    crumbs.append(crumb(name))
    return crumbs


def location(parents, inside):
    """
    The Location block (UX.144): the places above an item as links, and the cloud it sits inside. The galaxy is left out
    (UX.200: it is the galaxy anyway), and a star that shares the name of the system above it is listed once, as the system.

    Args:
        parents (list[dict]): The API's parent chain (`ref`, `kind`, `name`), galaxy first.
        inside (dict or None): The detail's `inside` (`type`, `id`, `name`), or `None` in open space.
    """
    chain = []
    for parent in parents:
        if parent["kind"] == "galaxy":
            continue
        if parent["kind"] == "star" and chain and chain[-1]["label"] == parent["name"]:
            continue
        chain.append({"label": parent["name"], "url": item_url(parent["ref"])})
    cloud = None if not inside else {"name": inside["name"], "url": own_page_url(inside["type"], inside["id"])}
    return {"chain": chain, "inside": cloud}


def sector_location(sector_id, sector_name, inside):
    """`location` for something whose only place above it is its sector (a system or a phenomenon): `sector_id` is
    `None` for one with no sector."""
    parents = [{"ref": "galaxy", "kind": "galaxy"}]
    if sector_id is not None:
        parents.append({"ref": f"sector:{sector_id}", "kind": "sector", "name": sector_name or "Sector"})
    return location(parents, inside)


CONTENT_LABELS = {
    "system": "Star systems", "nebula": "Nebulae", "supernova_remnant": "Supernova remnants",
    "asteroid_field": "Asteroid fields", "rogue_planet": "Rogue planets", "interstellar_comet": "Interstellar comets",
    "black_hole": "Black holes", "neutron_star": "Neutron stars",
}


def contents(listing):
    """The "What is inside it" block from `GET /api/phenomena/<type>/<id>/contents`, grouped by kind."""
    groups = []
    for kind, label in CONTENT_LABELS.items():
        count = listing["counts"].get(kind)
        if count:
            items = [{"name": item["name"], "url": own_page_url(kind, item["id"])}
                     for item in listing["items"] if item["kind"] == kind]
            groups.append({"label": label, "count": count, "entries": items})
    return {"total": listing["total"], "shown": len(listing["items"]), "groups": groups}


def _section(document, ref):
    """The description's section for `ref`, or `None`."""
    wanted = ref.lower()
    return next((piece for piece in sections(document) if piece["ref"] and piece["ref"].lower() == wanted), None)


def _ref_link(ref):
    """The link a description reference (`planet:ID`) gets: its item page."""
    return item_url(ref) if ref else None


BODY_REF_KINDS = ITEM_PAGE_KINDS
"""tuple: The kinds whose page this module draws (a reference to one is `<kind>:<printed ID>`)."""

TREE_PAGE = 50
"""int: Rows of a tree level the page and each "Show more" read at once."""


def tree_context():
    """The `TreeContext` of a page in this request: where a row links, where a row's children come from, and the class
    pages and habitability page its summary line links into."""
    return TreeContext(
        item_url=item_url, tree_url=lambda ref: page_url("object_tree", ref=ref), class_url=class_url,
        sprite_url=static_url("icons.svg"), help_url=page_url("habitability_levels"), static_url=static_url)


def tree_response(ref):
    """`GET /object/<ref>/tree[?after=]`: one page of an object's children as the HTML rows of a tree level
    (`{"html", "next"}`), for `static/objecttree.js`. A reference that names nothing is a 404, a bad `after` a 400."""
    try:
        kind, printed_id = objectref.parse_public(ref)
    except ValueError:
        abort(404)
    try:
        page = apiclient.get_object_children(db_name(), f"{kind}:{printed_id}", limit=TREE_PAGE,
                                             after=request.args.get("after") or None)
    except apiclient.NotFoundError:
        abort(404)
    except apiclient.ApiError as exc:
        if exc.status_code == 400:
            abort(400)
        raise
    response = jsonify({"html": nodes_html(page["items"], tree_context()), "next": page["next"]})
    response.headers["Cache-Control"] = "no-store"
    return response


def tree_panel(ref, name, page, loaded=None):
    """The "Belonging to it" tree of one object from its first children page, or `None` when it holds nothing. `loaded`
    is `objecttree.node_html`'s: the levels below already read, drawn open."""
    if not page["items"]:
        return None
    noun = CHILD_NOUNS.get(page["kind"], "children")
    return trusted_html(tree_html(page["items"], tree_context(), f"{noun.capitalize()} of {name}",
                                  token=page["next"], ref=ref, loaded=loaded))


def tree_has_marks(page, loaded=None):
    """Whether a tree drawn from `page` (and the levels `loaded` below it) carries a point-of-interest symbol, so the page shows
    the key to them (UX.233)."""
    levels = [page["items"], *(items for items, _ in (loaded or {}).values())]
    return any(has_poi_mark(items) for items in levels)


def system_tree(db, system_id, name):
    """The system page's tree: the system's stars (and whatever orbits the whole system, for a close pair), with each
    star's own bodies already open below it. `(panel, marks)`: the panel (`None` for a system with no stars or bodies) and whether
    its rows carry a point-of-interest symbol (UX.233)."""
    ref = f"system:{system_id}"
    page = apiclient.get_object_children(db, ref, limit=TREE_PAGE)
    loaded = {}
    for item in page["items"]:
        if item["kind"] == "star" and item["child_count"]:
            below = apiclient.get_object_children(db, item["ref"], limit=TREE_PAGE)
            loaded[item["ref"]] = (below["items"], below["next"])
    return tree_panel(ref, name, page, loaded), tree_has_marks(page, loaded)


def _listed(row):
    """Whether a page lists a sheet row: it holds something, or its absence is itself worth knowing (UX.188). A value that
    was never stored is not listed (UX.199)."""
    if not row.get("stored", True):    # UX.199: "not stored" is a gap in the data, not a fact about the object
        return False
    return not row.get("absent") or row.get("show_absent", False)


def _value_html(row):
    """A sheet row's value as HTML. A body's class letter is a link to what that class is (UX.218); a class the catalog
    does not know stays plain text."""
    html = trusted_html(inline_html(row["value"]))
    if row.get("key") == "planet_class" and row.get("stored", True):
        code = row["value"][0] if row["value"] and isinstance(row["value"][0], str) else row["text"]
        url = class_url("planet", code) or class_url("planet", code.lower())
        if url:
            return trusted_html(f'<a href="{esc(url)}">{esc(code)}</a>')
    return html


def sheet_view(sheet):
    """A property sheet (`GET /api/objects/<ref>/sheet`) as the page draws it: `[{"key", "label", "rows"}]` with each
    row's value as HTML (`describe.html.inline_html`, so a quantity reads in the reader's unit system). A row that holds
    nothing is left out unless its absence matters (`describe.sheet.SHOW_WHEN_ABSENT`); a group left empty goes too."""
    absent = set(sheet.get("missing") or ())
    groups = [{"key": group["key"], "label": group["label"],
               "rows": [{"key": row["key"], "label": row["label"], "html": _value_html(row), "stored": row["stored"],
                         "source": row.get("source")}
                        for row in group["rows"] if _listed(row)]}
              for group in sheet["groups"] if group["key"] not in absent]    # UX.203: a section of zeros is one sentence
    if "atmosphere" in absent:
        groups.insert(next((i for i, group in enumerate(groups) if group["key"] == "climate"), len(groups)),
                      {"key": "atmosphere", "label": "Atmosphere", "rows": [
                          {"key": "atmosphere_none", "label": "Atmosphere", "html": trusted_html("None"), "stored": True}]})
    return [group for group in groups if group["rows"]]


def bands_view(sheet):
    """
    The table of the seven bands of a planet or moon (GEN.1108) as the page draws it: `{"caption", "columns": [label], "rows": [{"label", "cells": [html]}]}`, the bands as columns
    (north pole to south pole, or day centre to night centre for a world that turns once an orbit) and the figures of each as rows, a cell with nothing in it a dash.
    `None` for a body with no bands.
    """
    table = sheet.get("bands")
    if not table:
        return None
    caption = ("The climate of the world in rings from the point under its star to the point away from it" if table["locked"]
               else "The climate of the world in seven bands of latitude, from the north pole to the south pole")
    return {"caption": caption, "columns": [band["label"] for band in table["bands"]],
            "rows": [{"label": row["label"], "cells": [trusted_html(inline_html(cell)) if cell else trusted_html("&mdash;") for cell in row["cells"]]}
                     for row in table["rows"]]}


def glance_view(sheet):
    """The rows of a page's "At a glance" panel (UX.188): the kind's `describe.sheet.GLANCE` rows in that order, each
    `{"label", "html"}`, and only those the page lists. Empty for a kind with no panel yet."""
    rows = {row["key"]: row for group in sheet["groups"] for row in group["rows"]}
    view = []
    for key in sheet_registry.GLANCE.get(sheet["kind"], ()):
        if key not in rows or not _listed(rows[key]):
            continue
        if key == "flux_variation":     # GEN.1135: a flag, only when the swing is very high; the number is on the full sheet
            swing = rows[key].get("raw")
            if swing is None or swing < sheet_registry.SWING_FLAG_MIN:
                continue
            view.append({"label": sheet_registry.SWING_FLAG_LABEL, "html": trusted_html(f"Very high, {_value_html(rows[key])} peak to peak")})
            continue
        view.append({"label": rows[key]["label"], "html": _value_html(rows[key])})
    return view


def life_air_view(db, kind, printed_id, row):
    """The Life and air chart of a planet's or moon's page (UX.237): `web.lib.life_air_chart.chart` of the body's air history, or `None` for a
    body with no living class or whose air the life did not change, or when the API cannot say."""
    if kind not in ("planet", "moon") or not row or not row.get("life_class_id"):
        return None
    try:
        answer = apiclient.get_life_air(db, kind, printed_id)
    except (apiclient.NotFoundError, apiclient.ApiError):
        return None
    grade = int(row.get("life_animal_grade") or 0)
    view = life_air_chart.chart(answer.get("life_air"), row["life_class_id"], grade)
    return None if view is None else {**view, "svg": trusted_html(view["svg"]), "legend": trusted_html(view["legend"]),
                                      "note": trusted_html(view["note"])}


def life_readings_view(db, kind, printed_id, row):
    """The three words for each kind of life on a planet or moon (UX.251, `GET /api/planets/<id>/life-readings`), or `[]` for a body with no
    living class or when the API cannot say."""
    if kind not in ("planet", "moon") or not row or not row.get("life_class_id"):
        return []
    try:
        return apiclient.get_life_readings(db, kind, printed_id).get("readings") or []
    except (apiclient.NotFoundError, apiclient.ApiError):
        return []


def _life_class_help_url(class_id):
    """The address of the help page of life class `class_id` ("3a"), at its own variant."""
    number = int("".join(ch for ch in class_id if ch.isdigit()))
    return page_url("life_help_class", number=number, _anchor=f"class-{class_id}")


def life_class_with_links(row, readings):
    """`planet_layout.life_class_view` with the class names of its glance and paragraph linked to their help pages, marked safe for the
    template (the view escapes every value it writes)."""
    view = planet_layout.life_class_view(row, readings, _life_class_help_url)
    if view is not None:
        view = {**view, "glance": trusted_html(view["glance"]) if view["glance"] else None,
                "paragraph": trusted_html(view["paragraph"]) if view["paragraph"] else None}
    return view


def planet_sections(name, row, sheet, species, piece, readings=()):
    """
    The sections of a planet's or moon's page in the order Boss set (UX.215 to UX.222), as template variables in
    `planet`: the summary sentences (UX.220), the Life and civilization paragraphs when there is something to say
    (UX.221), the PHI-6 section (the rating, radiation and survival equipment tables and the habitability rows), the
    Atmosphere section (its rows and the table of the air by mole fraction) and the Composition section (its rows and the
    table of the body by mass), and `sheet`, the property groups left for the Properties section (UX.219).
    """
    parts = planet_layout.split_blocks(piece["blocks"] if piece else [])
    values = {r["key"]: _value_html(r) for g in sheet["groups"] for r in g["rows"] if r.get("stored", True)}
    raw = {r["key"]: r.get("raw") for g in sheet["groups"] for r in g["rows"] if r.get("stored", True)}
    groups = {group["key"]: group for group in sheet_view(sheet)}

    def without(key, drop):
        group = groups.pop(key, None)
        return None if group is None else {**group, "rows": [r for r in group["rows"] if not drop(r["key"])]}

    atmosphere = without("atmosphere", lambda key: key.startswith("p_"))
    composition = without("composition", lambda key: key in planet_layout.MASS_SHARED_KEYS)
    phi6 = without("habitability", lambda key: False)
    return {
        "summary": [trusted_html(sentence) for sentence in planet_summary.summary_sentences(name, row or {}, values, raw)],
        "life": trusted_html(blocks_html(parts["life"], _ref_link)),
        "life_class": life_class_with_links(row, readings),
        "civilization": [trusted_html(sentence) for sentence in planet_layout.civilization_sentences(row, species)],
        "notes": trusted_html(blocks_html(parts["rest"], _ref_link)),
        "phi6_blocks": trusted_html(blocks_html(parts["phi6"], _ref_link)),
        "phi6": phi6 if phi6 and phi6["rows"] else None,
        "atmosphere": atmosphere if atmosphere and atmosphere["rows"] else None,
        "bands": bands_view(sheet),
        "gases": [{**gas, "pressure": trusted_html(gas["pressure"])} for gas in planet_layout.gas_table(sheet)],
        "composition": composition if composition and composition["rows"] else None,
        "mass": planet_layout.mass_table(sheet),
        "sheet": [group for group in groups.values() if group["rows"]],
    }


def sheet_context(sheet):
    """The template variables of a page's sheet: `sheet` (the groups), `glance` (the At a glance rows) and
    `hidden_rows` (how many rows hold nothing and are not listed, for the line that says so) and `absences` (the sentences
    for the sections left out because they hold only zeros, UX.203)."""
    listed = sum(1 for group in sheet["groups"] for row in group["rows"] if _listed(row))
    total = sum(len(group["rows"]) for group in sheet["groups"])
    notes = [sheet_registry.MISSING_SENTENCES[key] for key in sheet_registry.MISSING_SENTENCES
             if key in (sheet.get("missing") or ())]
    return {"sheet": sheet_view(sheet), "glance": glance_view(sheet), "hidden_rows": total - listed, "absences": notes}


def region_row(found):
    """
    The row that names the charted region an object lies in (UX.231): "in Vang Zarh Region". `found` is an API answer that carries
    `region` (`{"name", "parent"}`; a sector, a system or `GET /api/objects/<ref>`); `None` when it has none, which is a position in
    no charted region (a cell with no generated sector has no region), so the row is not drawn. The row's html is escaped text.
    """
    region = found.get("region") or {}
    if not region.get("name"):
        return None
    return {"label": "Region", "html": esc(f"in {region['name']}")}


def with_region(context, found):
    """A page's sheet context with the region the object lies in as the last row of its At a glance panel (UX.231), when it lies
    in a charted one. `found` is `GET /api/objects/<ref>`."""
    row = region_row(found)
    if row is None:
        return context
    return {**context, "glance": context["glance"] + [{"label": row["label"], "html": trusted_html(row["html"])}]}


def _bodies(detail):
    """`[(kind, row)]` for every star, planet, moon, belt and comet of a system's detail, each row with the `kind` it is
    and (a moon's) its planet."""
    found = [("star", star) for star in detail["stars"]]
    for planet in detail["planets"]:
        found.append(("planet", planet))
        found.extend(("moon", dict(moon, planet=planet)) for moon in planet.get("moons") or [])
    found += [("belt", belt) for belt in detail["belts"]] + [("comet", comet) for comet in detail["comets"]]
    return found


def belt_name(belt):
    """A belt's name: it has none stored, so it is numbered by its place out from the star."""
    return f"Asteroid belt {belt['orbital_index'] + 1}"


def _name_of(kind, row):
    return belt_name(row) if kind == "belt" else row["name"]


def _find(detail, kind, printed_id):
    """The detail's row for a body, or `None`."""
    wanted = printed_id.upper()
    return next((row for found_kind, row in _bodies(detail) if found_kind == kind and str(row["id"]).upper() == wanted), None)


def _ref(kind, row):
    return f"{kind}:{row['id']}"


def chain_places(detail, kind, row):
    """The places between a body's system and the body itself, as `[{"ref", "kind", "name"}]`: the star it orbits (none
    for a close pair's circumbinary body) and, for a moon, its planet."""
    host_id = (row["planet"] if kind == "moon" else row).get("star_id")
    places = []
    star = next((s for s in detail["stars"] if s["id"] == host_id), None) if kind != "star" else None
    if star is not None:
        places.append({"ref": _ref("star", star), "kind": "star", "name": star["name"]})
    if kind == "moon":
        places.append({"ref": _ref("planet", row["planet"]), "kind": "planet", "name": row["planet"]["name"]})
    return places


def neighbours(detail, kind, row):
    """`(previous, next)` of the body at its own level, each `{"ref", "kind", "name"}` or `None`: the planets that share
    its star in their order out, a planet's moons, the system's stars, the belts or the comets of its star."""
    if kind == "star":
        peers = detail["stars"]
    elif kind == "planet":
        peers = [p for p in detail["planets"] if p.get("star_id") == row.get("star_id")]
    elif kind == "moon":
        peers = row["planet"].get("moons") or []
    elif kind == "belt":
        peers = [b for b in detail["belts"] if b.get("star_id") == row.get("star_id")]
    else:
        peers = [c for c in detail["comets"] if c.get("star_id") == row.get("star_id")]
    index = next((i for i, peer in enumerate(peers) if peer["id"] == row["id"]), None)
    if index is None:
        return None, None

    def entry(peer):
        return None if peer is None else {"ref": _ref(kind, peer), "kind": kind, "name": _name_of(kind, peer),
                                          "url": item_url(_ref(kind, peer))}

    return (entry(peers[index - 1]) if index > 0 else None), (entry(peers[index + 1]) if index + 1 < len(peers) else None)


def hosted_facilities(db, system_id, kind, row):
    """The facilities stored on this body (`GET /api/systems/<id>/facilities`), as the list its page shows."""
    found = facilities_by_host(apiclient.get_system_facilities(db, system_id)).get((kind, row["id"]))
    return trusted_html(facilities_html(found, lambda facility: item_url(f"facility:{facility['id']}")))


def dominant_species(db, kind, row):
    """`{"name", "url"}` of a life world's dominant species, or `None` (also before population has made any)."""
    if kind != "planet" or not row.get("life_class_id") or not population_status()["species"]:
        return None
    try:
        species = apiclient.get_planet_species(db, row["id"])
    except (apiclient.ApiError, apiclient.NotFoundError):
        return None
    return {"name": species["name"], "url": page_url("species_page", species_id=species["id"])} if species else None


def _admin_context(kind, printed_id, system_id):
    """What the Admin menu of a body's page needs: `{"target"}` for an admin allowed to edit, else `None`."""
    if not edit_actions.can_edit(current_admin()):
        return None
    return {"target": f"{kind}:{printed_id}"}


def _edit_post(kind, printed_id, system_id):
    """An admin Regenerate or Delete POST for this body (`web/edit_actions.py`); after a delete it goes up to
    the body's system."""
    here = item_url(f"{kind}:{printed_id}")
    return edit_actions.handle_post({(kind, printed_id)}, here,
                                    after_delete={kind: page_url("system", system_id=system_id)})


def _free_moon_page(db, kind, printed_id, found, free_planet):
    """The page of a moon of a free planet (DB.94): its sheet and places, with no system to read a description,
    neighbours or admin actions from, and a link to the free planet instead of a system."""
    ref = f"{kind}:{printed_id}"
    name = found["name"]
    nav_ref = endpoint(kind, printed_id)
    return render_page(
        "item.html",
        title=name,
        section="systems",
        breadcrumbs=_crumbs(found, name),
        description=f"The {KIND_LABELS[kind].lower()} {name}: where it is and what it is like.",
        kind=kind,
        kind_label=KIND_LABELS[kind],
        habitability=found if kind in ("planet", "moon") else None,
        object_id=printed_id,
        facts=trusted_html(""),
        **with_region(sheet_context(apiclient.get_object_sheet(db, ref)), found),
        tree=tree_panel(ref, name, apiclient.get_object_children(db, ref, limit=TREE_PAGE)),
        location=location(found["parents"], None),
        owner_url=item_url(free_planet["ref"]),
        owner_name=esc(free_planet["name"]),
        bookmark=bookmark(kind, nav_ref, name, item_url(ref)),
        nav_links={"from": nav_url(origin=nav_ref), "to": nav_url(destination=nav_ref)},
        nearby_url=nearby_url(nav_ref),
        admin=current_admin(),
        object_ref=ref,
    )


def _facility_facts(facility, found):
    """The About text of a facility (UX.183): what it is built on, its type and affiliation, its description and its
    rating chips, as HTML."""
    lines = []
    host = found["parents"][-1] if found["parents"] and found["parents"][-1]["kind"] not in ("galaxy", "sector", "system") \
        else None
    if host is not None:
        lines.append(f'<p>Built on <a href="{esc(item_url(host["ref"]))}">{esc(host["name"])}</a>.</p>')
    kind_words = esc(facility["kind"].replace("-", " "))
    typed = facility_display.type_html(facility)
    lines.append(f"<p>A {kind_words}{' of the type ' + typed if typed else ''}"
                 f"{', affiliated with ' + esc(facility_display.affiliation_text(facility)) if facility.get('affiliation_name') else ''}.</p>")
    if facility.get("description"):
        lines.append(f"<p>{esc(facility['description'])}</p>")
    lines.append(facility_display.rating_chips(facility))
    return "".join(lines)


def _facility_page(db, printed_id, found):
    """The page of a facility (UX.183): where it stands, its stored properties and what its record says. It holds
    nothing, so it has no tree; a facility not in a system has no system to link to."""
    ref = f"facility:{printed_id}"
    name = found["name"]
    facility = apiclient.get_facility(db, printed_id)
    system = next((p for p in found["parents"] if p["kind"] == "system"), None)
    nav_ref = endpoint("facility", printed_id)
    return render_page(
        "item.html",
        title=name,
        section="systems",
        breadcrumbs=_crumbs(found, name),
        description=f"The facility {name}: where it stands, what it is and everything stored about it.",
        kind="facility",
        kind_label=KIND_LABELS["facility"],
        object_id=printed_id,
        facts=trusted_html(_facility_facts(facility, found)),
        **with_region(sheet_context(apiclient.get_object_sheet(db, ref)), found),
        tree="",
        location=location(found["parents"], None),
        system_url=page_url("system", system_id=system["ref"].partition(":")[2]) if system else None,
        system_name=esc(system["name"]) if system else None,
        bookmark=bookmark("facility", nav_ref, name, item_url(ref)),
        nav_links=None,    # a course names a system, a body or a phenomenon, not a facility
        nearby_url=None,
        admin=current_admin(),
        object_ref=ref,
    )


_CONTAINER_WORDS = {
    "globular_cluster": "a globular cluster, a dense ball of old stars",
    "superbubble": "a superbubble, a cavity that supernovae and stellar winds have blown in the gas",
}

def _container_contents_html(listing):
    """The "What is inside it" lists of a container (UX.184), one linked list per kind, as HTML."""
    shown = {}
    for item in listing["items"]:
        shown.setdefault(item["kind"], []).append(item)
    parts = []
    for kind, label in CONTENT_LABELS.items():
        count = listing["counts"].get(kind)
        if not count:
            continue
        links = "".join(f'<li><a href="{esc(own_page_url(kind, item["id"]))}">{esc(item["name"])}</a></li>'
                        for item in shown.get(kind, []))
        more = count - len(shown.get(kind, []))
        tail = f"<li><em>and {more:,} more</em></li>" if more > 0 else ""
        parts.append(f"<h3>{esc(label)} ({count:,})</h3><ul>{links}{tail}</ul>")
    if not parts:
        return "<p>Nothing the galaxy stores lies inside it.</p>"
    return "<h3>What is inside it</h3>" + "".join(parts)


def _container_facts(container, found, listing):
    """The About text of a cluster or a superbubble (UX.184): what it is, where its centre lies, a cluster's mass, age
    and metallicity, and what sits inside it, as HTML."""
    kind = found["kind"]
    lines = [f"<p>This is {_CONTAINER_WORDS[kind]}. Its sphere is {container['extent_pc']:,.0f} parsecs "
             "in radius.</p>"]
    centre = found["positions"]["galaxy_pc"]
    if centre:
        lines.append("<p>Its centre lies at (" + ", ".join(f"{value:,.0f}" for value in centre)
                     + ") parsecs in the galaxy frame.</p>")
    component = container.get("component")
    if component:
        facts = [f"{component['mass_solar']:,.0f} solar masses" if component.get("mass_solar") else None,
                 f"{component['age_gy']:.1f} billion years old" if component.get("age_gy") is not None else None,
                 f"[Fe/H] {component['metallicity_fe_h']:+.2f}" if component.get("metallicity_fe_h") is not None else None]
        lines.append("<p>" + esc(", ".join(fact for fact in facts if fact).capitalize()) + ".</p>")
    lines.append(_container_contents_html(listing))
    return "".join(lines)


def _container_page(db, kind, printed_id, found):
    """The page of a globular cluster or a superbubble (UX.184): where it lies, its stored properties and what sits
    inside it. It has no tree and is no endpoint of a course."""
    ref = f"{kind}:{printed_id}"
    container = apiclient.get_container(db, printed_id)
    listing = apiclient.get_container_contents(db, printed_id)
    name = found["name"] or f"{KIND_LABELS[kind]} {printed_id}"
    return render_page(
        "item.html",
        title=name,
        section="galaxy",
        breadcrumbs=_crumbs(found, name),
        description=f"The {KIND_LABELS[kind].lower()} {name}: where it lies, what it is and what sits inside it.",
        kind=kind,
        kind_label=KIND_LABELS[kind],
        object_id=printed_id,
        facts=trusted_html(_container_facts(container, found, listing)),
        **with_region(sheet_context(apiclient.get_object_sheet(db, ref)), found),
        tree="",
        location=location(found["parents"], None),
        bookmark=bookmark(kind, endpoint(kind, printed_id), name, item_url(ref)),
        nav_links=None,
        nearby_url=None,
        admin=current_admin(),
        object_ref=ref,
    )


def _moon_system_url(system_id, kind, printed_id, row):
    """MAP.308: the System Map's moon system map of a planet that has moons (the system page opens the diagram that holds
    a `?scene=`); a planet without moons, and anything but a planet, has none."""
    if kind != "planet" or row is None or not row.get("moons"):
        return None
    return page_url("system", system_id=system_id, scene=f"planet-{row['id']}")


def render_item(kind, printed_id):
    """The page of one star, planet, moon, belt, comet, facility, globular cluster or superbubble; a missing one is a 404."""
    ref = f"{kind}:{printed_id}"
    post = request.method == "POST"
    db = db_name()
    try:
        found = apiclient.get_object(db, ref)
    except apiclient.NotFoundError:
        abort(404)
    if kind == "facility":
        return _facility_page(db, printed_id, found)
    if kind in objectref.CONTAINER_KINDS:
        return _container_page(db, kind, printed_id, found)
    system = next((p for p in found["parents"] if p["kind"] == "system"), None)
    if system is None:
        free_planet = next((p for p in found["parents"] if p["kind"] == "rogue_planet"), None)
        if free_planet is None:
            abort(404)
        return _free_moon_page(db, kind, printed_id, found, free_planet)    # DB.94: no system, so no description
    system_id = system["ref"].partition(":")[2]
    if post and "edit_action" in request.form:
        return _edit_post(kind, printed_id, system_id)
    detail = apiclient.get_system(db, system_id)
    row = _find(detail, kind, printed_id)
    document = apiclient.get_system_description(db, system_id)
    piece = _section(document, ref)
    name = found["name"]
    label = KIND_LABELS[kind]
    parents = found["parents"]
    chain_parents = parents
    if row is not None:
        chain_parents = [p for p in parents if p["kind"] != "planet"] + chain_places(detail, kind, row)
    previous, following = neighbours(detail, kind, row) if row is not None else (None, None)
    sheet = apiclient.get_object_sheet(db, ref)
    species = dominant_species(db, kind, row) if row is not None else None
    children = apiclient.get_object_children(db, ref, limit=TREE_PAGE)
    tree = tree_panel(ref, name, children)
    admin = _admin_context(kind, printed_id, system_id)
    nav_ref = endpoint(kind, printed_id)
    return render_page(
        "item.html",
        title=name,
        section="systems",
        breadcrumbs=_crumbs({"parents": chain_parents}, name),
        description=f"The {label.lower()} {name}: where it is, what it is like and everything stored about it.",
        kind=kind,
        kind_label=label,
        habitability={**found, "life_class_id": (row or {}).get("life_class_id")} if kind in ("planet", "moon") else None,
        chips=trusted_html(chips_html(row, page_url("habitability_levels"))) if kind in ("planet", "moon") and row else "",
        facilities=hosted_facilities(db, system_id, kind, row) if row is not None else "",
        species=species,
        object_id=printed_id,
        facts=trusted_html(blocks_html(piece["blocks"], _ref_link) if piece else ""),
        example=bodyview.view(kind, row, detail["stars"]) if row is not None else None,
        life_air=life_air_view(db, kind, printed_id, row),
        **with_region(sheet_context(sheet), found),
        **({"planet": planet_sections(name, row, sheet, species, piece, life_readings_view(db, kind, printed_id, row))} if kind in ("planet", "moon") else {}),
        tree=tree,
        tree_marks=tree_has_marks(children),
        location=location(chain_parents, detail.get("inside")),
        previous=previous,
        following=following,
        system_url=page_url("system", system_id=system_id),
        system_name=esc(system["name"]),
        moon_system_url=_moon_system_url(system_id, kind, printed_id, row),
        bookmark=bookmark(kind, nav_ref, name, item_url(ref)),
        nav_links={"from": nav_url(origin=nav_ref), "to": nav_url(destination=nav_ref)},
        nearby_url=nearby_url(nav_ref),
        admin=current_admin(),
        edit=admin,
        object_ref=ref,
    )


def object_page_response(ref):
    """The response for `/object/<ref>`: the redirect to a page the site already has, the item page, or a 404."""
    try:
        kind, printed_id = objectref.parse_public(ref)
    except ValueError:
        abort(404)
    if kind in OWN_PAGE_KINDS:
        return redirect(own_page_url(kind, printed_id), code=302)
    return render_item(kind, printed_id)
