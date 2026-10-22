# planetgen/web/item_pages.py

"""
The dedicated page of an item (UX.143): `/object/<ref>` for any object reference (`planet:FE81000A2B-0000005-001`, a
bare system ID). Every kind a reference can name has a page, from one registry:

- a sector, a system and each phenomenon type already have a page of their own (`/sector`, `/system`, `/phenomenon`),
  so their reference goes there;
- a star, planet, moon, belt or comet gets the page here: its name, kind and Object ID, the chain of places it is in
  (galaxy, sector, system, the planet a moon circles), everything the system's description says about it (the facts,
  ratings, composition and prose, with units from the unit ladder), the bodies that belong to it, and a link back to its
  row in the system's list.

The facts are the description's section for the body (`GET /api/systems/<id>/description`, `describe/`), so the page
and the system's row never say different things. A reference that names nothing is a plain not-found page.
"""

from flask import abort, redirect

from planetgen.describe.html import blocks_html
from planetgen.describe.model import sections
from planetgen.galaxy import objectref
from planetgen.web.lib import apiclient
from planetgen.web.lib.fmt import esc

from planetgen.web.lib.systempage import body_anchor

from .helpers import crumb, db_name, page_url, render_page, trusted_html

OWN_PAGE_KINDS = ("sector", "system") + objectref.PHENOMENON_KINDS
"""tuple: The kinds whose reference opens a page the site already has."""

ITEM_PAGE_KINDS = objectref.BODY_KINDS
"""tuple: The kinds this module renders."""

KIND_LABELS = {"star": "Star", "planet": "Planet", "moon": "Moon", "belt": "Asteroid belt", "comet": "Comet"}


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
    """Home > Galaxy > sector > system > ... > this item, from the API's parent chain."""
    crumbs = []
    for parent in found["parents"]:
        if parent["kind"] == "galaxy":
            crumbs.append(crumb("Galaxy", "galaxy"))
        else:
            crumbs.append({"label": parent["name"], "url": item_url(parent["ref"])})
    crumbs.append(crumb(name))
    return crumbs


def _section(document, ref):
    """The description's section for `ref`, or `None`."""
    wanted = ref.lower()
    return next((piece for piece in sections(document) if piece["ref"] and piece["ref"].lower() == wanted), None)


def _ref_link(ref):
    """The link a description reference (`planet:ID`) gets: its item page."""
    return item_url(ref) if ref else None


def render_item(kind, printed_id):
    """The page of one star, planet, moon, belt or comet; a missing one is a 404."""
    ref = f"{kind}:{printed_id}"
    db = db_name()
    try:
        found = apiclient.get_object(db, ref)
    except apiclient.NotFoundError:
        abort(404)
    system = next((p for p in found["parents"] if p["kind"] == "system"), None)
    if system is None:
        abort(404)
    system_id = system["ref"].partition(":")[2]
    document = apiclient.get_system_description(db, system_id)
    piece = _section(document, ref)
    name = found["name"]
    label = KIND_LABELS[kind]
    children = [{"title": child["title"], "kind": KIND_LABELS.get(child["kind"], child["kind"].title()),
                 "url": _ref_link(child["ref"])} for child in (piece["children"] if piece else []) if child["ref"]]
    return render_page(
        "item.html",
        title=name,
        section="systems",
        breadcrumbs=_crumbs(found, name),
        description=f"The {label.lower()} {name}: where it is and what it is like.",
        kind_label=label,
        object_id=printed_id,
        facts=trusted_html(blocks_html(piece["blocks"], _ref_link) if piece else ""),
        children=children,
        system_url=page_url("system", system_id=system_id),
        row_url=page_url("system", _anchor=body_anchor(kind, printed_id), system_id=system_id),
        system_name=esc(system["name"]),
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
