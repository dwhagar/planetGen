# planetgen/web/belt_member_page.py

"""
The page of one body an asteroid belt keeps, `/object/belt:<ID>/member/<rank>` (GEN.1070): its picture (MAP.413: the look of `bodylook.member_look` painted by `static/bodyview.js`), designation, status, mass, size, makeup, orbit and the history of its belt. The data is
`GET /api/objects/<ref>/belt-members/<rank>` (`describe/belt_member_info.py`); this page lays it out.
"""

from flask import abort

from planetgen.galaxy import objectref, object_uid
from planetgen.web.lib import apiclient, bodyview

from . import bp
from .helpers import crumb, db_name, page_url, render_page
from .item_pages import _crumbs, item_url


def member_url(belt_ref, designation):
    """The page of the member a designation (`<belt ID>.<rank hex>`) names, in the belt `belt_ref`."""
    return page_url("belt_member_page", belt_ref=belt_ref, rank=designation.rpartition(object_uid.MEMBER_SEPARATOR)[2])


def _quantity(value):
    return value["text"] if value else ""


@bp.route("/object/<belt_ref>/member/<rank>")
def belt_member_page(belt_ref, rank):
    """One stored body of an asteroid belt; a reference that is not a belt, a rank the belt does not keep or one written otherwise than the designation writes it is a 404."""
    try:
        kind, printed_id = objectref.parse_public(belt_ref)
        object_uid.parse_rank(rank)
    except ValueError:
        abort(404)
    if kind != "belt":
        abort(404)
    ref, db = f"belt:{printed_id}", db_name()
    try:
        found = apiclient.get_object(db, ref)
        data = apiclient.get_belt_member(db, ref, rank)
    except apiclient.NotFoundError:
        abort(404)
    name = found["name"]
    crumbs = _crumbs(found, name)
    crumbs[-1] = {"label": name, "url": item_url(ref)}
    crumbs.append(crumb("Bodies and sizes", "belt_info_page", belt_ref=ref))
    crumbs.append(crumb(data["designation"]))
    member = data["member"]
    belt = data["belt"]
    return render_page(
        "belt_member.html",
        title=f"{data['designation']}: a body of {name}",
        section="systems",
        breadcrumbs=crumbs,
        description=f"{data['designation']}: one of the heaviest bodies of {name}, its mass, size, makeup and orbit.",
        designation=data["designation"],
        name=name,
        belt_url=item_url(ref),
        belt_info_url=page_url("belt_info_page", belt_ref=ref),
        belt=belt,
        member=member,
        composition_text=", ".join(f"{part} {member['composition'][part]['text']}" for part in ("rock", "metal", "water", "ice") if member["composition"][part]["value"] > 0.0),
        orbit=data["orbit"],
        share=data["share"],
        previous_url=member_url(ref, data["previous"]) if data["previous"] else None,
        next_url=member_url(ref, data["next"]) if data["next"] else None,
        previous=data["previous"],
        next=data["next"],
        quantity=_quantity,
        example=bodyview.member_view(data["designation"], member),
    )
