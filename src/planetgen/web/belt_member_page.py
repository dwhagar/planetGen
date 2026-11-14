# planetgen/web/belt_member_page.py

"""
The page of one body an asteroid belt keeps, `/object/belt:<ID>/member/<rank>` (GEN.1070): its picture (MAP.413: the look of `bodylook.member_look` painted by `static/bodyview.js`), designation, status, mass, size, makeup, orbit and the history of its belt. The data is
`GET /api/objects/<ref>/belt-members/<rank>` (`describe/belt_member_info.py`); this page lays it out.
"""

from flask import abort

from planetgen.galaxy import object_uid
from planetgen.web.lib import apiclient, bodyview

from . import bp
from .helpers import crumb, db_name, page_url, render_page
from .item_pages import item_url, linked_crumbs


def member_url(belt_ref, designation):
    """The page of the member a designation (`<belt ID>.<rank hex>`) names, in the belt `belt_ref`."""
    return page_url("belt_member_page", belt_ref=belt_ref, rank=designation.rpartition(object_uid.MEMBER_SEPARATOR)[2])


@bp.route("/object/<belt_ref>/member/<rank>")
def belt_member_page(belt_ref, rank):
    """One body of an asteroid belt at the class D size or larger (stored, or made when asked for, GEN.1674); a reference that is not a belt, a rank the belt does not keep or one written otherwise than the designation writes it is a 404."""
    ref, printed_id = bodyview.member_page_ref(belt_ref, rank, "belt")
    db = db_name()
    try:
        found = apiclient.get_object(db, ref)
        data = apiclient.get_belt_member(db, ref, rank)
    except apiclient.NotFoundError:
        abort(404)
    name = found["name"]
    crumbs = linked_crumbs(found, name, ref)
    crumbs.append(crumb("Bodies and sizes", "belt_info_page", belt_ref=ref))
    crumbs.append(crumb(data["designation"]))
    member = data["member"]
    belt = data["belt"]
    return render_page(
        "belt_member.html",
        title=f"{data['designation']}: a body of {name}",
        section="systems",
        breadcrumbs=crumbs,
        description=f"{data['designation']}: a body of {name} the size of a class D body or larger, its class, mass, size, makeup and orbit.",
        designation=data["designation"],
        name=name,
        belt_url=item_url(ref),
        belt_info_url=page_url("belt_info_page", belt_ref=ref),
        belt=belt,
        planet=data.get("planet"),
        member=member,
        composition_text=bodyview.member_composition_text(member),
        orbit=data["orbit"],
        share=data["share"],
        previous_url=member_url(ref, data["previous"]) if data["previous"] else None,
        next_url=member_url(ref, data["next"]) if data["next"] else None,
        previous=data["previous"],
        next=data["next"],
        quantity=bodyview.quantity_text,
        example=bodyview.member_view(data["designation"], member),
    )
