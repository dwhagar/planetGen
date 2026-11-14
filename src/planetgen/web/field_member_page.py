# planetgen/web/field_member_page.py

"""
The page of one body an interstellar asteroid field keeps, `/object/asteroid_field:<ID>/field-member/<rank>` (GEN.1105): its designation, status, mass, size and makeup, and where its field came from. The
data is `GET /api/objects/<ref>/field-members/<rank>` (`describe/field_member_info.py`); this page lays it out. The field's own page lists its bodies (`field_members_panel`).
"""

from flask import abort

from planetgen.galaxy import object_uid
from planetgen.web.lib import apiclient, bodyview

from . import bp
from .helpers import crumb, db_name, page_url, render_page


def member_url(field_ref, designation):
    """The page of the body a designation (`<field ID>.<rank hex>`) names, in the field `field_ref`."""
    rank = object_uid.parse_field_designation(designation)[1]
    return page_url("field_member_page", field_ref=field_ref, rank=f"{rank:0{object_uid.MEMBER_RANK_DIGITS}X}")


def field_members_panel(field_ref):
    """What the field's page shows of its bodies: `{"count", "rows"}` (each row's `designation`, `url`, `mass`, `diameter` and `status` text), or `None` for a field that keeps none."""
    try:
        data = apiclient.get_field_members(db_name(), field_ref)
    except (apiclient.ApiError, apiclient.NotFoundError):
        return None   # the list is extra: the page opens without it
    if not data["count"]:
        return None
    return {"count": data["count"], "rows": [
        {"designation": body["designation"], "url": member_url(field_ref, body["designation"]), "mass": bodyview.quantity_text(body["mass"]), "diameter": bodyview.quantity_text(body["diameter"]),
         "status": body["status"]["word"] or ""} for body in data["bodies"]]}


@bp.route("/object/<field_ref>/field-member/<rank>")
def field_member_page(field_ref, rank):
    """One stored body of an asteroid field; a reference that is not a field, a rank the field does not keep or one written otherwise than the designation writes it is a 404."""
    ref, printed_id = bodyview.member_page_ref(field_ref, rank, "asteroid_field")
    db = db_name()
    try:
        found = apiclient.get_object(db, ref)
        data = apiclient.get_field_member(db, ref, rank)
    except apiclient.NotFoundError:
        abort(404)
    name = found["name"]
    member = data["member"]
    field = data["field"]
    return render_page(
        "field_member.html",
        title=f"{data['designation']}: a body of {name}",
        section="phenomena",
        breadcrumbs=[crumb("Phenomena", "phenomena"), crumb(name, "phenomenon", phenomenon_type="asteroid_field", phenomenon_id=printed_id), crumb(data["designation"])],
        description=f"{data['designation']}: one of the heaviest bodies of {name}, its mass, size and makeup.",
        designation=data["designation"],
        name=name,
        field_url=page_url("phenomenon", phenomenon_type="asteroid_field", phenomenon_id=printed_id),
        field=field,
        member=member,
        composition_text=bodyview.member_composition_text(member),
        share=data["share"],
        previous_url=member_url(ref, data["previous"]) if data["previous"] else None,
        next_url=member_url(ref, data["next"]) if data["next"] else None,
        previous=data["previous"],
        next=data["next"],
        quantity=bodyview.quantity_text,
    )
