# planetgen/web/belt_info_page.py

"""
The information page of an asteroid belt, `/object/belt:<ID>/belt` (UX.241): the bodies the belt keeps one by one, by designation; the rest of its mass by size, with what it is made of;
and how crowded it is. The data is `GET /api/objects/<ref>/belt-page` (`describe/belt_info.py`); this page lays it out.
"""

from flask import abort

from planetgen.galaxy import objectref
from planetgen.web.lib import apiclient

from . import bp
from .helpers import crumb, db_name, render_page
from .item_pages import _crumbs, item_url


def _quantity(value):
    return value["text"] if value else ""


def _composition(shares):
    return ", ".join(f"{entry['analogue']} {entry['share']['text']}" for entry in shares)


@bp.route("/object/<belt_ref>/belt")
def belt_info_page(belt_ref):
    """The bodies, sizes and density of one asteroid belt; a reference that is not a belt, or names nothing, is a 404."""
    try:
        kind, printed_id = objectref.parse_public(belt_ref)
    except ValueError:
        abort(404)
    if kind != "belt":
        abort(404)
    ref, db = f"belt:{printed_id}", db_name()
    try:
        found = apiclient.get_object(db, ref)
        data = apiclient.get_belt_page(db, ref)
    except apiclient.NotFoundError:
        abort(404)
    name = found["name"]
    crumbs = _crumbs(found, name)
    crumbs[-1] = {"label": name, "url": item_url(ref)}
    crumbs.append(crumb("Bodies and sizes"))
    belt = data["belt"]
    return render_page(
        "belt_info.html",
        title=f"{name}: bodies and sizes",
        section="systems",
        breadcrumbs=crumbs,
        description=f"{name}: the bodies it keeps one by one, the sizes and contents of the rest of its mass, and how crowded it is.",
        name=name,
        belt_url=item_url(ref),
        belt=belt,
        belt_composition=_composition(belt["composition"]),
        bodies=[{**body, "composition_text": ", ".join(f"{part} {body['composition'][part]['text']}" for part in ("rock", "metal", "water", "ice")
                                                     if body["composition"][part]["value"] > 0.0)} for body in data["bodies"]],
        bins=[{**row, "composition_text": _composition(row["composition"])} for row in data["bins"]],
        density=data["density_factor"],
        totals=data["totals"],
        quantity=_quantity,
    )
