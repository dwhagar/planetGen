# planetgen/web/irregular_moon_pages.py

"""
The irregular moons of a planet, `/object/planet:<ID>/irregular-moons`, and the page of one, `/object/planet:<ID>/irregular/<rank>` (GEN.1086): the list by designation with each body's mass, size, makeup
and orbit about the planet, and a page per body with its picture, orbit and the capture window it came through. The data is `GET /api/objects/<ref>/irregular-moons` and `.../irregular-moons/<rank>`
(`describe/irregular_moon_info.py`); these pages lay it out. They follow the belt pages (`belt_info_page`, `belt_member_page`).
"""

from flask import abort

from planetgen.galaxy import objectref, object_uid
from planetgen.web.lib import apiclient, bodyview

from . import bp
from .helpers import crumb, db_name, page_url, render_page
from .item_pages import _crumbs, item_url

WINDOW_WORDS = {"gas": "the gas of the planet's disk", "encounter": "a close encounter between giants", "exchange": "an exchange with a binary"}
"""dict: The capture window (`captured_composition.WINDOWS`) in the words the page uses."""


def irregular_moons_url(planet_ref):
    """The list of the irregular moons of the planet `planet_ref`."""
    return page_url("irregular_moons_page", planet_ref=planet_ref)


def irregular_moon_url(planet_ref, rank):
    """The page of the irregular moon of the planet `planet_ref` whose designation rank is `rank` (a whole number)."""
    return page_url("irregular_moon_page", planet_ref=planet_ref, rank=f"{rank:0{object_uid.MEMBER_RANK_DIGITS}X}")


def _quantity(value):
    return value["text"] if value else ""


def _composition_text(entry):
    return ", ".join(f"{part} {entry['composition'][part]['text']}" for part in ("rock", "metal", "water", "ice") if entry["composition"][part]["value"] > 0.0)


def _planet_ref(planet_ref):
    try:
        kind, printed_id = objectref.parse_public(planet_ref)
    except ValueError:
        abort(404)
    if kind != "planet":
        abort(404)
    return f"planet:{printed_id}"


@bp.route("/object/<planet_ref>/irregular-moons")
def irregular_moons_page(planet_ref):
    """The irregular moons a planet keeps; a reference that is not a planet, or names nothing, is a 404."""
    ref, db = _planet_ref(planet_ref), db_name()
    try:
        found = apiclient.get_object(db, ref)
        data = apiclient.get_irregular_moons(db, ref)
    except apiclient.NotFoundError:
        abort(404)
    name = found["name"]
    crumbs = _crumbs(found, name)
    crumbs[-1] = {"label": name, "url": item_url(ref)}
    crumbs.append(crumb("Irregular moons"))
    return render_page(
        "irregular_moons.html",
        title=f"{name}: irregular moons",
        section="systems",
        breadcrumbs=crumbs,
        description=f"{name}: the small bodies it caught, their masses, sizes, makeup and orbits.",
        name=name,
        planet_url=item_url(ref),
        planet=data["planet"],
        not_stored=data["not_stored"],
        moons=[{**moon, "url": irregular_moon_url(ref, moon["rank"]), "composition_text": _composition_text(moon), "window_text": WINDOW_WORDS[moon["window"]]} for moon in data["moons"]],
        quantity=_quantity,
    )


@bp.route("/object/<planet_ref>/irregular/<rank>")
def irregular_moon_page(planet_ref, rank):
    """One irregular moon of a planet; a reference that is not a planet, a rank the planet does not keep or one written otherwise than a designation writes it is a 404."""
    ref, db = _planet_ref(planet_ref), db_name()
    try:
        object_uid.parse_rank(rank)
    except ValueError:
        abort(404)
    try:
        found = apiclient.get_object(db, ref)
        data = apiclient.get_irregular_moon(db, ref, rank)
    except apiclient.NotFoundError:
        abort(404)
    name = found["name"]
    crumbs = _crumbs(found, name)
    crumbs[-1] = {"label": name, "url": item_url(ref)}
    crumbs.append(crumb("Irregular moons", "irregular_moons_page", planet_ref=ref))
    crumbs.append(crumb(data["designation"]))
    moon = data["moon"]
    return render_page(
        "irregular_moon.html",
        title=f"{data['designation']}: an irregular moon of {name}",
        section="systems",
        breadcrumbs=crumbs,
        description=f"{data['designation']}: one of the irregular moons of {name}, its mass, size, makeup and orbit.",
        designation=data["designation"],
        name=name,
        planet_url=item_url(ref),
        list_url=irregular_moons_url(ref),
        planet=data["planet"],
        moon=moon,
        window_text=WINDOW_WORDS[moon["window"]],
        composition_text=_composition_text(moon),
        orbit=data["orbit"],
        previous_url=irregular_moon_url(ref, data["previous"]["rank"]) if data["previous"] else None,
        next_url=irregular_moon_url(ref, data["next"]["rank"]) if data["next"] else None,
        previous=data["previous"]["designation"] if data["previous"] else None,
        next=data["next"]["designation"] if data["next"] else None,
        quantity=_quantity,
        example=bodyview.member_view(data["designation"], moon),
    )
