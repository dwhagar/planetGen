# planetgen/web/lib/bodyview.py

"""
The Sensor view on an object's page (UX.171, UX.276): the look record (`web/maps/bodylook.py`, MAP.271) of a star,
planet or moon, handed to `static/bodyview.js` as `data-look`, which paints the same texture the system map paints
and shades it as a lit sphere on a canvas.

A belt or a comet has no look record yet, so its page has no example; a belt's member has one (`member_view`, MAP.413). The caption says how the picture was worked
out (`look_note`, a stand-in until the colour account of UX.158 / MAP.271 replaces it).
"""

from flask import abort

from planetgen.galaxy import objectref, object_uid
from planetgen.web.lib import sensor
from planetgen.web.maps.bodylook import look_json, look_note, member_look, planet_look, star_look


def view(kind, row, stars):
    """`{"look", "note", "estimate"}` for a star, planet or moon row of a system's detail (`stars` are the system's star rows, for the
    temperature of the star a body's air is seen under), or `None` for a kind the Sensor registry (`web/lib/sensor.py`) has no `sphere` painter
    for. `estimate` lists what the picture guessed (UX.276)."""
    entry = sensor.REGISTRY.get(kind)
    if entry is None or entry.painter != "sphere":
        return None
    if kind == "star":
        look = star_look(row)
    elif kind in ("planet", "moon"):
        planet = row["planet"] if kind == "moon" else row
        host = next((s for s in stars if s["id"] == planet.get("star_id")), stars[0] if stars else None)
        look = planet_look(dict(row, _star_temperature_k=host.get("temperature_k") if host else None), kind)
    else:
        return None
    return {"look": look_json(look), "note": look_note(look), "estimate": sensor.estimates_of(look)}


def member_view(designation, member):
    """
    `{"look", "note", "estimate"}` for the Sensor view of a belt member's page (MAP.413): the look record of `bodylook.member_look` over the page's `member` block
    (`physics/belt_page.stored_bodies`: quantities in the units the page shows), its id the body's designation so the page and the 3D view (MAP.414) draw the same
    body. Nothing is guessed, so `estimate` is empty.
    """
    value = lambda name: member[name]["value"]
    water, ice = member["composition"]["water"]["value"], member["composition"]["ice"]["value"]
    look = member_look({
        "id": designation, "taxon": member["taxon"], "mass_kg": value("mass"), "radius_km": value("diameter") / 2.0, "density_g_cm3": value("density") / 1000.0,
        "rock_fraction": member["composition"]["rock"]["value"], "metal_fraction": member["composition"]["metal"]["value"], "water_fraction": water, "ice_fraction": ice})
    return {"look": look_json(look), "note": look_note(look), "estimate": []}


MEMBER_COMPOSITION_PARTS = ("rock", "metal", "water", "ice")
"""tuple: The parts of a stored small body's composition, in the order its page lists them."""


def member_page_ref(public_ref, rank, kind):
    """`(api ref, printed id)` of the `kind` of object a member page's `public_ref` names, with `rank` checked to be written as a
    designation writes it; a 404 for a reference of another kind, one that does not parse or a rank written otherwise."""
    try:
        found, printed_id = objectref.parse_public(public_ref)
        object_uid.parse_rank(rank)
    except ValueError:
        abort(404)
    if found != kind:
        abort(404)
    return f"{kind}:{printed_id}", printed_id


def quantity_text(value):
    """The `text` of an API quantity (`{"value", "text", ...}`), or `""` for a missing one."""
    return value["text"] if value else ""


def member_composition_text(entry):
    """The parts of a stored small body's composition that are present, as `"rock 60 %, metal 10 %"`: the one line the
    belt member, field member and irregular moon pages all show."""
    return ", ".join(f"{part} {entry['composition'][part]['text']}" for part in MEMBER_COMPOSITION_PARTS
                     if entry["composition"][part]["value"] > 0.0)
