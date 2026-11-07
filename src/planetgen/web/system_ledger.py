# planetgen/web/system_ledger.py

"""
The system page's "Mass and barycenter" panel (GEN.1241): the mass of the whole system, where its center of mass lies against the primary star, how far the asteroid belts move that
star, the mass and place of each body, and the mass each belt's stored bodies and remainder add up to. The API (`GET /api/systems/<id>/bodies`, `api/system_bodies.py`) builds the
ledger from the typed tables; this builds the rows.
"""

import math

from planetgen.web.lib import apiclient
from planetgen.web.lib.fmt import format_number
from planetgen.web.lib.tabledisplay import format_body_distance, format_mass, format_star_radius

from .helpers import page_url

SHOWN_KINDS = ("star", "planet", "moon", "ring", "comet")
"""tuple: The kinds the table lists one row each; the bodies of a belt are summed in the belt table."""

_AU_KM = 149597870.7


def _link(entry):
    """The page of a body with a page of its own, else `None`."""
    ref = entry.get("ref")
    return page_url("object_page", ref=ref) if ref else None


def _share(part_km, radius_km):
    """A length as a share of the primary star's radius, as words."""
    if not radius_km or radius_km <= 0:
        return None
    ratio = part_km / radius_km
    if ratio < 0.001:
        return f"{ratio:.1e} of the star's radius"
    return f"{ratio:.3g} of the star's radius"


def panel(db, system_id):
    """
    `{summary, bodies, belts}` for the panel, or `None` for a system with no mass to put in a ledger (or one the API cannot read: the panel is extra and the page still opens).

    `summary` is a list of `(label, value, note)`; `bodies` rows are `{name, url, kind, parent, mass, distance}`; `belts` rows are `{name, url, mass, kept, extended, check}`.
    """
    try:
        data = apiclient.get_system_bodies(db, system_id)
    except (apiclient.ApiError, apiclient.NotFoundError):
        return None
    if not data["bodies"] or not data["total_mass_kg"]:
        return None
    radius = data.get("central_radius_km")
    names = {entry["ref"]: entry["name"] for entry in data["bodies"] if entry.get("ref")}
    summary = [("Mass of the system", format_mass(data["total_mass_kg"]), "every star, planet, moon, ring, comet and belt body, once"),
               ("Center of mass from the primary star", format_body_distance(data["central_offset_km"]), _share(data["central_offset_km"], radius))]
    if data["belts"]:
        summary.append(("Pull of the belts on the star", format_body_distance(data["belt_wobble_km"]) if data["belt_wobble_km"] else "none",
                        _share(data["belt_wobble_km"], radius) if data["belt_wobble_km"] else "the center of mass is the same with the belts left out"))
    bodies = []
    for entry in data["bodies"]:
        if entry["kind"] not in SHOWN_KINDS:
            continue
        bodies.append({"name": entry["name"], "url": _link(entry), "kind": entry["kind"].capitalize(), "parent": names.get(entry["parent"]) if entry["parent"] else None,
                       "mass": format_mass(entry["mass_kg"]) if entry["mass_kg"] else "–",
                       "distance": format_body_distance(entry["distance_from_barycenter_au"] * _AU_KM)})
    belts = []
    for belt in data["belts"]:
        residual = belt["residual"]
        belts.append({"name": belt["name"], "url": page_url("object_page", ref=belt["belt"]),
                      "mass": format_mass(belt["mass_kg"]) if belt["mass_kg"] else "–",
                      "kept": f"{belt['members']:,} bodies, {format_mass(belt['stored_members_kg'])}" if belt["members"] else "none stored",
                      "extended": format_mass(belt["extended_kg"]) if belt["extended_kg"] else "–",
                      "check": "–" if residual is None else ("adds up" if math.isclose(residual, 0.0, abs_tol=1e-9) else f"off by {format_number(residual * 100.0, '.2g')} percent")})
    return {"summary": summary, "bodies": bodies, "belts": belts}
