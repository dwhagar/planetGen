# planetgen/web/glance.py

"""
The "At a glance" rows of the pages that are not an item page (UX.188 to UX.194): a system's, and later a sector's and the
galaxy's. An item page reads its rows from the property sheet (`describe/sheet.py` `GLANCE`); a system or a sector has no
sheet, so its rows are worked out here from what its page already holds, in the order a captain asks: what is there,
where is the best place to stop, what is the hazard.

A row that would say none or 0 is left out unless its absence is the answer (a system with no planets says so).
"""

from planetgen.describe import model
from planetgen.describe.hierarchy import hierarchy_sentences
from planetgen.describe.html import inline_html
from planetgen.web import item_pages
from planetgen.web.helpers import trusted_html
from planetgen.web.lib import habitability_symbols
from planetgen.web.lib.fmt import esc

GOOD_TIER = 1
"""int: The worst overall PHI-6 tier counted as a habitable world: Blue (0) or Green (1)."""


def _count(word, plural, number):
    return f"{number:,} {word if number == 1 else plural}"


def _worlds(detail):
    """`[(kind, body)]` for every planet and moon of a system's detail."""
    worlds = []
    for planet in detail.get("planets") or []:
        worlds.append(("planet", planet))
        worlds.extend(("moon", moon) for moon in planet.get("moons") or [])
    return worlds


def system_rows(detail, facilities, symbols):
    """
    The rows of a system page's panel: `[{"label", "html"}]`.

    Args:
        detail (dict): `GET /api/systems/<id>`.
        facilities (dict or list): `GET /api/systems/<id>/facilities` (its `items`, or the list itself).
        symbols (callable): Draws the symbols of a scored body, `symbols(body, page=<its page's URL>)`
            (`web.habitability_symbol_row`).
    """
    rows = []
    stars = detail.get("stars") or []
    if stars:
        primary = next((s for s in stars if s.get("role") in ("single", "primary")), stars[0])
        make_up = hierarchy_sentences(
            [{"name": s.get("name") or "", "type": s.get("star_type"), "node_no": s.get("parent_node_no"), "slot": s.get("slot") or 0}
             for s in stars],
            [{"node_no": n["node_no"], "parent_node_no": n["parent_node_no"], "slot": n["slot"], "kind": n["kind"]}
             for n in detail.get("orbits") or []])
        rows.append({"label": "Stars", "html": esc(" ".join(make_up) if make_up else f"{_count('star', 'stars', len(stars))}, the main one {primary.get('star_type') or 'unclassified'}")})
    planets = detail.get("planets") or []
    rows.append({"label": "Planets", "html": esc(_count("planet", "planets", len(planets)) if planets else "none")})
    moons = sum(len(planet.get("moons") or []) for planet in planets)
    for label, number in (("Moons", moons), ("Asteroid belts", len(detail.get("belts") or [])),
                          ("Comets", len(detail.get("comets") or []))):
        if number:
            rows.append({"label": label, "html": esc(f"{number:,}")})
    scored = [(habitability_symbols.overall_tier(body), kind, body) for kind, body in _worlds(detail)
              if habitability_symbols.overall_tier(body) is not None]
    good = [item for item in scored if item[0] <= GOOD_TIER]
    rows.append({"label": "Habitable worlds", "html": esc(f"{len(good):,}" if good else "none")})
    if scored:
        tier, kind, body = min(scored, key=lambda item: item[0])
        url = item_pages.item_url(f"{kind}:{body['id']}")
        rows.append({"label": "Best world", "html": f'<a href="{esc(url)}">{esc(body["name"])}</a> {symbols(body, page=url)}'})
    held = facilities.get("items") if isinstance(facilities, dict) else facilities
    if held:
        rows.append({"label": "Facilities", "html": esc(f"{len(held):,}")})
    perimeters = [s["system_perimeter_km"] for s in stars if s.get("system_perimeter_km")]
    if perimeters:
        rows.append({"label": "System perimeter", "html": inline_html([model.qty("length", max(perimeters), "km")])})
    region = item_pages.region_row(detail)
    if region:
        rows.append(region)
    return [{"label": row["label"], "html": trusted_html(row["html"])} for row in rows]


HAZARD_TYPES = ("black_hole", "neutron_star", "quasar")
"""tuple: The standalone phenomena a ship must steer clear of (they have a keep-out radius): named in a sector's hazard row."""


def sector_rows(detail, facilities, type_labels):
    """
    The rows of a sector page's panel (UX.191): what the sector holds, the hazards in it, the facilities and where it lies.

    Args:
        detail (dict): `GET /api/sectors/<id>`.
        facilities (dict or list): `GET /api/sectors/<id>/facilities`.
        type_labels (dict): A phenomenon type to the words a page uses for it (plural not needed: it is counted).
    """
    rows = []
    systems = detail.get("systems") or []
    stars = detail.get("star_count")
    rows.append({"label": "Star systems", "html": esc(f"{len(systems):,}" + (f", {stars:,} stars" if stars else ""))})
    home = [item for item in detail.get("phenomena") or [] if item.get("home", True)]
    by_type = {}
    for item in home:
        by_type[item["type"]] = by_type.get(item["type"], 0) + 1
    for kind in sorted(by_type, key=lambda k: type_labels.get(k, k)):
        rows.append({"label": type_labels.get(kind, kind), "html": esc(f"{by_type[kind]:,}")})
    hazards = [item for item in home if item["type"] in HAZARD_TYPES]
    if hazards:
        nearest = min(hazards, key=lambda item: item.get("distance_ly") if item.get("distance_ly") is not None else 1e18)
        rows.append({"label": "Hazards", "html": esc(f"{len(hazards):,} compact object" + ("" if len(hazards) == 1 else "s")
                                                      + f", nearest {nearest['name']}")})
    held = facilities.get("items") if isinstance(facilities, dict) else facilities
    if held:
        rows.append({"label": "Facilities", "html": esc(f"{len(held):,}")})
    bubble = detail.get("superbubble")
    if bubble:
        rows.append({"label": "Inside", "html": esc(f"a superbubble {bubble['radius_pc'] * 2:g} pc across")})
    region = item_pages.region_row(detail)
    if region:
        rows.append(region)
    return [{"label": row["label"], "html": trusted_html(row["html"])} for row in rows]


def keep_out_row(found):
    """
    The keep-out row of a page whose object has one (UX.192): the radius a ship must not enter, or the note that it can
    pass through. `found` is `GET /api/objects/<ref>`; `None` for an object with neither.
    """
    keep_out = found.get("keep_out") or {}
    if keep_out.get("radius_km"):
        text = inline_html([model.qty("length", keep_out["radius_km"], "km")])
        return {"label": "Keep-out radius", "html": trusted_html(text)}
    if keep_out.get("basis") == "none" and keep_out.get("note"):
        return {"label": "Keep-out", "html": trusted_html(esc(keep_out["note"]))}
    return None


def galaxy_rows(sectors, containers):
    """
    The rows of the galaxy page's panel (UX.194): how much is charted and the big structures a course passes through.

    Args:
        sectors (list[dict]): `GET /api/galaxy/sectors`.
        containers (list[dict]): `GET /api/containers`.
    """
    rows = [{"label": "Charted sectors", "html": esc(f"{len(sectors):,}")},
            {"label": "Star systems", "html": esc(f"{sum(s.get('system_count') or 0 for s in sectors):,}")}]
    for kind, label in (("globular_cluster", "Globular clusters"), ("superbubble", "Superbubbles")):
        number = sum(1 for c in containers if c["kind"] == kind)
        if number:
            rows.append({"label": label, "html": esc(f"{number:,}")})
    return [{"label": row["label"], "html": trusted_html(row["html"])} for row in rows]
