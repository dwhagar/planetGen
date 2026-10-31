# planetgen/web/nearby_page.py

"""
The What's nearby page, `/nearby` (NAV.44): pick a place, type a distance in
parsecs, list what is there. The list comes from `GET /api/near`
(`planetgen.db.near`); this page only picks the place and renders the rows.

URL scheme (all GET, so every step is bookmarkable):

    /nearby                                  choose a place: a sector, then a system
    /nearby?sector=5                         ... a system in it
    /nearby?place=system:12                  ask for a distance
    /nearby?place=system:12&distance=10      the list
    /nearby?place=1.5,-2,0.25&distance=10    a point in the galaxy frame (pc)
    ...&kinds=planet&kinds=moon&page=2       a kind filter, a page
    ...&stars=3                              only systems with that many stars (3+, binary, ...; UX.239)

A place is an object reference (`planetgen.galaxy.objectref`; a bare number
is a system) or `x,y,z` in galaxy-frame parsecs.
"""

from urllib.parse import urlencode

from flask import request, url_for
from markupsafe import Markup

from planetgen.db import near as near_search
from planetgen.physics.units import pc_to_ly
from planetgen.util import star_count
from planetgen.web.lib import apiclient
from planetgen.web.lib.fmt import format_distance_ly
from planetgen.web.lib.pagination import PAGE_SIZE, parse_page, render_pagination

from . import bp
from .helpers import crumb, db_name, page_url, render_page
from .nav_page import _picker, _sector_id, _sector_options, _system_options, endpoint

DEFAULT_DISTANCE_PC = 10
"""int: The distance the form starts with."""

KIND_LABELS = {
    "system": "Star systems", "star": "Stars", "planet": "Planets", "moon": "Moons", "belt": "Asteroid belts",
    "comet": "Comets", "facility": "Facilities", "nebula": "Nebulae", "asteroid_field": "Asteroid fields",
    "black_hole": "Black holes", "neutron_star": "Neutron stars", "supernova_remnant": "Supernova remnants",
    "rogue_planet": "Rogue planets", "interstellar_comet": "Interstellar comets", "quasar": "Quasars",
}
"""dict: Each searchable kind's plural label."""

UNBUILT_LABELS = {"bright_star": "bright stars"}
"""dict: The label of a not-yet-generated kind that is not a searchable kind (the pre-placed bright stars, NAV.60)."""

KIND_NAMES = {
    "system": "Star system", "star": "Star", "planet": "Planet", "moon": "Moon", "belt": "Asteroid belt",
    "comet": "Comet", "facility": "Facility", "nebula": "Nebula", "asteroid_field": "Asteroid field",
    "black_hole": "Black hole", "neutron_star": "Neutron star", "supernova_remnant": "Supernova remnant",
    "rogue_planet": "Rogue planet", "interstellar_comet": "Interstellar comet", "quasar": "Quasar",
}
"""dict: Each searchable kind's singular label."""


def nearby_url(place=None, distance=None, **params):
    """`/nearby` with a place (an `endpoint(...)` string or `x,y,z`) and distance."""
    query = {"place": place, "distance": distance, **params}
    query = {name: value for name, value in query.items() if value not in (None, "")}
    return url_for("web.nearby") + (f"?{urlencode(query, safe=':,')}" if query else "")


def _row_url(row):
    """The page of one result."""
    kind = row["kind"]
    if row.get("ref"):        # UX.232: every list opens an object by the reference the API gave it, one builder
        from .item_pages import item_url        # imported here: item_pages reaches this module through the nav page
        return item_url(row["ref"])
    if kind in near_search.object_ref.PHENOMENON_KINDS:
        return page_url("phenomenon", phenomenon_type=kind, phenomenon_id=row["id"])
    if row["system_id"] is not None:
        return page_url("system", system_id=row["system_id"])
    parent = row["parent"]
    if parent and parent["kind"] == "sector":
        return page_url("sector", sector_id=parent["ref"].split(":")[1])
    return None


def _parent_link(row):
    parent = row["parent"]
    if not parent:
        return None
    kind, parent_id = parent["ref"].split(":")
    if kind == "sector":
        return {"name": parent["name"], "url": page_url("sector", sector_id=parent_id)}
    if kind == "system":
        return {"name": parent["name"], "url": page_url("system", system_id=parent_id)}
    return {"name": parent["name"], "url": page_url("system", system_id=row["system_id"])}


def _place_pickers(sector_raw):
    """The sector-then-system picker."""
    if not sector_raw:
        return [_picker("Choose a sector", "sector", "Sector", _sector_options(), {},
                        empty="No sectors have been charted yet.")]
    sector = apiclient.get_sector(db_name(), _sector_id(sector_raw))
    return [_picker(f"Choose a system in {sector['name']}", "place", "System",
                    _system_options(sector["systems"]), {}, start_over=nearby_url(),
                    empty="No systems are placed in this sector yet.")]


def _kinds(args):
    wanted = [kind for kind in args.getlist("kinds") if kind in near_search.SEARCH_KINDS]
    return list(dict.fromkeys(wanted))


def _stars(args):
    """The stars filter as typed when it is one (a count, a count and a plus, or a name), else empty: a hand-edited URL is ignored, not an error (UX.239)."""
    text = (args.get("stars") or "").strip()
    try:
        star_count.parse_star_count(text)
    except ValueError:
        return ""
    return text


def _error_field(message):
    """The form field a `/near` refusal is about (UX.124): "distance", else None (the page's own alert)."""
    return "distance" if "distance" in message else None


@bp.route("/nearby")
def nearby():
    """The What's nearby page; see the module docstring for its URL scheme."""
    args = request.args
    place_raw = (args.get("place") or "").strip()
    page = {"section": "nearby", "description": "List everything within a distance of a place.",
            "breadcrumbs": [crumb("What's nearby")]}
    if not place_raw:
        return render_page("nearby.html", title="What's nearby", pickers=_place_pickers(args.get("sector")), **page)

    if "," in place_raw:
        place_param = place_raw
        question = {"point": place_raw}
    else:
        try:
            kind, entity_id = near_search.object_ref.parse_public(place_raw)
        except ValueError:
            raise apiclient.NotFoundError(f"No such place: {place_raw!r}")
        place_param = endpoint(kind, entity_id)
        question = {"from": place_param}

    kinds = _kinds(args)
    stars = _stars(args)
    distance_raw = (args.get("distance") or "").strip()
    context = {"place_param": place_param, "distance_text": distance_raw or str(DEFAULT_DISTANCE_PC),
               "stars_text": stars, "star_choices": star_count.STAR_COUNT_CHOICES,
               "max_distance": near_search.MAX_DISTANCE_PC,
               "kind_options": [{"value": kind, "label": KIND_LABELS[kind], "checked": kind in kinds}
                                for kind in near_search.SEARCH_KINDS]}
    if not distance_raw:
        title = "What's nearby"
        return render_page("nearby.html", title=title, place_name=place_param, asking=True, **context, **page)

    page_number = parse_page(args.get("page"))
    params = {**question, "distance": distance_raw, "limit": PAGE_SIZE, "offset": (page_number - 1) * PAGE_SIZE}
    if kinds:
        params["kinds"] = ",".join(kinds)
    if stars:
        params["stars"] = stars
    try:
        result = apiclient.get_near(db_name(), params)
    except apiclient.NotFoundError:
        raise
    except apiclient.ApiError as exc:
        if exc.status_code != 400:
            raise
        message = str(exc).split("): ", 1)[-1]
        return render_page("nearby.html", title="What's nearby", place_name=place_param, asking=True,
                           error=message, error_field=_error_field(message), **context, **page)

    place = result["place"]
    rows = [{**row, "label": KIND_NAMES[row["kind"]], "url": _row_url(row), "parent_link": _parent_link(row),
             "distance_text": format_distance_ly(pc_to_ly(row["distance_pc"]))} for row in result["rows"]]
    keep = [("place", place_param), ("distance", distance_raw)] + [("kinds", kind) for kind in kinds] + ([("stars", stars)] if stars else [])
    pages = Markup(render_pagination(request.path, keep, "page", page_number, result["total"],
                                     anchor="nearby-results", label="Nearby pages"))
    not_generated = [{"label": UNBUILT_LABELS.get(kind) or KIND_LABELS.get(kind, kind.replace("-", " ").replace("_", " ")).lower(), "count": count}
                     for kind, count in sorted(result["not_generated"].items())]
    counts = [{"label": KIND_LABELS[kind], "count": count} for kind, count in
              sorted(result["by_kind"].items(), key=lambda item: near_search.SEARCH_KINDS.index(item[0]))]
    return render_page(
        "nearby.html", title=f"Within {distance_raw} pc of {place['name']}", place_name=place["name"],
        place_url=_place_url(place), rows=rows, total=result["total"], counts=counts, pages=pages,
        not_generated=not_generated, not_generated_total=result["not_generated_total"],
        systems_capped=result["systems_in_range_capped"], systems_in_range=result["systems_in_range"], **context,
        **{**page, "breadcrumbs": [crumb("What's nearby", "nearby"), crumb(place["name"])]})


def _place_url(place):
    if place["ref"] is None:
        return None
    kind, object_id = place["ref"].split(":")
    if kind in near_search.object_ref.PHENOMENON_KINDS:
        return page_url("phenomenon", phenomenon_type=kind, phenomenon_id=object_id)
    if kind == "sector":
        return page_url("sector", sector_id=object_id)
    if kind == "system":
        return page_url("system", system_id=object_id)
    return None
