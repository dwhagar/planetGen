# html/web/population_pages.py

"""
The population pages (POP.1 to POP.4), read from the population API
(`api/population.py`, schema v44):

- `/species`: every species, 50 a page, `?spacefaring=1|0` to filter.
- `/species/<id>`: one species: its homeworld, body plan, era and polity.
- `/polities`: every polity, 50 a page.
- `/polities/<id>`: one polity and the systems it owns, nearest its
  capital first.

None of them exists until a population pass has made species (or, for
the polity pages, polities): before that each is a 404 and the header
has no Species section (`helpers.population_status`).
"""

from flask import abort, request

import apiclient
from fmt import format_number
from pagination import fetch_page, parse_page

from . import bp
from .helpers import crumb, db_name, page_url, pager, population_status, render_page

SPACEFARING_FILTERS = (
    # (value of ?spacefaring=, label, API filter)
    (None, "All", None),
    ("1", "Spacefaring", True),
    ("0", "Not spacefaring", False),
)
"""tuple: The Species list's filter links."""


def _label(value):
    """`"technological_civilization"` -> `"Technological civilization"`."""
    return str(value).replace("_", " ").capitalize() if value else ""


def format_years(years):
    """A civilization's age in words: `"about 340 years"`, `"about 1.20 × 10⁴
    years"`, `"about 3.4 million years"`; `""` when unknown."""
    if years is None:
        return ""
    if years >= 1e9:
        return f"about {years / 1e9:.1f} billion years"
    if years >= 1e6:
        return f"about {years / 1e6:.1f} million years"
    rounded = round(years, -2) if years >= 1000 else round(years, -1)
    return f"about {format_number(rounded)} years"


def format_ly(value):
    return f"{format_number(value, ',.1f')} ly" if value is not None else ""


def _require(kind):
    """404 unless population data of `kind` (`"species"`/`"polities"`)
    exists."""
    if not population_status().get(kind):
        abort(404)


def _system_url(system_id):
    return page_url("system", system_id=system_id) if system_id is not None else None


def _species_row(item):
    return {
        "name": item["name"],
        "url": page_url("species_page", species_id=item["id"]),
        "homeworld": item.get("homeworld_name"),
        "system": item.get("system_name"),
        "system_url": _system_url(item.get("star_system_id")),
        "era": _label(item.get("era")) or _label(item.get("life_stage")),
        "spacefaring": item.get("spacefaring"),
        "polity": item.get("polity_name"),
        "polity_url": page_url("polity_page", polity_id=item["polity_id"]) if item.get("polity_id") else None,
    }


def _parse_spacefaring(raw):
    for value, _label_text, flag in SPACEFARING_FILTERS:
        if raw == value:
            return value, flag
    return None, None


@bp.route("/species")
def species():
    """Every species, by name."""
    _require("species")
    db = db_name()
    spacefaring, flag = _parse_spacefaring(request.args.get("spacefaring"))
    envelope, page = fetch_page(
        lambda limit, offset: apiclient.get_species_list(db, spacefaring=flag, limit=limit, offset=offset),
        parse_page(request.args.get("species_page")),
    )
    filters = [{
        "label": label,
        "url": page_url("species", spacefaring=value),
        "current": value == spacefaring,
    } for value, label, _flag in SPACEFARING_FILTERS]
    status = population_status()
    return render_page(
        "species_list.html",
        title="Species",
        section="species",
        breadcrumbs=[crumb("Species")],
        description="Every species this generated galaxy has, from life worlds to spacefaring civilizations.",
        rows=[_species_row(item) for item in envelope["items"]],
        total=envelope["total"],
        filters=filters,
        pager=pager("species_page", page, envelope["total"], anchor="species", label="Species pages",
                    keep={"spacefaring": spacefaring}),
        polities_url=page_url("polities") if status.get("polities") else None,
    )


@bp.route("/species/<int:species_id>")
def species_page(species_id):
    """One species."""
    _require("species")
    try:
        item = apiclient.get_species(db_name(), species_id)
    except apiclient.NotFoundError:
        abort(404)
    row = _species_row(item)
    facts = [
        ("Homeworld", item.get("homeworld_name")),
        ("Star system", item.get("system_name")),
        ("Most advanced life", _label(item.get("life_stage"))),
        ("Light-harvesting pigment", _label(item.get("life_chemical"))),
        ("Build", _label(item.get("build"))),
        ("Climate", _label(item.get("climate"))),
        ("Size", _label(item.get("size"))),
        ("Era", _label(item.get("era"))),
        ("Civilization age", format_years(item.get("civilization_age_years"))),
        ("Spacefaring", "Yes" if item.get("spacefaring") else "No"),
    ]
    return render_page(
        "species.html",
        title=item["name"],
        section="species",
        breadcrumbs=[crumb("Species", "species"), crumb(item["name"])],
        description=f"The {item['name']}, a species of this generated galaxy.",
        species=row,
        facts=[(label, text) for label, text in facts if text],
    )


def _polity_row(item):
    return {
        "name": item["name"],
        "url": page_url("polity_page", polity_id=item["id"]),
        "color": item.get("color"),
        "government": item.get("government"),
        "species": item.get("species_name"),
        "species_url": page_url("species_page", species_id=item["species_id"]) if item.get("species_id") else None,
        "capital": item.get("capital_name"),
        "capital_url": _system_url(item.get("capital_system_id")),
        "era": _label(item.get("era")),
        "systems": item.get("system_count") or 0,
        "reach": format_ly(item.get("reach_ly")),
    }


@bp.route("/polities")
def polities():
    """Every polity, by name."""
    _require("polities")
    db = db_name()
    envelope, page = fetch_page(
        lambda limit, offset: apiclient.get_polities(db, limit=limit, offset=offset),
        parse_page(request.args.get("polities_page")),
    )
    return render_page(
        "polities.html",
        title="Polities",
        section="species",
        breadcrumbs=[crumb("Species", "species"), crumb("Polities")],
        description="Every interstellar polity in this generated galaxy.",
        rows=[_polity_row(item) for item in envelope["items"]],
        total=envelope["total"],
        pager=pager("polities_page", page, envelope["total"], anchor="polities", label="Polity pages"),
        territories=population_status().get("territories"),
    )


@bp.route("/polities/<int:polity_id>")
def polity_page(polity_id):
    """One polity and the systems it owns, nearest its capital first."""
    _require("polities")
    db = db_name()
    page = parse_page(request.args.get("systems_page"))
    try:
        envelope, page = fetch_page(
            lambda limit, offset: _polity_envelope(db, polity_id, limit, offset), page,
        )
    except apiclient.NotFoundError:
        abort(404)
    item = envelope["polity"]
    systems = [{
        "name": system["name"],
        "url": _system_url(system["id"]),
        "distance": format_ly(system.get("distance_ly")),
    } for system in envelope["items"]]
    return render_page(
        "polity.html",
        title=item["name"],
        section="species",
        breadcrumbs=[crumb("Species", "species"), crumb("Polities", "polities"), crumb(item["name"])],
        description=f"{item['name']}, a polity of this generated galaxy, and the systems it holds.",
        polity=_polity_row(item),
        systems=systems,
        total=envelope["total"],
        pager=pager("systems_page", page, envelope["total"], anchor="polity-systems", label="System pages"),
        territories=population_status().get("territories"),
    )


def _polity_envelope(db, polity_id, limit, offset):
    """`GET /api/polities/<id>` as the pager's envelope: its systems are
    the `items`, its owned-system count the `total`."""
    item = apiclient.get_polity(db, polity_id, limit=limit, offset=offset)
    return {"items": item["systems"], "total": item.get("system_count") or 0, "polity": item}
