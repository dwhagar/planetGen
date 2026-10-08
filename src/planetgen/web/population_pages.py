# planetgen/web/population_pages.py

"""
The population pages (POP.1 to POP.4), read from the population API
(`api/population.py`, schema v44):

- `/species`: every species, a data table (UX.41) that sorts and filters
  by spacefaring and era.
- `/species/<id>`: one species: its homeworld, body plan, era and polity.
- `/polities`: every polity, a data table that sorts and filters by
  government and era.
- `/polities/<id>`: one polity and the systems it owns (a data table,
  nearest its capital first).

None of them exists until a population pass has made species (or, for
the polity pages, polities): before that each is a 404 and the header
has no Species section (`helpers.population_status`).
"""

from flask import abort, request

from planetgen.web.lib import apiclient
from planetgen.web.lib.fmt import format_number
from planetgen.web.lib.datatable import Column, Facet, Result, Table

from . import bp, tables
from .helpers import crumb, db_name, page_url, population_status, render_page


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


def _link_cell(text, href):
    return {"text": text, "href": href} if href else {"text": text}


def _species_cells(item):
    row = _species_row(item)
    return [
        _link_cell(row["name"], row["url"]),
        {"text": row["homeworld"] or ""},
        _link_cell(row["system"] or "", row["system_url"]),
        {"text": row["era"] or "\u2013"},
        {"text": "Yes" if row["spacefaring"] else "No"},
        _link_cell(row["polity"], row["polity_url"]) if row["polity"] else {"text": "\u2013"},
    ]


_SPACEFARING_LABELS = {"yes": "Spacefaring", "no": "Not spacefaring"}


def _one_flag(values):
    """The API's `spacefaring` flag for `yes`/`no` menu choices: one choice filters, none or both don't."""
    chosen = {value for value in values if value in _SPACEFARING_LABELS}
    return None if len(chosen) != 1 else chosen == {"yes"}


def _options(options, labels=None):
    return [{"value": option["value"], "label": (labels or {}).get(option["value"]) or _label(option["value"]),
             "count": option["count"]} for option in options]


def _species_load(state, limit, offset, want_facets):
    envelope = apiclient.get_species_list(
        db_name(), spacefaring=_one_flag(state.filters["spacefaring"]), limit=limit, offset=offset,
        sort=state.sort, descending=state.descending, eras=state.filters["era"], facets=want_facets)
    facets = None
    if want_facets:
        facets = {"spacefaring": _options(envelope["facets"]["spacefaring"], _SPACEFARING_LABELS),
                  "era": _options(envelope["facets"]["era"])}
    return Result([_species_cells(item) for item in envelope["items"]], envelope["total"], facets)


SPECIES_TABLE = tables.register(Table(
    "species", "Species",
    [Column("name", "Name"), Column("homeworld", "Homeworld"), Column("system", "System"), Column("era", "Era"),
     Column("spacefaring", "Spacefaring"), Column("polity", "Polity")],
    _species_load, facets=[Facet("spacefaring", "Spacefaring"), Facet("era", "Era")], prefix="species_",
    noun=("species", "species"),
))


@bp.route("/species")
def species():
    """Every species, as a data table."""
    _require("species")
    status = population_status()
    view = tables.render(SPECIES_TABLE, request.path, anchor="species")
    return render_page(
        "species_list.html",
        title="Species",
        section="species",
        breadcrumbs=[crumb("Species")],
        description="Every species this generated galaxy has, from life worlds to spacefaring civilizations.",
        table=view,
        total=view["total"],
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


def _polity_cells(item):
    row = _polity_row(item)
    return [
        {**_link_cell(row["name"], row["url"]), **({"swatch": row["color"]} if row["color"] else {})},
        _link_cell(row["species"] or "", row["species_url"]),
        {"text": row["government"] or ""},
        _link_cell(row["capital"] or "", row["capital_url"]),
        {"text": format_number(row["systems"])},
        {"text": row["reach"]},
    ]


def _polities_load(state, limit, offset, want_facets):
    envelope = apiclient.get_polities(
        db_name(), limit=limit, offset=offset, sort=state.sort, descending=state.descending,
        governments=state.filters["government"], eras=state.filters["era"], facets=want_facets)
    facets = None
    if want_facets:
        facets = {"government": _options(envelope["facets"]["government"]),
                  "era": _options(envelope["facets"]["era"])}
    return Result([_polity_cells(item) for item in envelope["items"]], envelope["total"], facets)


POLITIES_TABLE = tables.register(Table(
    "polities", "Polities",
    [Column("name", "Name"), Column("species", "Species"), Column("government", "Government"),
     Column("capital", "Capital"), Column("systems", "Systems"), Column("reach", "Reach")],
    _polities_load, facets=[Facet("government", "Government"), Facet("era", "Era")], prefix="polities_",
    noun=("polity", "polities"),
))


@bp.route("/polities")
def polities():
    """Every polity, as a data table."""
    _require("polities")
    view = tables.render(POLITIES_TABLE, request.path, anchor="polities")
    return render_page(
        "polities.html",
        title="Polities",
        section="species",
        breadcrumbs=[crumb("Species", "species"), crumb("Polities")],
        description="Every interstellar polity in this generated galaxy.",
        table=view,
        total=view["total"],
        territories=population_status().get("territories"),
    )


def _polity_systems_load(state, limit, offset, want_facets):
    """The systems a polity owns (its id is the page's, or the table route's `polity`)."""
    polity_id = request.view_args.get("polity_id") or request.args.get("polity", type=int)
    try:
        item = apiclient.get_polity(db_name(), polity_id, limit=limit, offset=offset, sort=state.sort,
                                    descending=state.descending)
    except apiclient.NotFoundError:
        abort(404)
    rows = [[_link_cell(system["name"], _system_url(system["id"])),
             {"text": format_ly(system.get("distance_ly"))}] for system in item["systems"]]
    return Result(rows, item.get("system_count") or 0, None)


POLITY_SYSTEMS_TABLE = tables.register(Table(
    "polity-systems", "Systems",
    [Column("name", "System"), Column("distance", "From the capital")],
    _polity_systems_load, prefix="systems_", default_sort="distance", noun=("system", "systems"),
))


@bp.route("/polities/<int:polity_id>")
def polity_page(polity_id):
    """One polity and the systems it owns, as a data table, nearest its capital first."""
    _require("polities")
    db = db_name()
    try:
        item = apiclient.get_polity(db, polity_id, limit=1, offset=0)
    except apiclient.NotFoundError:
        abort(404)
    view = tables.render(POLITY_SYSTEMS_TABLE, request.path, anchor="polity-systems", polity=polity_id)
    return render_page(
        "polity.html",
        title=item["name"],
        section="species",
        breadcrumbs=[crumb("Species", "species"), crumb("Polities", "polities"), crumb(item["name"])],
        description=f"{item['name']}, a polity of this generated galaxy, and the systems it holds.",
        polity=_polity_row(item),
        table=view,
        territories=population_status().get("territories"),
    )
