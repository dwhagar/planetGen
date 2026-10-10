"""
Population and politics read endpoints (POP.1 to POP.4, schema v44): species,
polities and territories. See docs/design/population-and-politics.md.

Kept in its own blueprint so the population work doesn't touch
`routes.py`; it shares that module's request-scoped connection and
pagination rules.
"""

from flask import Blueprint, jsonify, request

from planetgen.population import model

from .common import ApiError
from .routes import _paginate, _parse_sort, get_db

bp = Blueprint("population", __name__, url_prefix="/api")

TERRITORY_POINT_LIMIT = 20000
"""int: The most owned systems `/api/territories` returns in one call."""


def _spacefaring_filter(raw):
    """`?spacefaring=1|0|true|false`, or absent for no filter."""
    if raw is None:
        return None
    value = raw.strip().lower()
    if value in ("1", "true", "yes"):
        return True
    if value in ("0", "false", "no"):
        return False
    raise ApiError(f"spacefaring must be true or false, got {raw!r}")


@bp.route("/population")
def population_status():
    """`GET /api/population` -- `{generated, species, polities,
    territories}` booleans: which population data exists, so pages can
    hide themselves when there is none."""
    return jsonify(model.population_status(get_db()))


@bp.route("/species")
def species_list():
    """
    `GET /api/species` -- a page of species; `?spacefaring=` filters to
    (non-)spacefaring ones. The Species table (UX.41) also takes `sort`
    (`name` -- the default --, `homeworld`, `system`, `era`, `spacefaring`,
    `polity`) with `order=asc|desc`, the repeatable filter `era`, and
    `facets=1` for the option counts of its two menus (`facets`:
    `spacefaring` `yes`/`no`, `era`). `total` counts what passes the filters.
    """
    spacefaring = _spacefaring_filter(request.args.get("spacefaring"))
    sort, descending = _parse_sort(request.args, model.SPECIES_SORTS)
    eras = [v for v in request.args.getlist("era") if v]
    limit, offset = _paginate(request.args)
    db = get_db()
    body = {
        "items": model.list_species(db, spacefaring=spacefaring, limit=limit, offset=offset, sort=sort,
                                    descending=descending, eras=eras),
        "total": model.count_species(db, spacefaring=spacefaring, eras=eras),
        "limit": limit,
        "offset": offset,
    }
    if request.args.get("facets") == "1":
        body["facets"] = model.species_facets(db, spacefaring=spacefaring, eras=eras)
    return jsonify(body)


@bp.route("/species/<int:species_id>")
def species_detail(species_id):
    """`GET /api/species/<id>` -- one species."""
    found = model.species_detail(get_db(), species_id)
    if found is None:
        raise ApiError(f"no such species: {species_id}", status_code=404)
    return jsonify(found)


@bp.route("/planets/<uid:planet_id>/species")
def planet_species(planet_id):
    """`GET /api/planets/<id>/species` -- the dominant species of a life
    world; 404 when the planet has none."""
    found = model.species_on_planet(get_db(), planet_id)
    if found is None:
        raise ApiError(f"no species on planet {planet_id}", status_code=404)
    return jsonify(found)


@bp.route("/polities")
def polity_list():
    """
    `GET /api/polities` -- a page of polities. The Polities table (UX.41)
    takes `sort` (`name` -- the default --, `species`, `government`,
    `capital`, `systems`, `reach`) with `order=asc|desc`, the repeatable
    filters `government` and `era` (the species'), and `facets=1` for the
    option counts of those two menus (`facets`: `government`, `era`).
    `total` counts what passes the filters.
    """
    sort, descending = _parse_sort(request.args, model.POLITY_SORTS)
    filters = {"governments": [v for v in request.args.getlist("government") if v],
               "eras": [v for v in request.args.getlist("era") if v]}
    limit, offset = _paginate(request.args)
    db = get_db()
    body = {
        "items": model.list_polities(db, limit=limit, offset=offset, sort=sort, descending=descending, **filters),
        "total": model.count_polities(db, **filters),
        "limit": limit,
        "offset": offset,
    }
    if request.args.get("facets") == "1":
        body["facets"] = model.polity_facets(db, **filters)
    return jsonify(body)


@bp.route("/polities/<int:polity_id>")
def polity_detail(polity_id):
    """`GET /api/polities/<id>` -- one polity and a page of its systems,
    nearest the capital first (`limit`/`offset` page the systems; `sort`
    `name` or `distance` with `order=asc|desc` orders them)."""
    limit, offset = _paginate(request.args)
    sort, descending = _parse_sort(request.args, model.POLITY_SYSTEM_SORTS) if request.args.get("sort") \
        or request.args.get("order") else ("distance", False)
    found = model.polity_detail(get_db(), polity_id, limit=limit, offset=offset, sort=sort, descending=descending)
    if found is None:
        raise ApiError(f"no such polity: {polity_id}", status_code=404)
    return jsonify(found)


@bp.route("/systems/<uid:system_id>/owner")
def system_owner(system_id):
    """`GET /api/systems/<id>/owner` -- the polity that owns a system, or
    `{"owner": null}` when none does; 404 for an unknown system."""
    db = get_db()
    if db.execute("SELECT 1 FROM star_systems WHERE id = ?", (system_id,)).fetchone() is None:
        raise ApiError(f"no such system: {system_id}", status_code=404)
    return jsonify({"owner": model.system_owner(db, system_id)})


@bp.route("/territories")
def territories():
    """`GET /api/territories` -- owned systems with galaxy-frame positions
    (parsecs) and polity colors, for a territory overlay, plus every
    polity's capital and reach."""
    db = get_db()
    return jsonify({
        "points": model.territory_points(db, limit=TERRITORY_POINT_LIMIT),
        "polities": [
            {"id": polity_id, "capital_pc": list(capital) if capital else None, "reach_ly": reach}
            for polity_id, capital, reach in model.capital_positions(db)
        ],
    })
