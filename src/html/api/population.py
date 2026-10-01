"""
Population and politics read endpoints (POP.1 to POP.4, schema v44): species,
polities and territories. See docs/design/population-and-politics.md.

Kept in its own blueprint so the population work doesn't touch
`routes.py`; it shares that module's request-scoped connection and
pagination rules.
"""

from flask import Blueprint, jsonify, request

from stellarObjects import population

from .common import ApiError
from .routes import _paginate, get_db

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
    return jsonify(population.population_status(get_db()))


@bp.route("/species")
def species_list():
    """`GET /api/species` -- a page of species by name; `?spacefaring=`
    filters to (non-)spacefaring ones."""
    spacefaring = _spacefaring_filter(request.args.get("spacefaring"))
    limit, offset = _paginate(request.args)
    db = get_db()
    return jsonify({
        "items": population.list_species(db, spacefaring=spacefaring, limit=limit, offset=offset),
        "total": population.count_species(db, spacefaring=spacefaring),
        "limit": limit,
        "offset": offset,
    })


@bp.route("/species/<int:species_id>")
def species_detail(species_id):
    """`GET /api/species/<id>` -- one species."""
    found = population.species_detail(get_db(), species_id)
    if found is None:
        raise ApiError(f"no such species: {species_id}", status_code=404)
    return jsonify(found)


@bp.route("/planets/<int:planet_id>/species")
def planet_species(planet_id):
    """`GET /api/planets/<id>/species` -- the dominant species of a life
    world; 404 when the planet has none."""
    found = population.species_on_planet(get_db(), planet_id)
    if found is None:
        raise ApiError(f"no species on planet {planet_id}", status_code=404)
    return jsonify(found)


@bp.route("/polities")
def polity_list():
    """`GET /api/polities` -- a page of polities by name."""
    limit, offset = _paginate(request.args)
    db = get_db()
    return jsonify({
        "items": population.list_polities(db, limit=limit, offset=offset),
        "total": population.count_polities(db),
        "limit": limit,
        "offset": offset,
    })


@bp.route("/polities/<int:polity_id>")
def polity_detail(polity_id):
    """`GET /api/polities/<id>` -- one polity and a page of its systems,
    nearest the capital first (`limit`/`offset` page the systems)."""
    limit, offset = _paginate(request.args)
    found = population.polity_detail(get_db(), polity_id, limit=limit, offset=offset)
    if found is None:
        raise ApiError(f"no such polity: {polity_id}", status_code=404)
    return jsonify(found)


@bp.route("/systems/<int:system_id>/owner")
def system_owner(system_id):
    """`GET /api/systems/<id>/owner` -- the polity that owns a system, or
    `{"owner": null}` when none does; 404 for an unknown system."""
    db = get_db()
    if db.execute("SELECT 1 FROM star_systems WHERE id = ?", (system_id,)).fetchone() is None:
        raise ApiError(f"no such system: {system_id}", status_code=404)
    return jsonify({"owner": population.system_owner(db, system_id)})


@bp.route("/territories")
def territories():
    """`GET /api/territories` -- owned systems with galaxy-frame positions
    (parsecs) and polity colors, for a territory overlay, plus every
    polity's capital and reach."""
    db = get_db()
    return jsonify({
        "points": population.territory_points(db, limit=TERRITORY_POINT_LIMIT),
        "polities": [
            {"id": polity_id, "capital_pc": list(capital) if capital else None, "reach_ly": reach}
            for polity_id, capital, reach in population.capital_positions(db)
        ],
    })
