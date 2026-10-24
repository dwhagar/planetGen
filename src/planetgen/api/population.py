"""
Population and politics read endpoints (POP.1 to POP.4, schema v44): species,
polities and territories. See docs/design/habitability/population-and-politics.md.

Kept in its own blueprint so the population work doesn't touch
`routes.py`; it shares that module's request-scoped connection and
pagination rules.
"""

from flask import Blueprint, jsonify, request

from planetgen.db import species_fit
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


def _tech_ranges(args):
    """The `tech_min`/`tech_max` (tech level, 0 to 7) and `<domain>_min`/`<domain>_max` (a whole index, 0 to
    7) range filters of `/api/species`, as `{column: (low, high)}`."""
    ranges = {}
    for column, prefix in zip(model.TECH_RANGE_COLUMNS, ("tech",) + model.tech.DOMAINS):
        ends = []
        for end in ("min", "max"):
            raw = args.get(f"{prefix}_{end}")
            if raw in (None, ""):
                ends.append(None)
                continue
            try:
                value = float(raw)
            except ValueError:
                raise ApiError(f"{prefix}_{end} must be a number, got {raw!r}") from None
            if not 0 <= value <= model.tech.MAX_INDEX:
                raise ApiError(f"{prefix}_{end} must be from 0 to {model.tech.MAX_INDEX}, got {raw!r}")
            ends.append(value)
        if ends != [None, None]:
            ranges[column] = tuple(ends)
    return ranges


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
    Tech (POP.9): each species carries its six indices, `tech_level`,
    `tech_band`, `limiting_domain` and `lopsided`; `sort=tech` orders by
    level (those without one last), and `tech_min`/`tech_max` and
    `<domain>_min`/`<domain>_max` (`energy`, `materials`, `information`,
    `medical`, `propulsion`, `defense`) keep a range.
    """
    spacefaring = _spacefaring_filter(request.args.get("spacefaring"))
    sort, descending = _parse_sort(request.args, model.SPECIES_SORTS)
    eras = [v for v in request.args.getlist("era") if v]
    ranges = _tech_ranges(request.args)
    limit, offset = _paginate(request.args)
    db = get_db()
    body = {
        "items": model.list_species(db, spacefaring=spacefaring, limit=limit, offset=offset, sort=sort,
                                    descending=descending, eras=eras, tech_ranges=ranges),
        "total": model.count_species(db, spacefaring=spacefaring, eras=eras, tech_ranges=ranges),
        "limit": limit,
        "offset": offset,
    }
    if request.args.get("facets") == "1":
        body["facets"] = model.species_facets(db, spacefaring=spacefaring, eras=eras, tech_ranges=ranges)
    return jsonify(body)


@bp.route("/species/<int:species_id>")
def species_detail(species_id):
    """`GET /api/species/<id>` -- one species, with its environment `limits` (POP.12): every row of the human limits
    table with the species' own value, `changed` where it differs."""
    db = get_db()
    found = model.species_detail(db, species_id)
    if found is None:
        raise ApiError(f"no such species: {species_id}", status_code=404)
    found["limits"] = species_fit.limits_table(species_fit.limits_for_species(db, [species_id])[species_id])
    return jsonify(found)


def _habitable_species(kind, body_id):
    found = species_fit.species_for_body(get_db(), kind, body_id)
    if found is None:
        raise ApiError(f"no such {kind}: {body_id}", status_code=404)
    return jsonify(found)


@bp.route("/planets/<uid:planet_id>/habitable-species")
def planet_habitable_species(planet_id):
    """
    `GET /api/planets/<id>/habitable-species` -- the species that could live on this planet (POP.12): `unaided`
    (no Red habitability domain against the species' own limits) and `equipped` (could with the kit their tech level
    can make). `scored` is false for a gas giant, which has no surface.
    """
    return _habitable_species("planet", planet_id)


@bp.route("/moons/<uid:moon_id>/habitable-species")
def moon_habitable_species(moon_id):
    """`GET /api/moons/<id>/habitable-species` -- the same for a moon."""
    return _habitable_species("moon", moon_id)


@bp.route("/systems/<uid:system_id>/habitable-species")
def system_habitable_species(system_id):
    """
    `GET /api/systems/<id>/habitable-species` -- for each rocky planet and moon of a system that has any, the
    species that could live there (POP.12), as `bodies` of `{kind, id, name, unaided, equipped}`; the system page's
    "Species that could live here" panel.
    """
    found = species_fit.species_for_system(get_db(), system_id)
    if found is None:
        raise ApiError(f"no such system: {system_id}", status_code=404)
    return jsonify(found)


@bp.route("/species/<int:species_id>/habitable-worlds/systems/<uid:system_id>")
def species_habitable_in_system(species_id, system_id):
    """
    `GET /api/species/<id>/habitable-worlds/systems/<system>` -- the rocky planets and moons of one system that this
    species could live on (POP.12), as `unaided` and `equipped` lists of worlds with their verdicts.
    """
    found = species_fit.bodies_for_species(get_db(), species_id, system_id=system_id)
    if found is None:
        raise ApiError(f"no such species: {species_id}", status_code=404)
    return jsonify(found)


@bp.route("/species/<int:species_id>/habitable-worlds/sectors/<uid:sector_id>")
def species_habitable_in_sector(species_id, sector_id):
    """`GET /api/species/<id>/habitable-worlds/sectors/<sector>` -- the same over every system of one sector."""
    found = species_fit.bodies_for_species(get_db(), species_id, sector_id=sector_id)
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
