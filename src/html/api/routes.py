# html/api/routes.py

"""
JSON endpoints over the planetGen database.

The read endpoints each open a strictly-read-only connection
(`queryDb.open_readonly` -- same `file:...?mode=ro`-equivalent guarantee
the `queryDb.py` CLI already relies on, see that module's docstring for
how "read-only" is enforced by the configured account's grants). Listing
endpoints delegate straight to `queryDb.py`'s existing `list_sectors`/
`list_systems`/`systems_within_radius` functions rather than
re-implementing the same SQL a third time; the two detail endpoints go
through `stellarObjects._db.load_sector`/`load_star_system` for the full
nested object graph, serialized via each class's own `to_dict()` (Phase 1
serialization, `stellarObjects/serialization.py`).

Listing endpoints (`/sectors`, `/systems`) are paginated -- this project's
own roadmap (docs/TODO.md, Phase 4) plans galaxy-scale generation, so an
unbounded `SELECT *` here would eventually return an unbounded response.
`_paginate` centralizes parsing/validating `limit`/`offset` so both routes
apply the same defaults, cap, and error behavior.

The write endpoints (`POST`/`PATCH`/`DELETE` on `/sectors`/`/systems`,
near the bottom of this file) do real inserts/updates/deletes now, gated
behind admin authentication (`authz.require_admin`) and run against a
separate, write-capable database account -- see `docs/api.md`'s
"Write endpoints" section and `stellarObjects/adminAuth.py`.
"""

import os
import sys

from flask import Blueprint, current_app, g, jsonify, request

# generate.py lives at the repo root, two levels above src/html/api/ (this
# file) -- src/ itself is already on sys.path (see html/wsgi.py's own
# docstring), but the repo root isn't, so it's added here specifically for
# this import. Only `generate_sector_neighborhood_route` below needs it.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
import generate  # noqa: E402

from queryDb import (
    NO_SECTOR,
    NavUnavailable,
    SEARCH_TAG_FACETS,
    count_phenomena,
    count_sectors,
    count_systems,
    MAX_TILES_PER_REQUEST,
    galaxy_content_stamp,
    galaxy_density_shape,
    galaxy_placed_phenomena,
    galaxy_placed_sectors,
    galaxy_tiles,
    galaxy_view,
    list_phenomena,
    list_sectors,
    list_systems,
    nav_between,
    open_readonly,
    phenomenon_detail as query_phenomenon_detail,
    search as run_search,
    sector_detail as query_sector_detail,
    system_detail as query_system_detail,
    systems_within_radius,
)
from stellarObjects import _db
from stellarObjects._db import MySQLConfig, list_databases, resolve_database
from stellarObjects.config import SystemConfig
from stellarObjects.systemData import StarSystem
from stellarObjects.utils import ly_to_milliparsecs
from wikiClient import WikiClient, WikiClientAuthError, WikiClientPageExistsError, WikiClientRequestError

from .authz import audit, require_admin
from .common import ApiError, require_json_body
from .limiter import limiter

bp = Blueprint("api", __name__, url_prefix="/api")

DEFAULT_PAGE_LIMIT = 100
MAX_PAGE_LIMIT = 500

WRITE_RATE_LIMIT = "10 per minute"
"""str: Applied to every write route, on top of the app-wide default
limits (`config.Config.RATELIMIT_DEFAULT`) -- both apply together
(Flask-Limiter doesn't replace the default with a route's own `@limit`
unless told to), so a mutating endpoint is always throttled at least
this tightly regardless of how the global default is configured. Tighter
than any read endpoint's effective limit is today, on the general
principle that a request able to change data deserves more headroom
against abuse than one that only reads it."""


def _resolve_requested_db_config():
    """
    Resolves the connection config for this request: `current_app.config["MYSQL_CONFIG"]`
    (the process-wide default -- what every route used exclusively before
    this API learned about `?db=`) unless the request supplies a `db`
    query parameter, in which case it's validated against every schema on
    that same server whose name matches the configured prefix (see
    `stellarObjects._db.list_databases`/`resolve_database`) -- the same
    validation `html/lib/dbutil.resolve_db_name` has always applied for
    the CGI browser's own `?db=` picker, now shared by this API so a
    request can't select a schema this deployment never meant to expose.

    Returns:
        MySQLConfig: Ready to pass to `open_readonly`.

    Raises:
        ApiError: 404, if `db` is given but doesn't match a listed schema.
    """
    base_config = current_app.config["MYSQL_CONFIG"]
    requested = request.args.get("db")
    if requested is None:
        return base_config
    try:
        return resolve_database(base_config, requested)
    except ValueError as exc:
        raise ApiError(str(exc), status_code=404)


def get_db():
    """
    Returns the request-scoped read-only connection, opening one on first
    use (against `_resolve_requested_db_config`'s result -- the
    process-wide default database, or `?db=` when given). Reused for the
    lifetime of the request instead of one connection per query, then
    closed by `close_db` in the app's teardown handler.

    `queryDb.open_readonly` raises a bare `SystemExit` when the
    configured MySQL server is unreachable -- correct for its own CLI
    callers (a clean process exit with a message), but `SystemExit` is a
    `BaseException`, not an `Exception`: uncaught, it would propagate
    straight through Flask's request dispatch (and the WSGI worker
    handling it) instead of being turned into any HTTP response at all,
    from *every* route that calls this, not just `/health`'s own explicit
    check. Converted here into an `ApiError` (503 -- the same "database
    unreachable" status `/health` already wants to report) so it's just
    an ordinary exception from this point on: the app's registered
    `ApiError` handler turns it into the usual JSON error response for
    every other route, and `/health`'s own `except Exception` catches it
    directly.
    """
    if "db" not in g:
        try:
            g.db = open_readonly(_resolve_requested_db_config())
        except SystemExit as exc:
            raise ApiError(str(exc), status_code=503)
    return g.db


def close_db(exception=None):
    db = g.pop("db", None)
    if db is not None:
        db.close()


def _paginate(query_args):
    """
    Parses and validates the `limit`/`offset` query parameters shared by
    every listing endpoint.

    `limit` defaults to `DEFAULT_PAGE_LIMIT` and is capped at
    `MAX_PAGE_LIMIT` (silently clamped, not rejected -- a client asking for
    "too much" isn't an error, just more than this API will hand back in
    one response). `offset` defaults to 0. Both must be non-negative
    integers when given at all.

    Args:
        query_args (werkzeug.datastructures.MultiDict): `request.args`.

    Returns:
        tuple[int, int]: `(limit, offset)`.

    Raises:
        ApiError: If `limit`/`offset` is present but not a non-negative
            integer.
    """
    raw_limit = query_args.get("limit")
    raw_offset = query_args.get("offset")

    if raw_limit is None:
        limit = DEFAULT_PAGE_LIMIT
    else:
        try:
            limit = int(raw_limit)
        except ValueError:
            raise ApiError(f"limit must be an integer, got {raw_limit!r}")
        if limit < 1:
            raise ApiError("limit must be at least 1")
        limit = min(limit, MAX_PAGE_LIMIT)

    if raw_offset is None:
        offset = 0
    else:
        try:
            offset = int(raw_offset)
        except ValueError:
            raise ApiError(f"offset must be an integer, got {raw_offset!r}")
        if offset < 0:
            raise ApiError("offset must be at least 0")

    return limit, offset


def _parse_size_range(query_args, prefix):
    """
    Parses a `<prefix>_min_radius_km`/`<prefix>_max_radius_km` query
    parameter pair into a `(min_km, max_km)` tuple for `queryDb.search`'s
    `sizes` argument -- either bound may be omitted (`None`) for "no
    lower/upper limit"; giving neither returns `None` altogether (no size
    filter for this entity at all), the same "absent means inactive"
    convention every other search filter here already uses.

    Args:
        query_args (werkzeug.datastructures.MultiDict): `request.args`.
        prefix (str): `"star"`, `"planet"`, or `"moon"`.

    Returns:
        tuple[float or None, float or None] or None.

    Raises:
        ApiError: If either bound is present but not a valid non-negative
            number, or if both are given with min > max.
    """
    def _parse_bound(name):
        raw = query_args.get(name)
        if raw is None or raw.strip() == "":
            return None
        try:
            value = float(raw)
        except ValueError:
            raise ApiError(f"{name} must be a number, got {raw!r}")
        if value < 0:
            raise ApiError(f"{name} must be at least 0")
        return value

    min_km = _parse_bound(f"{prefix}_min_radius_km")
    max_km = _parse_bound(f"{prefix}_max_radius_km")
    if min_km is None and max_km is None:
        return None
    if min_km is not None and max_km is not None and min_km > max_km:
        raise ApiError(f"{prefix}_min_radius_km must not exceed {prefix}_max_radius_km")
    return (min_km, max_km)


@bp.route("/health")
@limiter.exempt
def health():
    """
    Liveness/readiness check for monitoring -- confirms the process is up
    and the configured database can actually be opened, not just that
    Flask is responding. Returns 503 (rather than letting the connection
    error propagate to a generic 500) so a monitor can tell "the API is
    running but its database is unreachable" apart from "the API itself is
    down".

    Also reports `schema_version`/`schema_current`: this API's read-only
    connection (`open_readonly`'s `ensure_schema=False`) never applies
    schema DDL itself, and even a read-write connection's `_ensure_schema`
    only runs `CREATE TABLE IF NOT EXISTS` -- a no-op against a table that
    already exists, so it never retroactively adds a migration's `ALTER
    TABLE` (e.g. schema v25's `idx_sectors_center`) to an existing
    database. Only `migrateDb.py` (run directly, or via `update.sh`/
    `install.sh`) actually advances an existing database's schema.
    Restarting this process alone -- a natural thing to try after pulling
    in a schema-fixing code change -- does *not* apply a pending
    migration, and previously wasn't surfaced anywhere: a stale schema
    silently kept e.g. the pre-v25 full-table-scan behind `GET
    /api/galaxy/view` that took the whole site down (`docs/apache-
    deployment.md`'s single `planetgen-api` process/thread pool serializes
    all API traffic, so one slow endpoint stalls every page). Surfaced
    here instead of only in `migrateDb.py`'s own output, so a live
    deployment that's fallen behind is visible without having to
    separately remember to go check.
    """
    try:
        get_db().execute("SELECT 1")
    except Exception as exc:
        return jsonify({"status": "error", "detail": str(exc)}), 503

    # A separate try/except from the reachability check above: a database
    # that answers `SELECT 1` fine but has never had `schema.sql`/
    # `migrateDb.py` applied to it at all has no `schema_migrations` table
    # yet either -- one step further back than "some migrations pending"
    # (schema_row would simply come back empty for that), not an
    # unreachable database. Reported the same way a stale-but-present
    # schema_migrations row is below (200, schema_current False), not
    # folded into the 503 case above.
    try:
        schema_row = get_db().execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()
        schema_version = schema_row["version"] if schema_row else None
    except Exception:
        schema_version = None

    body = {
        "status": "ok",
        "schema_version": schema_version,
        "schema_current": schema_version == _db.SCHEMA_VERSION,
    }
    if schema_version != _db.SCHEMA_VERSION:
        body["detail"] = (
            f"Database schema is at v{schema_version}, code expects v{_db.SCHEMA_VERSION} -- "
            f"run migrateDb.py (or update.sh/install.sh) against this database."
            if schema_version is not None else
            f"Database schema has not been initialized yet (no schema_migrations table), code expects "
            f"v{_db.SCHEMA_VERSION} -- run migrateDb.py (or update.sh/install.sh) against this database."
        )
    return jsonify(body)


@bp.route("/databases")
def databases():
    """
    Lists every MySQL schema on the configured server whose name matches
    this deployment's prefix (`stellarObjects._db.list_databases`), each
    with its size/last-modified stats plus a quick-glance sector/system
    count -- the data `html/index.py`'s database picker needs, and (via
    the sidenav's "Databases" link, shown on every page) whether that
    picker has anything to offer at all. Every other endpoint's own
    `?db=` selects among these same names (see `get_db`).
    """
    base_config = current_app.config["MYSQL_CONFIG"]
    entries = list_databases(base_config)
    items = []
    for entry in entries:
        conn = open_readonly(MySQLConfig(
            host=base_config.host, port=base_config.port,
            user=base_config.user, password=base_config.password, database=entry["name"],
        ))
        try:
            sector_count = conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"]
            system_count = conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"]
        except Exception:
            # A schema matching the configured prefix but missing this
            # project's own tables (e.g. mid-migration, or a stray
            # unrelated database sharing the prefix) shouldn't take down
            # the whole listing -- report it with unknown counts instead.
            sector_count = system_count = None
        finally:
            conn.close()
        items.append({**entry, "sector_count": sector_count, "system_count": system_count})
    return jsonify({"items": items})


@bp.route("/sectors")
def sectors():
    limit, offset = _paginate(request.args)
    db = get_db()
    return jsonify({
        "items": list_sectors(db, limit=limit, offset=offset),
        "total": count_sectors(db),
        "limit": limit,
        "offset": offset,
    })


@bp.route("/sectors/<int:sector_id>")
def sector_detail(sector_id):
    try:
        detail = query_sector_detail(get_db(), sector_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(detail)


def _parse_sector_id_filter(raw_sector_id):
    """
    Parses the `sector_id` query parameter `/api/systems` accepts: an
    integer (one sector), the literal `"none"` (standalone systems --
    `queryDb.NO_SECTOR`, `html/browse.py`'s own table of systems generated
    with no sector), or absent (no filter).

    Raises:
        ApiError: If present and neither `"none"` nor a valid integer.
    """
    if raw_sector_id is None:
        return None
    if raw_sector_id.strip().lower() == "none":
        return NO_SECTOR
    try:
        return int(raw_sector_id)
    except ValueError:
        raise ApiError(f"sector_id must be an integer or 'none', got {raw_sector_id!r}")


@bp.route("/systems")
def systems():
    star_type = request.args.get("star_type")
    sector_id = _parse_sector_id_filter(request.args.get("sector_id"))

    limit, offset = _paginate(request.args)
    db = get_db()
    rows = list_systems(db, star_type_prefix=star_type, sector_id=sector_id, limit=limit, offset=offset)
    return jsonify({
        "items": rows,
        "total": count_systems(db, star_type_prefix=star_type, sector_id=sector_id),
        "limit": limit,
        "offset": offset,
    })


@bp.route("/systems/<int:system_id>")
def system_detail(system_id):
    try:
        detail = query_system_detail(get_db(), system_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(detail)


@bp.route("/systems/<int:system_id>/near")
def systems_near(system_id):
    # This endpoint itself still returns bare JSON, just a list of ids and
    # distances -- but "no way to actually see the two points" is covered
    # now by /api/nav's web page (`html/nav.py`), which renders the
    # `navmap.py` top-down plot on top of the same course/route data.
    raw_radius = request.args.get("radius")
    if raw_radius is None:
        raise ApiError("radius query parameter is required")
    try:
        radius = float(raw_radius)
    except ValueError:
        raise ApiError(f"radius must be a number, got {raw_radius!r}")
    if radius <= 0:
        raise ApiError("radius must be greater than 0")

    try:
        matches = systems_within_radius(get_db(), system_id, radius)
    except SystemExit as exc:
        # systems_within_radius is shared with the queryDb.py CLI and raises
        # SystemExit (its CLI-appropriate error signal) for a missing/
        # unplaced system id -- caught here rather than changing its shared
        # behavior just for this one caller.
        return jsonify({"error": str(exc)}), 404
    return jsonify(matches)


def _parse_system_id_param(query_args, name):
    """
    Parses a required `<name>` query parameter as an id -- a
    `star_systems.id` when its matching `<name>_kind` is `"system"`
    (the default), or a phenomenon table's own row id when it's
    `"phenomenon"` (see `_parse_nav_endpoint_params`).

    Args:
        query_args (werkzeug.datastructures.MultiDict): `request.args`.
        name (str): The query parameter's name (`"from"` or `"to"`).

    Returns:
        int: The parsed id.

    Raises:
        ApiError: If the parameter is missing or not an integer.
    """
    raw_value = query_args.get(name)
    if raw_value is None:
        raise ApiError(f"{name} query parameter is required")
    try:
        return int(raw_value)
    except ValueError:
        raise ApiError(f"{name} must be an integer, got {raw_value!r}")


def _parse_nav_endpoint_params(query_args, name):
    """
    Parses one `/api/nav` endpoint's full reference: its id
    (`_parse_system_id_param`) plus `<name>_kind` (`"system"`, the
    default, or `"phenomenon"`) and, when it's a phenomenon,
    `<name>_type` (one of `queryDb._PHENOMENON_TYPE_TO_TABLE`'s keys --
    validated by `nav_between`/`_load_nav_phenomenon_endpoint` itself,
    via the same `ValueError` -> 404 handling `phenomenon()` above
    already uses for the same set of types, not re-validated here).

    Args:
        query_args (werkzeug.datastructures.MultiDict): `request.args`.
        name (str): `"from"` or `"to"`.

    Returns:
        tuple: `(id, kind, phenomenon_type)` -- `phenomenon_type` is
              `None` when `kind == "system"`.

    Raises:
        ApiError: If the id is missing/not an integer, `<name>_kind` is
            neither `"system"` nor `"phenomenon"`, or `kind ==
            "phenomenon"` with no `<name>_type` given.
    """
    endpoint_id = _parse_system_id_param(query_args, name)
    kind = query_args.get(f"{name}_kind", "system")
    if kind not in ("system", "phenomenon"):
        raise ApiError(f"{name}_kind must be 'system' or 'phenomenon', got {kind!r}")
    phenomenon_type = None
    if kind == "phenomenon":
        phenomenon_type = query_args.get(f"{name}_type")
        if not phenomenon_type:
            raise ApiError(f"{name}_type query parameter is required when {name}_kind is 'phenomenon'")
    return endpoint_id, kind, phenomenon_type


@bp.route("/nav")
def nav():
    """
    Course, distance, and optimal route between two endpoints -- each
    either a star system (the default) or a standalone phenomenon
    (`?from_kind=phenomenon&from_type=nebula&from=<id>`, and likewise for
    `to`) -- see `queryDb.nav_between` for the full availability rules
    and docs/api.md for the response shape.
    """
    from_id, from_kind, from_type = _parse_nav_endpoint_params(request.args, "from")
    to_id, to_kind, to_type = _parse_nav_endpoint_params(request.args, "to")

    try:
        result = nav_between(
            get_db(), from_id, to_id,
            from_kind=from_kind, to_kind=to_kind, from_type=from_type, to_type=to_type,
        )
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    except NavUnavailable as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify({
        "scope": result["scope"],
        "direct": result["direct"]._asdict(),
        "warp_times": [leg._asdict() for leg in result["warp_times"]],
        "origin_position": result["origin_position"],
        "destination_position": result["destination_position"],
        "route": _route_for_json(result["route"]),
    })


def _route_for_json(route):
    """
    `route["positions"]` can have both plain-int keys (a system hop) and
    `queryDb._phenomenon_nav_key`-shaped string keys (a phenomenon
    origin/destination) since phenomenon endpoints were added -- Flask's
    `jsonify` sorts dict keys by default, and comparing an int key against
    a string key mid-sort raises `TypeError: '<' not supported between
    instances of 'int' and 'str'` (confirmed directly: a system-to-
    phenomenon `/api/nav` request 500'd here before this fix). JSON object
    keys are always strings regardless, so this stringifies every key at
    this JSON-serialization boundary specifically, rather than changing
    `nav_between`'s own return shape -- its direct Python callers (e.g.
    `test_navigation.py`) still index a system hop's own position by its
    real int id.
    """
    if route is None:
        return None
    return {
        "path": route["path"],
        "distance_ly": route["distance_ly"],
        "positions": {str(node_id): position for node_id, position in route["positions"].items()},
    }


@bp.route("/galaxy/sectors")
def galaxy_sectors():
    """
    Every galaxy-placed sector (`sectors.center_x/y/z_pc` not NULL), with
    its live system count -- the data `html/galaxy.py`'s Galaxy Map plots.
    Not paginated: bounded by how much of the galaxy has actually been
    generated so far (see `docs/TODO.md`'s Phase 4 lazy-generation design),
    not by the addressable galaxy's own astronomical scale.
    """
    return jsonify({"items": galaxy_placed_sectors(get_db())})


@bp.route("/galaxy/phenomena")
def galaxy_phenomena():
    """
    Every galaxy-placed standalone phenomenon (`center_x/y/z_pc` not
    NULL, in any `queryDb._PHENOMENON_TABLES` table -- see `schema.sql`'s
    "v18"/"v21"/"v28" header notes) -- the phenomenon counterpart to `/api/galaxy/sectors`, plotted
    as small dots on the same `html/galaxy.py` Galaxy Map. Not paginated,
    for the same reason `/api/galaxy/sectors` isn't.
    """
    return jsonify({"items": galaxy_placed_phenomena(get_db())})


@bp.route("/galaxy/shape")
def galaxy_shape():
    """
    The galaxy's stored density-skeleton shape (`generate.py plan`'s own
    output, `queryDb.galaxy_density_shape`) -- the real spiral/disk/bulge
    model the Galaxy Map (`html/galaxy.py`) shades its "expected density"
    cloud from, for whatever space hasn't actually been generated yet.
    `"shape"` is `null` when the skeleton has never been built (the map
    then falls back to its own generic illustrative gradient).
    """
    return jsonify({"shape": galaxy_density_shape(get_db())})


MAX_GALAXY_VIEW_RADIUS_PC = 20000.0
"""float: Silently clamps an oversized `radius_pc` on `/galaxy/view` --
covers this project's own default Milky-Way-scale galaxy radius
(`program_constants.GALAXY_RADIUS_PC`, 15,000 pc) with headroom, while
still bounding how large a bounding-box scan `queryDb.
galaxy_sectors_in_view` ever has to run for one request. `planned`/
`density` are separately capped inside `queryDb.galaxy_view` itself
(`galaxyViewport.PLANNED_RADIUS_CAP_PC`/`DENSITY_SAMPLE_COUNT`) regardless
of this clamp."""


@bp.route("/galaxy/view")
def galaxy_view_route():
    """
    The interactive 3D Galaxy Map's live viewport query -- everywhere the
    flat overview map's `/galaxy/sectors` returns the whole galaxy's own
    placed sectors in one shot, this instead returns just what's near
    `cx`/`cy`/`cz` (galaxy-frame parsecs) within `radius_pc`, across all
    three content tiers `queryDb.galaxy_view` combines (placed/planned/
    density -- see that function's own docstring). Called repeatedly
    (debounced) as the 3D map's camera moves, via `html/galaxy_view.py`
    (the browser-facing CGI proxy for this route -- the browser itself
    never calls this API directly, same as every other page in `html/`).
    """
    try:
        cx = float(request.args["cx"])
        cy = float(request.args["cy"])
        cz = float(request.args["cz"])
        radius_pc = float(request.args["radius_pc"])
    except KeyError as exc:
        raise ApiError(f"{exc.args[0]} query parameter is required")
    except ValueError:
        raise ApiError("cx/cy/cz/radius_pc must all be numbers")
    if radius_pc <= 0:
        raise ApiError("radius_pc must be greater than 0")
    radius_pc = min(radius_pc, MAX_GALAXY_VIEW_RADIUS_PC)

    return jsonify(galaxy_view(get_db(), cx, cy, cz, radius_pc))



@bp.route("/galaxy/tiles")
def galaxy_tiles_route():
    """
    The 3D Galaxy Map's cube tiles -- `tiles` is a comma-separated list of
    `level/ix/iy/iz` keys (at most `MAX_TILES_PER_REQUEST`), `density` an
    optional single key to anchor a density cloud on. See
    `queryDb.galaxy_tiles` and `stellarObjects.galaxyViewport`'s "Cube
    tiles" section. Each tile's work is bounded, so unlike `/galaxy/view`
    no request can scan an unbounded region. Called by `html/
    galaxy_tiles.py`, which caches every tile on disk and only forwards
    the ones it doesn't already have.
    """
    tile_keys = [key for key in (request.args.get("tiles") or "").split(",") if key]
    density_key = request.args.get("density") or None
    if len(tile_keys) > MAX_TILES_PER_REQUEST:
        raise ApiError(f"at most {MAX_TILES_PER_REQUEST} tiles per request")
    try:
        return jsonify(galaxy_tiles(get_db(), tile_keys, density_key))
    except ValueError as exc:
        raise ApiError(str(exc))


@bp.route("/galaxy/stamp")
def galaxy_stamp_route():
    """
    `{"stamp": "<16 hex>"}` -- changes whenever the galaxy's tile contents
    could change (see `queryDb.galaxy_content_stamp`). Tile caches key on
    it, so a cached tile is never reused after new sectors are generated.
    """
    return jsonify({"stamp": galaxy_content_stamp(get_db())})

@bp.route("/phenomena")
def phenomena():
    """
    Every exotic phenomenon (nebula/asteroid field/black hole/neutron
    star/supernova remnant/rogue planet/interstellar comet), across every
    sector and regardless of galaxy placement -- the flat, paginated
    counterpart to `/api/galaxy/phenomena` (which only returns the
    galaxy-placed subset, for the Galaxy Map). `html/phenomena.py`'s own
    listing page.
    """
    limit, offset = _paginate(request.args)
    db = get_db()
    return jsonify({
        "items": list_phenomena(db, limit=limit, offset=offset),
        "total": count_phenomena(db),
        "limit": limit,
        "offset": offset,
    })


@bp.route("/phenomena/<phenomenon_type>/<int:phenomenon_id>")
def phenomenon(phenomenon_type, phenomenon_id):
    """
    One phenomenon's full detail -- `html/phenomenon.py`'s info page.
    `phenomenon_type` is one of `queryDb._PHENOMENON_TYPE_TO_TABLE`'s keys
    (`nebula`/`asteroid_field`/`black_hole`/`neutron_star`/
    `supernova_remnant`/`rogue_planet`/`interstellar_comet`); anything
    else, or an id that doesn't exist under it, is a 404.
    """
    try:
        detail = query_phenomenon_detail(get_db(), phenomenon_type, phenomenon_id)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404
    return jsonify(detail)


@bp.route("/search")
def search():
    """
    Faceted search over the chosen database: click-to-filter attribute
    tags (object type; star spectral/luminosity class; planet/moon class,
    body type, supported life chemistry; asteroid belt density), a
    min/max size range per entity, plus a per-entity name search -- see
    `queryDb.search`'s own docstring for the full request/response shape.
    `html/search.py` is a thin renderer over this endpoint's JSON.

    Query parameters: `sector_q`/`system_q`/`star_q`/`planet_q`/`moon_q`
    (name search terms); every name in `queryDb.SEARCH_TAG_FACETS`
    (repeatable, e.g. `?class=M&class=K` for two active Class tags); and
    `star_min_radius_km`/`star_max_radius_km` (likewise `planet_`/
    `moon_`) for a size range, in km -- either bound may be omitted.
    """
    args = request.args
    texts = {
        key: (args.get(key) or "").strip()
        for key in ("sector_q", "system_q", "star_q", "planet_q", "moon_q")
    }
    tags = {facet: {v.strip() for v in args.getlist(facet) if v.strip()} for facet in SEARCH_TAG_FACETS}
    tags["type"] &= {"star", "planet", "moon", "belt"}
    tags["body"] &= {"t", "g"}
    tags["moon_body"] &= {"t", "g"}
    sizes = {
        "star": _parse_size_range(args, "star"),
        "planet": _parse_size_range(args, "planet"),
        "moon": _parse_size_range(args, "moon"),
    }

    return jsonify(run_search(get_db(), texts, tags, sizes=sizes))


@bp.route("/wiki-config")
def wiki_config():
    """`GET /api/wiki-config` -- `{"wikijs": bool, "mediawiki": bool}`,
    whether each backend has a usable `base_url` plus credentials
    configured (see `config.py`'s `_wiki_config`) and so is offered as an
    "Upload to Wiki" target at all. Read by the CGI browser
    (`html/system.py`/`html/sector.py`) to decide which backend choice(s)
    to show its upload form -- never leaks any of `WIKI_CONFIG`'s actual
    credential values, only these two booleans."""
    wiki = current_app.config["WIKI_CONFIG"]
    return jsonify({"wikijs": wiki["wikijs"]["configured"], "mediawiki": wiki["mediawiki"]["configured"]})


# ---------------------------------------------------------------------
# Write endpoints. Every one requires an authenticated admin
# (`authz.require_admin`) whose credentials aren't still the seeded
# default (`fresh=True` -- see `authz.py`/`adminAuth.py`), runs against
# `WRITE_MYSQL_CONFIG` (a separate, less-privileged-than-`CREATE`
# account -- see `config.py`) rather than the request-scoped read-only
# `get_db()` every read route above uses, and records one
# `admin_audit_log` row (`authz.audit`) after it actually succeeds, never
# speculatively. See docs/api.md's "Write endpoints" section for the
# request/response shapes.
# ---------------------------------------------------------------------

SECTOR_FIELDS = {
    # field name -> (expected Python type(s), validator)
    "name": (str, lambda v: bool(v.strip())),
    "edge_ly": ((int, float), lambda v: v > 0),
}
"""dict: The `sectors` JSON object shape both `create_sector` and
`update_sector` validate against -- see docs/api.md's "Sectors" write
schema. Shared so the two can never validate the same field two
different ways."""

SECTOR_UPDATE_FIELDS = {
    **SECTOR_FIELDS,
    # `None` clears a manually-set/uploaded link back to "no page yet" --
    # see schema.sql's "v22" header note. `create_sector` deliberately
    # doesn't accept this (a brand-new, just-generated sector has never
    # been uploaded anywhere), so it's added only to `update_sector`'s own
    # allowed-fields set, not to SECTOR_FIELDS itself.
    "wiki_url": ((str, type(None)), lambda v: v is None or bool(v.strip())),
}
"""dict: `SECTOR_FIELDS` plus `update_sector`-only fields -- see
`_validate_sector_fields`'s `allowed` parameter."""


def _validate_sector_fields(body, required, allowed=SECTOR_FIELDS):
    """
    Checks `body` against `allowed`: every key in `required` must be
    present, and every key actually present (required or not) must match
    its expected type and pass its validator.

    Args:
        body (dict): The parsed request body.
        required (set[str]): Field names that must be present.
        allowed (dict): The field shape to validate against -- `SECTOR_FIELDS`
            for `create_sector`, `SECTOR_UPDATE_FIELDS` for `update_sector`
            (see that dict's own docstring for why they differ).

    Raises:
        ApiError: On a missing required field, an unrecognized field, a
            wrong-typed field, or one that fails its validator.
    """
    unknown = set(body) - set(allowed)
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")

    missing = required - set(body)
    if missing:
        raise ApiError(f"missing required field(s): {', '.join(sorted(missing))}")

    for field, (expected_type, is_valid) in allowed.items():
        if field not in body:
            continue
        value = body[field]
        if not isinstance(value, expected_type) or isinstance(value, bool) or not is_valid(value):
            raise ApiError(f"'{field}' is invalid: {value!r}")


def _resolve_requested_write_db_config():
    """
    Same as `_resolve_requested_db_config` above, but resolved against
    `WRITE_MYSQL_CONFIG` (see `config.py`) instead of the read-only
    `MYSQL_CONFIG` every read route uses -- every write route's `?db=`
    selection (or lack of one, falling back to the write account's own
    configured default database) needs the write-capable account's
    credentials, not the `SELECT`-only one.
    """
    base_config = current_app.config["WRITE_MYSQL_CONFIG"]
    requested = request.args.get("db")
    if requested is None:
        return base_config
    try:
        return resolve_database(base_config, requested)
    except ValueError as exc:
        raise ApiError(str(exc), status_code=404)


def _write_conn():
    """Opens a write-capable connection against the content database
    `?db=` (or the configured default) selects -- `ensure_schema=False`,
    same reasoning as `queryDb.open_readonly`/`_db.open_write`: the
    write-capable account has no `CREATE` grant (see `config.py`), so the
    schema must already exist."""
    return _db.open_write(_resolve_requested_write_db_config())


@bp.route("/sectors", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def create_sector():
    """`POST /api/sectors` `{"name": str, "edge_ly": number > 0}` (both
    required) -- inserts a new, empty `sectors` row (no systems -- use
    `sectorGen.py`/`galaxyGen.py` to generate a populated one) and
    returns its id."""
    body = require_json_body()
    _validate_sector_fields(body, required={"name", "edge_ly"})

    conn = _write_conn()
    try:
        with conn:
            cur = conn.execute(
                "INSERT INTO sectors (name, edge_mpc) VALUES (?, ?)",
                (body["name"], ly_to_milliparsecs(body["edge_ly"])),
            )
            sector_id = cur.lastrowid
    finally:
        conn.close()

    audit("sector.create", target=f"sector:{sector_id}", detail=f"name={body['name']!r} edge_ly={body['edge_ly']}")
    return jsonify({"id": sector_id, "name": body["name"], "edge_ly": body["edge_ly"]}), 201


@bp.route("/sectors/<int:sector_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def update_sector(sector_id):
    """`PATCH /api/sectors/<id>` -- any non-empty subset of `{"name": str,
    "edge_ly": number > 0, "wiki_url": str or null}`. `wiki_url` is the
    admin "manually set the wiki link" affordance (`html/admin.py`) --
    `null` clears it back to "no page yet" (see `schema.sql`'s "v22"
    header note); the same column is also written automatically by
    `POST /api/sectors/<id>/wiki` on a successful upload."""
    body = require_json_body()
    if not body:
        raise ApiError("body must include at least one field to update")
    _validate_sector_fields(body, required=set(), allowed=SECTOR_UPDATE_FIELDS)

    set_clauses, params = [], []
    if "name" in body:
        set_clauses.append("name = ?")
        params.append(body["name"])
    if "edge_ly" in body:
        set_clauses.append("edge_mpc = ?")
        params.append(ly_to_milliparsecs(body["edge_ly"]))
    if "wiki_url" in body:
        set_clauses.append("wiki_url = ?")
        params.append(body["wiki_url"])
    params.append(sector_id)

    conn = _write_conn()
    try:
        with conn:
            # Existence is checked explicitly rather than relying on the
            # UPDATE's own affected-row count: MySQL (and so pymysql)
            # reports 0 affected rows for a matched row whose new values
            # equal its current ones, which would otherwise misreport a
            # real, unchanged sector as "not found".
            if conn.execute("SELECT id FROM sectors WHERE id = ?", (sector_id,)).fetchone() is None:
                raise ApiError(f"no such sector: {sector_id}", status_code=404)
            conn.execute(f"UPDATE sectors SET {', '.join(set_clauses)} WHERE id = ?", tuple(params))
    finally:
        conn.close()

    audit("sector.update", target=f"sector:{sector_id}", detail=str(body))
    return jsonify({"status": "ok"})


@bp.route("/sectors/<int:sector_id>", methods=["DELETE"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def delete_sector(sector_id):
    """`DELETE /api/sectors/<id>` -- no request body. Its systems are
    detached, not deleted (`star_systems.sector_id` is `ON DELETE SET
    NULL`, per `schema.sql`), matching how a system can already exist
    with no sector at all."""
    conn = _write_conn()
    try:
        with conn:
            deleted = conn.execute("DELETE FROM sectors WHERE id = ?", (sector_id,)).rowcount > 0
    finally:
        conn.close()

    if not deleted:
        raise ApiError(f"no such sector: {sector_id}", status_code=404)
    audit("sector.delete", target=f"sector:{sector_id}")
    return jsonify({"status": "ok"})


@bp.route("/sectors/<int:sector_id>/generate-neighborhood", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def generate_sector_neighborhood_route(sector_id):
    """`POST /api/sectors/<id>/generate-neighborhood` -- generates every
    not-yet-generated sector within `radius_ly` (optional JSON body
    field; defaults to `program_constants.RANDOM_START_NEIGHBORHOOD_RADIUS_LY`,
    100 ly, the same "100 ly sphere" `generate.py galaxy`'s own
    random-start mode uses) of this already galaxy-placed sector -- see
    `generate.generate_sector_neighborhood`. **The default radius is
    genuinely large** (~2,000-3,000 candidate sector slots, confirmed by
    measurement -- see that function's own docstring), so this can run
    for minutes to hours, not seconds. Runs synchronously like every
    other write route regardless (there's no background job queue in
    this project to hand it off to) -- `apiclient.py`'s own caller uses a
    much longer timeout than its other calls for exactly this reason, but
    a production deployment's own reverse-proxy/gateway timeout (Apache,
    etc.) may still need raising for this one route to ever complete over
    HTTP at all."""
    body = request.get_json(silent=True) or {}
    radius_ly = body.get("radius_ly")
    if radius_ly is not None and (
        not isinstance(radius_ly, (int, float)) or isinstance(radius_ly, bool) or radius_ly <= 0
    ):
        raise ApiError(f"'radius_ly' is invalid: {radius_ly!r}")

    try:
        result = generate.generate_sector_neighborhood(
            sector_id, radius_ly=radius_ly, config=_resolve_requested_write_db_config(),
        )
    except ValueError as exc:
        raise ApiError(str(exc), status_code=404)
    except RuntimeError as exc:
        # The galaxy's density skeleton (`generate.py plan`) has never
        # been built -- generate_sector_neighborhood needs it to gate each
        # candidate slot's own generation on local stellar density.
        raise ApiError(str(exc), status_code=409)

    audit(
        "sector.generate_neighborhood", target=f"sector:{sector_id}",
        detail=f"radius_ly={radius_ly!r} generated={result['generated']}",
    )
    return jsonify(result)


SYSTEM_CONFIG_TRISTATE_FIELDS = {
    "habitable_world", "asteroid_belt", "large_star", "moons",
    "max_planets", "planets", "intelligent_life", "binary_system",
    "wide_binary",
}
"""set[str]: `SystemConfig` fields that are `bool` or `None` (`None` = let
the generator decide) -- see `config.py`'s own field docstrings."""

SYSTEM_CONFIG_ALLOWED_FIELDS = SYSTEM_CONFIG_TRISTATE_FIELDS | {"markdown", "star_type", "name", "age", "num_orbits"}
"""set[str]: Every recipe field `POST /api/systems` currently accepts --
`SystemConfig.SERIALIZABLE_FIELDS` minus `SLOTS` (a nested per-orbit
structure not validated in this pass; a request naming it is rejected as
an unrecognized field rather than silently accepted-but-ignored)."""


def _validate_system_config_body(body):
    """
    Validates a `POST /api/systems` recipe body against
    `SYSTEM_CONFIG_ALLOWED_FIELDS` -- deliberately stricter than
    `SystemConfig.from_dict`/`fields_from_dict` itself (which silently
    ignores an unrecognized key rather than rejecting it), since a typo
    in a request body should be a `400`, not a silently-ignored no-op.

    Raises:
        ApiError: On an unrecognized field or one with an invalid type/
            value.
    """
    unknown = set(body) - SYSTEM_CONFIG_ALLOWED_FIELDS
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")

    if "markdown" in body and not isinstance(body["markdown"], bool):
        raise ApiError("'markdown' must be a boolean")
    for field in SYSTEM_CONFIG_TRISTATE_FIELDS:
        if field in body and body[field] is not None and not isinstance(body[field], bool):
            raise ApiError(f"'{field}' must be a boolean or null")
    for field in ("star_type", "name"):
        if field in body and body[field] is not None and not isinstance(body[field], str):
            raise ApiError(f"'{field}' must be a string or null")
    if "age" in body and body["age"] not in (None, "young", "old"):
        raise ApiError("'age' must be 'young', 'old', or null")
    if "num_orbits" in body:
        value = body["num_orbits"]
        if value is not None and (not isinstance(value, int) or isinstance(value, bool) or value <= 0):
            raise ApiError("'num_orbits' must be a positive integer or null")


@bp.route("/systems", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def create_system():
    """
    `POST /api/systems` -- a generation "recipe" body (see
    `SYSTEM_CONFIG_ALLOWED_FIELDS`), the same shape `systemGen.py
    --system-file` already takes (`SystemConfig.from_dict`). Generates a
    new **standalone** system (no `sector_id` -- attaching a newly
    generated system to an existing sector needs that sector's own
    placement/Hill-sphere separation logic, tracked as follow-up work in
    docs/TODO.md, not implemented here) and persists it via the same
    `insert_star_system` every other entry point uses. Returns the new
    `star_systems.id`.
    """
    body = require_json_body()
    _validate_system_config_body(body)
    system_config = SystemConfig.from_dict(body)

    try:
        star_system = StarSystem(system_config=system_config)
    except Exception as exc:
        # Generation can reject an internally-inconsistent recipe (e.g. an
        # impossible num_orbits/planet-class combination) -- surfaced as a
        # 400 (the request body was the problem), not an unhandled 500.
        raise ApiError(f"system generation failed: {exc}")

    conn = _write_conn()
    try:
        with conn:
            system_id = _db.insert_star_system(conn, star_system, system_config)
    finally:
        conn.close()

    audit("system.create", target=f"system:{system_id}", detail=f"star_type={body.get('star_type')!r}")
    return jsonify({"id": system_id}), 201


SYSTEM_UPDATE_FIELDS = {
    "name": (str, lambda v: bool(v.strip())),
}
"""dict: `PATCH /api/systems/<id>`'s allowed fields -- metadata only
(a rename) in this pass; editing a system's generated content (stars/
planets/moons/belts) is out of scope here, same reasoning as `POST
/api/systems` accepting only a generation recipe rather than a
hand-edited object graph (see docs/api.md)."""


@bp.route("/systems/<int:system_id>", methods=["PATCH"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def update_system(system_id):
    """`PATCH /api/systems/<id>` `{"name": str}` -- renames a system."""
    body = require_json_body()
    if not body:
        raise ApiError("body must include at least one field to update")
    unknown = set(body) - set(SYSTEM_UPDATE_FIELDS)
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")
    for field, (expected_type, is_valid) in SYSTEM_UPDATE_FIELDS.items():
        if field not in body:
            continue
        value = body[field]
        if not isinstance(value, expected_type) or not is_valid(value):
            raise ApiError(f"'{field}' is invalid: {value!r}")

    conn = _write_conn()
    try:
        with conn:
            if conn.execute("SELECT id FROM star_systems WHERE id = ?", (system_id,)).fetchone() is None:
                raise ApiError(f"no such system: {system_id}", status_code=404)
            conn.execute("UPDATE star_systems SET name = ? WHERE id = ?", (body["name"], system_id))
    finally:
        conn.close()

    audit("system.update", target=f"system:{system_id}", detail=str(body))
    return jsonify({"status": "ok"})


@bp.route("/systems/<int:system_id>", methods=["DELETE"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def delete_system(system_id):
    """`DELETE /api/systems/<id>` -- no request body. Its stars/planets/
    moons/belts are deleted with it (all `ON DELETE CASCADE`, per
    `schema.sql`). Bumps the system's sector's `modified_at`, since a
    deleted row leaves no timestamp of its own behind (`schema.sql`'s
    "v27" header note)."""
    conn = _write_conn()
    try:
        with conn:
            row = conn.execute("SELECT sector_id FROM star_systems WHERE id = ?", (system_id,)).fetchone()
            deleted = conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,)).rowcount > 0
            if deleted:
                _db.touch_sector(conn, row["sector_id"])
    finally:
        conn.close()

    if not deleted:
        raise ApiError(f"no such system: {system_id}", status_code=404)
    audit("system.delete", target=f"system:{system_id}")
    return jsonify({"status": "ok"})


# ---------------------------------------------------------------------
# Wiki publishing (schema.sql's "v22" header note) -- both routes below
# share the same request shape and backend-resolution/error-mapping
# helpers, differing only in where their page content comes from
# (`star_systems.markdown_content`/`wikitext_content`, already rendered
# at generation time, vs. `_sector_wiki_content`'s on-the-fly summary --
# sectors have no persisted rendered page of their own) and which column
# the resulting URL is written back to.
# ---------------------------------------------------------------------

def _wiki_client_for(backend):
    """
    Builds a `wikiClient.WikiClient` for `backend`, from
    `current_app.config["WIKI_CONFIG"]` (see `config.py`'s `_wiki_config`).

    Args:
        backend (str): `"wikijs"` or `"mediawiki"`.

    Returns:
        wikiClient.WikiClient

    Raises:
        ApiError: 501 if `backend` has no `base_url`/credentials
            configured at all -- rather than trying, and failing, to
            reach an empty URL.
    """
    settings = current_app.config["WIKI_CONFIG"][backend]
    if not settings["configured"]:
        raise ApiError(f"wiki publishing is not configured for {backend!r}", status_code=501)

    if backend == "wikijs":
        return WikiClient(backend="wikijs", base_url=settings["base_url"], api_token=settings["api_token"])
    return WikiClient(
        backend="mediawiki", base_url=settings["base_url"],
        username=settings["username"], password=settings["password"],
    )


def _wiki_upload_request(body):
    """
    Validates a `POST .../wiki` request body.

    Args:
        body (dict): The parsed request body -- `{"backend": "wikijs" |
            "mediawiki", "path": str}`. `path` is the target page's
            path/slug for `"wikijs"` (which addresses a page separately
            from its title -- see `wikiClient/wikijs.py`) and required
            for it; `"mediawiki"` has no such separate concept (its
            title *is* its address, see `wikiClient/mediawiki.py`), so
            `path` is accepted but ignored for it.

    Returns:
        tuple[str, str or None]: `(backend, path)` -- `path` stripped, or
            `None` if not given/blank.

    Raises:
        ApiError: On an unrecognized field, a missing/invalid `backend`,
            or a missing/blank `path` when `backend == "wikijs"`.
    """
    unknown = set(body) - {"backend", "path"}
    if unknown:
        raise ApiError(f"unrecognized field(s): {', '.join(sorted(unknown))}")

    backend = body.get("backend")
    if backend not in ("wikijs", "mediawiki"):
        raise ApiError("'backend' must be 'wikijs' or 'mediawiki'")

    path = body.get("path")
    if path is not None and not isinstance(path, str):
        raise ApiError("'path' must be a string")
    path = path.strip() if path else None

    if backend == "wikijs" and not path:
        raise ApiError("'path' is required for the 'wikijs' backend")

    return backend, path


def _create_wiki_page(client, path, title, content):
    """Wraps `client.create_page`, mapping `wikiClient`'s own exception
    hierarchy onto this API's status codes -- shared by both routes below
    so neither duplicates the mapping.

    Raises:
        ApiError: 409 (a page already exists there) or 502 (auth
            rejected, or the wiki couldn't be reached/errored).
    """
    try:
        return client.create_page(path=path, title=title, content=content)
    except WikiClientPageExistsError as exc:
        raise ApiError(str(exc), status_code=409)
    except (WikiClientAuthError, WikiClientRequestError) as exc:
        raise ApiError(str(exc), status_code=502)


@bp.route("/systems/<int:system_id>/wiki", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def upload_system_wiki(system_id):
    """
    `POST /api/systems/<id>/wiki` `{"backend": "wikijs"|"mediawiki",
    "path": str}` -- publishes this system's already-generated page
    (`markdown_content` for `wikijs`, `wikitext_content` for `mediawiki`,
    see `_db.insert_star_system`) to the chosen wiki, then records the
    new page's URL on `star_systems.wikijs_url`/`mediawiki_url` (see
    `schema.sql`'s "v22" header note) so `html/system.py` can swap its
    description section for a link to it afterward.

    A real database read (the system's own name/content, via `get_db()`'s
    request-scoped read-only connection) but not a database *write* until
    after the wiki call succeeds -- unlike every other write route above,
    which mutate the local database directly, this one's real effect is
    external, so `_write_conn()` is only opened once there is a URL worth
    persisting.
    """
    body = require_json_body()
    backend, path = _wiki_upload_request(body)
    client = _wiki_client_for(backend)

    try:
        system = query_system_detail(get_db(), system_id)
    except ValueError:
        raise ApiError(f"no such system: {system_id}", status_code=404)

    content = system["markdown_content"] if backend == "wikijs" else system["wikitext_content"]
    page_path = path if backend == "wikijs" else system["name"]
    page = _create_wiki_page(client, page_path, system["name"], content)

    url_column = "wikijs_url" if backend == "wikijs" else "mediawiki_url"
    conn = _write_conn()
    try:
        with conn:
            conn.execute(f"UPDATE star_systems SET {url_column} = ? WHERE id = ?", (page.url, system_id))
    finally:
        conn.close()

    audit("system.wiki_upload", target=f"system:{system_id}", detail=f"backend={backend!r} path={page.path!r}")
    return jsonify({"id": page.id, "path": page.path, "title": page.title, "url": page.url}), 201


def _sector_wiki_content(sector):
    """
    Builds a sector-summary page's content, fresh, in both Markdown (for
    `wikijs`) and wikitext (for `mediawiki`) -- unlike a system, a sector
    has no persisted `markdown_content`/`wikitext_content` of its own
    (see `schema.sql`'s "v22" header note), so this is rendered on the
    fly from `sector`'s already-queried detail (`queryDb.sector_detail`'s
    shape) at upload time: the sector's name/size, and one table row per
    system placed in it, with the same star-type fallback
    `html/sector.py`'s own systems table already uses (a `'close'` pair's
    merged `binary_type` if it has one, else each of its own stars'
    `star_type` joined by `" / "`).

    Args:
        sector (dict): `queryDb.sector_detail`'s return shape.

    Returns:
        tuple[str, str]: `(markdown_content, wikitext_content)`.
    """
    def star_type_for(system):
        if system["is_binary"] and system.get("binary_type"):
            return system["binary_type"]
        return " / ".join(star["star_type"] for star in system["stars"]) if system["stars"] else ""

    systems = sector["systems"]

    md_rows = "\n".join(
        f'| {s["name"]} | {s["quadrant"] or ""} | {"Yes" if s["is_binary"] else "No"} | '
        f'{star_type_for(s)} | {s["location"] or ""} |'
        for s in systems
    ) or "| *(none)* | | | | |"
    markdown_content = (
        f'# {sector["name"]}\n\n'
        f'**Cube edge:** {sector["edge_ly"]:,.2f} ly  \n'
        f'**Systems:** {sector["system_count"]}\n\n'
        "## Systems\n\n"
        "| Name | Octant | Binary | Star type | Location |\n"
        "|---|---|---|---|---|\n"
        f"{md_rows}\n"
    )

    wiki_rows = "\n|-\n".join(
        f'| {s["name"]} || {s["quadrant"] or ""} || {"Yes" if s["is_binary"] else "No"} || '
        f'{star_type_for(s)} || {s["location"] or ""}'
        for s in systems
    ) or "| ''(none)'' ||  ||  ||  || "
    wikitext_content = (
        f'= {sector["name"]} =\n\n'
        f"'''Cube edge:''' {sector['edge_ly']:,.2f} ly\n\n"
        f"'''Systems:''' {sector['system_count']}\n\n"
        "== Systems ==\n\n"
        '{| class="wikitable"\n'
        "! Name !! Octant !! Binary !! Star type !! Location\n"
        "|-\n"
        f"{wiki_rows}\n"
        "|}\n"
    )
    return markdown_content, wikitext_content


@bp.route("/sectors/<int:sector_id>/wiki", methods=["POST"])
@limiter.limit(WRITE_RATE_LIMIT)
@require_admin(fresh=True)
def upload_sector_wiki(sector_id):
    """
    `POST /api/sectors/<id>/wiki` `{"backend": "wikijs"|"mediawiki",
    "path": str}` -- publishes a freshly generated sector-summary page
    (see `_sector_wiki_content`) to the chosen wiki, then records the new
    page's URL on `sectors.wiki_url` (a single column, not one per
    backend the way `star_systems` has -- see `schema.sql`'s "v22" header
    note) so `html/sector.py` can link to it afterward. The same column
    can also be set directly, without uploading anything, via
    `PATCH /api/sectors/<id>`'s own `wiki_url` field (`html/admin.py`'s
    manual-link admin section).
    """
    body = require_json_body()
    backend, path = _wiki_upload_request(body)
    client = _wiki_client_for(backend)

    try:
        sector = query_sector_detail(get_db(), sector_id)
    except ValueError:
        raise ApiError(f"no such sector: {sector_id}", status_code=404)

    markdown_content, wikitext_content = _sector_wiki_content(sector)
    content = markdown_content if backend == "wikijs" else wikitext_content
    page_path = path if backend == "wikijs" else sector["name"]
    page = _create_wiki_page(client, page_path, sector["name"], content)

    conn = _write_conn()
    try:
        with conn:
            conn.execute("UPDATE sectors SET wiki_url = ? WHERE id = ?", (page.url, sector_id))
    finally:
        conn.close()

    audit("sector.wiki_upload", target=f"sector:{sector_id}", detail=f"backend={backend!r} path={page.path!r}")
    return jsonify({"id": page.id, "path": page.path, "title": page.title, "url": page.url}), 201
