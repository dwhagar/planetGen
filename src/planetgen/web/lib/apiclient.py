# planetgen/web/lib/apiclient.py

"""
Client for the planetGen Flask API (`planetgen/api/`), used by every page in
`planetgen/web/` instead of querying the database directly (in-process there,
see `web/transport.py`; over HTTP anywhere else).

This is the concrete result of moving the interim web browser onto the
API: every page here used to open its own read-only MySQL connection
(`html/lib/dbutil.py`'s previous `open_readonly`/`resolve_db_name`) and
run its own bespoke SQL; now every page is a thin HTTP client over
`GET /api/...` (see `docs/api.md`), and the database-querying logic those
pages used to duplicate lives once, in `planetgen.db.query`, shared with the API
itself. `planetgen/web/lib/fmt.py` still holds the formatting-only helpers
(`esc`, `linkify_location`, `format_density`) that have nothing to do
with fetching data.

Stdlib only (`urllib.request`) -- same "nothing beyond a system Python 3"
deployment story this project's docs already describe for `html/`
(see `docs/html-interface.md`), now extended to reach the API process
over HTTP rather than a database socket.

The API must be reachable for this browser to work at all now -- see
`API_BASE_URL` below and `docs/deployment/apache.md` for how it's mounted
alongside `html/` in a real deployment.
"""

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request

from planetgen.util import log
from planetgen.util.appconfig import load_config

from planetgen.web.lib import pagecache

API_BASE_URL = os.environ.get("PLANETGEN_API_BASE_URL") or load_config()["api_base_url"]
"""str: Base URL of the planetGen API's `/api` mount point. Defaults to
the same host (see `examples/apache/planetgen.conf.example`) -- only used
by the HTTP transport, not by the Flask pages' in-process one --
override via the `PLANETGEN_API_BASE_URL` Apache `SetEnv` (or shell env,
for local testing), or `config.json`'s `api_base_url` (see
`docs/config.md`), when the API is deployed at a different host/port,
e.g. `http://127.0.0.1:5000/api` for `python src/html/wsgi.py`'s own dev
server."""

_TIMEOUT_SECONDS = 30
"""int: Was 15 -- confirmed too tight for GET /api/search specifically
(TimeoutError in production): that one endpoint always runs its full
facet+autocomplete query set up front regardless of whether any filter is
active (queryDb.search), which used to mean an unindexed full-table
scan/sort per query (see schema.sql's "v22" header note, which adds the
missing indexes -- the real fix). This wider margin is deliberately kept
as a second line of defense on top of that, not a replacement for it: even
an indexed query set can occasionally run long on a large, busy database,
and every other page here issues far fewer/cheaper queries per request, so
raising this shared constant costs them nothing in the common case."""


class NotFoundError(Exception):
    """Raised for a 404 response (an invalid `db`, sector/system id, etc.)
    -- callers turn this into a 404 page (`web/errors.py`)."""


class ApiError(Exception):
    """
    Raised for any other failure talking to the API: an unreachable
    process, a non-2xx/404 response, or an unparseable body.

    `status_code` (`None` for a connection-level failure -- the API never
    responded at all) lets a caller that needs to handle one specific
    status differently do so without string-matching the message -- e.g.
    `auth_me` treating 401 ("not logged in") as a normal, expected outcome
    rather than propagating it as a page-breaking error.
    """

    def __init__(self, message, status_code=None):
        super().__init__(message)
        self.status_code = status_code


# ---------------------------------------------------------------------
# Transport: how a request actually reaches the API.
#
# The default is HTTP (`_http_transport`). The Flask-served pages (`planetgen/web/`) run *inside* the API's own process, so
# they install an in-process transport with `set_transport` (see
# `web/transport.py`) that dispatches straight through the app's own
# routes -- no HTTP round trip to itself. Every typed wrapper below
# (`get_sectors`, `auth_me`, ...) is unchanged either way.
# ---------------------------------------------------------------------

class TransportResponse:
    """What a transport hands back for one request: the HTTP `status`
    code, the response body as text, every `Set-Cookie` header verbatim,
    and a `reason` phrase used as the error detail when the body has
    none."""

    __slots__ = ("status", "body", "set_cookie_headers", "reason")

    def __init__(self, status, body, set_cookie_headers=(), reason=""):
        self.status = status
        self.body = body
        self.set_cookie_headers = list(set_cookie_headers)
        self.reason = reason


class TransportUnreachable(Exception):
    """Raised by a transport when the API could not be reached at all
    (no HTTP status to report)."""


_transport = None
"""callable or None: The installed non-HTTP transport, see
`set_transport`."""


def set_transport(transport):
    """
    Installs a transport used for every API call instead of HTTP, or
    removes it again with `None`.

    Args:
        transport (callable or None): `transport(method, target, data,
            headers, timeout)`, where `target` is the path under `/api`
            with any query string (`"/sectors?db=x&limit=50"`), `data` the
            request body bytes (or `None`) and `headers` a dict. Returns a
            `TransportResponse`, or `None` to decline the call (e.g. no
            Flask request is active), in which case HTTP is used.
    """
    global _transport
    _transport = transport


_response_cache = None
"""callable or None: Returns the `pagecache.ResponseCache` to use for the
current call (or `None` for none); see `set_response_cache`."""


def set_response_cache(provider):
    """
    Lets `_request` serve public GETs from a cache (`planetgen/web/lib/pagecache.py`).

    Args:
        provider (callable or None): `provider()` -> a
            `pagecache.ResponseCache`, or `None` when no cache applies to
            this call (the Flask pages pass their app's own). `None`
            removes the hook.
    """
    global _response_cache
    _response_cache = provider


def _cache_db(params):
    """The `db` a GET's params name, or `None`."""
    if isinstance(params, dict):
        return params.get("db") or None
    for key, value in params or ():
        if key == "db":
            return value or None
    return None


def _http_transport(method, target, data, headers, timeout):
    """The default transport: a real HTTP request to `API_BASE_URL`."""
    request = urllib.request.Request(f"{API_BASE_URL}{target}", data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return TransportResponse(
                response.status, response.read().decode("utf-8"),
                response.headers.get_all("Set-Cookie") or [], response.reason,
            )
    except urllib.error.HTTPError as exc:
        try:
            body = exc.read().decode("utf-8", errors="replace")
        except Exception:
            body = ""
        return TransportResponse(exc.code, body, [], exc.reason)
    except urllib.error.URLError as exc:
        raise TransportUnreachable(exc.reason)


def _send(method, target, json_body=None, cookie_header=None, timeout=_TIMEOUT_SECONDS):
    """
    Sends one request through the installed transport (HTTP unless
    `set_transport` installed another that accepts the call) and maps
    the outcome onto this module's exceptions.

    Returns:
        TransportResponse: For any 2xx/3xx status.

    Raises:
        NotFoundError: On a 404 response.
        ApiError: On any other error status (with `.status_code` set), or
            if the API can't be reached.
    """
    data = None
    headers = {}
    if json_body is not None:
        data = json.dumps(json_body).encode("utf-8")
        headers["Content-Type"] = "application/json"
    if cookie_header:
        headers["Cookie"] = cookie_header

    start = time.perf_counter()
    response = None
    where = f"{API_BASE_URL}{target}"
    if _transport is not None:
        response = _transport(method, target, data, headers, timeout)
        if response is not None:
            where = f"(in-process) /api{target}"
    if response is None:
        try:
            response = _http_transport(method, target, data, headers, timeout)
        except TransportUnreachable as exc:
            _log_call(method, where, start, f"unreachable: {exc}")
            raise ApiError(f"Could not reach the planetGen API at {API_BASE_URL}: {exc}")

    if response.status >= 400:
        detail = _error_detail(response)
        _log_call(method, where, start, f"HTTP {response.status}: {detail}")
        if response.status == 404:
            raise NotFoundError(detail)
        raise ApiError(f"planetGen API error ({response.status}): {detail}", status_code=response.status)
    # Request bodies aren't logged: login and credential changes carry passwords.
    _log_call(method, where, start, f"HTTP {response.status}, {len(response.body)} bytes, "
                                    f"{len(response.set_cookie_headers)} Set-Cookie header(s), "
                                    f"{'a' if json_body is not None else 'no'} JSON request body")
    return response


def _parse_json(raw_body):
    try:
        return json.loads(raw_body)
    except ValueError as exc:
        raise ApiError(f"planetGen API returned an unparseable response: {exc}")


def _request(path, params=None):
    """
    Runs one `GET` against the API and returns the parsed JSON body.

    Args:
        path (str): The path under `API_BASE_URL`, e.g. `"/sectors/5"`.
        params (dict or list[tuple], optional): Query parameters --
            a plain `dict` for single-valued params, or a list of
            `(key, value)` pairs when a key repeats (e.g. search's tag
            facets). `None`/empty values are dropped from a `dict`; a
            list of pairs is used exactly as given.

    Returns:
        The parsed JSON body (usually a `dict`; a few endpoints return a bare `list`).

    Raises:
        NotFoundError: On a 404 response.
        ApiError: On any other non-2xx response, or if the API can't be
            reached/returns an unparseable body.
    """
    query = _build_query(params)
    target = f"{path}?{query}" if query else path
    cache = _response_cache() if _response_cache is not None and pagecache.is_cacheable(path) else None
    if cache is None:
        return _parse_json(_send("GET", target).body)
    db = _cache_db(params)
    body = cache.get(db, target)
    if body is None:
        generation = cache.generation
        body = _send("GET", target).body
        cache.put(db, target, body, generation)
    return _parse_json(body)


def _auth_request(method, path, json_body=None, cookie_header=None, timeout=_TIMEOUT_SECONDS):
    """
    Runs one JSON request against the API supporting any HTTP method and
    an optional request body/`Cookie` header -- the primitive every
    `auth_*` function below builds on.

    Args:
        method (str): `"POST"`, `"DELETE"`, etc.
        path (str): The path under `API_BASE_URL`, e.g. `"/auth/login"`
            (may carry its own query string).
        json_body (dict, optional): Sent as the request body
            (`Content-Type: application/json`) if given.
        cookie_header (str, optional): Forwarded as-is as the outgoing
            `Cookie` header -- callers pass the browser's own `Cookie`
            header verbatim (`request.headers.get("Cookie")`); this
            module never parses or constructs cookie values itself.
        timeout (float): Seconds to wait for a response over HTTP. A
            caller whose endpoint can legitimately run long (e.g.
            `generate_sector_neighborhood`) passes a larger value.

    Returns:
        tuple[dict or None, list[str]]: The parsed JSON body (`None` for
            an empty response), and every `Set-Cookie` response header
            verbatim (to relay back to the browser as-is).

    Raises:
        NotFoundError: On a 404 response.
        ApiError: On any other non-2xx response (with `.status_code` set
            -- see that class's docstring), or if the API can't be
            reached/returns an unparseable body.
    """
    response = _send(method, path, json_body=json_body, cookie_header=cookie_header, timeout=timeout)
    parsed_body = _parse_json(response.body) if response.body else None
    return parsed_body, response.set_cookie_headers


def _log_call(method, url, start, outcome):
    """One debug-log line per API call, attributed to the `apiclient`
    function that made it (`get_system`, `auth_me`, ...)."""
    log.debug(f"API {method} {url} -> {outcome} in {(time.perf_counter() - start) * 1000:.1f}ms", stacklevel=4)


def _build_query(params):
    if not params:
        return ""
    if isinstance(params, dict):
        pairs = [(key, value) for key, value in params.items() if value is not None and value != ""]
    else:
        pairs = [(key, value) for key, value in params if value is not None and value != ""]
    return urllib.parse.urlencode(pairs)


def _error_detail(response):
    """The API's own `{"error": ...}` message from an error response,
    else its raw body, else the status's reason phrase."""
    body = response.body
    try:
        parsed = json.loads(body)
    except ValueError:
        return body or response.reason
    if isinstance(parsed, dict):
        return parsed.get("error", body)
    return body


# ---------------------------------------------------------------------
# Typed wrappers -- one per endpoint a page needs, so a page never
# builds a `/api/...` path/query string by hand.
# ---------------------------------------------------------------------

def _require_db(db):
    """
    Every wrapper below except `list_databases` (which lists across the
    whole server, not one chosen schema) takes a `db` -- typically a
    page's `web.helpers.db_name()`, forwarded straight through.
    Omitting it wouldn't error (the API falls back to its own configured
    default database, see `api/routes.py`'s `get_db`), but a missing `db`
    is a bug in the caller, so it is raised eagerly instead, matching the old direct-database
    `dbutil.resolve_db_name`'s own "No database specified." check.

    Raises:
        NotFoundError: If `db` is empty/`None`.
    """
    if not db:
        raise NotFoundError("No database specified.")


def list_databases():
    """Returns `GET /api/databases`'s `items` list -- see that route's
    own docstring for the shape."""
    return _request("/databases")["items"]


def get_sectors(db, limit=None, offset=None, sort=None, descending=False, quadrants=(), facets=False):
    """Returns `GET /api/sectors`'s full paginated envelope
    (`items`/`total`/`limit`/`offset`, plus `facets` when asked). `sort`,
    `quadrants` and `facets` are the Sectors table's sort and filter (UX.41)."""
    _require_db(db)
    params = [("db", db), ("limit", limit), ("offset", offset), ("sort", sort),
              ("order", "desc" if descending else None)]
    params += [("quadrant", value) for value in quadrants] + [("facets", "1" if facets else None)]
    return _request("/sectors", [(key, value) for key, value in params if value is not None])


def get_sector(db, sector_id):
    """Returns `GET /api/sectors/<id>`'s detail dict -- see
    `queryDb.sector_detail`'s docstring for the shape."""
    _require_db(db)
    return _request(f"/sectors/{sector_id}", {"db": db})


def get_systems(db, star_type=None, sector_id=None, limit=None, offset=None, sort=None, descending=False,
                binary=None, placement=None, octants=(), facets=False):
    """
    Returns `GET /api/systems`'s full paginated envelope.

    Args:
        sector_id: An int, the literal string `"none"` (standalone
            systems), or `None` (no sector filter) -- passed straight
            through as the `sector_id` query parameter.
        sort, descending, binary (True/False), placement (`"sector"` or
            `"standalone"`), octants, facets: the Systems tables' sort and
            filters (UX.41).
    """
    _require_db(db)
    params = [("db", db), ("star_type", star_type), ("sector_id", sector_id), ("limit", limit),
              ("offset", offset), ("sort", sort), ("order", "desc" if descending else None),
              ("binary", None if binary is None else ("yes" if binary else "no")), ("placement", placement)]
    params += [("octant", value) for value in octants] + [("facets", "1" if facets else None)]
    return _request("/systems", [(key, value) for key, value in params if value is not None])


def get_system(db, system_id):
    """Returns `GET /api/systems/<id>`'s detail dict -- see
    `queryDb.system_detail`'s docstring for the shape."""
    _require_db(db)
    return _request(f"/systems/{system_id}", {"db": db})


def get_system_text(db, system_id, fmt):
    """Returns `GET /api/systems/<id>/text`'s `{id, format, content}` --
    the full wiki page as `fmt` (`"wikitext"` or `"markdown"`)."""
    _require_db(db)
    return _request(f"/systems/{system_id}/text", {"db": db, "format": fmt})


def get_system_sections(db, system_id):
    """Returns `GET /api/systems/<id>/sections` -- the page's Markdown split
    per body, see `planetgen.db.render.render_system_sections`."""
    _require_db(db)
    return _request(f"/systems/{system_id}/sections", {"db": db})


def get_system_scene(db, system_id):
    """Returns `GET /api/systems/<id>/scene` -- the 3D system view's data,
    see `planetgen.web.maps.systemscene`."""
    _require_db(db)
    return _request(f"/systems/{system_id}/scene", {"db": db})


def get_near(db, params):
    """Returns `GET /api/near` (NAV.43): everything within a distance of a
    place. `params` are the query parameters (`from` or `point`,
    `distance`, `kinds`, `limit`, `offset`)."""
    _require_db(db)
    return _request("/near", {"db": db, **params})


def get_nav(db, from_ref, to_ref):
    """Returns `GET /api/nav`'s course/route dict -- see `docs/api.md`'s
    "NAV" section for the full shape. `from_ref`/`to_ref` are object
    references (`planetgen.galaxy.objectref`): a system, a body in one, or
    a standalone phenomenon."""
    _require_db(db)
    return _request("/nav", {"db": db, "from": from_ref, "to": to_ref})


def get_galaxy_sectors(db):
    """Returns `GET /api/galaxy/sectors`'s `items` list (every
    galaxy-placed sector)."""
    _require_db(db)
    return _request("/galaxy/sectors", {"db": db})["items"]


def get_galaxy_phenomena(db):
    """Returns `GET /api/galaxy/phenomena`'s `items` list (every
    galaxy-placed nebula/asteroid field) -- see `queryDb.
    galaxy_placed_phenomena`'s docstring for the shape."""
    _require_db(db)
    return _request("/galaxy/phenomena", {"db": db})["items"]


def get_galaxy_shape(db):
    """Returns `GET /api/galaxy/shape`'s `shape` dict -- the galaxy's
    stored density-skeleton shape (`planetgen plan`'s output), or
    `None` if that skeleton has never been built."""
    _require_db(db)
    return _request("/galaxy/shape", {"db": db})["shape"]


def get_bright_star_status(db):
    """`GET /api/galaxy/shape`'s `bright_stars`: whether the bright-star
    scatter has run (`scattered`), its `min_luminosity_sol` and `seed`,
    and `default_min_luminosity_sol`."""
    _require_db(db)
    return _request("/galaxy/shape", {"db": db}).get("bright_stars")


def get_version_warning(db):
    """`GET /api/galaxy/shape`'s `version_warning` (DB.7): a sentence when
    the galaxy holds sectors generated by another PlanetGen version, else
    `None`."""
    _require_db(db)
    return _request("/galaxy/shape", {"db": db}).get("version_warning")


def get_bright_stars_in_cell(db, ring_index, layer_index, ring_slot_index):
    """`GET /api/galaxy/bright-stars`'s `items`: the pre-placed bright
    stars in one sector cell not yet built into a system, brightest
    first."""
    _require_db(db)
    return _request("/galaxy/bright-stars", {
        "db": db, "ring": ring_index, "layer": layer_index, "slot": ring_slot_index,
    })["items"]


def get_galaxy_tiles(db, tile_keys):
    """
    Returns `GET /api/galaxy/tiles`'s payload (`tiles`/`edge_pc`/
    `has_shape` -- see `queryDb.galaxy_tiles`). Callers go
    through `planetgen/web/lib/tilecache.py`, which only asks for tiles it hasn't
    already cached on disk.

    Args:
        db (str): The `?db=` value.
        tile_keys (list[str]): `level/ix/iy/iz` keys.
    """
    _require_db(db)
    return _request("/galaxy/tiles", {"db": db, "tiles": ",".join(tile_keys)})


def get_galaxy_stage(db, at=None):
    """Returns `GET /api/galaxy/stage`'s payload (`at`/`child_m`/
    `children`/`sectors` -- see `queryDb.galaxy_stage`). Callers go
    through `planetgen/web/lib/tilecache.py`, which caches each stage on disk."""
    _require_db(db)
    return _request("/galaxy/stage", {"db": db, "at": at})


def get_territories(db):
    """Returns `GET /api/territories`' payload (`points`: owned systems
    with galaxy-frame parsec positions and their polity's color;
    `polities`: each one's `id`, `capital_pc` and `reach_ly`)."""
    _require_db(db)
    return _request("/territories", {"db": db})


def get_polities(db, limit=None, offset=None, sort=None, descending=False, governments=(), eras=(),
                 facets=False):
    """Returns `GET /api/polities`' paginated envelope (`items`/`total`/
    `limit`/`offset`, plus `facets` when asked), polities by name unless
    `sort`ed; `governments` and `eras` filter (the Polities table, UX.41)."""
    _require_db(db)
    params = [("db", db), ("limit", limit), ("offset", offset), ("sort", sort),
              ("order", "desc" if descending else None)]
    params += [("government", value) for value in governments] + [("era", value) for value in eras]
    params.append(("facets", "1" if facets else None))
    return _request("/polities", [(key, value) for key, value in params if value is not None])


def get_polity(db, polity_id, limit=None, offset=None, sort=None, descending=False):
    """Returns `GET /api/polities/<id>`: one polity plus a page of its
    `systems` (nearest the capital first, or `sort`ed by `name`/`distance`).
    Raises `NotFoundError` for an unknown polity."""
    _require_db(db)
    params = [("db", db), ("limit", limit), ("offset", offset), ("sort", sort),
              ("order", "desc" if descending else None)]
    return _request(f"/polities/{int(polity_id)}", [(key, value) for key, value in params if value is not None])


def get_species_list(db, spacefaring=None, limit=None, offset=None, sort=None, descending=False, eras=(),
                     facets=False):
    """Returns `GET /api/species`' paginated envelope, species by name unless
    `sort`ed; `spacefaring` (`True`/`False`) and `eras` filter, `None` lists
    every one (the Species table, UX.41)."""
    _require_db(db)
    flag = None if spacefaring is None else ("1" if spacefaring else "0")
    params = [("db", db), ("spacefaring", flag), ("limit", limit), ("offset", offset), ("sort", sort),
              ("order", "desc" if descending else None)]
    params += [("era", value) for value in eras] + [("facets", "1" if facets else None)]
    return _request("/species", [(key, value) for key, value in params if value is not None])


def get_species(db, species_id):
    """Returns `GET /api/species/<id>`. Raises `NotFoundError` for an
    unknown species."""
    _require_db(db)
    return _request(f"/species/{int(species_id)}", {"db": db})


def get_planet_species(db, planet_id):
    """The dominant species whose homeworld is this planet
    (`GET /api/planets/<id>/species`), or `None` when it has none."""
    _require_db(db)
    try:
        return _request(f"/planets/{int(planet_id)}/species", {"db": db})
    except NotFoundError:
        return None


def get_system_owner(db, system_id):
    """The polity that owns a system (`GET /api/systems/<id>/owner`'s
    `owner`: `polity_id`, `polity_name`, `color`, `distance_ly`), or
    `None` when no polity does."""
    _require_db(db)
    return _request(f"/systems/{int(system_id)}/owner", {"db": db})["owner"]


POPULATION_NONE = {"generated": False, "species": False, "polities": False, "territories": False}
"""dict: `get_population_status`'s answer when there is nothing to show."""


def get_population_status(db):
    """
    What population data exists (`GET /api/population`): `generated`,
    `species`, `polities` and `territories` booleans. Before that endpoint
    existed it is worked out from the species and polity counts. Any
    failure (a database from before population, schema v44) means none.
    """
    _require_db(db)
    try:
        try:
            return dict(POPULATION_NONE, **_request("/population", {"db": db}))
        except NotFoundError:
            species = get_species_list(db, limit=1)["total"] > 0
            polities = get_polities(db, limit=1)["total"] > 0
            return {"generated": species, "species": species, "polities": polities, "territories": polities}
    except (ApiError, NotFoundError):
        return dict(POPULATION_NONE)


def get_galaxy_locate(db, q):
    """Returns `GET /api/galaxy/locate`'s `matches` (sectors and systems
    named like `q`, each with its sector address -- see
    `queryDb.galaxy_locate`)."""
    _require_db(db)
    return _request("/galaxy/locate", {"db": db, "q": q})["matches"]


def get_nebula_shape(db, nebula_id, lod="low"):
    """Returns `GET /api/nebulae/<id>/shape`'s payload: one nebula's mesh
    (`vertices` in nebula-radius units from its center, `faces`), at the
    `"low"` or `"full"` level of detail (GEN.75)."""
    _require_db(db)
    return _request(f"/nebulae/{int(nebula_id)}/shape", {"db": db, "lod": lod})


def get_nebula_surroundings(db, nebula_id):
    """Returns `GET /api/nebulae/<id>/surroundings`' payload: the brightest
    stars round one nebula, for its page's 3D view (MAP.105)."""
    _require_db(db)
    return _request(f"/nebulae/{int(nebula_id)}/surroundings", {"db": db})


def get_galaxy_changes(db, since=None):
    """Returns `GET /api/galaxy/changes`' payload (`stamp`/`state`/`full`/
    `tiles`/`stages`) -- which cube tiles and drill-down stages changed since `since`, an earlier call's
    `state` (see `queryDb.galaxy_changes`). `planetgen/web/lib/tilecache.py` uses it to
    refresh only the tiles an edit touched."""
    _require_db(db)
    return _request("/galaxy/changes", {"db": db, "since": since})

def get_phenomena(db, limit=None, offset=None, sort=None, descending=False, types=(), descriptors=(),
                  placed=None, facets=False):
    """Returns `GET /api/phenomena`'s full paginated envelope
    (`items`/`total`/`limit`/`offset`, plus `facets` when asked) -- see
    `queryDb.list_phenomena`'s docstring for the shape. `sort`, `types`,
    `descriptors` and `placed` (True/False) are the table's sort and
    filters (UX.41)."""
    _require_db(db)
    params = [("db", db), ("limit", limit), ("offset", offset), ("sort", sort),
              ("order", "desc" if descending else None)]
    params += [("type", value) for value in types] + [("descriptor", value) for value in descriptors]
    params += [("placed", None if placed is None else ("yes" if placed else "no")),
               ("facets", "1" if facets else None)]
    return _request("/phenomena", [(key, value) for key, value in params if value is not None])


def get_phenomenon(db, phenomenon_type, phenomenon_id):
    """Returns `GET /api/phenomena/<type>/<id>`'s detail dict -- see
    `queryDb.phenomenon_detail`'s docstring for the shape."""
    _require_db(db)
    return _request(f"/phenomena/{phenomenon_type}/{phenomenon_id}", {"db": db})


def get_object(db, ref):
    """Returns `GET /api/objects/<ref>`'s dict (kind, name, parent chain,
    sibling references, positions) -- see `queryDb.resolve_object`."""
    _require_db(db)
    return _request(f"/objects/{ref}", {"db": db})


def get_search(db, texts, tags, sizes=None, limit=None, offsets=None, panels=None):
    """
    Runs `GET /api/search` and returns its response dict -- see
    `queryDb.search`'s docstring for the full shape.

    Args:
        db (str): The `db` query parameter.
        texts (dict): `{"sector_q", "system_q", "star_q", "planet_q",
            "moon_q"}` -> search term.
        tags (dict): `{facet: iterable of value}` -- one entry per
            `queryDb.SEARCH_TAG_FACETS` name; a repeated query parameter
            per active value (e.g. `class=M&class=K`).
        sizes (dict, optional): `{"star", "planet", "moon"} -> (min_km,
            max_km)`, each bound `None` for "unbounded" -- sent as
            `<entity>_min_radius_km`/`<entity>_max_radius_km`, omitting
            either bound that's `None`. An absent key (or `sizes` itself
            being `None`) sends no size filter for that entity.
        limit (int, optional): Rows per result panel (the API's own
            default when `None`).
        offsets (dict, optional): `{panel: offset}` -- sent as
            `<panel>_offset`, one page per result panel.
        panels (iterable, optional): Run just these result panels.
    """
    _require_db(db)
    pairs = [("db", db)]
    pairs.extend((key, value) for key, value in texts.items() if value)
    for facet, values in tags.items():
        pairs.extend((facet, value) for value in values)
    for entity, size_range in (sizes or {}).items():
        if size_range is None:
            continue
        min_km, max_km = size_range
        if min_km is not None:
            pairs.append((f"{entity}_min_radius_km", min_km))
        if max_km is not None:
            pairs.append((f"{entity}_max_radius_km", max_km))
    if limit is not None:
        pairs.append(("limit", limit))
    for panel, offset in (offsets or {}).items():
        pairs.append((f"{panel}_offset", offset))
    if panels:
        pairs.append(("panels", ",".join(panels)))
    return _request("/search", pairs)


# ---------------------------------------------------------------------
# Admin auth -- see `docs/api.md`'s "Authentication" section. Every
# function below takes/returns the incoming/outgoing `Cookie`/`Set-Cookie`
# headers verbatim (see `_auth_request`'s docstring) -- the admin pages
# (`web/admin_pages.py`) are the only callers, and relay them to/from the
# browser.
# ---------------------------------------------------------------------

def auth_login(username, password, cookie_header=None):
    """
    `POST /api/auth/login`. `cookie_header` (the browser's own) carries
    the trusted-device cookie (SEC.22), if any.

    Returns:
        tuple[dict, list[str]]: `({"username", "must_change_credentials"},
            set_cookie_headers)` on success.

    Raises:
        ApiError: `status_code == 401` for a wrong username/password --
            callers should catch this specifically and show an inline
            "invalid username or password" message rather than a generic
            error page.
    """
    return _auth_request("POST", "/auth/login", json_body={"username": username, "password": password},
                         cookie_header=cookie_header)


def auth_login_totp(pending, code, cookie_header=None):
    """
    `POST /api/auth/login/totp`, the second step of a two-factor login
    (SEC.26). Returns `(identity, set_cookie_headers)` like `auth_login`.

    Raises:
        ApiError: `status_code == 401` for a wrong code or a stale
            `pending` (the message says which), 429 while locked.
    """
    return _auth_request("POST", "/auth/login/totp", json_body={"pending": pending, "code": code},
                         cookie_header=cookie_header)


def auth_totp_status(cookie_header):
    """`GET /api/auth/totp`: `{"enabled", "recovery_codes_left"}`."""
    return _auth_request("GET", "/auth/totp", cookie_header=cookie_header)[0]


def auth_totp_setup(cookie_header, current_password):
    """`POST /api/auth/totp/setup`: `{"secret", "uri", "qr_svg"}`."""
    return _auth_request("POST", "/auth/totp/setup", cookie_header=cookie_header,
                         json_body={"current_password": current_password})[0]


def auth_totp_confirm(cookie_header, code):
    """`POST /api/auth/totp/confirm`: `{"recovery_codes": [...]}`."""
    return _auth_request("POST", "/auth/totp/confirm", cookie_header=cookie_header, json_body={"code": code})[0]


def auth_totp_disable(cookie_header, current_password, code):
    """`POST /api/auth/totp/disable`. Returns the `Set-Cookie` headers to
    relay (a fresh trusted-device cookie: turning it off forgets every
    other device, TEST.46)."""
    return _auth_request("POST", "/auth/totp/disable", cookie_header=cookie_header,
                         json_body={"current_password": current_password, "code": code})[1]


def auth_logout(cookie_header):
    """`POST /api/auth/logout`. Returns the `Set-Cookie` headers to relay
    (clears the session cookie) -- a no-op, not an error, if the caller
    was already logged out (no valid cookie to begin with)."""
    try:
        _body, set_cookie_headers = _auth_request("POST", "/auth/logout", cookie_header=cookie_header)
        return set_cookie_headers
    except ApiError as exc:
        if exc.status_code == 401:
            return []
        raise


def auth_me(cookie_header):
    """
    `GET /api/auth/me` equivalent (this module has no bare-GET-with-cookie
    helper, so this goes through `_auth_request` too).

    Returns:
        dict or None: `{"username", "must_change_credentials"}`, or
            `None` if `cookie_header` names no valid session (a normal,
            expected "not logged in" outcome -- not an error).
    """
    try:
        body, _set_cookie_headers = _auth_request("GET", "/auth/me", cookie_header=cookie_header)
        return body
    except ApiError as exc:
        if exc.status_code == 401:
            return None
        raise


def auth_change_credentials(cookie_header, current_password, new_username, new_password):
    """
    `POST /api/auth/change-credentials`.

    Returns:
        tuple[dict, list[str]]: The updated identity and fresh session's
            `Set-Cookie` headers (the old session is invalidated server-
            side -- see `adminAuth.change_credentials`).

    Raises:
        ApiError: `status_code == 400` for a wrong current password, a
            taken username, or a `new_password` failing policy -- the
            message is safe to show the caller as-is (see
            `adminAuth.AuthError`).
    """
    return _auth_request(
        "POST", "/auth/change-credentials", cookie_header=cookie_header,
        json_body={
            "current_password": current_password,
            "new_username": new_username,
            "new_password": new_password,
        },
    )


def auth_list_api_keys(cookie_header):
    """`GET /api/auth/api-keys` equivalent -- returns the `items` list
    (label/timestamps only, never the key itself)."""
    body, _set_cookie_headers = _auth_request("GET", "/auth/api-keys", cookie_header=cookie_header)
    return body["items"]


def auth_create_api_key(cookie_header, label):
    """`POST /api/auth/api-keys` -- returns `{"id", "label", "key"}`; `key`
    is the raw key, shown this once (see `adminAuth.create_api_key`)."""
    body, _set_cookie_headers = _auth_request(
        "POST", "/auth/api-keys", cookie_header=cookie_header, json_body={"label": label},
    )
    return body


def auth_revoke_api_key(cookie_header, key_id):
    """`DELETE /api/auth/api-keys/<id>`."""
    _auth_request("DELETE", f"/auth/api-keys/{key_id}", cookie_header=cookie_header)


def admin_stats(cookie_header, db):
    """`GET /api/admin/stats?db=` -- server health and statistics about
    one database, for `html/adminstats.py`."""
    _require_db(db)
    body, _set_cookie_headers = _auth_request(
        "GET", f"/admin/stats?{_build_query({'db': db})}", cookie_header=cookie_header,
    )
    return body


def admin_login_failures(cookie_header, limit=None):
    """`GET /api/admin/login-failures?limit=` -- the newest refused
    sign-ins (`{"action", "username", "ip", "created_at"}` each), newest
    first."""
    query = _build_query({"limit": limit})
    path = f"/admin/login-failures?{query}" if query else "/admin/login-failures"
    body, _set_cookie_headers = _auth_request("GET", path, cookie_header=cookie_header)
    return body["items"]


def admin_lockouts(cookie_header):
    """`GET /api/admin/lockouts` -- `{"items": [{"scope", "subject",
    "retry_after", "locked_until", "level"}], "proxy_warning"}`."""
    body, _set_cookie_headers = _auth_request("GET", "/admin/lockouts", cookie_header=cookie_header)
    return body


def admin_naming_key(cookie_header, db):
    """`GET /api/admin/naming-key?db=` -- the galaxy's naming key and when it
    was drawn or last changed (GEN.70)."""
    _require_db(db)
    body, _set_cookie_headers = _auth_request(
        "GET", f"/admin/naming-key?{_build_query({'db': db})}", cookie_header=cookie_header)
    return body


def admin_set_naming_key(cookie_header, db, key=None, draw=False):
    """`POST /api/admin/naming-key?db=` -- sets the key (8 hex digits) or,
    with `draw`, a fresh random one (GEN.70)."""
    _require_db(db)
    body, _set_cookie_headers = _auth_request(
        "POST", f"/admin/naming-key?{_build_query({'db': db})}", cookie_header=cookie_header,
        json_body={"draw": True} if draw else {"key": key})
    return body


def admin_generation_stats(cookie_header):
    """`GET /api/admin/generation-stats` -- the server's generation speed
    per density bucket and each galaxy's size per star system (PERF.10)."""
    body, _set_cookie_headers = _auth_request("GET", "/admin/generation-stats", cookie_header=cookie_header)
    return body


def admin_work(cookie_header, limit=None, offset=None):
    """`GET /api/admin/work?limit=&offset=` -- the work queue's state and
    one page of job trees (ADM.10)."""
    query = _build_query({"limit": limit, "offset": offset})
    body, _set_cookie_headers = _auth_request("GET", f"/admin/work?{query}", cookie_header=cookie_header)
    return body


def admin_work_tree(cookie_header, node_id):
    """`GET /api/admin/work/<id>` -- the whole job tree holding a node."""
    body, _set_cookie_headers = _auth_request("GET", f"/admin/work/{urllib.parse.quote(str(node_id), safe='')}",
                                              cookie_header=cookie_header)
    return body


def admin_work_control(cookie_header, node_id, action):
    """`POST /api/admin/work/<id>/control` -- pause, resume or cancel a
    node and its subtree. Returns whether it took."""
    body, _set_cookie_headers = _auth_request("POST", f"/admin/work/{urllib.parse.quote(str(node_id), safe='')}/control",
                                              cookie_header=cookie_header, json_body={"action": action})
    return body["ok"]


def admin_work_delete(cookie_header, node_id):
    """`POST /api/admin/work/<id>/delete` -- deletes a finished tree."""
    body, _set_cookie_headers = _auth_request("POST", f"/admin/work/{urllib.parse.quote(str(node_id), safe='')}/delete",
                                              cookie_header=cookie_header)
    return body["ok"]


def admin_work_clear_finished(cookie_header):
    """`POST /api/admin/work/clear-finished` -- deletes every finished tree; returns how many."""
    body, _set_cookie_headers = _auth_request("POST", "/admin/work/clear-finished", cookie_header=cookie_header)
    return body["cleared"]


def admin_work_queue(cookie_header, action):
    """`POST /api/admin/work/queue` -- pauses or resumes the whole queue."""
    body, _set_cookie_headers = _auth_request("POST", "/admin/work/queue", cookie_header=cookie_header,
                                              json_body={"action": action})
    return body["ok"]


def admin_work_clear_lease(cookie_header):
    """`POST /api/admin/work/lease/clear` -- frees a stale lease; returns
    the holder it cleared, or `None`."""
    body, _set_cookie_headers = _auth_request("POST", "/admin/work/lease/clear", cookie_header=cookie_header)
    return body["cleared"]


def admin_lift_lockout(cookie_header, scope=None, subject=None, lift_all=False):
    """`POST /api/admin/lockouts/lift` -- one lockout, or every one with
    `lift_all`. Returns how many were lifted."""
    payload = {"all": True} if lift_all else {"scope": scope, "subject": subject}
    body, _set_cookie_headers = _auth_request("POST", "/admin/lockouts/lift", cookie_header=cookie_header,
                                              json_body=payload)
    return body["lifted"]


def admin_duplicate_names(cookie_header, db, limit=None, offset=None):
    """`GET /api/admin/duplicate-names?db=&limit=&offset=` -- one page of the
    names the uniqueness rules had to decorate, each with the rows that
    carry it."""
    _require_db(db)
    query = _build_query({"db": db, "limit": limit, "offset": offset})
    body, _set_cookie_headers = _auth_request("GET", f"/admin/duplicate-names?{query}", cookie_header=cookie_header)
    return body


# ---------------------------------------------------------------------
# Wiki publishing (schema.sql's "v22" header note) -- POST .../wiki, and
# GET /api/wiki-config so a page knows which backend(s), if any, to offer
# before even showing an "Upload to Wiki" form.
# ---------------------------------------------------------------------

def get_wiki_config():
    """Returns `GET /api/wiki-config`'s `{"wikijs": bool, "mediawiki":
    bool}` -- which backend(s) are configured deployment-wide, so a page
    knows whether to offer an "Upload to Wiki" choice at all (and which
    backend(s)) before rendering one."""
    return _request("/wiki-config")


def upload_system_to_wiki(cookie_header, db, system_id, backend, path=None):
    """
    `POST /api/systems/<id>/wiki` -- publishes the system's already-
    generated page to `backend`, and records the resulting URL on
    `star_systems.wikijs_url`/`mediawiki_url`.

    Args:
        backend (str): `"wikijs"` or `"mediawiki"`.
        path (str, optional): The target page's path/slug -- required for
            `"wikijs"` (which addresses a page separately from its title);
            ignored for `"mediawiki"` (its title -- the system's own name
            -- is its address).

    Returns:
        dict: The new page's `{"id", "path", "title", "url"}`, or, when the
            upload was queued and is still running (PERF.24),
            `{"status": "accepted", "job_id", "status_url"}`.

    Raises:
        ApiError: `status_code == 409` if a page already exists at the
            target path/title (`system.py` should show this as an inline
            "already uploaded" message, not a generic error);
            `status_code == 501` if `backend` isn't configured
            deployment-wide (see `get_wiki_config`).
    """
    _require_db(db)
    url_path = f"/systems/{system_id}/wiki?{_build_query({'db': db})}"
    body = {"backend": backend}
    if path:
        body["path"] = path
    result, _set_cookie_headers = _auth_request("POST", url_path, json_body=body, cookie_header=cookie_header)
    return result


def upload_sector_to_wiki(cookie_header, db, sector_id, backend, path=None):
    """`POST /api/sectors/<id>/wiki` -- same shape/errors as
    `upload_system_to_wiki`, but for a freshly generated sector-summary
    page (sectors have no persisted rendered page of their own), recording
    the result on `sectors.wiki_url`."""
    _require_db(db)
    url_path = f"/sectors/{sector_id}/wiki?{_build_query({'db': db})}"
    body = {"backend": backend}
    if path:
        body["path"] = path
    result, _set_cookie_headers = _auth_request("POST", url_path, json_body=body, cookie_header=cookie_header)
    return result


def admin_set_sector_wiki_url(cookie_header, db, sector_id, wiki_url):
    """`PATCH /api/sectors/<id>` `{"wiki_url": ...}` -- the admin "manually
    set the wiki link" affordance (`html/admin.py`), no upload involved.
    `wiki_url` may be `None`/empty to clear it back to "no page yet"."""
    _require_db(db)
    url_path = f"/sectors/{sector_id}?{_build_query({'db': db})}"
    body = {"wiki_url": wiki_url or None}
    _auth_request("PATCH", url_path, json_body=body, cookie_header=cookie_header)


# ---------------------------------------------------------------------
# Facilities (schema v42) -- see `docs/api.md`'s "Facilities" section.
# ---------------------------------------------------------------------

def get_system_facilities(db, system_id):
    """Returns `GET /api/systems/<id>/facilities`'s `items`: every
    facility in the system (`queryDb._facility_dict`'s shape)."""
    _require_db(db)
    return _request(f"/systems/{system_id}/facilities", {"db": db})["items"]


def get_sector_facilities(db, sector_id):
    """Returns `GET /api/sectors/<id>/facilities`'s `items`: the
    stand-alone facilities parked in the sector and those on its asteroid
    fields."""
    _require_db(db)
    return _request(f"/sectors/{sector_id}/facilities", {"db": db})["items"]


def get_facility_orbit(db, host_type, host_id, distance_km=None):
    """
    `GET /api/facilities/orbit` -- the orbit an orbital facility would get
    around a star, planet or moon, without saving anything.

    Returns:
        dict: `distance_km`, `period_years`, `orbital_speed_kms`.

    Raises:
        ApiError: `status_code == 400` with the API's reason (a distance
            inside the host, a host nothing orbits).
        NotFoundError: For an unknown host.
    """
    _require_db(db)
    return _request("/facilities/orbit", {
        "db": db, "host_type": host_type, "host_id": host_id, "distance_km": distance_km,
    })


def create_facility(cookie_header, db, body):
    """`POST /api/facilities` with `body` (`name`, `kind`, `placement`,
    `host_type`, `host_id`, plus optional `distance_km`/`description`).
    Returns the new facility's id. A rule the API refuses is an
    `ApiError` with `status_code == 400`; a missing host a
    `NotFoundError`."""
    _require_db(db)
    result, _set_cookie_headers = _auth_request(
        "POST", f"/facilities?{_build_query({'db': db})}", json_body=body, cookie_header=cookie_header,
    )
    return result["id"]


def delete_facility(cookie_header, db, facility_id):
    """`DELETE /api/facilities/<id>`. `NotFoundError` if it is already
    gone."""
    _require_db(db)
    _auth_request("DELETE", f"/facilities/{facility_id}?{_build_query({'db': db})}", cookie_header=cookie_header)


_NEIGHBORHOOD_GENERATION_TIMEOUT_SECONDS = 1800
"""float: `generate_sector_neighborhood` below reads the whole region's
density to count its candidate slots (and, with `estimate_only` off, to
refuse a run that won't fit), which on a large radius takes longer than a
quick CRUD call; `_TIMEOUT_SECONDS` alone would make a working request
look like a failure."""


def generate_sector_neighborhood(cookie_header, sector_id, radius_ly=None, estimate_only=False):
    """
    `POST /api/sectors/<id>/generate-neighborhood` -- the not-yet-generated
    sectors within `radius_ly` (`None` for the API's own default, 12 pc) of
    this already galaxy-placed sector. The admin-only "generate more
    sectors around this one" action on `sector.py`. Uses
    `_NEIGHBORHOOD_GENERATION_TIMEOUT_SECONDS` rather than this module's
    usual, much shorter timeout -- see that constant's own docstring.

    Returns:
        dict: With `estimate_only`, nothing is generated: the counts
            (`already_existed`/`candidates`/...) and the size and time
            `estimate` (PERF.3). Without it, the run is queued (PERF.24)
            and this is the `202` body: `job_id` and `status_url`
            (`GET /api/jobs/<id>`).

    Raises:
        ApiError: `status_code == 404` if the sector doesn't exist or was
            never placed in a galaxy; 507 when the database disk can't
            hold it.
    """
    json_body = {"radius_ly": radius_ly} if radius_ly is not None else {}
    if estimate_only:
        json_body["estimate_only"] = True
    body, _set_cookie_headers = _auth_request(
        "POST", f"/sectors/{sector_id}/generate-neighborhood", cookie_header=cookie_header,
        json_body=json_body,
        timeout=_NEIGHBORHOOD_GENERATION_TIMEOUT_SECONDS,
    )
    return body


_EDIT_TIMEOUT_SECONDS = 600
"""float: An admin edit (`admin_edit`) can regenerate a whole sector, which
takes longer than a quick CRUD call."""


def admin_edit(cookie_header, db, method, path, body=None):
    """
    One admin editing call (TODO ADM.1, `api/edits.py`, and the system
    `PATCH`/`DELETE`): `method` and `path` under `/api` (for example
    `"POST", "/planets/12/regenerate"`), with an optional JSON `body`.
    Returns the parsed answer (`summary`, `warnings`, ...). A refusal is an
    `ApiError` (409 when facilities would be lost or the edit can't be
    made); a missing target a `NotFoundError`.
    """
    _require_db(db)
    result, _set_cookie_headers = _auth_request(
        method, f"{path}?{_build_query({'db': db})}", json_body=body, cookie_header=cookie_header,
        timeout=_EDIT_TIMEOUT_SECONDS,
    )
    return result or {}
