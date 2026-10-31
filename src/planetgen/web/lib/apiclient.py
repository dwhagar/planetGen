# planetgen/web/lib/apiclient.py

"""
Client for the planetGen Flask API (`planetgen/api/`), used by every page in
`planetgen/web/` instead of querying the database directly (in-process there,
see `web/transport.py`; over HTTP anywhere else).

This is the concrete result of moving the interim web browser onto the
API: every page here used to open its own read-only MySQL connection
(`html/lib/dbutil.py`'s previous `open_readonly`/`resolve_db_name`) and
run its own bespoke SQL; now every page is a thin HTTP client over
`GET /api/...` (see `docs/reference/api.md`), and the database-querying logic those
pages used to duplicate lives once, in `planetgen.db.query`, shared with the API
itself. `planetgen/web/lib/fmt.py` still holds the formatting-only helpers
(`esc`, `linkify_location`, `format_density`) that have nothing to do
with fetching data.

Stdlib only (`urllib.request`) -- same "nothing beyond a system Python 3"
deployment story this project's docs already describe for `html/`
(see `docs/reference/html-interface.md`), now extended to reach the API process
over HTTP rather than a database socket.

The API must be reachable for this browser to work at all now -- see
`API_BASE_URL` below and `docs/deployment/apache.md` for how it's mounted
alongside `html/` in a real deployment.
"""

import json
import time
import urllib.error
import urllib.parse
import urllib.request

from planetgen.util import log
from planetgen.util.user_agent import user_agent
from planetgen.util.settings import get_settings

from planetgen.web.lib import pagecache, slowpages

API_BASE_URL = get_settings().api_base_url
"""str: Base URL of the planetGen API's `/api` mount point. Defaults to
the same host (see `examples/apache/planetgen.conf.example`) -- only used
by the HTTP transport, not by the Flask pages' in-process one --
override via the `PLANETGEN_API_BASE_URL` Apache `SetEnv` (or shell env,
for local testing), or `config.json`'s `api_base_url` (see
`docs/reference/config.md`), when the API is deployed at a different host/port,
e.g. `http://127.0.0.1:5000/api` for `python src/html/wsgi.py`'s own dev
server."""

_TIMEOUT_SECONDS = 30
"""int: Was 15 -- confirmed too tight for GET /api/search specifically
(TimeoutError in production): that one endpoint always runs its full
facet+autocomplete query set up front regardless of whether any filter is
active (db.query.search), which used to mean an unindexed full-table
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


cache_stats = {"hits": 0, "misses": 0}
"""dict: How many cached GETs were answered from the response cache and how many had to be fetched and were answered, in
this process (the warm-up reports them, PERF.104). A GET the API refused (a 404, a 400) built no answer and counts as neither:
the sector page of a sector with a grid address asks for its Profile, which a database with no planned galaxy refuses
every time (MAP.289)."""


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
    request = urllib.request.Request(f"{API_BASE_URL}{target}", data=data, headers={"User-Agent": user_agent(), **headers},
                                     method=method)
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
    if method == "GET":
        from planetgen.web import rebuild_gate
        rebuild_gate.check(target)    # a visitor's heavy read while the caches are rebuilt: CacheRebuilding, a 503
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
        slowpages.Pending: Only inside `slowpages.allow_pending()`: the answer is slow and a job is building it.
    """
    query = _build_query(params)
    target = f"{path}?{query}" if query else path
    return _parse_json(_body_of(path, target, _cache_db(params)))


def fetch_target(target):
    """
    Makes the cached GET `target` (a path under `/api` with its query) and returns its raw JSON body, never raising
    `slowpages.Pending`: a slow request is built here, in the caller, which is how a page job (PERF.105) builds
    the answer it stores in the page cache.
    """
    path, _, query = target.partition("?")
    db = next((value for key, value in urllib.parse.parse_qsl(query) if key == "db"), None) or None
    return _body_of(path, target, db, queueing=False)


def _body_of(path, target, db, queueing=True):
    """The raw JSON body of the GET `target`: from the response cache, else from the API (and stored, when it can
    be). With a cache, a miss on a page that opted in (`slowpages.allow_pending`) that the cost rule calls slow
    raises `slowpages.Pending` when `queueing`, instead of building it."""
    cache = _response_cache() if _response_cache is not None and pagecache.is_cacheable(path) else None
    if cache is None:
        return _send("GET", target).body
    body = cache.get(db, target)
    if body is not None:
        cache_stats["hits"] += 1
    if body is None:
        if queueing and slowpages.pending_allowed():
            if db is not None:
                body = slowpages.delivered(db, target)    # a job built it and the cache would not take it (UX.234)
                if body is not None:
                    return body
            pending = slowpages.pending_for(db, path, target)
            if pending is not None:
                raise pending
        token = cache.begin(db)
        started = time.perf_counter()
        body = _send("GET", target).body
        cache_stats["misses"] += 1    # an answer built; one the API refused (a 404, a 400) was not, and is asked again next time
        seconds = time.perf_counter() - started
        stored = cache.put(db, target, body, token)
        if db is not None:
            slowpages.record_duration(db, target, seconds)
            if not queueing and not stored:
                slowpages.deliver(db, target, body)
    return body


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


def get_sectors(db, limit=None, offset=None, sort=None, descending=False, quadrants=(), facets=False, after=None):
    """Returns `GET /api/sectors`'s full paginated envelope
    (`items`/`total`/`limit`/`offset`, plus `facets` when asked). `sort`,
    `quadrants` and `facets` are the Sectors table's sort and filter (UX.41); `after` is the sector ID the previous
    page ended on (PERF.137)."""
    _require_db(db)
    params = [("db", db), ("limit", limit), ("offset", offset), ("after", after), ("sort", sort),
              ("order", "desc" if descending else None)]
    params += [("quadrant", value) for value in quadrants] + [("facets", "1" if facets else None)]
    return _request("/sectors", [(key, value) for key, value in params if value is not None])


def get_sector(db, sector_id):
    """Returns `GET /api/sectors/<id>`'s detail dict -- see
    `db.query.sector_detail`'s docstring for the shape."""
    _require_db(db)
    return _request(f"/sectors/{sector_id}", {"db": db})


def get_systems(db, star_type=None, sector_id=None, limit=None, offset=None, sort=None, descending=False,
                binary=None, placement=None, octants=(), facets=False, after=None, start=None):
    """
    Returns `GET /api/systems`'s full paginated envelope.

    Args:
        sector_id: An int, the literal string `"none"` (standalone
            systems), or `None` (no sector filter) -- passed straight
            through as the `sector_id` query parameter.
        sort, descending, binary (True/False), placement (`"sector"` or
            `"standalone"`), octants, facets: the Systems tables' sort and
            filters (UX.41).
        after: The system ID the previous page ended on (PERF.75); with the name sort the page starts after
            it and `offset` is ignored.
        start: The jump-by-value box's text (PERF.75, the API's `from`): the list begins at the first name (name
            sort) or sector (sector sort) from it on.
    """
    _require_db(db)
    params = [("db", db), ("star_type", star_type), ("sector_id", sector_id), ("limit", limit),
              ("offset", offset), ("after", after), ("from", start or None), ("sort", sort), ("order", "desc" if descending else None),
              ("binary", None if binary is None else ("yes" if binary else "no")), ("placement", placement)]
    params += [("octant", value) for value in octants] + [("facets", "1" if facets else None)]
    return _request("/systems", [(key, value) for key, value in params if value is not None])


def get_uncharted_systems(db, limit=None, offset=None, sort=None, descending=False, after=None, start=None):
    """Returns `GET /api/uncharted-systems`'s envelope (UX.87): a page of
    the scattered stars no system was built around yet, with the `total`
    waiting. `sort` is one of `db.query.UNCHARTED_SYSTEM_SORTS`; `after` is the star id the previous page
    ended on (PERF.75); `start` is the jump-by-value box's text (PERF.135, the API's `from`): the list begins at the
    first star at or after that value in the sort."""
    _require_db(db)
    params = [("db", db), ("limit", limit), ("offset", offset), ("after", after), ("from", start or None), ("sort", sort),
              ("order", "desc" if descending else "asc" if sort else None)]
    return _request("/uncharted-systems", [(key, value) for key, value in params if value is not None])


def get_uncharted_system(db, bright_star_id):
    """Returns `GET /api/uncharted-systems/<id>`'s dict (UX.114) -- see `db.query.uncharted_system_detail`."""
    _require_db(db)
    return _request(f"/uncharted-systems/{int(bright_star_id)}", {"db": db})


def get_system(db, system_id):
    """Returns `GET /api/systems/<id>`'s detail dict -- see
    `db.query.system_detail`'s docstring for the shape."""
    _require_db(db)
    return _request(f"/systems/{system_id}", {"db": db})


def get_life_air(db, kind, body_id):
    """Returns `GET /api/planets/<id>/life-air` (or `moons`)'s `{id, life_class_id, life_air}`: the data of the Life and air chart,
    `life_air` being `None` when the body has nothing to chart (UX.237)."""
    _require_db(db)
    return _request(f"/{kind}s/{body_id}/life-air", {"db": db})


def get_life_readings(db, kind, body_id):
    """Returns `GET /api/planets/<id>/life-readings` (or `moons`)'s `{id, life_class_id, readings}`: the three words for each kind of life
    on the body, `readings` being empty when there is nothing to say (UX.251)."""
    _require_db(db)
    return _request(f"/{kind}s/{body_id}/life-readings", {"db": db})


def get_system_text(db, system_id, fmt):
    """Returns `GET /api/systems/<id>/text`'s `{id, format, content}` --
    the system's description as `fmt` (`"markdown"` or `"text"`)."""
    _require_db(db)
    return _request(f"/systems/{system_id}/text", {"db": db, "format": fmt})


def get_system_description(db, system_id):
    """Returns `GET /api/systems/<id>/description` -- the system's Document,
    see `planetgen.describe`."""
    _require_db(db)
    return _request(f"/systems/{system_id}/description", {"db": db})


def get_system_scene(db, system_id):
    """Returns `GET /api/systems/<id>/scene` -- the 3D system view's data,
    see `planetgen.web.maps.systemscene`."""
    _require_db(db)
    return _request(f"/systems/{system_id}/scene", {"db": db})


def get_table_estimates(db):
    """
    `GET /api/table-estimates` -- the database's own estimate of each table's row count, for the cost rule that
    decides which page requests are slow (`slowpages`). Never cached: it is asked about once a minute.

    Returns:
        dict: Table name -> estimated rows.
    """
    _require_db(db)
    query = _build_query({"db": db})
    return _parse_json(_send("GET", f"/table-estimates?{query}").body)["tables"]


def get_near(db, params):
    """Returns `GET /api/near` (NAV.43): everything within a distance of a
    place. `params` are the query parameters (`from` or `point`,
    `distance`, `kinds`, `limit`, `offset`)."""
    _require_db(db)
    return _request("/near", {"db": db, **params})


def get_nav(db, from_ref, to_ref, stay=None):
    """Returns `GET /api/nav`'s course/route dict -- see `docs/reference/api.md`'s
    "NAV" section for the full shape. `from_ref`/`to_ref` are object
    references (`planetgen.galaxy.objectref`): a system, a body in one, or
    a standalone phenomenon. `stay` is the minutes spent at each stop of
    the route (NAV.11; the API's default, 0, when `None`)."""
    _require_db(db)
    params = {"db": db, "from": from_ref, "to": to_ref}
    if stay is not None:
        params["stay"] = stay
    return _request("/nav", params)


def get_nav_chart(db, from_ref, to_ref, border=False):
    """Returns `GET /api/nav/chart` (NAV.48): the uncharted sectors that block the course between two
    objects -- `unknown_hops`, `cells`, `count`, `outside_galaxy`, `route_distance_ly`, `bypass` and
    `confirm_over`; `border` adds the uncharted cells next to them."""
    _require_db(db)
    params = {"db": db, "from": from_ref, "to": to_ref}
    if border:
        params["border"] = 1
    return _request("/nav/chart", params)


def get_galaxy_sectors(db):
    """Returns `GET /api/galaxy/sectors`'s `items` list (every
    galaxy-placed sector)."""
    _require_db(db)
    return _request("/galaxy/sectors", {"db": db})["items"]


def get_superbubbles(db):
    """Returns `GET /api/galaxy/superbubbles`'s `superbubbles` list (GEN.254): `{"x", "y", "z", "radius_pc"}` per
    cavity, galaxy-frame parsecs; empty when the galaxy has none to show (or the API can't say)."""
    _require_db(db)
    try:
        return _request("/galaxy/superbubbles", {"db": db})["superbubbles"]
    except ApiError:
        return []


def get_galaxy_made(db, since, until=None):
    """Returns `GET /api/galaxy/made` (`total` and `items`): the sectors
    created from `since` to `until`, Unix seconds (ADM.31)."""
    _require_db(db)
    params = {"db": db, "since": since}
    if until is not None:
        params["until"] = until
    return _request("/galaxy/made", params)


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


def get_shell_specs(db):
    """`GET /api/galaxy/shape`'s `shells` (ADM.28): the galaxy's shell
    count, reach and charted sectors, or `None` before a plan."""
    _require_db(db)
    return _request("/galaxy/shape", {"db": db}).get("shells")


def get_version_warning(db):
    """`GET /api/galaxy/shape`'s `version_warning` (DB.7): a sentence when
    the galaxy holds sectors generated by another PlanetGen version, else
    `None`."""
    _require_db(db)
    return _request("/galaxy/shape", {"db": db}).get("version_warning")


def get_bright_stars_in_cell(db, shell_index, band_index, wedge_index):
    """`GET /api/galaxy/bright-stars`'s `items`: the pre-placed bright
    stars in one sector cell not yet built into a system, brightest
    first."""
    _require_db(db)
    return _request("/galaxy/bright-stars", {
        "db": db, "shell": shell_index, "band": band_index, "wedge": wedge_index,
    })["items"]


def get_uncharted_sector(db, shell_index, band_index, wedge_index):
    """`GET /api/galaxy/uncharted`'s payload (MAP.162): one sector cell's
    place (`designation`, `center_pc`, `edge_pc`, `sector_id`) and what the
    scatters left in it (`stars`, `scattered`)."""
    _require_db(db)
    return _request("/galaxy/uncharted", {
        "db": db, "shell": shell_index, "band": band_index, "wedge": wedge_index,
    })


def get_galaxy_profile(db, shell_index, band_index, wedge_index):
    """`GET /api/galaxy/profile`'s payload (API.27, UX.113): the Sector Profile of one cell of the grid, charted or
    not, with `sector_id`, `unresolved` and the `estimate` of filling it."""
    _require_db(db)
    return _request("/galaxy/profile", {"db": db, "shell": shell_index, "band": band_index, "wedge": wedge_index})


def get_galaxy_nodes(db, node_keys):
    """
    Returns `GET /api/galaxy/nodes`' payload (`nodes`/`edge_pc`/`has_shape` -- see `db.nodes.galaxy_nodes`). Callers go
    through `planetgen/web/lib/tilecache.py`, which only asks for nodes it hasn't already cached on disk.

    Args:
        db (str): The `?db=` value.
        node_keys (list[str]): `m.shell.band.wedge` keys.
    """
    _require_db(db)
    return _request("/galaxy/nodes", {"db": db, "keys": ",".join(node_keys)})


def get_galaxy_stage(db, at=None):
    """Returns `GET /api/galaxy/stage`'s payload (`at`/`child_m`/
    `children`/`sectors` -- see `db.query.galaxy_stage`). Callers go
    through `planetgen/web/lib/tilecache.py`, which caches each stage on disk."""
    _require_db(db)
    return _request("/galaxy/stage", {"db": db, "at": at})


def get_galaxy_window(db, at, g):
    """Returns `GET /api/galaxy/window`'s payload (`at`/`g`/`child_m`/
    `cells`/`children` -- see `db.query.galaxy_window`): a window of cells
    bigger than a sector, read from the stored block totals in one query."""
    _require_db(db)
    return _request("/galaxy/window", {"db": db, "at": at, "g": g})


def get_galaxy_nebula_cover(db, at=None, g=None):
    """Returns `GET /api/galaxy/nebula-cover`'s payload (`at`/`child_m`/`cells` -- see `db.query.galaxy_nebula_cover`):
    the share of each cell of a stage (`at`) or window (`at` and `g`) inside a nebula. Read on every request, since
    it changes without a sector changing."""
    _require_db(db)
    return _request("/galaxy/nebula-cover", {"db": db, "at": at, "g": g})


def get_galaxy_charted_regions(db):
    """Returns `GET /api/galaxy/charted-regions`' payload (`regions`, each charted region with its node key, extent, share
    charted and name, by share charted; `edge_pc` -- see `db.region_stats.charted_regions`)."""
    _require_db(db)
    return _request("/galaxy/charted-regions", {"db": db})


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
                     facets=False, tech_ranges=None, after=None):
    """Returns `GET /api/species`' paginated envelope, species by name unless
    `sort`ed; `spacefaring` (`True`/`False`) and `eras` filter, `None` lists
    every one (the Species table, UX.41). `tech_ranges` maps a range name
    (`tech` or a domain such as `energy`) to `(low, high)`, either end
    `None` for open (POP.9). `after` is the id of the species the previous page ended on (PERF.205)."""
    _require_db(db)
    flag = None if spacefaring is None else ("1" if spacefaring else "0")
    params = [("db", db), ("spacefaring", flag), ("limit", limit), ("offset", offset), ("after", after), ("sort", sort),
              ("order", "desc" if descending else None)]
    params += [("era", value) for value in eras] + [("facets", "1" if facets else None)]
    for name, (low, high) in (tech_ranges or {}).items():
        params += [(f"{name}_min", low), (f"{name}_max", high)]
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
        return _request(f"/planets/{planet_id}/species", {"db": db})
    except NotFoundError:
        return None


def get_system_owner(db, system_id):
    """The polity that owns a system (`GET /api/systems/<id>/owner`'s
    `owner`: `polity_id`, `polity_name`, `color`, `distance_ly`), or
    `None` when no polity does."""
    _require_db(db)
    return _request(f"/systems/{system_id}/owner", {"db": db})["owner"]


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
    `db.query.galaxy_locate`)."""
    _require_db(db)
    return _request("/galaxy/locate", {"db": db, "q": q})["matches"]


def get_nebula_shape(db, nebula_id, lod="low"):
    """Returns `GET /api/nebulae/<id>/shape`'s payload: one nebula's mesh
    (`vertices` in nebula-radius units from its center, `faces`), at the
    `"low"` or `"full"` level of detail (GEN.75)."""
    _require_db(db)
    return _request(f"/nebulae/{nebula_id}/shape", {"db": db, "lod": lod})


def get_nebula_surroundings(db, nebula_id):
    """Returns `GET /api/nebulae/<id>/surroundings`' payload: the brightest
    stars round one nebula, for its page's 3D view (MAP.105)."""
    _require_db(db)
    return _request(f"/nebulae/{nebula_id}/surroundings", {"db": db})


def get_galaxy_changes(db, since=None):
    """Returns `GET /api/galaxy/changes`' payload (`stamp`/`state`/`full`/
    `stages`/`nodes`) -- which drill-down stages and nodes changed since `since`, an earlier call's
    `state` (see `db.query.galaxy_changes`). `planetgen/web/lib/tilecache.py` uses it to
    refresh only the nodes an edit touched."""
    _require_db(db)
    return _request("/galaxy/changes", {"db": db, "since": since})

def get_phenomena(db, limit=None, offset=None, sort=None, descending=False, types=(), descriptors=(),
                  placed=None, facets=False, after=None, start=None):
    """Returns `GET /api/phenomena`'s full paginated envelope
    (`items`/`total`/`limit`/`offset`, plus `facets` when asked) -- see
    `db.query.list_phenomena`'s docstring for the shape. `sort`, `types`,
    `descriptors` and `placed` (True/False) are the table's sort and
    filters (UX.41); `after` is the `scatter_id` the previous page ended on, inside the scattered tail (PERF.75);
    `start` is the jump-by-value box's scatter id (PERF.135, the API's `from`): the scattered tail from that id on."""
    _require_db(db)
    params = [("db", db), ("limit", limit), ("offset", offset), ("after", after), ("from", start or None), ("sort", sort),
              ("order", "desc" if descending else None)]
    params += [("type", value) for value in types] + [("descriptor", value) for value in descriptors]
    params += [("placed", None if placed is None else ("yes" if placed else "no")),
               ("facets", "1" if facets else None)]
    return _request("/phenomena", [(key, value) for key, value in params if value is not None])


def get_object_list(db, kind, limit=None, sort=None, descending=False, after=None, name=None, object_class=None,
                    sector_id=None, start=None):
    """Returns `GET /api/object-lists/<kind>`'s `items`/`next`/`total` -- one keyset page of the object index
    (`db.objectlist.list_objects`). `after` is the printed object ID the previous page ended on and `start` begins a
    name-sorted list at the first name from it on (UX.141)."""
    _require_db(db)
    params = [("db", db), ("limit", limit), ("sort", sort), ("order", "desc" if descending else None),
              ("after", after), ("name", name), ("class", object_class), ("sector_id", sector_id),
              ("from", start or None)]
    return _request(f"/object-lists/{kind}", [(key, value) for key, value in params if value is not None])


def get_phenomenon(db, phenomenon_type, phenomenon_id):
    """Returns `GET /api/phenomena/<type>/<id>`'s detail dict -- see
    `db.query.phenomenon_detail`'s docstring for the shape."""
    _require_db(db)
    return _request(f"/phenomena/{phenomenon_type}/{phenomenon_id}", {"db": db})


def get_phenomenon_contents(db, phenomenon_type, phenomenon_id, limit=None):
    """Returns `GET /api/phenomena/<type>/<id>/contents` -- what sits inside a nebula or supernova remnant
    (`db.cloudcontents.cloud_contents`)."""
    _require_db(db)
    return _request(f"/phenomena/{phenomenon_type}/{phenomenon_id}/contents",
                    [("db", db)] + ([("limit", limit)] if limit else []))


def get_uncharted_phenomenon(db, scatter_id):
    """Returns `GET /api/uncharted-phenomena/<id>`'s dict -- see `db.query.scattered_phenomenon_detail`."""
    _require_db(db)
    return _request(f"/uncharted-phenomena/{int(scatter_id)}", {"db": db})


def get_object(db, ref):
    """Returns `GET /api/objects/<ref>`'s dict (kind, name, parent chain,
    sibling references, positions) -- see `db.query.resolve_object`."""
    _require_db(db)
    return _request(f"/objects/{ref}", {"db": db})


def get_object_sheet(db, ref):
    """Returns `GET /api/objects/<ref>/sheet`'s dict (`kind`, `groups` of `{key, label, rows}`), the stored properties of
    a star, planet, moon, belt, comet or phenomenon -- see `describe/sheet.py`."""
    _require_db(db)
    return _request(f"/objects/{ref}/sheet", {"db": db})


def get_object_children(db, ref, limit=None, after=None):
    """Returns `GET /api/objects/<ref>/children`'s dict (`ref`, `kind`, `total`, `items`, `next`), one keyset page of what
    belongs to an object -- see `db/children.py`."""
    _require_db(db)
    return _request(f"/objects/{ref}/children", {"db": db, "limit": limit, "after": after})


def get_search(db, texts, tags, sizes=None, limit=None, offsets=None, panels=None, sorts=None, stars=None, ranges=None):
    """
    Runs `GET /api/search` and returns its response dict -- see
    `db.query.search`'s docstring for the full shape.

    Args:
        db (str): The `db` query parameter.
        texts (dict): `{"sector_q", "system_q", "star_q", "planet_q",
            "moon_q"}` -> search term.
        tags (dict): `{facet: iterable of value}` -- one entry per
            `db.query.SEARCH_TAG_FACETS` name; a repeated query parameter
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
        panels (iterable, optional): Run just these result panels; an empty one runs none, which still answers
            the facets and the name suggestions.
        sorts (dict, optional): `{"star", "planet", "moon"} -> sort key` (`-` in front for the best first), sent as
            `star_sort`/`planet_sort`/`moon_sort` (UX.214; the keys of the physical values, UX.226).
        stars (str, optional): The stars filter (UX.239), sent as `stars`: a count, a count and a plus, or a name.
        ranges (dict, optional): `{entity: {key: (least, most, unit)}}`, the range filters on the physical values (UX.226), sent as
            `<entity>_min_<key>`, `<entity>_max_<key>` and `<entity>_<key>_unit`; either bound `None` is left out, and a `None` unit
            is the key's own.
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
    for entity, sort in (sorts or {}).items():
        if sort:
            pairs.append((f"{entity}_sort", sort))
    if stars:
        pairs.append(("stars", stars))
    for entity, by_key in (ranges or {}).items():
        for key, (least, most, unit) in by_key.items():
            if least is not None:
                pairs.append((f"{entity}_min_{key}", least))
            if most is not None:
                pairs.append((f"{entity}_max_{key}", most))
            if unit:
                pairs.append((f"{entity}_{key}_unit", unit))
    if panels is not None:
        pairs.append(("panels", ",".join(panels) or "none"))
    return _request("/search", pairs)


# ---------------------------------------------------------------------
# Admin auth -- see `docs/reference/api.md`'s "Authentication" section. Every
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


def auth_confirm_password(cookie_header, current_password):
    """`POST /api/auth/confirm-password` (ADM.86): types the password again for the Owner-tier settings;
    `{"minutes": N}`. A wrong password is an `ApiError` (400), a locked-out one 429."""
    return _auth_request("POST", "/auth/confirm-password", cookie_header=cookie_header,
                         json_body={"current_password": current_password})[0]


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


def auth_avatar_set(cookie_header, data_base64):
    """`POST /api/auth/avatar` -- sets the account's picture from base64 text (SEC.36). A picture that can't be used
    raises `ApiError` (422) with the reason."""
    return _auth_request("POST", "/auth/avatar", cookie_header=cookie_header, json_body={"data": data_base64})[0]


def auth_avatar_remove(cookie_header):
    """`DELETE /api/auth/avatar` -- the account has no picture (SEC.36)."""
    return _auth_request("DELETE", "/auth/avatar", cookie_header=cookie_header)[0]


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
            side -- see `admin.auth.change_credentials`).

    Raises:
        ApiError: `status_code == 400` for a wrong current password, a
            taken username, or a `new_password` failing policy -- the
            message is safe to show the caller as-is (see
            `admin.auth.AuthError`).
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


def auth_create_api_key(cookie_header, label, scopes=None, expires_days=None):
    """`POST /api/auth/api-keys` -- returns `{"id", "label", "scopes",
    "expires_days", "key"}`; `key` is the raw key, shown this once (see
    `admin.auth.create_api_key`)."""
    payload = {"label": label}
    if scopes:
        payload["scopes"] = list(scopes)
    if expires_days is not None:
        payload["expires_days"] = expires_days
    body, _set_cookie_headers = _auth_request(
        "POST", "/auth/api-keys", cookie_header=cookie_header, json_body=payload,
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


def admin_settings(cookie_header):
    """`GET /api/admin/settings` -- every option with its help, limits, default, value and source, the file's
    hash, why it can't be written, and the options waiting for a restart (ADM.84)."""
    body, _set_cookie_headers = _auth_request("GET", "/admin/settings", cookie_header=cookie_header)
    return body


def admin_settings_plan(cookie_header, changes):
    """`POST /api/admin/settings/plan` -- the diff a save of `changes` (`{option: value or {"reset": true}}`)
    would make, or the error at each option; nothing is written (ADM.84)."""
    body, _set_cookie_headers = _auth_request("POST", "/admin/settings/plan", cookie_header=cookie_header,
                                              json_body={"changes": changes})
    return body


def admin_settings_save(cookie_header, changes, expect_hash=None):
    """`POST /api/admin/settings` -- checks and saves `changes`. A refused save is an answer with `saved` false
    and the error at each option; a file changed since `expect_hash` or one that can't be written raises
    `ApiError` (409) (ADM.84)."""
    body, _set_cookie_headers = _auth_request("POST", "/admin/settings", cookie_header=cookie_header,
                                              json_body={"changes": changes, "expect_hash": expect_hash})
    return body


def admin_images(cookie_header):
    """`GET /api/admin/images` -- each upload kind with the stored file in use (ADM.44)."""
    body, _set_cookie_headers = _auth_request("GET", "/admin/images", cookie_header=cookie_header)
    return body


def admin_image_upload(cookie_header, kind, data_base64):
    """`POST /api/admin/images/<kind>` -- stores the picture (base64 text) and saves the setting that names it;
    a picture that can't be used raises `ApiError` (422) with the reason (ADM.44)."""
    body, _set_cookie_headers = _auth_request("POST", f"/admin/images/{kind}", cookie_header=cookie_header,
                                              json_body={"data": data_base64})
    return body


def admin_image_remove(cookie_header, kind):
    """`DELETE /api/admin/images/<kind>` -- the built-in picture (or no icon) applies again (ADM.44)."""
    body, _set_cookie_headers = _auth_request("DELETE", f"/admin/images/{kind}", cookie_header=cookie_header)
    return body


def admin_backups(cookie_header):
    """`GET /api/admin/backups` -- the backup archives, newest first, with their parts (ADM.85)."""
    body, _set_cookie_headers = _auth_request("GET", "/admin/backups", cookie_header=cookie_header)
    return body


def admin_backup_take(cookie_header):
    """`POST /api/admin/backups` -- Back up now: one archive of every part (ADM.85)."""
    body, _set_cookie_headers = _auth_request("POST", "/admin/backups", cookie_header=cookie_header, json_body={})
    return body


def admin_backup(cookie_header, name):
    """`GET /api/admin/backups/<name>` -- one archive, whether it verifies, and the diff of its settings part
    against today's config.json (ADM.85). Raises `NotFoundError` for a name that isn't an archive."""
    body, _set_cookie_headers = _auth_request(
        "GET", f"/admin/backups/{urllib.parse.quote(name, safe='')}", cookie_header=cookie_header)
    return body


def admin_backup_restore(cookie_header, name, parts=None, expect_hash=None):
    """`POST /api/admin/backups/<name>/restore` -- restores the whole archive (`parts` None) or the named parts
    (ADM.85). A changed settings file or a backup that doesn't fit raises `ApiError` (409)."""
    body, _set_cookie_headers = _auth_request(
        "POST", f"/admin/backups/{urllib.parse.quote(name, safe='')}/restore", cookie_header=cookie_header,
        json_body={"parts": parts, "expect_hash": expect_hash})
    return body


def admin_settings_export(cookie_header):
    """`GET /api/admin/settings/export` -- the settings an admin may change, without secrets (ADM.85)."""
    body, _set_cookie_headers = _auth_request("GET", "/admin/settings/export", cookie_header=cookie_header)
    return body


def admin_settings_import(cookie_header, document, expect_hash=None, check_only=False):
    """`POST /api/admin/settings/import` -- checks (`check_only`) or saves an export (ADM.85)."""
    body, _set_cookie_headers = _auth_request(
        "POST", "/admin/settings/import", cookie_header=cookie_header,
        json_body={"document": document, "expect_hash": expect_hash, "check_only": check_only})
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


def admin_reset_generation_stats(cookie_header):
    """`POST /api/admin/generation-stats/reset` -- deletes every recorded
    generation rate (PERF.32). Returns how many rows went."""
    body, _set_cookie_headers = _auth_request("POST", "/admin/generation-stats/reset", cookie_header=cookie_header,
                                              json_body={})
    return body["deleted"]


def admin_work(cookie_header, limit=None, offset=None):
    """`GET /api/admin/work?limit=&offset=` -- the work queue's state and one page of runs (ADM.10, PERF.99)."""
    query = _build_query({"limit": limit, "offset": offset})
    body, _set_cookie_headers = _auth_request("GET", f"/admin/work?{query}", cookie_header=cookie_header)
    return body


def admin_work_warm(cookie_header, database):
    """`GET /api/admin/work/warm?database=` -- what the cache warm-up has done on `database`, phase by phase (PERF.141)."""
    body, _set_cookie_headers = _auth_request("GET", f"/admin/work/warm?{_build_query({'database': database})}",
                                              cookie_header=cookie_header)
    return body


def admin_work_run(cookie_header, run_id):
    """`GET /api/admin/work/<id>` -- one run with its stages and the last events of its log (PERF.99)."""
    body, _set_cookie_headers = _auth_request("GET", f"/admin/work/{urllib.parse.quote(str(run_id), safe='')}",
                                              cookie_header=cookie_header)
    return body


def admin_work_run_log(cookie_header, run_id):
    """`GET /api/admin/work/<id>/log` -- the run's whole log as text."""
    body, _set_cookie_headers = _auth_request("GET", f"/admin/work/{urllib.parse.quote(str(run_id), safe='')}/log",
                                              cookie_header=cookie_header)
    return body["text"]


def admin_work_control(cookie_header, run_id, action):
    """`POST /api/admin/work/<id>/control` -- cancels a live run or clears one that stopped answering. Returns whether it
    took."""
    body, _set_cookie_headers = _auth_request("POST", f"/admin/work/{urllib.parse.quote(str(run_id), safe='')}/control",
                                              cookie_header=cookie_header, json_body={"action": action})
    return body["ok"]


def admin_work_delete(cookie_header, run_id):
    """`POST /api/admin/work/<id>/delete` -- deletes a finished run's record."""
    body, _set_cookie_headers = _auth_request("POST", f"/admin/work/{urllib.parse.quote(str(run_id), safe='')}/delete",
                                              cookie_header=cookie_header)
    return body["ok"]


def admin_work_clear_finished(cookie_header):
    """`POST /api/admin/work/clear-finished` -- deletes every finished run's record; returns how many."""
    body, _set_cookie_headers = _auth_request("POST", "/admin/work/clear-finished", cookie_header=cookie_header)
    return body["cleared"]


def admin_work_queue(cookie_header, action):
    """`POST /api/admin/work/queue` -- pauses or resumes the whole queue."""
    body, _set_cookie_headers = _auth_request("POST", "/admin/work/queue", cookie_header=cookie_header,
                                              json_body={"action": action})
    return body["ok"]


def admin_work_stop_all(cookie_header):
    """`POST /api/admin/work/queue` `stop-all` -- stops everything and holds (PERF.186). Returns the summary
    (`runs`, `queued`, `job`, `killed`, `survivors`, `already`)."""
    body, _set_cookie_headers = _auth_request("POST", "/admin/work/queue", cookie_header=cookie_header,
                                              json_body={"action": "stop-all"}, timeout=120)
    return body["summary"]


def admin_work_workers(cookie_header, count=None, ludicrous=False):
    """`GET`, or with `count` / `ludicrous` `POST`, `/api/admin/work/workers` -- the bulk worker count (ADM.15):
    `{"wanted", "live", "ludicrous", "cpus"}`."""
    if count is None and not ludicrous:
        body, _set_cookie_headers = _auth_request("GET", "/admin/work/workers", cookie_header=cookie_header)
    else:
        body, _set_cookie_headers = _auth_request("POST", "/admin/work/workers", cookie_header=cookie_header,
                                                  json_body={"ludicrous": True} if ludicrous else {"count": count})
    return body


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
# Facilities (schema v42) -- see `docs/reference/api.md`'s "Facilities" section.
# ---------------------------------------------------------------------

def get_system_facilities(db, system_id):
    """Returns `GET /api/systems/<id>/facilities`'s `items`: every
    facility in the system (`db.query._facility_dict`'s shape)."""
    _require_db(db)
    return _request(f"/systems/{system_id}/facilities", {"db": db})["items"]


def get_system_habitable_species(db, system_id):
    """Returns `GET /api/systems/<id>/habitable-species`'s `bodies`: each rocky planet and moon some species could
    live on, with those species (POP.12)."""
    _require_db(db)
    return _request(f"/systems/{system_id}/habitable-species", {"db": db})["bodies"]


def get_sector_facilities(db, sector_id):
    """Returns `GET /api/sectors/<id>/facilities`'s `items`: the
    stand-alone facilities parked in the sector and those on its asteroid
    fields."""
    _require_db(db)
    return _request(f"/sectors/{sector_id}/facilities", {"db": db})["items"]


def get_containers(db, kind=None):
    """Returns `GET /api/containers`'s `containers`: every stored globular cluster and superbubble (or those of `kind`)."""
    _require_db(db)
    return _request("/containers", {"db": db, "kind": kind})["containers"]


def get_container(db, printed_id):
    """Returns `GET /api/containers/<ID>`'s dict -- one globular cluster or superbubble (`db.containers.get`'s shape)."""
    _require_db(db)
    return _request(f"/containers/{printed_id}", {"db": db})


def get_container_contents(db, printed_id, limit=None):
    """Returns `GET /api/containers/<ID>/contents`'s dict: `counts` by kind, `total` and the first `items`."""
    _require_db(db)
    return _request(f"/containers/{printed_id}/contents", {"db": db, "limit": limit})


def get_facility(db, facility_id):
    """Returns `GET /api/facilities/<id>`'s dict -- one facility (`db.query._facility_dict`'s shape)."""
    _require_db(db)
    return _request(f"/facilities/{facility_id}", {"db": db})


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


def get_facility_types(db, active_only=True):
    """`GET /api/facility-types`'s `items` (POP.10): every type with its `fields`, the active ones only unless
    `active_only` is off."""
    _require_db(db)
    return _request("/facility-types", {"db": db, "active": "1" if active_only else None})["items"]


def facility_type_write(method, cookie_header, db, path, body=None):
    """
    One admin write to `/api/facility-types<path>` (POP.10): `POST` a type (`path` empty), `PATCH` or `DELETE` one
    (`/<id>`), `POST` a field (`/<id>/fields`), `PATCH` or `DELETE` one (`/<id>/fields/<field_id>`) or `POST
    /examples`. Returns the parsed answer. A rule the API refuses is an `ApiError` with its status.
    """
    _require_db(db)
    result, _set_cookie_headers = _auth_request(
        method, f"/facility-types{path}?{_build_query({'db': db})}", json_body=body, cookie_header=cookie_header,
    )
    return result


def update_facility(cookie_header, db, facility_id, body):
    """`PATCH /api/facilities/<id>` with `body` (any of `name`, `description`, `type_id`, `affiliation_polity_id`,
    `affiliation_name`, `ratings`, `ratings_note`, `fields`). Returns the facility. A refusal is an `ApiError`
    with its status; a missing facility a `NotFoundError`."""
    _require_db(db)
    result, _set_cookie_headers = _auth_request(
        "PATCH", f"/facilities/{facility_id}?{_build_query({'db': db})}", json_body=body, cookie_header=cookie_header,
    )
    return result


def suggest_facility_ratings(cookie_header, db, facility_id):
    """`POST /api/facilities/<id>/suggest-ratings`: fills the ratings nobody set. Returns the dimensions filled."""
    _require_db(db)
    result, _set_cookie_headers = _auth_request(
        "POST", f"/facilities/{facility_id}/suggest-ratings?{_build_query({'db': db})}", cookie_header=cookie_header,
    )
    return result["filled"]


def delete_facility(cookie_header, db, facility_id):
    """`DELETE /api/facilities/<id>`. `NotFoundError` if it is already
    gone."""
    _require_db(db)
    _auth_request("DELETE", f"/facilities/{facility_id}?{_build_query({'db': db})}", cookie_header=cookie_header)


_NEIGHBORHOOD_GENERATION_TIMEOUT_SECONDS = 1800
"""float: `generate_sector_neighborhood` below reads the whole region's
density to count its candidate cells (and, with `estimate_only` off, to
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
