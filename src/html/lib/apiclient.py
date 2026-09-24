# html/lib/apiclient.py

"""
JSON HTTP client for the planetGen Flask API (`html/api/`), used by every
CGI script in `html/` instead of querying the database directly.

This is the concrete result of moving the interim web browser onto the
API: every page here used to open its own read-only MySQL connection
(`html/lib/dbutil.py`'s previous `open_readonly`/`resolve_db_name`) and
run its own bespoke SQL; now every page is a thin HTTP client over
`GET /api/...` (see `docs/api.md`), and the database-querying logic those
pages used to duplicate lives once, in `queryDb.py`, shared with the API
itself. `html/lib/fmt.py` still holds the formatting-only helpers
(`esc`, `linkify_location`, `format_density`) that have nothing to do
with fetching data.

Stdlib only (`urllib.request`) -- same "nothing beyond a system Python 3"
deployment story this project's docs already describe for `html/`
(see `docs/html-interface.md`), now extended to reach the API process
over HTTP rather than a database socket.

The API must be reachable for this browser to work at all now -- see
`API_BASE_URL` below and `docs/apache-deployment.md` for how it's mounted
alongside `html/` in a real deployment.
"""

import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

# stellarObjects/ lives at src/stellarObjects/ (src layout); this file is
# at src/html/lib/, two levels down from src/ -- add src/ to sys.path the
# same way every other html/ script already does (see e.g. nav.py) so
# `stellarObjects.appconfig` is importable here too.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from stellarObjects.appconfig import load_config  # noqa: E402

API_BASE_URL = os.environ.get("PLANETGEN_API_BASE_URL") or load_config()["api_base_url"]
"""str: Base URL of the planetGen API's `/api` mount point. Defaults to
the same host this CGI script itself runs on (see
`examples/apache/planetgen.conf.example`'s `WSGIScriptAlias /api`) --
override via the `PLANETGEN_API_BASE_URL` Apache `SetEnv` (or shell env,
for local testing), or `config.json`'s `api_base_url` (see
`docs/config.md`), when the API is deployed at a different host/port,
e.g. `http://127.0.0.1:5000/api` for `python src/html/wsgi.py`'s own dev
server running alongside a locally-invoked CGI script."""

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
    -- callers turn this into a 404 page, same as before this module
    existed (`html/lib/page.run` already handles it)."""


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
        The parsed JSON body (usually a `dict`, `/systems/<id>/near`
        returns a bare `list`).

    Raises:
        NotFoundError: On a 404 response.
        ApiError: On any other non-2xx response, or if the API can't be
            reached/returns an unparseable body.
    """
    url = f"{API_BASE_URL}{path}"
    query = _build_query(params)
    if query:
        url = f"{url}?{query}"

    try:
        with urllib.request.urlopen(url, timeout=_TIMEOUT_SECONDS) as response:
            body = response.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        detail = _error_detail(exc)
        if exc.code == 404:
            raise NotFoundError(detail)
        raise ApiError(f"planetGen API error ({exc.code}): {detail}", status_code=exc.code)
    except urllib.error.URLError as exc:
        raise ApiError(f"Could not reach the planetGen API at {API_BASE_URL}: {exc.reason}")

    try:
        return json.loads(body)
    except ValueError as exc:
        raise ApiError(f"planetGen API returned an unparseable response: {exc}")


def _auth_request(method, path, json_body=None, cookie_header=None, timeout=_TIMEOUT_SECONDS):
    """
    Runs one JSON request against the API supporting any HTTP method and
    an optional request body/`Cookie` header -- the primitive every
    `auth_*` function below builds on, kept separate from `_request`
    (GET-only, no body/cookie support) rather than complicating that
    function's simpler, far more common case.

    Args:
        method (str): `"POST"`, `"DELETE"`, etc.
        path (str): The path under `API_BASE_URL`, e.g. `"/auth/login"`.
        json_body (dict, optional): Sent as the request body
            (`Content-Type: application/json`) if given.
        cookie_header (str, optional): Forwarded as-is as the outgoing
            `Cookie` header -- callers pass the CGI request's own
            `HTTP_COOKIE` environment variable verbatim (see
            `page.incoming_cookie_header`); this module never parses or
            constructs cookie values itself, only relays them.
        timeout (float): Seconds to wait for a response. Defaults to
            `_TIMEOUT_SECONDS`, the same as every other request this
            module makes; a caller whose endpoint can legitimately run
            long (e.g. `generate_sector_neighborhood`'s own batch sector
            generation) passes a larger value explicitly instead of this
            module silently timing out a request that was still working.

    Returns:
        tuple[dict or None, list[str]]: The parsed JSON body (`None` for
            an empty response), and every `Set-Cookie` response header
            verbatim (for `page.py` to relay back to the browser as-is --
            this module never parses those either).

    Raises:
        NotFoundError: On a 404 response.
        ApiError: On any other non-2xx response (with `.status_code` set
            -- see that class's docstring), or if the API can't be
            reached/returns an unparseable body.
    """
    url = f"{API_BASE_URL}{path}"
    data = None
    headers = {}
    if json_body is not None:
        data = json.dumps(json_body).encode("utf-8")
        headers["Content-Type"] = "application/json"
    if cookie_header:
        headers["Cookie"] = cookie_header

    request = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            raw_body = response.read().decode("utf-8")
            set_cookie_headers = response.headers.get_all("Set-Cookie") or []
    except urllib.error.HTTPError as exc:
        detail = _error_detail(exc)
        if exc.code == 404:
            raise NotFoundError(detail)
        raise ApiError(f"planetGen API error ({exc.code}): {detail}", status_code=exc.code)
    except urllib.error.URLError as exc:
        raise ApiError(f"Could not reach the planetGen API at {API_BASE_URL}: {exc.reason}")

    parsed_body = None
    if raw_body:
        try:
            parsed_body = json.loads(raw_body)
        except ValueError as exc:
            raise ApiError(f"planetGen API returned an unparseable response: {exc}")
    return parsed_body, set_cookie_headers


def _build_query(params):
    if not params:
        return ""
    if isinstance(params, dict):
        pairs = [(key, value) for key, value in params.items() if value is not None and value != ""]
    else:
        pairs = [(key, value) for key, value in params if value is not None and value != ""]
    return urllib.parse.urlencode(pairs)


def _error_detail(http_error):
    try:
        body = http_error.read().decode("utf-8", errors="replace")
    except Exception:
        return http_error.reason
    try:
        return json.loads(body).get("error", body)
    except ValueError:
        return body or http_error.reason


# ---------------------------------------------------------------------
# Typed wrappers -- one per endpoint a CGI page needs, so a page never
# builds a `/api/...` path/query string by hand.
# ---------------------------------------------------------------------

def _require_db(db):
    """
    Every wrapper below except `list_databases` (which lists across the
    whole server, not one chosen schema) takes a `db` -- typically a
    page's own `db` (from `nav_params()`), forwarded straight through.
    Omitting it wouldn't error (the API falls back to its own configured
    default database, see `api/routes.py`'s `get_db`), but every one of
    these pages carries `db` in every link it renders (`page.post_link`'s
    hidden `db` field, e.g. on a `sector.py`/`system.py` link), so a
    missing `db` here would quietly build a page entirely out of a
    different database than the one every link on it claims to be
    showing. Raised eagerly instead, matching the old direct-database
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


def get_sectors(db, limit=None, offset=None):
    """Returns `GET /api/sectors`'s full paginated envelope
    (`items`/`total`/`limit`/`offset`)."""
    _require_db(db)
    return _request("/sectors", {"db": db, "limit": limit, "offset": offset})


def get_sector(db, sector_id):
    """Returns `GET /api/sectors/<id>`'s detail dict -- see
    `queryDb.sector_detail`'s docstring for the shape."""
    _require_db(db)
    return _request(f"/sectors/{sector_id}", {"db": db})


def get_systems(db, star_type=None, sector_id=None, limit=None, offset=None):
    """
    Returns `GET /api/systems`'s full paginated envelope.

    Args:
        sector_id: An int, the literal string `"none"` (standalone
            systems), or `None` (no sector filter) -- passed straight
            through as the `sector_id` query parameter.
    """
    _require_db(db)
    return _request("/systems", {
        "db": db, "star_type": star_type, "sector_id": sector_id, "limit": limit, "offset": offset,
    })


def get_system(db, system_id):
    """Returns `GET /api/systems/<id>`'s detail dict -- see
    `queryDb.system_detail`'s docstring for the shape."""
    _require_db(db)
    return _request(f"/systems/{system_id}", {"db": db})


def get_systems_near(db, system_id, radius):
    """Returns `GET /api/systems/<id>/near`'s bare list of
    `{id, name, distance_ly}`."""
    _require_db(db)
    return _request(f"/systems/{system_id}/near", {"db": db, "radius": radius})


def get_nav(db, from_id, to_id, from_kind="system", to_kind="system", from_type=None, to_type=None):
    """Returns `GET /api/nav`'s course/route dict -- see `docs/api.md`'s
    "NAV" section for the full shape. `from_kind`/`to_kind` default to
    `"system"`; pass `"phenomenon"` (plus the matching `from_type`/
    `to_type`, one of `queryDb._PHENOMENON_TYPE_TO_TABLE`'s keys) to
    route to/from a standalone nebula/asteroid field/black hole/neutron
    star instead."""
    _require_db(db)
    return _request("/nav", {
        "db": db, "from": from_id, "to": to_id,
        "from_kind": from_kind, "to_kind": to_kind,
        "from_type": from_type, "to_type": to_type,
    })


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
    stored density-skeleton shape (`generate.py plan`'s output), or
    `None` if that skeleton has never been built."""
    _require_db(db)
    return _request("/galaxy/shape", {"db": db})["shape"]


def get_galaxy_view(db, cx, cy, cz, radius_pc):
    """
    Returns `GET /api/galaxy/view`'s full live-viewport payload
    (`placed`/`planned`/`density`/`edge_pc`/`has_shape` -- see
    `queryDb.galaxy_view`'s docstring for the shape) for the interactive
    3D Galaxy Map's own camera position. Unlike every other `get_galaxy_*`
    function here, this is called from `html/galaxy_view.py` (the
    browser-facing JSON proxy the map's own client-side JS fetches from
    directly, on every camera move) rather than from a page's `handler()`
    at render time -- see that script's own docstring.

    Args:
        db (str): The `?db=` value.
        cx, cy, cz (float): The view center, galaxy-frame parsecs.
        radius_pc (float): The view radius, parsecs.
    """
    _require_db(db)
    return _request("/galaxy/view", {"db": db, "cx": cx, "cy": cy, "cz": cz, "radius_pc": radius_pc})



def get_galaxy_tiles(db, tile_keys, density_key=None):
    """
    Returns `GET /api/galaxy/tiles`'s payload (`tiles`/`density`/
    `edge_pc`/`has_shape` -- see `queryDb.galaxy_tiles`). Callers go
    through `lib/tilecache.py`, which only asks for tiles it hasn't
    already cached on disk.

    Args:
        db (str): The `?db=` value.
        tile_keys (list[str]): `level/ix/iy/iz` keys.
        density_key (str or None): A tile key to anchor a density cloud on.
    """
    _require_db(db)
    return _request("/galaxy/tiles", {
        "db": db, "tiles": ",".join(tile_keys), "density": density_key,
    })


def get_galaxy_stamp(db):
    """Returns `GET /api/galaxy/stamp`'s `stamp` -- the token tile caches
    key on (see `queryDb.galaxy_content_stamp`)."""
    _require_db(db)
    return _request("/galaxy/stamp", {"db": db})["stamp"]

def get_phenomena(db, limit=None, offset=None):
    """Returns `GET /api/phenomena`'s full paginated envelope
    (`items`/`total`/`limit`/`offset`) -- see `queryDb.list_phenomena`'s
    docstring for the shape."""
    _require_db(db)
    return _request("/phenomena", {"db": db, "limit": limit, "offset": offset})


def get_phenomenon(db, phenomenon_type, phenomenon_id):
    """Returns `GET /api/phenomena/<type>/<id>`'s detail dict -- see
    `queryDb.phenomenon_detail`'s docstring for the shape."""
    _require_db(db)
    return _request(f"/phenomena/{phenomenon_type}/{phenomenon_id}", {"db": db})


def get_search(db, texts, tags, sizes=None, limit=None, offsets=None):
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
    return _request("/search", pairs)


# ---------------------------------------------------------------------
# Admin auth -- see `docs/api.md`'s "Authentication" section. Every
# function below takes/returns the incoming/outgoing `Cookie`/`Set-Cookie`
# headers verbatim (see `_auth_request`'s docstring) -- `login.py`/
# `changecreds.py`/`admin.py` are the only callers, and relay them to/from
# the browser via `page.py`'s own cookie helpers.
# ---------------------------------------------------------------------

def auth_login(username, password):
    """
    `POST /api/auth/login`.

    Returns:
        tuple[dict, list[str]]: `({"username", "must_change_credentials"},
            set_cookie_headers)` on success.

    Raises:
        ApiError: `status_code == 401` for a wrong username/password --
            callers should catch this specifically and show an inline
            "invalid username or password" message rather than a generic
            error page.
    """
    return _auth_request("POST", "/auth/login", json_body={"username": username, "password": password})


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
        dict: The new page's `{"id", "path", "title", "url"}`.

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


_NEIGHBORHOOD_GENERATION_TIMEOUT_SECONDS = 1800
"""float: `generate_sector_neighborhood` below can legitimately run for a
very long time -- its default 100 ly radius holds on the order of
2,000-3,000 candidate sector slots (confirmed by measurement, see that
route's own docstring), each generated one at a time, synchronously (see
`routes.py`'s own note on this route having no background job queue to
hand off to) -- unlike every other quick CRUD call this module makes,
where `_TIMEOUT_SECONDS` alone would make a real, still-working request
look like a failure. Even this generous a timeout may not be enough for a
genuinely dense/large region -- there's no fully solving that without a
real job queue, which this project doesn't have; a caller triggering this
against an unfamiliar/large radius should pass a smaller `radius_ly`
first."""


def generate_sector_neighborhood(cookie_header, sector_id, radius_ly=None):
    """
    `POST /api/sectors/<id>/generate-neighborhood` -- generates every
    not-yet-generated sector within `radius_ly` (`None` for the API's own
    default, 100 ly) of this already galaxy-placed sector. The admin-only
    "generate more sectors around this one" action on `sector.py`. Uses
    `_NEIGHBORHOOD_GENERATION_TIMEOUT_SECONDS` rather than this module's
    usual, much shorter timeout -- see that constant's own docstring.

    Returns:
        dict: `generated`/`already_existed`/`candidates` -- see
            `generate.generate_sector_neighborhood`'s own docstring.

    Raises:
        ApiError: `status_code == 404` if the sector doesn't exist or was
            never placed in a galaxy.
    """
    body, _set_cookie_headers = _auth_request(
        "POST", f"/sectors/{sector_id}/generate-neighborhood", cookie_header=cookie_header,
        json_body={"radius_ly": radius_ly} if radius_ly is not None else {},
        timeout=_NEIGHBORHOOD_GENERATION_TIMEOUT_SECONDS,
    )
    return body
