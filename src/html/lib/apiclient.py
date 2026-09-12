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

_TIMEOUT_SECONDS = 15


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


def _auth_request(method, path, json_body=None, cookie_header=None):
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
        with urllib.request.urlopen(request, timeout=_TIMEOUT_SECONDS) as response:
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
    page's own `?db=` query parameter, forwarded straight through.
    Omitting it wouldn't error (the API falls back to its own configured
    default database, see `api/routes.py`'s `get_db`), but every one of
    these pages embeds `db` into the links it renders (`sector.py?db=...`,
    `system.py?db=...`, ...), so a missing `db` here would quietly build a
    page entirely out of a different database than the one every link on
    it claims to be showing. Raised eagerly instead, matching the old
    direct-database `dbutil.resolve_db_name`'s own "No database
    specified." check.

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


def get_nav(db, from_id, to_id):
    """Returns `GET /api/nav`'s course/route dict -- see `docs/api.md`'s
    "NAV" section for the full shape."""
    _require_db(db)
    return _request("/nav", {"db": db, "from": from_id, "to": to_id})


def get_galaxy_sectors(db):
    """Returns `GET /api/galaxy/sectors`'s `items` list (every
    galaxy-placed sector)."""
    _require_db(db)
    return _request("/galaxy/sectors", {"db": db})["items"]


def get_search(db, texts, tags, sizes=None):
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
