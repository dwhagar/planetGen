# html/lib/apiclient.py

"""
JSON HTTP client for the planetGen Flask API (`src/api/`), used by every
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
import urllib.error
import urllib.parse
import urllib.request

API_BASE_URL = os.environ.get("PLANETGEN_API_BASE_URL", "http://127.0.0.1/api")
"""str: Base URL of the planetGen API's `/api` mount point. Defaults to
the same host this CGI script itself runs on (see
`examples/apache/planetgen.conf.example`'s `WSGIScriptAlias /api`) --
override via the `PLANETGEN_API_BASE_URL` Apache `SetEnv` (or shell env,
for local testing) when the API is deployed at a different host/port,
e.g. `http://127.0.0.1:5000/api` for `python src/wsgi.py`'s own dev
server running alongside a locally-invoked CGI script."""

_TIMEOUT_SECONDS = 15


class NotFoundError(Exception):
    """Raised for a 404 response (an invalid `db`, sector/system id, etc.)
    -- callers turn this into a 404 page, same as before this module
    existed (`html/lib/page.run` already handles it)."""


class ApiError(Exception):
    """Raised for any other failure talking to the API: an unreachable
    process, a non-2xx/404 response, or an unparseable body."""


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
        raise ApiError(f"planetGen API error ({exc.code}): {detail}")
    except urllib.error.URLError as exc:
        raise ApiError(f"Could not reach the planetGen API at {API_BASE_URL}: {exc.reason}")

    try:
        return json.loads(body)
    except ValueError as exc:
        raise ApiError(f"planetGen API returned an unparseable response: {exc}")


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


def get_search(db, texts, tags):
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
    """
    _require_db(db)
    pairs = [("db", db)]
    pairs.extend((key, value) for key, value in texts.items() if value)
    for facet, values in tags.items():
        pairs.extend((facet, value) for value in values)
    return _request("/search", pairs)
