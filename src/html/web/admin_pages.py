# html/web/admin_pages.py

"""
The admin pages: `/login`, `/logout`, `/account` (change username/
password), `/admin` (API keys and a sector's manual wiki link) and
`/admin/stats` (server health and database statistics). They replace
`login.py`, `logout.py`, `changecreds.py`, `admin.py` and
`adminstats.py`.

Sessions work exactly as they did on CGI: the views call the same
`apiclient.auth_*` functions (in-process, see `transport.py`), and the
API's own `/api/auth/*` routes decide the session cookie. Whatever
`Set-Cookie` headers they send back are relayed onto the page's response
verbatim (`_relay`), so the cookie keeps the API's HttpOnly/Secure/
SameSite=Strict attributes and the login route keeps its own per-address
rate limit (`api.auth.LOGIN_RATE_LIMIT`; the transport passes the
visitor's address through).

Rules every view here follows:

- Every form that changes something is a POST carrying `csrf_field()`
  (checked app-wide by `csrf.protect`), answered with a 303 redirect to a
  GET page (POST-redirect-GET). What the POST did (a new API key, a
  message, an error) reaches that GET page through a one-shot signed
  cookie (`_flash`), never the URL. The one exception is a form that
  must be shown again with an inline error and the visitor's input
  (a wrong password on `/login` or `/account`): that answers the POST
  directly, and nothing changed.
- `GET /logout` changes nothing: it shows a button that POSTs.
- A visitor who isn't logged in goes to `/login?next=<this page>`;
  `next` is only ever followed when it is a local path (`safe_next`).
  An admin still on the seeded default credentials goes to `/account`
  first.
- Responses are `Cache-Control: no-store` (a new API key, admin data).
"""

import os
import re
import shutil
import time
from urllib.parse import urlsplit

from flask import current_app, make_response, redirect, request, url_for
from itsdangerous import BadSignature, URLSafeTimedSerializer

import apiclient
import tilecache
from pagination import fetch_page, page_slice, parse_page

from . import bp
from .helpers import crumb, current_admin, db_name, page_url, pager, render_page

# ---------------------------------------------------------------------
# Redirect targets
# ---------------------------------------------------------------------

_CONTROL_CHARS = re.compile(r"[\x00-\x20\x7f\\]")


def safe_next(value, default=None):
    """
    `value` if it is a local path on this site (`/admin/stats?names_page=2`),
    else `default` (the admin page when not given). Rejects anything a
    browser could resolve to another origin: absolute URLs, scheme-relative
    `//host`, backslashes (`/\\host`), whitespace and control characters.
    """
    if default is None:
        default = page_url("admin")
    if not value or not isinstance(value, str) or len(value) > 2000:
        return default
    if not value.startswith("/") or value.startswith("//") or _CONTROL_CHARS.search(value):
        return default
    parts = urlsplit(value)
    if parts.scheme or parts.netloc:
        return default
    return value


def _here():
    """This request's path and query, for a `?next=` back to it."""
    query = request.query_string.decode("utf-8", errors="replace")
    return f"{request.path}?{query}" if query else request.path


def _see_other(location):
    """A 303 (the POST-redirect-GET answer) to a local `location`."""
    return _no_store(redirect(location, code=303))


def _login_redirect():
    return _no_store(redirect(url_for("web.login", next=_here()), code=302))


def _no_store(response):
    response = make_response(response)
    response.headers["Cache-Control"] = "no-store"
    return response


def _relay(response, set_cookie_headers):
    """Adds the API's `Set-Cookie` headers to `response` unchanged."""
    for header in set_cookie_headers or ():
        response.headers.add("Set-Cookie", header)
    return response


def _cookie_header():
    return request.headers.get("Cookie")


def _require_admin(fresh=True):
    """
    The logged-in admin, or a redirect response: to `/login?next=...`
    when nobody is logged in, and (with `fresh`) to `/account?next=...`
    while the admin is still on the seeded default credentials.

    Returns:
        tuple: `(admin, None)` or `(None, response)`.
    """
    admin = current_admin()
    if admin is None:
        return None, _login_redirect()
    if fresh and admin.get("must_change_credentials"):
        return None, _no_store(redirect(url_for("web.account", next=_here()), code=302))
    return admin, None


_API_ERROR_PREFIX = re.compile(r"^planetGen API error \(\d+\): ")


def _api_message(exc):
    """The API's own message from an `ApiError`, without the
    `planetGen API error (400): ` prefix."""
    return _API_ERROR_PREFIX.sub("", str(exc))


# ---------------------------------------------------------------------
# One-shot messages across a redirect
# ---------------------------------------------------------------------

FLASH_COOKIE = "pg_flash"
_FLASH_MAX_AGE = 300


def _flash_serializer():
    return URLSafeTimedSerializer(current_app.config["SECRET_KEY"], salt="planetgen-web-flash")


def _flash(response, **data):
    """
    Attaches `data` (a new API key, a message, an error) to a redirect
    for the next GET to show once. Signed with `SECRET_KEY` so it can't
    be forged, HttpOnly, SameSite=Strict, Secure like the session cookie,
    and gone after five minutes or the next page view, whichever is
    first.
    """
    response.set_cookie(
        FLASH_COOKIE, _flash_serializer().dumps(data),
        max_age=_FLASH_MAX_AGE, httponly=True, samesite="Strict", path="/",
        secure=current_app.config.get("SESSION_COOKIE_SECURE", True),
    )
    return response


def _take_flash():
    """This request's flashed data (`{}` for none or a bad signature)."""
    raw = request.cookies.get(FLASH_COOKIE)
    if not raw:
        return {}
    try:
        data = _flash_serializer().loads(raw, max_age=_FLASH_MAX_AGE)
    except BadSignature:
        return {}
    return data if isinstance(data, dict) else {}


def _render(template, flashed=False, **kwargs):
    """`render_page` plus `no-store`, clearing the flash cookie if this
    page consumed one."""
    response = _no_store(render_page(template, **kwargs))
    if flashed and request.cookies.get(FLASH_COOKIE):
        response.delete_cookie(FLASH_COOKIE, path="/", samesite="Strict", httponly=True,
                               secure=current_app.config.get("SESSION_COOKIE_SECURE", True))
    return response


# ---------------------------------------------------------------------
# /login, /logout, /account
# ---------------------------------------------------------------------

def _login_page(next_url, error=None, username="", status=200):
    return _render(
        "login.html", title="Admin Login", section="login", breadcrumbs=[crumb("Login")],
        description="Sign in to administer this planetGen site.",
        error=error, username=username, next_url=next_url, status=status,
    )


@bp.route("/login", methods=["GET", "POST"])
def login():
    """The login form. A successful login goes to `next` (or `/admin`),
    by way of `/account` while the default credentials are in use."""
    next_url = safe_next(request.values.get("next"))
    if request.method == "GET":
        admin = current_admin()
        if admin is not None and not admin.get("must_change_credentials"):
            return _no_store(redirect(next_url, code=302))
        return _login_page(next_url)

    username = request.form.get("username", "")
    password = request.form.get("password", "")
    if not username.strip() or not password:
        return _login_page(next_url, error="Enter a username and password.", username=username)
    try:
        result, set_cookie_headers = apiclient.auth_login(username, password)
    except apiclient.ApiError as exc:
        if exc.status_code == 401:
            return _login_page(next_url, error="Invalid username or password.", username=username)
        if exc.status_code == 429:
            return _login_page(next_url, error="Too many login attempts. Wait a minute, then try again.",
                               username=username, status=429)
        if exc.status_code == 400:
            return _login_page(next_url, error=_api_message(exc), username=username)
        raise
    if result.get("must_change_credentials"):
        destination = url_for("web.account", next=next_url)
    else:
        destination = next_url
    return _relay(_see_other(destination), set_cookie_headers)


@bp.route("/logout", methods=["GET", "POST"])
def logout():
    """`GET` asks to confirm (it changes nothing); `POST` ends the
    session, clearing the cookie, and goes home."""
    if request.method == "GET":
        return _render(
            "logout.html", title="Log out", breadcrumbs=[crumb("Log out")], admin=current_admin(),
        )
    set_cookie_headers = apiclient.auth_logout(_cookie_header())
    return _relay(_see_other(url_for("web.index")), set_cookie_headers)


def _account_page(admin, next_url, error=None, new_username=None):
    return _render(
        "account.html", title="Change Credentials", section="account",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Account")],
        admin=admin, error=error, next_url=next_url,
        new_username=admin["username"] if new_username is None else new_username,
        forced=bool(admin.get("must_change_credentials")),
    )


@bp.route("/account", methods=["GET", "POST"])
def account():
    """Change the logged-in admin's username and password (forced while
    the seeded defaults are in use, voluntary otherwise)."""
    admin, bounce = _require_admin(fresh=False)
    if bounce is not None:
        return bounce
    next_url = safe_next(request.values.get("next"))
    if request.method == "GET":
        return _account_page(admin, next_url)

    new_username = request.form.get("new_username", "")
    try:
        _result, set_cookie_headers = apiclient.auth_change_credentials(
            _cookie_header(),
            current_password=request.form.get("current_password", ""),
            new_username=new_username,
            new_password=request.form.get("new_password", ""),
        )
    except apiclient.ApiError as exc:
        if exc.status_code in (400, 401):
            return _account_page(admin, next_url, error=_api_message(exc), new_username=new_username)
        raise
    response = _see_other(next_url)
    if next_url.split("?", 1)[0] == url_for("web.admin"):
        _flash(response, message="Your username and password were changed.")
    return _relay(response, set_cookie_headers)


# ---------------------------------------------------------------------
# /admin
# ---------------------------------------------------------------------

def _key_rows(keys):
    return [{
        "id": key["id"],
        "label": key["label"],
        "created_at": key.get("created_at") or "",
        "last_used_at": key.get("last_used_at") or "never",
        "revoked_at": key.get("revoked_at"),
    } for key in keys]


def _admin_action(cookie_header):
    """
    Runs one `POST /admin` form. Returns `(flash data, anchor)`.
    """
    action = request.form.get("action")
    try:
        if action == "create_key":
            label = request.form.get("label", "").strip()
            if not label:
                return {"error": "Label is required."}, "api-keys"
            new_key = apiclient.auth_create_api_key(cookie_header, label)
            return {"new_key": {"label": new_key["label"], "key": new_key["key"]}}, "new-key"
        if action == "revoke_key":
            try:
                key_id = int(request.form.get("key_id", ""))
            except ValueError:
                return {"error": "Invalid key id."}, "api-keys"
            apiclient.auth_revoke_api_key(cookie_header, key_id)
            return {"message": "API key revoked."}, "api-keys"
        if action == "set_sector_wiki_url":
            return _set_wiki_url(cookie_header), "sector-wiki-link"
    except apiclient.NotFoundError as exc:
        key = "wiki_error" if action == "set_sector_wiki_url" else "error"
        return {key: str(exc) or "Not found."}, "sector-wiki-link" if key == "wiki_error" else "api-keys"
    except apiclient.ApiError as exc:
        if exc.status_code is None or exc.status_code >= 500:
            raise
        key = "wiki_error" if action == "set_sector_wiki_url" else "error"
        return {key: _api_message(exc)}, "sector-wiki-link" if key == "wiki_error" else "api-keys"
    return {"error": "Unrecognized form action."}, None


def _set_wiki_url(cookie_header):
    sector_id_raw = request.form.get("sector_id", "").strip()
    wiki_url = request.form.get("wiki_url", "").strip() or None
    if not sector_id_raw:
        return {"wiki_error": "Sector ID is required."}
    try:
        sector_id = int(sector_id_raw)
    except ValueError:
        return {"wiki_error": "Invalid sector id."}
    apiclient.admin_set_sector_wiki_url(cookie_header, db_name(), sector_id, wiki_url)
    done = "cleared" if wiki_url is None else f"set to {wiki_url}"
    return {"wiki_message": f"Wiki link for sector {sector_id} {done}."}


@bp.route("/admin", methods=["GET", "POST"])
def admin():
    """API keys (list, create, revoke; `?keys_page=N`) and a sector's
    manual wiki link."""
    identity, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    keys_page = parse_page(request.args.get("keys_page"))

    if request.method == "POST":
        flash, anchor = _admin_action(cookie_header)
        keys_page = parse_page(request.form.get("keys_page"))
        return _flash(_see_other(page_url("admin", keys_page=keys_page if keys_page > 1 else None,
                                          _anchor=anchor)), **flash)

    keys = apiclient.auth_list_api_keys(cookie_header)
    page_keys, keys_page = page_slice(keys, keys_page)
    flashed = _take_flash()
    return _render(
        "admin.html", flashed=True, title="Admin", section="admin", breadcrumbs=[crumb("Admin")],
        admin=identity, keys=_key_rows(page_keys), keys_total=len(keys), keys_page=keys_page,
        keys_pager=pager("keys_page", keys_page, len(keys), anchor="api-keys", label="API key pages"),
        database=db_name(), new_key=flashed.get("new_key"), message=flashed.get("message"),
        error=flashed.get("error"), wiki_message=flashed.get("wiki_message"),
        wiki_error=flashed.get("wiki_error"),
    )


# ---------------------------------------------------------------------
# /admin/stats
# ---------------------------------------------------------------------

def format_bytes(value):
    if value is None:
        return "unknown"
    size = float(value)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024


def format_duration(seconds):
    if seconds is None:
        return "unknown"
    seconds = int(seconds)
    days, seconds = divmod(seconds, 86400)
    hours, seconds = divmod(seconds, 3600)
    minutes = seconds // 60
    if days:
        return f"{days}d {hours}h"
    if hours:
        return f"{hours}h {minutes}m"
    return f"{minutes}m"


def format_count(value):
    return "unknown" if value is None else f"{value:,}"


def tile_cache_info():
    """
    Where the galaxy tile cache lives, how much it holds against its
    budget, and how much room is left on that disk. Reads the directory
    without creating it (unlike `tilecache.cache_dir`).
    """
    configured = tilecache.configured_cache_dir()
    info = {
        "enabled": configured is not None,
        "dir": configured,
        "exists": False,
        "size_bytes": 0,
        "files": 0,
        "max_bytes": tilecache.max_cache_bytes(),
        "disk_free_bytes": None,
    }
    if configured is None or not os.path.isdir(configured):
        return info
    info["exists"] = True
    for dirpath, _dirnames, filenames in os.walk(configured):
        for name in filenames:
            try:
                info["size_bytes"] += os.stat(os.path.join(dirpath, name)).st_size
                info["files"] += 1
            except OSError:
                continue
    try:
        info["disk_free_bytes"] = shutil.disk_usage(configured).free
    except OSError:
        pass
    return info


def _health(stats, api_ms, cache):
    api = stats["api"]
    mysql = stats.get("mysql") or {}
    database = stats["database"]
    if database["reachable"] and database.get("schema_current"):
        status = "healthy"
    elif database["reachable"]:
        status = "schema out of date"
    else:
        status = "database unreachable"
    memory = api.get("memory")
    load = api.get("load_average")
    if not cache["enabled"]:
        cache_text, cache_dir, cache_tail = "disabled (max size is 0)", None, ""
    elif not cache["exists"]:
        cache_text, cache_dir, cache_tail = "", cache["dir"], " doesn't exist yet"
    else:
        free = format_bytes(cache["disk_free_bytes"]) if cache["disk_free_bytes"] is not None else "unknown"
        cache_text = (f"{format_bytes(cache['size_bytes'])} of {format_bytes(cache['max_bytes'])} "
                      f"({format_count(cache['files'])} files) in ")
        cache_dir, cache_tail = cache["dir"], f", {free} free on disk"
    return {
        "status": status,
        "rows": [
            ("API response", f"{api_ms:.0f} ms (queries took {stats.get('query_ms', 0):.0f} ms)"),
            ("API version", api["version"]),
            ("API uptime", format_duration(api["uptime_seconds"])),
            ("Python", api["python_version"]),
            ("Load average (1/5/15 min)", " / ".join(f"{value:.2f}" for value in load) if load else "unknown"),
            ("Memory", f"{format_bytes(memory['available_bytes'])} free of {format_bytes(memory['total_bytes'])}"
             if memory else "unknown"),
            ("MySQL version", mysql.get("version") or "unknown"),
            ("MySQL uptime", format_duration(mysql.get("uptime_seconds"))),
            ("MySQL connections", format_count(mysql.get("threads_connected"))),
        ],
        "cache_text": cache_text,
        "cache_dir": cache_dir,
        "cache_tail": cache_tail,
        "database_error": None if database["reachable"] else database.get("detail", ""),
    }


def _schema_text(database):
    version = database["schema_version"]
    if database["schema_current"]:
        return f"v{version} (current)"
    if version is None:
        return f"not initialized (code expects v{database['schema_expected']}); run migrateDb.py"
    return f"v{version}, code expects v{database['schema_expected']}; run migrateDb.py"


_LEVEL_LABELS = {
    "sector": "sector vs sector",
    "system": "system vs system or sector",
    "body": "planet or moon",
}


def _name_entry(row):
    """One "now named" row: its name, what it is, and a link (planets
    and moons link to their system)."""
    if row["kind"] == "sector":
        return {"name": row["name"], "url": page_url("sector", sector_id=row["id"]), "kind": "sector"}
    if row["kind"] == "system":
        return {"name": row["name"], "url": page_url("system", system_id=row["id"]), "kind": "system"}
    return {"name": row["name"], "kind": f"{row['kind']} in", "system_name": row["system_name"],
            "system_url": page_url("system", system_id=row["star_system_id"])}


def _names_panel(cookie_header, db):
    """One page of the duplicate-names list (`?names_page=N`)."""
    try:
        names, names_page = fetch_page(
            lambda limit, offset: apiclient.admin_duplicate_names(cookie_header, db, limit=limit, offset=offset),
            parse_page(request.args.get("names_page")),
        )
    except apiclient.ApiError as exc:
        return {"error": _api_message(exc), "items": [], "pager": ""}
    return {
        "error": None,
        "items": [{
            "base_name": item["base_name"],
            "levels": ", ".join(_LEVEL_LABELS.get(level, level) for level in item["levels"]),
            "rows": [_name_entry(row) for row in item["rows"]],
        } for item in names["items"]],
        "pager": pager("names_page", names_page, names["total"], anchor="duplicate-names",
                       label="Duplicate name pages"),
    }


@bp.route("/admin/stats")
def admin_stats():
    """Server health and statistics about the site's database, plus
    every name the uniqueness rules had to decorate (`?names_page=N`)."""
    _identity, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    db = db_name()
    started = time.perf_counter()
    stats = apiclient.admin_stats(cookie_header, db)
    api_ms = (time.perf_counter() - started) * 1000
    database = stats["database"]

    context = {"health": _health(stats, api_ms, tile_cache_info()), "database": database}
    if database["reachable"]:
        counts = database["counts"]
        collisions = database["name_collisions"]
        stamps = {row["table"]: row for row in database["timestamps"]}
        systems_stamp = stamps.get("star_systems", {})
        context.update(
            db_tiles=[
                ("Sectors", format_count(counts.get("sectors"))),
                ("Star systems", format_count(counts.get("star_systems"))),
                ("Size on disk", format_bytes(database["size_bytes"])),
                ("Names made unique", format_count(collisions.get("distinct_base_names"))),
            ],
            db_rows=[
                ("Schema", _schema_text(database)),
                ("Last system change", systems_stamp.get("last_modified_at") or "never"),
                ("Newest system", systems_stamp.get("newest_created_at") or "none"),
                ("Last sector change", stamps.get("sectors", {}).get("last_modified_at") or "never"),
            ],
            name_tiles=[
                ("Names made unique", format_count(collisions.get("distinct_base_names"))),
                ("Sector collisions", format_count(collisions.get("sector"))),
                ("System collisions", format_count(collisions.get("system"))),
                ("Planet/moon collisions", format_count(collisions.get("body"))),
            ],
            names=_names_panel(cookie_header, db),
            tables=[{
                "name": table["name"],
                "rows": format_count(table["approx_rows"]),
                "data": format_bytes(table["data_bytes"]),
                "indexes": format_bytes(table["index_bytes"]),
            } for table in database["tables"]],
        )
    return _render(
        "admin_stats.html", title="Server & Database Stats", section="admin_stats",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Stats")], **context,
    )
