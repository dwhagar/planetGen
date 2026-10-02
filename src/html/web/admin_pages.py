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
  cookie (`_flash`, scoped to `/admin`, the only page that reads it),
  never the URL. The one exception is a form that
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
from fmt import format_duration_seconds, format_number, utc_time_html
import tilecache
from pagination import fetch_page, page_slice, parse_page

from . import bp
from .helpers import crumb, current_admin, db_name, page_url, pager, render_page, trusted_html

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


def _flash_path():
    """The flash cookie's `Path`: `/admin`, the only page that reads it
    (`_take_flash`). Signed is not encrypted, and a flash can carry a new
    API key's raw value, so the browser must not send it with every
    request to the site (logs, other pages, a proxy's headers) -- only
    back to `/admin` (and the pages under it), which consumes it."""
    return url_for("web.admin")


def _flash_serializer():
    return URLSafeTimedSerializer(current_app.config["SECRET_KEY"], salt="planetgen-web-flash")


def _flash(response, **data):
    """
    Attaches `data` (a new API key, a message, an error) to a redirect
    for the next GET to show once. Signed with `SECRET_KEY` so it can't
    be forged, HttpOnly, SameSite=Strict, Secure like the session cookie,
    sent back only to `/admin` (`_flash_path`), and gone after five
    minutes or the next `/admin` view, whichever is first.
    """
    response.set_cookie(
        FLASH_COOKIE, _flash_serializer().dumps(data),
        max_age=_FLASH_MAX_AGE, httponly=True, samesite="Strict", path=_flash_path(),
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
        response.delete_cookie(FLASH_COOKIE, path=_flash_path(), samesite="Strict", httponly=True,
                               secure=current_app.config.get("SESSION_COOKIE_SECURE", True))
    return response


# ---------------------------------------------------------------------
# /login, /logout, /account
# ---------------------------------------------------------------------

def _login_page(next_url, error=None, username="", status=200, pending=None):
    return _render(
        "login.html", title="Admin Login", section="login", breadcrumbs=[crumb("Login")],
        description="Sign in to administer this planetGen site.",
        error=error, username=username, next_url=next_url, status=status, pending=pending,
    )


def _signed_in(result, set_cookie_headers, next_url):
    if result.get("must_change_credentials"):
        destination = url_for("web.account", next=next_url)
    else:
        destination = next_url
    return _relay(_see_other(destination), set_cookie_headers)


@bp.route("/login", methods=["GET", "POST"])
def login():
    """The login form. A successful login goes to `next` (or `/admin`),
    by way of `/account` while the default credentials are in use."""
    next_url = safe_next(request.values.get("next"))
    if request.method in ("GET", "HEAD"):  # HEAD is a GET without the body, never the POST branch
        admin = current_admin()
        if admin is not None and not admin.get("must_change_credentials"):
            return _no_store(redirect(next_url, code=302))
        return _login_page(next_url)

    username = request.form.get("username", "")
    password = request.form.get("password", "")
    if not username.strip() or not password:
        return _login_page(next_url, error="Enter a username and password.", username=username)
    try:
        result, set_cookie_headers = apiclient.auth_login(username, password, _cookie_header())
    except apiclient.ApiError as exc:
        if exc.status_code == 401:
            # 401, not 200, so the web server's access log shows a failed
            # login too (SEC.20).
            return _login_page(next_url, error="Invalid username or password.", username=username, status=401)
        if exc.status_code == 429:
            return _login_page(next_url, error=_too_many_message(exc), username=username, status=429)
        if exc.status_code == 400:
            return _login_page(next_url, error=_api_message(exc), username=username)
        raise
    if result.get("totp_required"):
        # Right password; now the authenticator code (SEC.26).
        return _login_page(next_url, username=username, pending=result["pending"])
    return _signed_in(result, set_cookie_headers, next_url)


@bp.route("/login/code", methods=["POST"])
def login_code():
    """The second step of a two-factor login: the authenticator (or
    recovery) code, with the signed `pending` value the first step gave."""
    next_url = safe_next(request.values.get("next"))
    pending = request.form.get("pending", "")
    code = request.form.get("code", "")
    username = request.form.get("username", "")
    if not pending:
        return _login_page(next_url, error="Enter your username and password first.", username=username)
    if not code.strip():
        return _login_page(next_url, error="Enter the code from your authenticator app.", username=username,
                           pending=pending)
    try:
        result, set_cookie_headers = apiclient.auth_login_totp(pending, code, _cookie_header())
    except apiclient.ApiError as exc:
        if exc.status_code == 401:
            message = _api_message(exc)
            if "password" in message:
                return _login_page(next_url, error=message[0].upper() + message[1:] + ".", username=username,
                                   status=401)
            return _login_page(next_url, error="That code isn't right. Try the newest one.", username=username,
                               pending=pending, status=401)
        if exc.status_code == 429:
            return _login_page(next_url, error=_too_many_message(exc), username=username, status=429)
        if exc.status_code == 400:
            return _login_page(next_url, error=_api_message(exc), username=username, pending=pending)
        raise
    return _signed_in(result, set_cookie_headers, next_url)


@bp.route("/logout", methods=["GET", "POST"])
def logout():
    """`GET` asks to confirm (it changes nothing); `POST` ends the
    session, clearing the cookie, and goes home."""
    if request.method in ("GET", "HEAD"):  # HEAD is a GET without the body, never the POST branch
        return _render(
            "logout.html", title="Log out", breadcrumbs=[crumb("Log out")], admin=current_admin(),
        )
    set_cookie_headers = apiclient.auth_logout(_cookie_header())
    return _relay(_see_other(url_for("web.index")), set_cookie_headers)


def _too_many_message(exc):
    """The message for a 429 from a password check: a lockout (per address
    or per username, `api/loginguard.py`) names its wait; the per-address
    rate limit doesn't."""
    message = _api_message(exc)
    if "try again in" in message:
        return message[0].upper() + message[1:].replace("; try", ". Try") + "."
    return "Too many login attempts. Wait a minute, then try again."


def _account_page(admin, next_url, error=None, new_username=None, status=200, totp=None):
    """`totp`: what the two-factor section shows -- `error`, `message`,
    `setup` (the API's setup response), `recovery_codes`."""
    totp = dict(totp or {})
    if not admin.get("must_change_credentials"):
        try:
            totp.setdefault("status", apiclient.auth_totp_status(_cookie_header()))
        except apiclient.ApiError:
            totp["status"] = None  # e.g. before update.sh made admin_totp
    if totp.get("setup", {}).get("qr_svg"):
        totp["qr"] = trusted_html(totp["setup"]["qr_svg"])
    return _render(
        "account.html", title="Change Credentials", section="account",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Account")],
        admin=admin, error=error, next_url=next_url,
        new_username=admin["username"] if new_username is None else new_username,
        forced=bool(admin.get("must_change_credentials")), status=status, totp=totp,
    )


@bp.route("/account/two-factor", methods=["POST"])
def account_two_factor():
    """Turn two-factor sign-in on (set up, then confirm with a code) or
    off (SEC.26)."""
    admin, bounce = _require_admin(fresh=True)
    if bounce is not None:
        return bounce
    next_url = safe_next(request.values.get("next"))
    action = request.form.get("action")
    cookies = _cookie_header()
    try:
        if action == "setup":
            setup = apiclient.auth_totp_setup(cookies, request.form.get("current_password", ""))
            return _account_page(admin, next_url, totp={"setup": setup})
        if action == "confirm":
            codes = apiclient.auth_totp_confirm(cookies, request.form.get("code", ""))["recovery_codes"]
            return _account_page(admin, next_url, totp={"recovery_codes": codes,
                                                        "message": "Two-factor sign-in is on."})
        if action == "disable":
            set_cookie_headers = apiclient.auth_totp_disable(cookies, request.form.get("current_password", ""),
                                                             request.form.get("code", ""))
            return _relay(_account_page(admin, next_url, totp={"message": "Two-factor sign-in is off."}),
                          set_cookie_headers)
    except apiclient.ApiError as exc:
        if exc.status_code == 429:
            return _account_page(admin, next_url, totp={"error": _too_many_message(exc)}, status=429)
        if exc.status_code in (400, 401):
            message = _api_message(exc)
            totp = {"error": message[0].upper() + message[1:] + "."}
            if action == "confirm":
                # Show the same QR code again (the pending secret is
                # unchanged); the form carried it back.
                secret = request.form.get("secret", "")
                if secret:
                    from stellarObjects import totp as totp_codes
                    uri = totp_codes.provisioning_uri(secret, admin["username"])
                    totp["setup"] = {"secret": secret, "uri": uri, "qr_svg": totp_codes.qr_svg(uri)}
            return _account_page(admin, next_url, totp=totp, status=400)
        raise
    return _account_page(admin, next_url, totp={"error": "Unknown action."}, status=400)


@bp.route("/account", methods=["GET", "POST"])
def account():
    """Change the logged-in admin's username and password (forced while
    the seeded defaults are in use, voluntary otherwise)."""
    admin, bounce = _require_admin(fresh=False)
    if bounce is not None:
        return bounce
    next_url = safe_next(request.values.get("next"))
    if request.method in ("GET", "HEAD"):  # HEAD is a GET without the body, never the POST branch
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
        if exc.status_code == 429:
            # A wrong current password counts as a failed login (SEC.23).
            return _account_page(admin, next_url, error=_too_many_message(exc), new_username=new_username,
                                 status=429)
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
        "last_used_at": key.get("last_used_at"),
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
    """An uptime or wait on the shared period ladder
    (`format_duration_seconds`, UX.14): "45 s", "12 minutes", "3.5 days";
    "unknown" for `None`."""
    if seconds is None:
        return "unknown"
    return format_duration_seconds(seconds)


def format_count(value):
    return "unknown" if value is None else format_number(value)


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
            ("Python", f"{api['python_version']} ({api['python_prefix']})"
             if api.get("python_prefix") else api["python_version"]),
            ("Libraries from", api.get("libraries_dir") or "unknown"),
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


def bright_star_text(database):
    """
    The Stats page's bright-star row: how many stars the plan's scatter
    pre-placed and how many of them a sector fill has built into a system
    (`database["bright_stars"]`, `adminStats.bright_star_counts`). Falls
    back to the `bright_stars` table's estimated row count when the API
    doesn't send the counts.
    """
    counts = database.get("bright_stars")
    if counts:
        if not counts["placed"]:
            return "none pre-placed"
        return (f"{format_count(counts['placed'])} placed: {format_count(counts['filled'])} built into "
                f"systems, {format_count(counts['unfilled'])} waiting for their sectors")
    row = next((table for table in database["tables"] if table["name"] == "bright_stars"), None)
    if row is None:
        return "not tracked by this schema"
    if not row["approx_rows"]:
        return "none pre-placed"
    return f"about {format_count(row['approx_rows'])} placed in all (estimate)"


def density_text(database):
    """
    The Stats page's density row (PERF.11): the decaying average of the
    systems sector fills got against the systems the galaxy model
    expected, and how many sectors are measured or backfilled
    (`database["sector_stats"]`, `adminStats.density_stats`).
    """
    stats = database.get("sector_stats")
    if not stats:
        return "not tracked by this schema"
    backfilled = f"{format_count(stats['backfilled'])} sectors backfilled to their own level"
    if stats["ratio"] is None:
        return f"no sector filled yet; {backfilled}"
    return (f"{stats['ratio']:.2f} of the expected systems (average of {format_count(stats['fills'])} fills); "
            f"{format_count(stats['measured'])} sectors measured, {backfilled}")


def _time_or(value, missing):
    """A stats time as local-time markup, or `missing` when there is none."""
    return trusted_html(utc_time_html(value)) if value else missing


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
}


def _name_entry(row):
    """One "now named" row: its name, what it is, and a link."""
    if row["kind"] == "sector":
        return {"name": row["name"], "url": page_url("sector", sector_id=row["id"]), "kind": "sector"}
    return {"name": row["name"], "url": page_url("system", system_id=row["id"]), "kind": "system"}


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


_FAILURE_LABELS = {
    "login.failed": "Wrong username or password",
    "login.locked": "Refused while locked",
    "password.failed": "Wrong current password (account page)",
    "totp.failed": "Wrong two-factor code",
}


def _failures_panel(cookie_header):
    """The newest refused sign-ins (SEC.20), for the stats page."""
    try:
        items = apiclient.admin_login_failures(cookie_header)
    except apiclient.ApiError as exc:
        return {"error": _api_message(exc), "items": []}
    return {"error": None, "items": [{
        "what": _FAILURE_LABELS.get(item["action"], item["action"]),
        "username": item["username"],
        "ip": item["ip"] or "unknown",
        "created_at": item["created_at"],
    } for item in items]}


def _lockouts_panel(cookie_header):
    """Every address and username locked right now (SEC.1, SEC.21)."""
    try:
        body = apiclient.admin_lockouts(cookie_header)
    except apiclient.ApiError as exc:
        return {"error": _api_message(exc), "items": [], "proxy_warning": False}
    return {"error": None, "proxy_warning": body.get("proxy_warning", False), "items": [{
        "scope": item["scope"],
        "what": "Address" if item["scope"] == "ip" else "Username",
        "subject": item["subject"],
        "wait": format_duration(item["retry_after"]),
        "level": item["level"],
    } for item in body["items"]]}


def _generation_panel(cookie_header, db):
    """PERF.10: this server's measured generation speed per density
    bucket, and this galaxy's size per star system."""
    try:
        body = apiclient.admin_generation_stats(cookie_header)
    except apiclient.ApiError as exc:
        return {"error": _api_message(exc), "rows": [], "size": None, "available": False}
    size = body["sizes"].get(db)
    return {
        "error": None,
        "available": body["available"],
        "size": None if size is None else (
            f"{format_bytes(size['bytes_per_system'])} per star system "
            f"({format_count(size['systems'])} systems, {format_bytes(size['total_bytes'])})"
        ),
        "rows": [{
            "what": "Sector fill" if row["kind"] == "sector" else "Bright-star layer",
            "density": f"{row['density_low']:.3g} to {row['density_high']:.3g}",
            "samples": format_count(row["samples"]),
            "per_task": f"{row['seconds_per_task']:.2f} s",
            "per_system": f"{row['seconds_per_system'] * 1000:.1f} ms",
            "systems": f"{row['systems_per_task']:.1f}",
        } for row in body["buckets"]],
    }


@bp.route("/admin/stats/lockouts", methods=["POST"])
def lift_lockout():
    """The Stats page's "Lift" buttons: one lockout, or all of them."""
    _identity, bounce = _require_admin()
    if bounce is not None:
        return bounce
    try:
        apiclient.admin_lift_lockout(
            _cookie_header(), scope=request.form.get("scope"), subject=request.form.get("subject"),
            lift_all=request.form.get("all") == "1",
        )
    except apiclient.ApiError as exc:
        if exc.status_code is None or exc.status_code >= 500:
            raise
    return _see_other(url_for("web.admin_stats", _anchor="lockouts"))


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

    context = {"health": _health(stats, api_ms, tile_cache_info()), "database": database,
               "failures": _failures_panel(cookie_header), "lockouts": _lockouts_panel(cookie_header),
               "generation": _generation_panel(cookie_header, db)}
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
                ("Last system change", _time_or(systems_stamp.get("last_modified_at"), "never")),
                ("Newest system", _time_or(systems_stamp.get("newest_created_at"), "none")),
                ("Last sector change", _time_or(stamps.get("sectors", {}).get("last_modified_at"), "never")),
                ("Bright stars", bright_star_text(database)),
                ("Sector density", density_text(database)),
            ],
            name_tiles=[
                ("Names made unique", format_count(collisions.get("distinct_base_names"))),
                ("Sector collisions", format_count(collisions.get("sector"))),
                ("System collisions", format_count(collisions.get("system"))),
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
