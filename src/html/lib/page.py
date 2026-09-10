# html/lib/page.py

"""
Tiny CGI response/layout helpers shared by every script in `html/`.

No templating engine -- these are plain scripts meant to run with nothing
beyond the standard library, so pages are built as f-strings against one
shared shell (`render`) and a common stylesheet (`static/style.css`,
served directly by Apache, not through CGI).
"""

import os
import sys
import traceback
from urllib.parse import parse_qs

from apiclient import ApiError, NotFoundError, auth_me
from fmt import esc


def query_params():
    """
    Parses the CGI `QUERY_STRING` environment variable.

    Returns:
        dict[str, str]: One value per key (the first, if a key repeats) --
                        every script here only ever needs single-valued
                        params (`db`, `id`, `format`).
    """
    raw = parse_qs(os.environ.get("QUERY_STRING", ""))
    return {key: values[0] for key, values in raw.items()}


def form_params():
    """
    Parses an `application/x-www-form-urlencoded` POST body from stdin
    (`CONTENT_LENGTH` bytes, per the CGI spec) -- the admin pages'
    (`login.py`/`changecreds.py`/`admin.py`) login/credential/API-key
    forms, the first POST handling anywhere in `html/` (every other page
    here is GET-only). Reads at most `CONTENT_LENGTH` bytes, never the
    whole of stdin, so a misbehaving/absent header can't make this block
    waiting for input that never arrives.

    Returns:
        dict[str, str]: One value per field (the first, if a field
                        repeats) -- same single-valued convention as
                        `query_params`.
    """
    try:
        length = int(os.environ.get("CONTENT_LENGTH", "0") or "0")
    except ValueError:
        length = 0
    raw = sys.stdin.buffer.read(length).decode("utf-8", errors="replace") if length > 0 else ""
    parsed = parse_qs(raw)
    return {key: values[0] for key, values in parsed.items()}


def incoming_cookie_header():
    """
    Returns the browser's own `Cookie` header for this request (the CGI
    `HTTP_COOKIE` environment variable), unparsed -- forwarded verbatim to
    the planetGen API by every `apiclient.auth_*` call that needs the
    session cookie. This module never inspects or parses cookie values
    itself, only relays them (see `apiclient._auth_request`'s docstring).

    Returns:
        str or None: The raw `Cookie` header value, or `None` if the
                     browser sent none.
    """
    return os.environ.get("HTTP_COOKIE")


def _write_security_headers():
    """
    Written by both `send_headers` and `redirect` -- `'self'` in the CSP
    is safe for every page this browser renders: `render()`'s shared shell
    (see below) only ever loads `static/style.css` and each page's own
    `static/*.js`, both same-origin, and no page here builds an inline
    `<script>`/`<style>` block from request- or database-derived content.
    """
    sys.stdout.write("X-Content-Type-Options: nosniff\r\n")
    sys.stdout.write("X-Frame-Options: DENY\r\n")
    sys.stdout.write("Referrer-Policy: no-referrer\r\n")
    sys.stdout.write("Content-Security-Policy: default-src 'self'\r\n")


def send_headers(status="200 OK", set_cookie_headers=None):
    """
    Writes the CGI response status + header block, terminated by the
    required blank line, then flushes -- must be called exactly once,
    before any body output.

    Args:
        status (str): e.g. `"200 OK"`, `"404 Not Found"`.
        set_cookie_headers (list[str], optional): Raw `Set-Cookie` header
            values to relay verbatim (from the planetGen API's own login/
            logout/change-credentials response -- see
            `apiclient._auth_request`) -- this module never constructs a
            cookie value itself, only relays what the API already set
            (`HttpOnly`/`Secure`/`SameSite` included).
    """
    # The response declares charset=utf-8 below; stdout's default encoding
    # is platform/locale-dependent (e.g. cp1252 on Windows) and would
    # otherwise raise UnicodeEncodeError on any non-ASCII generated name.
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stdout.write(f"Status: {status}\r\n")
    sys.stdout.write("Content-Type: text/html; charset=utf-8\r\n")
    _write_security_headers()
    for cookie in set_cookie_headers or []:
        sys.stdout.write(f"Set-Cookie: {cookie}\r\n")
    sys.stdout.write("\r\n")


def redirect(url, set_cookie_headers=None):
    """
    Writes a CGI 302 redirect response to stdout and nothing else -- must
    be called instead of (never alongside) `send_headers`/`render`, and
    before any other output, same as `send_headers`'s own requirement.

    Used by `index.py` to skip straight to `browse.py?db=...` when exactly
    one database exists, instead of rendering the picker table for a
    choice of one; also by `login.py`/`changecreds.py` to redirect after a
    successful POST, carrying the API's `Set-Cookie` response along (see
    `send_headers`'s own parameter of the same name).

    Args:
        url (str): The target URL, e.g. `"browse.py?db=planetgen.db"`.
                   Callers are responsible for URL-encoding any dynamic
                   piece of it themselves (see `urllib.parse.quote`) --
                   this just writes the header block, the same division of
                   responsibility `send_headers` already has for the
                   `Content-Type` header.
        set_cookie_headers (list[str], optional): See `send_headers`.
    """
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stdout.write("Status: 302 Found\r\n")
    sys.stdout.write(f"Location: {url}\r\n")
    _write_security_headers()
    for cookie in set_cookie_headers or []:
        sys.stdout.write(f"Set-Cookie: {cookie}\r\n")
    sys.stdout.write("\r\n")


def _sidenav_html():
    """
    Builds the site-wide vertical nav bar shown on the left of every page:
    a fixed set of "functions" rather than a breadcrumb -- Databases, then
    Search, Galaxy, Sectors, and Systems for the current `?db=` when one is
    present in the request's own query string. Sectors/Systems jump
    straight to the matching anchor on `browse.py` (`id="sectors"`/
    `id="standalone-systems"`) rather than duplicating that page's own
    listing here; Galaxy goes to `galaxy.py`, the galaxy-scale map of every
    sector that actually has a galaxy position (see that page's own
    docstring) -- a different view than `browse.py`'s flat sector list.

    The Databases item always links to `index.py?all=1` (never bare
    `index.py`) whenever at least one `.db` file exists -- `index.py` on its
    own redirects straight past the picker when there's exactly one
    database (see its own module docstring), which would otherwise make
    this link a silent no-op/bounce-back on the single-database deployments
    this project actually runs in production, leaving no way back to the
    database-info page at all. `all=1` tells `index.py` to render the
    picker table regardless of database count.

    `query_params()` (not a caller-supplied argument) is what lets this be
    computed uniformly from `render()` for every page -- `system.py` and
    `sector.py` otherwise have no way to reach `search.py`/`browse.py`
    without first going back through `browse.py` itself. `index.py` is the
    only script with no `db` in its query string, so it gets just the
    Databases item (or nothing, when no database exists at all).

    Returns:
        str: The `<nav class="sidenav">` element's inner HTML.
    """
    from apiclient import list_databases
    items = []
    try:
        has_databases = bool(list_databases())
    except ApiError:
        # The API itself is unreachable -- the page's own handler is
        # about to hit (or already hit) the same failure and render a
        # clear error page for it; this only decides whether the shared
        # nav chrome around that error page also tries (and fails) to
        # show a "Databases" link, so it fails quiet here instead of
        # taking the whole page shell down with it.
        has_databases = False
    if has_databases:
        items.append(("index.py?all=1", "Databases"))
    db_name = query_params().get("db")
    if db_name:
        db = esc(db_name)
        items.extend([
            (f"search.py?db={db}", "Search"),
            (f"galaxy.py?db={db}", "Galaxy"),
            (f"browse.py?db={db}#sectors", "Sectors"),
            (f"browse.py?db={db}#standalone-systems", "Systems"),
        ])

    # Admin/Login: one extra GET /api/auth/me per page render (same
    # "fails quiet, doesn't take the page down" treatment as the
    # Databases link above) -- an admin session is the exception, not the
    # common case, for this project's read-mostly browsing traffic, so
    # the extra request isn't worth caching/avoiding.
    admin = None
    try:
        admin = auth_me(incoming_cookie_header())
    except ApiError:
        admin = None
    if admin is None:
        items.append(("login.py", "Login"))
    else:
        items.append(("admin.py", "Admin"))
        items.append(("logout.py", "Logout"))

    return "".join(f'<a href="{href}" class="sidenav-item">{label}</a>' for href, label in items)


def render(title, body_html, status="200 OK", set_cookie_headers=None):
    """
    Sends headers and a complete HTML page (shared shell + `body_html`)
    to stdout.

    Args:
        title (str): Page title, also shown as the `<h1>`. Escaped here --
                     callers may build it from database content (a
                     filename, a system name), which is never trusted as
                     pre-sanitized HTML.
        body_html (str): Pre-built HTML for the page body (each caller is
                         responsible for escaping its own interpolated
                         values via `fmt.esc`).
        status (str): CGI status line -- `"200 OK"` unless the caller is
                      rendering an error page.
        set_cookie_headers (list[str], optional): See `send_headers`.
    """
    safe_title = esc(title)
    send_headers(status, set_cookie_headers=set_cookie_headers)
    sys.stdout.write(f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{safe_title} - planetGen</title>
<link rel="stylesheet" href="static/style.css">
</head>
<body>
<nav class="sidenav" aria-label="Main">{_sidenav_html()}</nav>
<header>
  <h1>{safe_title}</h1>
</header>
<main>
{body_html}
</main>
</body>
</html>
""")


def render_error(message, status="404 Not Found", raw=False):
    """
    Renders a minimal error page and exits.

    Args:
        message (str): The error text. HTML-escaped unless `raw=True` --
                       most callers pass user-influenced text (e.g. a
                       `NotFoundError` echoing back an invalid `?db=`
                       value), which must never be interpolated
                       unescaped into the page.
        raw (bool): Set only for the developer-only `PLANETGEN_DEBUG`
                   traceback dump, which is pre-wrapped in `<pre>`.
    """
    body = message if raw else esc(message)
    render("Error", f'<section class="panel"><p class="error">{body}</p></section>', status=status)
    sys.exit(0)


def run(handler):
    """
    Calls `handler()` (which returns `(title, body_html)`) and renders the
    result, turning a `NotFoundError` into a 404 page, an `ApiError`
    (the planetGen API unreachable, or itself erroring) into a 502 page,
    and any other exception into a generic 500 page instead of a raw
    traceback -- this is a public-facing script, so unhandled errors must
    not leak file paths or query text back to the browser.

    Args:
        handler (callable): Zero-argument function returning
                            `(title, body_html)`.
    """
    try:
        title, body_html = handler()
        render(title, body_html)
    except NotFoundError as exc:
        render_error(str(exc), status="404 Not Found")
    except ApiError as exc:
        traceback.print_exc(file=sys.stderr)
        render_error(str(exc), status="502 Bad Gateway")
    except Exception:
        # Always logged to stderr (Apache's error log); only echoed into
        # the page itself when PLANETGEN_DEBUG is set, since a public 500
        # page must not leak file paths or query text by default.
        traceback.print_exc(file=sys.stderr)
        if os.environ.get("PLANETGEN_DEBUG"):
            render_error(f"<pre>{esc(traceback.format_exc())}</pre>", status="500 Internal Server Error", raw=True)
        else:
            render_error("An unexpected error occurred.", status="500 Internal Server Error")
