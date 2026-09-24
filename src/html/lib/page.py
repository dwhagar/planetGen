# html/lib/page.py

"""
Tiny CGI response/layout helpers shared by every script in `html/`.

No templating engine -- these are plain scripts meant to run with nothing
beyond the standard library, so pages are built as f-strings against one
shared shell (`render`) and a common stylesheet (`static/style.css`,
served directly by Apache, not through CGI).

Every page-to-page link in `html/` posts its parameters as hidden form
fields (`fmt.post_link`) rather than putting them in a plain `<a href=
"page.py?db=...&id=...">`'s query string, so a database name, a sector/
system/phenomenon id, or a search term never ends up in the browser's own
address bar, browser history, or an outgoing `Referer` header. `nav_params`/
`nav_multi_params` below are what a page reads one of these back with.
This is a deliberate trade-off: it makes every page here un-bookmarkable
and un-shareable by URL, and the browser's own back/forward navigation
re-submits the last POST (the usual browser confirmation) rather than
just re-rendering a remembered URL. The one exception is a Galaxy Map/
Sector Map/NAV Map marker plotted inside an `<svg>`, where a `<form>`
can't nest -- see `fmt.data_nav_params` and `static/navform.js`.
"""

import json
import os
import sys
import time
import traceback
from urllib.parse import parse_qs

from apiclient import ApiError, NotFoundError, auth_me
from fmt import esc, post_link, static_url  # noqa: F401 -- post_link re-exported for `from page import post_link` callers

# stellarObjects/ lives at src/stellarObjects/ (src layout); this file is
# at src/html/lib/ -- add src/ to sys.path the same way apiclient.py
# already does, so `stellarObjects.appconfig` is importable here too.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from stellarObjects import log  # noqa: E402
from stellarObjects.appconfig import debug_enabled, load_config  # noqa: E402

_request_started = None
_SECRET_FIELD_WORDS = ("password", "token", "secret", "key", "cookie")


def _redact_fields(fields):
    """Form/query fields for the debug log, with anything credential-like
    (`password`, `new_password`, API key labels' secrets, ...) withheld."""
    return {name: ("<withheld>" if any(word in name.lower() for word in _SECRET_FIELD_WORDS) else values)
            for name, values in fields.items()}


def start_request_log():
    """
    Names this CGI process in the debug log after its script
    (`web/system.py`), keeps log output off stdout (it's the HTTP
    response), and records the incoming request. Called by `run`/
    `run_json`; safe to call more than once.
    """
    global _request_started
    if _request_started is not None:
        return
    _request_started = time.perf_counter()
    script = os.path.basename(os.environ.get("SCRIPT_FILENAME") or sys.argv[0] or "cgi")
    log.set_component(f"web/{script}")
    log.configure(log.NORMAL, console=False)
    if not log.debug_log_active():
        return
    env = os.environ
    log.debug(f"Request: {env.get('REQUEST_METHOD', 'GET')} {env.get('REQUEST_URI') or env.get('SCRIPT_NAME', script)} "
              f"from {env.get('REMOTE_ADDR', '?')} (user agent {env.get('HTTP_USER_AGENT', '?')!r}, referer "
              f"{env.get('HTTP_REFERER', '-')!r}, {'with' if env.get('HTTP_COOKIE') else 'no'} cookie)")
    query = parse_qs(env.get("QUERY_STRING", ""))
    if query:
        log.debug(f"Query parameters: {_redact_fields(query)}")
    if env.get("REQUEST_METHOD", "GET").upper() == "POST":
        log.debug(f"POST fields ({env.get('CONTENT_LENGTH', '0')} bytes): {_redact_fields(form_multi_params())}")


def _log_response(status, kind):
    start_request_log()
    if log.debug_log_active():
        elapsed = ""
        if _request_started is not None:
            elapsed = f" after {(time.perf_counter() - _request_started) * 1000:.1f}ms"
        log.debug(f"Response: {status} ({kind}){elapsed}", stacklevel=3)


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


_raw_body_cache = None


def _raw_post_body():
    """
    Reads and caches stdin's `CONTENT_LENGTH` bytes (per the CGI spec) --
    shared by `form_params`/`form_multi_params` so a page that reads the
    POST body more than once (e.g. `sector.py`'s own module-level wiki-
    upload check and its `handler()`'s admin-action check, both gated on
    "is this a POST") gets the same parsed body back each time instead of
    a second, already-exhausted read silently returning nothing. Safe as a
    process-global cache since a CGI process handles exactly one request
    before exiting.
    """
    global _raw_body_cache
    if _raw_body_cache is None:
        try:
            length = int(os.environ.get("CONTENT_LENGTH", "0") or "0")
        except ValueError:
            length = 0
        _raw_body_cache = sys.stdin.buffer.read(length).decode("utf-8", errors="replace") if length > 0 else ""
    return _raw_body_cache


def form_params():
    """
    Parses an `application/x-www-form-urlencoded` POST body (see
    `_raw_post_body`) -- every admin action form (`login.py`/
    `changecreds.py`/`admin.py`/wiki uploads) and, now, every plain
    navigational link in `html/` too (see `post_link`), since its
    parameters travel as hidden POST fields instead of a query string.

    Returns:
        dict[str, str]: One value per field (the first, if a field
                        repeats) -- same single-valued convention as
                        `query_params`.
    """
    parsed = parse_qs(_raw_post_body())
    return {key: values[0] for key, values in parsed.items()}


def form_multi_params():
    """
    Same as `form_params`, but keeps every value for a repeated field name
    (`search.py`'s tag-facet checkboxes-as-hidden-fields) instead of just
    the first.

    Returns:
        dict[str, list[str]]: Same shape `urllib.parse.parse_qs` returns.
    """
    return parse_qs(_raw_post_body())


def nav_params():
    """
    This request's own navigation parameters (`db`, `id`, and the like) --
    the POST body (`form_params`) for a POST request, which is what every
    navigational link in `html/` now submits (`post_link`) instead of a
    plain `<a href>`, so its parameters never appear in the browser's own
    address bar; falls back to the GET query string (`query_params`) for a
    request with nothing posted at all (e.g. a page rendered directly by
    another page's handler, like `index.py`'s direct call into
    `browse.handler`, with no request of its own to read from).

    Returns:
        dict[str, str]: Same single-valued shape as `query_params`/
                        `form_params`.
    """
    if os.environ.get("REQUEST_METHOD", "GET").upper() == "POST":
        return form_params()
    return query_params()


def nav_multi_params():
    """Multi-valued counterpart to `nav_params`, for a page with repeated
    field names (`search.py`'s tag facets) -- see `form_multi_params`."""
    if os.environ.get("REQUEST_METHOD", "GET").upper() == "POST":
        return form_multi_params()
    return parse_qs(os.environ.get("QUERY_STRING", ""))


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


# The Content-Security-Policy every HTML/JSON response from `html/` carries.
# `default-src 'self'` covers scripts, styles, images, fonts and fetch():
# the shell only loads same-origin `static/` files, the map pages build
# their canvases/WebGL textures in memory (no data:/blob: URLs), and
# nothing inline -- no `<script>` blocks, no `style="..."` attributes (JS
# setting `element.style` is not affected by CSP). The rest are the
# directives `default-src` does not fall back to: no `<base>` hijacking,
# forms only post back to this site, no framing (the modern form of
# X-Frame-Options, kept below for old browsers), and no plugins.
CONTENT_SECURITY_POLICY = ("default-src 'self'; base-uri 'self'; form-action 'self'; "
                           "frame-ancestors 'none'; object-src 'none'")

# The one place the HTML pages' security headers are decided. Apache's
# example vhost used to set the first three on every response too; it no
# longer does for pages (see examples/apache/planetgen.conf.example), so
# these work the same with or without that config. The JSON API sets its
# own stricter set in html/api/app.py.
SECURITY_HEADERS = (
    ("X-Content-Type-Options", "nosniff"),
    ("X-Frame-Options", "DENY"),
    ("Referrer-Policy", "no-referrer"),
    ("Content-Security-Policy", CONTENT_SECURITY_POLICY),
)


def _write_security_headers():
    """Written by `send_headers`, `send_json_headers` and `redirect` --
    see `SECURITY_HEADERS`."""
    for name, value in SECURITY_HEADERS:
        sys.stdout.write(f"{name}: {value}\r\n")


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
    _log_response(status, "HTML page")
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stdout.write(f"Status: {status}\r\n")
    sys.stdout.write("Content-Type: text/html; charset=utf-8\r\n")
    _write_security_headers()
    for cookie in set_cookie_headers or []:
        sys.stdout.write(f"Set-Cookie: {cookie}\r\n")
    sys.stdout.write("\r\n")


def send_json_headers(status="200 OK"):
    """
    Writes the CGI response status + header block for a JSON body
    (`Content-Type: application/json`) instead of `send_headers`'s own
    hardcoded `text/html` -- used only by a script that returns raw JSON
    directly to the browser rather than a rendered page (today, just
    `html/galaxy_tiles.py`, the interactive 3D Galaxy Map's own tile
    proxy: its client-side JS calls it directly via `fetch()`,
    unlike every other page here, which is rendered server-side and never
    fetched by the browser's own script).

    Args:
        status (str): CGI status line -- `"200 OK"` unless the caller is
                      reporting an error.
    """
    _log_response(status, "JSON")
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stdout.write(f"Status: {status}\r\n")
    sys.stdout.write("Content-Type: application/json; charset=utf-8\r\n")
    _write_security_headers()
    sys.stdout.write("\r\n")


def redirect(url, set_cookie_headers=None):
    """
    Writes a CGI 302 redirect response to stdout and nothing else -- must
    be called instead of (never alongside) `send_headers`/`render`, and
    before any other output, same as `send_headers`'s own requirement.

    Used by `login.py`/`changecreds.py` to redirect after a successful
    POST, carrying the API's `Set-Cookie` response along (see
    `send_headers`'s own parameter of the same name), and by `admin.py`
    to bounce an unauthenticated/stale-credentials visitor to `login.py`/
    `changecreds.py`. Every one of these targets is a bare page name with
    no params of its own to carry -- a redirect's `Location` URL is
    necessarily visible in the browser's own address bar, same as any
    other URL, which is exactly why `index.py` calls `browse.handler`
    in-process instead of redirecting to `browse.py?db=...` (see that
    page's own module docstring) and why every other page-to-page link in
    `html/` posts its params instead of putting them in a URL at all (see
    `post_link`).

    Args:
        url (str): The target URL, e.g. `"login.py"`. Callers are
                   responsible for URL-encoding any dynamic piece of it
                   themselves (see `urllib.parse.quote`) -- this just
                   writes the header block, the same division of
                   responsibility `send_headers` already has for the
                   `Content-Type` header.
        set_cookie_headers (list[str], optional): See `send_headers`.
    """
    _log_response("302 Found", f"redirect to {url}")
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
    a fixed set of "functions" rather than a breadcrumb -- Search, Galaxy,
    Sectors, Systems, Nav, and Phenomena for the current `db` when one is
    present in the request (`nav_params`, a hidden POST field for most
    requests now -- see that function). Sectors/Systems jump
    straight to the matching anchor on `browse.py` (`id="sectors"`/
    `id="standalone-systems"`) rather than duplicating that page's own
    listing here; Galaxy goes to `galaxy.py`, the galaxy-scale map of every
    sector that actually has a galaxy position (see that page's own
    docstring) -- a different view than `browse.py`'s flat sector list.
    Nav goes to `nav.py` with no `from=` -- that page's own sector-then-
    system picker is what lets NAV be reached from here rather than only
    ever via a specific system's own "Navigate from here" link. Phenomena
    goes to `phenomena.py`, the flat list of every exotic phenomenon (see
    that page's own docstring).

    No "Databases" item -- `index.py` no longer renders a picker table to
    link to at all (it renders `browse.py`'s own content directly, in
    place, for this deployment's one database -- see its own module
    docstring), so there was no longer a destination for one.

    `nav_params()` (not a caller-supplied argument) is what lets this be
    computed uniformly from `render()` for every page -- `system.py` and
    `sector.py` otherwise have no way to reach `search.py`/`browse.py`
    without first going back through `browse.py` itself. `index.py` is the
    only script with no `db` of its own, so it gets none of these (just
    Login/Admin below).

    Each item posts to its target instead of linking to it (see
    `post_link`), so `db` never shows up in the address bar just from
    using the sidenav. A URL fragment (`#sectors`/`#standalone-systems`)
    is never sent to the server either way, so `Sectors`/`Systems` can
    still target one directly on `post_link`'s own `action` -- the browser
    scrolls to it after the POST navigates there, same as a plain link.

    Returns:
        str: The `<nav class="sidenav">` element's inner HTML.
    """
    items = []
    db_name = nav_params().get("db")
    if db_name:
        items.extend([
            ("search.py", {"db": db_name}, "Search"),
            ("galaxy.py", {"db": db_name}, "Galaxy"),
            ("browse.py#sectors", {"db": db_name}, "Sectors"),
            ("browse.py#standalone-systems", {"db": db_name}, "Systems"),
            ("nav.py", {"db": db_name}, "Nav"),
            ("phenomena.py", {"db": db_name}, "Phenomena"),
        ])

    # Admin/Login: one extra GET /api/auth/me per page render (same
    # "fails quiet, doesn't take the page down" treatment as the
    # Databases link above) -- an admin session is the exception, not the
    # common case, for this project's read-mostly browsing traffic, so
    # the extra request isn't worth caching/avoiding.
    admin = None
    try:
        admin = auth_me(incoming_cookie_header())
    except (ApiError, NotFoundError):
        admin = None
    if admin is None:
        items.append(("login.py", {}, "Login"))
    else:
        items.append(("admin.py", {}, "Admin"))
        items.append(("adminstats.py", {"db": db_name} if db_name else {}, "Stats"))
        items.append(("logout.py", {}, "Logout"))

    return "".join(
        post_link(action, params, label, css_class="sidenav-item")
        for action, params, label in items
    ) + THEME_TOGGLE_HTML


# Light/dark/system switch at the bottom of the side nav. Starts `hidden`:
# `static/theme.js` reveals and wires it, so a browser without JavaScript
# never shows a button that does nothing (the page just follows the OS
# theme). The label names the current setting; `aria-pressed` isn't used
# because this is a three-way cycle, not an on/off switch.
THEME_TOGGLE_HTML = ('<button type="button" class="link-btn sidenav-item theme-toggle" '
                     'data-theme-toggle hidden>Theme: System</button>')

# Browser UI colour (address bar etc.) per OS theme -- `--bg` from
# static/style.css's light and dark token blocks. theme.js repoints both
# at one colour when the visitor picks a theme explicitly.
THEME_COLOR_LIGHT = "#f6f7fb"
THEME_COLOR_DARK = "#14151e"


def head_html(title, description=None):
    """
    The shared `<head>` contents (without the `<head>` tags themselves):
    charset, viewport, title, description, theme colours, favicon, the
    stylesheet, and the site-wide scripts. Split out of `render` so it
    can be tested on its own.

    `static/theme.js` is a small blocking script placed before the
    stylesheet on purpose: it only reads localStorage and sets
    `data-theme` on `<html>`, so the first paint already uses the chosen
    theme instead of flashing the OS one first. Everything else is
    `defer`/module.

    Args:
        title (str): The full, unescaped `<title>` text.
        description (str, optional): Page-specific `<meta name=
            "description">`; defaults to a site-wide one.

    Returns:
        str: HTML, every interpolated value escaped.
    """
    site_name = load_config()["site_name"]
    if not description:
        description = (f"{site_name}: a browsable, procedurally generated galaxy of sectors, "
                       f"star systems, planets and phenomena.")
    return f"""<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)}</title>
<meta name="description" content="{esc(description)}">
<meta name="color-scheme" content="light dark">
<meta name="theme-color" content="{THEME_COLOR_LIGHT}" media="(prefers-color-scheme: light)">
<meta name="theme-color" content="{THEME_COLOR_DARK}" media="(prefers-color-scheme: dark)">
<link rel="icon" type="image/svg+xml" href="{static_url('favicon.svg')}">
<script src="{static_url('theme.js')}"></script>
<link rel="stylesheet" href="{static_url('style.css')}">
<script src="{static_url('navform.js')}" defer></script>"""


def render(title, body_html, status="200 OK", set_cookie_headers=None, description=None):
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
        description (str, optional): `<meta name="description">` text for
            this page (see `head_html`); a site-wide default otherwise.
    """
    safe_title = esc(title)
    site_name = load_config()["site_name"]
    log.debug(f"Rendering page {title!r} ({len(body_html)} bytes of body HTML)")
    send_headers(status, set_cookie_headers=set_cookie_headers)
    sys.stdout.write(f"""<!doctype html>
<html lang="en">
<head>
{head_html(f"{title} - {site_name}", description)}
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
        raw (bool): Treat `message` as already-safe HTML.
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
    start_request_log()
    try:
        title, body_html = handler()
        render(title, body_html)
    except NotFoundError as exc:
        log.debug(f"Not found: {exc}")
        render_error(str(exc), status="404 Not Found")
    except ApiError as exc:
        traceback.print_exc(file=sys.stderr)
        log.exception(f"API error while building the page: {exc}")
        render_error(str(exc), status="502 Bad Gateway")
    except Exception:
        # Always logged to stderr (Apache's error log), and to the debug
        # log when debug is on -- never into the page itself, which is
        # public and must not leak file paths or query text, even when
        # someone has left debug on.
        traceback.print_exc(file=sys.stderr)
        log.exception("Unhandled exception while building the page")
        render_error(_unexpected_error_message(), status="500 Internal Server Error")


def _unexpected_error_message():
    if debug_enabled():
        return (f"An unexpected error occurred. The traceback is in the debug log "
                f"(process {os.getpid()}, {time.strftime('%Y-%m-%d %H:%M:%S')}).")
    return "An unexpected error occurred."


def render_json_error(message, status="400 Bad Request"):
    """`run_json`'s JSON counterpart to `render_error` -- writes `{"error":
    message}` instead of an HTML error page, and exits, same as that
    function."""
    send_json_headers(status)
    sys.stdout.write(json.dumps({"error": message}))
    sys.exit(0)


def run_json(handler):
    """
    `run`'s JSON counterpart: calls `handler()` (which returns a
    JSON-serializable value, usually a `dict`) and writes it as the whole
    response body, with the same exception handling `run` gives an HTML
    page (`NotFoundError` -> 404, `ApiError` -> 502, anything else -> 500)
    but via `render_json_error` instead of a rendered error page -- used
    by a script whose only job is a browser `fetch()` target (today, just
    `html/galaxy_tiles.py`), never a page a person navigates to directly.

    Args:
        handler (callable): Zero-argument function returning a
                            JSON-serializable value.
    """
    start_request_log()
    try:
        payload = handler()
        send_json_headers()
        body = json.dumps(payload)
        log.debug(f"JSON response body: {len(body)} bytes")
        sys.stdout.write(body)
    except NotFoundError as exc:
        log.debug(f"Not found: {exc}")
        render_json_error(str(exc), status="404 Not Found")
    except ApiError as exc:
        traceback.print_exc(file=sys.stderr)
        log.exception(f"API error while building the response: {exc}")
        render_json_error(str(exc), status="502 Bad Gateway")
    except Exception:
        traceback.print_exc(file=sys.stderr)
        log.exception("Unhandled exception while building the response")
        render_json_error(_unexpected_error_message(), status="500 Internal Server Error")


# Under a real CGI server, start logging before the page's own module-level
# code runs (some pages make API calls or redirect before ever reaching
# `run`).
if os.environ.get("GATEWAY_INTERFACE"):
    start_request_log()
