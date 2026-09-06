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

from dbutil import NotFoundError


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


def send_headers(status="200 OK"):
    """
    Writes the CGI response status + header block, terminated by the
    required blank line, then flushes -- must be called exactly once,
    before any body output.

    Args:
        status (str): e.g. `"200 OK"`, `"404 Not Found"`.
    """
    # The response declares charset=utf-8 below; stdout's default encoding
    # is platform/locale-dependent (e.g. cp1252 on Windows) and would
    # otherwise raise UnicodeEncodeError on any non-ASCII generated name.
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stdout.write(f"Status: {status}\r\n")
    sys.stdout.write("Content-Type: text/html; charset=utf-8\r\n")
    sys.stdout.write("\r\n")


NAV = '<a href="index.py">Databases</a>'


def render(title, body_html, status="200 OK"):
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
                         values via `dbutil.esc`).
        status (str): CGI status line -- `"200 OK"` unless the caller is
                      rendering an error page.
    """
    from dbutil import esc
    safe_title = esc(title)
    send_headers(status)
    sys.stdout.write(f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{safe_title} - planetGen</title>
<link rel="stylesheet" href="static/style.css">
</head>
<body>
<header>
  <h1>{safe_title}</h1>
  <nav>{NAV}</nav>
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
    from dbutil import esc
    body = message if raw else esc(message)
    render("Error", f'<section class="panel"><p class="error">{body}</p></section>', status=status)
    sys.exit(0)


def run(handler):
    """
    Calls `handler()` (which returns `(title, body_html)`) and renders the
    result, turning a `NotFoundError` into a 404 page and any other
    exception into a generic 500 page instead of a raw traceback -- this
    is a public-facing script, so unhandled errors must not leak file
    paths or query text back to the browser.

    Args:
        handler (callable): Zero-argument function returning
                            `(title, body_html)`.
    """
    try:
        title, body_html = handler()
        render(title, body_html)
    except NotFoundError as exc:
        render_error(str(exc), status="404 Not Found")
    except Exception:
        # Always logged to stderr (Apache's error log); only echoed into
        # the page itself when PLANETGEN_DEBUG is set, since a public 500
        # page must not leak file paths or query text by default.
        traceback.print_exc(file=sys.stderr)
        if os.environ.get("PLANETGEN_DEBUG"):
            from dbutil import esc
            render_error(f"<pre>{esc(traceback.format_exc())}</pre>", status="500 Internal Server Error", raw=True)
        else:
            render_error("An unexpected error occurred.", status="500 Internal Server Error")
