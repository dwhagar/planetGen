# tests/test_page_shell.py

"""
The shared CGI page shell (`html/lib/page.py`'s `head_html`/`render`/
`SECURITY_HEADERS`) and the versioned static URLs (`html/lib/fmt.py`'s
`static_url`). Same `sys.path` setup as `test_pagination.py`; `auth_me`
is replaced so `render` never needs a running API.
"""

import io
import os
import re
import sys

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_HTML_DIR = os.path.join(_SRC_DIR, "html")
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

import pytest  # noqa: E402

import fmt  # noqa: E402
import page  # noqa: E402
from mdconvert import markdown_to_html  # noqa: E402
from stellarObjects._version import __version__  # noqa: E402


class _Stdout(io.StringIO):
    def reconfigure(self, **kwargs):
        pass


def _capture(call):
    """Runs `call` with stdout swapped for a buffer (inside the test body:
    pytest's own capture resets sys.stdout between fixture setup and the
    test call) and returns what it wrote."""
    out = _Stdout()
    real = sys.stdout
    sys.stdout = out
    try:
        call()
    finally:
        sys.stdout = real
    return out.getvalue()


# --- static_url ---------------------------------------------------------------

def test_static_url_appends_package_version():
    assert fmt.static_url("style.css") == f"static/style.css?v={__version__}"
    assert fmt.STATIC_VERSION == __version__


def test_static_url_keeps_subdirectories():
    assert fmt.static_url("vendor/three.module.min.js") == f"static/vendor/three.module.min.js?v={__version__}"


def test_static_url_quotes_an_odd_version(monkeypatch):
    monkeypatch.setattr(fmt, "STATIC_VERSION", "1.0+local build")
    assert fmt.static_url("a.js") == "static/a.js?v=1.0%2Blocal%20build"


def test_no_page_links_an_unversioned_static_file():
    """Every `<script src>`/`<link href>` into static/ goes through
    static_url; a bare "static/x" URL would be cached for a year by the
    Apache example without ever being refreshed."""
    offenders = []
    for folder in (_HTML_DIR, os.path.join(_HTML_DIR, "lib")):
        for name in os.listdir(folder):
            if name.endswith(".py"):
                text = open(os.path.join(folder, name), encoding="utf-8").read()
                offenders += [f"{name}: {m}" for m in re.findall(r'(?:src|href)="static/[^"]*"', text)]
    assert offenders == []


def test_js_modules_import_siblings_with_their_own_version():
    """A static `import ... from "./x.js"` would drop `?v=`, so the sibling
    would be cached unversioned (and loaded twice if a page also linked
    it by its versioned URL)."""
    static_dir = os.path.join(_HTML_DIR, "static")
    for name in os.listdir(static_dir):
        if name.endswith(".js"):
            text = open(os.path.join(static_dir, name), encoding="utf-8").read()
            assert not re.search(r'^\s*import\s[^(]*from\s', text, re.MULTILINE), name
            for spec in re.findall(r'await import\(`([^`]*)`\)', text):
                assert spec.endswith("${VERSION_QUERY}"), (name, spec)


# --- head and security headers ------------------------------------------------

def test_head_has_viewport_description_theme_colors_and_favicon():
    head = page.head_html("Sector: A & B - planetGen")
    assert '<meta name="viewport" content="width=device-width, initial-scale=1">' in head
    assert "<title>Sector: A &amp; B - planetGen</title>" in head
    assert '<meta name="description" content="' in head
    assert f'content="{page.THEME_COLOR_LIGHT}" media="(prefers-color-scheme: light)"' in head
    assert f'content="{page.THEME_COLOR_DARK}" media="(prefers-color-scheme: dark)"' in head
    assert f'<link rel="icon" type="image/svg+xml" href="static/favicon.svg?v={__version__}">' in head


def test_head_description_is_escaped():
    head = page.head_html("t", description='Say "hi" <b>')
    assert '<meta name="description" content="Say &quot;hi&quot; &lt;b&gt;">' in head


def test_head_loads_theme_script_blocking_before_the_stylesheet():
    head = page.head_html("t")
    theme = f'<script src="static/theme.js?v={__version__}"></script>'
    css = f'<link rel="stylesheet" href="static/style.css?v={__version__}">'
    assert theme in head and css in head
    assert head.index(theme) < head.index(css)
    assert f'<script src="static/navform.js?v={__version__}" defer></script>' in head
    assert "<script>" not in head  # nothing inline: the CSP forbids it


def test_csp_has_the_hardening_directives():
    directives = {d.strip().split(" ", 1)[0]: d.strip().split(" ", 1)[1]
                  for d in page.CONTENT_SECURITY_POLICY.split(";")}
    assert directives == {
        "default-src": "'self'",
        "base-uri": "'self'",
        "form-action": "'self'",
        "frame-ancestors": "'none'",
        "object-src": "'none'",
    }


@pytest.mark.parametrize("send", [
    lambda: page.send_headers(),
    lambda: page.send_json_headers(),
    lambda: page.redirect("login.py"),
])
def test_every_response_kind_carries_the_security_headers(send):
    headers = _capture(send).split("\r\n\r\n", 1)[0]
    for name, value in page.SECURITY_HEADERS:
        assert f"\r\n{name}: {value}" in headers
        assert headers.count(f"{name}:") == 1


def test_render_includes_head_and_hidden_theme_toggle(monkeypatch):
    monkeypatch.setattr(page, "auth_me", lambda cookie: None)
    monkeypatch.setenv("REQUEST_METHOD", "GET")
    monkeypatch.setenv("QUERY_STRING", "db=planetgen")
    body = _capture(lambda: page.render("Hello", "<p>body</p>", description="A test page")).split("\r\n\r\n", 1)[1]
    assert body.startswith("<!doctype html>")
    assert '<meta name="description" content="A test page">' in body
    assert "width=device-width" in body
    nav = body[body.index('<nav class="sidenav"'):body.index("</nav>")]
    assert page.THEME_TOGGLE_HTML in nav
    assert "data-theme-toggle hidden" in page.THEME_TOGGLE_HTML
    assert 'type="button"' in page.THEME_TOGGLE_HTML


# --- stylesheet and static files ----------------------------------------------

def _block(css, opener):
    start = css.index(opener) + len(opener)
    return css[start:css.index("}", start)]


def _tokens(block):
    return dict(re.findall(r"(--[\w-]+):\s*([^;]+);", block))


def test_explicit_dark_theme_matches_the_os_dark_theme():
    css = open(os.path.join(_HTML_DIR, "static", "style.css"), encoding="utf-8").read()
    os_dark = _tokens(_block(css, ':root:not([data-theme="light"]) {'))
    explicit = _tokens(_block(css, ':root[data-theme="dark"] {'))
    assert os_dark and os_dark == explicit
    # Same colours the <meta name="theme-color"> tags use.
    assert os_dark["--bg"] == page.THEME_COLOR_DARK
    assert _tokens(_block(css, ":root {"))["--bg"] == page.THEME_COLOR_LIGHT


def test_favicon_and_theme_script_exist():
    static_dir = os.path.join(_HTML_DIR, "static")
    assert open(os.path.join(static_dir, "favicon.svg"), encoding="utf-8").read().lstrip().startswith("<svg")
    theme = open(os.path.join(static_dir, "theme.js"), encoding="utf-8").read()
    assert '"planetgen-theme"' in theme and "data-theme-toggle" in theme


def test_markdown_tables_scroll_inside_their_own_box():
    html = markdown_to_html("| a | b |\n|---|---|\n| 1 | 2 |\n")
    assert '<div class="table-scroll"><table>' in html
    assert "</table></div>" in html
