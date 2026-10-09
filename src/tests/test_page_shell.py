# tests/test_page_shell.py

"""
Site-wide page pieces: the versioned static URLs (`planetgen/web/lib/fmt.py`'s
`static_url`), the pages' security headers (`web.SECURITY_HEADERS`), and
the stylesheet's theme tokens. The page shell itself (`base.html`) is
covered by `test_web_pages.py`. Same `sys.path` setup as
`test_pagination.py`.
"""

import os
import re

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_HTML_DIR = os.path.join(_SRC_DIR, "html")

from planetgen.web.lib import fmt  # noqa: E402
from planetgen import web  # noqa: E402
from planetgen.web.lib.mdconvert import markdown_to_html  # noqa: E402
from planetgen._version import __version__  # noqa: E402


# --- static_url ---------------------------------------------------------------

def test_static_url_appends_package_version_and_the_static_fingerprint():
    fingerprint = fmt._static_fingerprint()
    assert re.fullmatch(r"[0-9a-f]{8}", fingerprint)
    assert fmt.static_url("style.css") == f"static/style.css?v={__version__}-{fingerprint}"
    assert fmt.STATIC_VERSION == __version__


def test_static_url_keeps_subdirectories():
    url = fmt.static_url("vendor/three.module.min.js")
    assert url.startswith("static/vendor/three.module.min.js?v=" + __version__ + "-")


def test_static_url_quotes_an_odd_version(monkeypatch):
    monkeypatch.setattr(fmt, "STATIC_VERSION", "1.0+local build")
    monkeypatch.setattr(fmt, "_static_fingerprint", lambda: "abcd1234")
    assert fmt.static_url("a.js") == "static/a.js?v=1.0%2Blocal%20build-abcd1234"


def test_a_changed_static_file_gets_a_new_url_even_within_one_release(tmp_path, monkeypatch):
    """Scripts deployed before the release stamp lands (or an Apache that
    wasn't restarted) must not be served under the URL their previous
    content was cached at for a year: two modules from two versions meet,
    an import fails, and the maps stay blank."""
    static = tmp_path / "static"
    static.mkdir()
    (static / "a.js").write_text("export const a = 1;\n")
    monkeypatch.setattr(fmt, "_STATIC_DIR", str(static))
    monkeypatch.setattr(fmt, "STATIC_FINGERPRINT_TTL_SECONDS", 0.0)
    before = fmt.static_url("a.js")
    assert fmt.static_url("a.js") == before
    (static / "a.js").write_text("export const a = 2; // changed\n")
    assert fmt.static_url("a.js") != before
    (static / "new.js").write_text("export {};\n")
    assert fmt.static_url("a.js") != before


def test_no_page_links_an_unversioned_static_file():
    """Every `<script src>`/`<link href>` into static/ goes through
    static_url; a bare "static/x" URL would be cached for a year by the
    Apache example without ever being refreshed."""
    offenders = []
    web_dir = os.path.join(_SRC_DIR, "planetgen", "web")
    for folder in (os.path.join(web_dir, "lib"), os.path.join(web_dir, "maps"), web_dir,
                   os.path.join(web_dir, "templates"),
                   os.path.join(web_dir, "templates", "partials")):
        for name in os.listdir(folder):
            if name.endswith((".py", ".html")):
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


def test_vendored_shoelace_files_import_only_files_that_exist():
    """UX.40: scripts/vendor_shoelace.py copies a component with the chunks it
    imports; a missing one would break every page's menus at load time."""
    root = os.path.join(_HTML_DIR, "static", "vendor", "shoelace")
    imports = re.compile(r'(?:from|import)\s*"(\.[^"]+)"|import\("(\.[^"]+)"\)')
    checked = 0
    for folder, _dirs, names in os.walk(root):
        for name in names:
            if name.endswith(".js"):
                text = open(os.path.join(folder, name), encoding="utf-8").read()
                for match in imports.finditer(text):
                    target = os.path.normpath(os.path.join(folder, match.group(1) or match.group(2)))
                    assert os.path.exists(target), (name, target)
                    checked += 1
    assert checked > 0
    components = open(os.path.join(_HTML_DIR, "static", "components.js"), encoding="utf-8").read()
    for name in re.findall(r'"([a-z-]+/[a-z-]+)",', components):
        assert os.path.exists(os.path.join(root, "components", f"{name}.js")), name


# --- security headers -----------------------------------------------------------

def test_csp_has_the_hardening_directives():
    directives = {d.strip().split(" ", 1)[0]: d.strip().split(" ", 1)[1]
                  for d in web.CONTENT_SECURITY_POLICY.split(";")}
    assert directives == {
        "default-src": "'self'",
        "base-uri": "'self'",
        "form-action": "'self'",
        "frame-ancestors": "'none'",
        "object-src": "'none'",
    }
    assert ("Content-Security-Policy", web.CONTENT_SECURITY_POLICY) in web.SECURITY_HEADERS


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
    # Same colours base.html's <meta name="theme-color"> tags use.
    base = open(os.path.join(_SRC_DIR, "planetgen", "web", "templates", "base.html"), encoding="utf-8").read()
    assert f'content="{os_dark["--bg"]}" media="(prefers-color-scheme: dark)"' in base
    light = _tokens(_block(css, ":root {"))["--bg"]
    assert f'content="{light}" media="(prefers-color-scheme: light)"' in base


def test_favicon_and_theme_script_exist():
    static_dir = os.path.join(_HTML_DIR, "static")
    assert open(os.path.join(static_dir, "favicon.svg"), encoding="utf-8").read().lstrip().startswith("<svg")
    theme = open(os.path.join(static_dir, "theme.js"), encoding="utf-8").read()
    assert '"planetgen-theme"' in theme and "data-theme-toggle" in theme


def test_markdown_tables_scroll_inside_their_own_box():
    html = markdown_to_html("| a | b |\n|---|---|\n| 1 | 2 |\n")
    assert '<div class="table-scroll" tabindex="0"><table>' in html
    assert "</table></div>" in html


def _stylesheet():
    return open(os.path.join(_HTML_DIR, "static", "style.css"), encoding="utf-8").read()


def test_one_shared_rule_spaces_every_button_group():
    """UX.16: one gap token, larger on touch screens, used by a rule that
    matches any box holding two or more buttons."""
    css = _stylesheet()
    assert re.search(r":root \{[^}]*--btn-gap: 0\.5rem;", css)
    assert re.search(r"@media \(pointer: coarse\) \{\s*:root \{\s*--btn-gap: 0\.75rem;", css)
    rule = re.search(r":is\(div, p, span, form[^{]*\{([^}]*)\}", css)
    assert rule, "the shared button-gap rule is missing"
    assert "~ :is(button, .btn, .starmap-btn)" in rule.group(0)
    assert "gap: var(--btn-gap);" in rule.group(1) and "flex-wrap: wrap;" in rule.group(1)


def test_data_sits_beside_the_render_when_there_is_room():
    """UX.15: container queries (not viewport ones) put the System Map's
    info panel and a phenomenon's data beside the render."""
    css = _stylesheet()
    assert "container: sysmap-layout / inline-size;" in css
    assert "@container sysmap-layout (min-width: 46rem)" in css
    assert "container: object-view / inline-size;" in css
    assert "@container object-view (min-width: 50rem)" in css
