# tests/test_web_icons_under_csp_browser.py

"""
UX.166: the site's icons and habitability symbols in a real browser under the site's own Content-Security-Policy, and a
second policy header on top of it (the way an Apache in front of the site adds one).

The sprite (`static/icons.svg`) is a separate document, and Firefox refuses an inline style inside it under
`default-src 'self'` (Chrome did not), which left every icon a solid black shape on Boss's Firefox. The symbols' fill and
outline are in `static/style.css` now. Every symbol and life icon is drawn here in each browser that is installed, in
both themes, from a plain static server with both policy headers, and each must show its own colour: pixels of the
tier colour (shaded by the symbol's opacity over the page) in at least `MIN_SHARE` of its box. A black shape on the dark
page, the old fault, shows none.

Skipped without Playwright, Pillow or the browser: Chromium is the one `python -m playwright install chromium` fetches,
and Firefox is skipped until `python -m playwright install firefox` has been run.
"""

import functools
import http.server
import io
import os
import shutil
import threading

import pytest

pytest.importorskip("playwright.sync_api")
Image = pytest.importorskip("PIL.Image")

from playwright import sync_api  # noqa: E402

from planetgen.physics.habitability_world import TIER_NAMES  # noqa: E402
from planetgen.web import CONTENT_SECURITY_POLICY, STATIC_DIR  # noqa: E402
from planetgen.web.lib import habitability_symbols as symbols  # noqa: E402

BROWSERS = ("chromium", "firefox")
"""tuple: The browsers checked, each skipped when it is not installed."""

THEMES = ("light", "dark")
"""tuple: The `data-theme` values checked."""

SECOND_POLICY = "default-src 'self'; style-src 'self'; img-src 'self' data:"
"""str: The extra policy header an Apache in front of the site adds; the browser applies both and the stricter wins."""

MIN_SHARE = 0.12
"""float: The least share of a symbol's box that must be its own colour (a trefoil or a ring has about a third)."""

COLOUR_TOLERANCE = 14
"""int: How far each channel of a pixel may be from the expected colour and count as it."""

SCALE = 3
"""int: The device scale factor of the screenshots, for solid interiors."""

PAGE = """<!doctype html>
<html lang="en" data-theme="{theme}"><head><meta charset="utf-8"><title>symbols</title>
<link rel="stylesheet" href="style.css"></head>
<body><p class="hab-symbols" id="row">{row}</p><p id="life">{life}</p>
<p><svg class="icon" id="pencil" role="img" aria-label="edit"><use href="icons.svg#edit"></use></svg></p></body></html>
"""


class _Handler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Content-Security-Policy", CONTENT_SECURITY_POLICY)
        self.send_header("Content-Security-Policy", SECOND_POLICY)
        super().end_headers()

    def log_message(self, *args):  # noqa: D401 -- keep the test output quiet
        pass


@pytest.fixture(scope="module")
def site(tmp_path_factory):
    """The base URL of a static server of `style.css`, `icons.svg` and one page per theme."""
    root = tmp_path_factory.mktemp("icons-under-csp")
    for name in ("style.css", "icons.svg"):
        shutil.copy(os.path.join(STATIC_DIR, name), root / name)
    row = "".join(symbols.symbol_html(name, tier, "icons.svg")
                  for name in symbols.SYMBOL_NAMES for tier in range(len(TIER_NAMES)))
    life = "".join(symbols.life_html(name, 5, "icons.svg") for name in symbols.LIFE_NAMES)
    for theme in THEMES:
        (root / f"{theme}.html").write_text(PAGE.format(theme=theme, row=row, life=life), encoding="utf-8")
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(_Handler, directory=str(root)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_address[1]}"
    finally:
        server.shutdown()
        thread.join()


@pytest.fixture(scope="module", params=BROWSERS)
def browser(request):
    try:
        playwright = sync_api.sync_playwright().start()
    except Exception as exc:  # noqa: BLE001 -- any failure here means "no browser"
        pytest.skip(f"Playwright could not start ({exc})")
    try:
        launched = getattr(playwright, request.param).launch()
    except Exception as exc:  # noqa: BLE001
        playwright.stop()
        pytest.skip(f"No {request.param} for Playwright ({str(exc).splitlines()[0]}); run "
                    f"`python -m playwright install {request.param}`.")
    try:
        yield launched
    finally:
        launched.close()
        playwright.stop()


def _parse_rgb(css_colour):
    """`(r, g, b)` of a computed `rgb(...)` or `rgba(...)` colour."""
    parts = css_colour[css_colour.index("(") + 1:css_colour.index(")")].replace("/", ",").replace(" ", ",").split(",")
    return tuple(float(part) for part in parts if part)[:3]


def _blend(colour, background, opacity):
    return tuple(opacity * c + (1.0 - opacity) * b for c, b in zip(colour, background))


def _own_colour_share(png, expected):
    """The share of a screenshot's pixels within `COLOUR_TOLERANCE` of the colour `expected`."""
    image = Image.open(io.BytesIO(png)).convert("RGB")
    pixels = list(image.getdata())
    near = sum(all(abs(c - e) <= COLOUR_TOLERANCE for c, e in zip(pixel, expected)) for pixel in pixels)
    return near / len(pixels)


@pytest.mark.parametrize("theme", THEMES)
def test_every_symbol_and_icon_shows_its_own_colour_under_both_policy_headers(browser, site, theme):
    page = browser.new_page(device_scale_factor=SCALE, viewport={"width": 900, "height": 500})
    try:
        page.goto(f"{site}/{theme}.html", wait_until="load")
        page.wait_for_timeout(300)
        background = _parse_rgb(page.evaluate("getComputedStyle(document.body).backgroundColor"))
        targets = list(page.locator("#row svg.hab-symbol").all()) + list(page.locator("#life svg.hab-life-shape").all())
        targets.append(page.locator("#pencil"))
        assert len(targets) == len(symbols.SYMBOL_NAMES) * len(TIER_NAMES) + len(symbols.LIFE_NAMES) + 1
        for target in targets:
            label = target.get_attribute("aria-label") or target.get_attribute("class")
            colour = _parse_rgb(target.evaluate("el => getComputedStyle(el).color"))
            opacity = float(target.evaluate("el => getComputedStyle(el).opacity"))
            pencil = target.get_attribute("id") == "pencil"
            share = _own_colour_share(target.screenshot(), _blend(colour, background, opacity))
            floor = 0.03 if pencil else MIN_SHARE  # the pencil is an outline: a thin stroke
            assert share >= floor, f"{label} in {theme}: {share:.3f} of its box is its own colour (want {floor})"
    finally:
        page.close()
