# tests/test_logfollow_browser.py

"""
The admin log boxes follow the end of their log (static/logfollow.js): scrolled to the last line on load, kept there as
lines are added, left alone once the reader scrolls up, and following again when they return to the end. Needs
Playwright with Chromium; skipped without it. The script is tested on a bare page, so no site or database is needed.

Run with: pytest tests/test_logfollow_browser.py
"""

from pathlib import Path

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

import os

_LAUNCH = {"executable_path": os.environ["PLAYWRIGHT_CHROMIUM"]} if os.environ.get("PLAYWRIGHT_CHROMIUM") else {}
CSS = Path(__file__).resolve().parent.parent / "html" / "static" / "style.css"
SCRIPT = Path(__file__).resolve().parent.parent / "html" / "static" / "logfollow.js"
LINES = "\n".join(f"line {i}" for i in range(200))
PAGE = f"""<!doctype html><title>t</title>
<pre id="log" data-follow-end style="height:100px;overflow:auto;white-space:pre">{LINES}</pre>"""


@pytest.fixture(scope="module")
def browser():
    with sync_api.sync_playwright() as p:
        try:
            b = p.chromium.launch(**_LAUNCH)
        except Exception as exc:  # no browser installed
            pytest.skip(f"Chromium is not available: {exc}")
        yield b
        b.close()




def _loaded(browser):
    page = browser.new_page()
    page.route("http://logtest.invalid/", lambda route: route.fulfill(content_type="text/html", body=PAGE))
    page.goto("http://logtest.invalid/")
    page.add_script_tag(path=str(SCRIPT))
    return page


def _gap(page):
    return page.evaluate("(() => { const p = document.getElementById('log'); return p.scrollHeight - p.scrollTop - p.clientHeight; })()")


def test_starts_at_the_end_and_follows_new_lines(browser):
    page = _loaded(browser)
    assert page.evaluate("document.getElementById('log').scrollTop") > 0
    assert _gap(page) < 2
    page.evaluate("document.getElementById('log').append('\\nmore\\nmore\\nmore')")
    page.wait_for_function("(() => { const p = document.getElementById('log'); return p.scrollHeight - p.scrollTop - p.clientHeight < 2; })()")
    page.close()


def test_a_reader_who_scrolled_up_is_left_alone_until_they_return(browser):
    page = _loaded(browser)
    page.evaluate("document.getElementById('log').scrollTop = 50")
    page.wait_for_function("document.getElementById('log').dataset.following === 'false'")
    page.evaluate("document.getElementById('log').append('\\nnew line')")
    page.wait_for_timeout(100)
    assert page.evaluate("document.getElementById('log').scrollTop") == 50
    page.evaluate("document.getElementById('log').scrollTop = document.getElementById('log').scrollHeight")
    page.wait_for_function("document.getElementById('log').dataset.following === 'true'")
    page.evaluate("document.getElementById('log').append('\\nlast line')")
    page.wait_for_function("(() => { const p = document.getElementById('log'); return p.scrollHeight - p.scrollTop - p.clientHeight < 2; })()")
    page.close()


def _luminance(rgb):
    def lin(v):
        v /= 255
        return v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4
    r, g, b = (lin(v) for v in rgb)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _numbers(color):
    return [int(float(v)) for v in color[color.index("(") + 1:color.index(")")].replace("/", " ").replace(",", " ").split()[:3]]


@pytest.mark.parametrize("scheme", ["light", "dark"])
def test_each_level_has_its_own_colour_readable_on_the_log_box(browser, scheme):
    """UX.284: error red, warning amber, info blue, debug green, each 4.5:1 or better on the box, light and dark."""
    page = browser.new_page(color_scheme=scheme)
    levels = ("debug", "info", "notice", "warning", "error", "critical")
    spans = "\n".join(f'<span class="log-{lv}" id="{lv}">2026-11-07 10:00:00 EST {lv.upper()} x</span>' for lv in levels)
    page.route("http://logtest.invalid/", lambda route: route.fulfill(content_type="text/html", body=f'<!doctype html><pre class="job-log">{spans}</pre>'))
    page.goto("http://logtest.invalid/")
    page.add_style_tag(path=str(CSS))
    background = _numbers(page.evaluate("getComputedStyle(document.querySelector('.job-log')).backgroundColor"))
    colors = {lv: _numbers(page.evaluate(f"getComputedStyle(document.getElementById('{lv}')).color")) for lv in levels}
    page.close()
    for lv, rgb in colors.items():
        a, b = sorted((_luminance(rgb), _luminance(background)), reverse=True)
        assert (a + 0.05) / (b + 0.05) >= 4.5, f"{lv} in {scheme}: {rgb} on {background}"
    red, amber, blue, green = colors["error"], colors["warning"], colors["info"], colors["debug"]
    assert red[0] > red[1] + 40 and red[0] > red[2] + 40
    assert amber[0] > amber[2] + 60 and amber[1] > amber[2] + 40
    assert blue[2] > blue[0] + 30 and blue[2] > blue[1]
    assert green[1] > green[0] + 30 and green[1] > green[2] + 30
