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
