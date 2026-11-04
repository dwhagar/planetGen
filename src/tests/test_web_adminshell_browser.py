# tests/test_web_adminshell_browser.py

"""
The admin shell (ADM.92) in a real browser: the rail, the search box and the density switch pass the same
accessibility and phone-fit checks as every other page, Ctrl+K opens the palette and Enter on a result opens it, the
`g` then `s` keys go to Settings, and the density choice is remembered. Needs Playwright with Chromium and the MySQL
test server, like test_web_a11y.py; skipped without them.

Run with: pytest tests/test_web_adminshell_browser.py
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_a11y import (  # noqa: E402,F401 -- the browser and axe fixtures, and the checks
    SCHEMES, VIEWPORTS, _axe_blocking, _overflow, axe_source, browser,
)
from tests.test_web_settings_browser import _page, config, server  # noqa: E402,F401


@pytest.mark.parametrize("scheme", SCHEMES)
@pytest.mark.parametrize("viewport", sorted(VIEWPORTS))
def test_the_shell_is_accessible_and_fits(browser, axe_source, server, config, viewport, scheme):
    context, page, problems, base_url = _page(browser, server, viewport, scheme)
    try:
        for path in ("/admin", "/admin/queue", "/admin/search?q=cache"):
            assert page.goto(base_url + path, wait_until="load").status == 200
            assert page.locator(".admin-rail a[aria-current='page']").count() == 1
            problems += [f"{path}: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]
    finally:
        context.close()
    assert not problems, f"at {viewport}/{scheme}:\n  " + "\n  ".join(problems)


def test_ctrl_k_opens_the_palette_and_enter_opens_a_result(browser, server, config):
    context, page, problems, base_url = _page(browser, server)
    try:
        page.goto(base_url + "/admin", wait_until="load")
        page.keyboard.press("Control+k")
        page.wait_for_selector("dialog.admin-palette[open] li[role=option]")
        page.keyboard.type("cache size")
        page.wait_for_selector("dialog.admin-palette[open] li[aria-selected=true]")
        page.keyboard.press("Enter")
        page.wait_for_url("**/admin/settings**")
        page.keyboard.press("Control+k")
        page.keyboard.press("Escape")
        assert page.locator("dialog.admin-palette[open]").count() == 0
    finally:
        context.close()
    assert not problems, "\n".join(problems)


def test_g_then_s_goes_to_settings_but_not_while_typing(browser, server, config):
    context, page, problems, base_url = _page(browser, server)
    try:
        page.goto(base_url + "/admin", wait_until="load")
        page.fill("#admin-find-q", "gs")
        assert page.url.endswith("/admin")
        page.locator("body").click()
        page.keyboard.press("g")
        page.keyboard.press("s")
        page.wait_for_url("**/admin/settings")
        page.keyboard.press("?")
        page.wait_for_selector("dialog.admin-keys[open]")
    finally:
        context.close()
    assert not problems, "\n".join(problems)


def test_the_density_choice_is_remembered(browser, server, config):
    context, page, problems, base_url = _page(browser, server)
    try:
        page.goto(base_url + "/admin", wait_until="load")
        assert page.is_checked("[data-admin-density]")
        assert page.get_attribute("body", "data-density") is None
        page.uncheck("[data-admin-density]")
        assert page.get_attribute("body", "data-density") == "roomy"
        page.goto(base_url + "/admin/queue", wait_until="load")
        assert not page.is_checked("[data-admin-density]")
        assert page.get_attribute("body", "data-density") == "roomy"
    finally:
        context.close()
    assert not problems, "\n".join(problems)


@pytest.mark.parametrize("scheme", SCHEMES)
def test_the_shell_stays_in_view_under_the_header_while_a_long_admin_page_scrolls(browser, server, config, scheme):
    """UX.142: the rail, the search and the density switch are not scrolled away; a phone keeps them in the flow."""
    context, page, problems, base_url = _page(browser, server, "desktop", scheme)
    try:
        page.set_viewport_size({"width": 1280, "height": 400 + 200})
        page.goto(base_url + "/admin/settings", wait_until="load")
        page.wait_for_selector(".admin-shell")
        assert page.evaluate("document.documentElement.scrollHeight") > page.evaluate("innerHeight") + 200, "needs a long page"
        page.evaluate("window.scrollTo(0, document.documentElement.scrollHeight)")
        page.wait_for_timeout(200)
        box = page.evaluate("""() => { const s = document.querySelector('.admin-shell').getBoundingClientRect();
          const h = document.querySelector('.site-header').getBoundingClientRect();
          return {shellTop: s.top, shellBottom: s.bottom, headerBottom: h.bottom, win: innerHeight}; }""")
        assert box["shellTop"] == pytest.approx(box["headerBottom"], abs=2), box
        assert 0 < box["shellBottom"] < box["win"], box
        page.set_viewport_size({"width": 390, "height": 844})
        page.wait_for_timeout(200)
        assert page.evaluate("getComputedStyle(document.querySelector('.admin-shell')).position") == "static"
    finally:
        context.close()
    assert not problems, "\n".join(problems)
