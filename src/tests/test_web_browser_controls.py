"""
UX.21: every control does something, and the maps' open menus hold no overlap.

Two checks over every GET page of the `web` blueprint (the pages and
fixtures of `test_web_a11y.py`), in headless Chromium:

- every plain button on a page (not the maps' own buttons, which
  `test_web_browser_maps.py` covers one by one, and not a submit button)
  changes something when it is clicked: the DOM, the address or a
  request. A button that does nothing is a dead control.
- with one of the maps' popover menus (Menu, Bookmarks, the history
  panel) open at each width, no two controls inside that menu overlap.
  The menu dropping over the page's own rows is by design (TEST.106), so
  only pairs inside the menu count.

Skipped without Playwright, Chromium or the MySQL test server, like
`test_web_browser_layout.py`.
"""

import pytest

pytest.importorskip("playwright.sync_api")

from planetgen.api.authz import SESSION_COOKIE_NAME  # noqa: E402

from tests.test_web_a11y import (  # noqa: E402,F401 -- fixtures
    PAGE_ENDPOINTS,
    admin_token,
    base_url,
    browser,
    page_targets,
    sample_job,
    sample_job_tree,
    sample_params,
    site_app,
    site_db,
)
from tests.test_web_browser_layout import CONTROLS, POPOVER_MENUS, WIDTHS  # noqa: E402

PLAIN_BUTTONS = ("css:light=button:not([type=submit]):not([disabled]):not(.starmap-btn), "
                 "sl-button:not([disabled]):not([type=submit])")

WATCH = """() => {
  window.__changes = 0;
  new MutationObserver((records) => { window.__changes += records.length; })
    .observe(document, { subtree: true, attributes: true, childList: true, characterData: true });
}"""

# The controls inside an open popover menu that overlap one another.
MENU_OVERLAPS = """([menu, selector]) => {
  const root = document.querySelectorAll(menu)[0];
  const boxes = [...root.querySelectorAll(selector)].filter((el) => el.checkVisibility()).map((el) => ({ el, r: el.getBoundingClientRect() }))
    .filter(({ r }) => r.width > 1 && r.height > 1);
  const found = [];
  for (let i = 0; i < boxes.length; i++) {
    for (let j = i + 1; j < boxes.length; j++) {
      const a = boxes[i], b = boxes[j];
      if (a.el.contains(b.el) || b.el.contains(a.el)) continue;
      const w = Math.min(a.r.right, b.r.right) - Math.max(a.r.left, b.r.left);
      const h = Math.min(a.r.bottom, b.r.bottom) - Math.max(a.r.top, b.r.top);
      if (w > 0.5 && h > 0.5) found.push(`${a.el.textContent.trim().slice(0, 24)} / ${b.el.textContent.trim().slice(0, 24)} ${Math.round(w)}x${Math.round(h)}`);
    }
  }
  return found;
}"""


def _context(browser, base_url, admin_token, needs_admin, width):
    context = browser.new_context(viewport={"width": width, "height": 900}, device_scale_factor=1,
                                  reduced_motion="reduce")
    if needs_admin:
        context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
    return context


@pytest.mark.parametrize("endpoint", PAGE_ENDPOINTS)
def test_every_plain_button_does_something(browser, base_url, page_targets, admin_token, endpoint):
    path, needs_admin = page_targets[endpoint]
    if path is None:
        pytest.skip(needs_admin)
    context = _context(browser, base_url, admin_token, needs_admin, 1280)
    page = context.new_page()
    dead = []
    try:
        page.goto(base_url + path, wait_until="load")
        page.wait_for_timeout(500)
        count = page.locator(PLAIN_BUTTONS).count()
        for index in range(count):
            page.goto(base_url + path, wait_until="load")
            page.wait_for_timeout(300)
            button = page.locator(PLAIN_BUTTONS).nth(index)
            if not button.is_visible():
                continue
            name = button.evaluate("""(e) => (e.getAttribute("aria-label") || e.textContent || e.className || e.tagName).trim().slice(0, 40)""")
            page.evaluate(WATCH)
            requests = []
            record = lambda request: requests.append(request)
            page.on("request", record)
            before = page.url
            button.click()
            page.wait_for_timeout(300)
            page.remove_listener("request", record)
            if not page.evaluate("() => window.__changes") and page.url == before and not requests:
                dead.append(name)
    finally:
        context.close()
    assert not dead, f"{path}: buttons that do nothing: {dead}"


@pytest.mark.parametrize("width", WIDTHS)
@pytest.mark.parametrize("endpoint", PAGE_ENDPOINTS)
def test_open_map_menus_hold_no_overlap(browser, base_url, page_targets, admin_token, endpoint, width):
    path, needs_admin = page_targets[endpoint]
    if path is None:
        pytest.skip(needs_admin)
    context = _context(browser, base_url, admin_token, needs_admin, width)
    page = context.new_page()
    problems = []
    try:
        page.goto(base_url + path, wait_until="load")
        page.wait_for_timeout(500)
        for index in range(page.locator(POPOVER_MENUS).count()):
            menu = page.locator(POPOVER_MENUS).nth(index)
            if not menu.is_visible():
                continue
            menu.locator("summary").first.click()
            page.wait_for_timeout(150)
            problems += page.evaluate(
                """([popovers, index, selector]) => {
                  const menu = document.querySelectorAll(popovers)[index];
                  menu.setAttribute("data-audit", "1");
                  return (""" + MENU_OVERLAPS + """)(['[data-audit]', selector]);
                }""", [POPOVER_MENUS, index, CONTROLS])
            # Closed in the page, not by a click: under load the open panel can still be
            # settling over the button, and the overlap was already measured above.
            menu.evaluate("m => { m.open = false; }")
    finally:
        context.close()
    assert not problems, f"{path} at {width} px: {sorted(set(problems))}"
