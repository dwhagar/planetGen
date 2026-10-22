"""
UX.145: the Admin menu (and the Navigate menu beside it) opens at its button and inside the visible window.

The root cause was the action row's size container (`container-type: inline-size`, UX.27), which is also a containing block
for the menus in it: Shoelace placed them against the page and they opened as far from their button as the row is from
the top of the window, off the screen on a short window. The Admin menu is one shared partial (`templates/partials/admin_menu.html`), drawn on the sector, system and phenomenon
pages. Loaded in headless Chromium as the admin, at the window sizes of a laptop, a small laptop, a tablet and a phone
and in both themes, opened on every page that draws it: the open menu (and the submenu of a group) lies inside the
window, and every entry can be reached (the menu scrolls inside itself when it is taller than the space left).
Skipped without Playwright, Chromium or the MySQL test server, like `test_web_a11y.py`.
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

VIEWPORTS = {"laptop": (1280, 720), "small laptop": (1366, 768), "tablet": (768, 1024), "phone": (390, 844),
             "short": (1280, 330)}
"""dict: The window sizes UX.145 names, as (width, height), and a short one where the menu has to open above its button."""

GAP_PX = 16.0
"""float: How far a menu may be from its button (Shoelace's distance is 4)."""

SLACK_PX = 1.0
"""float: How far a box may lie outside the window and still count as inside it (rounding)."""


@pytest.fixture(scope="module")
def menu_pages(site_app, page_targets, admin_token):
    """`endpoint -> path` of every page that draws an Admin menu for the admin."""
    admin = site_app.test_client()
    admin.set_cookie(SESSION_COOKIE_NAME, admin_token)
    pages = {}
    for endpoint, (path, _needs_admin) in sorted(page_targets.items()):
        if path is None:
            continue
        response = admin.get(path)
        if response.status_code == 200 and b'class="admin-menu"' in response.data:
            pages[endpoint] = path
    assert pages, "no page draws an Admin menu: the pages this test is about were not found"
    return pages


def _attached(panel, trigger):
    """Whether `panel` hangs just below `trigger` or stands just above it."""
    below = abs(panel["y"] - (trigger["y"] + trigger["height"])) <= GAP_PX
    above = abs((panel["y"] + panel["height"]) - trigger["y"]) <= GAP_PX
    return below or above


def _outside(box, width, height):
    return [side for side, over in (("left", -box["x"]), ("top", -box["y"]), ("right", box["x"] + box["width"] - width),
                                    ("bottom", box["y"] + box["height"] - height)) if over > SLACK_PX]


@pytest.mark.parametrize("scheme", ["light", "dark"])
@pytest.mark.parametrize("viewport", sorted(VIEWPORTS))
def test_the_admin_menu_opens_inside_the_window_on_every_page_that_has_it(
        browser, base_url, menu_pages, admin_token, viewport, scheme):
    width, height = VIEWPORTS[viewport]
    context = browser.new_context(viewport={"width": width, "height": height}, color_scheme=scheme,
                                  reduced_motion="reduce")
    context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
    problems = []
    try:
        page = context.new_page()
        for endpoint, path in menu_pages.items():
            page.goto(base_url + path, wait_until="load")
            page.locator("sl-dropdown.admin-menu:defined").first.wait_for()
            trigger = page.locator("sl-dropdown.admin-menu >> [slot=trigger]").first
            trigger_box = trigger.bounding_box()
            trigger.click()
            menu = page.locator("sl-dropdown.admin-menu >> sl-menu").first
            menu.wait_for(state="visible")
            popup = page.locator("sl-dropdown.admin-menu >> .dropdown__panel").first
            if not _attached(popup.bounding_box(), trigger_box):
                problems.append(f"{endpoint}: the panel {popup.bounding_box()} is not at its button {trigger_box}")
            for name, box in (("menu", menu.bounding_box()), ("panel", popup.bounding_box())):
                if box is None:
                    problems.append(f"{endpoint}: the {name} has no box")
                    continue
                outside = _outside(box, width, height)
                if outside:
                    problems.append(f"{endpoint}: the {name} {box} lies outside the {width}x{height} window "
                                    f"({', '.join(outside)})")
            # Every entry can be reached: scrolling it into view of its menu puts it inside the window.
            for item in page.locator("sl-dropdown.admin-menu > sl-menu > sl-menu-item, "
                                     "sl-dropdown.admin-menu sl-menu > sl-menu-item").all():
                item.scroll_into_view_if_needed()
                box = item.bounding_box()
                if box is None or _outside(box, width, height):
                    problems.append(f"{endpoint}: the entry {item.inner_text()!r} cannot be reached ({box})")
            page.keyboard.press("Escape")
            navigate = page.locator(".nav-menu >> [slot=trigger]").first
            if navigate.count() and navigate.is_visible():
                trigger_box = navigate.bounding_box()
                navigate.click()
                panel = page.locator(".nav-menu >> .dropdown__panel").first
                panel.wait_for(state="visible")
                box = panel.bounding_box()
                if not _attached(box, trigger_box) or _outside(box, width, height):
                    problems.append(f"{endpoint}: the Navigate panel {box} is not at its button {trigger_box} "
                                    f"inside the {width}x{height} window")
                page.keyboard.press("Escape")
    finally:
        context.close()
    assert not problems, "\n".join(problems)
