"""
UX.27, UX.201: the system page's buttons (the Nav menu, Show on Galaxy Map, Bookmark) sit on one row at every width. The
Nav menu holds Navigate from here, Navigate to here and Nearby locations.

Skipped without Playwright, Chromium or the MySQL test server, like
`test_web_a11y.py`.
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
    sample_run,
    sample_params,
    site_app,
    site_db,
)

ROW = ".page-actions.nav-row"
ROW_BUTTONS = ROW + " > a.btn-small:visible, " + ROW + " > sl-dropdown.nav-menu:visible, " + ROW + " > .btn-bookmark:visible"

# (viewport width, buttons on the row): the Nav menu, Show on Galaxy Map and Bookmark at every width.
WIDTHS = [(1280, 3), (820, 3), (600, 3), (390, 3)]


def _load(browser, base_url, page_targets, admin_token, width):
    path, _needs_admin = page_targets["web.system"]
    context = browser.new_context(viewport={"width": width, "height": 900}, reduced_motion="reduce")
    context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
    page = context.new_page()
    page.goto(base_url + path, wait_until="load")
    page.locator("sl-dropdown.nav-menu:defined").wait_for(state="attached")
    return context, page


@pytest.mark.parametrize("width,count", WIDTHS)
def test_system_buttons_share_one_row(browser, base_url, page_targets, admin_token, width, count):
    context, page = _load(browser, base_url, page_targets, admin_token, width)
    try:
        # The bookmark button is shown by a script a moment after load.
        page.locator(ROW + " > .btn-bookmark").wait_for(state="visible")
        boxes = [el.bounding_box() for el in page.locator(ROW_BUTTONS).all()]
        assert len(boxes) == count, boxes
        centers = [box["y"] + box["height"] / 2 for box in boxes]
        assert max(centers) - min(centers) < 2, f"buttons on more than one row: {boxes}"
        for first, second in zip(boxes, boxes[1:]):
            assert first["x"] + first["width"] <= second["x"] + 0.5, f"overlap: {first} {second}"
        assert boxes[-1]["x"] + boxes[-1]["width"] <= width, boxes
        assert page.locator(ROW + " .nav-menu").is_visible()
    finally:
        context.close()


def test_the_nav_menu_lists_navigate_from_to_and_nearby(browser, base_url, page_targets, admin_token):
    context, page = _load(browser, base_url, page_targets, admin_token, 390)
    try:
        page.locator(ROW + " .nav-menu > [slot=trigger]").click()
        names = [name.strip() for name in page.locator(ROW + " .nav-menu sl-menu-item").all_inner_texts()]
        assert names[:2] == ["Navigate from here", "Navigate to here"], names
        assert set(names[2:]) <= {"Nearby locations"}, names
        page.get_by_role("menuitem", name="Navigate to here").click()
        page.wait_for_url(lambda url: "/nav" in url and "to=" in url)
    finally:
        context.close()
