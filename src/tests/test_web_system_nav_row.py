"""
UX.27: the system page's buttons (Navigate from here, Navigate to here, Show
on Galaxy Map, Bookmark) sit on one row at every width; where the two
navigate buttons don't fit they fold into one Navigate menu (From here / To
here), chosen by a container query, not a device check.

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
    sample_job_tree,
    sample_params,
    site_app,
    site_db,
)

ROW = ".page-actions.nav-row"
ROW_BUTTONS = ROW + " > a.btn-small:visible, " + ROW + " > sl-dropdown.nav-menu:visible, " + ROW + " > .btn-bookmark:visible"

# (viewport width, buttons on the row): four from 38 rem (608 px) of room,
# three (Navigate menu, Show on Galaxy Map, Bookmark) down to 26 rem
# (416 px), then two (Navigate menu, Bookmark). The page pads 32-48 px.
WIDTHS = [(1280, 4), (820, 4), (600, 3), (390, 2)]


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
        assert page.locator(ROW + " .nav-menu").is_visible() == (count < 4)
        assert page.locator(ROW + " a.nav-wide").first.is_visible() == (count == 4)
        assert page.locator(ROW + " a.nav-galaxy").is_visible() == (count >= 3)
    finally:
        context.close()


@pytest.mark.parametrize("item,link", [("To here", "a.nav-wide >> nth=1"), ("From here", "a.nav-wide >> nth=0"),
                                       ("Show on Galaxy Map", "a.nav-galaxy")])
def test_the_navigate_menu_goes_where_the_wide_buttons_go(browser, base_url, page_targets, admin_token, item, link):
    context, page = _load(browser, base_url, page_targets, admin_token, 390)
    try:
        href = page.locator(ROW + " " + link).evaluate("el => el.getAttribute('href')")
        page.locator(ROW + " .nav-menu > [slot=trigger]").click()
        page.get_by_role("menuitem", name=item).click()
        # The Galaxy Map link redirects (to /galaxy?...), so match its start only.
        page.wait_for_url(lambda url: url.startswith(base_url + (href if item != "Show on Galaxy Map" else "/galaxy")))
    finally:
        context.close()
