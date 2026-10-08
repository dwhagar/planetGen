"""
UX.26 / UX.31: the admin edit actions are one Edit button that opens a menu,
and each pick opens its confirm step in a dialog on top of the page.

Loaded in headless Chromium as the admin on a system page (a row of Edit
menus), with the mouse and with the keyboard: the menu opens, Regenerate
opens its dialog, Escape and Cancel close it, and nothing is submitted.
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

FIRST_MENU = "sl-dropdown.edit-menu >> nth=0"
FIRST_TRIGGER = FIRST_MENU + " >> [slot=trigger]"


@pytest.fixture
def system_page(browser, base_url, page_targets, admin_token):
    path, _needs_admin = page_targets["web.system"]
    context = browser.new_context(viewport={"width": 1280, "height": 900}, reduced_motion="reduce")
    context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.goto(base_url + path, wait_until="load")
    # Locators, not wait_for_function: the page's CSP forbids evaluating strings.
    page.locator("sl-dropdown.edit-menu:defined").first.wait_for()
    page.locator("sl-dialog.edit-dialog:defined").first.wait_for(state="attached")
    try:
        yield page, base_url + path
    finally:
        context.close()
    assert not errors, errors


def test_edit_menu_opens_a_dialog_with_the_mouse(system_page):
    page, url = system_page
    page.locator(FIRST_TRIGGER).click()
    page.get_by_role("menuitem", name="Regenerate").first.click()
    dialog = page.locator("sl-dialog.edit-dialog[open]")
    dialog.wait_for()
    assert dialog.count() == 1
    assert "Regenerate" in dialog.get_attribute("label")
    assert dialog.locator("input[name=edit_action]").input_value() == "regenerate"
    # Cancel closes it; the page was not submitted.
    dialog.locator("sl-button[data-dialog-close]").click()
    page.locator("sl-dialog.edit-dialog[open]").wait_for(state="detached")
    assert page.url == url


def test_edit_menu_works_from_the_keyboard_and_escape_closes_the_dialog(system_page):
    page, url = system_page
    page.evaluate("document.querySelector('sl-dropdown.edit-menu > sl-button').focus()")
    page.keyboard.press("Enter")
    page.locator(FIRST_MENU + " >> sl-menu").wait_for()
    page.keyboard.press("ArrowDown")
    assert page.evaluate("document.activeElement.tagName") == "SL-MENU-ITEM"
    page.keyboard.press("Enter")
    page.locator("sl-dialog.edit-dialog[open]").wait_for()
    page.keyboard.press("Escape")
    page.locator("sl-dialog.edit-dialog[open]").wait_for(state="detached")
    assert page.url == url


def test_the_delete_dialog_has_a_danger_button_and_the_forms_post_target(system_page):
    page, _url = system_page
    page.locator(FIRST_TRIGGER).click()
    page.get_by_role("menuitem", name="Delete").first.click()
    dialog = page.locator("sl-dialog.edit-dialog[open]")
    dialog.wait_for()
    assert dialog.locator("sl-button[variant=danger]").count() == 1
    assert dialog.locator("form input[name=edit_action]").input_value() == "delete"


@pytest.fixture
def sector_page(browser, base_url, page_targets, admin_token):
    path, _needs_admin = page_targets["web.sector"]
    context = browser.new_context(viewport={"width": 1280, "height": 900}, reduced_motion="reduce")
    context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.goto(base_url + path, wait_until="load")
    page.locator("sl-dropdown.admin-menu:defined").first.wait_for()
    page.locator("sl-dialog.edit-dialog:defined").first.wait_for(state="attached")
    try:
        yield page, base_url + path
    finally:
        context.close()
    assert not errors, errors


def test_the_sector_pages_admin_actions_are_one_menu_that_opens_dialogs(sector_page):
    """UX.26: no admin form is inline on the page; the Admin menu's items
    open them in a dialog, and Cancel and Escape close it without posting."""
    page, url = sector_page
    inline = page.locator("section#admin-sector form:visible")
    assert inline.count() == 0
    page.locator("sl-dropdown.admin-menu >> [slot=trigger]").click()
    page.get_by_role("menuitem", name="Generate neighborhood").click()
    dialog = page.locator("sl-dialog#admin-neighborhood[open]")
    dialog.wait_for()
    assert dialog.locator("form input[name=action]").input_value() == "generate_neighborhood"
    dialog.locator("sl-button[data-dialog-close]").click()
    page.locator("sl-dialog#admin-neighborhood[open]").wait_for(state="detached")
    page.locator("sl-dropdown.admin-menu >> [slot=trigger]").click()
    page.get_by_role("menuitem", name="Generate neighborhood").click()
    page.locator("sl-dialog#admin-neighborhood[open]").wait_for()
    page.keyboard.press("Escape")
    page.locator("sl-dialog#admin-neighborhood[open]").wait_for(state="detached")
    assert page.url == url
