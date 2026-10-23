# tests/test_web_ui_fields_browser.py

"""
The one form field (UX.124) in a real browser: a form of each kind (a text box, a number box, a select, a checkbox, the
datatable's "Starting at" box, covered in test_web_browser_tables.py) is filled and submitted on the search, nearby and system-list pages in light and dark and
at phone width, the pages pass the accessibility and phone-fit checks with the fields drawn, a validation error shows on
the right field, and with script off each form still posts through its plain twin. Needs Playwright with Chromium and
the MySQL test server, like test_web_a11y.py; skipped without them.

Run with: pytest tests/test_web_ui_fields_browser.py
"""

from urllib.parse import parse_qs, urlparse

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_a11y import (  # noqa: E402,F401 -- the browser, axe and site fixtures, and the checks
    SCHEMES, VIEWPORTS, _axe_blocking, _overflow, axe_source, base_url, browser, site_app, site_db,
)


def _context(browser, viewport="desktop", scheme="light", scripts=True):
    context = browser.new_context(viewport=VIEWPORTS[viewport], color_scheme=scheme, reduced_motion="reduce",
                                  java_script_enabled=scripts)
    page = context.new_page()
    problems = []
    page.on("console", lambda msg: msg.type == "error" and problems.append(f"console error: {msg.text}"))
    page.on("pageerror", lambda exc: problems.append(f"uncaught exception: {exc}"))
    return context, page, problems


def _query(page):
    return parse_qs(urlparse(page.url).query)


def _wait_defined(page, selector):
    """The field is drawn once Shoelace has defined it: wait for its inner control."""
    page.locator(f"{selector}:defined").wait_for(state="visible", timeout=20000)


@pytest.mark.parametrize("scheme", SCHEMES)
@pytest.mark.parametrize("viewport", sorted(VIEWPORTS))
def test_each_kind_of_field_is_filled_and_submitted(browser, axe_source, base_url, viewport, scheme):
    context, page, problems = _context(browser, viewport, scheme)
    try:
        # A text box and number boxes (the search page), with the page's own checks.
        page.goto(f"{base_url}/search", wait_until="load")
        _wait_defined(page, "sl-input[name='q']")
        page.locator("details.search-more > summary").click()
        problems += [f"/search: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]
        page.locator("sl-input[name='q'] input").fill("Sol")
        page.locator("sl-input[name='star_min_radius_km'] input").fill("1200")
        page.locator("#search-filters-form button[type=submit]").click()
        page.wait_for_url("**/search?*")
        assert _query(page)["q"] == ["Sol"] and _query(page)["star_min_radius_km"] == ["1200"]

        # A select, then a number box and checkboxes (nearby).
        page.goto(f"{base_url}/nearby", wait_until="load")
        _wait_defined(page, "sl-select[name='sector']")
        problems += [f"/nearby: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]
        page.locator("sl-select[name='sector']").click()
        page.locator("sl-select[name='sector'] sl-option").first.click()
        page.locator("main form button[type=submit]").first.click()
        page.wait_for_url("**/nearby?*")
        assert _query(page)["sector"]
        page.goto(f"{base_url}/nearby?place=0,0,0", wait_until="load")
        _wait_defined(page, "sl-input[name='distance']")
        problems += [f"/nearby asking: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]
        page.locator("sl-input[name='distance'] input").fill("12")
        page.locator("sl-checkbox[name='kinds']").first.click()
        page.locator("main form button[type=submit]").click()
        page.wait_for_url("**distance=12**")
        assert _query(page)["distance"] == ["12"] and len(_query(page)["kinds"]) == 1
    finally:
        context.close()
    assert not problems, f"at {viewport}/{scheme}:\n  " + "\n  ".join(problems)


def test_a_validation_error_shows_on_the_field_it_is_about(browser, axe_source, base_url):
    context, page, problems = _context(browser)
    try:
        page.goto(f"{base_url}/nearby?place=0,0,0&distance=500", wait_until="load")
        _wait_defined(page, "sl-input[name='distance']")
        error = page.locator("sl-input[name='distance'] .field-error")
        error.wait_for(state="visible")
        assert "largest distance" in error.inner_text()
        assert page.locator("sl-input[name='distance']").get_attribute("data-invalid") is not None
        assert page.locator("sl-checkbox .field-error").count() == 0
        assert not _axe_blocking(page, axe_source)
    finally:
        context.close()
    assert not problems, "\n".join(problems)


def test_without_script_each_form_posts_through_its_plain_twin(browser, base_url):
    context, page, problems = _context(browser, scripts=False)
    try:
        page.goto(f"{base_url}/search", wait_until="load")
        assert not page.locator("sl-input[name='q']").is_visible()
        page.get_by_label("Name", exact=True).first.fill("Sol")
        page.locator("#search-filters-form button[type=submit]").click()
        page.wait_for_url("**/search?*")
        assert _query(page)["q"] == ["Sol"]

        page.goto(f"{base_url}/nearby?place=0,0,0", wait_until="load")
        page.get_by_label("Distance in parsecs (up to 50)").fill("7")
        page.get_by_label("Star system").first.check()
        page.locator("main form button[type=submit]").click()
        page.wait_for_url("**distance=7**")
        assert _query(page)["distance"] == ["7"] and _query(page)["kinds"]

        page.goto(f"{base_url}/nearby", wait_until="load")
        page.get_by_role("combobox").first.select_option(index=0)
        page.locator("main form button[type=submit]").first.click()
        page.wait_for_url("**/nearby?*")
        assert _query(page)["sector"]
    finally:
        context.close()
    assert not problems, "\n".join(problems)
