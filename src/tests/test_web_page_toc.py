"""
UX.195, UX.196, UX.204: every page has one h1 and no skipped level, an object page carries an "On this page" list built
from its headings, each entry points at a heading that exists, and following one into a collapsed section opens it.

Skipped without Playwright, Chromium or the MySQL test server, like `test_web_a11y.py`.
"""

import pytest

pytest.importorskip("playwright.sync_api")

from planetgen.api.authz import SESSION_COOKIE_NAME  # noqa: E402

from tests.test_web_a11y import (  # noqa: E402,F401 -- fixtures
    HEADING_AUDIT_JS,
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


def test_the_heading_audit_finds_a_skipped_level_a_second_h1_and_an_unheaded_table(browser):
    page = browser.new_page()
    try:
        page.set_content("<h1>a</h1><h3>skip</h3><section><table><tr><td>x</td></tr></table></section><h1>two</h1>")
        found = page.evaluate(HEADING_AUDIT_JS)
        assert any("one h1" in problem for problem in found)
        assert any("skips a level" in problem for problem in found)
        assert any("no heading" in problem for problem in found)
        page.set_content("<h1>a</h1><section><h2>T</h2><table><tr><td>x</td></tr></table></section>")
        assert page.evaluate(HEADING_AUDIT_JS) == []
    finally:
        page.close()


def test_an_object_page_has_a_contents_list_whose_links_reach_their_headings(browser, base_url, page_targets, admin_token):
    path, _needs_admin = page_targets["web.object_page"]
    context = browser.new_context(viewport={"width": 1280, "height": 900}, reduced_motion="reduce")
    context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
    page = context.new_page()
    try:
        page.goto(base_url + path, wait_until="load")
        toc = page.locator("nav.page-toc[aria-label='On this page']")
        toc.wait_for(state="attached")
        hrefs = toc.locator("a").evaluate_all("els => els.map(a => a.getAttribute('href'))")
        assert len(hrefs) >= 3
        for href in hrefs:
            assert page.locator(f"[id='{href[1:]}']").count() == 1, href
        # a link to a heading inside a collapsed property group opens that group
        inside = page.evaluate_handle("""() => [...document.querySelectorAll("nav.page-toc a")]
            .find(a => document.getElementById(a.getAttribute("href").slice(1)).closest("details.sheet-group"))""")
        assert inside.as_element() is not None, "no contents entry reaches a property group"
        target = inside.as_element().get_attribute("href")[1:]
        assert page.evaluate("id => document.getElementById(id).closest('details').open", target) is False
        inside.as_element().click()
        assert page.evaluate("id => document.getElementById(id).closest('details').open", target) is True
    finally:
        context.close()
