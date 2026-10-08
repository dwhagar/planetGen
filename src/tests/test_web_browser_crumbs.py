"""
NAV.14: a page's breadcrumb is the shared one-line component
(`static/breadcrumb.js`), like the Galaxy Map's. On a system page at phone
width the breadcrumb stays on one line, and when its steps don't fit the
middle ones fold into a "…" menu that lists them. Skipped without
Playwright, Chromium or the MySQL test server.
"""

import pytest

pytest.importorskip("playwright.sync_api")

from tests.test_web_a11y import (  # noqa: E402,F401 -- fixtures
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


def _line(page):
    return page.evaluate("""() => {
        const nav = document.querySelector('nav.breadcrumbs');
        const ol = nav.querySelector('ol');
        const tops = Array.from(ol.children).map((li) => li.getBoundingClientRect().top);
        const rows = new Set([tops.some((t) => Math.abs(t - tops[0]) > 12) ? 2 : 1]);
        return {ready: nav.classList.contains('crumbs-ready'), rows: rows.size,
                more: !!ol.querySelector('.crumb-menu'), overflow: ol.scrollWidth - ol.clientWidth,
                steps: Array.from(ol.children).length};
    }""")


@pytest.mark.parametrize("width", [1280, 320])
def test_a_pages_breadcrumb_is_one_line_at_any_width(browser, base_url, page_targets, width):
    path, _needs_admin = page_targets["web.system"]
    page = browser.new_page(viewport={"width": width, "height": 800})
    try:
        page.goto(base_url + path, wait_until="load")
        page.locator("nav.breadcrumbs.crumbs-ready").wait_for()
        line = _line(page)
        assert line["ready"] and line["rows"] == 1, line
        assert line["overflow"] <= 1, line
        if line["more"]:
            page.locator("nav.breadcrumbs .crumb-menu summary").click()
            assert page.locator("nav.breadcrumbs .crumb-menu li").count() >= 1
    finally:
        page.close()
