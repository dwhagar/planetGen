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


def _wait_until(page, condition, seconds=20):
    """Polls `condition` until it holds (the page's scripts run on animation frames, which a loaded machine delays)."""
    for _ in range(seconds * 10):
        if condition():
            return
        page.wait_for_timeout(100)
    raise AssertionError("the page did not reach the expected state")


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
        _wait_until(page, lambda: page.evaluate("id => document.getElementById(id).closest('details').open", target) is True)
    finally:
        context.close()


OWN_CONTENTS_JS = """() => {
    const main = document.getElementById("main");
    if (!main) return {auto: 0, own: []};
    const own = [...main.querySelectorAll("nav, [data-own-toc]")].filter(el => !el.closest(".page-toc")).filter(el => {
        const anchors = [...el.querySelectorAll("a[href^='#']")].filter(a => a.getAttribute("href").length > 1 && document.getElementById(a.getAttribute("href").slice(1)));
        return anchors.length >= 3;
    }).map(el => el.getAttribute("aria-label") || el.className || el.localName);
    return {auto: main.querySelectorAll(".page-toc").length, own: own};
}"""
"""A page's own contents lists (UX.235): a nav in <main> (or one flagged `data-own-toc`) with three or more links to
headings or sections on the same page, next to how many automatic "On this page" bars the page shows."""

PAGETOC_SOURCE = __import__("pathlib").Path(__file__).resolve().parents[1] / "html" / "static" / "pagetoc.js"


def _toc_page(browser, body):
    page = browser.new_page()
    page.set_content("<main id='main'><h1>Title</h1>" + body + "</main>")
    page.add_script_tag(path=str(PAGETOC_SOURCE))
    return page


HEADED = "<section><h2 id='a'>A</h2></section><section><h2 id='b'>B</h2></section><section><h2 id='c'>C</h2></section>"


def test_the_automatic_contents_bar_skips_a_page_that_has_its_own_contents_list_and_stays_on_one_without(browser):
    own = _toc_page(browser, "<nav aria-label='Groups'><a href='#a'>A</a> <a href='#b'>B</a> <a href='#c'>C</a></nav>" + HEADED)
    try:
        assert own.locator("nav.page-toc").count() == 0
    finally:
        own.close()
    flagged = _toc_page(browser, "<div data-own-toc><a href='#a'>A</a> <a href='#b'>B</a> <a href='#c'>C</a></div>" + HEADED)
    try:
        assert flagged.locator("nav.page-toc").count() == 0
    finally:
        flagged.close()
    # a nav that does not point inside the page (a site menu, the neighbours of an object) is not a contents list
    plain = _toc_page(browser, "<nav aria-label='Neighbours'><a href='/x'>X</a> <a href='/y'>Y</a> <a href='/z'>Z</a></nav>" + HEADED)
    try:
        assert plain.locator("nav.page-toc a").count() == 3
    finally:
        plain.close()
    bare = _toc_page(browser, HEADED)
    try:
        assert bare.locator("nav.page-toc a").count() == 3
    finally:
        bare.close()


def test_the_wide_map_page_gets_no_contents_menu_so_its_map_keeps_the_page_width(browser):
    """TEST.442: the Galaxy Map's page (main column `app-main-wide`) shows no "On this page" bar, which on a wide window is a
    column down the left that narrowed its map past the window (UX.215); a page without that class still gets one."""
    for classes, bars in (("app-main app-main-wide", 0), ("app-main", 1)):
        page = browser.new_page()
        try:
            page.set_content(f"<main id='main' class='{classes}'><h1>Title</h1>" + HEADED + "</main>")
            page.add_script_tag(path=str(PAGETOC_SOURCE))
            assert page.locator("nav.page-toc").count() == bars, classes
            assert page.evaluate("document.getElementById('main').classList.contains('has-page-toc')") is bool(bars)
        finally:
            page.close()


def test_no_page_shows_two_contents_lists(browser, base_url, page_targets, admin_token):
    """UX.235: over every page of the site, the automatic bar and a page's own contents list never both show."""
    context = browser.new_context(viewport={"width": 1280, "height": 900}, reduced_motion="reduce")
    context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
    page = context.new_page()
    two = []
    try:
        for endpoint, (path, _needs_admin) in sorted(page_targets.items()):
            if path is None:
                continue
            page.goto(base_url + path, wait_until="load")
            found = page.evaluate(OWN_CONTENTS_JS)
            if found["auto"] and found["own"]:
                two.append(f"{endpoint} {path}: automatic bar and {found['own']}")
        assert not two, "\n".join(two)
    finally:
        context.close()

def test_the_contents_menu_of_an_object_page_is_grouped_and_never_scrolls_sideways(browser, base_url, page_targets, admin_token):
    """The contents menu of an object page lists sections with their subsections nested under them: a column down the left on a
    wide window, a list under the title on a phone; neither scrolls sideways."""
    path, _needs_admin = page_targets["web.object_page"]
    for width, column in ((1280, True), (390, False)):
        context = browser.new_context(viewport={"width": width, "height": 900}, reduced_motion="reduce")
        context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
        page = context.new_page()
        try:
            page.goto(base_url + path, wait_until="load")
            toc = page.locator("nav.page-toc[aria-label='On this page']")
            toc.wait_for(state="attached")
            assert toc.locator("li ul li a").count() > 0, (width, "no subsections nested under their sections")
            menu = toc.bounding_box()
            first = page.locator("main .panel").first.bounding_box()
            if column:
                assert menu["x"] + menu["width"] <= first["x"] + 1, (width, menu, first)
            else:
                assert menu["y"] + menu["height"] <= first["y"] + 1, (width, menu, first)
            assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth + 1"), "the page scrolls sideways"
            assert page.evaluate("n => n.scrollWidth <= n.clientWidth + 1", toc.element_handle()), "the menu scrolls sideways"
        finally:
            context.close()


def test_the_contents_menu_is_one_row_when_it_fits_and_the_grouped_list_when_it_does_not(browser):
    """The links go in one row across the top only if all of them fit without scrolling; otherwise the grouped list is used."""
    few = "<h2>One</h2><h3>Sub</h3><h2>Two</h2>"
    many = "".join(f"<h2>Section number {i}</h2><h3>Subsection {i}</h3>" for i in range(30))
    for body, strip in ((few, True), (many, False)):
        page = browser.new_page(viewport={"width": 1280, "height": 800})
        try:
            page.set_content("<main id='main' class='app-main'><h1>Title</h1>" + body + "</main>")
            page.add_script_tag(path=str(PAGETOC_SOURCE))
            assert page.evaluate("document.getElementById('main').classList.contains('toc-strip')") is strip
            assert page.evaluate("n => n.scrollWidth <= n.clientWidth + 1", page.locator("nav.page-toc").element_handle())
        finally:
            page.close()
