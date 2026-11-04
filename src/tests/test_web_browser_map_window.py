"""
MAP.350: the Galaxy Map window fills the browser window: the map is always 95% of the window's height (less a held header), recomputed on
resize, the page rests with the map centred in the window, and once the reader scrolls the page it is left where they put it. A short window
keeps the least height.

Skipped without Playwright, Chromium or WebGL, like the other map browser tests.
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_browser_fixture_maps import GALAXY_CANVAS, _open_galaxy, _settle, gl_browser, map_site, page  # noqa: E402,F401 -- fixtures

FILL_SHARE = 0.95
MIN_HEIGHT_PX = 320
RECT = """() => { const r = document.querySelector('.galaxymap3d-panel .starmap-viewport').getBoundingClientRect();
    const h = document.querySelector('header'); const sticky = h && ['fixed', 'sticky'].includes(getComputedStyle(h).position);
    return {top: r.top, bottom: r.bottom, height: r.height, inner: innerHeight, scrollY: scrollY, covered: sticky ? h.getBoundingClientRect().bottom : 0}; }"""


def _rested(page, expected_inner):
    page.wait_for_function(f"() => innerHeight === {expected_inner}")
    last = None
    for _ in range(30):
        page.wait_for_timeout(200)
        now = page.evaluate(RECT)
        if last and abs(now["top"] - last["top"]) < 0.5 and abs(now["height"] - last["height"]) < 0.5:
            return now
        last = now
    return last


@pytest.mark.parametrize("size", [(1920, 950), (1280, 720)])
def test_the_map_is_ninety_five_percent_of_the_window(page, map_site, size):
    page.set_viewport_size({"width": size[0], "height": size[1]})
    _open_galaxy(page, map_site)
    _settle(page)
    rect = _rested(page, size[1])
    assert rect["covered"] > 0, "the fixture's page holds its header at the top of the window"
    assert abs(rect["height"] - FILL_SHARE * (size[1] - rect["covered"])) <= 1, rect
    margin = (size[1] - rect["covered"] - rect["height"]) / 2
    assert abs(rect["top"] - (rect["covered"] + margin)) <= 2, rect
    assert abs(rect["bottom"] - (size[1] - margin)) <= 3, rect
    assert page.locator(GALAXY_CANVAS).bounding_box()["height"] > rect["height"] - 2, "the canvas fills the window"


def test_a_short_window_keeps_the_least_height_and_the_page_still_rests_with_the_map_at_the_top(page, map_site):
    page.set_viewport_size({"width": 1200, "height": 330})
    _open_galaxy(page, map_site)
    _settle(page)
    rect = _rested(page, 330)
    assert rect["height"] == pytest.approx(MIN_HEIGHT_PX, abs=1), rect
    assert rect["top"] >= rect["covered"] - 2, rect


def test_a_page_the_reader_has_scrolled_is_left_where_they_put_it(page, map_site):
    page.set_viewport_size({"width": 1600, "height": 900})
    _open_galaxy(page, map_site)
    _settle(page)
    _rested(page, 900)
    page.evaluate("() => window.scrollTo(0, 0)")
    page.wait_for_timeout(300)
    page.set_viewport_size({"width": 1600, "height": 880})
    rect = _rested(page, 880)
    assert rect["scrollY"] == 0 and rect["top"] > rect["covered"] + 50, rect
    assert abs(rect["height"] - FILL_SHARE * (880 - rect["covered"])) <= 1, rect
