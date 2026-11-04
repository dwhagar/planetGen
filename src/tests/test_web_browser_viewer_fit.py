"""
UX.142 (slice 1): the viewers of fixed shape stay inside the browser window. The nav chart (5:3) and the System Map diagram (1:1) used to be as tall as
their column was wide, so a wide window pushed them past the bottom of it; `boxfit.js` (`mapfit.fitWidth`) now caps their width so their height is what
the window leaves under the held header, as the Galaxy Map, the Sector Map and the 3D system view already are (MAP.350, MAP.359). They keep their shape
(the nav chart's labels are placed in percent of it), sit in the middle of their column, and re-fit when the window is resized, with no reload.

Sizes: 1280x720, 1366x768, 1440x900, 1920x1080, a 768x1024 tablet and a 390x844 phone, both orientations. Skipped without Playwright or Chromium.
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_browser_fixture_maps import gl_browser, page  # noqa: E402,F401 -- fixtures
from tests.test_web_browser_maps import _system_with_planets  # noqa: E402
from tests.test_web_a11y import admin_token, base_url, sample_job, sample_run, sample_params, site_app, site_db  # noqa: E402,F401 -- fixtures

FILL_SHARE = 0.95
MIN_HEIGHT_PX = 320
SIZES = [(1280, 720), (1366, 768), (1440, 900), (1920, 1080), (768, 1024), (1024, 768), (390, 844), (844, 390)]
SHAPES = {".navmap-viewport": 5 / 3, "#sysmap-viewport": 1.0}
RECT = """(sel) => { const e = document.querySelector(sel); const r = e.getBoundingClientRect();
  const h = document.querySelector('header'); const held = h && ['fixed', 'sticky'].includes(getComputedStyle(h).position);
  const col = e.parentElement.getBoundingClientRect();
  return {top: r.top + scrollY, left: r.left, width: r.width, height: r.height, win: innerHeight,
          covered: held ? h.getBoundingClientRect().bottom : 0, colLeft: col.left, colWidth: col.width}; }"""


def _limit(rect):
    """The tallest the box may be: 95% of what the window leaves under the held header, never under the least height."""
    return max(MIN_HEIGHT_PX, FILL_SHARE * (rect["win"] - rect["covered"]))


def _nav(page, base_url, site_app):
    items = site_app.test_client().get("/api/systems?limit=100").get_json()["items"]
    page.goto(f"{base_url}/nav?from={items[0]['id']}&to={items[1]['id']}", wait_until="load")
    page.wait_for_selector(".navmap-viewport")


def _check(rect, selector):
    assert rect["height"] <= _limit(rect) + 1.5, rect
    assert rect["width"] / rect["height"] == pytest.approx(SHAPES[selector], rel=0.02), rect      # the labels are in percent of this shape
    middle = rect["left"] + rect["width"] / 2
    assert middle == pytest.approx(rect["colLeft"] + rect["colWidth"] / 2, abs=2), rect          # in the middle of its column


@pytest.mark.parametrize("size", SIZES)
def test_the_nav_chart_is_never_taller_than_the_window_leaves(page, base_url, site_app, size):
    page.set_viewport_size({"width": size[0], "height": size[1]})
    _nav(page, base_url, site_app)
    page.wait_for_timeout(300)
    rect = page.evaluate(RECT, ".navmap-viewport")
    _check(rect, ".navmap-viewport")


@pytest.mark.parametrize("scheme", ["light", "dark"])
def test_the_nav_chart_re_fits_when_the_window_is_resized_in_either_theme(page, base_url, site_app, scheme):
    page.emulate_media(color_scheme=scheme)
    page.set_viewport_size({"width": 1600, "height": 1000})
    _nav(page, base_url, site_app)
    page.wait_for_timeout(300)
    wide = page.evaluate(RECT, ".navmap-viewport")
    page.set_viewport_size({"width": 1600, "height": 700})
    page.wait_for_function("(h) => document.querySelector('.navmap-viewport').getBoundingClientRect().height < h", arg=wide["height"] - 5)
    short = page.evaluate(RECT, ".navmap-viewport")
    _check(short, ".navmap-viewport")
    assert short["height"] < wide["height"]
    page.set_viewport_size({"width": 1600, "height": 1000})
    page.wait_for_function("(h) => document.querySelector('.navmap-viewport').getBoundingClientRect().height > h", arg=short["height"] + 5)


@pytest.mark.parametrize("size", [(1280, 720), (1920, 1080), (390, 844), (844, 390)])
def test_the_system_diagram_is_never_taller_than_the_window_leaves(page, base_url, site_app, size):
    page.set_viewport_size({"width": size[0], "height": size[1]})
    system_id, _scene = _system_with_planets(site_app.test_client())
    page.goto(f"{base_url}/system/{system_id}?view=3d", wait_until="load")
    page.wait_for_selector("#sysmap-view-diagram", state="visible")
    page.click("#sysmap-view-diagram")
    page.wait_for_selector("#sysmap-viewport", state="visible")
    page.wait_for_timeout(500)
    rect = page.evaluate(RECT, "#sysmap-viewport")
    _check(rect, "#sysmap-viewport")


def test_the_system_diagram_fits_when_it_is_shown_after_the_page_loaded_in_3d(page, base_url, site_app):
    """The diagram is hidden behind the Diagram/3D switch when the page loads; the box is fitted when the reader brings it up."""
    page.set_viewport_size({"width": 1920, "height": 1080})
    system_id, _scene = _system_with_planets(site_app.test_client())
    page.goto(f"{base_url}/system/{system_id}?view=3d", wait_until="load")
    page.wait_for_selector("#sysview3d-viewport", state="visible")
    page.click("#sysmap-view-diagram")
    page.wait_for_selector("#sysmap-viewport", state="visible")
    page.wait_for_timeout(500)
    _check(page.evaluate(RECT, "#sysmap-viewport"), "#sysmap-viewport")
