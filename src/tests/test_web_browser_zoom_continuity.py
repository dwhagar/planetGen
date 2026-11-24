"""
MAP.257: one continuous scroll from the whole galaxy to a sector and out again, in a real browser on the fixture site: every wheel tick takes the camera the same steady factor nearer or
farther, however the stage under it changes (a sector opening, a sector or a system closing), so leaving an opened sector no longer snaps the camera to the whole galaxy.
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.map_site_support import SECTOR_ID, SECTORS, _sector_center_pc, map_site  # noqa: E402,F401 -- fixture
from tests.test_web_browser_fixture_maps import GALAXY_CANVAS, _open, gl_browser, page  # noqa: E402,F401 -- fixtures

CAMERA = "(() => { const k = document.querySelector('#galaxymap3d-canvas').galaxyCamera(); return k ? k.dist : null; })()"
SCREEN = "(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)"
WHEEL_PX = 300
"""int: One wheel tick, pixels."""
STEADY_STEP = 1.66
"""float: The most a tick may move the camera out (the wheel's own step is e to the half, 1.649)."""


def _dist(page):
    return page.evaluate(CAMERA)


def test_leaving_an_opened_sector_carries_the_zoom_on_by_the_wheels_own_step(page, map_site):
    center = _sector_center_pc(*next(s[1:4] for s in SECTORS if s[0] == SECTOR_ID))
    _open(page, f"{map_site}/galaxy", GALAXY_CANVAS)
    page.wait_for_function("() => { const c = document.querySelector('#galaxymap3d-canvas'); return !c.galaxyReady || c.galaxyReady(); }", timeout=60000)
    page.wait_for_timeout(1500)
    box = page.locator(GALAXY_CANVAS).bounding_box()
    middle = (box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    for _ in range(60):
        spot = page.evaluate(SCREEN, list(center))
        aim = (spot["x"], spot["y"]) if spot and 0 < spot["x"] < 1000 and 0 < spot["y"] < 800 else middle
        page.mouse.move(*aim)
        page.mouse.wheel(0, -WHEEL_PX)
        page.wait_for_timeout(150)
        if "open=1" in page.url:
            break
    else:
        raise AssertionError(f"60 wheel ticks did not open the sector: {page.url}")
    page.wait_for_timeout(400)
    farthest = _dist(page)
    for _ in range(8):                       # a little nearer in the opened sector, as a visitor reading it would
        page.mouse.move(*middle)
        page.mouse.wheel(0, -WHEEL_PX)
        page.wait_for_timeout(150)
    previous = _dist(page)
    assert previous < farthest, (previous, farthest)
    closed_once = False
    for tick in range(120):                  # out until the camera is back at the whole galaxy
        page.mouse.move(*middle)
        page.mouse.wheel(0, WHEEL_PX)
        page.wait_for_timeout(150)
        now = _dist(page)
        assert now <= previous * STEADY_STEP, f"tick {tick}: the camera went from {previous:.6g} pc to {now:.6g} pc in one wheel tick (url {page.url})"
        closed_once = closed_once or "open=1" not in page.url
        previous = now
        if now > 30000:
            break
    assert closed_once, "the sector never closed on the way out"
    assert previous > 30000, f"the camera did not get back to the whole galaxy: {previous} pc"
