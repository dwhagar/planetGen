"""
The Galaxy Map's mini map (MAP.75): a second, locked camera on the map's own renderer, drawn in the bottom right corner of the canvas.
It stays at its fit distance when the wheel zooms the main view, follows the direction the main camera looks from, and the Sector Map
has none. Skipped without Playwright, Chromium or MySQL, like the other map browser tests.
"""

import math

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_browser_maps import (  # noqa: E402,F401 -- fixtures
    GL_ARGS,
    VIEWPORT,
    admin_token,
    base_url,
    gl_browser,
    page,
    sample_job,
    sample_run,
    sample_params,
    site_app,
    site_db,
)
from tests.test_web_browser_sector_open import _ready  # noqa: E402
from tests.test_web_browser_sector_reach import CANVAS, _at_rest, _wheel  # noqa: E402

CANVAS_JS = "document.querySelector('#galaxymap3d-canvas')"


def _mini(page):
    return page.evaluate(f"() => {CANVAS_JS}.galaxyMiniMap()")


def test_the_mini_map_sits_in_the_corner_and_the_wheel_does_not_move_its_camera(page, base_url):
    page.goto(f"{base_url}/galaxy", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    page.wait_for_timeout(500)
    before = _mini(page)
    assert before and before["rect"], "the Galaxy Map has no mini map"
    box = page.locator(CANVAS).bounding_box()
    rect = before["rect"]
    assert rect["x"] + rect["side"] <= box["width"] and rect["y"] >= 0, rect
    assert rect["x"] > box["width"] / 2 and rect["side"] >= 96
    assert math.isclose(math.hypot(*before["eye"]), before["distance"], rel_tol=1e-6)
    # The pixels of the square are not the map behind it: the panel and the outline are drawn there.
    shot = page.screenshot(clip={"x": box["x"] + rect["x"], "y": box["y"] + box["height"] - rect["y"] - rect["side"], "width": rect["side"], "height": rect["side"]})
    assert len(set(shot[1000:5000:7])) > 4, "the mini map's square is blank"
    _wheel(page, (box["x"] + box["width"] / 3, box["y"] + box["height"] / 3), -240)
    _at_rest(page)
    after = _mini(page)
    assert math.isclose(after["distance"], before["distance"], rel_tol=1e-9), "the wheel zoomed the mini map"
    assert math.isclose(math.hypot(*after["eye"]), before["distance"], rel_tol=1e-6)
