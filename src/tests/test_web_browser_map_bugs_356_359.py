"""
Boss's four map bugs of 10-30 (MAP.356 to MAP.359): the Key is a small box in the map's corner; every charted region is a point of interest
with a link that goes to it (regions too close to tell apart are one point); the sector card's buttons stay in view on a small map; and the
sector page's map and the 3D system view fill the browser window as the Galaxy Map does.

Skipped without Playwright, Chromium or WebGL, like the other map browser tests.
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_browser_fixture_maps import (  # noqa: E402,F401 -- fixtures
    GALAXY_CANVAS, SECTOR_CANVAS, _open, _open_galaxy, far_site, _open_sector, _poi_layers, _settle, gl_browser, map_site, page,
)
from tests.test_web_browser_maps import _system_with_planets  # noqa: E402
from tests.test_web_browser_sector_open import SCENE_WAIT_MS, _designation, _sector, _settled, moving_page  # noqa: E402,F401
from tests.test_web_a11y import admin_token, base_url, sample_job, sample_run, sample_params, site_app, site_db  # noqa: E402,F401 -- fixtures

FILL_SHARE = 0.95
KEY = "#galaxymap3d-layer-legend"
CANVAS = "document.querySelector('#galaxymap3d-canvas')"


def _view_box(page, selector=".starmap-viewport"):
    return page.evaluate("(sel) => { const r = document.querySelector(sel).getBoundingClientRect(); return {width: r.width, height: r.height, top: r.top}; }", selector)


def test_the_key_is_a_small_box_in_the_corner_of_the_map(page, map_site):
    """MAP.356: the Key is one column of short lines in a box no wider than two fifths of the map (13 rem at most) and no taller than
    two fifths of it, in the map's lower left, so it never covers much of the map."""
    page.set_viewport_size({"width": 1280, "height": 800})
    _open_galaxy(page, map_site)
    page.wait_for_function(f"() => {CANVAS}.__mapLayers.keys().some((k) => k.drawn)", timeout=20000)
    page.wait_for_timeout(1000)
    got = page.evaluate("""() => {
        const key = document.querySelector('#galaxymap3d-layer-legend').getBoundingClientRect();
        const view = document.querySelector('.galaxymap3d-panel .starmap-viewport').getBoundingClientRect();
        const lines = getComputedStyle(document.querySelector('#galaxymap3d-layer-legend .map-key-lines'));
        return {key: [key.width, key.height], view: [view.width, view.height], direction: lines.flexDirection,
                left: key.left - view.left, bottom: view.bottom - key.bottom};
    }""")
    assert got["direction"] == "column", got
    assert got["key"][0] <= 0.4 * got["view"][0] + 1 and got["key"][0] <= 13 * 16 + 2, got
    assert got["key"][1] <= 0.4 * got["view"][1] + 1, got
    assert got["key"][0] * got["key"][1] <= 0.08 * got["view"][0] * got["view"][1], got
    assert got["left"] < 30 and got["bottom"] < 120, got


POIS = "() => document.querySelector('#galaxymap3d-canvas').galaxyRegionPois()"


def test_a_charted_region_is_a_point_of_interest_with_a_link_that_goes_to_it(page, map_site):
    """MAP.357: from the galactic view each region with generated sectors has its point of interest, inside the map, and its link flies the
    camera to the region (the camera comes closer and stays within the region's reach)."""
    page.set_viewport_size({"width": 1280, "height": 800})
    _open_galaxy(page, map_site)
    page.wait_for_function(f"() => {CANVAS}.__mapRegions && {CANVAS}.__mapRegions.drawn() > 0", timeout=20000)
    page.wait_for_function(f"() => {CANVAS}.galaxyRegionPois().length > 0", timeout=20000)
    pois = page.evaluate(POIS)
    box = _view_box(page, "#galaxymap3d-canvas")
    for poi in pois:
        assert 0 <= poi["x"] <= box["width"] and 0 <= poi["y"] <= box["height"], poi
        assert poi["count"] == len(poi["keys"]) >= 1, poi
    # Every region with sectors is in exactly one point: the points hold all the regions the layer lists.
    listed = page.evaluate(f"() => {CANVAS}.galaxyRegionSpots().length")
    assert sum(p["count"] for p in pois) <= listed
    before = page.evaluate(f"() => {CANVAS}.galaxyCamera().dist")
    page.locator(".map-region-poi").first.click()
    page.wait_for_function(f"() => !{CANVAS}.galaxyMoving()", timeout=60000)
    _settle(page)
    after = page.evaluate(f"() => {CANVAS}.galaxyCamera().dist")
    assert after < before, (before, after)


def test_the_goto_box_lists_the_charted_places_when_it_is_empty(page, map_site):
    """MAP.357: a blank Go to box offers the regions with the most generated sectors, each a link."""
    page.set_viewport_size({"width": 1280, "height": 800})
    _open_galaxy(page, map_site)
    page.wait_for_function(f"() => {CANVAS}.galaxyRegionPois().length > 0", timeout=20000)
    page.locator("#galaxymap3d-address-input").evaluate("i => { i.closest('form').hidden = false; i.focus(); }")
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-matches').hidden", timeout=10000)
    text = page.locator("#galaxymap3d-matches").inner_text()
    assert "charted places" in text.lower(), text


@pytest.mark.parametrize("width,height", [(1000, 380), (800, 340)])
def test_the_sector_card_keeps_its_buttons_in_view_on_a_small_map(page, map_site, width, height):
    """MAP.358: on a small map the card still holds its name, values and buttons (Expand and the rest) inside the map and the card, and a
    button can be reached."""
    page.set_viewport_size({"width": width, "height": height})
    _open_sector(page, map_site)
    page.locator(".starmap-sr-list button").first.dispatch_event("click")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    page.wait_for_timeout(600)
    got = page.evaluate("""() => {
        const info = document.querySelector('#galaxymap3d-info');
        const view = document.querySelector('#galaxymap3d-canvas').closest('.starmap-viewport').getBoundingClientRect();
        const box = info.getBoundingClientRect();
        const buttons = [...info.querySelectorAll(':scope > .map-info-actions > *, :scope .map-info-expand')].filter((e) => e.getClientRects().length);
        return {buttons: buttons.map((b) => { const r = b.getBoundingClientRect(); return [b.textContent.trim(), r.top >= box.top - 1 && r.bottom <= box.bottom + 1 && r.top >= view.top - 1 && r.bottom <= view.bottom + 1]; }),
                height: box.height, view: view.height};
    }""")
    assert got["buttons"] and all(inside for _, inside in got["buttons"]), got


@pytest.mark.parametrize("size", [(1280, 800), (390, 760)])
def test_the_sector_page_map_fills_the_window_as_the_galaxy_map_does(page, map_site, size):
    """MAP.359: the sector page's map is as tall as the Galaxy Map is for the browser window (95 percent of it less a held header)."""
    page.set_viewport_size({"width": size[0], "height": size[1]})
    _open_sector(page, map_site)
    page.wait_for_timeout(800)
    box = _view_box(page)
    covered = page.evaluate("() => { const h = document.querySelector('header'); return h && ['fixed', 'sticky'].includes(getComputedStyle(h).position) ? h.getBoundingClientRect().bottom : 0; }")
    assert box["height"] == pytest.approx(max(320, FILL_SHARE * (size[1] - covered)), abs=2), (box, covered)


def test_the_3d_system_view_fills_the_window_as_the_galaxy_map_does(page, base_url, site_app):
    """MAP.359: the 3D system view is as tall as the Galaxy Map is for the browser window."""
    page.set_viewport_size({"width": 1280, "height": 800})
    system_id, _scene = _system_with_planets(site_app.test_client())
    page.goto(f"{base_url}/system/{system_id}", wait_until="load")
    page.evaluate("window.localStorage.clear()")
    page.goto(f"{base_url}/system/{system_id}", wait_until="load")
    page.wait_for_selector('#sysview3d[data-ready="true"]')
    page.wait_for_timeout(800)
    box = _view_box(page, ".sysview-viewport")
    covered = page.evaluate("() => { const h = document.querySelector('header'); return h && ['fixed', 'sticky'].includes(getComputedStyle(h).position) ? h.getBoundingClientRect().bottom : 0; }")
    assert box["height"] == pytest.approx(max(320, FILL_SHARE * (800 - covered)), abs=2), (box, covered)


def test_an_opened_sector_fills_the_window_to_within_two_percent(moving_page, base_url, site_app, sample_params):
    """MAP.360: the sector opened in place is fitted to the window with a 2 percent margin (its corners reach nearly to the edge
    of the shorter side)."""
    detail = _sector(site_app, sample_params["sector_id"])
    moving_page.goto(f"{base_url}/galaxy?sector={_designation(detail)}&open=1", wait_until="load")
    moving_page.wait_for_selector(GALAXY_CANVAS, state="attached")
    centre = (detail["center_x_pc"], detail["center_y_pc"], detail["center_z_pc"])
    _settled(moving_page, centre)
    moving_page.wait_for_function("() => { const o = document.querySelector('#galaxymap3d-canvas').galaxyOpenSector(); return !!o && o.share >= 0.999; }", timeout=SCENE_WAIT_MS)
    moving_page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=60000)
    moving_page.wait_for_timeout(500)
    share = moving_page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxySectorFitShare()")
    assert share and 0.85 <= max(share["x"], share["y"]) <= 1.0 and min(share["x"], share["y"]) <= 1.0, share


def test_superbubbles_are_off_until_the_menu_button_shows_them(page, far_site):
    """MAP.362: the superbubbles are not drawn when the map opens; the Menu's Superbubbles button (and `?show=superbubbles`) draws them and the
    button turns them off again."""
    _open_galaxy(page, far_site, "")
    look = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyBubbleLook()")
    assert look and look["count"] == 1 and look["shown"] is False, look
    page.click("#galaxymap3d-menu > summary")
    button = page.locator('#galaxymap3d-kinds button[data-layer="superbubbles"], #galaxymap3d-controls button[data-layer="superbubbles"]').first
    button.click()
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyBubbleLook().shown") is True
    button.click()
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyBubbleLook().shown") is False
    _open_galaxy(page, far_site, "?show=superbubbles")
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyBubbleLook().shown") is True
