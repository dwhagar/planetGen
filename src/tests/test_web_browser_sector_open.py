"""
Opening a sector from the Galaxy Map, on the real generated database, with the
flight playing (MAP.209, MAP.210).

The other map browser tests run with reduced motion, where a flight takes no
time. Here the page keeps its animations, so the camera really flies from the
sector's window into the sector; the move must then be one continuous view (no
page load, one URL, the panel beside the map never empty) and what is drawn in
the opened sector must be selectable: a click on a star's projected position
selects it and fills the panel with its facts.

Runs against `test_web_a11y.py`'s small generated database like
`test_web_browser_maps.py`; skipped without Playwright, Chromium or MySQL.
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from planetgen.web.sector_page import sector_designation  # noqa: E402
from tests.test_web_browser_maps import (  # noqa: E402,F401 -- fixtures
    GL_ARGS,
    VIEWPORT,
    admin_token,
    base_url,
    gl_browser,
    sample_job,
    sample_job_tree,
    sample_params,
    site_app,
    site_db,
)

CANVAS = "#galaxymap3d-canvas"
INFO = "#galaxymap3d-info"

WATCH = """
window.__sameDocument = true;
window.__panel = [];
new MutationObserver(() => {
  const info = document.querySelector('#galaxymap3d-info');
  window.__panel.push([location.search, info ? info.innerText.trim().slice(0, 60) : null]);
}).observe(document, {childList: true, subtree: true, characterData: true});
"""
"""str: Run in the page before anything else: marks the document (a reload makes a new one without the mark)
and records the panel's text each time anything in the page changes."""


@pytest.fixture()
def moving_page(gl_browser):
    """A page whose animations play (reduced motion left off), unlike `page`."""
    context = gl_browser.new_context(viewport=VIEWPORT, device_scale_factor=1, reduced_motion="no-preference",
                                     color_scheme="dark")
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.add_init_script(WATCH)
    page.errors = errors
    try:
        yield page
    finally:
        context.close()
    assert not errors, f"uncaught exceptions: {errors}"


def _ready(page):
    """Until the map has drawn what its URL names and the camera holds still."""
    page.wait_for_selector("#galaxymap3d-steps [data-steps-panel] li", state="attached")
    page.wait_for_function("""() => {
        const canvas = document.querySelector("#galaxymap3d-canvas");
        return !canvas || !canvas.galaxyReady || canvas.galaxyReady();
    }""", timeout=30000)
    last = None
    for _ in range(60):
        now = page.evaluate("() => (document.querySelector('#galaxymap3d-scale') || {}).innerHTML || ''")
        if now == last:
            return
        last = now
        page.wait_for_timeout(250)


def _sector(site_app, sector_id):
    detail = site_app.test_client().get(f"/api/sectors/{sector_id}").get_json()
    assert detail["placed"], detail
    return detail


def _designation(detail):
    """The designation the Galaxy Map names the sector by, worked out in the test process from the sector's detail.
    It used to be read from the redirect of the "Show on Galaxy Map" link, a request that waited behind the page's
    own tile requests and timed out at 30 s under load (TEST.195); that redirect is tested in test_web_sector_nav."""
    return sector_designation(detail)


def _scene(site_app, sector_id):
    """The sector's scene JSON, asked of the app in the test process. Through the served site it waited behind the
    page's own requests and the software renderer's load and timed out at 30 s in the gate (TEST.272)."""
    response = site_app.test_client().get(f"/sector/{sector_id}/scene")
    assert response.status_code == 200, response.status_code
    return response.get_json()


def _star_spots(page, site_app, sector_id):
    """Where the opened sector's stars fall on the page: [(name, (x, y))], from the sector's scene JSON."""
    scene = _scene(site_app, sector_id)
    unit = scene["halfEdgePc"] / (scene["sceneHalfPx"] or 160)
    cx, cy, cz = scene["centerPc"]
    box = page.locator(CANVAS).bounding_box()
    spots = []
    for star in scene["stars"]:
        if not star.get("endpoint", "").startswith("system:"):
            continue
        point = (cx + star["x"] * unit, cy - star["y"] * unit, cz + star["z"] * unit)  # the scene flips y
        spot = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(point))
        if spot and box["x"] + 4 < spot["x"] < box["x"] + box["width"] - 4 and box["y"] + 4 < spot["y"] < box["y"] + box["height"] - 4:
            spots.append((star["name"], (spot["x"], spot["y"])))
    return spots


def _picks_a_star(page, spots):
    """Clicks each star's spot until the panel names a star; returns its name."""
    info = page.locator(INFO)
    names = {name for name, _spot in spots}
    for name, (x, y) in spots:
        page.mouse.move(x, y)
        page.wait_for_timeout(100)
        page.mouse.click(x, y)
        page.wait_for_timeout(200)
        heading = info.locator("h3").inner_text() if info.locator("h3").count() else None
        if heading in names:
            return heading
    return None


def test_opening_a_sector_is_one_continuous_view_with_a_filled_panel(moving_page, base_url, sample_params, site_app):
    """MAP.209: a click on a charted sector (the flight playing) opens it in place: the document is the same one
    (no reload), the address changes once to `?sector=...&open=1` and stays, the panel never goes empty and ends
    describing the opened sector."""
    page = moving_page
    sector_id = sample_params["sector_id"]
    detail = _sector(site_app, sector_id)
    code = _designation(detail)
    page.goto(f"{base_url}/galaxy?sector={code}", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    page.evaluate("window.__panel.length = 0")
    point = (detail["center_x_pc"], detail["center_y_pc"], detail["center_z_pc"])
    spot = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(point))
    assert spot, "the sector's centre is not on the map"
    page.mouse.move(spot["x"], spot["y"])
    page.wait_for_timeout(200)
    page.mouse.click(spot["x"], spot["y"])
    page.wait_for_function("() => location.search.includes('open=1')", timeout=20000)
    _ready(page)
    assert page.evaluate("window.__sameDocument") is True, "the page reloaded while the sector opened"
    assert page.url.endswith(f"?sector={code}&open=1"), page.url
    panel = page.evaluate("window.__panel")
    assert all(text for _query, text in panel if text is not None), f"the panel went empty: {panel}"
    heading = page.locator(f"{INFO} h3").inner_text()
    assert heading == detail["name"], heading
    assert page.locator(INFO).inner_text().strip(), "the panel is empty"
    assert any(page.eval_on_selector_all(f"{INFO} a", "els => els.map(e => e.textContent)")), "the panel offers no link"


def test_a_star_in_an_opened_sector_is_selectable(moving_page, base_url, sample_params, site_app):
    """MAP.210: in the opened sector (loaded straight from its address, as a reload does, with the flight
    playing) a click on a star's projected position selects it: the panel names the star and links its system."""
    page = moving_page
    sector_id = sample_params["sector_id"]
    detail = _sector(site_app, sector_id)
    code = _designation(detail)
    page.goto(f"{base_url}/galaxy?sector={code}&open=1", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    spots = _star_spots(page, site_app, sector_id)
    assert spots, "no star of the sector is on the map"
    name = _picks_a_star(page, spots)
    assert name, f"no click on a star of {len(spots)} picked one: {page.locator(INFO).inner_text()!r}"
    assert page.locator(INFO).locator("a", has_text="View system").count() == 1
    assert page.evaluate("window.__sameDocument") is True
