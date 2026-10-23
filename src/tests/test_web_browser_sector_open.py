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

from planetgen.util import draw  # noqa: E402
from planetgen.web.sector_page import sector_designation  # noqa: E402
from tests.bughunt_support import forced_system_config, mysql_argv, run_cli  # noqa: E402
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


def _settled(page, centre):
    """`_ready`, and then until `centre` (a point of the scene in parsecs, the opened sector's centre) stays where it is
    on the screen for a second, with the canvas scrolled into the window. The scale text `_ready` reads is the same
    while a flight is still on (a sector far from the core flies longer), and a click on a moving map is ignored."""
    _ready(page)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    held, steady = None, 0
    for _ in range(100):
        spot = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(centre))
        now = None if spot is None else (round(spot["x"]), round(spot["y"]))
        steady = steady + 1 if now == held else 0
        held = now
        if steady >= 5:
            return
        page.wait_for_timeout(200)


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


def _star_points(site_app, sector_id):
    """The opened sector's stars: [(name, (x, y, z))], the point each is drawn at in the map's frame (parsecs), from the
    sector's scene JSON."""
    scene = _scene(site_app, sector_id)
    unit = scene["halfEdgePc"] / (scene["sceneHalfPx"] or 160)
    cx, cy, cz = scene["centerPc"]
    return [(star["name"], (cx + star["x"] * unit, cy - star["y"] * unit, cz + star["z"] * unit))  # the scene flips y
            for star in scene["stars"] if star.get("endpoint", "").startswith("system:")]


def _spot(page, point):
    """Where `point` falls on the page now, if that is inside the canvas and the window: (x, y) or None."""
    spot = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(point))
    if not spot:
        return None
    box = page.locator(CANVAS).bounding_box()
    window = page.viewport_size
    left, top = max(box["x"], 0) + 4, max(box["y"], 0) + 4
    right, bottom = min(box["x"] + box["width"], window["width"]) - 4, min(box["y"] + box["height"], window["height"]) - 4
    return (spot["x"], spot["y"]) if left < spot["x"] < right and top < spot["y"] < bottom else None


PICK_WAIT_MS = 3000
"""int: How long a click on a star may take to name it in the panel."""

PICK_TRIES = 4
"""int: How many times a star is clicked, where it is drawn then, before the next star is tried."""


def _steady_spot(page, point, samples=4):
    """Where `point` falls on the page once it has held the same place for `samples` looks 150 ms apart (rounded to a
    pixel), or its last place if it keeps moving: after the flight ends the map can still settle, and on a loaded
    machine a star moved 100 px between being worked out and being clicked (TEST.291)."""
    held, steady, spot = None, 0, None
    for _ in range(40):
        spot = _spot(page, point)
        now = None if spot is None else (round(spot[0]), round(spot[1]))
        steady = steady + 1 if now is not None and now == held else 0
        held = now
        if steady >= samples:
            break
        page.wait_for_timeout(150)
    return spot


def _picks_a_star(page, points):
    """Clicks each star where it is drawn (worked out again at each click, once it holds still) until the panel names a
    star; returns its name."""
    info = page.locator(INFO)
    names = {name for name, _point in points}
    for _name, point in points:
        for _ in range(PICK_TRIES):
            spot = _steady_spot(page, point)
            if not spot:
                break
            x, y = spot
            page.mouse.move(x, y)
            page.wait_for_timeout(100)
            page.mouse.click(x, y)
            # The panel names the star once the page has handled the click; on a loaded machine that takes longer than
            # the 200 ms this used to wait, and the click was then taken for a miss (TEST.291).
            try:
                page.wait_for_function("""([selector, names]) => {
                    const heading = document.querySelector(selector + " h3");
                    return !!heading && names.includes(heading.innerText.trim());
                }""", arg=[INFO, sorted(names)], timeout=PICK_WAIT_MS)
            except sync_api.TimeoutError:
                continue
            return info.locator("h3").inner_text()
    return None


def test_opening_a_sector_is_one_continuous_view_with_a_filled_panel(moving_page, base_url, sample_params, site_app):
    """MAP.209: a click on a charted sector (the flight playing) opens it in place: the document is the same one
    (no reload), the address changes once to `?sector=...&open=1` and stays, the panel never goes empty and ends
    describing the opened sector."""
    page = moving_page
    sector_id = sample_params["sector_id"]
    detail = _sector(site_app, sector_id)
    code = _designation(detail)
    point = (detail["center_x_pc"], detail["center_y_pc"], detail["center_z_pc"])
    page.goto(f"{base_url}/galaxy?sector={code}", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _settled(page, point)
    page.evaluate("window.__panel.length = 0")
    spot = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(point))
    assert spot, "the sector's centre is not on the map"
    page.mouse.move(spot["x"], spot["y"])
    page.wait_for_timeout(200)
    page.mouse.click(spot["x"], spot["y"])
    page.wait_for_function("() => location.search.includes('open=1')", timeout=20000)
    _settled(page, point)
    assert page.evaluate("window.__sameDocument") is True, "the page reloaded while the sector opened"
    assert page.url.endswith(f"?sector={code}&open=1"), page.url
    panel = page.evaluate("window.__panel")
    assert all(text for _query, text in panel if text is not None), f"the panel went empty: {panel}"
    heading = page.locator(f"{INFO} h3").inner_text()
    assert heading == detail["name"], heading
    assert page.locator(INFO).inner_text().strip(), "the panel is empty"
    assert any(page.eval_on_selector_all(f"{INFO} a", "els => els.map(e => e.textContent)")), "the panel offers no link"


def _selects_a_star(page, base_url, site_app, detail):
    """Opens `detail`'s sector from its address and clicks its stars where they are drawn until one is selected."""
    sector_id = detail["id"]
    page.goto(f"{base_url}/galaxy?sector={_designation(detail)}&open=1", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _settled(page, (detail["center_x_pc"], detail["center_y_pc"], detail["center_z_pc"]))
    points = _star_points(site_app, sector_id)
    assert points, "no star of the sector is on the map"
    name = _picks_a_star(page, points)
    assert name, f"no click on a star of {len(points)} picked one: {page.locator(INFO).inner_text()!r}"
    assert page.locator(INFO).locator("a", has_text="View system").count() == 1
    assert page.evaluate("window.__sameDocument") is True


def test_a_star_in_an_opened_sector_is_selectable(moving_page, base_url, sample_params, site_app):
    """MAP.210: in the opened sector (loaded straight from its address, as a reload does, with the flight
    playing) a click on a star's projected position selects it: the panel names the star and links its system."""
    _selects_a_star(moving_page, base_url, site_app, _sector(site_app, sample_params["sector_id"]))


FAR_SHELL = 2000
"""int: The shell of the far sector, 8,000 pc from the core (shells are one sector edge, 4 pc, deep)."""


@pytest.fixture(scope="module")
def far_sector(site_db, site_app):
    """The detail of a sector generated 8 kpc from the core, where real sectors lie (the sample sector is at the
    origin): the map's frame holds large numbers and the camera flies a long way to it."""
    target = mysql_argv(site_db)
    state = draw.getstate()
    draw.set_run_seed(20261023)
    try:
        with forced_system_config(PLANETS=True):
            run_cli("galaxy", ["--shell", str(FAR_SHELL), "--band", "0", "--wedge", "0", "--quiet"] + target)
    finally:
        draw.setstate(state)
    client = site_app.test_client()
    for item in client.get("/api/sectors?limit=100").get_json()["items"]:
        detail = client.get(f"/api/sectors/{item['id']}").get_json()
        if detail["shell_index"] == FAR_SHELL:
            assert detail["placed"], detail
            return detail
    raise AssertionError(f"no sector of shell {FAR_SHELL} was generated")


def test_a_star_in_a_sector_far_from_the_core_is_selectable(moving_page, base_url, far_sector, site_app):
    """MAP.210: the same click on a star of a sector 8 kpc out, the distance of most of the galaxy."""
    assert far_sector["center_x_pc"] > 7000, far_sector
    _selects_a_star(moving_page, base_url, site_app, far_sector)
