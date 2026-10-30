"""
The maps in a real browser: every map button does something (TEST.55; the
sector page's map, now the Galaxy Map's engine, is in
`test_web_browser_fixture_maps.py`), the System Map's selection and measuring (TEST.58's browser half: it is
SVG laid out by the browser, which a fake DOM can't stand in for) and
the Galaxy Map's drill-down by clicks (TEST.59).

The pages run against `test_web_a11y.py`'s small generated database in
headless Chromium with SwiftShader WebGL, so the 3D maps really draw: a
button "changes the view" when the map's canvas (or the diagram's
viewBox) is different afterwards. Skipped without Playwright, Chromium
or the MySQL test server, like `test_web_a11y.py`.
"""

import re
import time
from urllib.parse import parse_qs, urlparse

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_a11y import (  # noqa: E402,F401 -- fixtures
    admin_token,
    base_url,
    sample_job,
    sample_run,
    sample_params,
    site_app,
    site_db,
)

VIEWPORT = {"width": 800, "height": 600}
"""dict: The default window. Reading a software-WebGL canvas back costs time in proportion to its pixels (2.5 s
at 1280x900, 1.1 s at 800x600, TEST.139); a test that needs a wide or tall window says so with
`@pytest.mark.viewport(width, height)`."""
GL_ARGS = ["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"]


@pytest.fixture(scope="module")
def gl_browser():
    try:
        playwright = sync_api.sync_playwright().start()
    except Exception as exc:  # noqa: BLE001 -- any failure here means "no browser"
        pytest.skip(f"Playwright could not start ({exc})")
    try:
        chromium = playwright.chromium.launch(args=GL_ARGS)
    except Exception as exc:  # noqa: BLE001
        playwright.stop()
        pytest.skip(f"No Chromium for Playwright ({str(exc).splitlines()[0]})")
    try:
        yield chromium
    finally:
        chromium.close()
        playwright.stop()


@pytest.fixture()
def page(gl_browser, request):
    marker = request.node.get_closest_marker("viewport")
    viewport = {"width": marker.args[0], "height": marker.args[1]} if marker else VIEWPORT
    context = gl_browser.new_context(viewport=viewport, device_scale_factor=1, reduced_motion="reduce",
                                     color_scheme="dark")
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.errors = errors
    try:
        yield page
    finally:
        context.close()
    assert not errors, f"uncaught exceptions: {errors}"


def _open(page, url, ready):
    response = page.goto(url, wait_until="load")
    assert response.status == 200, f"{url} answered {response.status}"
    page.wait_for_selector(ready, state="attached")
    if page.locator("canvas").count() and not page.evaluate("!!document.createElement('canvas').getContext('webgl2')"):
        pytest.skip("no WebGL in this Chromium")
    page.wait_for_timeout(300)


def _shot(page, selector):
    page.wait_for_timeout(150)
    return page.locator(selector).screenshot()


# --- TEST.55, TEST.58: the System Map ------------------------------------------------

def _active_scene(page):
    return page.evaluate("""() => {
        const s = document.querySelector('#sysmap-root .sysmap-svg:not(.sysmap-orbits-layer):not(.sysmap-hidden)');
        return s ? s.dataset.scene : null;
    }""")


PUBLIC_ID = r"[0-9A-F]{14}-[0-9A-F]{7}-[0-9A-F]{3}"
"""str: An object's public ID (API.23): 80 bits as `<birth sector>-<serial>-<body>`, not a row number. TEST.125:
these tests matched `\\d+` and so found no moon on every page, long after the generator was fine."""


def _click_body(marker):
    """Click a System Map body on its dot. A marker's `<g>` also holds its
    label, which the layout can push to one side, so the middle of the
    whole group can be empty map that takes the click (TEST.91)."""
    marker.locator("circle.sysmap-body-fill").click()


def _expand_info(page):
    """Presses Expand on the info card (MAP.347) when it has one that is closed: what is behind it is drawn only then."""
    button = page.locator("#sysmap-info .map-info-expand")
    if button.count() and button.is_visible() and button.get_attribute("aria-expanded") != "true":
        button.click()


@pytest.mark.viewport(1280, 900)  # the planet's info card (MAP.332 sizes it to its content) covers the moons in a window of 800x600
def test_system_map_selection_drill_and_measure(page, base_url, site_app):
    client = site_app.test_client()
    systems = client.get("/api/systems?limit=100").get_json()["items"]
    target = None
    for system in systems:
        resp = client.get(f"/system/{system['id']}")
        # The whole-system scene draws planets only for a single-star system (a system of several stars draws its
        # stars there, and each star's planets in a scene of its own), so the planet with moons must be in that scene.
        whole = re.search(rb'<svg class="sysmap-svg" data-scene="system".*?</svg>', resp.data, re.S)
        if whole and re.search(rb'data-kind="planet"[^>]*data-scene="', whole.group(0)):
            target = system["id"]
            break
    assert target, "no generated single-star system has a planet with moons"
    # The map's own ready mark (its click handlers are on), not a fixed
    # wait: under load the script can still be loading after a second (TEST.104).
    # The diagram, as a visitor who chose it sees it (3D is the default).
    page.add_init_script("try { window.localStorage.setItem('planetgen.systemView', 'diagram'); } catch (e) {}")
    _open(page, f"{base_url}/system/{target}", '#sysmap-root[data-ready="true"]')
    assert _active_scene(page) == "system"

    # A body shows its details.
    body = page.locator('#sysmap-root .sysmap-svg[data-scene="system"] [data-kind="star"]').first
    _click_body(body)
    assert page.locator("#sysmap-info h3").count() == 1
    # MAP.254: every body's panel has the View button to its own page.
    view = page.locator("#sysmap-info .map-info-actions a", has_text="View star")
    view.wait_for(state="attached")
    assert re.fullmatch(rf"/object/star:{PUBLIC_ID}", view.get_attribute("href")), view.get_attribute("href")

    # A planet with moons opens its moons; the crumb's button comes back.
    planet = page.locator('#sysmap-root .sysmap-svg[data-scene="system"] [data-kind="planet"][data-scene]').first
    moons = planet.get_attribute("data-scene")
    _click_body(planet)
    assert _active_scene(page) == moons
    assert "Moons of" in page.locator("#sysmap-crumb").inner_text()
    # MAP.92: the planet and its moons show their radius and mass (behind Expand, MAP.347).
    _expand_info(page)
    info = page.locator("#sysmap-info").inner_text()
    assert "Radius" in info and "Mass" in info, info
    _click_body(page.locator(f'#sysmap-root .sysmap-svg[data-scene="{moons}"] [data-kind="moon"]').first)
    _expand_info(page)
    info = page.locator("#sysmap-info").inner_text()
    assert "Orbits" in info and "Radius" in info and "Mass" in info, info
    # NAV.50: a moon can start or end a course.
    page.wait_for_selector("#sysmap-info .map-info-actions a", state="attached")  # `.all()` does not wait for them
    links = {link.inner_text(): link.get_attribute("href") for link in page.locator("#sysmap-info .map-info-actions a").all()}
    assert re.fullmatch(rf"/nav\?from=moon:{PUBLIC_ID}", links.get("Start Here", "")), links
    assert re.fullmatch(rf"/nav\?to=moon:{PUBLIC_ID}", links.get("End Here", "")), links
    assert re.fullmatch(rf"/object/moon:{PUBLIC_ID}", links.get("View moon →", "")), links
    page.click("#sysmap-crumb .sysmap-back-btn")
    assert _active_scene(page) == "system"
    assert page.locator("#sysmap-crumb button").count() == 0

    # Measure: the button toggles, two bodies give a distance and a path.
    measure = page.locator("#sysmap-measure-btn")
    measure.click()
    assert measure.get_attribute("aria-pressed") == "true"
    assert "Click two" in page.locator("#sysmap-info").inner_text()
    # The star and the planet with moons: a system may have only one
    # planet (TEST.91).
    _click_body(body)
    _click_body(planet)
    assert page.locator(".sysmap-measure-selected").count() == 2
    assert page.locator(".sysmap-measure-path").count() == 1
    assert re.search(r"\d", page.locator("#sysmap-info").inner_text())
    assert _active_scene(page) == "system", "in measure mode a planet with moons doesn't open them"
    measure.click()
    assert measure.get_attribute("aria-pressed") == "false"
    assert page.locator(".sysmap-measure-path").count() == 0


# --- TEST.55: the phenomenon diagram -------------------------------------------------

def _view_box(page):
    return [float(v) for v in page.get_attribute("#phenomenonmap-svg", "viewBox").split()]


def _phenomena(site_app):
    """Every phenomenon in the database (all pages of the list: the
    rogue planets alone run past one)."""
    client = site_app.test_client()
    found = []
    while True:
        items = client.get(f"/api/phenomena?limit=100&offset={len(found)}").get_json()["items"]
        found += items
        if len(items) < 100:
            return found


def _diagram(page, base_url, phenomenon):
    """Opens a phenomenon's page; True when it has the AU diagram (a
    rogue planet or comet has its own body view instead)."""
    _open(page, f"{base_url}/phenomenon/{phenomenon['type']}/{phenomenon['id']}", "main")
    return page.locator("#phenomenonmap-svg").count() == 1


def test_phenomenon_diagram_zoom_in_and_reset(page, base_url, site_app):
    seen = 0
    for phenomenon in [p for p in _phenomena(site_app) if p["type"] in ("nebula", "supernova_remnant")]:
        if not _diagram(page, base_url, phenomenon):
            continue
        seen += 1
        start = _view_box(page)
        scale = page.locator("#phenomenonmap-scale").inner_text()
        page.click('#phenomenonmap-controls [data-action="zoom-in"]')
        zoomed = _view_box(page)
        assert zoomed[2] < start[2], f"{phenomenon['type']}: + didn't zoom in"
        assert page.locator("#phenomenonmap-scale").inner_text() != scale
        page.click('#phenomenonmap-controls [data-action="zoom-out"]')
        assert _view_box(page)[2] > zoomed[2], f"{phenomenon['type']}: - after + didn't zoom out"
        page.click('#phenomenonmap-controls [data-action="zoom-in"]')
        page.click('#phenomenonmap-controls [data-action="reset"]')
        assert _view_box(page) == start, f"{phenomenon['type']}: Reset view"
    assert seen >= 1, "no phenomenon with an AU diagram in the database"


def test_phenomenon_diagram_minus_zooms_out_from_the_start(page, base_url, site_app):
    """UX.38: "-" works from the first view, even on a remnant too big to
    open inside the 1 ly zoom-out limit (a smaller object's view always did)."""
    seen = 0
    for phenomenon in [p for p in _phenomena(site_app) if p["type"] == "supernova_remnant"]:
        if not _diagram(page, base_url, phenomenon):
            continue
        seen += 1
        start = _view_box(page)
        page.click('#phenomenonmap-controls [data-action="zoom-out"]')
        assert _view_box(page)[2] > start[2], "- does nothing on the remnant diagram"
    assert seen >= 1, "no supernova remnant with an AU diagram in the database"


def test_nebula_page_draws_its_shape_in_3d(page, base_url, site_app):
    """MAP.105: a nebula's page swaps its still drawing for the 3D view of
    its shape, with the endpoints answering."""
    nebula = next(p for p in _phenomena(site_app) if p["type"] == "nebula")
    _open(page, f"{base_url}/phenomenon/nebula/{nebula['id']}", "main")
    assert page.locator("#phenomenonmap-svg").count() == 0
    page.locator("#nebulaview-canvas[data-nebula-view='ready']").wait_for(state="attached", timeout=15000)
    # MAP.138: the view can be put back about the nebula's center.
    assert page.locator("#nebulaview-recenter").is_visible()
    page.locator("#nebulaview-recenter").click()
    client = site_app.test_client()
    assert client.get(f"/galaxy/nebula/{nebula['id']}/shape?lod=full").status_code == 200
    assert client.get(f"/galaxy/nebula/{nebula['id']}/surroundings").status_code == 200


# --- TEST.55, TEST.59: the Galaxy Map ------------------------------------------------

def _crumbs(page):
    """Every step to here, galaxy first: the Steps menu's list, which
    holds them all however much of the breadcrumb line is folded into
    "…" (MAP.93, MAP.94)."""
    steps = page.eval_on_selector_all("#galaxymap3d-steps [data-steps-panel] li",
                                      "els => els.map(e => e.textContent.trim())")
    return steps[1:] if steps[:1] == ["Home"] else steps  # UX.59: the trail starts at Home


def _query(page):
    """The address as the tests read it: its query without the camera it carries once the view is still
    (MAP.228), led by "win=m.shell.band.wedge" naming the window on screen (the stage is no part of the address)."""
    parts = [part for part in urlparse(page.url).query.split("&") if part and not part.startswith("cam=")]
    win = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyStage()")["win"]
    if win and not any(part.startswith("sector=") for part in parts):
        centre = win["centre"]
        parts.insert(0, f"win={centre['m']}.{centre['shell']}.{centre['band']}.{centre['wedge']}")
    return "&".join(parts)


def _state(page):
    return _query(page), _crumbs(page)


def _moving(page):
    """Whether the Galaxy Map's camera is still on its way to a stage, waiting for the stage's data or flying: the
    trail follows the camera (MAP.227), a click opens the new stage's address at once and the camera leaves for it
    only when its data has loaded, and the software-drawn map's frames can be more than 250 ms apart, so the address
    and the trail can hold still for that long and then move to the cell the camera lands in (TEST.223)."""
    return page.evaluate("""() => {
        const canvas = document.querySelector('#galaxymap3d-canvas');
        return !!(canvas && canvas.galaxyMoving && canvas.galaxyMoving());
    }""")


def _wait_settled(page, before=None, timeout_s=90.0):
    """Waits until the URL and the steps have held still for 250 ms with the camera not on its way
    (and, given `before`, moved off it). TEST.89: a fixed 250 ms wait
    read the map mid-change under a loaded parallel run. The limit is a ceiling (the wait returns the moment the
    map has settled): the software-drawn map took longer than 15 s to answer a click on a loaded machine, and the
    caller then read the unchanged state and reported "did nothing" (TEST.172, TEST.202)."""
    deadline = time.monotonic() + timeout_s
    page.wait_for_timeout(100)
    last, since = _state(page), time.monotonic()
    while time.monotonic() < deadline:
        page.wait_for_timeout(50)
        now = _state(page)
        if now != last or _moving(page):
            last, since = now, time.monotonic()
        elif (before is None or now != before) and time.monotonic() - since >= 0.25:
            return


def _open_galaxy(page, base_url, query=""):
    _open(page, f"{base_url}/galaxy{query}", "#galaxymap3d-canvas")
    page.wait_for_selector("#galaxymap3d-steps [data-steps-panel] li", state="attached")
    _wait_settled(page)


WINDOW_SIZES = [243, 81, 27, 9, 3, 1]
"""The sizes of the windows on the way down from the whole galaxy to a sector, round the core's cells (the first window
the galaxy opens is the 81-block one)."""


def _descend(page):
    """Opens the window one size down, round the core's cell, as a wheel run toward it hands off to the next (a click
    on a sector block took it down before MAP.297 removed the blocks from the Galaxy Map). Returns the window's
    size."""
    before = _state(page)
    win = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyStage()")["win"]
    size = WINDOW_SIZES[1] if not win else WINDOW_SIZES[WINDOW_SIZES.index(win["centre"]["m"]) + 1]
    cell = {"m": size, "shell": 0, "band": 0, "wedge": 0}
    page.evaluate("(cell) => document.querySelector('#galaxymap3d-canvas').galaxyOpenWindow(cell)", cell)
    _wait_settled(page, before)
    assert _state(page) != before, f"opening the window of {size} did nothing"
    return size


@pytest.mark.viewport(1280, 900)
def test_galaxy_map_drill_down_by_windows(page, base_url):
    _open_galaxy(page, base_url)
    assert _crumbs(page) == ["Galaxy"]
    assert _query(page) == ""
    steps = [("", ["Galaxy"])]
    reached_sector = False
    for _ in range(16):
        _descend(page)
        crumbs = _crumbs(page)
        query = _query(page)
        if crumbs[-1].startswith("Sector "):
            reached_sector = True
            break
        # MAP.150: each click opens a window centred on the cell clicked, one size down.
        assert crumbs[-1].startswith("Block "), crumbs
        assert re.fullmatch(r"win=(81|27|9|3|1)\.\d+\.-?\d+\.\d+", query), query
        assert query != steps[-1][0], "the URL names the new stage"
        steps.append((query, crumbs))
    assert reached_sector, f"never reached a sector: {steps[-1]}"
    kinds = " ".join(q for q, _ in steps)
    for marker in ("win=81.", "win=27.", "win=9.", "win=3."):
        assert marker in kinds, f"no {marker} stage on the way down: {[q for q, _ in steps]}"
    assert not re.search(r"[=,]r\d|p=", kinds), kinds

    # Back and Forward (the browser's) walk the same stages, in order.
    # (Picking the sector may or may not add an entry of its own, so the
    # walk is matched against the stages seen rather than counted.)
    steps.append((_query(page), _crumbs(page)))
    queries = [q for q, _ in steps]
    at = len(steps) - 1
    for _ in range(2):  # (three would reach the galaxy, where the map's own Back is off)
        before = _state(page)
        page.go_back()
        _wait_settled(page, before)
        query = _query(page)
        assert query in queries[:at], f"Back went to {query!r}, not an earlier stage of {queries[:at]}"
        at = max(n for n in range(at) if queries[n] == query)
        # The trail follows the camera (MAP.227): give a flight that a loaded run stalls time to land.
        for _ in range(100):
            if _crumbs(page) == steps[at][1]:
                break
            page.wait_for_timeout(100)
        assert _crumbs(page) == steps[at][1]
    before = _state(page)
    page.go_forward()
    _wait_settled(page, before)
    forward = _query(page)
    assert forward in queries[at + 1:], f"Forward went to {forward!r}"
    before = _state(page)
    page.go_back()
    _wait_settled(page, before)
    assert _query(page) == queries[at]

    # The map's own Back, Forward and Up buttons.
    here = _query(page)
    before = _state(page)
    page.click('#galaxymap3d-controls [data-action="back"]')
    _wait_settled(page, before)
    back = _query(page)
    assert back != here and back in queries
    before = _state(page)
    page.click('#galaxymap3d-controls [data-action="forward"]')
    _wait_settled(page, before)
    assert _query(page) == here
    crumbs = _crumbs(page)
    before = _state(page)
    page.click('#galaxymap3d-controls [data-action="up"]')
    _wait_settled(page, before)
    # The trail follows the camera (MAP.227): the window above is centred where the camera lands, so its steps need
    # not be the ones the window below was inside, but there is one fewer.
    assert len(_crumbs(page)) == len(crumbs) - 1 and _crumbs(page)[0] == crumbs[0], "Up drops the last breadcrumb"

    # An address opens the view it was copied at directly (MAP.228): the camera is in it, the stage worked out.
    page.wait_for_timeout(600)
    deep_url, deep_crumbs = page.url, _crumbs(page)
    assert "cam=" in deep_url, deep_url
    _open(page, deep_url, "#galaxymap3d-canvas")
    page.wait_for_selector("#galaxymap3d-steps [data-steps-panel] li", state="attached")
    _wait_settled(page)
    assert _crumbs(page) == deep_crumbs

    # A breadcrumb button goes back up; Reset goes home.
    more = page.locator("#galaxymap3d-crumbs .crumb-menu")
    if more.count():  # the second step is folded into "…" (MAP.93)
        more.locator("summary").click()
        more.locator("button").first.click()
    else:
        page.locator("#galaxymap3d-crumbs button.crumb").nth(1).click()
    _wait_settled(page)
    # The crumb is the 243-block above where the walk ended, which the clicks may have carried past the first
    # block's edge: the 3-ary ladder takes more steps than the old one did, each centred on the cell clicked.
    crumbs = _crumbs(page)
    assert len(crumbs) == 2 and crumbs[0] == "Galaxy" and crumbs[1].endswith("(243)"), crumbs
    before = _state(page)
    page.click('#galaxymap3d-controls [data-action="reset"]')
    _wait_settled(page, before)
    assert _crumbs(page) == ["Galaxy"]
    assert _query(page) == ""


def _scale_text(page):
    """The scale line: its label and its bar's width (a zoom between two
    nice lengths only changes the bar)."""
    return page.locator("#galaxymap3d-scale").inner_html()


def test_galaxy_map_free_camera_from_a_window_down(page, base_url):
    _open_galaxy(page, base_url)
    reset_view = page.locator('#galaxymap3d-controls [data-action="reset-view"]')
    assert not reset_view.is_disabled(), "the whole galaxy turns (MAP.85)"
    _descend(page)
    assert parse_qs(_query(page))["win"][0].startswith("81."), "the galaxy's first pick is a window of 81-blocks"

    canvas = "#galaxymap3d-canvas"
    box = page.locator(canvas).bounding_box()
    middle = (box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    first = _shot(page, canvas)
    scale = _scale_text(page)
    stage = (_query(page), _crumbs(page))

    page.mouse.move(*middle)
    page.mouse.down()
    page.mouse.move(middle[0] + 120, middle[1] + 40, steps=6)
    page.mouse.up()
    turned = _shot(page, canvas)
    assert turned != first, "dragging didn't turn the view"
    assert (_query(page), _crumbs(page)) == stage, "a drag picks nothing"

    page.mouse.move(*middle)
    page.mouse.wheel(0, 400)
    page.wait_for_timeout(150)
    assert _scale_text(page) != scale, "the wheel didn't zoom"

    page.mouse.move(*middle)
    page.mouse.down(button="right")
    page.mouse.move(middle[0] - 150, middle[1], steps=6)
    page.mouse.up(button="right")
    assert _shot(page, canvas) != turned, "right-drag didn't move the view"

    page.click("#galaxymap3d-menu > summary")
    reset_view.click()
    _wait_settled(page)
    assert _scale_text(page) == scale, "Reset view is back at the stage's own zoom"
    assert (_query(page), _crumbs(page)) == stage


@pytest.mark.viewport(1280, 900)
def test_galaxy_map_buttons(page, base_url):
    _open_galaxy(page, base_url)
    controls = page.locator("#galaxymap3d-controls")
    actions = controls.locator("[data-action]").evaluate_all("els => els.map(e => e.dataset.action)")
    canvas = "#galaxymap3d-canvas"
    for action in ("back", "forward", "current", "up"):
        assert page.locator(f'#galaxymap3d-controls [data-action="{action}"]').is_disabled(), f"{action} at the galaxy"
    assert "wedges" not in actions, "no Wedges button (MAP.85)"
    page.click("#galaxymap3d-menu > summary")

    if "charted-only" in actions:
        only = controls.locator('[data-action="charted-only"]')
        before = _shot(page, canvas)
        only.click()
        assert only.get_attribute("aria-pressed") == "true"
        assert _shot(page, canvas) != before, "Charted only didn't change the map"
        only.click()
        assert only.get_attribute("aria-pressed") == "false"

    if "territories" in actions:
        territories = controls.locator('[data-action="territories"]')
        territories.click()
        page.wait_for_selector("#galaxymap3d-territories h3")
        assert page.locator("#galaxymap3d-territories").is_visible()
        assert territories.get_attribute("aria-pressed") == "true"
        territories.click()
        assert not page.locator("#galaxymap3d-territories").is_visible()

    # One step down: Back and Up work; Back then enables Forward. (The Menu's
    # star filters make it tall enough to cover the map, so it is closed.)
    page.click("#galaxymap3d-menu > summary")
    _descend(page)
    for action in ("back", "up"):
        assert not controls.locator(f'[data-action="{action}"]').is_disabled()
    controls.locator('[data-action="back"]').click()
    _wait_settled(page)
    assert _crumbs(page) == ["Galaxy"]
    assert not controls.locator('[data-action="forward"]').is_disabled()


def test_galaxy_map_forward_to_current_jumps_to_the_newest_view(page, base_url):
    """MAP.95: three levels in, back two, then one press of Current is
    back at the deepest view; it is disabled once there."""
    _open_galaxy(page, base_url)
    controls = page.locator("#galaxymap3d-controls")
    steps = page.locator("#galaxymap3d-steps")
    current = steps.locator('[data-action="current"]')
    for _ in range(3):
        _descend(page)
    deepest = (_query(page), _crumbs(page))
    assert current.is_disabled(), "already at the newest view"
    for _ in range(2):
        controls.locator('[data-action="back"]').click()
        _wait_settled(page)
    assert (_query(page), _crumbs(page)) != deepest
    steps.locator("summary").click()
    assert not current.is_disabled()
    before = _state(page)
    current.click()
    _wait_settled(page, before)
    assert (_query(page), _crumbs(page)) == deepest
    assert current.is_disabled()


# --- NAV.40: picking a course end from a bookmark ---------------------------------

def _seed_bookmarks(page, entries):
    page.evaluate("""(entries) => {
        const db = document.querySelector("[data-bookmark-db]").getAttribute("data-bookmark-db");
        localStorage.setItem("planetgen.bookmarks." + db, JSON.stringify(entries));
    }""", entries)


def _course_ends(page):
    params = parse_qs(urlparse(page.url).query)
    return params.get("from", [None])[0], params.get("to", [None])[0]


def test_nav_ends_picked_from_bookmarks_on_every_page(page, base_url, sample_params):
    """NAV.40: while picking a start or destination, a bookmark on the NAV
    page, the Galaxy Map or the sector page sets that end and keeps the
    other; on the course page a bookmark replaces either end."""
    start, dest = f"system:{sample_params['from_id']}", f"system:{sample_params['to_id']}"
    assert start != dest
    systems = page.request.get(f"{base_url}/api/systems?sector_id={sample_params['sector_id']}&limit=100").json()["items"]
    other = next(f"system:{s['id']}" for s in systems if f"system:{s['id']}" not in (start, dest))
    entries = [
        {"name": "Bookmarked destination", "kind": "system", "value": dest, "url": f"/system/{sample_params['to_id']}", "created": 1},
        {"name": "Bookmarked other", "kind": "system", "value": other, "url": "/system/" + other.split(":")[1], "created": 2},
    ]
    _open(page, f"{base_url}/nav?from={start}", "body")
    _seed_bookmarks(page, entries)

    # The NAV page's destination step.
    _open(page, f"{base_url}/nav?from={start}", "form[data-bookmarks-nav]")
    form = page.locator("form[data-bookmarks-nav]")
    form.locator("select").select_option(dest)
    with page.expect_navigation():
        form.locator("button[type=submit]").click()
    assert _course_ends(page) == (start, dest)
    assert page.locator("#course-heading").count() == 1

    # The course page: a bookmark as the new start, then the new
    # destination, each keeping the other end.
    group = page.locator("[data-bookmarks-nav-group]")
    assert group.is_visible()
    new_start = group.locator("form[data-pick=from]")
    new_start.locator("select").select_option(other)
    with page.expect_navigation():
        new_start.locator("button[type=submit]").click()
    assert _course_ends(page) == (other, dest)
    new_dest = page.locator("[data-bookmarks-nav-group] form[data-pick=to]")
    assert new_dest.locator(f"option[value='{other}']").count() == 0, "not the start already chosen"
    new_dest.locator("select").select_option(dest)
    with page.expect_navigation():
        new_dest.locator("button[type=submit]").click()
    assert _course_ends(page) == (other, dest)

    # The Galaxy Map's menu while picking a destination.
    _open(page, f"{base_url}/galaxy?pick=to&from={start}", "#galaxymap3d-canvas")
    page.locator("[data-bookmarks-menu] summary").click()
    with page.expect_navigation():
        page.locator("[data-bookmarks-panel] a", has_text="Bookmarked destination").click()
    assert _course_ends(page) == (start, dest)

    # The sector page's menu while picking a destination.
    _open(page, f"{base_url}/sector/{sample_params['sector_id']}?pick=to&from={start}", "#galaxymap3d-canvas")
    page.locator(".pick-bookmarks summary").click()
    with page.expect_navigation():
        page.locator(".pick-bookmarks a", has_text="Bookmarked destination").click()
    assert _course_ends(page) == (start, dest)


# --- MAP.72 to MAP.74: the 3D system view ------------------------------------------------

def _system_with_planets(client):
    for system in client.get("/api/systems?limit=100").get_json()["items"]:
        scene = client.get(f"/api/systems/{system['id']}/scene").get_json()
        if scene["planets"]:
            return system["id"], scene
    pytest.fail("no generated system has a planet")


def _canvas_pixels(page):
    return page.locator("#sysview3d-canvas").screenshot()


def test_system_page_names_its_bodies_through_the_label_engine_with_one_switch(page, base_url, site_app):
    """MAP.231: the 3D view's names are the label engine's (the same `map-label` elements as the Galaxy Map's, none
    over another), the diagram's names are placed by the engine too, and the Menu's Body names switch hides both."""
    system_id, scene = _system_with_planets(site_app.test_client())
    page.goto(f"{base_url}/system/{system_id}?view=3d", wait_until="load")
    page.wait_for_selector('#sysview3d[data-ready="true"]')
    page.wait_for_timeout(500)
    labels = page.evaluate("() => [...document.querySelectorAll('#sysview3d-labels .map-label')].map((e) => [e.textContent, e.dataset.kind])")
    names = {text for text, _kind in labels}
    known = {body["name"] for body in scene["stars"] + scene["planets"] + scene["comets"]} \
        | {moon["name"] for planet in scene["planets"] for moon in planet["moons"]}
    assert names and names <= known, (names - known, names)
    assert not page.locator("#sysview3d-labels .sysview-label").count(), "no second label code"
    rects = page.evaluate("() => [...document.querySelectorAll('#sysview3d-labels .map-label')].map((e) => { const r = e.getBoundingClientRect(); return [r.left, r.top, r.right, r.bottom]; })")
    for i, a in enumerate(rects):
        for b in rects[i + 1:]:
            assert not (a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]), (a, b)
    page.click("#sysmap-view-diagram")
    assert page.locator("text.sysmap-label").count() > 0
    page.click("#sysmap-menu > summary")
    toggle = page.locator('#sysmap-layers [data-label="body"]')
    toggle.click()
    assert page.locator("#sysmap-root.sysmap-hide-labels").count() == 1
    assert "label-body" in page.url
    assert not page.locator("text.sysmap-label").first.is_visible()


@pytest.mark.fresh_database
@pytest.mark.parametrize("wide", [False, True])
def test_system_map_diagram_draws_each_star_beyond_the_pair_with_its_planets(page, mysql_config, wide):
    """MAP.290: a system of three on the 2D System Map shows its third star as a marker that opens a scene of its own
    with the planets that orbit it, and Back returns to the whole system."""
    from planetgen.api.config import Config
    from planetgen.db import store
    from planetgen.galaxy.sector import SpaceSector
    from planetgen.web.app import create_app
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    from tests.publicids import pid
    from tests.served import served

    for seed in range(400):    # a system of three stars of the asked pair with planets round a further star
        cfg = SystemConfig()
        cfg.STARS = 3
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(body.body_type != "a" for bodies in system.member_planets for body in bodies) and (system.binary_type == "wide") == wide:
            break
    else:
        pytest.fail("no three-star system with planets round a further star in 400 draws")
    sector = SpaceSector("Member Map Sector", edge_ly=10.0)
    sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
    sector_id = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
        member = conn.execute("SELECT id, name FROM stars WHERE star_system_id = ? AND role = 'member'",
                              (system_id,)).fetchone()
        member_planets = [r["name"] for r in conn.execute(
            "SELECT name FROM planets WHERE star_id = ?", (member["id"],)).fetchall()]
    finally:
        conn.close()
    assert member_planets

    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    app = create_app(RealConfig)
    app.testing = True
    with served(app) as base_url:
        _open(page, f"{base_url}/system/{pid('system', system_id)}", '#sysmap-root[data-ready="true"]')
        page.evaluate("window.localStorage.setItem('planetgen.systemView', 'diagram')")
        page.goto(f"{base_url}/system/{pid('system', system_id)}", wait_until="load")
        page.wait_for_selector('#sysmap-root[data-ready="true"]', state="attached")
        if page.locator("#sysmap-diagram").is_hidden():
            page.click("#sysmap-view-diagram")
        assert _active_scene(page) == "system"
        marker = page.locator(f'#sysmap-root .sysmap-svg:not(.sysmap-hidden) [data-kind="star"][data-name="{member["name"]}"]')
        assert marker.count() == 1 and marker.get_attribute("data-scene").startswith("star-")
        marker.dispatch_event("click")
        assert _active_scene(page).startswith("star-")
        shown = page.evaluate("""() => [...document.querySelectorAll('#sysmap-root .sysmap-svg:not(.sysmap-orbits-layer):not(.sysmap-hidden) .sysmap-planet')]
            .map((e) => e.dataset.name)""")
        assert sorted(shown) == sorted(member_planets), (shown, member_planets)
        assert page.locator(".sysmap-crumb-current").inner_text() == f"Planets of {member['name']}"
        page.click(".sysmap-back-btn")
        assert _active_scene(page) == "system"


def test_system_page_3d_view_draws_switches_scale_and_keeps_the_diagram(page, base_url, site_app):
    system_id, scene = _system_with_planets(site_app.test_client())
    _open(page, f"{base_url}/system/{system_id}", '#sysmap-root[data-ready="true"]')
    # 3D is the default with no choice remembered.
    page.evaluate("window.localStorage.clear()")
    page.goto(f"{base_url}/system/{system_id}", wait_until="load")
    page.wait_for_selector('#sysview3d[data-ready="true"]')
    assert page.locator("#sysmap-diagram").is_hidden() and page.locator("#sysview3d").is_visible()
    assert "view=3d" in page.url
    page.wait_for_timeout(500)
    compressed = _canvas_pixels(page)
    assert len(set(compressed[100:3000])) > 20, "the canvas drew something"
    assert "Compressed" in page.locator("#sysview3d-scale-note").inner_text()

    # Every body is in the list for screen readers, and picking one shows it.
    expected = 1 + len(scene["planets"]) + sum(len(p["moons"]) for p in scene["planets"]) + len(scene["comets"]) \
        + (len(scene["stars"]) - 1)
    page.click(".sysview-list summary")
    assert page.locator("#sysview3d-list button").count() == expected
    page.locator("#sysview3d-list button").nth(1 if expected > 1 else 0).click()
    assert page.locator("#sysmap-info h3").count() == 1
    assert "object=" in page.url
    # NAV.50: the picked body can start or end a course.
    page.wait_for_selector("#sysmap-info .map-info-actions a", state="attached")  # `.all()` does not wait for them
    hrefs = [link.get_attribute("href") for link in page.locator("#sysmap-info .map-info-actions a").all()]
    assert any(re.fullmatch(rf"/nav\?to=[a-z]+:{PUBLIC_ID}", href or "") for href in hrefs), hrefs

    # MAP.138: Center moves the view's center to the picked body.
    before_center = _canvas_pixels(page)
    page.click("#sysview3d-recenter")
    page.wait_for_timeout(300)
    assert _canvas_pixels(page) != before_center, "Center changes the view"

    # The scale changes the picture and its note.
    page.select_option("#sysview3d-scale", "true")
    page.wait_for_timeout(500)
    assert "True scale" in page.locator("#sysview3d-scale-note").inner_text()
    assert _canvas_pixels(page) != compressed

    # Time: pause freezes the label, play resumes.
    page.click("#sysview3d-play")
    assert page.locator("#sysview3d-rate").inner_text() == "Paused"
    page.click("#sysview3d-faster")
    assert page.locator("#sysview3d-rate").inner_text() != "Paused"
    page.click("#sysview3d-now")
    assert page.locator("#sysview3d-rate").inner_text() == "Time: real time (1×)"

    # Reload on the address: back in 3D with the body selected.
    selected = page.locator("#sysmap-info h3").inner_text()
    page.reload(wait_until="load")
    page.wait_for_selector('#sysview3d[data-ready="true"]')
    page.wait_for_timeout(300)
    assert page.locator("#sysmap-info h3").inner_text() == selected

    # The diagram is still there.
    page.click("#sysmap-view-diagram")
    assert page.locator("#sysmap-diagram").is_visible() and page.locator("#sysview3d").is_hidden()
    assert "view=3d" not in page.url
    # The diagram, once chosen, is what opens next.
    page.goto(f"{base_url}/system/{system_id}", wait_until="load")
    page.wait_for_selector('#sysmap-root[data-ready="true"]')
    assert page.locator("#sysview3d").is_hidden() and page.locator("#sysmap-diagram").is_visible()


def test_system_page_3d_view_paints_its_bodies_from_the_look_record(page, base_url, site_app):
    """MAP.287: the 3D system view takes the body look record the System Map takes: a star, planet and moon is drawn
    from its `look` (seeded, composition-led), and once the painters have run a star and a solid world show a painted
    texture of their own seed, where the view used to tint every body by its class letter."""
    system_id, scene = _system_with_planets(site_app.test_client())
    assert all("color" not in body for body in scene["stars"] + scene["planets"])
    star = scene["stars"][0]
    _open(page, f"{base_url}/system/{system_id}", '#sysmap-root[data-ready="true"]')
    page.evaluate("window.localStorage.clear()")
    page.goto(f"{base_url}/system/{system_id}", wait_until="load")
    page.wait_for_selector('#sysview3d[data-ready="true"]')
    # A body is painted once it is big enough on screen to show it: fly to the star, then to a planet.
    for body in [star, next(p for p in scene["planets"] if p["look"]["body"] == "solid")]:
        page.evaluate("(ref) => document.querySelector('#sysview3d-canvas').systemFlyTo(ref)", body["ref"])
        page.wait_for_function(
            """(ref) => document.querySelector('#sysview3d-canvas').systemBodies()
                .some((b) => b.ref === ref && b.texture)""", arg=body["ref"], timeout=20000)
        shown = next(b for b in page.evaluate("document.querySelector('#sysview3d-canvas').systemBodies()")
                     if b["ref"] == body["ref"])
        assert shown["texture"].split(":")[0] == str(body["look"]["seed"]), (body["ref"], shown)


def test_system_page_3d_view_draws_a_triple_with_every_star_clear_of_the_others(page, base_url, site_app):
    """MAP.213: a system of three stars is drawn node by node: the scene's tree puts the inner pair about its centre and the
    third star about the whole, each on its own scale, so the three stars sit apart (none at the middle, on top of another)."""
    import json
    import math

    from tests.test_body_positions import triple_scene

    system_id, scene = _system_with_planets(site_app.test_client())
    tree = triple_scene()["tree"]
    refs = {"star:1": "star:9001", "star:2": "star:9002", "star:3": "star:9003"}
    scene["tree"] = {"nodes": tree["nodes"], "origin": None,
                     "members": [dict(member, ref=refs[member["ref"]]) for member in tree["members"]]}
    stars = []
    for number, ref in enumerate(refs.values(), start=1):
        star = dict(scene["stars"][0], ref=ref, id=9000 + number, name=f"Triple {number}", role="member" if number == 3 else "primary",
                    orbit=None, position_km=[0.0, 0.0, 0.0])
        stars.append(star)
    scene["stars"] = stars
    scene["epoch_unix"] = tree["nodes"][0]["epoch_unix"]
    for planet in scene["planets"]:
        planet["orbit"]["around"] = "barycenter"
    body = json.dumps(scene)
    page.route(re.compile(r"/system/[^/?]+/scene"),
               lambda route: route.fulfill(status=200, content_type="application/json", body=body))
    _open(page, f"{base_url}/system/{system_id}?view=3d", '#sysview3d[data-ready="true"]')
    page.wait_for_timeout(500)
    placed = {b["ref"]: b for b in page.evaluate("document.querySelector('#sysview3d-canvas').systemView.placed()")}
    assert set(refs.values()) <= set(placed)
    for a, b in (("star:9001", "star:9002"), ("star:9001", "star:9003"), ("star:9002", "star:9003")):
        gap = math.dist([placed[a][axis] for axis in "xyz"], [placed[b][axis] for axis in "xyz"])
        assert gap > placed[a]["radius"] + placed[b]["radius"], (a, b, gap)


def _system_with_moons(client):
    for system in client.get("/api/systems?limit=100").get_json()["items"]:
        scene = client.get(f"/api/systems/{system['id']}/scene").get_json()
        if any(planet["moons"] for planet in scene["planets"]):
            return system["id"], scene
    pytest.fail("no generated system has a planet with moons")


def _fit_half_field(page):
    """Half the smaller field of view of the 3D system view's camera now, in radians."""
    return page.evaluate("""() => {
        const view = document.querySelector('#sysview3d-canvas').systemView;
        const canvas = document.querySelector('#sysview3d-canvas');
        const tan = Math.tan(view.camera.fov * Math.PI / 360);
        return Math.atan(Math.min(tan, tan * canvas.clientWidth / canvas.clientHeight));
    }""")


def test_flying_to_a_planet_with_moons_fills_the_view_with_its_moon_system(page, base_url, site_app):
    """MAP.319: flying to a planet that has moons puts the camera where the sphere holding the planet and its
    outermost moon orbit fills the smaller side of the view (8 percent to spare), and a planet with no moon is just
    come close to."""
    system_id, scene = _system_with_moons(site_app.test_client())
    page.goto(f"{base_url}/system/{system_id}?view=3d", wait_until="load")
    page.wait_for_selector('#sysview3d[data-ready="true"]')
    planet = next(p for p in scene["planets"] if p["moons"])
    page.evaluate("(ref) => document.querySelector('#sysview3d-canvas').systemFlyTo(ref)", planet["ref"])
    page.wait_for_function("() => !document.querySelector('#sysview3d-canvas').systemView.flying()", timeout=20000)
    got = page.evaluate("""(ref) => { const v = document.querySelector('#sysview3d-canvas').systemView;
        return {distance: v.distance(), reach: v.reachOf(ref), framing: v.framing(v.reachOf(ref))}; }""", planet["ref"])
    assert got["reach"] and got["reach"] > 0, got
    assert got["distance"] == pytest.approx(got["framing"], rel=1e-6), got
    # The sphere of that radius, seen from that far, takes 1/1.08 of half the smaller field.
    import math
    seen = math.asin(got["reach"] / got["distance"])
    assert seen == pytest.approx(_fit_half_field(page) / 1.08, rel=0.02), (seen, got)
    bare = [p for p in scene["planets"] if not p["moons"]]
    if bare:
        page.evaluate("(ref) => document.querySelector('#sysview3d-canvas').systemFlyTo(ref)", bare[0]["ref"])
        page.wait_for_function("() => !document.querySelector('#sysview3d-canvas').systemView.flying()", timeout=20000)
        assert page.evaluate("(ref) => document.querySelector('#sysview3d-canvas').systemView.reachOf(ref)", bare[0]["ref"]) is None


def test_the_system_view_opens_with_every_orbit_in_the_window(page, base_url, site_app):
    """MAP.319: a star system's view starts at the distance where the sphere holding every body and orbit line out to
    the outermost apoapsis fills the view."""
    system_id, scene = _system_with_planets(site_app.test_client())
    page.goto(f"{base_url}/system/{system_id}?view=3d", wait_until="load")
    page.wait_for_selector('#sysview3d[data-ready="true"]')
    got = page.evaluate("""() => { const v = document.querySelector('#sysview3d-canvas').systemView;
        return {distance: v.distance(), reach: v.reachOf(null)}; }""")
    import math
    assert got["reach"] > 0, got
    seen = math.asin(min(1.0, got["reach"] / got["distance"]))
    assert seen <= _fit_half_field(page) / 1.08 + 1e-6, (seen, got)


def _moon_scene_facts(page, scene_id):
    return page.evaluate("""(id) => {
        const scene = document.querySelector('#sysmap-root .sysmap-svg:not(.sysmap-orbits-layer)[data-scene="' + id + '"]');
        const orbits = document.querySelector('#sysmap-root .sysmap-orbits-layer[data-scene="' + id + '"]');
        const circle = (g) => g.querySelector('circle.sysmap-body-fill');
        const self = scene.querySelector('[data-self="true"]');
        const rect = (el) => el.getBoundingClientRect();
        const planet = rect(circle(self));
        return {
            hidden: scene.classList.contains('sysmap-hidden'),
            paths: orbits.querySelectorAll('path.sysmap-orbit').length,
            circles: orbits.querySelectorAll('circle.sysmap-orbit').length,
            moons: scene.querySelectorAll('[data-kind="moon"]').length,
            perKm: parseFloat(scene.dataset.pxperkm), enlarge: parseFloat(scene.dataset.enlarge),
            note: scene.dataset.scaleNote, twin: scene.dataset.scaleTwin,
            planetPx: planet.width / 2, look: !!self.dataset.look,
            scale: scene.getBoundingClientRect().width / scene.viewBox.baseVal.width,
        };
    }""", scene_id)


def test_a_planet_with_moons_has_its_own_moon_system_map(page, base_url, site_app):
    """MAP.308: a planet with moons opens a map of them at one linear scale: each moon on its stored orbit's ellipse,
    the planet drawn from its look record, a note of how much bodies are enlarged and a True scale switch that draws
    them at their true size; the planet's page links straight to it, and the link is in the URL."""
    client = site_app.test_client()
    system_id, scene = _system_with_moons(client)
    planet = next(p for p in scene["planets"] if p["moons"])
    # A planet whose bodies are drawn bigger than true has a True scale twin; find one if the generated systems have it.
    twinned = None
    for system in client.get("/api/systems?limit=100").get_json()["items"]:
        found = re.search(rb'data-scene="(planet-[^"]+)"(?=[^>]*data-scale-twin)', client.get(f"/system/{system['id']}").data)
        if found:
            twinned = (system["id"], found.group(1).decode())
            break
    if twinned:
        system_id, twin_scene = twinned
        scene = client.get(f"/api/systems/{system_id}/scene").get_json()
        planet = next(p for p in scene["planets"] if p["moons"] and f"planet-{p['id']}" == twin_scene)
    page.add_init_script("try { window.localStorage.setItem('planetgen.systemView', '3d'); } catch (e) {}")
    link = page.request.get(f"{base_url}/object/{planet['ref']}")
    assert link.status == 200
    moon_href = re.search(r'href="([^"]*scene=planet-[^"]+)">moon system map', link.text())
    assert moon_href, "the planet's page offers its moon system map"
    _open(page, base_url + moon_href.group(1).replace("&amp;", "&"), '#sysmap-root[data-ready="true"]')
    default_id = f"planet-{planet['id']}"
    assert _active_scene(page) == default_id
    facts = _moon_scene_facts(page, default_id)
    assert not facts["hidden"] and facts["look"], facts
    assert facts["paths"] == len(planet["moons"]) and facts["circles"] == 0 and facts["moons"] == len(planet["moons"]), facts
    assert page.locator("#sysview3d").is_hidden(), "a link to a moon system map opens the diagram, not the 3D view"
    crumb = page.locator("#sysmap-crumb")
    assert "Moons of " + planet["name"] in crumb.inner_text()
    assert facts["note"]
    if not twinned:
        assert facts["twin"] is None and facts["enlarge"] == pytest.approx(1.0, abs=0.06), facts
        assert crumb.get_by_role("button", name="True scale").count() == 0
        return
    assert facts["twin"] == default_id + "-true"
    # The True scale switch: the other scene, the URL follows, the planet is drawn to the scale of the orbits.
    crumb.get_by_role("button", name="True scale").click()
    assert _active_scene(page) == default_id + "-true"
    assert f"scene={default_id}-true" in page.url
    true = _moon_scene_facts(page, default_id + "-true")
    assert true["perKm"] == pytest.approx(facts["perKm"], rel=1e-4), (facts, true)
    assert true["enlarge"] == pytest.approx(1.0) and "true size" in true["note"], true
    radius_km = planet["radius_km"]
    assert true["planetPx"] / true["scale"] == pytest.approx(max(radius_km * true["perKm"], 2.0), abs=1.0), (true, radius_km)
    assert facts["planetPx"] / facts["scale"] >= true["planetPx"] / true["scale"] - 1.0
    crumb.get_by_role("button", name="Enlarge bodies").click()
    assert _active_scene(page) == default_id
    crumb.get_by_role("button", name="Back to system").click()
    assert _active_scene(page) == "system" and "scene=" not in page.url


def test_a_planet_without_moons_offers_no_moon_system_map(page, base_url, site_app):
    """MAP.308: a planet without moons has no moon system map, and its page does not offer one."""
    client = site_app.test_client()
    for system in client.get("/api/systems?limit=100").get_json()["items"]:
        scene = client.get(f"/api/systems/{system['id']}/scene").get_json()
        bare = [p for p in scene["planets"] if not p["moons"]]
        if bare:
            break
    else:
        pytest.skip("every generated planet has a moon")
    response = page.request.get(f"{base_url}/object/{bare[0]['ref']}")
    assert response.status == 200 and "moon system map" not in response.text() and "scene=planet-" not in response.text()
    page.goto(f"{base_url}/system/{system['id']}?scene=planet-{bare[0]['id']}", wait_until="load")
    page.wait_for_selector('#sysmap-root[data-ready="true"]', state="attached")


def _inside_window(page, selector):
    box = page.evaluate("""(sel) => {
        const r = document.querySelector(sel).getBoundingClientRect();
        return {top: r.top, bottom: r.bottom, left: r.left, right: r.right,
                w: window.innerWidth, h: window.innerHeight};
    }""", selector)
    return box, (box["top"] >= -1 and box["left"] >= -1 and box["bottom"] <= box["h"] + 1
                 and box["right"] <= box["w"] + 1)


@pytest.mark.parametrize("size", [{"width": 1280, "height": 420}, {"width": 400, "height": 520}])
@pytest.mark.parametrize("menu, panel", [("#galaxymap3d-menu", "#galaxymap3d-menu .galaxy-menu-panel"),
                                         ("#galaxymap3d-steps", "#galaxymap3d-steps .galaxy-steps-panel")])
def test_galaxy_map_menus_open_inside_the_window(page, base_url, menu, panel, size):
    """UX.85: a button menu opens where it can be seen, in a short or a narrow window, with the
    controls scrolled to the bottom edge, where a panel hung below its button ran off the screen."""
    page.set_viewport_size(size)
    _open_galaxy(page, base_url)
    page.evaluate("document.querySelector('#galaxymap3d-controls').scrollIntoView({block: 'end'})")
    page.click(f"{menu} > summary")
    page.wait_for_selector(f"{menu}[data-placed]")        # menuplace.js has put the panel where it goes
    box, inside = _inside_window(page, panel)
    assert inside, f"{menu} opened out of sight in a {size['width']}x{size['height']} window: {box}"


@pytest.mark.parametrize("size", [{"width": 1280, "height": 800}, {"width": 400, "height": 800}])
def test_galaxy_map_menu_filters_sit_side_by_side(page, base_url, size):
    """UX.86: the Menu's filters take their natural width and wrap in a row, not one full-width
    button each, and the Menu's buttons are the smaller size."""
    page.set_viewport_size(size)
    _open_galaxy(page, base_url)
    page.click("#galaxymap3d-menu > summary")
    page.wait_for_selector("#galaxymap3d-menu[data-placed]")
    chips = page.evaluate("""() => [...document.querySelectorAll('#galaxymap3d-kinds .starmap-toggle')]
        .map(b => { const r = b.getBoundingClientRect(); return {w: r.width, top: Math.round(r.top)}; })""")
    panel = page.evaluate("document.querySelector('#galaxymap3d-menu .galaxy-menu-panel').getBoundingClientRect().width")
    assert len(chips) >= 7
    # Natural width, not one full-width button each: a long label (a 70-chip Menu has a few) may pass half the panel,
    # none takes the row, and nearly all stay under half (TEST.418).
    assert all(c["w"] < panel * 0.75 for c in chips), chips
    assert sum(c["w"] < panel / 2 for c in chips) >= 0.9 * len(chips), chips
    assert len({c["top"] for c in chips}) < len(chips) / 2, "toggles should share rows"
    assert page.locator("#galaxymap3d-show > summary").is_visible()




# --- MAP.328: the bodies turn slowly on their axes ---------------------------------

SPIN_SAMPLE = """async (wait) => {
    const read = () => document.querySelector('#sysmap-spheres-canvas').__spins();
    const first = read();
    await new Promise((resolve) => setTimeout(resolve, wait));
    return {first: first, second: read()};
}"""


def _single_star_system_with_planets(client):
    """A generated system of one star whose planets store a rotation period: the whole-system scene then draws planets
    (a system of several stars draws only its stars there)."""
    for system in client.get("/api/systems?limit=100").get_json()["items"]:
        scene = client.get(f"/api/systems/{system['id']}/scene").get_json()
        if len(scene["stars"]) == 1 and sum(1 for p in scene["planets"] if p["look"].get("spinH")) >= 2:
            return system["id"], scene
    pytest.fail("no generated system has one star and two planets with a stored rotation period")


def _open_spinning_diagram(page, base_url, site_app, reduced="no-preference"):
    system_id, scene = _single_star_system_with_planets(site_app.test_client())
    page.emulate_media(reduced_motion=reduced)
    _open(page, f"{base_url}/system/{system_id}", "#sysmap-root")
    _show_diagram(page)
    return system_id, scene


def _show_diagram(page):
    if page.locator("#sysmap-view-diagram").count():
        page.click("#sysmap-view-diagram")
    page.wait_for_function("() => { const c = document.querySelector('#sysmap-spheres-canvas');"
                           " return c && c.__spins && c.__spins().length > 0; }")


def test_the_system_map_spheres_turn_at_the_rate_their_stored_period_gives_on_one_shared_law(page, base_url, site_app):
    """MAP.328: each sphere turns at the rate spin.js gives its stored rotation period, every one by the same elapsed
    spin time, no faster than one turn in 20 s nor slower than one in 10 minutes, a longer period never faster."""
    import math
    _open_spinning_diagram(page, base_url, site_app)
    got = page.evaluate(SPIN_SAMPLE, 800)
    before = {s["name"]: s for s in got["first"] if s["spinH"]}
    after = {s["name"]: s for s in got["second"] if s["spinH"]}
    assert len(before) >= 2, "planets with a stored period are drawn"
    elapsed = {name: (after[name]["angle"] - before[name]["angle"]) / before[name]["rate"] for name in before}
    assert all(0.3 < e < 3.0 for e in elapsed.values()), elapsed
    assert max(elapsed.values()) - min(elapsed.values()) < 1e-6, f"the bodies turned by different spin times: {elapsed}"
    for name, one in before.items():
        assert 2 * math.pi / 600 - 1e-9 <= one["rate"] <= 2 * math.pi / 20 + 1e-9, (name, one["rate"])
    ordered = sorted(before.values(), key=lambda one: one["spinH"])
    for shorter, longer in zip(ordered, ordered[1:]):
        assert shorter["rate"] >= longer["rate"] - 1e-12, (shorter, longer)


def test_the_spin_switch_stops_the_bodies_and_the_choice_is_kept(page, base_url, site_app):
    """MAP.328: the Menu's Spin switch holds every body still (and a reload keeps it so); turning it on again sets them
    turning."""
    system_id, _ = _open_spinning_diagram(page, base_url, site_app)
    moving = page.evaluate(SPIN_SAMPLE, 500)
    assert any(b["angle"] != a["angle"] for a, b in zip(moving["first"], moving["second"])), "bodies turn by default"
    page.click("#sysmap-menu > summary")
    button = page.locator('#sysmap-layers [data-setting="spin"]')
    assert button.get_attribute("aria-pressed") == "true"
    button.click()
    assert button.get_attribute("aria-pressed") == "false"
    page.wait_for_timeout(200)
    held = page.evaluate(SPIN_SAMPLE, 500)
    assert [b["angle"] for b in held["second"]] == [a["angle"] for a in held["first"]], "switched off, nothing turns"
    assert page.evaluate("() => localStorage.getItem('planetgen.spin')") == "off"
    assert "spin" not in page.url, "a standing preference, not part of the address"
    page.reload(wait_until="load")
    _show_diagram(page)
    still = page.evaluate(SPIN_SAMPLE, 500)
    assert [b["angle"] for b in still["second"]] == [a["angle"] for a in still["first"]], "the choice survives a reload"
    page.click("#sysmap-menu > summary")
    page.locator('#sysmap-layers [data-setting="spin"]').click()
    again = page.evaluate(SPIN_SAMPLE, 500)
    assert any(b["angle"] != a["angle"] for a, b in zip(again["first"], again["second"]))


def test_reduced_motion_holds_the_bodies_still(page, base_url, site_app):
    """MAP.328: with the reader's reduced-motion setting on, no body turns although the Spin switch is on."""
    _open_spinning_diagram(page, base_url, site_app, reduced="reduce")
    held = page.evaluate(SPIN_SAMPLE, 600)
    assert [b["angle"] for b in held["second"]] == [a["angle"] for a in held["first"]]
    assert page.evaluate("() => localStorage.getItem('planetgen.spin')") is None, "the switch itself was not touched"


def test_the_3d_system_view_turns_each_body_about_its_stored_axis_at_its_rate(page, base_url, site_app):
    """MAP.328: in the 3D view a body's mesh turns by the angle spin.js gives for the elapsed spin time, about its stored
    axis; a body with no stored period stays where it is."""
    import math
    system_id, scene = _system_with_moons(site_app.test_client())
    page.emulate_media(reduced_motion="no-preference")
    page.goto(f"{base_url}/system/{system_id}?view=3d", wait_until="load")
    page.wait_for_selector('#sysview3d[data-ready="true"]')
    got = page.evaluate("""async () => {
        const view = document.querySelector('#sysview3d-canvas').systemView;
        const first = view.spins();
        await new Promise((resolve) => setTimeout(resolve, 900));
        return {first: first, second: view.spins()};
    }""")
    checked = 0
    for a, b in zip(got["first"], got["second"]):
        elapsed = b["at"] - a["at"]
        assert elapsed > 0.2, elapsed
        dot = abs(sum(x * y for x, y in zip(a["quaternion"], b["quaternion"])))
        turned = 2 * math.acos(min(1.0, dot))
        if a["locked"]:
            assert turned < 0.02, (a["ref"], turned)      # UX.244: a locked body keeps its face toward what it orbits
            continue
        assert turned == pytest.approx(a["rate"] * elapsed, abs=2e-3), (a["ref"], turned, a["rate"], elapsed)
        checked += a["rate"] > 0
    assert checked >= 2, ("bodies with a stored period turn", [(a["ref"], a["rate"], a["locked"]) for a in got["first"]])


@pytest.mark.viewport(1280, 900)
def test_galaxy_map_steps_to_the_neighbouring_region_at_its_own_level(page, base_url):
    """MAP.121: the Menu's step buttons and the keys move the window to the cell beside it at the same level, in one
    camera move (the address carries it, the trail follows the new centre), and back the other way, without going
    up a level."""
    _open_galaxy(page, base_url)
    size = _descend(page)
    stage = lambda: page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyStage()")["win"]  # noqa: E731
    centre = lambda: {k: stage()["centre"][k] for k in ("m", "shell", "band", "wedge")}  # noqa: E731
    start = centre()
    assert start["m"] == size

    def enabled():
        return page.eval_on_selector_all(
            "#galaxymap3d-controls [data-action=step]",
            "els => Object.fromEntries(els.map(e => [e.dataset.step, !e.disabled]))")

    ways = enabled()
    assert set(ways) == {"north", "south", "east", "west", "in", "out"}
    # The core's own cell has nothing nearer the core, and some way out.
    assert not ways["in"] and ways["out"], ways

    # Out, with the Menu's button: the window moves to the next shell at the same size, and the address and the
    # trail follow it (no Up: the window is not a level higher or lower).
    crumbs, before = len(_crumbs(page)), _state(page)
    page.click("#galaxymap3d-menu > summary")
    page.click('#galaxymap3d-controls [data-step="out"]')
    _wait_settled(page, before)
    moved = centre()
    assert moved["m"] == size and moved["shell"] == start["shell"] + 1, (start, moved)
    assert _query(page) != before[0] and "cam=" in page.url
    assert abs(len(_crumbs(page)) - crumbs) <= 1, "the trail is the position's, not a level deeper"

    # In again, with the key: back to where it started (the core's cell is the only one a step in from shell 1 reaches).
    before = _state(page)
    page.locator("#galaxymap3d-canvas").focus()
    page.keyboard.press("[")
    _wait_settled(page, before)
    assert centre() == start

    # A step the geometry does not allow does nothing: nothing is nearer the core than the core.
    before = _state(page)
    page.keyboard.press("[")
    page.wait_for_timeout(300)
    assert _state(page) == before and centre() == start

    # Sideways: from the shell outside the core, round the ring: a step east then west lands on the same cell.
    page.keyboard.press("]")
    _wait_settled(page)
    ring = centre()
    sideways = enabled()
    if sideways["east"]:
        before = _state(page)
        page.keyboard.press("Alt+Shift+ArrowRight")
        _wait_settled(page, before)
        east = centre()
        assert east["m"] == size and east["shell"] == ring["shell"] and east != ring, (ring, east)
        before = _state(page)
        page.keyboard.press("Alt+Shift+ArrowLeft")
        _wait_settled(page, before)
        assert centre() == ring
    else:
        assert ways["out"]


def test_a_locked_moon_in_the_system_map_keeps_its_day_side_toward_its_planet(page, base_url, site_app):
    """UX.245: in a planet's moon scene a tidally locked moon is turned to face the planet (the scene's own marker) instead of
    turning by the display law: the direction it faces is the one from the moon's marker to the planet's."""
    client = site_app.test_client()
    found = None
    for system in client.get("/api/systems?limit=100").get_json()["items"]:
        scene = client.get(f"/api/systems/{system['id']}/scene").get_json()
        for planet in scene["planets"]:
            if any(moon["look"].get("locked") for moon in planet.get("moons", [])):
                found = (system["id"], planet["id"])
                break
        if found:
            break
    if not found:
        pytest.skip("no generated moon is tidally locked")
    system_id, planet_id = found
    page.emulate_media(reduced_motion="no-preference")
    _open(page, f"{base_url}/system/{system_id}?scene=planet-{planet_id}", "#sysmap-root")
    page.wait_for_function("() => { const c = document.querySelector('#sysmap-spheres-canvas');"
                           " return c && c.__spins && c.__spins().some((s) => s.kind === 'moon' && s.locked && s.toward); }")
    spins = page.evaluate("() => document.querySelector('#sysmap-spheres-canvas').__spins()")
    parent = next(s for s in spins if s["isSelf"])
    for moon in (s for s in spins if s["kind"] == "moon" and s["locked"] and s["toward"]):
        want = (parent["cx"] - moon["cx"], -(parent["cy"] - moon["cy"]))
        assert moon["toward"][0] == pytest.approx(want[0], abs=1e-6) and moon["toward"][1] == pytest.approx(want[1], abs=1e-6), (moon, want)
