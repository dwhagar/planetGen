"""
The maps in a real browser: every map button does something (TEST.55),
the System Map's selection and measuring (TEST.58's browser half: it is
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
    sample_job_tree,
    sample_params,
    site_app,
    site_db,
)

VIEWPORT = {"width": 1280, "height": 900}
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
def page(gl_browser):
    context = gl_browser.new_context(viewport=VIEWPORT, device_scale_factor=1, reduced_motion="reduce",
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


# --- TEST.55: the Sector Map ---------------------------------------------------------

def _scale(page, selector):
    return page.locator(selector).inner_text()


def test_sector_map_buttons_change_the_view(page, base_url, sample_params):
    _open(page, f"{base_url}/sector/{sample_params['sector_id']}", "#starmap-canvas")
    canvas = "#starmap-canvas"
    first = _shot(page, canvas)
    scale = _scale(page, "#starmap-scale-label")
    actions = page.eval_on_selector_all("#starmap-controls [data-action]", "els => els.map(e => e.dataset.action)")
    assert {"zoom-in", "zoom-out", "reset"} <= set(actions)

    page.click('#starmap-controls [data-action="zoom-in"]')
    zoomed = _shot(page, canvas)
    assert zoomed != first, "+ didn't change the Sector Map"

    page.click('#starmap-controls [data-action="zoom-out"]')
    page.click('#starmap-controls [data-action="zoom-out"]')
    out = _shot(page, canvas)
    assert out != zoomed, "- didn't change the Sector Map"

    page.click('#starmap-controls [data-action="reset"]')
    assert _shot(page, canvas) != out, "Reset view didn't change the Sector Map"
    assert _scale(page, "#starmap-scale-label") == scale, "Reset view is back at the opening scale"

    if "toggle-rogue-markers" in actions:
        toggle = page.locator('#starmap-controls [data-action="toggle-rogue-markers"]')
        before = _shot(page, canvas)
        assert toggle.get_attribute("aria-pressed") == "false", "off by default (MAP.83)"
        toggle.click()
        assert toggle.get_attribute("aria-pressed") == "true"
        assert _shot(page, canvas) != before


def test_sector_map_drag_and_click(page, base_url, sample_params):
    _open(page, f"{base_url}/sector/{sample_params['sector_id']}", "#starmap-canvas")
    canvas = page.locator("#starmap-canvas")
    box = canvas.bounding_box()
    first = _shot(page, "#starmap-canvas")
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    page.mouse.down()
    page.mouse.move(box["x"] + box["width"] / 2 + 80, box["y"] + box["height"] / 2 + 20, steps=5)
    page.mouse.up()
    assert _shot(page, "#starmap-canvas") != first, "dragging didn't turn the Sector Map"
    # The screen-reader list selects too (a click on the canvas is the
    # node tests' job, tests/js/sectormap.test.mjs).
    first_entry = page.locator(".starmap-sr-list button").first
    name = first_entry.inner_text()
    first_entry.evaluate("b => b.click()")
    assert page.locator("#starmap-info h3").inner_text() == name


# --- TEST.55, TEST.58: the System Map ------------------------------------------------

def _active_scene(page):
    return page.evaluate("""() => {
        const s = document.querySelector('#sysmap-root .sysmap-svg:not(.sysmap-orbits-layer):not(.sysmap-hidden)');
        return s ? s.dataset.scene : null;
    }""")


def test_system_map_selection_drill_and_measure(page, base_url, site_app):
    client = site_app.test_client()
    systems = client.get("/api/systems?limit=100").get_json()["items"]
    target = None
    for system in systems:
        resp = client.get(f"/system/{system['id']}")
        if b'data-scene="' in resp.data and re.search(rb'data-kind="planet"[^>]*data-scene="', resp.data):
            target = system["id"]
            break
    assert target, "no generated system has a planet with moons"
    _open(page, f"{base_url}/system/{target}", "#sysmap-root")
    assert _active_scene(page) == "system"

    # A body shows its details.
    body = page.locator('#sysmap-root .sysmap-svg[data-scene="system"] [data-kind="star"]').first
    body.click()
    assert page.locator("#sysmap-info h3").count() == 1

    # A planet with moons opens its moons; the crumb's button comes back.
    planet = page.locator('#sysmap-root .sysmap-svg[data-scene="system"] [data-kind="planet"][data-scene]').first
    moons = planet.get_attribute("data-scene")
    planet.click()
    assert _active_scene(page) == moons
    assert "Moons of" in page.locator("#sysmap-crumb").inner_text()
    # MAP.92: the planet and its moons show their radius and mass.
    info = page.locator("#sysmap-info").inner_text()
    assert "Radius" in info and "Mass" in info, info
    page.locator(f'#sysmap-root .sysmap-svg[data-scene="{moons}"] [data-kind="moon"]').first.click()
    info = page.locator("#sysmap-info").inner_text()
    assert "Orbits" in info and "Radius" in info and "Mass" in info, info
    page.click("#sysmap-crumb .sysmap-back-btn")
    assert _active_scene(page) == "system"
    assert page.locator("#sysmap-crumb button").count() == 0

    # Measure: the button toggles, two bodies give a distance and a path.
    measure = page.locator("#sysmap-measure-btn")
    measure.click()
    assert measure.get_attribute("aria-pressed") == "true"
    assert "Click two" in page.locator("#sysmap-info").inner_text()
    markers = page.locator('#sysmap-root .sysmap-svg[data-scene="system"] [data-kind="planet"]')
    assert markers.count() >= 2
    markers.nth(0).click()
    markers.nth(1).click()
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


@pytest.mark.xfail(strict=True, reason=(
    "UX.21: a nebula or remnant about half a light-year across or more opens at the 1 ly "
    "zoom-out limit (planetgen/web/maps/phenomenonmap.py), so - does nothing"))
def test_phenomenon_diagram_minus_zooms_out_from_the_start(page, base_url, site_app):
    """A large nebula's "-" (a smaller object's, whose view opens inside
    the limit, works: tests/js/mapzoom.test.mjs and the "- after +" step
    above)."""
    nebula = next(p for p in _phenomena(site_app) if p["type"] == "nebula")
    assert _diagram(page, base_url, nebula)
    start = _view_box(page)
    assert start[2] >= float(page.get_attribute("#phenomenonmap-svg", "data-max-view-size")) * (1 - 1e-6), \
        "the generated nebula is too small to show the bug"
    page.click('#phenomenonmap-controls [data-action="zoom-out"]')
    assert _view_box(page)[2] > start[2], "- does nothing on the nebula diagram"


# --- TEST.55, TEST.59: the Galaxy Map ------------------------------------------------

def _crumbs(page):
    """Every step to here, galaxy first: the Steps menu's list, which
    holds them all however much of the breadcrumb line is folded into
    "…" (MAP.93, MAP.94)."""
    return page.eval_on_selector_all("#galaxymap3d-steps [data-steps-panel] li",
                                     "els => els.map(e => e.textContent.trim())")


def _query(page):
    return urlparse(page.url).query


def _state(page):
    return page.url, _crumbs(page)


def _wait_settled(page, before=None, timeout_s=15.0):
    """Waits until the URL and the steps have held still for 250 ms
    (and, given `before`, moved off it). TEST.89: a fixed 250 ms wait
    read the map mid-change under a loaded parallel run."""
    deadline = time.monotonic() + timeout_s
    page.wait_for_timeout(100)
    last, since = _state(page), time.monotonic()
    while time.monotonic() < deadline:
        page.wait_for_timeout(50)
        now = _state(page)
        if now != last:
            last, since = now, time.monotonic()
        elif (before is None or now != before) and time.monotonic() - since >= 0.25:
            return


def _hover_choice(page, wanted=None):
    """Moves the pointer over the map until the tooltip names a choice
    (matching `wanted`, a regex); returns (x, y, text) or None."""
    box = page.locator("#galaxymap3d-canvas").bounding_box()
    tooltip = page.locator("#galaxymap3d-tooltip")
    for fy in [i / 24 for i in range(2, 23)]:
        for fx in [i / 24 for i in range(2, 23)]:
            x, y = box["x"] + fx * box["width"], box["y"] + fy * box["height"]
            page.mouse.move(x, y)
            if tooltip.is_visible():
                text = tooltip.inner_text()
                if text and (wanted is None or re.search(wanted, text)):
                    return x, y, text
    return None


def _open_galaxy(page, base_url, query=""):
    _open(page, f"{base_url}/galaxy{query}", "#galaxymap3d-canvas")
    page.wait_for_selector("#galaxymap3d-steps [data-steps-panel] li", state="attached")
    _wait_settled(page)


def _click_choice(page, wanted=None):
    found = _hover_choice(page, wanted)
    assert found, f"nothing to click on the Galaxy Map ({wanted}); at {_crumbs(page)}"
    x, y, text = found
    before = _state(page)
    page.mouse.click(x, y)
    _wait_settled(page, before)
    assert _state(page) != before, f"clicking {text!r} did nothing"
    return text


def _label_of(tooltip):
    """The breadcrumb label a tooltip's choice leads to: its text up to
    the first ': ' or ', ' (a sector's own label has commas in its
    numbers, "Sector 3,837·14·12,287")."""
    return re.split(r"[:,] ", tooltip, maxsplit=1)[0].strip()


def test_galaxy_map_drill_down_by_clicks(page, base_url):
    _open_galaxy(page, base_url)
    assert _crumbs(page) == ["Galaxy"]
    assert _query(page) == ""
    steps = [("", ["Galaxy"])]
    reached_sector = False
    for _ in range(16):
        text = _click_choice(page)
        crumbs = _crumbs(page)
        query = _query(page)
        assert crumbs[-1] == _label_of(text), f"clicked {text!r}, breadcrumb now {crumbs}"
        if crumbs[-1].startswith("Sector "):
            reached_sector = True
            break
        assert len(crumbs) == len(steps[-1][1]) + 1 or crumbs[-1].startswith("Block "), crumbs
        assert re.fullmatch(r"(at=\d+\.\d+\.\d+\.-?\d+)?(&?p=[asL0-9.~,-]+)?", query), query
        assert query != steps[-1][0], "the URL names the new stage"
        steps.append((query, crumbs))
    assert reached_sector, f"never reached a sector: {steps[-1]}"
    kinds = " ".join(q for q, _ in steps)
    # MAP.56: arc, slab, segment, slab, ... down to a sector; no regions.
    for marker in ("p=a", "L", "at=243.", "at=27.", "at=3."):
        assert marker in kinds, f"no {marker} stage on the way down: {[q for q, _ in steps]}"
    assert not re.search(r"[=,]r\d", kinds), kinds

    # Back and Forward (the browser's) walk the same stages, in order.
    # (Picking the sector may or may not add an entry of its own, so the
    # walk is matched against the stages seen rather than counted.)
    steps.append((_query(page), _crumbs(page)))
    queries = [q for q, _ in steps]
    at = len(steps) - 1
    for _ in range(3):
        before = _state(page)
        page.go_back()
        _wait_settled(page, before)
        query = _query(page)
        assert query in queries[:at], f"Back went to {query!r}, not an earlier stage of {queries[:at]}"
        at = max(n for n in range(at) if queries[n] == query)
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
    assert _crumbs(page) == crumbs[:-1], "Up drops the last breadcrumb"

    # A stage's URL opens it directly.
    deep_query, deep_crumbs = steps[-2]
    _open_galaxy(page, base_url, "?" + deep_query)
    assert _crumbs(page) == deep_crumbs

    # A breadcrumb button goes back up; Reset goes home.
    more = page.locator("#galaxymap3d-crumbs .galaxy-crumb-menu")
    if more.count():  # the second step is folded into "…" (MAP.93)
        more.locator("summary").click()
        more.locator("button").first.click()
    else:
        page.locator("#galaxymap3d-crumbs button.galaxy-crumb").nth(1).click()
    _wait_settled(page)
    assert _crumbs(page) == steps[1][1]
    before = _state(page)
    page.click('#galaxymap3d-controls [data-action="reset"]')
    _wait_settled(page, before)
    assert _crumbs(page) == ["Galaxy"]
    assert _query(page) == ""


def _scale_text(page):
    """The scale line: its label and its bar's width (a zoom between two
    nice lengths only changes the bar)."""
    return page.locator("#galaxymap3d-scale").inner_html()


def test_galaxy_map_free_camera_from_an_arc_down(page, base_url):
    _open_galaxy(page, base_url)
    reset_view = page.locator('#galaxymap3d-controls [data-action="reset-view"]')
    assert not reset_view.is_disabled(), "the whole galaxy turns (MAP.85)"
    _click_choice(page)
    assert parse_qs(_query(page))["p"][0].startswith("a"), "the galaxy's first pick is an arc"

    canvas = "#galaxymap3d-canvas"
    box = page.locator(canvas).bounding_box()
    middle = (box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    first = _shot(page, canvas)
    scale = _scale_text(page)
    stage = (page.url, _crumbs(page))

    page.mouse.move(*middle)
    page.mouse.down()
    page.mouse.move(middle[0] + 120, middle[1] + 40, steps=6)
    page.mouse.up()
    turned = _shot(page, canvas)
    assert turned != first, "dragging didn't turn the view"
    assert (page.url, _crumbs(page)) == stage, "a drag picks nothing"

    page.mouse.move(*middle)
    page.mouse.wheel(0, 400)
    page.wait_for_timeout(150)
    assert _scale_text(page) != scale, "the wheel didn't zoom"

    page.mouse.move(*middle)
    page.mouse.down(button="right")
    page.mouse.move(middle[0] - 150, middle[1], steps=6)
    page.mouse.up(button="right")
    assert _shot(page, canvas) != turned, "right-drag didn't move the view"

    page.click("#galaxymap3d-menu summary")
    reset_view.click()
    _wait_settled(page)
    assert _scale_text(page) == scale, "Reset view is back at the stage's own zoom"
    assert (page.url, _crumbs(page)) == stage


def test_galaxy_map_buttons(page, base_url):
    _open_galaxy(page, base_url)
    controls = page.locator("#galaxymap3d-controls")
    actions = controls.locator("[data-action]").evaluate_all("els => els.map(e => e.dataset.action)")
    canvas = "#galaxymap3d-canvas"
    for action in ("back", "forward", "up"):
        assert controls.locator(f'[data-action="{action}"]').is_disabled(), f"{action} at the galaxy"
    assert "wedges" not in actions, "no Wedges button (MAP.85)"
    page.click("#galaxymap3d-menu summary")

    if "generated-only" in actions:
        only = controls.locator('[data-action="generated-only"]')
        before = _shot(page, canvas)
        only.click()
        assert only.get_attribute("aria-pressed") == "true"
        assert _shot(page, canvas) != before, "Generated only didn't change the map"
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

    # One step down: Back and Up work; Back then enables Forward.
    _click_choice(page)
    for action in ("back", "up"):
        assert not controls.locator(f'[data-action="{action}"]').is_disabled()
    controls.locator('[data-action="back"]').click()
    _wait_settled(page)
    assert _crumbs(page) == ["Galaxy"]
    assert not controls.locator('[data-action="forward"]').is_disabled()


# --- NAV.40: picking a course end from a bookmark ---------------------------------

def _seed_bookmarks(page, entries):
    page.evaluate("""(entries) => {
        const db = document.querySelector("[data-bookmark-db]").getAttribute("data-bookmark-db");
        localStorage.setItem("planetgen.bookmarks." + db, JSON.stringify(entries));
    }""", entries)


def _course_ends(page):
    params = parse_qs(_query(page))
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
    _open(page, f"{base_url}/sector/{sample_params['sector_id']}?pick=to&from={start}", "#starmap-canvas")
    page.locator(".pick-bookmarks summary").click()
    with page.expect_navigation():
        page.locator(".pick-bookmarks a", has_text="Bookmarked destination").click()
    assert _course_ends(page) == (start, dest)
