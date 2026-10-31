"""
The Galaxy Map and the Sector Map in a real browser with no database
(TEST.70): picking, hover, keys, Back/Forward and URL state, bookmarks and
the scale line, pinned before the shared map code (MAP.61: MAP.63's
helpers, MAP.64's camera and input controller) moves them, so a refactor
can't change what they do unnoticed.

The pages come from `tests/map_site_support.py`'s fixture site (the real
Flask views over fixture data) in headless Chromium with SwiftShader
WebGL, so they run wherever Playwright and Chromium are installed, MySQL
or not. `test_web_browser_maps.py` covers the same maps against a
generated database.
"""

import time
import io
import math
import re
from urllib.parse import parse_qs, urlparse

import pytest

sync_api = pytest.importorskip("playwright.sync_api")
Image = pytest.importorskip("PIL.Image")

from planetgen.galaxy.spherical import provisional_sector_designation  # noqa: E402
from tests.map_site_support import (  # noqa: E402,F401 -- fixtures
    DENSE_CENTER_PC, DENSE_HALO_RADIUS_PC, DENSE_POINTS, LIFE_FIELD, EDGE_PC, DENSE_STARS, DENSE_WORLDS, HALO_CENTER_PC, HALO_POINTS, HALO_STARS, PHENOMENA, POINT_FIELD,
    BARREN_FIELD, BESIDE_ID, BESIDE_SYSTEMS, GIANT_FIELD, SECTOR_ID, SECTORS, STAR_FIELD_CENTER_PC, SYSTEMS, UNCHARTED_CELL, WORLD_FIELD, _sector_center_pc,
    dense_region_query, dense_site, HOLE_FIELD, FAR_BUBBLE,
    far_site, halo_site, holes_site,
    map_site,
)

PRIME_CODE = provisional_sector_designation(1, 0, 1)
"""The designation of the fixture's Fixture Prime (shell 1, band 0, wedge 1)."""
BESIDE_CODE = provisional_sector_designation(1, 0, 2)
"""The designation of Fixture Beside, the sector next to it (MAP.329)."""
UNCHARTED_CODE = provisional_sector_designation(*UNCHARTED_CELL)
"""The designation of the fixture's uncharted cell."""

VIEWPORT = {"width": 1000, "height": 800}
"""dict: The default window. Reading a software-WebGL canvas back costs time in proportion to its pixels (2.5 s
at 1280x900, 1.1 s at 800x600, TEST.139), and most of these tests take several screenshots; a test that needs
a wide or tall window says so with `@pytest.mark.viewport(width, height)`. It is 1000 wide and 800 high because the
page's header takes the top 230 px or more and a map in a 600 px window ran off its bottom (TEST.394)."""
HIDE_NAMES = "label-system,label-star,label-object"
"""str: The address's hide list that takes the names off the map. A name takes the click or the hover that lands on it
(MAP.298: the map picks a name before the star under it), so a test that points at a star or a marker and is not about
the names is marked `names_hidden` and opens its maps with them hidden (`_open`); `_raw_query` leaves that token out."""
GL_ARGS = ["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"]
SECTOR_CANVAS = "#galaxymap3d-canvas"
GALAXY_CANVAS = "#galaxymap3d-canvas"
SYSTEM_NAMES = {name for _id, name, *_rest in SYSTEMS}


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
    page.names_hidden = request.node.get_closest_marker("names_hidden") is not None
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    try:
        yield page
    finally:
        context.close()
    assert not errors, f"uncaught exceptions: {errors}"

DOUBLE_CLICK_GAP_MS = 400
"""int: Longer than the map's double-click window (350 ms), so two clicks on one spot stay two clicks."""


def _open(page, url, ready):
    if getattr(page, "names_hidden", False) and "hide=" not in url:
        url += ("&" if "?" in url else "?") + "hide=" + HIDE_NAMES
    response = page.goto(url, wait_until="load")
    assert response.status == 200, f"{url} answered {response.status}"
    page.wait_for_selector(ready, state="attached")
    if not page.evaluate("!!document.createElement('canvas').getContext('webgl2')"):
        pytest.skip("no WebGL in this Chromium")
    page.wait_for_timeout(300)


def _shot(page, selector):
    """The picture of `selector` without the names laid over it (the label engine's, MAP.231, fade and move), the ruler's words
    (MAP.336: their antialiasing differs by a level between a canvas with the focus and one without), the Key (MAP.234: its lines
    follow what is drawn, so they change as the view settles) and the hover preview (MAP.348)."""
    page.wait_for_timeout(150)
    covers = ".map-labels, .starmap-ruler-stops, .starmap-altitude, #galaxymap3d-layer-legend, .map-hover-preview"
    page.evaluate("(q) => document.querySelectorAll(q).forEach((e) => { e.style.visibility = 'hidden'; })", covers)
    try:
        return page.locator(selector).screenshot()
    finally:
        page.evaluate("(q) => document.querySelectorAll(q).forEach((e) => { e.style.visibility = ''; })", covers)


# --- The Sector Map -------------------------------------------------------------------

def _open_sector(page, base, query=""):
    _open(page, f"{base}/sector/{SECTOR_ID}{query}", SECTOR_CANVAS)
    page.wait_for_selector(".starmap-sr-list button", state="attached")


def _selection_displays(page):
    """The ids of the selection displays drawn on the page right now: a panel or box that is on the screen and holds
    a picked thing's heading (MAP.256: one, never two)."""
    return page.evaluate("""() => [...document.querySelectorAll('.starmap-info, .map-infobox')]
        .filter((e) => e.getClientRects().length && e.querySelector('h3')).map((e) => e.id)""")


def _expand_card(page, wait=True):
    """Presses Expand on the info card (MAP.347), where it has one and is not open: what is behind it is drawn only then. The button
    shows once the card has something behind it (the region line and the nearest systems arrive a moment after the card), so by
    default this waits up to 20 seconds for it; `wait=False` presses it only if it is there."""
    button = page.locator("#galaxymap3d-info .map-info-expand")
    if wait:
        try:
            button.wait_for(state="visible", timeout=20000)
        except sync_api.TimeoutError:
            return
    if button.count() and button.is_visible() and button.get_attribute("aria-expanded") != "true":
        button.evaluate("b => b.click()")    # a Menu left open lies over the card's corner


def _info_title(page):
    heading = page.locator("#galaxymap3d-info h3")
    return heading.inner_text() if heading.count() else None


def _bright_spots(page, selector, limit=40, white=False):
    """Screen points (page coordinates) of the brightest small spots on a
    canvas, brightest first: where its points of light are drawn (with
    `white`, only near-white ones: a star's whitened core, not the Galaxy
    Map's lit blocks)."""
    # The map's names (MAP.268) are text over the canvas, not points of light: leave them out of the picture.
    page.evaluate("() => document.querySelectorAll('.map-labels').forEach((e) => { e.style.visibility = 'hidden'; })")
    try:
        image = Image.open(io.BytesIO(_shot(page, selector))).convert("RGB")
    finally:
        page.evaluate("() => document.querySelectorAll('.map-labels').forEach((e) => { e.style.visibility = ''; })")
    box = page.locator(selector).bounding_box()  # after the shot, which may scroll
    width, height = image.size
    pixels = image.load()
    background = pixels[2, 2]

    def contrast(x, y):
        return sum(abs(a - b) for a, b in zip(pixels[x, y], background))

    spots = []
    for y in range(2, height - 2):
        for x in range(2, width - 2):
            here = contrast(x, y)
            if here < 120 or (white and min(pixels[x, y]) < 225):
                continue
            if all(here >= contrast(x + dx, y + dy) for dx in (-2, 0, 2) for dy in (-2, 0, 2)):
                spots.append((here, x, y))
    spots.sort(reverse=True)
    kept = []
    for _value, x, y in spots:
        if all(abs(x - kx) + abs(y - ky) > 8 for kx, ky in kept):
            kept.append((x, y))
        if len(kept) >= limit:
            break
    return [(box["x"] + x * box["width"] / width, box["y"] + y * box["height"] / height) for x, y in kept]


def _click_a_star(page):
    """Clicks points of light until one opens a star system's details;
    returns (x, y, name)."""
    for x, y in _bright_spots(page, SECTOR_CANVAS):
        page.mouse.click(x, y)
        name = _info_title(page)
        if name in SYSTEM_NAMES:
            return x, y, name
    pytest.fail("no click on the Sector Map picked a star")


def test_sector_map_click_picks_a_star(page, map_site):
    _open_sector(page, map_site)
    assert _info_title(page) is None
    before = _shot(page, SECTOR_CANVAS)
    _x, _y, name = _click_a_star(page)
    system_id = next(i for i, n, *_rest in SYSTEMS if n == name)
    info = page.locator("#galaxymap3d-info")
    assert "Star type" in info.inner_text()
    view = info.locator("a", has_text="View system")
    assert view.count() == 1 and view.get_attribute("href").endswith(f"/system/{system_id}")
    nav = [b.inner_text() for b in info.locator("button").all()]
    assert "Start Here" in nav and "End Here" in nav
    assert _shot(page, SECTOR_CANVAS) != before, "the picked star gets a highlight ring"
    # MAP.256: the selection is shown once, over the map window, with a Close button; there is no second box.
    assert _selection_displays(page) == ["galaxymap3d-info"]
    assert page.locator("#galaxymap3d-infobox").count() == 0
    assert info.locator("h3").inner_text() == name
    info.locator("button", has_text="Close").click()
    assert _selection_displays(page) == []


@pytest.mark.names_hidden
def test_sector_map_hover_names_what_is_under_the_pointer(page, map_site):
    """MAP.65: the Sector Map shows a tooltip for what is under the pointer
    (it had none before the shared picking layer), and a picked system's
    panel has a ☆ Bookmark button."""
    _open_sector(page, map_site)
    x, y, name = _click_a_star(page)
    tip = page.locator("#galaxymap3d-tooltip")
    box = page.locator(SECTOR_CANVAS).bounding_box()
    page.mouse.move(x, y)
    assert tip.is_visible()
    assert tip.inner_text().startswith(name), tip.inner_text()
    star = page.locator("#galaxymap3d-info .map-info-bookmark")
    assert star.inner_text() == "☆" and star.is_enabled()
    page.mouse.move(box["x"] - 20, box["y"] - 20)
    assert not tip.is_visible(), "leaving the map hides it"


def test_sector_map_click_on_empty_space_keeps_the_selection(page, map_site):
    _open_sector(page, map_site)
    _x, _y, name = _click_a_star(page)
    box = page.locator(SECTOR_CANVAS).bounding_box()
    page.mouse.click(box["x"] + 3, box["y"] + 3)
    assert _info_title(page) == name


def test_sector_map_drag_turns_and_picks_nothing(page, map_site):
    _open_sector(page, map_site)
    x, y, name = _click_a_star(page)
    page.locator("#galaxymap3d-info").evaluate("panel => panel.textContent = ''")
    before = _shot(page, SECTOR_CANVAS)
    page.mouse.move(x, y)
    page.mouse.down()
    page.mouse.move(x + 60, y + 15, steps=6)
    page.mouse.up()
    assert _info_title(page) is None, "a drag that starts on a star doesn't pick it"
    assert _shot(page, SECTOR_CANVAS) != before, "dragging turns the map"


def test_sector_map_arrow_keys_turn_the_view_and_keep_the_scale(page, map_site):
    _open_sector(page, map_site)
    _settle(page)    # the opening flight (MAP.333) lands after the page is ready: the scale reads 10 pc before it and 1 pc after
    canvas = page.locator(SECTOR_CANVAS)
    scale = page.locator("#galaxymap3d-scale .starmap-scale-label").inner_text()
    canvas.focus()
    for key in ("ArrowLeft", "ArrowRight", "ArrowUp", "ArrowDown"):
        before = _shot(page, SECTOR_CANVAS)
        page.keyboard.press(key)
        assert _shot(page, SECTOR_CANVAS) != before, f"{key} didn't turn the Sector Map"
    assert page.locator("#galaxymap3d-scale .starmap-scale-label").inner_text() == scale, "turning doesn't change the scale"
    _settle(page)    # the last turn and the nodes it asked for have landed: a frame still changing is not the key's doing (TEST.426)
    before = _shot(page, SECTOR_CANVAS)
    for _ in range(20):    # the turn eases and the sector's scene may still be arriving: the picture is taken once it holds still
        now = _shot(page, SECTOR_CANVAS)
        if now == before:
            break
        before = now
    page.keyboard.press("a")
    assert _shot(page, SECTOR_CANVAS) == before, "other keys do nothing"


def _scale_line(page):
    """(label, bar width in px) of the map's scale line."""
    width = page.locator("#galaxymap3d-scale .starmap-scale-bar").evaluate("bar => parseFloat(bar.style.width)")
    return page.locator("#galaxymap3d-scale .starmap-scale-label").inner_text(), width


def _pc(label):
    """The parsecs of a scale label, "<n> sectors · <distance>"."""
    number, unit = re.search(r"· ([\d.,]+) (\w+)", label).groups()
    return float(number.replace(",", "")) * {"pc": 1, "kpc": 1000, "ly": 1 / 3.26156, "mly": 1000 / 3.26156}[unit]


def test_sector_map_scale_line_follows_the_zoom(page, map_site):
    _open_sector(page, map_site)
    _settle(page)  # TEST.458: the opening zoom is the one the page's flight in lands at, not a frame of it
    label, width = _scale_line(page)
    assert re.search(r"sectors? · [\d.,]+ (pc|kpc|ly|mly)", label), label
    assert 30 <= width <= 160, "the scale line is about 70 px long"
    ly_per_px = _pc(label) / width

    zoom_in = page.locator('#galaxymap3d-controls [data-action="zoom-in"]')
    zoom_out = page.locator('#galaxymap3d-controls [data-action="zoom-out"]')
    zoom_in.click()
    in_label, in_width = _scale_line(page)
    assert _pc(in_label) / in_width < ly_per_px, "+ shows fewer light-years per pixel"

    page.locator(SECTOR_CANVAS).scroll_into_view_if_needed()
    box = page.locator(SECTOR_CANVAS).bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, min(box["y"] + box["height"] / 2, VIEWPORT["height"] - 100))
    page.mouse.wheel(0, 300)
    page.wait_for_timeout(300)
    wheel_label, wheel_width = _scale_line(page)
    assert _pc(wheel_label) / wheel_width > _pc(in_label) / in_width, "scrolling down zooms out"

    for _ in range(30):
        zoom_in.click()
    closest = _scale_line(page)
    zoom_in.click()
    assert _scale_line(page) == closest, "zoom in stops at its limit"
    for _ in range(30):
        zoom_out.click()
    farthest = _scale_line(page)
    zoom_out.click()
    assert _scale_line(page) == farthest, "zoom out stops at its limit"

    page.locator("#galaxymap3d-menu > summary").click()
    page.locator('#galaxymap3d-controls [data-action="reset-view"]').click()
    for _ in range(40):  # the view eases back; a loaded machine takes longer than a fixed pause
        page.wait_for_timeout(250)
        if _scale_line(page) == (label, width):
            break
    assert _scale_line(page) == (label, width), "Reset view goes back to the opening zoom"


def test_reset_view_after_the_minus_button_has_run_out_of_range_goes_back(page, map_site):
    """TEST.458: pressing minus past the end of the sector's range leaves a distance goal being chased; Reset view ended it no more than
    a fly-to did, so the chase took the view back out as soon as the reset's flight was over and the button seemed to do nothing."""
    _open_sector(page, map_site)
    _settle(page)
    opening = _scale_line(page)
    zoom_out = page.locator('#galaxymap3d-controls [data-action="zoom-out"]')
    for _ in range(30):
        zoom_out.click()
    assert _scale_line(page) != opening
    page.locator("#galaxymap3d-menu > summary").click()
    page.locator('#galaxymap3d-controls [data-action="reset-view"]').click()
    for _ in range(40):
        page.wait_for_timeout(250)
        if _scale_line(page) == opening:
            break
    page.wait_for_timeout(1500)
    assert _scale_line(page) == opening, "Reset view goes back to the opening zoom and stays there"


def test_sector_map_show_on_map_selects_a_rogue_planet(page, map_site):
    # The Contents table folds rogue planets into one row; filtered to them it lists each.
    _open_sector(page, map_site, "?contents_type=Rogue+Planet")
    rogue = next(p for p in PHENOMENA if p["type"] == "rogue_planet")
    button = page.locator(f'[data-map-target="rogue_planet:{rogue["id"]}"]')
    button.wait_for(state="attached")
    assert button.count() == 1
    button.click()
    assert _info_title(page) == rogue["name"]
    assert page.evaluate("document.activeElement.id") == "galaxymap3d-canvas"
    view = page.locator("#galaxymap3d-info a", has_text="View phenomenon")
    assert view.get_attribute("href").endswith(f"/phenomenon/rogue_planet/{rogue['id']}")


def test_sector_map_has_no_rogue_planet_marker_button(page, map_site):
    """MAP.184: the rogue-only "Mark rogue planets" button is gone; the Menu's per-kind Mark covers it."""
    _open_sector(page, map_site)
    assert page.locator('#galaxymap3d-controls [data-action="toggle-rogue-markers"]').count() == 0
    assert "Mark rogue planets" not in page.locator("#galaxymap3d-controls").inner_text()


def _open_menu(page):
    if not page.locator("#galaxymap3d-menu").evaluate("m => m.open"):
        page.locator("#galaxymap3d-menu > summary").click()


POI_LAYERS = ("marker-inhabited", "marker-habitable", "marker-uninhabitable", "marker-gas-giant", "marker-asteroid-field")
POI_DEFAULT_ON = ("marker-inhabited", "marker-habitable", "marker-gas-giant")
"""The Menu's point-of-interest layers, and the ones that are on until the visitor says otherwise (MAP.303)."""


def _poi_layers(page, *on):
    """Opens the Menu and sets the point-of-interest layers so that exactly `on` are on."""
    _open_menu(page)
    for layer in POI_LAYERS:
        button = page.locator(f'#galaxymap3d-kinds button[data-layer="{layer}"]')
        if (button.get_attribute("aria-pressed") == "true") != (layer in on):
            button.click()


def test_sector_map_kinds_can_be_hidden_one_by_one_and_the_choice_is_kept(page, map_site):
    """MAP.79: a toggle per kind of object in the Menu; a hidden kind is gone
    from the map, can't be hovered or picked, and its choice is in the URL."""
    _open_sector(page, map_site)
    _x, _y, name = _click_a_star(page)
    page.locator("#galaxymap3d-info").evaluate("panel => panel.textContent = ''")
    before = _shot(page, SECTOR_CANVAS)
    _open_menu(page)
    buttons = page.locator("#galaxymap3d-kinds button[data-kind]")
    names = buttons.all_inner_texts()
    assert "Stars" in names and "Rogue planets" in names, names
    assert [b.get_attribute("aria-pressed") for b in buttons.all()] == [
        "false" if name == "Rogue planets" else "true" for name in names], "all on but the rogue planets (MAP.137)"
    stars = page.locator('#galaxymap3d-kinds button[data-kind="star"]')
    stars.click()
    assert stars.get_attribute("aria-pressed") == "false"
    assert "hide=star" in page.url
    page.locator("#galaxymap3d-menu > summary").click()   # the open Menu covers the spot of the star
    page.mouse.click(_x, _y)
    assert _info_title(page) != name, "a hidden star can't be picked"
    assert _shot(page, SECTOR_CANVAS) != before, "hiding the stars changes the map"
    assert "Stars" not in [b.inner_text() for b in page.locator("#galaxymap3d-info button").all()]
    # A reload keeps it.
    page.reload()
    page.wait_for_selector(".starmap-sr-list button", state="attached")
    page.wait_for_selector('#galaxymap3d-kinds button[data-kind="star"]', state="attached")
    _open_menu(page)
    assert page.locator('#galaxymap3d-kinds button[data-kind="star"]').get_attribute("aria-pressed") == "false"
    page.locator('#galaxymap3d-kinds button[data-kind="star"]').click()
    assert "hide" not in page.url


def test_sector_map_a_registered_layer_has_a_menu_button_a_url_key_and_a_legend_line(page, map_site):
    """MAP.217: a layer declared once in the shared registry (static/maplayers.js) shows in the Menu beside
    the kinds, writes its off state into the URL's hide list and draws its legend line while on."""
    _open_sector(page, map_site)
    page.evaluate("""() => {
        const layers = document.getElementById("galaxymap3d-canvas").__mapLayers;
        layers.register({id: "testLayer", label: "Test layer", modes: ["galaxy", "sector"], group: "Show on the map",
            urlKey: "hide", urlKeyOn: "show", legend: {text: "marks the test layer", swatch: "#f90"}});
        layers.refresh();
    }""")
    _open_menu(page)
    button = page.locator('#galaxymap3d-kinds button[data-layer="testLayer"]')
    assert button.count() == 1 and button.get_attribute("aria-pressed") == "true"
    legend = page.locator("#galaxymap3d-layer-legend")
    assert "marks the test layer" in legend.inner_text()
    button.click()
    assert button.get_attribute("aria-pressed") == "false"
    assert "hide=testLayer" in page.url
    assert "marks the test layer" not in legend.inner_text(), "an off layer has no legend line"
    button.click()
    assert "hide" not in page.url
    assert "marks the test layer" in legend.inner_text()
    # The kinds the sector stage owns still work beside it.
    assert page.locator('#galaxymap3d-kinds button[data-kind="star"]').get_attribute("aria-pressed") == "true"


def test_sector_map_show_on_map_shows_a_hidden_kind_again(page, map_site):
    _open_sector(page, map_site, "?contents_type=Rogue+Planet")
    _open_menu(page)
    rogue = next(p for p in PHENOMENA if p["type"] == "rogue_planet")
    toggle = page.locator('#galaxymap3d-kinds button[data-kind="roguePlanet"]')
    assert toggle.get_attribute("aria-pressed") == "false", "rogue planets begin hidden (MAP.137)"
    page.locator("#galaxymap3d-menu > summary").click()  # closed, so its panel isn't over the table
    button = page.locator(f'[data-map-target="rogue_planet:{rogue["id"]}"]')
    button.click()
    assert _info_title(page) == rogue["name"]
    assert toggle.get_attribute("aria-pressed") == "true", "selecting it shows its kind again"
    assert "show=roguePlanet" in page.url and "hide" not in page.url


def test_sector_map_screen_reader_list_selects(page, map_site):
    _open_sector(page, map_site)
    names = page.locator(".starmap-sr-list button").all_inner_texts()
    assert SYSTEM_NAMES <= set(names)
    assert {p["name"] for p in PHENOMENA} <= set(names)
    page.locator(".starmap-sr-list button", has_text="Fixture Pulsar").evaluate("b => b.click()")
    assert _info_title(page) == "Fixture Pulsar"


def test_sector_page_map_keeps_to_its_sector(page, map_site):
    """MAP.68: the sector page's map is the Galaxy Map's engine locked to the
    sector: no breadcrumb, steps or history of its own, and Escape and Home
    don't leave the sector."""
    _open_sector(page, map_site)
    for gone in ("#galaxymap3d-crumbs", "#galaxymap3d-slabs", "#galaxymap3d-address", '[data-action="up"]',
                 '[data-action="back"]', '[data-action="reset"]'):
        assert page.locator(gone).count() == 0, gone
    url = _raw_query(page)
    first = _shot(page, SECTOR_CANVAS)
    page.locator(SECTOR_CANVAS).focus()
    for key in ("Escape", "Backspace", "Home"):
        page.keyboard.press(key)
    assert _raw_query(page) == url
    assert _shot(page, SECTOR_CANVAS) == first, "no key leaves the sector"
    _x, _y, name = _click_a_star(page)
    assert name in SYSTEM_NAMES and _raw_query(page) == url, "picking a star adds no history entry"


# --- The Galaxy Map -------------------------------------------------------------------

def _crumbs(page):
    """Every step to here, galaxy first: the Steps menu's list, which
    holds them all however much of the breadcrumb line is folded into
    "…" (MAP.93, MAP.94)."""
    steps = page.eval_on_selector_all("#galaxymap3d-steps [data-steps-panel] li",
                                      "els => els.map(e => e.textContent.trim())")
    return steps[1:] if steps[:1] == ["Home"] else steps  # UX.59: the trail starts at Home


def _raw_query(page):
    """The address's query without the camera it carries once the view is still (MAP.228)."""
    return "&".join(part for part in urlparse(page.url).query.split("&")
                    if part and not part.startswith("cam=") and part != "hide=" + HIDE_NAMES)


def _stage(page):
    return page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyStage()")


def _query(page):
    """The address as these tests read it: its query without the camera, led by "win=m.shell.band.wedge" (with
    "&g=" when not 9) naming the window on screen when no sector is in it. The stage is no part of the address
    since MAP.228 (the camera is); the window is read from the page."""
    query = _raw_query(page)
    win = _stage(page)["win"]
    if "sector=" in query or not win:
        return query
    centre = win["centre"]
    token = f"win={centre['m']}.{centre['shell']}.{centre['band']}.{centre['wedge']}" + (f"&g={win['g']}" if win["g"] != 9 else "")
    return token + ("&" + query if query else "")


_stage_query = _query


def _open_window(page, base, m, shell, band, wedge, g=9):
    """Opens the Galaxy Map and takes the window of `g` cells round a cell, as a click on it does: the address
    holds the camera now and no stage token (MAP.228)."""
    _open_galaxy(page, base)
    page.evaluate("([cell, g]) => document.querySelector('#galaxymap3d-canvas').galaxyOpenWindow(cell, g)",
                  [{"m": m, "shell": shell, "band": band, "wedge": wedge}, g])
    _settle(page)


def _open_star_field(page, base):
    """The window round the fixture's star field, then close in on the field, as far in as its systems stand apart on the screen (marks
    closer than `MERGE_PX` are drawn as one, MAP.305)."""
    _open_window(page, base, 27, 27, 0, 0, g=3)
    page.evaluate("([c, r]) => document.querySelector('#galaxymap3d-canvas').galaxyGoTo({center: c, radius: r})",
                  [list(STAR_FIELD_CENTER_PC), 7])
    _settle(page)


def _settle(page):
    """Waits out a stage's flight: until the scale bar (its width follows
    the camera) holds still, a few seconds at most on a busy machine."""
    page.wait_for_timeout(250)
    # The map has drawn the stage it names (its data loads after the
    # breadcrumb moves), then the scale bar holds still.
    page.wait_for_function("""() => {
        const canvas = document.querySelector("#galaxymap3d-canvas");
        return !canvas || !canvas.galaxyReady || canvas.galaxyReady();
    }""", timeout=20000)
    last = None
    for _ in range(40):
        now = page.evaluate("() => (document.querySelector('#galaxymap3d-scale') || {}).innerHTML || ''")
        if now == last:
            return
        last = now
        page.wait_for_timeout(150)


# Added to every document the page loads from here on (the map reloads itself when it opens a sector from the address).
NO_WALL_ICONS = """{
    // A constructed sheet: the page's content security policy refuses a style element of ours.
    const sheet = new CSSStyleSheet();
    sheet.replaceSync(".map-wall-icons { display: none !important; }");
    document.adoptedStyleSheets = [...document.adoptedStyleSheets, sheet];
}"""


def _pin_window(page, height=434):
    """The map's window as it was before it filled the browser window (MAP.350), for the tests whose numbers were worked out on it: pinned to
    `height` pixels, the height it had at the fixture's viewport. Call it before the page opens."""
    page.add_init_script("""{
    const sheet = new CSSStyleSheet();
    sheet.replaceSync(".galaxymap3d-panel .starmap-viewport { height: %dpx !important; }");
    document.adoptedStyleSheets = [...document.adoptedStyleSheets, sheet];
}""" % height)


def _open_galaxy(page, base, query="", wall_icons=False):
    """Opens the Galaxy Map; the icons on the walls of an opened sector (MAP.327) are left out unless `wall_icons`, as
    they sit over the points the tests look for."""
    if not wall_icons:
        page.add_init_script(NO_WALL_ICONS)
    _open(page, f"{base}/galaxy{query}", GALAXY_CANVAS)
    page.wait_for_selector("#galaxymap3d-steps [data-steps-panel] li", state="attached")
    _settle(page)


HOVER_SCAN = """([wanted, steps]) => {
    const canvas = document.querySelector("#galaxymap3d-canvas");
    const tip = document.querySelector("#galaxymap3d-tooltip");
    const box = canvas.getBoundingClientRect();
    const re = wanted ? new RegExp(wanted) : null;
    for (let j = 2; j < steps - 1; j++) {
        for (let i = 2; i < steps - 1; i++) {
            const x = box.left + (i / steps) * box.width, y = box.top + (j / steps) * box.height;
            // Only where a real click can land.
            if (y < 2 || y > window.innerHeight - 2) continue;
            canvas.dispatchEvent(new PointerEvent("pointermove", {clientX: x, clientY: y, bubbles: true, pointerType: "mouse"}));
            const text = tip.hidden ? "" : tip.textContent;
            if (text && (!re || re.test(text))) return [x, y, text];
        }
    }
    return null;
}"""


NAMED_SPOTS = """([spots, names]) => {
    const canvas = document.querySelector("#galaxymap3d-canvas");
    const tip = document.querySelector("#galaxymap3d-tooltip");
    const found = [];
    for (const [x, y] of spots) {
        if (document.elementFromPoint(x, y) !== canvas) continue;  // a real pointer there would not be on the map
        canvas.dispatchEvent(new PointerEvent("pointermove", {clientX: x, clientY: y, bubbles: true, pointerType: "mouse"}));
        const text = tip.hidden ? "" : tip.textContent;
        if (names.some((name) => text.startsWith(name))) found.push([x, y, text]);
    }
    return found;
}"""
"""str: Of the screen points given, those where a hover names one of the given systems: [x, y, tooltip]."""


HOVER_RETRIES = 20


CHARTED_SCAN = """([point, heights, wanted]) => {
    const canvas = document.querySelector("#galaxymap3d-canvas");
    const tip = document.querySelector("#galaxymap3d-tooltip");
    const re = new RegExp(wanted);
    for (const height of heights) {
        const spot = canvas.galaxyScreenOf(point[0], point[1], height);
        if (!spot || document.elementFromPoint(spot.x, spot.y) !== canvas) continue;  // a real pointer there would not be on the map
        canvas.dispatchEvent(new PointerEvent("pointermove", {clientX: spot.x, clientY: spot.y, bubbles: true, pointerType: "mouse"}));
        const text = tip.hidden ? "" : tip.textContent;
        if (re.test(text)) return [spot.x, spot.y, text];
    }
    return null;
}"""
"""str: The first of `heights` (scene z, parsecs) above the scene point (x, y) where a hover names a choice
matching `wanted`: [x, y, tooltip], or None. One round trip for all of them."""


def _charted_spot(page):
    """Where the fixture's charted space (Fixture Prime) falls on the page now:
    a click or a hover there lands on the cell holding generated sectors, at
    every size of window (a block of 243 sectors is too big to find them by
    scanning: the window it opens is round the point clicked)."""
    from planetgen.galaxy.spherical import sector_position_pc
    shell, band, wedge = next((sh, b, w) for _id, sh, b, w, name in SECTORS if name == "Fixture Prime")
    x, y, z = sector_position_pc(shell, band, wedge, 4.0)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    # What the pointer meets is the first block of its ray, the top face of the
    # cell under the sector in a slanted view: aim a little under the top of
    # the cell at each size (1, 3, 9, 27, 81 and 243 sectors tall), the midplane too.
    heights = [z, z + 2 - 0.01, z + 6 - 0.01, z + 18 - 0.01, z + 54 - 0.01, z + 162 - 0.01, z + 486 - 0.01,
               z - 2 + 0.01, z - 6 + 0.01, z - 18 + 0.01, z - 54 + 0.01, z - 162 + 0.01, z - 486 + 0.01]  # and the bottoms, seen from under the plane
    for _ in range(HOVER_RETRIES):
        found = page.evaluate(CHARTED_SCAN, [[x, y], heights, GENERATED_CHOICE])
        if found:
            # The scan used synthetic pointer events, which are quick; confirm with the real pointer (TEST.109).
            page.mouse.move(found[0], found[1])
            page.wait_for_timeout(120)
            tip = page.locator("#galaxymap3d-tooltip")
            text = tip.inner_text() if tip.is_visible() else ""
            if re.search(GENERATED_CHOICE, text):
                return found[0], found[1], text
        page.wait_for_timeout(250)
    return _charted_by_key(page)


class _KeySpot(tuple):
    """A choice taken with the keyboard: `(x, y, text)` of the map's middle and the tooltip, standing for a
    choice the pointer cannot reach."""


CYCLE_TO = """([wanted, limit]) => {
    const canvas = document.querySelector('#galaxymap3d-canvas');
    const tip = document.querySelector('#galaxymap3d-tooltip');
    const pattern = new RegExp(wanted);
    canvas.focus();
    for (let n = 0; n < limit; n++) {
        canvas.dispatchEvent(new KeyboardEvent('keydown', {key: 'ArrowRight', bubbles: true, cancelable: true}));
        const text = tip.hidden ? '' : tip.textContent;
        if (pattern.test(text)) return text;
    }
    return null;
}"""


def _charted_by_key(page):
    """The charted choice when no ray of the pointer reaches it: the blocks are shells, bands and wedges that
    tile the space round the sector, so a charted one inside a window is behind others. The keyboard walks
    the window's choices in turn, hidden ones too; the choice is left selected for `_take`."""
    text = page.evaluate(CYCLE_TO, [GENERATED_CHOICE, 2000])
    if not text:
        return None
    box = page.locator(SECTOR_CANVAS).bounding_box()
    return _KeySpot((box["x"] + box["width"] / 2, box["y"] + box["height"] / 2, text))


def _take(page, found):
    """Clicks the choice `found` names, or presses Enter on one the keyboard chose."""
    if isinstance(found, _KeySpot):
        page.focus(SECTOR_CANVAS)
        page.keyboard.press("Enter")
    else:
        page.mouse.click(found[0], found[1])


def _screen_of(page, point):
    """Where a point of the scene (parsecs) falls on the page now, or None. A card open over the map is closed first (MAP.332: a card
    too tall for the room above the zoom controls now stands beside them, over more of the map, and its More takes the pointer)."""
    close = page.locator("#galaxymap3d-info .map-info-close")
    if close.count() and close.is_visible():
        close.click()
        page.wait_for_timeout(100)
    spot = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(point))
    return (spot["x"], spot["y"]) if spot else None


def _hover_choice(page, wanted=None, steps=32, retries=HOVER_RETRIES):
    """Moves the pointer over the map (as pointer events, for speed) until
    the tooltip names a choice matching `wanted` (a regex, JavaScript
    syntax); returns (x, y, text) or None, the pointer left there. The
    map ignores the pointer while a stage's flight runs, which takes longer
    on a busy machine, so the scan is tried again for a few seconds. A scan
    of the whole map costs about 2.5 s, so a test asking whether a choice
    exists at all (where "no" is an answer) settles the map first and gives
    `retries=1`."""
    if wanted == GENERATED_CHOICE:
        return _charted_spot(page)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    found = page.evaluate(HOVER_SCAN, [wanted, steps])
    for _ in range(retries):
        if found:
            break
        page.wait_for_timeout(250)
        found = page.evaluate(HOVER_SCAN, [wanted, steps])
    if found:
        page.mouse.move(found[0], found[1])
    return found


def _click_steady(page, wanted, steps=48):
    """Clicks where the tooltip names `wanted`, once that has held for a moment
    with the real pointer there (TEST.109): a scan with synthetic pointer events
    can land on a spot the map is still redrawing under (its contents table and
    shapes settle after load), where a real click picks something else. Scans
    again until the pointer, moved to the spot, still reads `wanted`."""
    for _ in range(20):
        found = _hover_choice(page, wanted, steps=steps)
        assert found, f"no spot on the map names {wanted}"
        page.mouse.move(found[0], found[1])
        page.wait_for_timeout(200)
        page.mouse.move(found[0] + 1, found[1])
        page.mouse.move(found[0], found[1])
        tooltip = page.evaluate("() => { const tip = document.querySelector('#galaxymap3d-tooltip'); return tip.hidden ? '' : tip.textContent; }")
        if re.search(wanted, tooltip):
            page.mouse.click(found[0], found[1])
            return found
        page.wait_for_timeout(250)
    pytest.fail(f"the pointer never steadily read {wanted}")


def _click_choice(page, wanted=None):
    found = _hover_choice(page, wanted)
    assert found, f"nothing to click on the Galaxy Map ({wanted}); at {_crumbs(page)}"
    x, y, text = found
    before = (_query(page), _crumbs(page))
    if re.match(r"Sector .*, charted$", text):
        # A generated sector opens in the map itself (MAP.66).
        _take(page, found)
        page.wait_for_function("() => location.search.includes('open=1')", timeout=15000)
        _settle(page)
        return text
    _take(page, found)
    _settle(page)
    assert (_query(page), _crumbs(page)) != before, f"clicking {text!r} did nothing"
    return text


def _on_galaxy(page):
    return urlparse(page.url).path == "/galaxy"


def _label_of(tooltip):
    """The breadcrumb label a tooltip's choice leads to (see
    test_web_browser_maps.py)."""
    return re.split(r"[:,] ", tooltip, maxsplit=1)[0].strip()


GENERATED_CHOICE = r"(^|, )([1-9][\d,.]*\S* (of .+ )?sectors )?charted$"
"""A tooltip for a choice holding a generated sector (the fixture's
sectors, near the core)."""


PRIME_WINDOWS = [(243, 0, 0, 0), (81, 0, 0, 0), (27, 0, 0, 0), (9, 0, 0, 0), (3, 0, 0, 0), (1, 1, 0, 1)]
"""The cells of the windows on the way down from the whole galaxy to Fixture Prime (shell 1, band 0, wedge 1), as the
wheel carries a zoom toward it (the fixture's sectors lie near the core): [243, 81, 27, 9, 3, a sector]."""


def _descend(page):
    """Opens the window one size down on the way to Fixture Prime (from the galaxy the first is the 81-block window),
    as a wheel run toward it hands off to the next
    window (the sector blocks that a click took down are gone from the Galaxy Map, MAP.297). Returns the window's
    size, or None when the last window is open already."""
    win = _stage(page)["win"]
    size = win["centre"]["m"] if win else None
    sizes = [cell[0] for cell in PRIME_WINDOWS]
    if size == sizes[-1]:
        return None
    m, shell, band, wedge = PRIME_WINDOWS[1 if size is None else sizes.index(size) + 1]  # the galaxy's first window is 81 blocks
    page.evaluate("([cell]) => document.querySelector('#galaxymap3d-canvas').galaxyOpenWindow(cell)",
                  [{"m": m, "shell": shell, "band": band, "wedge": wedge}])
    _settle(page)
    return m


def _walk_down(page, stages=24):
    """Opens the windows down from the galaxy to the one round Fixture Prime, up to `stages` of them; returns the
    stages seen, [(query, crumbs)]."""
    steps = [(_stage_query(page), _crumbs(page))]
    for _ in range(stages):
        if _descend(page) is None:
            break
        steps.append((_stage_query(page), _crumbs(page)))
    return steps


def test_galaxy_map_hover_names_an_object_and_leaving_hides_it(page, map_site):
    _open_galaxy(page, map_site)
    found = _hover_choice(page)
    assert found, "hovering the whole galaxy shows no object"
    box = page.locator(GALAXY_CANVAS).bounding_box()
    page.mouse.move(box["x"] + box["width"] + 200, box["y"] + 10)
    page.wait_for_timeout(100)
    assert not page.locator("#galaxymap3d-tooltip").is_visible(), "leaving the map hides the tooltip"


def test_galaxy_map_window_stage_zooms_in_and_back_out(page, map_site):
    """MAP.151: a window stage (?win=) shows `g` cells around a cell; a click
    opens the window one size down, centred on the cell clicked, and Up comes
    back out, the breadcrumb following. (The click on a block is gone, MAP.297: the wheel and the crumbs go down.)"""
    _open_window(page, map_site, 243, 0, 0, 0)
    assert _stage_query(page) == "win=243.0.0.0"
    assert _crumbs(page) == ["Galaxy", "Block 0·0·0 (243)"], _crumbs(page)
    _descend(page)
    assert re.fullmatch(r"win=81\.\d+\.-?\d+\.\d+", _stage_query(page)), _stage_query(page)
    assert len(_crumbs(page)) == 3 and _crumbs(page)[-1].endswith("(81)"), _crumbs(page)
    page.locator(GALAXY_CANVAS).focus()
    page.keyboard.press("Escape")
    _settle(page)
    assert re.fullmatch(r"win=243\.\d+\.-?\d+\.\d+", _stage_query(page)), _stage_query(page)
    assert len(_crumbs(page)) == 2
    page.keyboard.press("Escape")
    _settle(page)
    assert _crumbs(page) == ["Galaxy"] and _query(page) == ""


def _as_admin(page, base):
    """Signs the page in as the fixture's admin: any session cookie does (`FixtureApi.auth_me`)."""
    from planetgen.api.authz import SESSION_COOKIE_NAME
    page.context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": "fixture", "url": base}])


def test_the_block_menu_offers_generate_to_an_admin_and_starts_the_blocks_job(page, map_site):
    """UX.127: the info card of a block carries "Generate this block" for an admin, which posts the Generate page's block mode
    with the block's key (the job `planetgen sectors --block` runs) and lands on the page that shows its progress; a visitor
    who is not signed in has no such button."""
    _open_window(page, map_site, 243, 0, 0, 0)
    assert page.locator("#galaxymap3d-info .starmap-generate").count() == 0, "a visitor is offered no Generate"
    _as_admin(page, map_site)
    _open_window(page, map_site, 243, 0, 0, 0)
    _expand_card(page)    # the Generate button is behind Expand
    page.wait_for_selector("#galaxymap3d-info .starmap-generate-block")
    posted = []

    def start(route):
        posted.append(route.request.post_data)
        route.fulfill(status=200, content_type="text/html", body="<!doctype html><title>Generate</title><p>Job started</p>")

    page.route("**/admin/generate", start)
    page.once("dialog", lambda dialog: dialog.accept())    # a block this size can run for hours: it asks first
    button = page.locator("#galaxymap3d-info .starmap-generate-block button", has_text="Generate this block")
    assert button.is_enabled()
    button.click()
    page.wait_for_function("() => document.body.textContent.includes('Job started')")
    (body,) = posted
    fields = dict(pair.split("=", 1) for pair in body.split("&"))
    assert fields["action"] == "galaxy" and fields["mode"] == "block" and fields["block_key"] == "243.0.0.0", fields
    assert fields["whole_block"] == "1" and fields["csrf_token"], fields


def test_the_block_generate_is_disabled_with_its_reason_while_a_job_runs(page, map_site, monkeypatch):
    """UX.127: while another run holds the Generate page, the block's button is disabled and says which job and links to it."""
    from planetgen.web import jobs
    monkeypatch.setattr(jobs, "active_job", lambda root=None: {"title": "Generate sectors around Fixture Prime"})
    _as_admin(page, map_site)
    _open_window(page, map_site, 243, 0, 0, 0)
    _expand_card(page)
    page.wait_for_selector("#galaxymap3d-info .starmap-generate-block")
    button = page.locator("#galaxymap3d-info .starmap-generate-block button", has_text="Generate this block")
    assert button.is_disabled() and button.get_attribute("title") == "Another job is running"
    hint = page.locator("#galaxymap3d-info .starmap-generate-block .hint")
    assert "Generate sectors around Fixture Prime is still running" in hint.inner_text()
    assert hint.locator("a").get_attribute("href").endswith("/admin/generate#current-job")


SAFE_EXTENT_PC = 10000
"""float: The farthest, in parsecs, a layer's buffer values may be from the object's own origin (`static/camerarelative.js`'s
SAFE_EXTENT_PC): float32 holds 10 kpc to 0.001 pc, a tenth of a pixel at the closest zoom of a layer."""

BUFFER_EXTENTS = "() => document.querySelector('#galaxymap3d-canvas').galaxyBufferExtents()"


@pytest.mark.parametrize("scheme", ("dark", "light"))
def test_a_superbubble_has_a_crisp_rim_brighter_than_its_fill_and_a_fill_fainter_than_the_old_circles(page, far_site, scheme):
    """MAP.236: a superbubble is its sphere, a faint fill (under the 0.14 of the three circles it was) and a rim where the sphere meets
    the line of sight that stands out from the background more than the fill does, in the dark and the light theme."""
    page.emulate_media(color_scheme=scheme)
    _open_galaxy(page, far_site, "?show=superbubbles")
    x, y, z, radius = FAR_BUBBLE
    page.evaluate("([c, r]) => document.querySelector('#galaxymap3d-canvas').galaxyGoTo({center: c, radius: r})",
                  [[x, y, z], radius * 1.25])
    _settle(page)
    page.wait_for_timeout(500)
    look = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyBubbleLook()")
    assert look and look["count"] == 1, look
    assert look["fillAlpha"] < 0.14, look
    assert look["rimAlpha"] > 3 * look["fillAlpha"], look
    spot = _screen_of(page, (x, y, z))
    assert spot
    box = page.locator(GALAXY_CANVAS).bounding_box()
    image = Image.open(io.BytesIO(_shot(page, GALAXY_CANVAS))).convert("RGB")
    row = int(spot[1] - box["y"])
    middle = int(spot[0] - box["x"])
    luminance = [sum(image.getpixel((column, row))) for column in range(image.width)]
    background = sum(luminance[2:12]) / 10
    shown = [column for column, value in enumerate(luminance) if abs(value - background) > 3 and 12 < column < image.width - 12]
    assert shown, "no bubble on the row through its middle"
    left, right = min(shown), max(shown)
    assert right - left > 100, (left, right)
    # Halfway out from the middle on the left the fill alone is drawn (the view's cross sits at the middle); at the edge the rim is.
    fill = sum(luminance[(left + middle) // 2 - 3:(left + middle) // 2 + 4]) / 7
    rim = max(luminance[left:left + 6], key=lambda value: abs(value - background))
    assert abs(fill - background) > 3, (fill, background)
    assert abs(rim - background) > abs(fill - background) + 6, (rim, fill, background)


def test_layers_far_from_the_centre_are_built_relative_to_origins_of_their_own(page, far_site):
    """MAP.240: a superbubble 250 kpc out, a polity's dots and capital there, and the sectors a job made 150 kpc out keep every buffer
    value within 10 kpc of the origin of the object that holds it (the object carries the origin as its position, in double
    precision), where a float32 of the absolute place is 0.016 pc off at 250 kpc."""
    _open_galaxy(page, far_site, "?made=0,9999999999&show=superbubbles")
    page.evaluate("() => document.querySelector(\"[data-action='territories']\").click()")   # it sits in the closed Menu
    page.wait_for_function(f"() => ({BUFFER_EXTENTS})().some((e) => e.layer === 'territory dots')", timeout=20000)
    extents = page.evaluate(BUFFER_EXTENTS)
    layers = {e["layer"] for e in extents}
    assert {"superbubbles", "made rings", "territory dots"} <= layers, layers
    for e in extents:
        assert e["extentPc"] <= SAFE_EXTENT_PC, e
    far = {e["layer"] for e in extents if max(abs(v) for v in e["origin"]) > 100000}
    assert {"superbubbles", "territory dots"} <= far or "made rings" in far, (far, extents)
    assert {"superbubbles", "territory dots"} <= far, "the bubble and the dots are really out there"
    assert sum(1 for e in extents if e["layer"] == "territory dots") >= 2, "dots spread over several kpc are split into cubes"


def test_the_blocks_and_their_outlines_are_built_relative_to_origins_of_their_own(page, map_site):
    """MAP.240: the block prisms of a window, the outline and shade of the block under the pointer and the frame of an opened
    sector hold offsets from the middle of what they draw, not its place in the galaxy."""
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    extents = page.evaluate(BUFFER_EXTENTS)
    prisms = [e for e in extents if e["layer"] == "block prisms"]
    assert prisms, extents
    assert all(e["extentPc"] <= SAFE_EXTENT_PC for e in extents), extents
    assert any(max(abs(v) for v in e["origin"]) > 50 for e in prisms), "the prisms sit away from the galaxy's middle, on an origin"


OPEN_PRIME = f"?sector={PRIME_CODE}&open=1"
"""The fixture's Fixture Prime (shell 1, band 0, wedge 1) opened in place."""


@pytest.mark.names_hidden
def test_galaxy_map_opens_a_sector_in_place(page, map_site):
    """MAP.66: a generated sector opens in the Galaxy Map itself, drawn by
    the Sector Map's scene: a star can be hovered and picked there (its
    panel, tooltip and NAV links as on the sector page), Up leaves the
    sector with it still selected, and Back and a reload come back to it."""
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    assert page.locator("#galaxymap3d-info h3").inner_text() == "Fixture Prime"
    page.wait_for_timeout(500)
    tip = page.locator("#galaxymap3d-tooltip")
    info = page.locator("#galaxymap3d-info")
    picked = None
    for x, y in _bright_spots(page, GALAXY_CANVAS):
        page.mouse.move(x, y)
        page.wait_for_timeout(120)
        hovered = tip.inner_text() if tip.is_visible() else ""
        page.mouse.click(x, y)
        heading = info.locator("h3").inner_text() if info.locator("h3").count() else None
        if heading in SYSTEM_NAMES:
            picked = (heading, hovered)
            break
    assert picked, "no click on the opened sector picked a star"
    name, hovered = picked
    assert name in hovered, f"the tooltip for {name}: {hovered!r}"
    system_id = next(i for i, n, *_rest in SYSTEMS if n == name)
    assert info.locator("a", has_text="View system").get_attribute("href").endswith(f"/system/{system_id}")
    nav = [b.inner_text() for b in info.locator("button").all()]
    assert "Start Here" in nav and "End Here" in nav

    assert _on_galaxy(page) and _query(page) == OPEN_PRIME[1:]
    page.click('#galaxymap3d-controls [data-action="up"]')
    _settle(page)
    assert _query(page) == f"sector={PRIME_CODE}", _query(page)
    page.go_back()
    _settle(page)
    assert _query(page) == OPEN_PRIME[1:]
    page.reload()
    page.wait_for_selector("#galaxymap3d-steps [data-steps-panel] li", state="attached")
    _settle(page)
    assert _query(page) == OPEN_PRIME[1:]
    assert any(_bright_spots(page, GALAXY_CANVAS)), "the reloaded sector draws its stars"


FAR_CODE = provisional_sector_designation(40, 0, 0)
"""A cell far outside every generated sector's 9-block."""


def _window_size(page):
    win = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyStage()")["win"]
    return win["centre"]["m"] if win else None


@pytest.mark.names_hidden
def test_galaxy_map_infobox_names_the_region_and_flies_to_it(page, map_site):
    """MAP.315, MAP.349: a pick inside a charted region says "Region: <name>" and "In: <family>", each line with its own Go to that flies
    the camera to that block (the region's 9-block, then its 81-block family), and a pick in no charted region has no such lines."""
    _open_galaxy(page, map_site, OPEN_PRIME)
    info = page.locator("#galaxymap3d-info")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    assert info.locator("h3").inner_text() == "Fixture Prime"
    line = info.locator(".map-info-region")
    _expand_card(page)    # the region lines are behind Expand
    line.wait_for(state="visible", timeout=20000)
    places = info.locator(".map-info-region .map-info-place")
    assert [t.replace("\nGo to", " Go to") for t in places.all_inner_texts()] == ["Region: Fixture Beez Region Go to", "In: Fixture Reach Go to"], places.all_inner_texts()
    # The family first: its 81-block, then the region itself from the same card, framed by its extent (MAP.319).
    places.filter(has_text="Fixture Reach").locator(".map-info-go").click()
    _settle(page)
    assert _window_size(page) == 81, _query(page)
    _open_galaxy(page, map_site, OPEN_PRIME)
    _expand_card(page)
    info.locator(".map-info-place", has_text="Fixture Beez Region").locator(".map-info-go").click()
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=30000)
    _settle(page)
    from planetgen.galaxy.spherical import sector_position_pc
    from tests.map_site_support import REGION_RADIUS_PC
    after = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
    centre = sector_position_pc(*SECTORS[0][1:4], 4.0)
    assert max(abs(a - b) for a, b in zip(after["target"], centre)) < 0.05, (after["target"], centre)
    size = page.evaluate("() => { const c = document.querySelector('#galaxymap3d-canvas'); return [c.clientWidth, c.clientHeight]; }")
    tan = math.tan(math.radians(page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyFov()")) / 2)
    half = math.atan(min(tan, tan * size[0] / size[1]))
    assert math.asin(min(1.0, REGION_RADIUS_PC / after["dist"])) <= half + 1e-6, (after["dist"], REGION_RADIUS_PC)
    assert after["dist"] == pytest.approx(REGION_RADIUS_PC / math.sin(half) * 1.08, rel=0.05) or after["dist"] > REGION_RADIUS_PC / math.sin(half)
    assert page.evaluate(LINES)["chartedLines"] == 0, "no sector box is drawn once the flight to the region is over"
    # Where nothing is charted the block has no name on the wire and the card has no line.
    _open_galaxy(page, map_site, f"?sector={FAR_CODE}&open=1")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    page.wait_for_timeout(1500)
    assert not info.locator(".map-info-region").is_visible()


OPEN_UNCHARTED = f"?sector={UNCHARTED_CODE}&open=1"
"""The fixture's ungenerated shell 1, band 0, wedge 0, where the scatters left objects, opened in place."""


@pytest.mark.viewport(1280, 900)
def test_galaxy_map_opens_an_uncharted_sector_with_its_scattered_objects_and_marks_it(page, map_site):
    """MAP.162: a sector nothing was generated in opens in the Galaxy Map from what the scatters left in it,
    carrying an "uncharted" mark in its crumb, its info panel and a badge in the map's frame; its objects
    pick (a scattered phenomenon with its uncharted page, none with course buttons) and Up leaves it selected without the mark."""
    _open_galaxy(page, map_site, OPEN_UNCHARTED)
    page.wait_for_selector(".sector-uncharted-badge", state="attached")
    badge = page.locator(".sector-uncharted-badge")
    assert "Uncharted sector" in badge.inner_text() and "3 objects" in badge.inner_text()
    assert _crumbs(page)[-1].endswith("(uncharted)"), _crumbs(page)
    info = page.locator("#galaxymap3d-info")
    assert info.locator("h3").inner_text() == "Uncharted sector"
    assert "Uncharted: not charted yet" in info.inner_text()
    assert _query(page) == OPEN_UNCHARTED[1:]
    page.wait_for_timeout(500)
    from planetgen.galaxy.spherical import sector_position_pc
    cx, cy, cz = sector_position_pc(*UNCHARTED_CELL, 4.0)
    picked = None
    # Its bright star, black hole and nebula, each by where it sits.
    for point in ((cx + 0.5, cy + 0.4, cz), (cx - 0.6, cy, cz + 0.3), (cx, cy - 0.7, cz - 0.2)):
        spot = _screen_of(page, point)
        assert spot, point
        page.mouse.click(*spot)
        page.wait_for_timeout(150)
        heading = info.locator("h3").inner_text() if info.locator("h3").count() else None
        if heading and "(uncharted)" in heading:
            picked = heading
            break
    assert picked, "no click on the uncharted sector picked one of its objects"
    assert "its sector is not charted yet" in info.inner_text()
    # MAP.211: a scattered phenomenon has an uncharted page to view; a bright star waiting for its sector has none.
    waiting_star = picked.endswith(" star (uncharted)") and not picked.startswith("Hypervelocity")
    view = info.locator("a", has_text="View phenomenon")
    if waiting_star:
        assert view.count() == 0, "a waiting bright star has no page to view"
    else:
        assert view.count() == 1
        assert "/uncharted-phenomenon/" in view.get_attribute("href")
    assert info.locator("button", has_text="Start Here").count() == 0
    page.click('#galaxymap3d-controls [data-action="up"]')
    _settle(page)
    assert _query(page) == f"sector={UNCHARTED_CODE}", _query(page)
    assert page.locator(".sector-uncharted-badge").count() == 0
    assert not _crumbs(page)[-1].endswith("(uncharted)")


@pytest.mark.viewport(1280, 900)
def test_a_binary_picks_as_one_system_in_an_opened_sector(page, map_site):
    """MAP.136: every dot of an opened sector, a binary's companion too,
    picks a system by its own name; a star is never picked on its own. The names are off: a click on a name flies
    to its item (MAP.298), and a flight would move the spots this test tries."""
    _open_galaxy(page, map_site, OPEN_PRIME + "&hide=label-system,label-star,label-object")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    page.wait_for_timeout(500)
    info = page.locator("#galaxymap3d-info")
    headings = set()
    # Only dots whose tooltip names a system are clicked: a click on the
    # opened cell's outline, among the brightest spots, would open another sector.
    # The 400 spots are tried with pointer events sent from the page (a real pointer move to each cost
    # 0.8 s under load, 340 s in all, TEST.139 stage 4); the few that name a system are then clicked for real.
    for x, y, text in page.evaluate(NAMED_SPOTS, [_bright_spots(page, GALAXY_CANVAS, limit=400), sorted(SYSTEM_NAMES)]):
        page.mouse.move(x, y)
        page.mouse.click(x, y)
        page.wait_for_timeout(100)
        if info.locator("h3").count():
            headings.add(info.locator("h3").inner_text())
    assert "Twin Lamps" in headings, headings
    assert headings <= SYSTEM_NAMES | {"Fixture Prime"}, headings


def _camera(page):
    page.wait_for_function("() => { const c = document.querySelector('#galaxymap3d-canvas').galaxyCamera(); return c && !c.flying; }")
    return page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")


@pytest.mark.names_hidden
def test_c_recenters_the_view_on_the_selection_and_the_mouse_turns_about_it(page, map_site):
    """MAP.138: C (or Menu's Center on selection) moves the center the view
    turns and zooms about to the selected star, keeping the angle and zoom;
    a mouse drag then turns about it and leaves it where it is."""
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    page.wait_for_timeout(500)
    _x, _y, _name = _click_a_star(page)
    before = _camera(page)
    page.locator(GALAXY_CANVAS).focus()
    page.keyboard.press("c")
    after = _camera(page)
    assert after["target"] != before["target"], "the center moved to the star"
    assert abs(after["tilt"] - before["tilt"]) < 1e-6 and abs(after["zoom"] - before["zoom"]) < 1e-9
    box = page.locator(GALAXY_CANVAS).bounding_box()
    cx, cy = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
    page.mouse.move(cx, cy)
    page.mouse.down()
    page.mouse.move(cx + 80, cy + 40, steps=5)
    page.mouse.up()
    turned = _camera(page)
    assert turned["tilt"] != after["tilt"] and turned["target"] == after["target"], (after, turned)
    # With nothing selected, C goes back to the middle of the view.
    page.keyboard.press("Escape")
    page.keyboard.press("c")
    again = _camera(page)
    assert again["target"] != after["target"]


def _on_line_of_sight(camera, point, tolerance=1e-3):
    """Whether `point` lies on the line from the camera's eye through its center (the cross)."""
    ray = [t - e for t, e in zip(camera["target"], camera["eye"])]
    toward = [p - e for p, e in zip(point, camera["eye"])]
    cross = [ray[1] * toward[2] - ray[2] * toward[1], ray[2] * toward[0] - ray[0] * toward[2],
             ray[0] * toward[1] - ray[1] * toward[0]]
    return math.sqrt(sum(c * c for c in cross)) <= tolerance * math.sqrt(sum(r * r for r in ray)) * math.sqrt(sum(t * t for t in toward))


def test_wheel_zoom_turns_about_the_centre_and_a_shift_click_chooses_it(page, map_site):
    """MAP.199: the mouse turns the view about the cross in the middle of the
    map, and a wheel with the pointer on the cross zooms about it; a
    Shift-click makes the spot clicked the new center, keeping the angle and
    zoom, and the cross stays in the middle. (MAP.225: with the pointer off
    the cross the wheel zooms toward what is under it instead.)"""
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    page.wait_for_timeout(500)
    canvas = page.locator(GALAXY_CANVAS)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    page.wait_for_timeout(200)
    box = canvas.bounding_box()
    cx, cy = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
    cross = page.locator("#galaxymap3d-centre")
    assert cross.is_visible()
    cross_box = cross.bounding_box()
    assert abs(cross_box["x"] + cross_box["width"] / 2 - cx) < 2 and abs(cross_box["y"] + cross_box["height"] / 2 - cy) < 2
    before = _camera(page)
    page.mouse.move(cx, cy)
    page.mouse.wheel(0, -300)
    page.wait_for_timeout(300)
    zoomed = _camera(page)
    assert zoomed["zoom"] != before["zoom"], "the wheel zoomed"
    assert _on_line_of_sight(zoomed, before["target"]), "the pointer was on the cross, the view stayed on its line of sight"
    page.mouse.wheel(0, 300)
    page.wait_for_timeout(300)
    assert _on_line_of_sight(_camera(page), before["target"])
    pre = _camera(page)
    spot_x, spot_y = cx + box["width"] / 4, cy - box["height"] / 6
    page.mouse.move(spot_x, spot_y)
    page.keyboard.down("Shift")
    page.mouse.click(spot_x, spot_y)
    page.keyboard.up("Shift")
    moved = _camera(page)
    assert moved["target"] != before["target"], "the spot clicked became the center"
    assert moved["zoom"] == pytest.approx(pre["zoom"]) and abs(moved["tilt"] - pre["tilt"]) < 1e-6
    page.mouse.move(cx, cy)
    page.mouse.wheel(0, -300)
    page.wait_for_timeout(300)
    assert _on_line_of_sight(_camera(page), moved["target"])


def test_wheel_in_zooms_toward_the_star_under_the_pointer(page, map_site):
    """MAP.225: with the pointer off the cross the wheel zooms toward the point
    under it: the center moves that way, the star under the pointer stays
    under the pointer, and a wheel out brings the view back out about the
    same point."""
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    page.wait_for_timeout(500)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    page.wait_for_timeout(200)
    box = page.locator(GALAXY_CANVAS).bounding_box()
    cx, cy = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
    spots = [(x, y) for x, y in _bright_spots(page, GALAXY_CANVAS, white=True)
             if math.hypot(x - cx, y - cy) > box["width"] / 8 and 40 < y < VIEWPORT["height"] - 40]
    assert spots, "the Galaxy Map shows a star away from its middle"
    x, y = spots[0]
    before = _camera(page)
    page.mouse.move(x, y)
    page.mouse.wheel(0, -200)
    page.wait_for_timeout(400)
    zoomed = _camera(page)
    assert zoomed["zoom"] < before["zoom"], "the wheel zoomed in"
    assert zoomed["target"] != pytest.approx(before["target"], abs=1e-9), "the center moved toward the pointer"
    near = [(sx, sy) for sx, sy in _bright_spots(page, GALAXY_CANVAS, white=True) if math.hypot(sx - x, sy - y) <= 6]
    assert near, "the star under the pointer is still under it"


def test_the_header_and_every_crumb_follow_the_cell_the_camera_is_in(page, map_site):
    """MAP.227: the breadcrumb and the line under it come from where the view's middle is, not from the stage the
    map was opened at: moving the center into the next sector names that sector, and a crumb above it flies out."""
    from planetgen.galaxy.spherical import sector_position_pc
    from tests.map_site_support import EDGE_PC

    _open_galaxy(page, map_site, f"?sector={PRIME_CODE}")
    where = page.locator("#galaxymap3d-where")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-where') && !document.querySelector('#galaxymap3d-where').hidden")
    first = where.inner_text()
    crumbs = _crumbs(page)
    assert first.startswith("In Sector ") and crumbs[0] == "Galaxy" and crumbs[-1].startswith("Sector "), (first, crumbs)

    point = sector_position_pc(1, 0, 2, EDGE_PC)
    page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyCentreOn(...p)", list(point))
    page.wait_for_function("(was) => document.querySelector('#galaxymap3d-where').textContent !== was", arg=first)
    _settle(page)
    moved = where.inner_text()
    assert moved.startswith("In Sector ") and moved != first
    after = _crumbs(page)
    assert after[:-1] == crumbs[:-1] or len(after) == len(crumbs), (crumbs, after)
    assert after[-1] != crumbs[-1] and after[-1].split(" (")[0] == "Sector " + moved[len("In Sector "):], (moved, after)

    # The first crumb flies out to the whole galaxy.
    page.locator("#galaxymap3d-steps summary").click()
    page.locator("#galaxymap3d-steps [data-steps-panel] li", has_text="Galaxy").first.click()
    page.wait_for_function("() => document.querySelectorAll('#galaxymap3d-steps [data-steps-panel] li').length <= 2")
    _settle(page)
    assert _crumbs(page) == ["Galaxy"], _crumbs(page)


def test_on_map_zoom_slider_and_direction_star(page, map_site):
    """MAP.200: a zoom slider and a direction star on the map. The slider and
    the wheel move together; the arrows move the center (the readout and the
    cross show where it is), Alt and the arrow keys do too, and the middle
    button puts the center back at the middle of the region viewed."""
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    page.wait_for_timeout(500)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    start = _camera(page)
    slider = page.locator("#galaxymap3d-zoom-slider")
    assert slider.is_visible() and page.locator(".starmap-dirstar button").count() == 5
    value = int(slider.input_value())
    slider.fill(str(min(1000, value + 300)))
    zoomed = _camera(page)
    assert zoomed["zoom"] < start["zoom"], "the slider zooms in"
    assert zoomed["target"] == pytest.approx(start["target"], abs=1e-6), "zooming keeps the center"
    page.locator(GALAXY_CANVAS).hover()
    page.mouse.wheel(0, 200)
    page.wait_for_timeout(300)
    assert int(slider.input_value()) < value + 300, "the wheel moves the slider"
    for name in ("right", "up"):
        before = _camera(page)
        page.click(f'.starmap-dir-{name}')
        moved = _camera(page)
        assert moved["target"] != before["target"], f"the {name} arrow moves the center"
        assert abs(moved["tilt"] - before["tilt"]) < 1e-6
    assert page.locator(".starmap-centre-readout").is_visible()
    assert "is-moving" in page.locator("#galaxymap3d-centre").get_attribute("class")
    page.locator(GALAXY_CANVAS).focus()
    before = _camera(page)
    page.keyboard.press("Alt+ArrowLeft")
    assert _camera(page)["target"] != before["target"], "Alt and an arrow move the center"
    page.click(".starmap-dir-reset")
    reset = _camera(page)
    assert reset["target"] == pytest.approx(start["target"], abs=1e-6), "the middle button returns to the region's middle"


def test_galaxy_map_back_forward_and_url_state(page, map_site):
    _open_galaxy(page, map_site)
    steps = _walk_down(page, stages=3)
    queries = [q for q, _ in steps]
    page.go_back()
    _settle(page)
    assert _stage_query(page) in queries[:-1]
    at = queries.index(_stage_query(page))
    assert _crumbs(page) == steps[at][1]
    page.go_forward()
    _settle(page)
    assert _stage_query(page) in queries[at + 1:]

    # The map's own Back, Forward and Up.
    here, crumbs = _stage_query(page), _crumbs(page)
    page.click('#galaxymap3d-controls [data-action="back"]')
    _settle(page)
    assert _stage_query(page) != here
    page.click('#galaxymap3d-controls [data-action="forward"]')
    _settle(page)
    assert _stage_query(page) == here
    page.click('#galaxymap3d-controls [data-action="up"]')
    _settle(page)
    assert _crumbs(page) == crumbs[:-1]

    # An address opens the view it was copied at directly (MAP.228); Whole galaxy goes home.
    page.wait_for_timeout(600)
    deep_url, deep_crumbs = page.url, _crumbs(page)
    assert "cam=" in deep_url, deep_url
    _open(page, deep_url, GALAXY_CANVAS)
    page.wait_for_selector("#galaxymap3d-steps [data-steps-panel] li", state="attached")
    _settle(page)
    assert _crumbs(page) == deep_crumbs
    page.click('#galaxymap3d-controls [data-action="reset"]')
    _settle(page)
    assert _crumbs(page) == ["Galaxy"] and _stage_query(page) == ""


def test_galaxy_map_keys_go_up_and_home(page, map_site):
    """Escape goes up a window and Home to the whole galaxy; the arrows no longer mark a block (MAP.297: no sector
    boxes to move among on the Galaxy Map)."""
    _open_galaxy(page, map_site)
    canvas = page.locator(GALAXY_CANVAS)
    canvas.focus()
    page.keyboard.press("ArrowRight")
    assert not page.locator("#galaxymap3d-tooltip").is_visible(), "an arrow key marks no choice"
    page.keyboard.press("Enter")
    _settle(page)
    assert _crumbs(page) == ["Galaxy"], "Enter with nothing selected goes nowhere"
    _descend(page)
    _descend(page)
    one_down = _crumbs(page)
    assert len(one_down) > 2
    page.keyboard.press("Escape")
    _settle(page)
    assert len(_crumbs(page)) == len(one_down) - 1
    page.keyboard.press("Home")
    _settle(page)
    assert _crumbs(page) == ["Galaxy"]


def _galaxy_scale(page):
    """The scale line: its label and its bar's width (a zoom between two
    nice lengths only changes the bar)."""
    return page.locator("#galaxymap3d-scale").inner_html()


def test_galaxy_map_scale_is_one_line(page, map_site):
    """MAP.60: just the bar and its length, side by side."""
    _open_galaxy(page, map_site)
    scale = page.locator("#galaxymap3d-scale")
    assert scale.locator(".starmap-scale-bar").count() == 1
    labels = scale.locator(".starmap-scale-label")
    assert labels.count() == 1
    assert re.fullmatch(r"[\d.,]+ sectors? · .+", labels.inner_text()), labels.inner_text()
    assert "px" not in scale.inner_text() and "block" not in scale.inner_text()
    bar, label = scale.locator(".starmap-scale-bar").bounding_box(), labels.bounding_box()
    assert abs((bar["y"] + bar["height"] / 2) - (label["y"] + label["height"] / 2)) < label["height"], "one line"


def test_galaxy_map_scale_line_follows_the_zoom_on_a_free_stage(page, map_site):
    _open_galaxy(page, map_site)
    whole = _galaxy_scale(page)
    assert whole, "the whole galaxy has a scale line"
    reset_view = page.locator('#galaxymap3d-controls [data-action="reset-view"]')
    _descend(page)
    assert not reset_view.is_disabled(), f"no free stage by {_crumbs(page)}"
    stage_scale = _galaxy_scale(page)
    assert stage_scale != whole, "a closer stage has a shorter scale"
    box = page.locator(GALAXY_CANVAS).bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    page.mouse.wheel(0, -400)
    _settle(page)
    assert _galaxy_scale(page) != stage_scale, "the wheel zooms and the scale follows"
    page.click('#galaxymap3d-menu > summary')
    reset_view.click()
    _settle(page)
    assert _galaxy_scale(page) == stage_scale


def _bookmarks(page):
    return page.evaluate("""() => {
        const key = Object.keys(localStorage).find(k => k.startsWith('planetgen.bookmarks.'));
        return key ? JSON.parse(localStorage.getItem(key)) : [];
    }""")


def test_galaxy_map_bookmark_star_saves_the_stage_and_the_menu_opens_it(page, map_site):
    _open_galaxy(page, map_site)
    _descend(page)
    stage_query = _query(page)
    star = page.locator("#galaxymap3d-crumbs .galaxy-bookmark")
    assert star.count() == 1
    star.click()
    saved = _bookmarks(page)
    assert len(saved) == 1 and saved[0]["kind"] == "stage" and "cam=" in saved[0]["value"], saved
    assert star.get_attribute("aria-pressed") == "true"

    page.click('#galaxymap3d-controls [data-action="reset"]')
    _settle(page)
    menu = page.locator("[data-bookmarks-menu]")
    menu.locator("summary").click()
    links = menu.locator("[data-bookmarks-panel] a")
    assert links.count() == 1
    assert menu.locator("kbd").first.inner_text() == "1"
    links.first.click()
    page.wait_for_load_state("load")
    _settle(page)
    assert _query(page) == stage_query


def test_galaxy_map_bookmark_keys_are_1_to_9_while_the_map_has_focus(page, map_site):
    """MAP.81: plain 1 to 9 (Ctrl+1 to 9 are the browser's tab keys),
    only while the map has focus."""
    _open_galaxy(page, map_site)
    _descend(page)
    stage_query = _query(page)
    page.locator("#galaxymap3d-crumbs .galaxy-bookmark").click()
    page.click('#galaxymap3d-controls [data-action="reset"]')
    _settle(page)
    page.locator(GALAXY_CANVAS).focus()
    for key in ("Control+1", "2"):
        page.keyboard.press(key)
        _settle(page)
        assert _crumbs(page) == ["Galaxy"], f"{key} does nothing"
    page.locator("h1").first.evaluate("h => { h.tabIndex = -1; h.focus(); }")
    page.keyboard.press("1")
    _settle(page)
    assert _crumbs(page) == ["Galaxy"], "1 does nothing when the map doesn't have focus"
    page.locator("#galaxymap3d-address-input").evaluate("i => { i.closest('form').hidden = false; i.focus(); }")
    page.keyboard.press("1")
    _settle(page)
    assert _crumbs(page) == ["Galaxy"], "typing 1 in the address box is just typing"
    page.locator(GALAXY_CANVAS).focus()
    with page.expect_navigation():
        page.keyboard.press("1")
    _settle(page)
    assert _query(page) == stage_query


def test_sector_page_bookmark_shows_in_the_galaxy_menu(page, map_site):
    _open_sector(page, map_site)
    toggle = page.locator("[data-bookmark-toggle]")
    toggle.click()
    saved = _bookmarks(page)
    assert len(saved) == 1 and saved[0]["kind"] == "sector" and saved[0]["name"] == "Fixture Prime"
    _open_galaxy(page, map_site)
    menu = page.locator("[data-bookmarks-menu]")
    menu.locator("summary").click()
    assert "Fixture Prime" in menu.locator("[data-bookmarks-panel]").inner_text()


# --- NAV.30: no link out of a course being picked -----------------------------------

def test_sector_map_pick_mode_offers_only_the_pick_button(page, map_site):
    _open_sector(page, map_site, "?pick=to")
    _x, _y, name = _click_a_star(page)
    info = page.locator("#galaxymap3d-info")
    buttons = [b.inner_text() for b in info.locator(".map-info-actions button:not(.map-info-bookmark):not(.map-info-close):not(.map-info-expand)").all()]
    assert info.locator("a", has_text="View").count() == 0
    assert buttons == ["End Here"], buttons
    page.locator(".starmap-sr-list button", has_text="Fixture Pulsar").evaluate("b => b.click()")
    assert info.locator("a", has_text="View").count() == 0
    assert [b.inner_text() for b in info.locator(".map-info-actions button:not(.map-info-bookmark):not(.map-info-close):not(.map-info-expand)").all()] == ["End Here"]


def test_start_here_keeps_the_view_and_end_here_opens_the_course(page, map_site):
    """NAV.29 and NAV.33: "Start Here" on a star keeps the user on the map at
    the same view, with the start in the URL and a banner asking for the
    destination; "End Here" on another object then opens the NAV page with
    both ends."""
    _open_sector(page, map_site)
    _x, _y, first = _click_a_star(page)
    start_id = next(i for i, n, *_rest in SYSTEMS if n == first)
    url_before = page.url
    _expand_card(page)    # Start Here is behind Expand
    page.locator("#galaxymap3d-info button", has_text="Start Here").click()
    assert page.url.split("?")[0] == url_before.split("?")[0], "no page change"
    assert parse_qs(urlparse(page.url).query) == {"pick": ["to"], "from": [f"system:{start_id}"]}, page.url
    banner = page.locator(".pick-banner")
    assert "Choosing a destination" in banner.inner_text() and first in banner.inner_text()
    other = next(n for n in SYSTEM_NAMES if n != first)
    page.locator(".starmap-sr-list button", has_text=other).first.evaluate("b => b.click()")
    assert _info_title(page) == other
    end_id = next(i for i, n, *_rest in SYSTEMS if n == other)
    link = page.locator("#galaxymap3d-info a", has_text="End Here")
    assert link.count() == 1 and page.locator("#galaxymap3d-info a", has_text="View").count() == 0
    query = parse_qs(urlparse(link.get_attribute("href")).query)
    assert query == {"from": [f"system:{start_id}"], "to": [f"system:{end_id}"]}, query


def test_cancel_ends_a_pick_begun_on_the_map(page, map_site):
    _open_sector(page, map_site)
    _click_a_star(page)
    _expand_card(page)
    page.locator("#galaxymap3d-info button", has_text="Start Here").click()
    page.locator(".pick-banner button", has_text="Cancel").click()
    assert page.locator(".pick-banner").count() == 0
    assert "pick" not in parse_qs(urlparse(page.url).query)


def _click_galaxy_stars(page, count=5):
    """Clicks the brightest fixture stars on the Galaxy Map, each on a freshly settled view, where the map
    draws them now (they are points of light a pixel or two wide until the camera is about a sector away,
    MAP.208, too small to find in a screenshot by their brightness); the info panel's heading after each click."""
    from tests.map_site_support import STAR_FIELD

    headings = []
    box = page.locator(GALAXY_CANVAS).bounding_box()
    for _ in range(count):
        # A black hole, neutron star or quasar does take the click (MAP.80): leave the stars beside one alone.
        points = [spot for spot in (_screen_of(page, (p["x"], p["y"], p["z"])) for p in POINT_FIELD) if spot]
        spot = None
        for star in sorted(STAR_FIELD, key=lambda star: -star["luminosity_sol"]):
            here = _screen_of(page, (star["x"], star["y"], star["z"]))
            if (here and box["x"] + 20 < here[0] < box["x"] + box["width"] - 20
                    and box["y"] + 20 < here[1] < box["y"] + box["height"] - 20
                    and all(abs(here[0] - x) + abs(here[1] - y) > 24 for x, y in points)):
                spot = here
                break
        assert spot, "no stars drawn to click"
        page.mouse.click(*spot)
        page.wait_for_timeout(500)
        _settle(page)
        heading = page.locator("#galaxymap3d-info h3")
        headings.append(heading.inner_text() if heading.count() else None)
    return headings


def test_sector_map_draws_its_nebula_from_the_shape_mesh(page, map_site):
    """MAP.103: the Sector Map's nebula (the fixture's "Fixture Veil") is
    drawn from its shape mesh once that arrives, in place of the sphere."""
    shape_requests = []
    page.on("request", lambda request: shape_requests.append(request.url) if "/shape" in request.url else None)
    _open_sector(page, map_site)
    page.locator("#galaxymap3d-canvas[data-sector-nebula-meshes='1']").wait_for(state="attached", timeout=15000)
    assert any("/galaxy/nebula/31/shape?lod=low" in url for url in shape_requests), shape_requests


def test_galaxy_map_draws_a_nebula_from_its_shape_mesh(page, map_site):
    """MAP.103: a nebula big enough on screen is drawn from its shape mesh
    (fetched from the site), not only as a sprite; far out it stays a sprite."""
    shape_requests = []
    page.on("request", lambda request: shape_requests.append(request.url) if "/shape" in request.url else None)
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    # A locator, not wait_for_function: the page's CSP forbids evaluating strings.
    page.locator("#galaxymap3d-canvas:not([data-nebula-meshes='0'])").wait_for(state="attached", timeout=15000)
    assert any("/galaxy/nebula/801/shape?lod=low" in url for url in shape_requests), shape_requests


def test_galaxy_map_shades_a_nebula_over_unfilled_sectors(page, map_site):
    """MAP.104: a nebula is drawn over space where no sector is generated
    (the fixture's generated sectors are all in the core, thousands of
    parsecs away), in a close view of the block as well as a far one."""
    _open_window(page, map_site, 3, 250, 0, 0, g=3)
    page.locator("#galaxymap3d-canvas:not([data-nebula-meshes='0'])").wait_for(state="attached", timeout=15000)   # MAP.264: more than one, nothing is clipped now


LABELS = "() => (document.querySelector('#galaxymap3d-canvas').__mapLabels.shown())"


def test_sector_map_names_its_systems_through_the_label_engine_with_the_system_names_switch(page, map_site):
    """MAP.231: a sector's own page names its stars through the label engine (DOM names, none over another, never
    picking), with no text sprites in the scene; the System names switch (the address's `hide`) takes them away."""
    _open_sector(page, map_site)
    page.wait_for_function(f"() => ({LABELS})().some((label) => label.kind === 'system')")
    shown = page.evaluate(LABELS)
    assert all(label["kind"] in ("system", "sector") for label in shown), shown
    rects = page.evaluate("() => [...document.querySelectorAll('.map-label')].map((e) => { const r = e.getBoundingClientRect(); return [r.left, r.top, r.right, r.bottom]; })")
    for i, a in enumerate(rects):
        for b in rects[i + 1:]:
            assert not (a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]), (a, b)
    _open_sector(page, map_site, "?hide=label-system")
    page.wait_for_timeout(600)
    assert not [label for label in page.evaluate(LABELS) if label["kind"] == "system"]


@pytest.mark.viewport(1000, 900)
def test_galaxy_map_names_a_nebula_where_it_is_large_enough_and_not_beyond(page, map_site):
    """MAP.268 (MAP.229, MAP.233): the label engine names the fixture nebula through the Labels group of the Menu:
    in the window it fills, off with the toggle (which the address keeps), and not from the whole galaxy, where it
    is far too small to read."""
    _open_galaxy(page, map_site)
    assert not [label for label in page.evaluate(LABELS) if label["kind"] == "nebula"], "far out nothing is named"
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_function(f"() => ({LABELS})().some((label) => label.kind === 'nebula' && label.text === 'Fixture Cloud')")
    labels = page.locator("#galaxymap3d-canvas").evaluate("c => c.parentElement.querySelectorAll('.map-label').length")
    assert labels > 0
    page.click("#galaxymap3d-menu > summary")
    toggle = page.locator('#galaxymap3d-kinds [data-label="nebula"]')
    assert toggle.inner_text() == "Nebula names"
    toggle.click()
    page.wait_for_function(f"() => !({LABELS})().some((label) => label.kind === 'nebula')")
    assert "label-nebula" in parse_qs(urlparse(page.url).query).get("hide", [""])[0]
    toggle.click()
    page.wait_for_function(f"() => ({LABELS})().some((label) => label.kind === 'nebula')")


def test_galaxy_map_names_the_stars_of_the_open_sector_and_the_toggle_takes_them_away(page, map_site):
    """MAP.268 (MAP.230): in a sector opened in place the named stars are labelled, brightest first, and the System
    names toggle removes them (and its address entry brings it back)."""
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    page.wait_for_function(f"() => ({LABELS})().some((label) => label.kind === 'system')")
    named = [label["text"] for label in page.evaluate(LABELS) if label["kind"] == "system"]
    assert named and set(named) <= SYSTEM_NAMES, named
    page.click("#galaxymap3d-menu > summary")
    page.click('#galaxymap3d-kinds [data-label="system"]')
    page.wait_for_function(f"() => !({LABELS})().some((label) => label.kind === 'system')")
    assert "label-system" in parse_qs(urlparse(page.url).query).get("hide", [""])[0]


def test_galaxy_map_names_a_generated_sector_in_its_window(page, map_site):
    """MAP.268 (MAP.230): in a window of sectors close enough to read, a generated sector shows its name."""
    _open_window(page, map_site, 1, 1, 0, 1, g=3)
    page.wait_for_function(f"() => ({LABELS})().some((label) => label.kind === 'sector' && label.text === 'Fixture Prime')")


@pytest.mark.viewport(1000, 900)
def test_galaxy_map_labels_never_take_the_pick_and_stay_inside_the_map(page, map_site):
    """MAP.229: the labels are decoration over the map: they do not catch the pointer, none overlap and every one sits
    inside the map's window."""
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_function(f"() => ({LABELS})().length > 0")
    boxes = page.evaluate("""() => {
        const root = document.querySelector('.map-labels');
        const frame = root.getBoundingClientRect();
        return {
            events: getComputedStyle(root).pointerEvents,
            inside: [...root.children].every((e) => { const r = e.getBoundingClientRect(); return r.left >= frame.left - 1 && r.right <= frame.right + 1 && r.top >= frame.top - 1 && r.bottom <= frame.bottom + 1; }),
            rects: [...root.children].map((e) => { const r = e.getBoundingClientRect(); return [r.left, r.top, r.right, r.bottom]; }),
        };
    }""")
    assert boxes["events"] == "none" and boxes["inside"]
    rects = boxes["rects"]
    for i, a in enumerate(rects):
        for b in rects[i + 1:]:
            assert not (a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]), (a, b)


# --- NAV.40: bookmarks keep a course pick ---------------------------------------------

SEED_BOOKMARKS = """(entries) => {
    const db = document.querySelector("[data-bookmark-db]").getAttribute("data-bookmark-db");
    localStorage.setItem("planetgen.bookmarks." + db, JSON.stringify(entries));
}"""

PICK_BOOKMARKS = [
    {"name": "Saved view", "kind": "stage", "value": "/galaxy?cam=3003.436762,5.294542628,-0.6473350564,151.006709,0.325344516,-0.3247703979,-0.6274082299,0.6285173408", "created": 1},
    {"name": "Middling Sun", "kind": "system", "value": "system:702", "url": "/system/702", "created": 2},
    {"name": "Fixture Second", "kind": "sector", "value": "8·0·3", "url": "/sector/8", "sectorId": 8, "created": 3},
]


def _menu_links(page, menu):
    menu.locator("summary").click()
    return dict(menu.locator("[data-bookmarks-panel] a").evaluate_all(
        "els => els.map(a => [a.querySelector('.bookmarks-name').textContent, a.getAttribute('href')])"))


def _assert_pick_links(links):
    assert links["Middling Sun"] == "/nav?from=system:701&to=system:702", links
    assert links["Fixture Second"] == "/sector/8?pick=to&from=system:701", links
    view = parse_qs(urlparse(links["Saved view"]).query)
    assert view == {"cam": ["3003.436762,5.294542628,-0.6473350564,151.006709,0.325344516,-0.3247703979,-0.6274082299,0.6285173408"], "pick": ["to"], "from": ["system:701"]}, links


def test_bookmarks_keep_a_course_pick_on_the_galaxy_map(page, map_site):
    """NAV.40: while a destination is picked, the Galaxy Map's bookmarks
    (menu and keys) keep the pick and the start already chosen, and the
    map's own URLs keep it as it moves."""
    _open_galaxy(page, map_site, "?pick=to&from=system:701")
    page.evaluate(SEED_BOOKMARKS, PICK_BOOKMARKS)
    _open_galaxy(page, map_site, "?pick=to&from=system:701")
    _assert_pick_links(_menu_links(page, page.locator("[data-bookmarks-menu]")))
    page.keyboard.press("Escape")

    # Key 1 opens the saved view, still picking.
    page.locator(GALAXY_CANVAS).focus()
    with page.expect_navigation():
        page.keyboard.press("1")
    _settle(page)
    params = parse_qs(urlparse(page.url).query)
    assert params["pick"] == ["to"] and params["from"] == ["system:701"] and "cam" in params
    assert page.locator(".pick-banner").count() == 1
    # A step on the map keeps the pick in the URL.
    shown = _crumbs(page)
    page.click('#galaxymap3d-controls [data-action="up"]')
    _settle(page)
    assert len(_crumbs(page)) < len(shown), (shown, _crumbs(page))
    params = parse_qs(_raw_query(page))
    assert params.get("pick") == ["to"] and params.get("from") == ["system:701"], params


def test_bookmarks_keep_a_course_pick_on_the_sector_page(page, map_site):
    """NAV.40: the sector page in pick mode offers bookmarks, each keeping
    the pick."""
    _open_sector(page, map_site)
    page.evaluate(SEED_BOOKMARKS, PICK_BOOKMARKS)
    _open_sector(page, map_site, "?pick=to&from=system:701")
    menu = page.locator(".pick-bookmarks[data-bookmarks-menu]")
    assert menu.count() == 1
    _assert_pick_links(_menu_links(page, menu))
    _open_sector(page, map_site)
    assert page.locator(".pick-bookmarks").count() == 0, "no extra menu outside pick mode"


# --- MAP.93, MAP.94: the breadcrumb on one line, and the phone's Steps button --------

CRUMB_LINE = """() => {
    const ol = document.querySelector("#galaxymap3d-crumbs ol");
    const items = Array.from(ol.children);
    const mids = items.map(li => { const r = li.getBoundingClientRect(); return (r.top + r.bottom) / 2; });
    const here = ol.querySelector("[aria-current]");
    return {cut: !!here && here.scrollWidth > here.clientWidth + 1,
            tops: Math.max(...mids) - Math.min(...mids) < 6 ? 1 : 2, overflow: ol.scrollWidth - ol.clientWidth, count: items.length,
            more: !!ol.querySelector(".crumb-menu"), visible: ol.offsetParent !== null,
            starTop: Math.round(document.querySelector("#galaxymap3d-crumbs .galaxy-bookmark").getBoundingClientRect().top),
            olTop: Math.round(ol.getBoundingClientRect().top), olHeight: ol.getBoundingClientRect().height};
}"""


def _open_deep_galaxy(page, map_site):
    """The Galaxy Map at the sector level of a fixture sector (many steps)."""
    _open_sector(page, map_site)
    href = page.locator(".page-actions a", has_text="Show on Galaxy Map").get_attribute("href")
    _open_galaxy(page, map_site, href[len("/galaxy"):])
    return _crumbs(page)


def _crumb_line_settled(page):
    """CRUMB_LINE once the breadcrumb has refolded after a resize (it
    does so a moment after the window changes, longer on a busy machine)."""
    line = page.evaluate(CRUMB_LINE)
    for _ in range(30):
        if line["tops"] == 1 and line["overflow"] <= 1 and not line["cut"]:
            break
        page.wait_for_timeout(100)
        line = page.evaluate(CRUMB_LINE)
    return line


def test_galaxy_breadcrumb_stays_on_one_line_folding_its_middle(page, map_site):
    """MAP.93: at any width the breadcrumb is one line: the first step,
    "…" for the steps that don't fit, the last ones and the current one;
    it fits again on resize, and "…" lists the folded steps."""
    page.set_viewport_size({"width": 1280, "height": 900})
    full = _open_deep_galaxy(page, map_site)
    assert len(full) >= 5, full
    for width in (1280, 900, 700, 600):
        page.set_viewport_size({"width": width, "height": 900})
        line = _crumb_line_settled(page)
        assert line["tops"] == 1 and line["overflow"] <= 1, (width, line)
        assert line["olHeight"] < 40, (width, line)
        assert not line["cut"], f"the current step is cut short at {width} px while steps could fold: {line}"
        shown = _crumbs(page)
        assert shown[0] == "Galaxy" and shown[-1] == full[-1], (width, shown)
        if line["more"]:
            assert line["count"] < len(full) + 1
    # At 600 px the middle is folded; "…" lists exactly those steps.
    line = page.evaluate(CRUMB_LINE)
    assert line["more"], "a deep breadcrumb at 600 px folds its middle"
    more = page.locator("#galaxymap3d-crumbs .crumb-menu")
    more.locator("summary").click()
    folded = more.locator("li").all_inner_texts()
    line_items = page.locator("#galaxymap3d-crumbs > ol > li").all_inner_texts()
    kept = line_items[2:]
    assert line_items[0] == "Home" and folded, line_items
    assert [line_items[0]] + folded + kept == ["Home"] + full, (line_items, folded, full)
    more.locator("button").first.click()
    _settle(page)
    assert _crumbs(page) == full[:1], "the first folded step goes there"
    assert not page.evaluate(CRUMB_LINE)["more"] or page.evaluate(CRUMB_LINE)["overflow"] <= 1
    # Wider again: everything shown, no "…".
    page.set_viewport_size({"width": 1280, "height": 900})
    page.wait_for_timeout(200)
    assert not page.evaluate(CRUMB_LINE)["more"] or page.evaluate(CRUMB_LINE)["overflow"] <= 1


def test_galaxy_phone_steps_button_between_back_and_forward(page, map_site):
    """MAP.94, UX.59: at phone width the breadcrumb line shows only the
    current level beside the star, and a labelled Steps button between Back
    and Forward lists every step with the current one marked; Whole galaxy
    stays beside the arrows."""
    page.set_viewport_size({"width": 390, "height": 844})
    full = _open_deep_galaxy(page, map_site)
    line = page.evaluate(CRUMB_LINE)
    assert line["visible"] and line["count"] == 1 or page.locator(
        "#galaxymap3d-crumbs > ol > li:visible").count() == 1, "only the current level on a phone"
    assert page.locator("#galaxymap3d-crumbs > ol > li:visible").inner_text() == full[-1]
    order = page.locator("#galaxymap3d-controls > *").evaluate_all(
        "els => els.map(e => e.dataset.action || e.id || e.className)")
    assert order.index("back") + 1 == order.index("galaxymap3d-steps") == order.index("forward") - 1, order
    assert page.locator('#galaxymap3d-controls [data-action="reset"]').is_visible()
    steps = page.locator("#galaxymap3d-steps")
    assert steps.locator("summary").inner_text().strip().endswith("Steps"), "labelled"
    steps.locator("summary").click()
    items = steps.locator("li").all_inner_texts()
    assert items == ["Home"] + full, (items, full)
    assert steps.locator("[aria-current]").inner_text() == full[-1]
    steps.locator("button", has_text="Galaxy").first.click()
    _settle(page)
    assert _crumbs(page) == ["Galaxy"]
    assert not steps.evaluate("d => d.open"), "taking a step closes the menu"
    # MAP.95, UX.58: Jump to newest is the last item of the Steps menu.
    page.locator('#galaxymap3d-controls [data-action="reset"]').click()
    _settle(page)
    current = steps.locator('[data-action="current"]')
    steps.locator("summary").click()
    assert current.is_visible() and current.is_disabled()
    steps.locator("summary").click()
    page.locator('#galaxymap3d-controls [data-action="back"]').click()
    _settle(page)
    steps.locator("summary").click()
    assert not current.is_disabled()
    assert current.bounding_box()["x"] + current.bounding_box()["width"] <= 390
    current.click()
    _settle(page)
    assert _crumbs(page) == ["Galaxy"] and current.is_disabled()
    assert not steps.evaluate("d => d.open"), "taking a step closes the menu"
    page.set_viewport_size({"width": 1280, "height": 900})
    page.wait_for_timeout(200)
    assert steps.is_visible() and page.evaluate(CRUMB_LINE)["visible"], "Steps shows at every width"


# --- NAV.31: hover lights every choice while picking a course -------------------------

CHOICE_SCAN = """() => {
    const canvas = document.querySelector("#galaxymap3d-canvas");
    const tip = document.querySelector("#galaxymap3d-tooltip");
    const box = canvas.getBoundingClientRect();
    const seen = new Set();
    for (let j = 2; j < 22; j++) {
        for (let i = 2; i < 22; i++) {
            canvas.dispatchEvent(new PointerEvent("pointermove", {clientX: box.left + (i / 24) * box.width,
                clientY: box.top + (j / 24) * box.height, bubbles: true, pointerType: "mouse"}));
            if (!tip.hidden) seen.add(tip.textContent.replace(/ \\(nothing charted here to pick\\)$/, ""));
        }
    }
    return Array.from(seen).sort();
}"""


def _choices_seen(page):
    seen = page.evaluate(CHOICE_SCAN)
    for _ in range(HOVER_RETRIES):  # the map ignores the pointer during a flight
        if seen:
            break
        page.wait_for_timeout(250)
        seen = page.evaluate(CHOICE_SCAN)
    return seen


def test_galaxy_map_hover_lights_every_choice_while_picking_a_course(page, map_site):
    """NAV.31: picking a NAV end, hovering lights the arc (and at each
    later stage the choice) under the pointer, generated or not (outside pick mode there are no blocks to light,
    MAP.297); only an empty one can't be taken."""
    _open_galaxy(page, map_site)
    _descend(page)
    stages = ["", "?" + _query(page)]
    for query in stages:
        _open_galaxy(page, map_site, query + ("&" if query else "?") + "pick=to&from=system:701")
        picking = _choices_seen(page)
        assert len(picking) > 1 and any("charted" in text for text in picking), (query, picking)
    # An empty arc lights but a click there stays put.
    _open_galaxy(page, map_site, "?pick=to&from=system:701")
    found = _hover_choice(page, r"nothing charted here to pick\)$")
    assert found, "no empty arc to hover"
    before = (_query(page), _crumbs(page))
    page.mouse.click(found[0], found[1])
    _settle(page)
    assert (_query(page), _crumbs(page)) == before


def test_galaxy_map_charted_only_dims_uncharted_objects_by_their_place_in_the_cut_order(page, halo_site):
    """MAP.252: with "Charted only" on, an uncharted star keeps a share of its light by the cut order (a star of a
    satellite galaxy the least, then the halo, then the ramp from small and dim to large and bright), a charted one
    keeps all, and the shader that applies the share and the distance rule compiled: the scene drew."""
    _open_galaxy(page, halo_site)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyUnchartedKeeps().entries.length > 0")
    page.click("#galaxymap3d-menu > summary")
    page.click('[data-action="charted-only"]')
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyUnchartedKeeps().on")
    keeps = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyUnchartedKeeps()")
    entries = keeps["entries"]
    by_id = {entry["id"]: entry for entry in entries}
    halo = [by_id[12000 + i] for i in range(40) if 12000 + i in by_id]
    assert halo, "the loose group is drawn"
    satellite = [e["keep"] for e in halo if e["tier"] == 2]
    halo_tier = [e["keep"] for e in halo if e["tier"] == 1]
    ramp = [e["keep"] for e in halo if e["tier"] == 0]
    assert satellite and halo_tier and ramp, (len(satellite), len(halo_tier), len(ramp))
    assert max(satellite) < min(halo_tier) and max(halo_tier) < min(ramp) and max(ramp) <= 0.5 + 1e-6, (satellite, halo_tier, ramp)
    brighter = sorted(((e["keep"], by_id[e["id"]]) for e in halo if e["tier"] == 0), key=lambda pair: pair[0])
    assert brighter[0][0] < brighter[-1][0], "the ramp runs from small and dim to large and bright"
    assert all(e["keep"] == 1 for e in entries if e["id"] < 12000), "a charted star keeps all of its light"
    assert 0.0 < keeps["far"] < 1.0


def test_galaxy_map_charted_only_works_while_picking_a_course(page, map_site):
    """NAV.32: choosing a NAV end, "Charted only" is a working toggle as in
    browsing, and a choice holding nothing generated still can't be taken."""
    _open_galaxy(page, map_site, "?pick=to&from=system:701")
    only = page.locator('#galaxymap3d-controls [data-action="charted-only"]')
    assert only.is_enabled() and only.get_attribute("aria-pressed") == "false"
    page.click("#galaxymap3d-menu > summary")
    only.click()
    assert only.get_attribute("aria-pressed") == "true"
    only.click()
    assert only.get_attribute("aria-pressed") == "false"
    found = _hover_choice(page, r"nothing charted here to pick\)$")
    assert found, "no empty arc to hover"
    before = (_query(page), _crumbs(page))
    page.mouse.click(found[0], found[1])
    _settle(page)
    assert (_query(page), _crumbs(page)) == before


FRAME = "() => document.querySelector('#galaxymap3d-canvas').galaxyFrame()"


def _framed(frame, least=0.75):
    """The stage's blocks all on the map, and filling much of it one way."""
    assert frame, "no frame"
    assert frame["left"] >= -1 and frame["top"] >= -1, frame
    assert frame["right"] <= frame["width"] + 1 and frame["bottom"] <= frame["height"] + 1, frame
    fill = max((frame["right"] - frame["left"]) / frame["width"], (frame["bottom"] - frame["top"]) / frame["height"])
    assert fill >= least, f"the stage fills only {fill:.0%} of the map: {frame}"
    return fill


@pytest.mark.parametrize("window", [None, (243, 0, 0, 0), (27, 27, 0, 0)])
def test_galaxy_map_frames_the_whole_stage_at_any_size_and_turn(page, map_site, window):
    """MAP.53 and MAP.78: every stage (the galaxy, a window of blocks, a window of smaller ones) is
    framed whole on the map, fitted to the map's own size, so a bigger
    window shows it bigger; it refits when the window is resized, and
    stays whole while Shift and the arrow keys turn it about its middle."""
    if window:
        _open_window(page, map_site, *window)
        # MAP.337: the map draws the stage the page loaded on; the window's own is what a load of the address gives.
        page.goto(page.url, wait_until="load")
        page.wait_for_selector(GALAXY_CANVAS, state="attached")
        _settle(page)
    else:
        _open_galaxy(page, map_site)
    if not window:
        # The galaxy opens zoomed in on the charted space (MAP.124); Re-center
        # is its whole fit.
        page.locator("#galaxymap3d-menu > summary").click()
        page.click('#galaxymap3d-controls [data-action="reset-view"]')
        _settle(page)
    small = page.evaluate(FRAME)
    _framed(small)
    page.set_viewport_size({"width": 1700, "height": 1300})
    page.wait_for_timeout(400)
    big = page.evaluate(FRAME)
    _framed(big)
    assert big["right"] - big["left"] > (small["right"] - small["left"]) * 1.15, (small, big)
    canvas = page.locator(GALAXY_CANVAS)
    canvas.focus()
    for key in ["Shift+ArrowLeft"] * 9 + ["Shift+ArrowDown"] * 4:
        page.keyboard.press(key)
    turned = page.evaluate(FRAME)
    # Turned about its middle, the picture may sit off center, but whole.
    _framed(turned, 0.5)
    assert abs(turned["left"] - big["left"]) + abs(turned["top"] - big["top"]) > 2, "Shift and the arrows turn the view"


LINES = "() => document.querySelector('#galaxymap3d-canvas').galaxyLines()"


def test_galaxy_map_draws_block_lines_only_in_a_window(page, map_site):
    """MAP.77: the whole galaxy draws no lines between its blocks; a window draws the divisions between its blocks,
    at every size down the ladder, nearly invisible at rest (MAP.202) and never lit by the pointer (MAP.297)."""
    _open_galaxy(page, map_site)
    seen = []
    for _ in range(12):
        lines = page.evaluate(LINES)
        seen.append((bool(parse_qs(_query(page)).get("win")), lines["blockEdges"]))
        if _descend(page) is None:
            break
    assert seen[0] == (False, 0), seen
    windows = [edges for is_window, edges in seen if is_window]
    assert len(windows) >= 2 and all(0 < edges <= 0.15 for edges in windows), ("nearly invisible at rest (MAP.202)", seen)
    page.mouse.move(*[v / 2 for v in (VIEWPORT["width"], VIEWPORT["height"])])
    page.wait_for_timeout(300)
    assert 0 < page.evaluate(LINES)["blockEdges"] <= 0.15, "the pointer over the map lights no block"


def open_dense_region(page, base, dist_pc=30.0):
    """TEST.368: the Galaxy Map of `dense_site` looking at its dense region from `dist_pc` parsecs, settled."""
    _open_galaxy(page, base, dense_region_query(dist_pc))
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyStarFill().count > 0", timeout=20000)
    # The page's header may leave part of the map below the window: bring the whole map on screen.
    page.evaluate("() => document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    _settle(page)


def test_dense_fixture_region_is_there_and_can_be_zoomed_into(page, dense_site):
    """TEST.368: the dense fixture galaxy holds a few thousand stars packed round one place, and a wheel scroll over it
    takes the view in until a single star of its cluster is the nearest thing to the middle of the screen."""
    _pin_window(page)
    page.set_viewport_size(VIEWPORT)
    open_dense_region(page, dense_site)
    filled = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyStarFill()")
    assert filled["count"] >= 1000, filled
    assert len(DENSE_STARS) >= 2500 and DENSE_POINTS and DENSE_WORLDS
    start = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
    page.locator(GALAXY_CANVAS).hover()
    for _ in range(12):
        page.mouse.wheel(0, -400)
        page.wait_for_timeout(100)
    _settle(page)
    end = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
    assert end["dist"] < start["dist"] / 3, (start["dist"], end["dist"])
    # MAP.300: the wheel goes on in on a star of the cluster (the one under the middle of the screen), so the cluster's centre may be
    # off the screen by then; what is under the middle is a star of the cluster.
    box = page.locator(GALAXY_CANVAS).bounding_box()
    middle = page.evaluate("([x, y]) => document.querySelector('#galaxymap3d-canvas').galaxyNearestThing(x, y)",
                           [box["x"] + box["width"] / 2, box["y"] + box["height"] / 2])
    assert middle["kind"] == "star" and math.dist(middle["point"], DENSE_CENTER_PC) < DENSE_HALO_RADIUS_PC, (middle, DENSE_CENTER_PC)


def _click_star(page, star):
    box = page.locator(GALAXY_CANVAS).bounding_box()
    page.mouse.move(box["x"] + star["x"], box["y"] + star["y"])
    page.mouse.click(box["x"] + star["x"], box["y"] + star["y"])


def _panel_state(page):
    return page.evaluate("() => document.querySelector('#galaxymap3d-info').dataset.state")


@pytest.mark.names_hidden
def test_a_click_in_the_dense_region_selects_a_star_and_never_a_sector_box(page, dense_site):
    """MAP.297, MAP.304: on the dense fixture region a click on a drawn star selects it (its panel shows in the map
    window), a click on nothing selects nothing and puts the panel away, and no sector box is outlined or shaded,
    hovered or selected."""
    page.set_viewport_size(VIEWPORT)
    open_dense_region(page, dense_site, 14.0)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyPickableStars().length > 20", timeout=20000)
    stars = [s for s in page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyPickableStars()") if s["visible"]]
    assert len(stars) > 20, len(stars)
    # A star with room round it: the click lands on it and on nothing else.
    chosen = next(s for s in stars if all(abs(s["x"] - o["x"]) + abs(s["y"] - o["y"]) > 30 for o in stars if o is not s))
    _click_star(page, chosen)
    # (A click waits a moment in case it is the first of two, MAP.226.)
    wanted = chosen["name"] or "Star"
    page.wait_for_function("(wanted) => { const h = document.querySelector('#galaxymap3d-info h3'); return h && h.textContent === wanted; }",
                           arg=wanted)
    assert page.evaluate("() => document.querySelector('#galaxymap3d-info').dataset.state") == "selected"
    lines = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyLines()")
    assert lines["hovered"]["lines"] == 0 and lines["hovered"]["shade"] == 0, lines

    # Nothing under the pointer: nothing is selected and the panel goes.
    box = page.locator(GALAXY_CANVAS).bounding_box()
    gaps = [(x, y) for x in range(40, 700, 53) for y in range(60, 520, 47)
            if all(abs(x - o["x"]) + abs(y - o["y"]) > 40 for o in stars)]
    # A spot the map itself is under: its panels (the selection's, the Menu) cover part of it.
    bare = page.evaluate("""(spots) => spots.filter(([x, y]) => {
        const found = document.elementFromPoint(x, y);
        return found && found.id === 'galaxymap3d-canvas';
    })""", [[box["x"] + x, box["y"] + y] for x, y in gaps])
    assert bare, "no spot of the map is free of its panels"
    page.mouse.click(*bare[0])
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info').dataset.state !== 'selected'")
    lines = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyLines()")
    assert lines["hovered"]["lines"] == 0, lines


def _visible_stars(page):
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyPickableStars().length > 20", timeout=20000)
    return [s for s in page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyPickableStars()") if s["visible"]]


def _free_of_panels(page, spots):
    """The spots (page coordinates) where the map itself, and none of its panels, is under the point."""
    return page.evaluate("""(spots) => spots.filter(([x, y]) => {
        const found = document.elementFromPoint(x, y);
        return found && found.id === 'galaxymap3d-canvas';
    })""", spots)


def test_a_selected_star_can_be_bookmarked(page, dense_site):
    """MAP.122: the panel of a picked star has a ☆ that saves its system as a bookmark (kind "system"), and a second press
    takes it out again."""
    page.set_viewport_size(VIEWPORT)
    open_dense_region(page, dense_site, 14.0)
    stars = _visible_stars(page)
    chosen = next(s for s in stars if all(abs(s["x"] - o["x"]) + abs(s["y"] - o["y"]) > 30 for o in stars if o is not s))
    _click_star(page, chosen)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info').dataset.state === 'selected'")
    star = page.locator("#galaxymap3d-info .map-info-bookmark")
    assert star.inner_text() == "☆" and star.is_enabled()
    star.click()
    saved = _bookmarks(page)
    assert len(saved) == 1 and saved[0]["kind"] == "system" and saved[0]["value"].startswith("system:"), saved
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info .map-info-bookmark').getAttribute('aria-pressed') === 'true'")
    page.once("dialog", lambda dialog: dialog.accept())   # taking a bookmark out asks first
    star.click()
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info .map-info-bookmark').getAttribute('aria-pressed') === 'false'")
    assert _bookmarks(page) == []


def _close_pair_on_the_map(page):
    """MAP.300: two fixture stars of the dense region under a tenth of a parsec apart that the map draws, in view, as stars of their own
    at the start (the nearest pair of the bright ones): `(a, b)`."""
    by_id = {star["id"]: star for star in DENSE_STARS}
    drawn = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyPickableStars()")
    stars = sorted((by_id[star["id"]] for star in drawn if star["visible"] and star["id"] in by_id and by_id[star["id"]]["luminosity_sol"] >= 1.0),
                   key=lambda star: star["x"])
    best = None
    for index, first in enumerate(stars):
        for second in stars[index + 1:index + 150]:
            if second["x"] - first["x"] > 0.1:
                break
            gap = math.dist((first["x"], first["y"], first["z"]), (second["x"], second["y"], second["z"]))
            if 0.02 < gap < 0.1 and (best is None or gap < best[0]):
                best = (gap, first, second)
    assert best, "the dense region has no two bright stars in view under a tenth of a parsec apart"
    return best[1], best[2]


def _spot(page, star):
    """The place of `star` on the page (client pixels): where a pointer is on it."""
    return page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", [star["x"], star["y"], star["z"]])


def _under(page, star):
    """What the map finds under the pointer placed on `star`'s spot: `galaxyNearestThing`'s `{kind, point}`."""
    spot = _spot(page, star)
    return page.evaluate("([x, y]) => document.querySelector('#galaxymap3d-canvas').galaxyNearestThing(x, y)", [spot["x"], spot["y"]])


def _wheel_in_on(page, star, steps, until=None, delta=-150):
    """Scrolls in with the pointer on `star`, each step, until `until(camera)` holds or `steps` are used; returns the camera."""
    for _ in range(steps):
        camera = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
        if until and until(camera):
            return camera
        spot = _spot(page, star)
        thing = page.evaluate("([x, y]) => document.querySelector('#galaxymap3d-canvas').galaxyNearestThing(x, y)", [spot["x"], spot["y"]])
        page.mouse.move(spot["x"], spot["y"])
        page.mouse.wheel(0, delta)
        page.wait_for_timeout(350)
        _settle(page)
    return page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")


def test_the_wheel_comes_close_enough_to_tell_two_stars_under_a_tenth_of_a_parsec_apart(page, dense_site):
    """MAP.300: scrolling in on a star of the dense region brings the camera down to where the star's neighbour, under a tenth of a parsec
    away, is a star of its own: each is the thing under the pointer on its own spot, tens of pixels apart, and the wheel goes on to the
    least distance the map allows (a thousandth of a parsec, `FIELD_FLOOR_PC`)."""
    page.set_viewport_size(VIEWPORT)
    open_dense_region(page, dense_site, 3.0)
    first, second = _close_pair_on_the_map(page)
    gap = math.dist((first["x"], first["y"], first["z"]), (second["x"], second["y"], second["z"]))
    assert gap < 0.1

    def apart(_camera):
        one, other = _spot(page, first), _spot(page, second)
        return math.hypot(one["x"] - other["x"], one["y"] - other["y"])

    # In on the first star, the pointer on it at each step, until the second stands clear of it.
    camera = _wheel_in_on(page, first, 40, until=lambda camera: apart(camera) >= 40)
    assert apart(camera) >= 40, (apart(camera), camera["dist"])
    # The map lists its stars to a coarser precision than the fixture's: the point is the star's, to a hundredth of a parsec at worst.
    for star, other in ((first, second), (second, first)):
        found = _under(page, star)
        near = math.dist(found["point"], (star["x"], star["y"], star["z"]))
        assert found["kind"] == "star" and near < 0.01 and near < math.dist(found["point"], (other["x"], other["y"], other["z"])), (found, star)
    # And on: the camera goes down to the least the map allows, and the first star is still the thing under the pointer on it.
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyZoomLimits()")["nearestPc"] <= 0.002
    camera = _wheel_in_on(page, first, 60, until=lambda camera: camera["dist"] <= 0.0011)
    assert camera["dist"] <= 0.0011, camera["dist"]
    assert math.dist(camera["eye"], (first["x"], first["y"], first["z"])) < gap, (camera["eye"], first)


def test_a_click_between_close_stars_offers_the_ones_within_reach(page, dense_site):
    """MAP.122: a click that has several stars inside the pointer's reach (10 px for a mouse) selects the nearest and lists
    the others under its panel; choosing another from the list shows that star instead."""
    page.set_viewport_size(VIEWPORT)
    open_dense_region(page, dense_site, 14.0)
    stars = _visible_stars(page)
    # Two stars 3 to 9 px apart, so a click on the first has the second inside the 10 px reach and still picks the first.
    def apart(a, b):
        return math.hypot(a["x"] - b["x"], a["y"] - b["y"])
    pairs = [(a, b) for a in stars for b in stars if a["id"] < b["id"] and 3 < apart(a, b) < 9
             and not any(c is not a and c is not b and apart(c, a) < apart(a, b) for c in stars)]
    assert pairs, "the dense region has no two stars a few pixels apart"
    # A pair that is not under a panel of the map.
    box = page.locator(GALAXY_CANVAS).bounding_box()
    free = set(map(tuple, _free_of_panels(page, [[box["x"] + a["x"], box["y"] + a["y"]] for a, _b in pairs])))
    first, second = next((a, b) for a, b in pairs if (box["x"] + a["x"], box["y"] + a["y"]) in free)
    page.mouse.click(box["x"] + first["x"], box["y"] + first["y"])
    _expand_card(page)    # the list is behind Expand
    page.wait_for_selector("#galaxymap3d-info .map-star-choices button", timeout=5000)
    buttons = page.locator("#galaxymap3d-info .map-star-choices button")
    assert buttons.count() >= 2
    pressed = page.locator("#galaxymap3d-info .map-star-choices button[aria-pressed='true']")
    assert pressed.count() == 1 and pressed.get_attribute("data-star-choice") == str(first["id"])
    other = page.locator(f"#galaxymap3d-info .map-star-choices button[data-star-choice='{second['id']}']")
    assert other.count() == 1
    other.click()
    page.wait_for_selector(f"#galaxymap3d-info .map-star-choices button[data-star-choice='{second['id']}'][aria-pressed='true']")
    assert page.locator("#galaxymap3d-info .map-star-choices button[aria-pressed='true']").count() == 1


def test_a_lone_star_has_no_choice_list(page, dense_site):
    """MAP.122: with only one star within reach there is nothing to choose, so its panel carries no list."""
    page.set_viewport_size(VIEWPORT)
    open_dense_region(page, dense_site, 14.0)
    stars = _visible_stars(page)
    chosen = next(s for s in stars if all(abs(s["x"] - o["x"]) + abs(s["y"] - o["y"]) > 40 for o in stars if o is not s))
    _click_star(page, chosen)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info').dataset.state === 'selected'")
    assert page.locator("#galaxymap3d-info .map-star-choices").count() == 0


def test_a_double_click_on_a_star_flies_to_it(page, dense_site):
    """MAP.304: a double click on a drawn star flies the camera there: its place becomes the middle of the view and the view
    comes in to a sector and a half."""
    page.set_viewport_size(VIEWPORT)
    open_dense_region(page, dense_site, 14.0)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyPickableStars().length > 20", timeout=20000)
    stars = [s for s in page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyPickableStars()") if s["visible"]]
    chosen = next(s for s in stars if all(abs(s["x"] - o["x"]) + abs(s["y"] - o["y"]) > 30 for o in stars if o is not s)
                  and abs(s["x"] - 400) + abs(s["y"] - 300) > 120)
    before = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
    box = page.locator(GALAXY_CANVAS).bounding_box()
    page.mouse.dblclick(box["x"] + chosen["x"], box["y"] + chosen["y"])
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=20000)
    _settle(page)
    after = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
    star = next(s for s in DENSE_STARS if s["id"] == chosen["id"])
    target = after["target"]
    assert max(abs(target[0] - star["x"]), abs(target[1] - star["y"]), abs(target[2] - star["z"])) < 0.05, (target, star)
    assert after["dist"] <= before["dist"] + 1e-6, (before["dist"], after["dist"])


def test_stars_glow_brighter_they_do_not_grow_bigger_and_their_colour_follows_their_temperature(page, dense_site):
    """MAP.241: two stars of very different luminosity are drawn the same size, within a narrow band, but the bright one
    glows harder (and its core is as opaque or more), a larger star has a slightly larger core that stays small, and a
    hotter star is bluer than a cooler one."""
    page.set_viewport_size(VIEWPORT)
    open_dense_region(page, dense_site, 14.0)
    looks = [look for look in page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyStarLooks()")
             if not look["phenomenon"] and look["luminosity"] and look["temperature"]]
    assert len(looks) > 100, len(looks)
    dim = min(looks, key=lambda look: look["luminosity"])
    bright = max(looks, key=lambda look: look["luminosity"])
    assert bright["luminosity"] > 1000 * dim["luminosity"], (dim["luminosity"], bright["luminosity"])
    assert abs(bright["sizePx"] - dim["sizePx"]) < 3, (dim["sizePx"], bright["sizePx"])
    assert bright["glow"] > 2 * dim["glow"], (dim["glow"], bright["glow"])
    assert bright["bright"] >= dim["bright"]
    assert all(7 <= look["sizePx"] <= 13 and look["corePx"] <= 2.5 for look in looks), "every star stays small"
    hot = max(looks, key=lambda look: look["temperature"])
    cool = min(looks, key=lambda look: look["temperature"])
    assert hot["temperature"] > cool["temperature"] + 1000
    assert hot["color"][2] - hot["color"][0] > cool["color"][2] - cool["color"][0], (hot["color"], cool["color"])


@pytest.mark.parametrize("radius_pc", [4.0, 9.0])
def test_fly_to_fit_puts_the_sphere_in_the_window_with_no_sector_box(page, dense_site, radius_pc):
    """MAP.319: flyToFit(sphere) puts the camera where the sphere fills the smaller side of the view with 8 percent to
    spare (or as near as the view's zoom allows), over the sphere's middle, and no block of sectors is drawn once the
    flight is over."""
    page.set_viewport_size(VIEWPORT)
    open_dense_region(page, dense_site, 40.0)
    before = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
    center = [before["target"][0] + 3.0, before["target"][1] - 2.0, before["target"][2]]
    page.evaluate("([center, radius]) => document.querySelector('#galaxymap3d-canvas').galaxyGoTo({center: center, radius: radius})",
                  [center, radius_pc])
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=20000)
    _settle(page)
    after = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
    size = page.evaluate("() => { const c = document.querySelector('#galaxymap3d-canvas'); return [c.clientWidth, c.clientHeight]; }")
    fov = math.radians(page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyFov()"))
    tan = math.tan(fov / 2)
    half = math.atan(min(tan, tan * size[0] / size[1]))
    wanted = max(radius_pc, 2.0) / math.sin(half) * 1.08
    assert max(abs(a - b) for a, b in zip(after["target"], center)) < 0.05, (after["target"], center)
    assert after["dist"] == pytest.approx(wanted, rel=0.02) or after["dist"] > wanted, (after["dist"], wanted)
    # The sphere is in the window: its edge, seen from the camera, lies inside half the smaller field.
    assert math.asin(min(1.0, radius_pc / after["dist"])) <= half + 1e-6, (after["dist"], radius_pc)
    lines = page.evaluate(LINES)
    assert lines["chartedLines"] == 0 and lines["hovered"]["lines"] == 0, lines


LABELS = "() => document.querySelector('#galaxymap3d-canvas').__mapLabels.shown()"


def _label_spot(page, label):
    """The page point at the middle of a placed name. A card open over the map is closed first (MAP.332: a tall card's More takes the pointer)."""
    close = page.locator("#galaxymap3d-info .map-info-close")
    if close.count() and close.is_visible():
        close.click()
        page.wait_for_timeout(100)
    box = page.locator(GALAXY_CANVAS).bounding_box()
    return box["x"] + label["left"] + label["width"] / 2, box["y"] + label["top"] + label["height"] / 2


def _on_the_canvas(page, spot):
    """Whether the map's canvas is what the page shows at `spot` (the info panel lies over part of the map)."""
    return page.evaluate("([x, y]) => { const e = document.elementFromPoint(x, y); return !!e && e.id === 'galaxymap3d-canvas'; }", list(spot))


LABEL_RECTS = """() => [...document.querySelector('.map-labels').children].filter((e) => e.style.opacity !== '0').map((e) => {
    const r = e.getBoundingClientRect(); return {text: e.textContent, kind: e.dataset.kind, left: r.left, top: r.top, right: r.right, bottom: r.bottom};
})"""


def _label_gaps(rects):
    """The smallest clear space between two names' boxes (negative where they overlap), or None for fewer than two."""
    gaps = []
    for i, a in enumerate(rects):
        for b in rects[i + 1:]:
            gaps.append(max(b["left"] - a["right"], a["left"] - b["right"], b["top"] - a["bottom"], a["top"] - b["bottom"]))
    return min(gaps) if gaps else None


def _named_stars_by_brightness(page):
    """The dense region's named stars that fall inside the map now, brightest first: [(luminosity, name, x, y)] in page pixels."""
    box = page.locator(GALAXY_CANVAS).bounding_box()
    found = []
    for star in DENSE_STARS:
        if not star["name"]:
            continue
        spot = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", [star["x"], star["y"], star["z"]])
        if spot and box["x"] + 20 <= spot["x"] <= box["x"] + box["width"] - 20 and box["y"] + 20 <= spot["y"] <= box["y"] + box["height"] - 20:
            found.append((star["luminosity_sol"], star["name"], spot["x"], spot["y"]))
    return sorted(found, reverse=True)


@pytest.mark.viewport(1000, 900)
@pytest.mark.parametrize("dist_pc, buffer_px", ((6.0, 4), (14.0, 4), (30.0, 4), (14.0, 12)))
def test_names_in_the_dense_region_keep_a_buffer_and_the_brightest_in_a_crowd_keeps_its_name(page, dense_site, dist_pc, buffer_px):
    """MAP.299: in the dense fixture region no two names overlap or stand closer than the buffer (4 px, the `map.labelBuffer`
    preference, here also 12), and where names must be cut the brighter star keeps its name: the brightest named star in view is
    named, and a bright star that lost its name lost it to a name that is no dimmer (another kind of name, or a star at least as bright)."""
    if buffer_px != 4:
        page.add_init_script(f"try {{ localStorage.setItem('planetgen.pref.map.labelBuffer', '{buffer_px}'); }} catch (e) {{}}")
    open_dense_region(page, dense_site, dist_pc)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapLabels.shown().length > 10", timeout=20000)
    page.wait_for_timeout(2500)
    rects = page.evaluate(LABEL_RECTS)
    assert len(rects) > 10, rects
    gap = _label_gaps(rects)
    assert gap >= buffer_px - 0.5, (gap, buffer_px)
    brightness = {s["name"]: s["luminosity_sol"] for s in DENSE_STARS if s["name"]}
    shown = {r["text"] for r in rects}
    ranked = _named_stars_by_brightness(page)
    assert ranked, "no named star of the dense region is in view"
    assert ranked[0][1] in shown, ("the brightest named star in view lost its name", ranked[0], sorted(shown))
    for luminosity, name, x, y in ranked[:10]:
        if name in shown:
            continue
        # Where its name would stand (a name of this length is about 74 px by 16, above the star). A bright star that lost its name lost it
        # to a name in the way that outranks it: one of another kind (a sector, a world) or a star at least as bright, never only to dimmer stars.
        left, top = x - 37 - buffer_px, y - 22 - buffer_px
        in_the_way = [rect for rect in rects
                      if rect["left"] < left + 74 + 2 * buffer_px and left < rect["right"]
                      and rect["top"] < top + 16 + 2 * buffer_px and top < rect["bottom"]]
        justified = [rect for rect in in_the_way if rect["text"] not in brightness or brightness[rect["text"]] >= luminosity]
        assert justified, (name, luminosity, [(r["text"], r["kind"], brightness.get(r["text"])) for r in in_the_way])


def test_names_stand_on_stars_and_other_objects_and_a_click_on_one_flies_to_its_item(page, dense_site):
    """MAP.298: the dense region's named stars and its named objects (rogue planets, a neutron star, a black hole)
    carry their names on the map, and a click on a name selects its item (the panel names it) and flies the camera there."""
    page.set_viewport_size(VIEWPORT)
    open_dense_region(page, dense_site, 14.0)
    page.wait_for_function("""() => {
        const kinds = new Set(document.querySelector('#galaxymap3d-canvas').__mapLabels.shown().map(l => l.kind));
        return kinds.has('star') && kinds.has('object');
    }""", timeout=20000)
    shown = page.evaluate(LABELS)
    stars = {s["name"]: s for s in DENSE_STARS if s["name"]}
    objects = {p["name"]: p for p in DENSE_POINTS}
    for kind, known in (("star", stars), ("object", objects)):
        label = next(l for l in shown if l["kind"] == kind and l["text"] in known and _on_the_canvas(page, _label_spot(page, l)))
        item = known[label["text"]]
        page.mouse.move(*_label_spot(page, label))
        page.wait_for_function("(id) => !!document.querySelector('.map-label[data-hot]')", arg=label["id"])
        page.mouse.click(*_label_spot(page, label))
        page.wait_for_function("(name) => { const h = document.querySelector('#galaxymap3d-info h3'); return h && h.textContent === name; }",
                               arg=label["text"])
        page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=20000)
        _settle(page)
        target = page.evaluate(CAMERA)["target"]
        assert max(abs(target[0] - item["x"]), abs(target[1] - item["y"]), abs(target[2] - item["z"])) < 0.05, (label["text"], target, item)
        # Back to the start for the next kind.
        open_dense_region(page, dense_site, 14.0)
        page.wait_for_function("""() => {
            const kinds = new Set(document.querySelector('#galaxymap3d-canvas').__mapLabels.shown().map(l => l.kind));
            return kinds.has('star') && kinds.has('object');
        }""", timeout=20000)
        shown = page.evaluate(LABELS)


def test_galaxy_map_grid_never_adds_up_to_more_than_the_cap_while_zooming_in(page, map_site):
    """MAP.188: the blocks' faces and lines are drawn at a share that keeps the cells along any
    ray from adding up to more than the cap (50%), at every stage and at the deepest zoom."""
    problems = []
    page.on("console", lambda message: problems.append(message.text) if message.type == "error" else None)
    _open_galaxy(page, map_site)
    reading = "() => document.querySelector('#galaxymap3d-canvas').galaxyGridOpacity()"
    readings = []
    for _ in range(12):
        readings.append(page.evaluate(reading))
        if _descend(page) is None:
            break
    page.locator(GALAXY_CANVAS).hover()
    for _ in range(14):
        page.mouse.wheel(0, -400)
        page.wait_for_timeout(120)
        readings.append(page.evaluate(reading))
    assert len(readings) > 4
    assert not [text for text in problems if "shader" in text.lower() or "WebGLProgram" in text], problems
    for entry in readings:
        assert 0 < entry["scale"] <= 1, entry
        assert entry["worst"] <= entry["cap"] + 1e-6, entry
    deepest = readings[-1]
    assert deepest["worst"] <= 0.5 + 1e-6, deepest


def test_galaxy_map_charted_only_outlines_the_charted_blocks(page, map_site):
    """MAP.111: "Charted only" outlines the blocks holding generated
    sectors on the picked arc, and takes the outlines away when it is off."""
    _open_galaxy(page, map_site)
    _descend(page)
    assert page.evaluate(LINES)["chartedLines"] == 0
    page.click("#galaxymap3d-menu > summary")
    page.click('[data-action="charted-only"]')
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyLines().chartedLines > 0")
    page.click('[data-action="charted-only"]')
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyLines().chartedLines === 0")


def test_galaxy_map_draws_no_blocks_unless_the_cell_grid_or_a_fill_is_on(page, map_site):
    """MAP.265: the blocks (translucent prisms) are not drawn by default, in the whole galaxy or in a window; the
    Menu's Cell grid draws them, as do Charted only and a Color by fill, and turning each off takes them away."""
    _open_galaxy(page, map_site)
    drawn = "() => document.querySelector('#galaxymap3d-canvas').galaxyLines().blocksDrawn"
    assert page.evaluate(drawn) == 0
    _descend(page)
    assert parse_qs(_query(page)).get("win"), "a window is open"
    assert page.evaluate(drawn) == 0
    page.click("#galaxymap3d-menu > summary")
    for action in ("cell-grid", "charted-only"):
        page.click(f'[data-action="{action}"]')
        page.wait_for_function(f"() => ({drawn})() > 0")
        assert page.get_attribute(f'[data-action="{action}"]', "aria-pressed") == "true"
        page.click(f'[data-action="{action}"]')
        page.wait_for_function(f"() => ({drawn})() === 0")


CAMERA = "() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()"


def test_galaxy_map_turns_past_edge_on_and_under_the_plane(page, map_site):
    """MAP.96: the view turns any way by any amount: past 80 degrees,
    through edge-on and under the plane, and the window stays framed."""
    _open_window(page, map_site, 243, 0, 0, 0)
    canvas = page.locator(GALAXY_CANVAS)
    canvas.focus()
    tilts = []
    for _ in range(24):
        page.keyboard.press("Shift+ArrowDown")
        tilts.append(page.evaluate(CAMERA)["tilt"])
    assert max(tilts) > 170, tilts
    assert any(80 < t < 100 for t in tilts), "turned through edge-on"
    assert page.evaluate(CAMERA)["tilt"] > 120, "seen from under the plane"
    _framed(page.evaluate(FRAME), 0.4)


BELOW_THE_PLANE = r"latitude -[\d.]+°–(?:-[\d.]+|0\.0)°"
"""str: A tooltip naming a block that lies wholly under the galactic plane (its latitude range ends at or below 0)."""


def test_galaxy_map_camera_presets_at_each_zoom_step(page, map_site):
    """MAP.97: the whole galaxy is seen straight down, a window of blocks
    at the isometric slant, going in and coming back out; a manual turn
    holds only within its step."""
    iso = math.degrees(math.atan(math.sqrt(2)))
    _open_galaxy(page, map_site)
    seen = []
    for _ in range(12):
        tilt = page.evaluate(CAMERA)["tilt"]
        windowed = bool(_stage(page)["win"])
        seen.append((windowed, round(tilt, 1)))
        assert tilt == pytest.approx(iso if windowed else 0, abs=0.5), (tilt, seen)
        if windowed:
            # A manual turn doesn't carry over to the next step.
            page.locator(GALAXY_CANVAS).focus()
            page.keyboard.press("Shift+ArrowDown")
        if _descend(page) is None:
            break
    assert seen[0][0] is False and sum(1 for windowed, _ in seen if windowed) >= 2, seen
    # Back out: each step flies back to its own preset.
    _open_window(page, map_site, 27, 27, 0, 0)
    assert page.evaluate(CAMERA)["tilt"] == pytest.approx(iso, abs=0.5)
    page.locator(GALAXY_CANVAS).focus()
    for _ in range(3):  # 27, 81, 243 and the galaxy
        page.keyboard.press("Escape")
        _settle(page)
    assert _crumbs(page) == ["Galaxy"]
    assert page.evaluate(CAMERA)["tilt"] == pytest.approx(0, abs=0.5)


# --- MAP.80: every star of the sector at sector zoom ------------------------------------

def test_galaxy_map_fetches_the_sector_nodes_around_a_sector_it_shows(page, map_site):
    """MAP.80, MAP.261: shown a sector, the map also fetches the finest nodes around it (they list every star and
    point phenomenon), not just the view's coarser ones."""
    urls = []
    page.on("request", lambda request: urls.append(request.url) if "/galaxy/nodes" in request.url else None)
    # Straight to the sector: a visit to its page first would leave the
    # nodes in the browser's cache.
    _open_galaxy(page, map_site, f"?sector={PRIME_CODE}")
    page.wait_for_timeout(2500)
    requests = [[int(key.split(".")[0]) for key in parse_qs(urlparse(url).query)["keys"][0].split(",")]
                for url in urls]
    assert requests, "the map asked for nodes"
    assert any(1 in levels for levels in requests), requests
    assert any(min(levels) > 1 for levels in requests), requests


def _point_colored_spots(page, selector, limit=80):
    """Screen points of the spots no star could be: colors off the
    blackbody line (a star's red, green and blue always run in order), as
    the black holes, neutron stars and quasars are drawn (MAP.80)."""
    image = Image.open(io.BytesIO(_shot(page, selector))).convert("RGB")
    box = page.locator(selector).bounding_box()
    width, height = image.size
    pixels = image.load()
    spots = []
    for y in range(2, height - 2):
        for x in range(2, width - 2):
            r, g, b = pixels[x, y]
            off = max(g - max(r, b), min(r, b) - g)  # green above or below both others
            if off > 20 and max(r, g, b) > 150:
                spots.append((off, x, y))
            elif b > 200 and r < 110 and g < 110:
                # MAP.164: a neutron star's deep blue, darker than any star's blue-white.
                spots.append((b - r - g, x, y))
    spots.sort(reverse=True)
    kept = []
    for _off, x, y in spots:
        if all(abs(x - kx) + abs(y - ky) > 8 for kx, ky in kept):
            kept.append((x, y))
        if len(kept) >= limit:
            break
    return [(box["x"] + x * box["width"] / width, box["y"] + y * box["height"] / height) for x, y in kept]


def test_galaxy_map_draws_point_phenomena_that_link_to_their_pages(page, map_site):
    """MAP.80: black holes, neutron stars and quasars are drawn among the
    stars in their own colors, and a click on one names it and links to
    its phenomenon page."""
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_timeout(1000)
    _poi_layers(page)    # a crowd's mark of the points of interest would take the click meant for the pulsar
    page.locator("#galaxymap3d-menu > summary").click()
    # One step in, so the glow of the bright stars around them (every star of
    # a node comes in together, MAP.152) leaves the dark blue neutron star clear.
    # The wheel zooms toward what is under the pointer (MAP.225): Shift-click the
    # neutron star to center the view on it, then step in with the pointer on the center.
    pulsar = next(point for point in POINT_FIELD if point["type"] == "neutron_star")
    spot = _screen_of(page, (pulsar["x"], pulsar["y"], pulsar["z"]))
    page.keyboard.down("Shift")
    page.mouse.click(*spot)
    page.keyboard.up("Shift")
    _camera(page)
    middle = page.locator(GALAXY_CANVAS).bounding_box()
    page.mouse.move(middle["x"] + middle["width"] / 2, middle["y"] + middle["height"] / 2)
    # A step in is as far as the nearest thing along the ray allows (MAP.334), so on the wide map the three points still stand within a
    # few pixels of one another after one: steps in until each is 40 px from the next, which a click can tell apart (TEST.445).
    for _ in range(8):
        page.mouse.wheel(0, -400)
        page.wait_for_timeout(600)
        _settle(page)
        page.wait_for_timeout(1500)
        spots = [_screen_of(page, (point["x"], point["y"], point["z"])) for point in POINT_FIELD]
        if all(spots) and min(math.dist(a, b) for i, a in enumerate(spots) for b in spots[i + 1:]) >= 40:
            break
    names = {point["name"]: point for point in POINT_FIELD}
    found = {}
    url = page.url
    # They are points of light now (MAP.208), too small for a colour scan to find reliably: click
    # where the map draws each.
    for point in POINT_FIELD:
        spot = _screen_of(page, (point["x"], point["y"], point["z"]))
        assert spot, point["name"]
        # The pulsar sits inside the fixture nebula over the field, whose core takes a click on the line
        # through it first, and a click beside it can land on another point of the map (TEST.459): the pointer is moved over the
        # pixels around the spot until the tooltip names the point, and only that click is made.
        tip = page.locator("#galaxymap3d-tooltip")
        for dx, dy in ((0, 0), (4, 0), (-4, 0), (0, 4), (0, -4), (4, 4), (-4, -4), (4, -4), (-4, 4)):
            page.mouse.move(spot[0] + dx, spot[1] + dy)
            page.wait_for_timeout(100)
            if not (tip.is_visible() and tip.inner_text().startswith(point["name"])):
                continue
            page.wait_for_timeout(DOUBLE_CLICK_GAP_MS)   # clicks this close in time would be a double click (MAP.226)
            page.mouse.click(spot[0] + dx, spot[1] + dy)
            heading = page.locator("#galaxymap3d-info h3")
            if heading.count() and heading.inner_text() in names:
                found_point = names[heading.inner_text()]
                link = page.locator("#galaxymap3d-info a", has_text="View phenomenon")
                found[found_point["type"]] = (link.get_attribute("href"), page.locator("#galaxymap3d-info").inner_text())
                break
        if page.url != url:
            page.go_back()
            _settle(page)
    assert set(found) == {"black_hole", "neutron_star", "quasar"}, found
    for kind, (href, text) in found.items():
        point = next(p for p in POINT_FIELD if p["type"] == kind)
        assert href and str(point["id"]) in href, (kind, href)
    assert "Millisecond pulsar" in found["neutron_star"][1]


MARKERS = """() => {
    const layer = document.getElementById("galaxymap3d-canvas").__mapMarkers;
    return layer ? layer.entries().map((e) => ({kind: e.kind.id, name: e.item.name, regionKey: e.item.regionKey || null, count: e.item.count || null, pulled: !!e.pulled, members: (e.item.members || []).map((m) => m.name), x: e.x, y: e.y, z: e.z})) : null;
}"""


def _markers_until(page, done, seconds=20):
    """The marker layer's entries once `done(entries)` holds: the layer redraws on the next frame, which a busy machine
    renders late."""
    deadline = time.monotonic() + seconds
    while True:
        found = page.evaluate(MARKERS)
        if done(found) or time.monotonic() > deadline:
            return found
        page.wait_for_timeout(100)


def _hole_names(page):
    """The names of the black holes the marks stand for (a crowd's mark, MAP.340, for each of its members), sorted."""
    names = []
    for entry in page.evaluate(MARKERS):
        if entry["kind"] == "marker-black-hole":
            names.extend(entry["members"] or [entry["name"]])
    return sorted(names)


def _hole_view(page, base, density=None, dist_pc=10.0):
    """The Galaxy Map of `holes_site` looking at its star field from `dist_pc` parsecs (a 10 pc zoom), the density preference set to
    `density` pixels first when given; returns the black hole marks drawn once they have settled, and the view's facts."""
    if density is not None:
        page.add_init_script(f"try {{ localStorage.setItem('planetgen.pref.map.blackHoleDensity', '{density}'); }} catch (e) {{}}")
    x, y, z = STAR_FIELD_CENTER_PC
    _open_galaxy(page, base, f"?cam={x:.6f},{y:.6f},{z:.6f},{dist_pc:.6f},0,0,0,1")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapMarkers.entries().length > 0", timeout=20000)
    last = None
    for _ in range(30):
        names = _hole_names(page)
        if names == last:
            break
        last = names
        page.wait_for_timeout(300)
    view = page.evaluate("""() => { const c = document.querySelector('#galaxymap3d-canvas');
        return {dist: c.galaxyCamera().dist, width: c.clientWidth, height: c.clientHeight}; }""")
    return last, view


def _hole_rule(hole, dist):
    """MAP.339's rule in Python: within 3 view distances, or heavy enough for how far (5 solar masses at the reach, rising as the square)."""
    cx, cy, cz = STAR_FIELD_CENTER_PC
    r = math.dist((hole["x"], hole["y"], hole["z"]), (cx, cy, cz))
    return r <= 3 * dist or hole["mass_solar"] >= 5 * (r / (3 * dist)) ** 2


def test_black_holes_are_the_near_ones_and_the_heavy_enough_at_a_10_pc_zoom_and_the_density_key_sets_the_cap(page, holes_site):
    """MAP.339: at a 10 pc zoom of a galaxy listing 183 black holes, the marks are those within 3 view distances of what the camera
    looks at and those heavy enough for how far they are, no more than one per 120 by 120 pixels of the screen (the default of the
    `map.blackHoleDensity` preference), and a smaller side in the key lets more of the near ones show."""
    _pin_window(page)
    page.set_viewport_size(VIEWPORT)
    drawn, view = _hole_view(page, holes_site)
    by_name = {hole["name"]: hole for hole in HOLE_FIELD}
    accepted = {name for name, hole in by_name.items() if _hole_rule(hole, view["dist"])}
    assert drawn and set(drawn) <= accepted, sorted(set(drawn) - accepted)
    assert not [name for name in drawn if name.startswith("Far Hole")], "no far light hole is drawn"
    cap = int(view["width"] * view["height"] // (120 * 120))
    assert len(drawn) <= cap, (len(drawn), cap)
    near = [name for name in drawn if name.startswith("Near Hole")]
    assert len(near) >= 1 and len([n for n in accepted if n.startswith("Near Hole")]) > cap - 3, "the near holes fill the cap"
    assert {"Heavy Hole 1", "Heavy Hole 2"} <= set(drawn) or len(drawn) == cap, "a rare heavy hole follows the same rule, ranked by mass over distance squared"
    # A smaller square per mark lets more of the near ones show.
    page.evaluate("() => localStorage.setItem('planetgen.pref.map.blackHoleDensity', '50')")
    page.reload(wait_until="load")
    _settle(page)
    denser, _view = _hole_view_after_reload(page)
    assert len(denser) > len(drawn) and set(denser) <= accepted, (len(denser), len(drawn))
    assert len(denser) <= int(view["width"] * view["height"] // (50 * 50))


def _hole_view_after_reload(page):
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapMarkers.entries().length > 0", timeout=20000)
    last = None
    for _ in range(30):
        names = _hole_names(page)
        if names == last:
            break
        last = names
        page.wait_for_timeout(300)
    return last, None


def test_crowded_black_holes_are_one_mark_with_a_count_and_a_list_and_split_on_zoom_in(page, holes_site):
    """MAP.340: from 200 pc the thirty near black holes of `holes_site` lie within a few pixels of each other, so they are drawn as one
    mark with their count (a badge on the glyph), its card lists them each with its mass and has a Go to them button, and after
    zooming in the crowd has split into more, smaller crowds or single marks."""
    page.set_viewport_size(VIEWPORT)
    drawn, view = _hole_view(page, holes_site, density=30, dist_pc=200.0)

    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-black-hole" and e["count"] for e in entries))
    crowds = [e for e in found if e["kind"] == "marker-black-hole" and e["count"]]
    assert crowds, found
    crowd = max(crowds, key=lambda e: e["count"])
    assert crowd["count"] >= 5, crowds
    button = page.locator(".starmap-sr-markers button", has_text=f"{crowd['count']} black holes")
    assert button.count() >= 1
    button.first.evaluate("b => b.click()")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    assert page.locator("#galaxymap3d-info h3").inner_text() == f"{crowd['count']} black holes"
    members = page.locator("#galaxymap3d-info .map-info-crowd-system")
    assert members.count() == crowd["count"]
    assert "Near Hole" in members.first.inner_text() or "Heavy Hole" in members.first.inner_text() or "Far Hole" in members.first.inner_text()
    zoom = page.locator("#galaxymap3d-info button", has_text="Go to them")
    assert zoom.count() == 1
    zoom.click()
    _settle(page)
    after = _markers_until(page, lambda entries: not any(e["kind"] == "marker-black-hole" and e["count"] and e["count"] >= crowd["count"] for e in entries))
    holes_after = [e for e in after if e["kind"] == "marker-black-hole"]
    assert not any(e["count"] and e["count"] >= crowd["count"] for e in holes_after), "the crowd has split"
    assert len(holes_after) >= 2, "the holes of the crowd are drawn as smaller crowds and single marks"


def _spot_clear_of_names(page, spot, tooltip_text):
    """A pixel within a mark's reach that names it on hover: a name label drawn over the mark takes the pointer first (the map's
    label engine places names by the room the canvas has, so a mark may lie under one), and the mark's own spot a few pixels
    above or to the side is clear of it."""
    for dx, dy in ((0, 0), (0, -6), (6, -6), (-6, -6), (0, -9), (8, 0), (-8, 0), (0, 6), (0, 9)):
        page.mouse.move(spot[0] + dx - 30, spot[1] + dy - 30)
        page.mouse.move(spot[0] + dx, spot[1] + dy)
        page.wait_for_timeout(120)
        if tooltip_text in page.locator("#galaxymap3d-tooltip").inner_text():
            return (spot[0] + dx, spot[1] + dy)
    raise AssertionError(f"no pixel near {spot} names {tooltip_text!r} on hover")


def test_galaxy_map_black_hole_markers_have_a_toggle_a_legend_line_and_can_be_picked(page, map_site):
    """MAP.219: the marker layer draws a glyph over each black hole in view, with a Menu toggle that rides in the
    address, a legend line while on, and a click on the glyph shows the black hole's card."""
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    _poi_layers(page)    # a crowd of a black hole with a world's mark has no name of its own (MAP.305), which the long runs hit
    page.locator("#galaxymap3d-menu > summary").click()
    holes = {p["name"]: p for p in POINT_FIELD if p["type"] == "black_hole"}
    drawn = _markers_until(page, lambda found: {e["name"] for e in found if e["kind"] == "marker-black-hole"} >= set(holes))
    assert {e["name"] for e in drawn if e["kind"] == "marker-black-hole"} >= set(holes), drawn
    legend = page.locator("#galaxymap3d-layer-legend")
    assert "Black holes" in legend.inner_text()
    button = page.locator('#galaxymap3d-kinds button[data-layer="marker-black-hole"]')
    _open_menu(page)
    assert button.count() == 1 and button.get_attribute("aria-pressed") == "true"
    # The keyboard and screen readers get a button per marker.
    page.locator(".starmap-sr-markers button", has_text="Fixture Maw").evaluate("b => b.click()")
    assert page.locator("#galaxymap3d-info h3").inner_text() == "Fixture Maw"
    page.keyboard.press("Escape")
    close = page.locator("#galaxymap3d-info .map-info-close")
    if close.is_visible():   # MAP.256: the card over the map covers the points under it until it is closed
        close.click()
    # A click on the glyph names the black hole.
    maw = holes["Fixture Maw"]
    if page.evaluate("document.querySelector('#galaxymap3d-menu').open"):   # the Menu's panel would cover the glyph (TEST.394)
        page.locator("#galaxymap3d-menu > summary").click()
    page.wait_for_timeout(600)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    spot = _screen_of(page, (maw["x"], maw["y"], maw["z"]))
    assert spot
    spot = _spot_clear_of_names(page, spot, "Fixture Maw")
    page.mouse.click(*spot)
    assert page.locator("#galaxymap3d-info h3").inner_text() == "Fixture Maw"
    # Off: no markers, no legend line, and the address says so.
    _open_menu(page)
    button.click()
    assert "marker-black-hole" in parse_qs(urlparse(page.url).query).get("hide", [""])[0].split(",")
    assert [e for e in _markers_until(page, lambda found: not any(e["kind"] == "marker-black-hole" for e in found))
            if e["kind"] == "marker-black-hole"] == []
    assert "Black holes" not in legend.inner_text()
    # A reload keeps it off; on again brings them back.
    page.reload()
    page.wait_for_selector("#galaxymap3d-steps [data-steps-panel] li", state="attached")
    _settle(page)
    page.wait_for_timeout(800)
    assert [e for e in page.evaluate(MARKERS) if e["kind"] == "marker-black-hole"] == []
    _open_menu(page)
    page.locator('#galaxymap3d-kinds button[data-layer="marker-black-hole"]').click()
    assert {e["name"] for e in _markers_until(page, lambda found: bool(found))} >= set(holes)
    assert "marker-black-hole" not in parse_qs(urlparse(page.url).query).get("hide", [""])[0].split(",")


GLOW = "() => { const g = document.querySelector('#galaxymap3d-canvas').__mapGlow; return g && {count: g.count(), drawn: g.drawn(), visible: g.isVisible()}; }"


def test_galaxy_map_glow_shows_the_unlisted_light_and_the_menu_switches_it(page, map_site):
    """MAP.266: a node's light beyond the stars it lists is a soft glow, drawn from the first view, a layer of the Menu
    that is on until the visitor turns it off (and `hide=glow` in the address says so)."""
    _open_galaxy(page, map_site)
    page.wait_for_function("() => { const g = document.querySelector('#galaxymap3d-canvas').__mapGlow; return g && g.drawn() > 0; }",
                           timeout=20000)
    state = page.evaluate(GLOW)
    assert state["count"] > 0 and state["visible"], state
    button = page.locator('#galaxymap3d-kinds button[data-layer="glow"]')
    _open_menu(page)
    assert button.count() == 1 and button.get_attribute("aria-pressed") == "true"
    button.click()
    assert "hide=glow" in page.url
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapGlow.isVisible() === false")
    _open_menu(page)
    button.click()
    assert "hide=glow" not in page.url
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapGlow.isVisible() === true")


def test_the_glow_fades_out_as_the_stars_resolve_and_never_adds_up(page, dense_site):
    """MAP.324: from the whole galaxy the unresolved light is a soft glow, but looking into a crowd of stars it is gone
    (it flooded the map white and hid the stars), and glows that overlap take the stronger, never the sum."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, dense_site)
    page.wait_for_function("() => { const g = document.querySelector('#galaxymap3d-canvas').__mapGlow; return g && g.drawn() > 0; }",
                           timeout=20000)
    far = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').__mapGlow.strongest()")
    assert 0.05 < far <= 0.4 + 1e-9, far
    # Blending: the maximum of the premultiplied glows (THREE.MaxEquation is 104), not the sum.
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').__mapGlow.mesh.material.blendEquation") == 104
    for dist in (14.0, 3.0):
        open_dense_region(page, dense_site, dist)
        page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapGlow.count() > 0", timeout=20000)
        _settle(page)
        try:
            page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapGlow.strongest() < 0.01", timeout=20000)
        except Exception:
            raise AssertionError((dist, page.evaluate("() => [document.querySelector('#galaxymap3d-canvas').__mapGlow.strongest(), document.querySelector('#galaxymap3d-canvas').galaxyCamera().dist]")))


NEBULAE_STATE = ("() => { const g = document.querySelector('#galaxymap3d-canvas').__mapNebulae; return g && {count: g.count(), "
                 "drawn: g.drawn(), visible: g.isVisible(), blobs: g.blobs().map(b => ({key: b.key, family: b.family, count: b.count, "
                 "ring: b.ring, color: b.color}))}; }")


def test_nebula_regions_are_soft_groups_in_the_family_colour_and_the_menu_switches_each_family(page, map_site):
    """MAP.176: the groups of nebulae the drawn nodes list (`nebula_regions`) are soft sprites at the galactic view, each in the
    colour of its dominant family with a ring when it holds more than one nebula, a layer of the Menu for each family (on until the
    visitor says otherwise, `hide=nebulae-emission` in the address), with a legend line each and one that says what a ring is."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, map_site)
    page.wait_for_function("() => { const g = document.querySelector('#galaxymap3d-canvas').__mapNebulae; return g && g.drawn() > 0; }",
                           timeout=20000)
    state = page.evaluate(NEBULAE_STATE)
    assert state["visible"] and state["count"] >= 1, state
    families = {blob["family"] for blob in state["blobs"]}
    assert families <= {"emission", "reflection", "diffuse"} and len(families) >= 2, families
    assert any(blob["ring"] for blob in state["blobs"]), "a group of several nebulae has a ring"
    by_family = {}
    for blob in state["blobs"]:
        by_family.setdefault(blob["family"], set()).add(tuple(blob["color"]))
    assert all(len(colors) == 1 for colors in by_family.values()) and len({c for v in by_family.values() for c in v}) == len(by_family)
    legend = page.locator("#galaxymap3d-layer-legend")
    assert legend.locator('[data-layer="nebulae-emission"]').inner_text() == "Emission nebulae"
    assert "zoom in to separate" in legend.locator('[data-layer="nebulae-groups"]').inner_text()
    _open_menu(page)
    button = page.locator('#galaxymap3d-kinds button[data-layer="nebulae-emission"]')
    assert button.count() == 1 and button.get_attribute("aria-pressed") == "true"
    assert page.locator('#galaxymap3d-kinds button[data-layer="nebulae-groups"]').count() == 0
    button.click()
    assert "hide=nebulae-emission" in page.url
    page.wait_for_function("""() => !document.querySelector('#galaxymap3d-canvas').__mapNebulae.blobs().some(b => b.family === 'emission')""")
    assert legend.locator('[data-layer="nebulae-emission"]').count() == 0
    assert page.evaluate(NEBULAE_STATE)["drawn"] > 0, "the other families stay"
    _open_menu(page)
    button.click()
    assert "hide=nebulae-emission" not in page.url
    page.wait_for_function("""() => document.querySelector('#galaxymap3d-canvas').__mapNebulae.blobs().some(b => b.family === 'emission')""")


def test_at_most_six_hundred_nebula_regions_are_drawn_out_of_three_thousand_candidates(page, map_site):
    """MAP.176: given 3,000 groups the layer draws the strongest 600 and no more."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, map_site)
    page.wait_for_function("() => { const g = document.querySelector('#galaxymap3d-canvas').__mapNebulae; return g && g.drawn() > 0; }",
                           timeout=20000)
    drawn = page.evaluate("""() => {
        const canvas = document.querySelector('#galaxymap3d-canvas');
        const layer = canvas.__mapNebulae;
        const here = canvas.galaxyCamera().target;
        const span = canvas.galaxyCamera().dist;
        let seed = 7;
        const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
        const list = [];
        for (let i = 0; i < 3000; i++) {
            list.push({key: 'x' + i, block: layer.blobs()[0].block, ring: i % 2 === 0, color: [0.9, 0.4, 0.5],
                       x: here[0] + (random() - 0.5) * span, y: here[1] + (random() - 0.5) * span, z: here[2] + (random() - 0.5) * span,
                       radiusPc: span * (0.03 + random() * 0.05), count: 2 + Math.floor(random() * 90), family: 'emission'});
        }
        layer.setBlobs(list, [0, 0, 0]);
        return list.length;
    }""")
    assert drawn == 3000
    page.wait_for_timeout(500)
    seen = page.evaluate("() => { const g = document.querySelector('#galaxymap3d-canvas').__mapNebulae; return [g.count(), g.drawn()]; }")
    assert seen[0] == 3000 and 0 < seen[1] <= 600, seen


REGION_SPOTS = "() => document.querySelector('#galaxymap3d-canvas').galaxyRegionSpots()"
REGION_KEYS = ["9.0.0.0", "9.1.0.0"]
"""The fixture's two charted regions (the 9-blocks holding its sectors), by cell."""


def _region_pixels(page, spot, half=6):
    """The picture of the canvas about the region's middle `spot` (canvas pixels), as a flat list of RGB tuples."""
    image = Image.open(io.BytesIO(_shot(page, GALAXY_CANVAS))).convert("RGB")
    width, height = image.size
    box = page.locator(GALAXY_CANVAS).bounding_box()
    scale = width / box["width"]
    x, y = int(spot["x"] * scale), int(spot["y"] * scale)
    pixels = image.load()
    return [pixels[i, j] for i in range(max(0, x - half), min(width, x + half)) for j in range(max(0, y - half), min(height, y + half))]


def test_the_charted_regions_are_soft_spheres_at_galactic_scale_colored_by_share_with_no_sector_box(page, map_site):
    """MAP.316: at the whole-galaxy view every charted region (a 9-block with a generated sector, however little of it) is a soft
    sphere at its extent, colored by its share charted, drawn by a layer of the Menu; no sector box is drawn with it."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, map_site)
    page.wait_for_function("() => { const g = document.querySelector('#galaxymap3d-canvas').__mapRegions; return g && g.drawn() > 0; }", timeout=20000)
    spots = page.evaluate(REGION_SPOTS)
    assert sorted(spot["key"] for spot in spots) == REGION_KEYS
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').__mapRegions.count()") == 2
    assert all(spot["alpha"] > 0.15 and spot["radiusPx"] >= 5 for spot in spots), spots
    by_key = {spot["key"]: spot for spot in spots}
    assert by_key["9.0.0.0"]["share"] > by_key["9.1.0.0"]["share"] > 0    # Fixture Beside (MAP.329) is the second sector of the first
    assert by_key["9.1.0.0"]["color"] != by_key["9.0.0.0"]["color"], "the colour follows the share"
    assert page.evaluate(LINES)["chartedLines"] == 0, "no sector box at galactic scale"
    # The layer is what is drawn there: the picture about a region changes when the Menu's Charted regions is switched off.
    visible = next(spot for spot in spots if spot["visible"])
    with_regions = _region_pixels(page, visible)
    _open_menu(page)
    button = page.locator('#galaxymap3d-kinds button[data-layer="regions"]')
    assert button.count() == 1 and button.get_attribute("aria-pressed") == "true"
    button.click()
    assert "hide=regions" in page.url
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapRegions.isVisible() === false")
    assert _region_pixels(page, visible) != with_regions
    _open_menu(page)
    button.click()
    assert "hide=regions" not in page.url
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapRegions.isVisible() === true")


def test_a_charted_region_with_an_inhabited_world_carries_an_inhabited_token_that_opens_its_card(page, map_site):
    """MAP.321 c: at galactic scale a region whose list entry has inhabited worlds carries the inhabited symbol at its middle
    (the Menu's Inhabited worlds layer, on by default), in the colour of its best world's tier; picking it shows the region's card with
    the inhabited line, and a region without any has no token."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, map_site)
    page.wait_for_function("() => { const g = document.querySelector('#galaxymap3d-canvas').__mapRegions; return g && g.drawn() > 0; }", timeout=20000)
    found = _markers_until(page, lambda entries: any(e["regionKey"] for e in entries))
    tokens = [e for e in found if e["regionKey"]]
    assert [(e["kind"], e["regionKey"]) for e in tokens] == [("marker-inhabited", "9.0.0.0")], found
    assert tokens[0]["name"] == "Fixture Beez Region", tokens
    spot = next(s for s in page.evaluate(REGION_SPOTS) if s["key"] == "9.0.0.0")
    assert spot["visible"]
    box = page.locator(GALAXY_CANVAS).bounding_box()
    where = (box["x"] + spot["x"], box["y"] + spot["y"])
    assert _on_the_canvas(page, where)
    page.mouse.move(where[0] + 1, where[1])
    page.mouse.move(*where)
    page.wait_for_timeout(300)
    assert "Fixture Beez Region, 3 inhabited worlds in 2 systems, best Blue" in page.locator("#galaxymap3d-tooltip").inner_text()
    page.mouse.click(*where)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    info = page.locator("#galaxymap3d-info").inner_text()
    assert "3 inhabited worlds in 2 systems, best Blue" in info, info
    # Switching the Inhabited worlds layer off takes the token, and the regions layer off takes it with its sphere.
    _open_menu(page)
    page.locator('#galaxymap3d-kinds button[data-layer="marker-inhabited"]').click()
    _markers_until(page, lambda entries: not any(e["regionKey"] for e in entries))
    assert not [e for e in page.evaluate(MARKERS) if e["regionKey"]]
    page.locator('#galaxymap3d-kinds button[data-layer="marker-inhabited"]').click()
    _markers_until(page, lambda entries: any(e["regionKey"] for e in entries))
    page.locator('#galaxymap3d-kinds button[data-layer="regions"]').click()
    _markers_until(page, lambda entries: not any(e["regionKey"] for e in entries))
    assert not [e for e in page.evaluate(MARKERS) if e["regionKey"]]


def test_the_regions_give_way_to_the_stars_as_the_camera_comes_in(page, map_site):
    """MAP.316: a region's sphere is strongest from far and fades as the camera comes in (the exact curve is the unit test's),
    and with a sector opened in place the layer is out of the picture."""
    _open_galaxy(page, map_site)
    page.wait_for_function("() => { const g = document.querySelector('#galaxymap3d-canvas').__mapRegions; return g && g.drawn() > 0; }", timeout=20000)
    spot = lambda: next(s for s in page.evaluate(REGION_SPOTS) if s["key"] == "9.0.0.0")
    far = spot()
    _open_window(page, map_site, 9, 0, 0, 0, g=3)
    page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyGoTo({center: [2.83, 2.83, 0], radius: 2})")
    _settle(page)
    near = spot()
    assert near["radiusPx"] > 4 * far["radiusPx"], (far, near)
    assert near["alpha"] < far["alpha"], (far, near)
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapRegions.isVisible() === false", timeout=20000)


REGION_FRAME = "() => { const f = document.querySelector('#galaxymap3d-canvas').__mapRegionFrame; return f && {framed: f.state().framed, positions: f.positions()}; }"


def test_the_galactic_view_frames_the_region_it_is_in_with_a_fixed_wireframe_and_dims_nothing(page, map_site):
    """MAP.327: in the galactic view a locked wireframe stands round the charted region the view is in (the pointer moving
    changes nothing of it), another region is framed once the view is in that, the stars are not dimmed, and opening a
    sector hands the framing to the sector."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, map_site)
    _open_window(page, map_site, 9, 1, 0, 0, g=3)
    page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyGoTo({center: [28.88, 21.5, 0], radius: 12})")
    _settle(page)
    page.wait_for_function("() => { const f = document.querySelector('#galaxymap3d-canvas').__mapRegionFrame; return f && f.state().framed; }", timeout=20000)
    first = page.evaluate(REGION_FRAME)
    assert len(first["positions"]) > 0
    # MAP.331: the legend says what the outline is while it is drawn.
    page.wait_for_function("() => document.getElementById('galaxymap3d-layer-legend').textContent.includes('Outline:')", timeout=20000)
    assert "(9 sectors a side)" in page.locator("#galaxymap3d-layer-legend").text_content()
    assert page.locator("#galaxymap3d-layer-legend .map-frame-swatch").count() == 1
    assert page.evaluate(SECTOR_VIEW)["outsideDim"] == 1, "nothing is dimmed in the galactic view"
    page.mouse.move(300, 420)
    page.mouse.move(520, 480)
    page.wait_for_timeout(300)
    assert page.evaluate(REGION_FRAME) == first
    page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyGoTo({center: [2.83, 2.83, 0], radius: 12})")
    _settle(page)
    page.wait_for_function("(old) => JSON.stringify(document.querySelector('#galaxymap3d-canvas').__mapRegionFrame.positions()) !== old",
                           arg=__import__("json").dumps(first["positions"]), timeout=20000)
    assert page.evaluate(REGION_FRAME)["framed"]
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapRegionFrame.state().framed === false", timeout=20000)
    assert "Outline:" not in page.locator("#galaxymap3d-layer-legend").text_content(), "no region outline, so no legend line for it"
    assert page.evaluate(SECTOR_FRAME)["framed"]


def test_the_menu_lists_the_charted_regions_by_share_or_name_and_a_click_flies_to_one(page, map_site):
    """MAP.316: the Menu's Charted regions list is the always-available way to known space: by share charted (most first) or by
    name, a button each that flies to the region with the camera framing its extent."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, map_site)
    page.wait_for_function("() => { const l = document.querySelector('#galaxymap3d-canvas').__mapRegionList; return l && l.order().length === 2; }", timeout=20000)
    _open_menu(page)
    details = page.locator("#galaxymap3d-regions")
    assert details.is_visible()
    assert "(2)" in details.locator("summary").inner_text()
    details.locator("summary").click()
    order = lambda: page.evaluate("() => document.querySelector('#galaxymap3d-canvas').__mapRegionList.order()")
    by_share = order()
    assert [text.split(",")[0] for text in by_share] == ["Fixture Beez Region", "Fixture Zarh Region"], by_share
    assert all(text.endswith("under 1 percent charted") for text in by_share), by_share
    details.locator("[data-regions-sort]").select_option("name")
    assert [text.split(",")[0] for text in order()] == ["Fixture Beez Region", "Fixture Zarh Region"]
    details.locator("[data-regions-sort]").select_option("share")
    assert [text.split(",")[0] for text in order()] == ["Fixture Beez Region", "Fixture Zarh Region"]
    details.locator('button[data-region-key="9.0.0.0"]').click()
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=30000)
    _settle(page)
    from tests.map_site_support import REGION_RADIUS_PC
    after = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
    centre = next(spot for spot in [{"x": 2.8284271247461903, "y": 2.82842712474619, "z": 0.0}])
    assert max(abs(after["target"][i] - centre["xyz"[i]]) for i in range(3)) < 0.05, after["target"]
    assert after["dist"] < 6 * REGION_RADIUS_PC, after["dist"]


def test_a_click_on_a_region_shows_its_card_and_a_double_click_flies_to_it(page, map_site):
    """MAP.316: a region's sphere takes the click (its card: how much is charted, a Go to region button) and the double click."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, map_site)
    page.wait_for_function("() => { const g = document.querySelector('#galaxymap3d-canvas').__mapRegions; return g && g.drawn() > 0; }", timeout=20000)
    _poi_layers(page)    # the inhabited token of the other region lies over this one's middle
    _open_menu(page)
    page.locator("#galaxymap3d-menu > summary").click()
    box = page.locator(GALAXY_CANVAS).bounding_box()
    spot = next(s for s in page.evaluate(REGION_SPOTS) if s["key"] == "9.1.0.0" and s["visible"])
    where = (box["x"] + spot["x"], box["y"] + spot["y"])
    assert _on_the_canvas(page, where)
    page.mouse.move(*where)
    page.wait_for_timeout(300)
    assert "Fixture Zarh Region" in page.locator("#galaxymap3d-tooltip").inner_text()
    page.mouse.click(*where)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    assert page.locator("#galaxymap3d-info h3").inner_text() == "Fixture Zarh Region"
    assert "under 1 percent charted" in page.locator("#galaxymap3d-info").inner_text()
    page.locator("#galaxymap3d-info button", has_text="Go to region").click()
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=30000)
    _settle(page)
    after = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
    assert math.dist(after["target"], [28.876, 21.498, 0.0]) < 0.1, after["target"]


def test_galaxy_map_opens_with_the_glow_off_when_the_address_says_so(page, map_site):
    _open_galaxy(page, map_site, "?hide=glow")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapGlow !== undefined")
    assert page.evaluate(GLOW)["visible"] is False


def test_galaxy_map_habitable_world_markers_are_on_by_default_and_name_the_system(page, map_site):
    """MAP.220, MAP.303: systems with a habitable world get a marker of their own kind, on until the Menu turns it off (and
    `hide=` in the address says so), a legend line while on, and a click on one shows how many worlds it holds; the
    inhabited, habitable and gas giant layers begin on and the uninhabitable and asteroid field ones off."""
    _open_star_field(page, map_site)
    # MAP.333: the flight lands where it was asked, so frame the black hole (15 pc out) as well.
    page.evaluate("([c, r]) => document.querySelector('#galaxymap3d-canvas').galaxyGoTo({center: c, radius: r})",
                  [list(STAR_FIELD_CENTER_PC), 16])
    _settle(page)
    page.wait_for_timeout(1000)
    wanted = {"marker-black-hole", *POI_DEFAULT_ON}
    found = _markers_until(page, lambda entries: wanted <= {e["kind"] for e in entries})
    assert {e["kind"] for e in found} == wanted, found
    assert "show=" not in page.url and "hide=marker" not in page.url
    systems = {s["name"]: s for s in WORLD_FIELD}
    button = page.locator('#galaxymap3d-kinds button[data-layer="marker-habitable"]')
    _open_menu(page)
    assert button.count() == 1
    for layer in POI_LAYERS:
        pressed = page.locator(f'#galaxymap3d-kinds button[data-layer="{layer}"]').get_attribute("aria-pressed")
        assert pressed == ("true" if layer in POI_DEFAULT_ON else "false"), layer
    button.click()
    page.wait_for_timeout(400)
    assert "hide=marker-habitable" in page.url
    _markers_until(page, lambda entries: not any(e["kind"] == "marker-habitable" for e in entries))
    assert not [e for e in page.evaluate(MARKERS) if e["kind"] == "marker-habitable"]
    _poi_layers(page, "marker-habitable")
    assert "marker-habitable" not in page.url and "hide=marker-inhabited,marker-gas-giant" in page.url, page.url
    drawn = {e["name"] for e in _markers_until(page, lambda found: any(e["kind"] == "marker-habitable" for e in found))
             if e["kind"] == "marker-habitable"}
    assert drawn == set(systems), drawn
    assert "Habitable worlds" in page.locator("#galaxymap3d-layer-legend").inner_text()
    page.locator("#galaxymap3d-menu > summary").click()
    page.wait_for_timeout(300)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    eden = systems["Fixture Eden"]
    spot = _screen_of(page, (eden["x"], eden["y"], eden["z"]))
    assert spot
    page.mouse.move(*spot)
    page.mouse.move(spot[0] + 1, spot[1])
    page.wait_for_timeout(200)
    assert "Fixture Eden, 2 habitable worlds, best Green" in page.locator("#galaxymap3d-tooltip").inner_text()
    page.mouse.click(*spot)
    assert page.locator("#galaxymap3d-info h3").inner_text() == "Fixture Eden"
    assert "2" in page.locator("#galaxymap3d-info").inner_text()


def test_galaxy_map_marks_gas_giants_and_a_system_with_two_reasons_once_and_lists_them_in_its_infobox(page, map_site):
    """MAP.305: a gas giant is a reason for a mark of its own, which the Menu can turn on alone; a system with a gas giant and a
    habitable world is drawn once however many reasons are on, a system with neither has no mark, and the infobox of a mark lists
    every reason its system has."""
    _open_star_field(page, map_site)
    page.wait_for_timeout(1000)
    giants = {s["name"] for s in GIANT_FIELD}
    worlds = {s["name"] for s in WORLD_FIELD}
    button = page.locator('#galaxymap3d-kinds button[data-layer="marker-gas-giant"]')
    _poi_layers(page, "marker-gas-giant")
    assert button.count() == 1 and button.get_attribute("aria-pressed") == "true"
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-gas-giant" for e in entries))
    assert {e["name"] for e in found if e["kind"] == "marker-gas-giant"} == giants
    assert not [e for e in found if e["kind"] == "marker-habitable"]
    page.locator('#galaxymap3d-kinds button[data-layer="marker-habitable"]').click()
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-habitable" for e in entries))
    marks = [e["name"] for e in found if e["kind"] in ("marker-habitable", "marker-gas-giant")]
    assert sorted(marks) == sorted(giants | worlds), marks    # Fixture Eden once, though it has both reasons
    page.locator("#galaxymap3d-menu > summary").click()
    page.wait_for_timeout(300)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    eden = next(s for s in WORLD_FIELD if s["name"] == "Fixture Eden")
    spot = _screen_of(page, (eden["x"], eden["y"], eden["z"]))
    assert spot
    page.mouse.move(*spot)
    page.mouse.move(spot[0] + 1, spot[1])
    page.mouse.click(*spot)
    page.wait_for_function("() => (document.querySelector('#galaxymap3d-info h3') || {}).textContent === 'Fixture Eden'", timeout=10000)
    _expand_card(page)    # all the reasons are behind Expand
    info = page.locator("#galaxymap3d-info").inner_text()
    assert "2 habitable worlds, best Green" in info and "1 gas giant" in info, info


def _fly_to_system(page, system):
    """Opens the window round the sector a fixture system is in, so the nodes of its own sector are in (the class markers come from
    nodes of 3 shells or fewer, MAP.220), and waits until the flight ends."""
    from planetgen.galaxy.spherical import sector_address_at

    shell, band, wedge = sector_address_at((system["x"], system["y"], system["z"]), EDGE_PC)
    page.evaluate("([cell, g]) => document.querySelector('#galaxymap3d-canvas').galaxyOpenWindow(cell, g)",
                  [{"m": 1, "shell": shell, "band": band, "wedge": wedge}, 1])
    _settle(page)


def test_galaxy_map_class_markers_are_off_until_chosen_one_class_at_a_time_and_come_in_with_the_sector(page, map_site):
    """MAP.220: every planet class has a marker of its own in the Menu's World classes group, all off until the visitor picks them;
    choosing a class marks the systems that hold a world of it once the map is close enough to have their sector (its glyph carries the
    class letter, its hover text counts them), a system a world reason that is on already marks is marked once, and turning the class
    off takes its marks away."""
    _open_star_field(page, map_site)
    page.wait_for_timeout(1000)
    _open_menu(page)
    group = page.locator("#galaxymap3d-kinds")
    assert "World classes" in group.inner_text()
    assert group.locator('button[data-layer^="marker-class-"]').count() == 20
    assert all(group.locator(f'button[data-layer="marker-class-{letter}"]').get_attribute("aria-pressed") == "false" for letter in "KOGJ")
    eden = next(s for s in WORLD_FIELD if s["name"] == "Fixture Eden")
    rock = next(s for s in BARREN_FIELD if s["name"] == "Fixture Rock")
    _fly_to_system(page, eden)
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-habitable" for e in entries))
    assert not [e for e in found if e["kind"].startswith("marker-class-")], "no class is marked until one is chosen"
    group.locator('button[data-layer="marker-class-K"]').click()
    group.locator('button[data-layer="marker-class-G"]').click()
    found = _markers_until(page, lambda entries: any(e["kind"].startswith("marker-class-") for e in entries), seconds=8)
    assert not [e for e in found if e["kind"].startswith("marker-class-")], "Fixture Eden is already marked by its habitable worlds"
    assert "show=marker-class-G,marker-class-K" in page.url or "show=marker-class-K,marker-class-G" in page.url, page.url
    _poi_layers(page)    # no world reasons on: a system the classes list is named by its class
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-class-K" for e in entries))
    assert [e["name"] for e in found if e["kind"] == "marker-class-K"] == ["Fixture Eden"], found
    assert "Class K worlds" in page.locator("#galaxymap3d-layer-legend").inner_text()
    page.locator("#galaxymap3d-menu > summary").click()
    page.wait_for_timeout(300)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    spot = _screen_of(page, (eden["x"], eden["y"], eden["z"]))
    assert spot
    page.mouse.move(*spot)
    page.mouse.move(spot[0] + 1, spot[1])
    page.wait_for_timeout(200)
    assert "Fixture Eden, 2 class K worlds" in page.locator("#galaxymap3d-tooltip").inner_text()
    page.mouse.click(*spot)
    _expand_card(page)
    info = page.locator("#galaxymap3d-info").inner_text()
    assert "2 class K worlds" in info and "1 class O world" in info, info
    _fly_to_system(page, rock)
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-class-G" for e in entries))
    assert [e["name"] for e in found if e["kind"] == "marker-class-G"] == ["Fixture Rock"], found
    _open_menu(page)
    group.locator('button[data-layer="marker-class-G"]').click()
    found = _markers_until(page, lambda entries: not any(e["kind"] == "marker-class-G" for e in entries))
    assert not [e for e in found if e["kind"] == "marker-class-G"]


def test_galaxy_map_life_class_markers_are_off_until_chosen_and_mark_the_systems_where_the_class_lives(page, map_site):
    """MAP.351: the Menu has a Life classes group with a marker for each of the 15 life classes, all off until the visitor picks them;
    choosing one marks the systems with a world where that class lives once the map is close enough to have their sector, and the hover
    text counts the worlds."""
    _open_star_field(page, map_site)
    page.wait_for_timeout(1000)
    _open_menu(page)
    group = page.locator("#galaxymap3d-kinds")
    assert "Life classes" in group.inner_text()
    assert group.locator('button[data-layer^="marker-lifeclass-"]').count() == 15
    assert all(group.locator(f'button[data-layer="marker-lifeclass-{class_id}"]').get_attribute("aria-pressed") == "false"
               for class_id in ("1a", "2", "6", "9"))
    dawn = next(s for s in LIFE_FIELD if s["name"] == "Fixture Dawn")
    ember = next(s for s in LIFE_FIELD if s["name"] == "Fixture Ember")
    _fly_to_system(page, dawn)
    found = _markers_until(page, lambda entries: bool(entries))
    assert not [e for e in found if e["kind"].startswith("marker-lifeclass-")], "no life class is marked until one is chosen"
    group.locator('button[data-layer="marker-lifeclass-2"]').click()
    group.locator('button[data-layer="marker-lifeclass-6"]').click()
    assert "show=marker-lifeclass-2,marker-lifeclass-6" in page.url or "show=marker-lifeclass-6,marker-lifeclass-2" in page.url, page.url
    _poi_layers(page)    # no world reasons on: a system the life classes list is named by its class
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-lifeclass-2" for e in entries))
    assert [e["name"] for e in found if e["kind"] == "marker-lifeclass-2"] == ["Fixture Dawn"], found
    assert "Class 2 life" in page.locator("#galaxymap3d-layer-legend").inner_text()
    _fly_to_system(page, ember)
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-lifeclass-6" for e in entries))
    assert [e["name"] for e in found if e["kind"] == "marker-lifeclass-6"] == ["Fixture Ember"], found
    _open_menu(page)
    group.locator('button[data-layer="marker-lifeclass-6"]').click()
    found = _markers_until(page, lambda entries: not any(e["kind"] == "marker-lifeclass-6" for e in entries))
    assert not [e for e in found if e["kind"] == "marker-lifeclass-6"]


def test_galaxy_map_life_filters_show_only_the_systems_with_life_sentient_life_or_a_civilization(page, map_site):
    """MAP.311: the Menu's Points of interest group has an optional marker for worlds with life, with sentient life and with a civilization,
    each off until turned on; each shows only the systems that have it, a system with several is drawn once by the strongest sign, and
    the card lists the reasons."""
    _open_star_field(page, map_site)
    page.wait_for_timeout(1000)
    _open_menu(page)
    layers = {"life": "marker-life", "sentient": "marker-sentient", "civilization": "marker-civilization"}
    for layer in layers.values():
        assert page.locator(f'#galaxymap3d-kinds button[data-layer="{layer}"]').get_attribute("aria-pressed") == "false", layer
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-habitable" for e in entries))
    assert not [e for e in found if e["kind"] in layers.values()], "none of them is on until it is chosen"
    named = {reason: {s["name"] for s in LIFE_FIELD if s["reason"] == reason} for reason in layers}
    for reason, layer in layers.items():
        _poi_layers(page)    # no other reason on
        button = page.locator(f'#galaxymap3d-kinds button[data-layer="{layer}"]')
        button.click()
        found = _markers_until(page, lambda entries: {e["name"] for e in entries if e["kind"] == layer} == named[reason])
        assert {e["name"] for e in found if e["kind"] == layer} == named[reason], (reason, found)
        assert {e["name"] for e in found if e["kind"] not in (layer, "marker-black-hole")} == set(), found
        button.click()
    _poi_layers(page)
    page.locator('#galaxymap3d-kinds button[data-layer="marker-life"]').click()
    page.locator('#galaxymap3d-kinds button[data-layer="marker-sentient"]').click()
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-sentient" for e in entries) and
                           len([e for e in entries if e["kind"] in ("marker-life", "marker-sentient")]) == 4)
    marks = [(e["name"], e["kind"]) for e in found if e["kind"] in ("marker-life", "marker-sentient")]
    assert sorted(marks) == [("Fixture Eden", "marker-life"), ("Fixture Jove", "marker-sentient"), ("Fixture Rest", "marker-sentient"),
                             ("Fixture Rock", "marker-life")], marks
    page.locator("#galaxymap3d-menu > summary").click()
    page.wait_for_timeout(300)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    rest = next(s for s in LIFE_FIELD if s["name"] == "Fixture Rest")
    spot = _screen_of(page, (rest["x"], rest["y"], rest["z"]))
    assert spot
    page.mouse.move(*spot)
    page.mouse.move(spot[0] + 1, spot[1])
    page.mouse.click(*spot)
    _expand_card(page)
    info = page.locator("#galaxymap3d-info").inner_text()
    assert "1 world with a civilization" in info and "1 world with sentient life" in info and "1 world with life" in info, info


def test_galaxy_map_civilizations_can_be_chosen_by_technology_band(page, map_site):
    """MAP.344: the Menu has a Civilization technology group with a marker for each technology band, off until chosen; a band shows only the
    systems whose most advanced species is in it, a system with no species drawn yet is only in the Civilizations marker, and the card
    names the band."""
    _open_star_field(page, map_site)
    page.wait_for_timeout(1000)
    _open_menu(page)
    group = page.locator("#galaxymap3d-kinds")
    assert "Civilization technology" in group.inner_text()
    bands = group.locator('button[data-layer^="marker-tech-"]')
    assert bands.count() == 8
    assert all(button.get_attribute("aria-pressed") == "false" for button in bands.all())
    assert [button.inner_text().strip() for button in bands.all()][:3] == ["Primal civilizations", "Agrarian civilizations", "Industrial civilizations"]
    _poi_layers(page)    # no other reason on
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-black-hole" for e in entries), seconds=2)
    assert not [e for e in found if e["kind"].startswith("marker-tech-")], "no band is marked until one is chosen"
    group.locator('button[data-layer="marker-tech-industrial"]').click()
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-tech-industrial" for e in entries))
    assert [e["name"] for e in found if e["kind"].startswith("marker-tech-")] == ["Fixture Rest"], found
    assert "show=marker-tech-industrial" in page.url, page.url
    group.locator('button[data-layer="marker-tech-gravitic"]').click()
    found = _markers_until(page, lambda entries: {e["name"] for e in entries if e["kind"].startswith("marker-tech-")} == {"Fixture Rest", "Fixture Dawn"})
    assert {(e["kind"], e["name"]) for e in found if e["kind"].startswith("marker-tech-")} == {
        ("marker-tech-industrial", "Fixture Rest"), ("marker-tech-gravitic", "Fixture Dawn")}, found
    page.locator('#galaxymap3d-kinds button[data-layer="marker-civilization"]').click()
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-civilization" for e in entries))
    assert {(e["kind"], e["name"]) for e in found if e["kind"] in ("marker-civilization",) or e["kind"].startswith("marker-tech-")} == {
        ("marker-tech-industrial", "Fixture Rest"), ("marker-tech-gravitic", "Fixture Dawn"), ("marker-civilization", "Fixture Ember")}, found
    page.locator("#galaxymap3d-menu > summary").click()
    page.wait_for_timeout(300)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    dawn = next(s for s in LIFE_FIELD if s["name"] == "Fixture Dawn")
    spot = _screen_of(page, (dawn["x"], dawn["y"], dawn["z"]))
    assert spot
    page.mouse.move(*spot)
    page.mouse.move(spot[0] + 1, spot[1])
    page.mouse.click(*spot)
    info = page.locator("#galaxymap3d-info").inner_text()
    assert "2 worlds with a civilization, most advanced Gravitic" in info, info


def test_the_default_map_marks_a_green_and_a_yellow_world_with_the_reason_listed_and_a_red_one_not(page, map_site):
    """MAP.303 acceptance: with no Menu choice made, the map marks the star of a Green world and of a Yellow one (the Yellow one's
    system is also inhabited, so its mark is the inhabited symbol), each card lists the habitable reason with its best tier, and the
    system whose only world is Red (below the Yellow floor) has no mark."""
    _open_star_field(page, map_site)
    page.wait_for_timeout(1000)
    found = _markers_until(page, lambda entries: {"Fixture Eden", "Fixture Rest"} <= {e["name"] for e in entries})
    marks = {e["name"]: e["kind"] for e in found if e["name"]}
    assert marks["Fixture Eden"] == "marker-habitable" and marks["Fixture Rest"] == "marker-inhabited", marks
    assert "Fixture Rock" not in marks, "a system with only a Red world has no mark"
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    for name, wanted in (("Fixture Eden", "2 habitable worlds, best Green"), ("Fixture Rest", "1 habitable world, best Yellow")):
        mark = next(e for e in found if e["name"] == name)
        spot = _screen_of(page, (mark["x"], mark["y"], mark["z"]))
        assert spot
        page.mouse.move(*spot)
        page.mouse.move(spot[0] + 1, spot[1])
        page.mouse.click(*spot)
        assert page.locator("#galaxymap3d-info h3").inner_text() == name
        _expand_card(page)
        assert wanted in page.locator("#galaxymap3d-info").inner_text()


def test_the_points_of_interest_menu_group_shows_pulls_forward_or_hides_each_kind_and_remembers_it(page, map_site):
    """MAP.306: the Menu's Points of interest group has a row for each reason, whose second button pulls it forward (drawn over the
    other marks, named first); a crowd takes its strongest member's state; hiding a kind takes its marks away; one row sets every
    kind at once; the address (`front=`, `hide=`) and the browser keep the choice."""
    _open_star_field(page, map_site)
    page.wait_for_timeout(1000)
    _open_menu(page)
    kinds = page.locator("#galaxymap3d-kinds")
    assert "Points of interest" in kinds.inner_text()
    forward = {layer: kinds.locator(f'button[data-forward="{layer}"]') for layer in POI_LAYERS}
    assert all(button.count() == 1 for button in forward.values())
    assert all(forward[layer].get_attribute("aria-pressed") == "false" for layer in POI_LAYERS)
    assert forward["marker-uninhabitable"].is_disabled(), "a hidden row cannot be pulled forward by its own button"
    # Pull the gas giants forward: their marks are drawn in front and named first.
    forward["marker-gas-giant"].click()
    assert "front=marker-gas-giant" in page.url
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-gas-giant" and e["pulled"] for e in entries))
    assert {e["name"] for e in found if e["kind"] == "marker-gas-giant" and e["pulled"]} == {"Fixture Jove", "Fixture Eden"}, found
    assert not [e for e in found if e["kind"] != "marker-gas-giant" and e["pulled"]]
    page.wait_for_function(f"() => ({LABELS})().some((l) => l.kind === 'world' && l.text === 'Fixture Jove')", timeout=20000)
    # Hide the habitable worlds: none of their marks is drawn.
    kinds.locator('button[data-layer="marker-habitable"]').click()
    found = _markers_until(page, lambda entries: not any(e["kind"] == "marker-habitable" for e in entries))
    assert not [e for e in found if e["kind"] == "marker-habitable"]
    assert "hide=marker-habitable" in page.url
    # A reload keeps the choices, from the address ...
    url = page.url
    page.reload()
    page.wait_for_selector("#galaxymap3d-steps [data-steps-panel] li", state="attached")
    _settle(page)
    _open_menu(page)
    assert page.url == url
    assert page.locator('#galaxymap3d-kinds button[data-forward="marker-gas-giant"]').get_attribute("aria-pressed") == "true"
    assert page.locator('#galaxymap3d-kinds button[data-layer="marker-habitable"]').get_attribute("aria-pressed") == "false"
    # ... and a fresh address, from the browser.
    _open(page, f"{map_site}/galaxy", GALAXY_CANVAS)
    page.wait_for_selector("#galaxymap3d-steps [data-steps-panel] li", state="attached")
    _open_menu(page)
    assert page.locator('#galaxymap3d-kinds button[data-forward="marker-gas-giant"]').get_attribute("aria-pressed") == "true"
    assert page.locator('#galaxymap3d-kinds button[data-layer="marker-habitable"]').get_attribute("aria-pressed") == "false"
    # The one row sets every kind: all hidden, all shown.
    page.locator('#galaxymap3d-kinds button[data-bulk="hidden"][data-group="Points of interest"]').click()
    assert all(page.locator(f'#galaxymap3d-kinds button[data-layer="{layer}"]').get_attribute("aria-pressed") == "false" for layer in POI_LAYERS)
    page.locator('#galaxymap3d-kinds button[data-bulk="shown"][data-group="Points of interest"]').click()
    assert all(page.locator(f'#galaxymap3d-kinds button[data-layer="{layer}"]').get_attribute("aria-pressed") == "true" for layer in POI_LAYERS)
    assert all(forward[layer].get_attribute("aria-pressed") == "false" for layer in POI_LAYERS)
    page.locator('#galaxymap3d-kinds button[data-bulk="forward"][data-group="Points of interest"]').click()
    assert all(forward[layer].get_attribute("aria-pressed") == "true" for layer in POI_LAYERS)


def test_a_crowd_of_marks_takes_the_state_of_its_strongest_kind(page, map_site):
    """MAP.306: from far out the fixture's points of interest merge into one mark; with the gas giants pulled forward the crowd is a
    gas giant's symbol, drawn in front, and with none pulled it is the best reason's symbol, not in front."""
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_timeout(1000)
    crowd = lambda found: next((e for e in found if e["count"]), None)
    found = _markers_until(page, lambda entries: crowd(entries) is not None)
    assert crowd(found)["kind"] == "marker-inhabited" and not crowd(found)["pulled"]
    _open_menu(page)
    page.locator('#galaxymap3d-kinds button[data-forward="marker-gas-giant"]').click()
    found = _markers_until(page, lambda entries: crowd(entries) and crowd(entries)["pulled"])
    assert crowd(found)["kind"] == "marker-gas-giant" and crowd(found)["pulled"], found


def test_galaxy_map_marks_each_system_with_its_best_worlds_symbol(page, map_site):
    """MAP.321: with every world reason on a system has one mark, the symbol of its best world (an inhabited world before a habitable world before a barren
    one before a gas giant before a belt), and each reason has a button of its own (the barren and belt ones off until asked for)."""
    _open_star_field(page, map_site)
    page.wait_for_timeout(1000)
    ids = ("marker-inhabited", "marker-habitable", "marker-uninhabitable", "marker-gas-giant", "marker-asteroid-field")
    _open_menu(page)
    for layer in ids:
        button = page.locator(f'#galaxymap3d-kinds button[data-layer="{layer}"]')
        assert button.count() == 1 and button.get_attribute("aria-pressed") == ("true" if layer in POI_DEFAULT_ON else "false"), layer
    _poi_layers(page, *ids)
    found = _markers_until(page, lambda entries: {e["kind"] for e in entries} >= set(ids))
    marks = {e["name"]: e["kind"] for e in found if e["kind"] in ids}
    assert sorted(e["name"] for e in found if e["kind"] in ids) == sorted(marks), "no system is marked twice"
    assert marks == {
        "Fixture Eden": "marker-habitable", "Fixture Rest": "marker-inhabited", "Fixture Rock": "marker-uninhabitable",
        "Fixture Jove": "marker-gas-giant", "Fixture Belt": "marker-asteroid-field"}, marks
    # Turn the habitable reason off: Eden's best remaining world is its barren one, Rest keeps its inhabited one.
    page.locator('#galaxymap3d-kinds button[data-layer="marker-habitable"]').click()
    found = _markers_until(page, lambda entries: not any(e["kind"] == "marker-habitable" for e in entries))
    marks = {e["name"]: e["kind"] for e in found if e["kind"] in ids}
    assert marks == {
        "Fixture Eden": "marker-uninhabitable", "Fixture Rest": "marker-inhabited", "Fixture Rock": "marker-uninhabitable",
        "Fixture Jove": "marker-gas-giant", "Fixture Belt": "marker-asteroid-field"}, marks


def test_a_name_is_only_offered_for_a_habitable_world_marker_that_is_drawn_and_a_sparse_view_names_them_all(page, map_site):
    """MAP.317: the habitable-world markers are the inhabited slice of the label budget: a sparse view names every one it
    draws that has room (the fixture's two lie close, so the one with more worlds keeps its name), and none while the marker
    layer is off or the World names switch is."""
    _open_star_field(page, map_site)
    page.wait_for_timeout(1000)
    systems = {"Fixture Eden"}
    _poi_layers(page)
    page.wait_for_function(f"() => !({LABELS})().some((label) => label.kind === 'world')")
    button = page.locator('#galaxymap3d-kinds button[data-layer="marker-habitable"]')
    button.click()
    _markers_until(page, lambda found: any(e["kind"] == "marker-habitable" for e in found))
    page.wait_for_function(f"([names]) => {{ const got = new Set(({LABELS})().filter((l) => l.kind === 'world').map((l) => l.text)); return names.every((n) => got.has(n)); }}",
                           arg=[sorted(systems)], timeout=20000)
    # The switch for the names takes them away while the markers stay.
    page.locator('#galaxymap3d-kinds [data-label="world"]').click()
    page.wait_for_function(f"() => !({LABELS})().some((label) => label.kind === 'world')")
    assert any(e["kind"] == "marker-habitable" for e in page.evaluate(MARKERS))
    page.locator('#galaxymap3d-kinds [data-label="world"]').click()
    page.wait_for_function(f"() => ({LABELS})().some((label) => label.kind === 'world')")
    # No marker drawn, no name offered.
    button.click()
    page.wait_for_function(f"() => !({LABELS})().some((label) => label.kind === 'world')")


def test_crowded_marks_merge_into_one_with_a_count_whose_card_lists_every_system_and_zooming_in_splits_it(page, map_site):
    """MAP.305: from far out the fixture's five points of interest lie too close to tell apart, so they are drawn as fewer marks, a
    crowd showing its count; a click on a crowd opens a card that lists every system in it with its reasons, and zooming in on it
    spreads the marks until each system has its own."""
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_timeout(1000)
    ids = ("marker-inhabited", "marker-habitable", "marker-uninhabitable", "marker-gas-giant", "marker-asteroid-field")
    systems = {"Fixture Eden", "Fixture Rest", "Fixture Rock", "Fixture Jove", "Fixture Belt"}
    _poi_layers(page, *ids)

    def marks(found):
        return [e for e in found if e["kind"] in ids]

    def represented(found):
        return sorted(name for e in marks(found) for name in (e["members"] if e["count"] else [e["name"]]))

    found = _markers_until(page, lambda entries: len(represented(entries)) == 5)
    assert represented(found) == sorted(systems), found
    crowds = [e for e in marks(found) if e["count"]]
    assert crowds and len(marks(found)) < 5, "the close systems are drawn as fewer marks"
    crowd = max(crowds, key=lambda e: e["count"])
    page.locator("#galaxymap3d-menu > summary").click()
    page.wait_for_timeout(300)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    spot = _screen_of(page, (crowd["x"], crowd["y"], crowd["z"]))
    assert spot
    page.mouse.move(*spot)
    page.mouse.move(spot[0] + 1, spot[1])
    page.wait_for_timeout(200)
    assert f"{crowd['count']} points of interest" in page.locator("#galaxymap3d-tooltip").inner_text()
    page.mouse.click(*spot)
    assert page.locator("#galaxymap3d-info h3").inner_text() == f"{crowd['count']} points of interest"
    _expand_card(page)    # the members are behind Expand
    info = page.locator("#galaxymap3d-info").inner_text()
    assert all(name in info for name in crowd["members"]), info
    page.locator("#galaxymap3d-info button", has_text="Go to them").click()
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=30000)
    _settle(page)
    after = _markers_until(page, lambda entries: not any(e["count"] and e["count"] >= crowd["count"] and e["members"] == crowd["members"] for e in entries))
    assert represented(after) == sorted(systems)
    assert not any(e["members"] == crowd["members"] for e in marks(after)), "zooming in splits the crowd"


def test_an_inhabited_world_marker_is_offered_a_name_for_as_long_as_it_is_drawn(page, map_site):
    """MAP.321 d: an inhabited token is a candidate for a name of the inhabited slice (the budget's order of the two is the unit test's):
    with only the inhabited layer on the fixture's inhabited world is named and the habitable one that lies close is not, with both on
    both are, and with both off no world name is offered for a token that is not drawn."""
    _open_star_field(page, map_site)
    page.wait_for_timeout(1000)
    worlds = f"() => ({LABELS})().filter((l) => l.kind === 'world').map((l) => l.text).sort()"
    _poi_layers(page)
    page.locator('#galaxymap3d-kinds button[data-layer="marker-inhabited"]').click()
    _markers_until(page, lambda found: any(e["kind"] == "marker-inhabited" for e in found))
    page.wait_for_function(f"() => ({LABELS})().some((l) => l.kind === 'world' && l.text === 'Fixture Rest')", timeout=20000)
    assert page.evaluate(worlds) == ["Fixture Rest"]
    page.locator('#galaxymap3d-kinds button[data-layer="marker-habitable"]').click()
    page.wait_for_function(f"() => ({LABELS})().some((l) => l.kind === 'world' && l.text === 'Fixture Eden')", timeout=20000)
    assert page.evaluate(worlds) == ["Fixture Eden", "Fixture Rest"]
    page.locator('#galaxymap3d-kinds button[data-layer="marker-habitable"]').click()
    page.locator('#galaxymap3d-kinds button[data-layer="marker-inhabited"]').click()
    page.wait_for_function(f"() => !({LABELS})().some((l) => l.kind === 'world')", timeout=20000)


SECTOR_FRAME = "() => { const f = document.querySelector('#galaxymap3d-canvas').__mapSectorFrame; return f && f.state(); }"
SECTOR_VIEW = "() => document.querySelector('#galaxymap3d-canvas').galaxySectorView()"
SECTOR_FRAME_LINES = "() => document.querySelector('#galaxymap3d-canvas').__mapSectorFrame.positions()"


def test_the_opened_sector_is_framed_and_dimmed_around_with_a_labelled_icon_on_each_wall_that_flies_to_the_neighbour(page, map_site):
    """MAP.327: with a sector opened in place the space outside it is dimmed (no glow), a wireframe round its cell stays
    where the cell is, and each wall facing another sector has one labelled icon at its centre whose click arrives there."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, map_site, wall_icons=True)
    assert page.evaluate(SECTOR_VIEW)["outsideDim"] == 1 and page.evaluate(SECTOR_FRAME)["framed"] is False
    _open_galaxy(page, map_site, OPEN_PRIME, wall_icons=True)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    page.wait_for_function(f"() => {{ const f = ({SECTOR_FRAME})(); return f && f.framed && f.icons.some((i) => i.shown); }}", timeout=20000)
    _settle(page)
    view = page.evaluate(SECTOR_VIEW)
    assert view == {"outsideDim": 0.3, "glowVisible": False, "opened": True}, view
    state = page.evaluate(SECTOR_FRAME)
    names = [icon["name"] for icon in state["icons"]]
    assert names and len(set(names)) == len(names) and all(name.startswith("Sector ") for name in names), names
    # Each shown icon stands at the centre of its wall, the icon's own button what the page shows there.
    shown = [icon for icon in state["icons"] if icon["shown"]]
    assert shown
    for icon in shown:
        spot = _screen_of(page, tuple(icon["wall"]))
        assert spot
        on = page.evaluate("([x, y]) => { const e = document.elementFromPoint(x, y); const b = e && e.closest('.map-wall-icon'); return b && b.dataset.dir; }", [spot[0], spot[1]])
        assert on == icon["dir"], (icon["dir"], on)
    # The frame is fixed: the pointer moving over the map changes nothing of it.
    before = page.evaluate(SECTOR_FRAME_LINES)
    page.mouse.move(300, 300)
    page.mouse.move(420, 360)
    page.wait_for_timeout(300)
    assert page.evaluate(SECTOR_FRAME_LINES) == before
    # A click on an icon arrives at that neighbour: the address names it, the frame is round it.
    target = shown[0]
    designation = target["name"].replace("Sector ", "")
    page.locator(f'.map-wall-icon[data-dir="{target["dir"]}"]').click()
    page.wait_for_function("([d]) => location.search.includes('sector=' + d)", arg=[designation], timeout=20000)
    _settle(page)
    assert page.evaluate(SECTOR_FRAME_LINES) != before or not page.evaluate(SECTOR_VIEW)["opened"]


OPACITIES = "() => document.querySelector('#galaxymap3d-canvas').__mapMarkers.opacities()"


def test_markers_outside_an_opened_sector_dim_as_the_stars_do(page, holes_site):
    """MAP.341: with a sector opened in place the stars outside it dim to 0.3 (MAP.327), and so does every marker outside it, black holes
    and the rest; with none opened the markers are whole, and closing the sector gives them back."""
    _pin_window(page)
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, holes_site)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapMarkers.entries().length > 0", timeout=20000)
    whole = page.evaluate(OPACITIES)
    assert whole and all(mark["opacity"] == 1 for mark in whole), whole
    _open_galaxy(page, holes_site, OPEN_PRIME)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    _settle(page)
    assert page.evaluate(SECTOR_VIEW)["outsideDim"] == 0.3
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapMarkers.entries().length > 0", timeout=20000)
    dimmed = page.evaluate(OPACITIES)
    assert dimmed and all(mark["opacity"] == 0.3 for mark in dimmed), dimmed
    assert any(mark["kind"] == "marker-black-hole" for mark in dimmed), dimmed


def _sector_in_address(page):
    return parse_qs(urlparse(page.url).query).get("sector", [None])[0]


CENTRE_NOW = """() => { const c = document.querySelector("#galaxymap3d-canvas").galaxyStage().win.centre;
    return [c.m, c.shell, c.band, c.wedge].join("."); }"""


def _sector_centre(page):
    """The sector the view is on, `m.shell.band.wedge` (m is 1), as the window it sits in names it."""
    return page.evaluate(CENTRE_NOW)


def _wait_off(page, was):
    page.wait_for_function(f"([was]) => ({CENTRE_NOW})() !== was", arg=[was], timeout=20000)
    _settle(page)
    return _sector_centre(page)


def test_the_step_buttons_and_keys_go_from_an_opened_sector_to_the_sector_beside_it_and_back(page, map_site):
    """MAP.121: with a sector opened in place the Menu's step buttons and the keys take the view to the sector across
    the wall in that direction (the step the Galaxy Map's windows use and MAP.327's wall icons call), in one camera move
    with the address naming it; a button is off where the wall faces nothing, and a sideways step back lands where it
    left."""
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    _settle(page)
    start = _sector_centre(page)
    assert start.startswith("1."), start
    assert _sector_in_address(page) == PRIME_CODE

    def enabled():
        return page.eval_on_selector_all(
            "#galaxymap3d-controls [data-action=step]",
            "els => Object.fromEntries(els.map(e => [e.dataset.step, !e.disabled]))")

    ways = enabled()
    assert set(ways) == {"north", "south", "east", "west", "in", "out"}
    assert ways["out"] and any(ways.values()), ways

    # Out, with the Menu's button: the sector across the outer wall, the address naming it.
    page.click("#galaxymap3d-menu > summary")
    page.click('#galaxymap3d-controls [data-step="out"]')
    outer = _wait_off(page, start)
    assert outer.split(".")[1] == str(int(start.split(".")[1]) + 1), (start, outer)
    assert _sector_in_address(page) not in (None, PRIME_CODE)

    # In again, with the key: the sector across the inner wall (not always the one it came from: the rings differ in how
    # many sectors they hold, so a step out and back need not land on the same cell).
    page.locator(GALAXY_CANVAS).focus()
    page.keyboard.press("[")
    here = _wait_off(page, outer)
    assert here.split(".")[1] == start.split(".")[1], (start, outer, here)

    # Sideways with the arrow keys: a step east then west lands on the same sector.
    if enabled()["east"]:
        page.keyboard.press("Alt+Shift+ArrowRight")
        east = _wait_off(page, here)
        assert east.split(".")[1] == here.split(".")[1], (here, east)
        page.keyboard.press("Alt+Shift+ArrowLeft")
        assert _wait_off(page, east) == here


def test_galaxy_map_point_phenomena_hover_and_offer_nav_links(page, map_site):
    """MAP.65: a black hole, neutron star or quasar shows a tooltip with its
    name under the pointer, and its panel offers Start Here, End Here
    and a ☆, as the Sector Map's does."""
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_timeout(1000)
    names = {point["name"]: point for point in POINT_FIELD}
    tip = page.locator("#galaxymap3d-tooltip")
    close = page.locator("#galaxymap3d-info .map-info-close")
    if close.is_visible():   # MAP.256: the panel over the map covers the points until it is closed
        close.click()
    for x, y in _point_colored_spots(page, GALAXY_CANVAS):
        page.mouse.move(x, y)
        text = tip.inner_text() if tip.is_visible() else ""
        name = next((n for n in names if text.startswith(n)), None)
        if not name:
            continue
        point = names[name]
        page.mouse.click(x, y)
        info = page.locator("#galaxymap3d-info")
        assert info.locator("h3").inner_text() == name
        endpoint = f'{point["type"]}:{point["id"]}'
        _expand_card(page)    # the course buttons are behind Expand
        assert [b.inner_text() for b in info.locator(".map-info-actions button:not(.map-info-bookmark):not(.map-info-close):not(.map-info-expand)").all()] == ["Start Here", "End Here"]
        assert info.locator(".map-info-bookmark").inner_text() == "☆"
        info.locator("button", has_text="Start Here").click()
        assert parse_qs(urlparse(page.url).query)["from"] == [endpoint], page.url
        return
    pytest.fail("no point phenomenon showed a tooltip under the pointer")


def test_galaxy_map_pick_mode_clicks_down_through_every_stage(page, map_site):
    """NAV.46: choosing a NAV end, the arcs, slabs and blocks are clicked to
    zoom in, stage after stage as outside pick mode, down to a sector."""
    _open_galaxy(page, map_site, "?pick=to&from=system:701")
    seen = [_crumbs(page)]
    for _ in range(8):
        text = _click_choice(page, GENERATED_CHOICE)
        if re.match(r"Sector .*, charted$", text):
            break
        assert _crumbs(page) != seen[-1], (text, _crumbs(page))
        seen.append(_crumbs(page))
    assert len(seen) >= 3, seen


def test_galaxy_map_stars_are_drawn_relative_to_the_eye(page, map_site):
    """MAP.102, MAP.267: each star is a high and a low float32 and the shader subtracts the camera's own pair, so
    zoomed in the offset the GPU works with is small and exact (to a float32 of the offset, not of the position)
    rather than hundreds or thousands of parsecs out."""
    _open_window(page, map_site, 3, 250, 0, 0)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyStarFrame().count > 0")
    drawn = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyStarFrame()")
    # Measured from the camera, no star is as far out as the galaxy's edge, and the pair puts each where it is.
    assert drawn["far"] < 20000, drawn
    assert drawn["error"] < 1e-3 * max(drawn["far"], 1.0), drawn


def test_galaxy_map_stars_are_filled_node_by_node_and_never_rebuilt(page, map_site):
    """MAP.267: a node landing writes its own stars at the end of the one star buffer and a node leaving is cut out by
    moving the tail into its place, so across a run of zoom steps the stars written are a small share of what a rebuild
    at every step would have written, and a step that changes no node writes nothing."""
    _open_galaxy(page, map_site)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyStarFrame().count > 0")
    fill = "() => document.querySelector('#galaxymap3d-canvas').galaxyStarFill()"
    start = page.evaluate(fill)
    assert start["capacity"] >= start["count"] > 0 and start["last"] is not None, start
    rebuilt = 0
    box = page.locator(GALAXY_CANVAS).bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    previous = start["syncs"]
    for _ in range(14):
        page.mouse.wheel(0, -300)
        page.wait_for_timeout(250)
        now = page.evaluate(fill)
        if now["syncs"] > previous:
            rebuilt += now["count"]
            previous = now["syncs"]
    end = page.evaluate(fill)
    assert end["syncs"] > start["syncs"], "the zoom drew new nodes"
    written = end["written"] - start["written"]
    assert written < 0.75 * rebuilt, (written, rebuilt, end)
    # A sync that keeps every node holds stars where they were: nothing written or moved.
    page.wait_for_timeout(600)
    settled = page.evaluate(fill)
    page.wait_for_timeout(600)
    again = page.evaluate(fill)
    assert (again["written"], again["moved"]) == (settled["written"], settled["moved"])


def test_galaxy_map_stars_fade_in_with_the_zoom(page, map_site):
    """MAP.153, MAP.261: each star comes in with the weight of its node, so the opacity on screen changes by a
    small share over a zoom step instead of the whole set arriving at once."""
    problems = []
    page.on("console", lambda message: problems.append(message.text) if message.type == "error" else None)
    _open_galaxy(page, map_site)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyStarFrame().count > 0")
    page.wait_for_timeout(500)
    assert not [text for text in problems if "shader" in text.lower() or "WebGLProgram" in text], problems
    opacity = "(r) => document.querySelector('#galaxymap3d-canvas').galaxyStarOpacity(r)"
    here = page.evaluate(opacity, 1.0)
    camera = here["camera"]
    assert here["count"] > 0 and here["camera"] > 0, here
    # Out to 1.4 times the camera's distance: further out the cut would hand the stars over to coarser nodes that
    # this view has not fetched, so the loaded nodes alone fade out there.
    steps = [camera * 2 ** (k / 8) for k in range(-8, 5)]
    sums = [page.evaluate(opacity, r)["sum"] for r in steps]
    assert all(0 <= s <= here["count"] + 1e-6 for s in sums), sums
    # A star listed in a node and in the nodes under it keeps its whole share as the weights hand over (the unit test
    # `starfade.test.mjs` bounds the same for a field whose parents drop stars their children list): `count` is the
    # entries (a star once for every node that lists it), and the opacity on screen is the stars, whatever the step.
    assert sums[0] > 0 and all(abs(s - sums[0]) <= 0.02 * sums[0] for s in sums), f"{here['count']} entries: {sums}"


def test_galaxy_map_near_field_thins_what_stands_between_the_camera_and_the_focus(page, map_site):
    """MAP.149: the shaders get the camera-to-focus distance; a point on the
    line of sight, nearer than the focus, is faded to the tube's minimum
    and one close to the camera to nothing, while the focus itself and a
    point off to the side stay whole."""
    problems = []
    page.on("console", lambda message: problems.append(message.text) if message.type == "error" else None)
    _open_galaxy(page, map_site)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyStarFrame().count > 0")
    page.wait_for_timeout(500)
    assert not [text for text in problems if "shader" in text.lower() or "WebGLProgram" in text], problems
    read = "(p) => document.querySelector('#galaxymap3d-canvas').galaxyNearField(...p)"
    near = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyNearField()")
    camera, target, distance = near["camera"], near["target"], near["distance"]
    assert distance > 0, near
    assert abs(sum((a - b) ** 2 for a, b in zip(camera, target)) ** 0.5 - distance) < 1e-6 * distance, near

    def along(t):
        return [c + t * (f - c) for c, f in zip(camera, target)]

    assert page.evaluate(read, along(1.0))["share"] > 0.999
    assert page.evaluate(read, along(1.5))["share"] > 0.999
    assert page.evaluate(read, along(0.05))["share"] < 1e-6
    assert abs(page.evaluate(read, along(0.6))["share"] - 0.25) < 0.02


def test_galaxy_map_stars_are_never_thinned_by_the_near_field(page, map_site):
    """MAP.187: the near-field window thins the blocks only. The star shader holds no near-field
    term, so a star in front of the camera or on the line to the focus is drawn whole at any zoom,
    while a point on that line is still thinned for the blocks."""
    problems = []
    page.on("console", lambda message: problems.append(message.text) if message.type == "error" else None)
    _open_galaxy(page, map_site)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyStarFrame().count > 0")
    page.wait_for_timeout(500)
    assert not [text for text in problems if "shader" in text.lower() or "WebGLProgram" in text], problems
    shader = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyStarShader()")
    assert "nearField" not in shader["vertex"]
    assert "nearD" not in shader["uniforms"] and "nearR" not in shader["uniforms"]
    near = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyNearField()")
    midway = [c + 0.6 * (f - c) for c, f in zip(near["camera"], near["target"])]
    share = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyNearField(...p)", midway)["share"]
    assert share < 1, "the blocks are still thinned on the line to the focus"


def test_galaxy_map_stars_follow_the_apparent_magnitude_law(page, map_site):
    """MAP.148: the stars are drawn against a limiting magnitude the page
    works out from the stars on screen; the shaders compile, the limit
    settles on a number, and the stars the map holds show (the fixture has
    far fewer than the target, so the limit sits past the dimmest star
    and none is thinned away)."""
    problems = []
    page.on("console", lambda message: problems.append(message.text) if message.type == "error" else None)
    _open_galaxy(page, map_site)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyStarFrame().count > 0")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyStarMagnitude().limit < 999")
    page.wait_for_timeout(800)
    assert not [text for text in problems if "shader" in text.lower() or "WebGLProgram" in text], problems
    info = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyStarMagnitude()")
    assert info["target"] in (8000, 20000), info["target"]
    assert info["count"] > 0 and info["onScreen"] <= info["count"], info
    assert abs(info["limit"] - info["goal"]) < 0.05, "the limit has eased to its goal"
    assert info["limit"] > -5, info["limit"]
    shown = [share for magnitude, share in info["stars"] if share > 0]
    assert shown, "some stars show"
    assert all(0 <= share <= 1 for _, share in info["stars"])
    # Fewer stars than the target on screen, so the count thins none of them.
    assert info["onScreen"] < info["target"]


def _select_the_nebula(page):
    nebula = next(p for p in PHENOMENA if p["type"] == "nebula")
    # The screen-reader list holds a button for every object.
    page.locator("#galaxymap3d-canvas ~ .starmap-sr-list button, .starmap-sr-list button",
                 has_text=nebula["name"]).first.dispatch_event("click")
    assert _info_title(page) == nebula["name"]
    return nebula


def test_sector_map_escape_clears_the_selection(page, map_site):
    """MAP.113: Escape clears what is selected, ahead of anything else."""
    _open_sector(page, map_site)
    resting = _info_title(page)
    nebula = _select_the_nebula(page)
    assert resting != nebula["name"]
    page.locator(SECTOR_CANVAS).focus()
    page.keyboard.press("Escape")
    assert _info_title(page) == resting
    assert page.url.split("?")[0].endswith(f"/sector/{SECTOR_ID}"), page.url


@pytest.mark.names_hidden
@pytest.mark.viewport(1000, 900)  # the contents list (UX.215) takes height from the map in 1000x800, where no spot on it names the nebula
def test_sector_map_click_on_the_selected_nebula_clears_it(page, map_site):
    """MAP.113: a nebula can be unselected by clicking it again, which is all
    there is to click when it covers the sector."""
    _open_sector(page, map_site)
    resting = _info_title(page)
    nebula = next(p for p in PHENOMENA if p["type"] == "nebula")
    found = _click_steady(page, re.escape(nebula["name"]))
    assert _info_title(page) == nebula["name"]
    page.wait_for_timeout(DOUBLE_CLICK_GAP_MS)   # a second click this soon would be a double click (MAP.226)
    page.mouse.click(found[0], found[1])
    assert _info_title(page) == resting
    page.wait_for_timeout(DOUBLE_CLICK_GAP_MS)
    page.mouse.click(found[0], found[1])
    assert _info_title(page) == nebula["name"], "and a click selects it again"


def _url_settled(page, quiet_ms=900):
    """The page's URL once it has stopped changing. The map writes the camera into the address 350 ms after the
    camera last moved, in its frame loop, so on a loaded machine the scale bar can hold still (`_settle`) before
    the address names where the camera ended (TEST.126)."""
    last, held = None, 0
    for _ in range(200):
        if page.url == last:
            held += 100
            if held >= quiet_ms:
                return last
        else:
            last, held = page.url, 0
        page.wait_for_timeout(100)
    return page.url


def _crumbs_on_a_fresh_load(page):
    """The steps a new page shows for this page's URL."""
    _url_settled(page)
    other = page.context.new_page()
    seen = []
    other.on("pageerror", lambda error: seen.append(f"page error: {error}"))
    other.on("console", lambda message: seen.append(f"console {message.type}: {message.text}")
             if message.type in ("error", "warning") else None)
    try:
        other.goto(page.url)
        try:
            other.wait_for_selector("#galaxymap3d-steps [data-steps-panel] li", state="attached")
        except sync_api.TimeoutError:
            # TEST.287: a fresh load of this address showed no steps in 120 s on a loaded machine, once in two runs.
            body = other.evaluate("() => (document.body ? document.body.innerText : '').slice(0, 600)")
            raise AssertionError(
                f"no steps on a fresh load of {page.url}: now at {other.url}, titled {other.title()!r}; "
                f"the page said {body!r}; it logged {seen[-8:]}") from None
        _settle(other)
        return _crumbs(other)
    finally:
        other.close()


def test_galaxy_breadcrumb_always_matches_the_view_its_url_names(page, map_site):
    """MAP.106: walking down to a sector, out again, Home, Back and Forward,
    the steps shown are the ones a fresh load of the same URL shows."""
    _open_galaxy(page, map_site)
    checked = []

    def in_sync(label):
        _settle(page)
        _url_settled(page)
        live = _crumbs(page)
        assert live == _crumbs_on_a_fresh_load(page), (label, page.url, live)
        checked.append(label)

    _walk_down(page)
    in_sync("in the sector")
    assert _crumbs(page)[-1].startswith("Sector"), _crumbs(page)
    canvas = page.locator("#galaxymap3d-canvas")
    for key in ("Escape", "Home"):
        canvas.focus()
        page.keyboard.press(key)
        in_sync(key)
    for action in ("back", "forward"):
        button = page.locator(f'#galaxymap3d-controls [data-action="{action}"]')
        if button.is_enabled():
            button.click()
            in_sync(action)
    assert len(checked) >= 4, checked


# --- UX.50, UX.57, UX.58, UX.60, UX.62: the map pages' controls and layout ---------------

def test_the_map_help_is_a_menu_item_that_opens_a_dialog(page, map_site):
    """UX.50: the how-to text is behind Map help in the Menu, on the Galaxy
    Map and on the Sector Map; nothing of it sits under the map."""
    for opener, name in ((lambda: _open_galaxy(page, map_site), "Galaxy"),
                         (lambda: _open_sector(page, map_site), "Sector")):
        opener()
        assert page.locator("#map > .panel-header .hint").count() == 0, name
        dialog = page.locator("#galaxymap3d-help")
        assert not dialog.evaluate("d => d.open"), name
        # The dialog is a Shoelace element; before it is upgraded it has no show() and the item does
        # nothing, which a loaded machine reaches before the click (TEST.126).
        page.wait_for_function("() => typeof document.querySelector('#galaxymap3d-help').show === 'function'")
        _open_menu(page)
        page.locator('#galaxymap3d-controls [data-action="map-help"]').click()
        page.wait_for_function("() => document.querySelector('#galaxymap3d-help').open")
        assert "Drag to turn" in dialog.inner_text() and "Re-center" in dialog.inner_text(), name
        assert not page.locator("#galaxymap3d-menu").evaluate("m => m.open"), "the Menu closes behind it"
        dialog.locator("[data-dialog-close]").click()
        page.wait_for_function("() => !document.querySelector('#galaxymap3d-help').open")


def test_the_galaxy_toolbar_says_whole_galaxy_and_the_menu_re_center(page, map_site):
    """UX.57: no two controls called Reset."""
    _open_galaxy(page, map_site)
    assert page.locator('#galaxymap3d-controls [data-action="reset"]').inner_text() == "Whole galaxy"
    assert page.locator('#galaxymap3d-controls [data-action="reset-view"]').text_content() == "Re-center"


@pytest.mark.viewport(1280, 900)
def test_the_map_card_goes_wide_but_the_title_keeps_the_page_edge(page, map_site):
    """UX.62: the Galaxy page's title and crumb sit where every other page's do. A page with a contents menu (UX.215) is wider and
    puts its title beside the menu, so the comparison is with a page that has none (TEST.445)."""
    _open_galaxy(page, map_site)
    galaxy = page.evaluate("() => document.querySelector('h1').getBoundingClientRect().left")
    card = page.evaluate("() => document.querySelector('#map').getBoundingClientRect().width")
    assert not page.locator("nav.page-toc").count()
    page.goto(f"{map_site}/", wait_until="load")
    assert not page.locator("nav.page-toc").count(), "the page compared with has a contents menu"
    other = page.evaluate("() => document.querySelector('h1').getBoundingClientRect().left")
    assert abs(galaxy - other) < 1, (galaxy, other)
    assert card > 1100, card


def test_the_map_crumb_replaces_the_page_crumb_and_starts_at_home(page, map_site):
    """UX.59: one trail on the Galaxy page: Home, Galaxy, the levels, the star."""
    _open_galaxy(page, map_site)
    assert not page.locator("nav.breadcrumbs").is_visible()
    main = page.locator("main > #galaxymap3d-crumbs")
    assert main.count() == 1 and main.is_visible()
    first = page.locator("#galaxymap3d-crumbs > ol > li").first
    assert first.inner_text() == "Home" and first.locator("a").get_attribute("href") == "/"
    assert page.locator("#galaxymap3d-crumbs .galaxy-bookmark").is_visible()


# --- MAP.125: a system opened in place ------------------------------------------------

def _galaxy_system(page):
    return page.evaluate("() => { const s = document.querySelector('#galaxymap3d-canvas').galaxySystem(); "
                         "return s && {id: s.id, name: s.name, mode: s.mode, refs: s.refs, flying: s.flying, "
                         "following: s.following, distance: s.distance, selected: s.selected}; }")


def _wait_system(page, present=True):
    page.wait_for_function(
        "(present) => { const c = document.querySelector('#galaxymap3d-canvas'); "
        "const s = c.galaxySystem && c.galaxySystem(); return present ? !!s && !s.flying : !s; }",
        arg=present, timeout=20000)
    _settle(page)


def _open_system_in_place(page, map_site, name="Middling Sun"):
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_selector(".starmap-sr-list button", state="attached")
    page.locator(".starmap-sr-list button", has_text=name).first.evaluate("b => b.click()")
    page.locator("#galaxymap3d-info button", has_text="Open system here").click()
    _wait_system(page)


def test_galaxy_map_names_the_bodies_of_an_open_system_down_to_a_planet(page, map_site):
    """MAP.268 (MAP.231): the bodies of a system opened in place are named by the label engine once they are a few
    pixels across: flying to a planet names it (and the labels never overlap), the Body names toggle takes them away,
    and no stage is involved from the sector to the planet."""
    _open_system_in_place(page, map_site)
    button = page.locator(".starmap-sr-system button", has_text="(planet)").first
    planet = button.inner_text().split(" (")[0]
    button.evaluate("b => b.click()")
    page.locator("#galaxymap3d-info button", has_text="Go to").click()
    _wait_system(page)
    page.wait_for_function(f"([name]) => ({LABELS})().some((label) => label.kind === 'body' && label.text === name)", arg=[planet])
    page.click("#galaxymap3d-menu > summary")
    page.click('#galaxymap3d-kinds [data-label="body"]')
    page.wait_for_function(f"() => !({LABELS})().some((label) => label.kind === 'body')")


def test_the_dots_and_belts_of_an_open_system_are_built_in_the_systems_own_frame(page, map_site):
    """MAP.240: a system opened in place is built in its own units inside a group at the system, so its body dots and belt particles hold
    offsets from the system (in parsecs, well inside the bound), not the system's place in the galaxy."""
    _open_system_in_place(page, map_site)
    extents = page.evaluate(BUFFER_EXTENTS)
    layers = {e["layer"] for e in extents}
    assert "system dots" in layers, layers
    for e in extents:
        if e["layer"].startswith("system"):
            assert e["extentPc"] <= SAFE_EXTENT_PC and e["extentPc"] < 1.0, e
            assert max(abs(v) for v in e["origin"]) > 0.5, "the group sits at the system's place in the galaxy"


def test_the_bodies_of_an_open_system_carry_their_own_symbols_by_the_same_layers(page, map_site):
    """MAP.321: with a system open in place its planets, moons and belt carry their own symbol (an inhabited Green world, a gas
    giant, a barren moon, an asteroid belt) in the layers the Menu has on, the symbols ride the bodies, a layer that is off draws
    none, and leaving the system takes them away."""
    _open_system_in_place(page, map_site)
    ids = ("marker-inhabited", "marker-habitable", "marker-uninhabitable", "marker-gas-giant", "marker-asteroid-field")
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').__bodyMarkersOn()") is True
    _poi_layers(page, *ids)
    found = _markers_until(page, lambda entries: {e["kind"] for e in entries} >= {"marker-inhabited", "marker-gas-giant",
                                                                                  "marker-uninhabitable", "marker-asteroid-field"})
    kinds = sorted(e["kind"] for e in found if e["kind"] in ids)
    assert kinds.count("marker-inhabited") == 1 and kinds.count("marker-gas-giant") == 1, kinds
    assert kinds.count("marker-asteroid-field") == 1 and kinds.count("marker-uninhabitable") >= 1, kinds
    assert "marker-habitable" not in kinds, "the Green world is inhabited, so it carries the inhabited symbol alone"
    inhabited = next(e for e in found if e["kind"] == "marker-inhabited")
    assert inhabited["name"].endswith(" I"), inhabited
    # MAP.321 d: the inhabited body's name is offered in the inhabited slice of the label budget, the rest stay in the uninhabited one
    slices = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').__mapBodyNames().map((l) => [l.text, l.slice.domain, !!l.slice.inhabited])")
    assert [s for s in slices if s[1] == "inhabited"] == [[inhabited["name"], "inhabited", True]], slices
    assert any(s[1] == "uninhabited" for s in slices), slices
    # the symbol rides its body: where the marker is is where the body is
    refs = _galaxy_system(page)["refs"]
    assert "planet:1" in refs
    spot = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxySystem().screen('planet:1')")
    assert spot is not None
    # a layer that is off draws none of its symbols
    page.locator('#galaxymap3d-kinds button[data-layer="marker-gas-giant"]').click()
    found = _markers_until(page, lambda entries: not any(e["kind"] == "marker-gas-giant" for e in entries))
    assert not [e for e in found if e["kind"] == "marker-gas-giant"]
    assert any(e["kind"] == "marker-inhabited" for e in found)
    # leaving the system hands the symbols back to the stars
    page.click('#galaxymap3d-controls [data-action="up"]')
    _wait_system(page, present=False)
    found = _markers_until(page, lambda entries: not any(e["kind"] == "marker-inhabited" and e["name"].endswith(" I") for e in entries))
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').__bodyMarkersOn()") is False
    assert not [e for e in found if e["kind"] == "marker-inhabited" and e["name"].endswith(" I")]


def test_the_bodies_of_an_open_system_carry_the_life_and_class_symbols_when_those_layers_are_on(page, map_site):
    """MAP.343: with a system open in place its bodies carry the symbols of the reasons to look (life, sentient life, a civilization)
    and of the planet classes when the Menu has them on (they are off until chosen), the strongest reason that is on names a body's
    symbol, and a layer that is off draws none."""
    _open_system_in_place(page, map_site)
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').__bodyMarkersOn()") is True
    _poi_layers(page)    # no world reason on
    found = _markers_until(page, lambda entries: not [e for e in entries if e["kind"] != "marker-black-hole"])
    assert not [e for e in found if e["kind"] != "marker-black-hole"], "the life and class layers are off until chosen"
    kinds = page.locator("#galaxymap3d-kinds")

    def named(layer):
        entries = _markers_until(page, lambda found: any(e["kind"] == layer for e in found))
        return sorted(e["name"] for e in entries if e["kind"] == layer)

    kinds.locator('button[data-layer="marker-life"]').click()
    life = named("marker-life")
    assert len(life) == 2 and life[0].endswith(" I") and life[1].endswith(" I a"), life    # the first planet and its first moon
    kinds.locator('button[data-layer="marker-sentient"]').click()
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-sentient" for e in entries))
    sentient = [e["name"] for e in found if e["kind"] == "marker-sentient"]
    life_left = [e["name"] for e in found if e["kind"] == "marker-life"]
    assert len(sentient) == 1 and sentient[0].endswith(" I") and sentient[0] not in life_left, (sentient, life_left)
    kinds.locator('button[data-layer="marker-civilization"]').click()
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-civilization" for e in entries))
    civ = [e["name"] for e in found if e["kind"] == "marker-civilization"]
    assert civ == sentient and not [e for e in found if e["kind"] == "marker-sentient"], (civ, found)
    kinds.locator('button[data-layer="marker-class-F"]').click()
    kinds.locator('button[data-layer="marker-class-J"]').click()
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-class-J" for e in entries))
    assert [e["name"] for e in found if e["kind"] == "marker-class-J"] and [e["name"] for e in found if e["kind"] == "marker-class-J"][0].endswith(" II")
    assert not [e for e in found if e["kind"] == "marker-class-F"], "the first planet is a civilization, which is the stronger reason"
    kinds.locator('button[data-layer="marker-life"]').click()
    kinds.locator('button[data-layer="marker-sentient"]').click()
    kinds.locator('button[data-layer="marker-civilization"]').click()
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-class-F" for e in entries))
    assert [e["name"] for e in found if e["kind"] == "marker-class-F"][0].endswith(" I")
    assert not [e for e in found if e["kind"] in ("marker-life", "marker-sentient", "marker-civilization")]


def _nearest_that_way(here, others, direction, cone_deg=60.0):
    """The rule of static/systemstep.js worked out independently: the nearest of `others` (name, x, y, z) from `here` in
    the cone about `direction` in the galaxy's frame at `here`, or None."""
    x, y, z = here
    plane = math.hypot(x, y)
    r = math.hypot(plane, z)
    east = (-y / plane, x / plane, 0.0)
    south = (z * x / (plane * r), z * y / (plane * r), -plane / r)
    out = (x / r, y / r, z / r)
    axis = {"out": out, "in": tuple(-v for v in out), "east": east, "west": tuple(-v for v in east),
            "south": south, "north": tuple(-v for v in south)}[direction]
    best = None
    for name, ox, oy, oz in others:
        offset = (ox - x, oy - y, oz - z)
        distance = math.dist(offset, (0.0, 0.0, 0.0))
        if distance > 0 and sum(a * b for a, b in zip(offset, axis)) / distance >= math.cos(math.radians(cone_deg)):
            if best is None or distance < best[0]:
                best = (distance, name)
    return best[1] if best else None


def _system_places():
    """Every fixture system's place in the galaxy's frame, pc, by name: a sector's scene lies in the sector's own frame, turned about the
    z axis to the sector's azimuth (the Fixture Prime wedge 1 at 45 degrees), its stored x and y in mpc from the sector's centre."""
    places = {}
    for sector_id, systems in ((SECTOR_ID, SYSTEMS), (BESIDE_ID, BESIDE_SYSTEMS)):
        centre = _sector_center_pc(*next(s[1:4] for s in SECTORS if s[0] == sector_id))
        azimuth = math.atan2(centre[1], centre[0])
        for _i, name, x, y, z, _s in systems:
            places[name] = (centre[0] + (x * math.cos(azimuth) - y * math.sin(azimuth)) / 1000.0,
                            centre[1] + (x * math.sin(azimuth) + y * math.cos(azimuth)) / 1000.0, centre[2] + z / 1000.0)
    return places


def _step_buttons(page):
    return page.eval_on_selector_all(
        "#galaxymap3d-controls [data-action=step]",
        "els => Object.fromEntries(els.map(e => [e.dataset.step, !e.disabled]))")


def test_the_step_buttons_and_keys_go_from_an_open_system_to_the_nearest_system_that_way(page, map_site):
    """MAP.121, MAP.329: with a system opened in place the step buttons and keys take the view to the nearest system lying that
    way (north, south, east, west, out, in, in the galaxy's frame at the system), keeping the scale, across the sector's edge
    too once the sectors round it are read; a button is off where no system lies that way."""
    _open_system_in_place(page, map_site)
    places = _system_places()
    here = places["Middling Sun"]
    others = [(name, *place) for name, place in places.items() if name != "Middling Sun"]
    expected = {d: _nearest_that_way(here, others, d) for d in ("north", "south", "east", "west", "in", "out")}
    assert any(expected.values()), expected
    # The sectors round the system are read after it opens: a step goes where it should once they are.
    page.wait_for_function("(want) => Object.entries(want).every(([d, name]) => document.querySelector('#galaxymap3d-canvas').galaxyStepTo(d) === name)",
                           arg=expected, timeout=20000)
    assert _step_buttons(page) == {d: name is not None for d, name in expected.items()}, (_step_buttons(page), expected)
    direction = next(d for d, name in expected.items() if name and name in {n for _i, n, *_r in SYSTEMS})
    was = _galaxy_system(page)["id"]
    page.click("#galaxymap3d-menu > summary")
    page.click(f'#galaxymap3d-controls [data-step="{direction}"]')
    page.wait_for_function("([was]) => { const s = document.querySelector('#galaxymap3d-canvas').galaxySystem(); return !!s && s.id !== was && !s.flying; }",
                           arg=[was], timeout=20000)
    _settle(page)
    state = _galaxy_system(page)
    assert state["name"] == expected[direction] and state["mode"] == "compressed", (direction, expected, state)
    assert f"&system={state['id']}" in _query(page)
    assert _crumbs(page)[-1] == expected[direction], _crumbs(page)


def test_a_step_across_the_sector_edge_opens_the_system_beside_in_its_own_sector(page, map_site):
    """MAP.329: from Middling Sun, whose nearest system to the east is in the next sector (Across The Edge, 0.6 pc away), the east
    step opens that sector beside the first and the system in it, keeping the scale, each with its history entry; Back returns
    to the system it came from."""
    _open_system_in_place(page, map_site)
    places = _system_places()
    others = [(name, *place) for name, place in places.items() if name != "Middling Sun"]
    assert _nearest_that_way(places["Middling Sun"], others, "east") == "Across The Edge"
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyStepTo('east') === 'Across The Edge'", timeout=20000)
    was = _galaxy_system(page)["id"]
    page.click("#galaxymap3d-menu > summary")
    page.click('#galaxymap3d-controls [data-step="east"]')
    page.wait_for_function("([was]) => { const s = document.querySelector('#galaxymap3d-canvas').galaxySystem(); return !!s && s.id !== was && !s.flying; }",
                           arg=[was], timeout=30000)
    _settle(page)
    state = _galaxy_system(page)
    assert state["name"] == "Across The Edge" and state["mode"] == "compressed", state
    assert _query(page) == f"sector={BESIDE_CODE}&open=1&system={state['id']}", _query(page)
    crumbs = _crumbs(page)
    assert crumbs[-1] == "Across The Edge" and crumbs[-2].startswith("Sector"), crumbs
    page.go_back()    # the system's entry, then the sector's own
    page.wait_for_function("() => !new URL(location.href).searchParams.get('system')", timeout=30000)
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=30000)
    _settle(page)
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxySystem()", timeout=30000)
    page.go_back()
    page.wait_for_function("([was]) => { const s = document.querySelector('#galaxymap3d-canvas').galaxySystem(); return !!s && s.id === was; }",
                           arg=[was], timeout=30000)


@pytest.mark.parametrize("gap", [-1, 0, 30, 60, 90, 150, 400, 1000])
def test_two_quick_presses_of_back_across_a_sector_change_land_on_the_system_the_second_one_asks_for(page, map_site, gap):
    """MAP.352: Back pressed twice with no wait, from Middling Sun reached by stepping west from Across The Edge (the sector beside),
    goes to the first sector's own entry and then to Across The Edge in the sector beside: the second press, a system in a different
    sector, is not lost to the first one's move, and the map ends on that system."""
    _open_system_in_place(page, map_site)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyStepTo('east') === 'Across The Edge'", timeout=20000)
    page.click("#galaxymap3d-menu > summary")
    for direction, name in (("east", "Across The Edge"), ("west", "Middling Sun")):
        page.wait_for_function("([d, n]) => document.querySelector('#galaxymap3d-canvas').galaxyStepTo(d) === n", arg=[direction, name], timeout=30000)
        before = _galaxy_system(page)["id"]
        page.click(f'#galaxymap3d-controls [data-step="{direction}"]')
        page.wait_for_function("([was]) => { const s = document.querySelector('#galaxymap3d-canvas').galaxySystem(); return !!s && s.id !== was && !s.flying; }",
                               arg=[before], timeout=30000)
        _settle(page)
        assert _galaxy_system(page)["name"] == name
    beside = next(i for i, n, *_r in BESIDE_SYSTEMS if n == "Across The Edge")
    if gap < 0:
        page.go_back()
        page.go_back()
    else:
        page.evaluate("(gap) => new Promise((done) => { const back = document.querySelector('#galaxymap3d-controls [data-action=back]'); back.click(); setTimeout(() => { back.click(); done(); }, gap); })", gap)
    page.wait_for_function("([wanted]) => { const s = document.querySelector('#galaxymap3d-canvas').galaxySystem(); return !!s && s.id === wanted && !s.flying; }",
                           arg=[beside], timeout=40000)
    _settle(page)
    assert _galaxy_system(page)["id"] == beside
    assert f"system={beside}" in _query(page) and f"sector={BESIDE_CODE}" in _query(page), _query(page)


def test_a_close_pairs_two_orbits_are_drawn_each_at_its_share(page, map_site):
    """MAP.136: both stars of a close pair circle the barycenter on their own
    line, the lighter one on the wider (the fixture pair's mass fraction is
    0.25, so its second star keeps three times the first's distance)."""
    _open_system_in_place(page, map_site)
    radii = page.evaluate("() => { const s = document.querySelector('#galaxymap3d-canvas').galaxySystem(); "
                          "return ['star:1', 'star:2'].map((r) => s.trailRadius(r)); }")
    assert all(r and r > 0 for r in radii), radii
    assert 1.5 < radii[1] / radii[0] < 6, radii


def test_galaxy_map_opens_a_system_in_place_and_zooms_to_a_moon(page, map_site):
    """MAP.125: a star of the open sector opens its system in the map itself:
    its bodies are picked and flown to (down to a moon, nine orders of
    magnitude inside the sector's fit), the scale can be true, Up returns to
    the sector, and Back, Forward and a reload keep the system."""
    _open_system_in_place(page, map_site)
    system_id = next(i for i, n, *_rest in SYSTEMS if n == "Middling Sun")
    assert _query(page) == f"sector={PRIME_CODE}&open=1&system={system_id}"
    state = _galaxy_system(page)
    assert state["id"] == system_id and state["mode"] == "compressed"
    assert {"star:1", "star:2", "planet:1", "planet:2", "moon:1", "comet:1", "comet:2"} <= set(state["refs"])
    crumbs = _crumbs(page)
    assert crumbs[-1] == "Middling Sun" and crumbs[-2].startswith("Sector"), crumbs
    fit = state["distance"]

    # A planet picked from the list shows its panel; Go to goes in close.
    page.locator(".starmap-sr-system button", has_text="(planet)").first.evaluate("b => b.click()")
    assert "Go to" in [b.inner_text() for b in page.locator("#galaxymap3d-info button").all()]
    # MAP.126: the picked body's path is drawn full and the others' fade back.
    picked = _galaxy_system(page)["selected"]
    opacity = page.evaluate("(refs) => { const s = document.querySelector('#galaxymap3d-canvas').galaxySystem(); "
                            "return refs.map((r) => s.trailOpacity(r)); }", [picked, "planet:1", "planet:2"])
    assert opacity[0] == 1
    others = [o for ref, o in zip([picked, "planet:1", "planet:2"], opacity) if ref != picked and o is not None]
    assert others and all(o < 0.5 for o in others), opacity
    page.locator("#galaxymap3d-info button", has_text="Go to").click()
    _wait_system(page)
    near = _galaxy_system(page)["distance"]
    assert near < fit / 5, (near, fit)
    page.locator(".starmap-sr-system button", has_text="(moon)").first.evaluate("b => b.click()")
    page.locator("#galaxymap3d-info button", has_text="Go to").click()
    _wait_system(page)
    assert _galaxy_system(page)["distance"] < near
    assert _galaxy_system(page)["following"] == "moon:1"
    assert _query(page).endswith("&object=moon:1"), _query(page)

    # The scale can be true; the same bodies are still there.
    page.locator("#galaxymap3d-info button", has_text="True scale").click()
    _wait_system(page)
    assert _galaxy_system(page)["mode"] == "true"
    assert "star:1" in _galaxy_system(page)["refs"]

    # Up goes to the sector; Back to the system; a reload keeps it.
    page.click('#galaxymap3d-controls [data-action="up"]')
    _wait_system(page, present=False)
    assert _query(page) == OPEN_PRIME[1:]
    page.go_back()
    _wait_system(page)
    page.reload()
    page.wait_for_selector(".starmap-sr-list button", state="attached")
    _wait_system(page)
    assert _galaxy_system(page)["id"] == system_id
    assert f"&system={system_id}" in _query(page)
    # The body the address names is picked and flown to again.
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxySystem().following === 'moon:1'", timeout=20000)
    assert _galaxy_system(page)["selected"] == "moon:1"


def test_galaxy_map_wheel_carries_the_zoom_into_a_system_and_back(page, map_site):
    """MAP.125: zooming right in on an opened sector opens the star picked
    (else the one nearest the middle of the view); zooming right out of the system comes back to
    the sector."""
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_selector(".starmap-sr-list button", state="attached")
    page.locator(".starmap-sr-list button", has_text="Middling Sun").first.evaluate("b => b.click()")
    box = page.locator("#galaxymap3d-canvas").bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    for _ in range(40):
        page.mouse.wheel(0, -400)
        page.wait_for_timeout(30)
    _wait_system(page)
    state = _galaxy_system(page)
    assert state["id"] and f"&system={state['id']}" in _query(page)
    page.wait_for_timeout(1000)
    # One gesture carries across one stage (MAP.152): the scroll stops at the sector, however it runs on.
    assert _wheel_until(page, 400, lambda: not (_galaxy_system(page) or {}).get("id"), ticks=80), "the wheel never came back out of the system"
    _wait_system(page, present=False)
    assert _query(page) == OPEN_PRIME[1:]


def test_galaxy_map_color_by_switches_the_fill_with_a_legend_and_keeps_it_in_the_url(page, map_site):
    """MAP.131: the Menu's Color by recolours the blocks by one statistic,
    shows what the colours mean, keeps the choice in the address and goes
    back to the default."""
    _open_galaxy(page, map_site, f"?sector={PRIME_CODE}")
    legend = page.locator("#galaxymap3d-legend")
    assert legend.is_hidden()
    page.locator("#galaxymap3d-menu > summary").click()
    for mode in ("density", "age", "luminosity", "stars", "habitable"):
        page.select_option("#galaxymap3d-color-by", mode)
        _settle(page)
        assert f"color={mode}" in _query(page)
        assert legend.is_visible(), mode
        assert legend.locator(".galaxy-legend-low").inner_text() != ""
        assert legend.locator(".galaxy-legend-title").inner_text() != ""
    # Reload on the address: the select and the mode follow it.
    page.reload()
    page.wait_for_selector("#galaxymap3d-steps [data-steps-panel] li", state="attached")
    assert page.input_value("#galaxymap3d-color-by") == "habitable"
    legend.wait_for(state="visible", timeout=15000)  # drawn once the nodes are back; under load that is not instant (TEST.126)
    assert legend.is_visible()
    page.locator("#galaxymap3d-menu > summary").click()
    page.select_option("#galaxymap3d-color-by", "default")
    _settle(page)
    assert "color=" not in _query(page)
    assert legend.is_hidden()


def test_galaxy_map_opens_zoomed_in_on_the_charted_space(page, map_site):
    """MAP.124: with no view in the address the galaxy opens as close as its
    zoom range allows with every charted sector in view; Reset goes back to
    the whole galaxy and the wheel still zooms out to it."""
    _open_galaxy(page, map_site, "")
    zoom = page.evaluate(CAMERA)["zoom"]
    assert 0.5 <= zoom < 1, zoom
    page.locator("#galaxymap3d-menu > summary").click()
    page.click('#galaxymap3d-controls [data-action="reset-view"]')
    _settle(page)
    assert page.evaluate(CAMERA)["zoom"] == pytest.approx(1)
    # A view the address names is not moved.
    _open_window(page, map_site, 3, 250, 0, 0)
    # MAP.337: the window is a camera move, so its fit is the distance the camera rests at (it is a ratio of the whole galaxy's fit now).
    rest = page.evaluate(CAMERA)["dist"]
    assert rest < 1e4, rest
    page.wait_for_timeout(1000)
    assert page.evaluate(CAMERA)["dist"] == pytest.approx(rest)


def test_galaxy_map_star_classes_and_the_luminosity_floor_hide_stars_and_ride_in_the_url(page, map_site):
    """MAP.123: the Menu's star classes hide the stars of a class, the
    luminosity slider hides the dimmer ones, both go to the screen-reader
    list and the address, and a reload keeps them."""
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_selector(".starmap-sr-list button", state="attached")

    def listed():
        return page.evaluate("""() => Array.from(document.querySelectorAll('.starmap-sr-list li'))
            .filter((li) => !li.hidden && li.querySelector('button') && !li.querySelector('button').hidden
                    && !li.hidden).length""")

    page.locator("#galaxymap3d-menu > summary").click()
    classes = page.locator("#galaxymap3d-kinds [data-star-class]")
    assert classes.count() >= 1
    before = listed()
    first = classes.first
    letter = first.get_attribute("data-star-class")
    first.click()
    assert first.get_attribute("aria-pressed") == "false"
    assert f"stars={letter}" in _query(page)
    assert listed() < before
    first.click()
    assert "stars=" not in _query(page)
    assert listed() == before

    # The floor: the top of the slider hides every star.
    page.evaluate("""() => { const s = document.getElementById('galaxymap3d-lum'); s.value = s.max;
                              s.dispatchEvent(new Event('input', {bubbles: true})); }""")
    assert "lum=" in _query(page)
    assert listed() < before
    page.reload()
    page.wait_for_selector(".starmap-sr-list button", state="attached")
    assert "lum=" in _query(page)
    assert listed() < before


def test_galaxy_map_bookmark_keeps_the_view_filters(page, map_site):
    """MAP.123: a bookmark of a stage view keeps what the map shows (hidden
    kinds, color), so opening it restores them."""
    _open_galaxy(page, map_site)
    page.locator("#galaxymap3d-menu > summary").click()
    page.select_option("#galaxymap3d-color-by", "age")
    _settle(page)
    page.locator("#galaxymap3d-crumbs .galaxy-bookmark").click()
    saved = _bookmarks(page)
    assert len(saved) == 1 and saved[0]["kind"] == "stage" and "color=age" in saved[0]["value"]


def test_sector_map_highlight_draws_a_phenomenon_kind_larger_and_keeps_it_in_the_url(page, map_site):
    """MAP.123: the Menu's Highlight draws a kind of phenomenon larger and
    brighter, keeps the choice in the address, and a reload restores it."""
    _open_sector(page, map_site)
    _open_menu(page)
    page.wait_for_selector('#galaxymap3d-kinds [data-highlight-kind]', state="attached")
    assert page.locator('#galaxymap3d-kinds [data-highlight-kind="star"]').count() == 0
    toggle = page.locator('#galaxymap3d-kinds [data-highlight-kind="roguePlanet"]')
    assert toggle.get_attribute("aria-pressed") == "false"
    # MAP.137: rogue planets start hidden, and a hidden kind has nothing to highlight; show them first.
    page.locator('#galaxymap3d-kinds button[data-kind="roguePlanet"]').click()
    page.locator("#galaxymap3d-info").evaluate("panel => panel.textContent = ''")
    # The shots are taken with the Menu closed: its panel hangs over the map.
    page.locator("#galaxymap3d-menu").evaluate("m => { m.open = false; }")
    before = _shot(page, SECTOR_CANVAS)
    page.locator("#galaxymap3d-menu").evaluate("m => { m.open = true; }")
    toggle.click()
    assert toggle.get_attribute("aria-pressed") == "true"
    assert "mark=roguePlanet" in page.url
    page.locator("#galaxymap3d-menu").evaluate("m => { m.open = false; }")
    assert _shot(page, SECTOR_CANVAS) != before, "a highlight changes the map"
    page.reload()
    page.wait_for_selector('#galaxymap3d-kinds [data-highlight-kind="roguePlanet"]', state="attached")
    _open_menu(page)
    assert page.locator('#galaxymap3d-kinds [data-highlight-kind="roguePlanet"]').get_attribute("aria-pressed") == "true"
    page.locator('#galaxymap3d-kinds [data-highlight-kind="roguePlanet"]').click()
    assert "mark" not in page.url


def test_galaxy_map_star_filters_apply_to_the_stars_drawn_at_galaxy_scale(page, map_site):
    """MAP.123: with no sector open the Menu still offers the star classes and
    the luminosity floor, and they thin the stars the nodes draw."""
    _open_galaxy(page, map_site)
    _settle(page)
    count = "() => document.querySelector('#galaxymap3d-canvas').galaxyStarFrame().count"
    before = page.evaluate(count)
    assert before > 0
    page.locator("#galaxymap3d-menu > summary").click()
    page.wait_for_selector("#galaxymap3d-kinds [data-star-class]", state="attached")
    assert page.locator("#galaxymap3d-kinds [data-star-class]").count() == 8
    page.evaluate("""() => { const s = document.getElementById('galaxymap3d-lum'); s.value = s.max;
                              s.dispatchEvent(new Event('input', {bubbles: true})); }""")
    assert page.evaluate(count) < before
    page.evaluate("""() => { const s = document.getElementById('galaxymap3d-lum'); s.value = 0;
                              s.dispatchEvent(new Event('input', {bubbles: true})); }""")
    assert page.evaluate(count) == before
    # MAP.166: the galaxy view is never complete, so floor 0 names the dimmest
    # star it carries rather than "every star".
    label = page.locator("label.galaxy-lum span").inner_text()
    assert "every star" not in label and "and up" in label, label
    for button in page.locator("#galaxymap3d-kinds [data-star-class]").all():
        button.click()
    assert page.evaluate(count) < before


# --- MAP.158: nodes kept in IndexedDB, a gentler prefetch ---------------------------------

STORED_NODES = """() => new Promise(resolve => {
    const open = indexedDB.open("planetgen-nodes");
    open.onerror = () => resolve(-1);
    open.onsuccess = () => {
        const count = open.result.transaction("nodes").objectStore("nodes").count();
        count.onsuccess = () => resolve(count.result);
        count.onerror = () => resolve(-1);
    };
})"""


def _node_requests(page):
    """A list that fills with (seconds since now, [node keys]) for each /galaxy/nodes request."""
    import time
    began = time.monotonic()
    found = []

    def note(request):
        if "/galaxy/nodes" in request.url:
            keys = parse_qs(urlparse(request.url).query)["keys"][0].split(",")
            found.append((time.monotonic() - began, keys))
    page.on("request", note)
    return found


def test_galaxy_map_keeps_its_nodes_in_indexeddb_for_a_revisit(page, map_site):
    """MAP.158: the nodes fetched on a visit are stored in IndexedDB (not
    localStorage), so a reload finds them and asks the server for none of them."""
    first = _node_requests(page)
    _open_galaxy(page, map_site, f"?sector={PRIME_CODE}")
    page.wait_for_timeout(2500)
    asked = {key for _at, keys in first for key in keys}
    assert asked, "the first visit fetched nodes"
    # The nodes are written to IndexedDB after they arrive; under load the fixed wait above is not always enough (TEST.126).
    deadline = time.monotonic() + 20
    while page.evaluate(STORED_NODES) < len(asked) and time.monotonic() < deadline:
        page.wait_for_timeout(250)
    assert page.evaluate(STORED_NODES) >= len(asked)
    assert not page.evaluate("Object.keys(localStorage).filter(k => k.startsWith('planetgen:node:') && !k.endsWith('@stamp'))"), \
        "no node is kept in localStorage any more"

    again = _node_requests(page)
    page.reload(wait_until="load")
    page.wait_for_selector(GALAXY_CANVAS, state="attached")
    _settle(page)
    page.wait_for_timeout(2500)
    assert not asked & {key for _at, keys in again for key in keys}, "a revisit found the nodes it kept"


def test_galaxy_map_prefetches_only_the_way_in_and_only_when_still(page, map_site):
    """MAP.158: after the view's own nodes are in and the camera has been still
    a moment, the nodes one zoom step IN are fetched; nothing for the way out."""
    import math

    from planetgen.galaxy.nodecut import node_cut
    from planetgen.galaxy.nodecut import MAX_NODES
    from planetgen.web.maps.galaxymap3d import CAMERA_FOV_DEG, camera_position_pc, node_bounds_pc, view_radius_bounds
    from tests.map_site_support import EDGE_PC, galaxy_shape

    requests = _node_requests(page)
    _open_galaxy(page, map_site, f"?sector={PRIME_CODE}")
    page.wait_for_timeout(3500)
    near = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyNearField()")
    target, camera = near["target"], near["camera"]
    bounds = node_bounds_pc(EDGE_PC, galaxy_shape())
    height = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').clientHeight")
    focal = height / (2 * math.tan(math.radians(CAMERA_FOV_DEG) / 2))

    def cut_at(eye, scale=1.0, look=target):
        eye = [look[i] + (eye[i] - look[i]) * scale for i in range(3)]
        return set(node_cut(eye, look, EDGE_PC, focal_px=focal, bounds_pc=bounds, max_nodes=MAX_NODES).keys)

    asked = {key for _at, keys in requests for key in keys}
    here, in_one, out_one = cut_at(camera), cut_at(camera, 1.0 / 3.0), cut_at(camera, 3.0)
    # The page also loads the opening view (the whole galaxy from its farthest camera) before it settles on the sector.
    centre = (0.0, 0.0, 0.0)
    opening = cut_at(camera_position_pc(view_radius_bounds(EDGE_PC, galaxy_shape())[1], centre), look=centre)
    way_out = out_one - here - in_one - opening
    assert way_out, "the way out needs nodes of its own for this test to mean anything"
    assert not asked & way_out, f"no zoom-out prefetch: {sorted(asked & way_out)}"
    assert asked <= here | in_one | opening, sorted(asked - (here | in_one | opening))


def test_galaxy_map_does_not_prefetch_on_a_data_saving_connection(page, map_site):
    """MAP.158: `navigator.connection.saveData` turns the prefetch off: a visit that settles asks for the nodes one
    zoom step in (a request marked `X-Map-Prefetch`), and a data saving one does not."""
    def visit(save_data):
        context = page.context.browser.new_context(viewport=VIEWPORT, device_scale_factor=1, reduced_motion="reduce")
        try:
            tab = context.new_page()
            if save_data:
                tab.add_init_script(
                    "Object.defineProperty(navigator, 'connection', {value: {saveData: true}, configurable: true})")
            requests, prefetches = _node_requests(tab), []
            tab.on("request", lambda request: prefetches.append(request.url)
                   if "/galaxy/nodes" in request.url and request.headers.get("x-map-prefetch") else None)
            _open(tab, f"{map_site}/galaxy", GALAXY_CANVAS)
            # The view's own nodes first (batches can pause for seconds on a busy machine), then the prefetch, which
            # waits for the camera to be still; a data saving visit gets the same time and no prefetch.
            deadline = time.monotonic() + 45
            seen, quiet_since = 0, time.monotonic()
            while time.monotonic() < deadline and not prefetches:
                tab.wait_for_timeout(500)
                if len(requests) != seen:
                    seen, quiet_since = len(requests), time.monotonic()
                elif save_data and requests and time.monotonic() - quiet_since > 8.0:
                    break
            return prefetches
        finally:
            context.close()

    assert visit(False), "a visit that settles prefetches the way in"
    assert not visit(True), "a data saving connection prefetches nothing"


# --- MAP.161: the map's modules are preloaded, each fetched once ---------------------------

def test_galaxy_map_fetches_each_of_its_modules_once_with_the_hints_it_was_given(page, map_site):
    """MAP.161: the page's `modulepreload` hints name the same URLs the script's own import()
    asks for, so no module is fetched twice, and every hinted module is one the page used."""
    urls = []
    page.on("request", lambda request: urls.append(request.url) if "/static/" in request.url and ".js" in request.url else None)
    _open_galaxy(page, map_site)
    hinted = page.evaluate("[...document.querySelectorAll('link[rel=modulepreload]')].map(l => l.href)")
    assert len(hinted) >= 20
    repeated = sorted({url for url in urls if urls.count(url) > 1})
    assert not repeated, f"fetched more than once: {repeated}"
    assert set(hinted) <= set(urls)


def _wheel_until(page, deltas, changed, ticks=80):
    """Turns the wheel (as real wheel events, a pause between) until `changed()` holds; returns whether it did."""
    for _ in range(ticks):
        page.mouse.wheel(0, deltas)
        page.wait_for_timeout(30)
        if changed():
            return True
    return False


def test_galaxy_map_wheel_carries_the_zoom_down_through_the_windows_and_back_up(page, map_site):
    """MAP.152, MAP.337: zoomed right in on a window the wheel goes on in, the camera reaching the zoom level of the block under it (the
    window one size down), and zoomed right out it comes back to the level it came from; the stage the map draws is not switched."""
    _open_window(page, map_site, 243, 0, 0, 0)
    box = page.locator(GALAXY_CANVAS).bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    start = page.evaluate(STAGE_M)
    assert _wheel_until(page, -400, lambda: page.evaluate(STAGE_M) != start), "the wheel never carried on into a block"
    inside = page.evaluate(STAGE_M)
    assert inside == 81, f"one level down is the 81-shell window, not {inside!r}"
    _settle(page)
    page.wait_for_timeout(1000)
    assert _wheel_until(page, 400, lambda: page.evaluate(STAGE_M) != inside, ticks=120), "the wheel never came back out"
    _settle(page)
    assert page.evaluate(STAGE_M) == 243, "back out on the parent level, the 243-shell blocks"


def test_nothing_is_clipped_when_you_zoom_into_a_window_or_a_sector(page, map_site):
    """MAP.264: the surroundings stay in view: in a window of sectors, and with a sector opened in place, the galaxy's
    stars are all still drawn but those inside the opened sector's own cell (which that sector draws itself)."""
    _open_window(page, map_site, 1, 1, 0, 1, g=3)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyHiddenStars().total > 0")
    stars = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyHiddenStars()")
    assert stars["hidden"] == 0, f"a window hid {stars['hidden']} of its {stars['total']} stars"
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyHiddenStars().total > 0")
    stars = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyHiddenStars()")
    assert 0 < stars["total"] - stars["hidden"], f"the neighbours are gone: {stars}"


def test_one_wheel_scroll_runs_from_the_galaxy_to_a_sector_with_no_stop_and_no_swoop(page, map_site):
    """MAP.263: the camera is the only state. From the whole galaxy, turned to an angle of the reader's own, one
    unbroken wheel scroll carries the view down through every window to the sectors: the camera is never flown to a
    preset (never `flying`), its turn is never reset, and the distance falls at every tick (no stop at a stage's end)."""
    _open_galaxy(page, map_site)
    box = page.locator(GALAXY_CANVAS).bounding_box()
    centre = (box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    page.mouse.move(*centre)
    page.mouse.down()
    page.mouse.move(centre[0], centre[1] + 60, steps=6)
    page.mouse.up()
    start = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
    assert start["tilt"] > 5, "the reader's own turn is not the preset's"
    seen = []
    for _ in range(260):
        page.mouse.wheel(0, -400)
        page.wait_for_timeout(30)
        camera = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
        assert not camera["flying"], "a flight to a preset"
        seen.append(camera)
        win = parse_qs(_query(page)).get("win", [""])[0]
        if win.startswith("1.") or "open=1" in _query(page):
            break
    else:
        pytest.fail(f"the scroll never reached the sectors: {seen[-1]}")
    assert all(abs(c["tilt"] - start["tilt"]) < 0.5 for c in seen), "the turn was reset on the way"
    # The tick that carries across a stage is the only one that does not come closer: never two in a row.
    stalls = [n for n in range(1, len(seen)) if seen[n]["dist"] >= seen[n - 1]["dist"] * (1 - 1e-9)]
    assert not any(b - a == 1 for a, b in zip(stalls, stalls[1:])), f"the zoom stopped at ticks {stalls} of {len(seen)}"
    assert seen[-1]["dist"] < start["dist"] / 50, (start["dist"], seen[-1]["dist"])


def test_the_wheel_keeps_the_view_over_the_disc_on_the_way_down_to_a_sector(page, map_site):
    """MAP.292: the blocks are not drawn (MAP.265), so the point the wheel zooms toward is not the top face of the
    invisible column of them (kiloparsecs above the disc, where no star is): scrolling in with the pointer anywhere
    over the galaxy carries the view down to a sector with its centre still in the plane the stars lie in."""
    for fraction in (0.5, 0.2, 0.8):
        _open_galaxy(page, map_site)
        page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
        box = page.locator(GALAXY_CANVAS).bounding_box()
        page.mouse.move(box["x"] + box["width"] * fraction, box["y"] + box["height"] * 0.5)
        targets = []

        def reached():
            camera = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")
            targets.append(camera["target"][2])
            return parse_qs(_query(page)).get("win", [""])[0].startswith("1.")

        assert _wheel_until(page, -400, reached, ticks=260), f"the scroll never reached a sector (pointer at {fraction})"
        assert max(abs(z) for z in targets) < 600, (fraction, targets)


# --- The halo beyond the outline (MAP.243) -----------------------------------------------

def _tooltip(page):
    return page.evaluate("() => { const tip = document.querySelector('#galaxymap3d-tooltip'); return tip.hidden ? '' : tip.textContent; }")


def _hover_spot(page, point):
    """Moves the pointer over the point of the scene `point` and returns where it is on the page, once the tooltip
    names something (the map takes a moment after its nodes arrive)."""
    for _ in range(40):
        spot = _screen_of(page, point)
        assert spot, f"{point} is not on the map"
        page.mouse.move(spot[0] + 1, spot[1] + 1)
        page.mouse.move(*spot)
        if _tooltip(page):
            return spot
        page.wait_for_timeout(250)
    pytest.fail(f"nothing under the pointer at {point}")


@pytest.mark.names_hidden
def test_a_black_hole_in_the_halo_is_picked_on_the_whole_galaxy_view(page, halo_site):
    """MAP.243: beyond the galaxy's outline no block is under the pointer; a phenomenon there is picked on the screen
    like one inside a window, its details shown."""
    _open_galaxy(page, halo_site)
    black_hole = HALO_POINTS[0]
    point = (black_hole["x"], black_hole["y"], black_hole["z"])
    assert page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyOutsideOutline(...p)", list(point))
    spot = _hover_spot(page, point)
    assert "Fixture Outlier" in _tooltip(page)
    page.mouse.click(*spot)
    page.wait_for_function("() => { const h = document.querySelector('#galaxymap3d-info h3'); return h && h.textContent.includes('Fixture Outlier'); }")


@pytest.mark.names_hidden
def test_a_star_in_the_halo_is_picked_and_the_wheel_opens_the_window_round_it(page, halo_site):
    """MAP.243: a star of a loose group beyond the outline is picked (its position and brightness shown), and the
    wheel run toward it hands off to a window round the cell that holds it, where there was no block to open before.
    (The black hole markers, MAP.219, are off: on this view a marker's reach covers the star beside it.)"""
    _open_galaxy(page, halo_site, "?hide=marker-black-hole")
    star = HALO_STARS[0]
    point = (star["x"], star["y"], star["z"])
    assert page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyOutsideOutline(...p)", list(point))
    spot = _hover_spot(page, point)
    assert "Star beyond the galaxy's outline" in _tooltip(page)
    page.mouse.click(*spot)
    page.wait_for_function("() => { const h = document.querySelector('#galaxymap3d-info h3'); return h && h.textContent.includes('Star beyond'); }")
    assert _crumbs(page) == ["Galaxy"]
    # The wheel toward the star hands off to the window round the cell that holds it, and on, level by level. Once
    # the star is out of the visible part of the page the wheel goes on at the middle of what is shown.
    viewport = page.viewport_size
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    for _ in range(24):
        spot = _screen_of(page, point)
        if not spot or not (0 <= spot[0] < viewport["width"] and 0 <= spot[1] < viewport["height"]):
            spot = page.evaluate("() => { const r = document.querySelector('#galaxymap3d-canvas').getBoundingClientRect(); return [r.left + r.width / 2, Math.min(r.top + r.height / 2, innerHeight - 10)]; }")
        page.mouse.move(*spot)
        page.mouse.wheel(0, -400)
        page.wait_for_timeout(1000)  # a moment for the next window to arrive
        if len(_crumbs(page)) >= 6:
            break
    _settle(page)
    crumbs = _crumbs(page)
    assert len(crumbs) >= 6 and crumbs[0] == "Galaxy", crumbs
    camera = _camera(page)
    # The blocks are not drawn (MAP.265, MAP.292), so once the star is off the screen the wheel goes on toward the
    # focus plane at the middle of what is shown, not a block face beside the star: a little wider than it was.
    assert math.dist(camera["target"], point) < 2500.0, (camera["target"], point)


@pytest.mark.names_hidden
@pytest.mark.viewport(1280, 900)
def test_a_selection_is_shown_once_and_over_the_map(page, map_site):
    """MAP.256: an opened sector (selected on opening) and a system picked in it each show their information in one
    place, inside the map's window, not also in a box beside it; Close takes it away and the next pick brings it back."""
    _check_a_selection_is_a_popup(page, map_site, pick=True)


@pytest.mark.viewport(390, 844)
def test_a_selection_is_a_popup_over_the_map_on_a_phone_too(page, map_site):
    """MAP.288: there is no panel beside or under the map at any width: on a phone the selection is the same popup."""
    _check_a_selection_is_a_popup(page, map_site)


def _check_a_selection_is_a_popup(page, map_site, pick=False):
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    info = page.locator("#galaxymap3d-info")
    canvas = page.locator(GALAXY_CANVAS).bounding_box()

    def inside_the_window():
        box = info.bounding_box()
        return (box["x"] >= canvas["x"] - 1 and box["y"] >= canvas["y"] - 1
                and box["x"] + box["width"] <= canvas["x"] + canvas["width"] + 1
                and box["y"] + box["height"] <= canvas["y"] + canvas["height"] + 1)

    assert _selection_displays(page) == ["galaxymap3d-info"] and inside_the_window()
    assert info.locator("h3").inner_text() == "Fixture Prime"
    info.locator("button", has_text="Close").click()   # the popup is drawn into the map's picture: clear it to look for stars
    assert _selection_displays(page) == []
    if not pick:
        return
    page.wait_for_timeout(500)
    picked = None
    for x, y in _bright_spots(page, GALAXY_CANVAS):
        page.mouse.click(x, y)
        heading = info.locator("h3").inner_text() if info.locator("h3").count() else None
        if heading in SYSTEM_NAMES:
            picked = heading
            break
    assert picked, "no click on the opened sector picked a star"
    assert _selection_displays(page) == ["galaxymap3d-info"] and inside_the_window()
    info.locator("button", has_text="Close").click()
    assert _selection_displays(page) == []


GEOMETRY = """() => {
    const box = (sel) => { const r = document.querySelector(sel).getBoundingClientRect();
        return {left: r.left, top: r.top, width: r.width, height: r.height, bottom: r.bottom}; };
    const canvas = document.querySelector("#galaxymap3d-canvas");
    return {
        view: box(".galaxymap3d-panel .starmap-viewport"), canvas: box("#galaxymap3d-canvas"),
        controls: box("#galaxymap3d-controls"), column: box(".galaxymap3d-panel .starmap-layout"),
        buffer: [canvas.width, canvas.height], dpr: devicePixelRatio, scrollWidth: document.documentElement.scrollWidth,
        inner: [innerWidth, innerHeight],
    };
}"""


@pytest.mark.parametrize("width,height", [(1000, 900), (1280, 720), (390, 700)])
def test_galaxy_map_window_fills_the_page_width_and_the_height_the_window_leaves(page, map_site, width, height):
    """MAP.291, MAP.350: the map's window is as wide as its column and as tall as the browser window less the header and 20 px
    above and below (at least 320 px), with no fixed aspect ratio, and the canvas's drawing buffer follows it."""
    page.set_viewport_size({"width": width, "height": height})
    _open_galaxy(page, map_site)
    g = page.evaluate(GEOMETRY)
    assert g["view"]["width"] == pytest.approx(g["column"]["width"], abs=1), g
    assert g["canvas"]["width"] == pytest.approx(g["view"]["width"], abs=1) and g["canvas"]["height"] == pytest.approx(g["view"]["height"], abs=1)
    assert g["buffer"] == [round(g["canvas"]["width"] * g["dpr"]), round(g["canvas"]["height"] * g["dpr"])], g
    covered = page.evaluate("() => document.querySelector('header').getBoundingClientRect().bottom")
    assert g["view"]["height"] == pytest.approx(max(320, height - covered - 40), abs=1), (g, covered)
    assert g["scrollWidth"] <= width, g


def test_galaxy_map_window_follows_the_browser_window_when_it_is_resized(page, map_site):
    """MAP.291: making the window taller makes the map taller (and wider, wider); the canvas follows."""
    page.set_viewport_size({"width": 900, "height": 800})
    _open_galaxy(page, map_site)
    first = page.evaluate(GEOMETRY)
    page.set_viewport_size({"width": 1200, "height": 1000})
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').clientWidth > 1000")
    second = page.evaluate(GEOMETRY)
    assert second["view"]["width"] > first["view"]["width"] + 200, (first, second)
    assert second["view"]["height"] > first["view"]["height"] + 100, (first, second)
    assert second["buffer"] == [round(second["canvas"]["width"] * second["dpr"]), round(second["canvas"]["height"] * second["dpr"])]


# --- MAP.295: the info popup follows the map's window ---------------------------------

INFO_BOX = """() => {
    const info = document.querySelector("#galaxymap3d-info");
    const view = document.querySelector(".galaxymap3d-panel .starmap-viewport").getBoundingClientRect();
    const r = info.getBoundingClientRect();
    const rest = info.querySelector(".map-info-rest");
    return {
        inside: r.left >= view.left - 1 && r.right <= view.right + 1 && r.top >= view.top - 1 && r.bottom <= view.bottom + 1,
        fits: info.scrollHeight <= info.clientHeight + 1 && info.scrollWidth <= info.clientWidth + 1,
        scroll: !!rest && !rest.hidden && rest.scrollHeight > rest.clientHeight + 1,
        width: r.width, height: r.height, view: [view.width, view.height], fontSize: parseFloat(getComputedStyle(info).fontSize),
    };
}"""

LONG_CARD = """(rows) => {
    const panel = document.querySelector("#galaxymap3d-canvas").__mapInfo;
    panel.show({title: "A system with a very long name that goes on and on and on", fields: Array.from({length: rows},
        (_, n) => ["Row " + n, "A value with a good deal of text so that it has to wrap onto a second line " + n])});
}"""


@pytest.mark.parametrize("width,height", [(800, 600), (1280, 900), (390, 700), (600, 480)])
def test_galaxy_map_info_popup_stays_inside_the_map_window_at_every_size(page, map_site, width, height):
    """MAP.295: the popup never leaves the map's window and holds all it says (smaller text first, then a scroll
    that the pointer can reach), for a short card and for one far too long for the window."""
    page.set_viewport_size({"width": width, "height": height})
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_timeout(500)
    page.locator(".starmap-sr-markers button", has_text="Fixture Maw").evaluate("b => b.click()")
    short = page.evaluate(INFO_BOX)
    assert short["inside"] and (short["fits"] or short["scroll"]), short
    if width >= 1280:
        assert short["fits"], short   # a window this size holds a black hole's card whole
    page.evaluate(LONG_CARD, 40)
    long = page.evaluate(INFO_BOX)
    assert long["inside"], long
    assert long["fits"] or long["scroll"], long
    assert long["width"] <= long["view"][0] - 12 and long["height"] <= long["view"][1], long


def test_galaxy_map_info_popup_grows_with_the_map_window(page, map_site):
    """MAP.295: a wider map gives the popup more room (up to a limit) and the text stays at a reading size (MAP.325: 13 to 16 px)."""
    widths = []
    for width, height in ((700, 600), (1400, 900)):
        page.set_viewport_size({"width": width, "height": height})
        _open_window(page, map_site, 27, 27, 0, 0, g=3)
        page.wait_for_timeout(500)
        page.locator(".starmap-sr-markers button", has_text="Fixture Maw").evaluate("b => b.click()")
        box = page.evaluate(INFO_BOX)
        widths.append(box["width"])
        assert 13 <= box["fontSize"] <= 16, box
    assert widths[1] > widths[0] + 20, widths


def test_galaxy_map_info_popup_scrolls_with_the_wheel_and_hands_the_rest_to_the_map(page, map_site):
    """MAP.295, MAP.347: over an expanded popup that has more to show the wheel scrolls what is behind Expand; at the end of its
    scroll (or when it has none) the same turn zooms the map under it, so a popup never takes the zoom away."""
    page.set_viewport_size({"width": 800, "height": 600})
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_timeout(500)
    page.locator(".starmap-sr-markers button", has_text="Fixture Maw").evaluate("b => b.click()")
    page.evaluate(LONG_CARD, 40)
    _expand_card(page)
    assert page.evaluate(INFO_BOX)["scroll"]
    box = page.locator("#galaxymap3d-info .map-info-rest").bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    page.mouse.wheel(0, 200)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info .map-info-rest').scrollTop > 0")
    page.evaluate("() => { document.querySelector('#galaxymap3d-info .map-info-rest').scrollTop = 0; }")
    camera = "() => JSON.stringify(document.querySelector('#galaxymap3d-canvas').galaxyCamera())"
    before = page.evaluate(camera)
    for _ in range(12):
        page.mouse.wheel(0, -400)
        page.wait_for_timeout(60)
    assert page.evaluate(camera) != before, "the wheel at the top of the popup did not reach the map"


def _open_habitable_card(page, base, name, width=1280, height=900):
    """The Galaxy Map's card for one fixture system with a habitable world, opened from the marker's button."""
    page.set_viewport_size({"width": width, "height": height})
    _open_star_field(page, base)
    page.wait_for_timeout(500)
    _poi_layers(page, "marker-habitable")
    _markers_until(page, lambda found: any(e["kind"] == "marker-habitable" for e in found))
    page.locator(".starmap-sr-markers button", has_text=name).evaluate("b => b.click()")
    assert page.locator("#galaxymap3d-info h3").inner_text() == name


def test_galaxy_map_system_card_lists_the_nearest_systems_each_a_link(page, map_site):
    """MAP.296: under its position a system's card names the systems nearest it, nearest first with their distances,
    each name a link to that system's page; and the card still fits the map's window with the list in it."""
    _open_habitable_card(page, map_site, "Fixture Eden")
    _expand_card(page)    # the list is behind Expand
    near = page.locator("#galaxymap3d-info .map-info-nearest-list li a")
    near.first.wait_for()
    names = near.all_inner_texts()
    assert names[0] == "Fixture Rest" and len(names) <= 3, names
    assert near.first.get_attribute("href") == "/system/9002"
    assert "ly" in page.locator("#galaxymap3d-info .map-info-nearest-list li").first.inner_text()
    assert page.locator("#galaxymap3d-info .map-info-nearest h4").inner_text() == "Nearest systems"
    order = page.evaluate("""() => { const info = document.querySelector('#galaxymap3d-info');
        return [...info.querySelectorAll('dt, h4, .map-info-actions')].map((e) => e.tagName === 'DT' ? e.textContent : e.tagName); }""")
    assert order.index("P") < order.index("Position x, y, z") < order.index("H4"), order
    box = page.evaluate(INFO_BOX)
    assert box["inside"] and (box["fits"] or box["scroll"]), box


CARD_CLEARANCE = """() => {
    const box = (e) => { const r = e.getBoundingClientRect(); return {l: r.left, t: r.top, r: r.right, b: r.bottom}; };
    const info = document.querySelector("#galaxymap3d-info");
    const pad = document.querySelector("#galaxymap3d-navpad");
    const a = box(info), b = box(pad);
    const alpha = (getComputedStyle(info).backgroundColor.match(/[\\d.]+/g) || []).map(Number)[3];
    return {
        overlap: a.l < b.r && b.l < a.r && a.t < b.b && b.t < a.b,
        font: parseFloat(getComputedStyle(info).fontSize), alpha: alpha === undefined ? 1 : alpha,
        more: document.querySelectorAll("#galaxymap3d-info details.map-info-more").length,
        open: document.querySelectorAll("#galaxymap3d-info details.map-info-more[open]").length,
    };
}"""


@pytest.mark.parametrize("width,height", [(1280, 900), (800, 600), (390, 700)])
def test_galaxy_map_card_is_reading_size_translucent_and_clear_of_the_zoom_controls(page, map_site, width, height):
    """MAP.325: a card is set at a reading size, its background a darker see-through panel, its rarely needed facts
    under a More that starts collapsed, and it never lies over the zoom controls in the corner."""
    _open_habitable_card(page, map_site, "Fixture Eden", width, height)
    seen = page.evaluate(CARD_CLEARANCE)
    assert not seen["overlap"], seen
    assert 13 <= seen["font"] <= 16, seen
    assert 0.8 <= seen["alpha"] < 1.0, seen
    assert seen["more"] == 1 and seen["open"] == 0, seen
    assert page.locator("#galaxymap3d-info details.map-info-more summary").inner_text() == "More"
    page.evaluate(LONG_CARD, 40)
    assert not page.evaluate(CARD_CLEARANCE)["overlap"], "a long card stays clear of the controls too"


STACK_BOXES = """() => {
    const box = (e) => { const r = e.getBoundingClientRect(); return {l: r.left, t: r.top, r: r.right, b: r.bottom}; };
    const legend = document.querySelector("#galaxymap3d-layer-legend");
    const scale = document.querySelector("#galaxymap3d-scale");
    const pad = document.querySelector("#galaxymap3d-navpad");
    const view = document.querySelector(".galaxymap3d-panel .starmap-viewport").getBoundingClientRect();
    const a = box(legend), s = box(scale), p = box(pad);
    const hit = (x, y) => x.l < y.r && y.l < x.r && x.t < y.b && y.t < x.b;
    return {legendOverScale: hit(a, s), legendAboveScale: a.b <= s.t + 1, scaleOverPad: hit(s, p), legendOverPad: hit(a, p),
            inside: a.l >= view.left - 1 && s.b <= view.bottom + 1 && s.r <= view.right + 1};
}"""


@pytest.mark.parametrize("width,height", [(1280, 900), (700, 600), (390, 700)])
def test_galaxy_map_legend_sits_above_the_scale_and_covers_neither_it_nor_the_controls(page, map_site, width, height):
    """MAP.323: the layer legend and the scale line are one stack in the lower left, the legend above the scale, at
    any window size and with several legend lines."""
    page.set_viewport_size({"width": width, "height": height})
    _open_sector(page, map_site)
    page.evaluate("""() => {
        const layers = document.getElementById("galaxymap3d-canvas").__mapLayers;
        ["one", "two", "three"].forEach((id) => layers.register({id: "legend-" + id, label: "Layer " + id,
            modes: ["galaxy", "sector"], group: "Show on the map", urlKey: "hide", urlKeyOn: "show",
            legend: {text: "marks the layer called " + id + " on the map", swatch: "#f90"}}));
        layers.refresh();
    }""")
    page.wait_for_timeout(300)
    assert page.locator("#galaxymap3d-layer-legend").is_visible()
    seen = page.evaluate(STACK_BOXES)
    assert seen["legendAboveScale"] and not seen["legendOverScale"], seen
    assert not seen["scaleOverPad"] and not seen["legendOverPad"], seen
    assert seen["inside"], seen


def test_galaxy_map_system_nearest_route_answers_json(page, map_site):
    """MAP.296: `/galaxy/system/<id>/nearest` is the card's data: the nearest neighbours with a page link each, a 404
    for a system that is not there."""
    ok = page.request.get(map_site + "/galaxy/system/9001/nearest")
    assert ok.status == 200
    body = ok.json()
    names = [n["name"] for n in body["nearest"]]
    assert names[0] == "Fixture Rest" and len(names) <= 3 and body["nearest"][0]["url"] == "/system/9002", names
    assert page.request.get(map_site + "/galaxy/system/424242/nearest").status == 404


def test_galaxy_map_names_the_heaviest_sectors_first_by_their_total_mass(page, map_site):
    """MAP.322: the sectors offered a name are ranked by the total mass of what they hold (the stored `mass_sum`), the heaviest first, so
    a crowded region names its heaviest sectors and a sparse one names them all. The fixture's three sectors in this window hold
    1.2e32, 3e31 and 0 kg; the label engine places them in that order."""
    _open_window(page, map_site, 1, 9, 1, 4, g=3)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').__mapSectorNames().length === 3")
    offered = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').__mapSectorNames().map(l => [l.text, l.slice.massKg])")
    assert dict(offered) == {"Fixture Second": 1.2e32, "Fixture High": 3e31, "Fixture Shell Nine": 0}, offered
    page.wait_for_function(f"() => ({LABELS})().filter((label) => label.kind === 'sector').length === 3")
    placed = [label["text"] for label in page.evaluate(LABELS) if label["kind"] == "sector"]
    assert placed == ["Fixture Second", "Fixture High", "Fixture Shell Nine"], placed


BLOCK_LIKE_CARD = """() => {
    const panel = document.querySelector("#galaxymap3d-canvas").__mapInfo;
    const region = document.createElement("p");
    region.className = "map-info-region";
    region.textContent = "In the Fixture Region";
    panel.show({title: "Block 10·1·63 (9 sectors a side)",
        fields: [["Sectors", "469"], ["Charted", "0"], ["Distance from core", "376 pc (1,226 ly)"]],
        more: [["Shells", "90-98"], ["Row", "1 of 12"], ["Wedge", "63 of 130"]],
        after: [region],
        buttons: [{label: "Fly into this block", onClick: function () {}}]});
}"""

CARD_REACH = """() => {
    const info = document.querySelector("#galaxymap3d-info");
    const view = document.querySelector(".galaxymap3d-panel .starmap-viewport").getBoundingClientRect();
    const rest = info.querySelector(".map-info-rest");
    const box = info.getBoundingClientRect();
    const button = [...info.querySelectorAll("button")].find((b) => b.textContent.includes("Fly into this block"));
    const region = info.querySelector(".map-info-region");
    const content = rest.scrollHeight;
    rest.scrollTop = rest.scrollHeight;
    const reach = (el) => { const r = el.getBoundingClientRect(); return r.top >= box.top - 1 && r.bottom <= box.bottom + 1; };
    return {content: content, height: box.height, view: view.height, scrolls: info.scrollHeight > info.clientHeight + 1,
            restScrolls: rest.scrollHeight > rest.clientHeight + 1, inside: box.top >= view.top - 1 && box.bottom <= view.bottom + 1,
            buttonReachable: !!button && reach(button), regionReachable: !!region && reach(region)};
}"""


@pytest.mark.parametrize("width,height", [(1000, 560), (800, 520), (390, 480)])
def test_galaxy_map_info_card_sizes_to_its_content_and_keeps_its_buttons_and_region_line_reachable_in_a_short_window(page, map_site, width, height):
    """MAP.332, MAP.347: the expanded card stays inside the map window, what is behind Expand scrolls inside the card, never the card, and
    its region line and its buttons are never left below the fold of a short window."""
    page.set_viewport_size({"width": width, "height": height})
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_timeout(500)
    page.evaluate(BLOCK_LIKE_CARD)
    _expand_card(page)
    page.evaluate("() => { document.querySelector('#galaxymap3d-info .map-info-more').open = true; }")   # the card fits itself again
    page.wait_for_timeout(300)
    got = page.evaluate(CARD_REACH)
    assert got["buttonReachable"] and got["regionReachable"] and got["inside"] and not got["scrolls"], got


# --- MAP.333: zoom Z0: moves are held by the global limits, not the stage's range -------------------------------------

CAMERA_NOW = "() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()"
STAGE_M = "() => { const s = document.querySelector('#galaxymap3d-canvas').galaxyZoomLevel(); return s.win ? s.win.centre.m : null; }"
"""str: The block size of the zoom level the camera is at (MAP.337: the camera's distance alone, not a stage the map has switched to)."""


def _slider_to(page, share, seconds=7.0):
    """Sets the zoom slider as a drag to `share` would (an input event) and waits for the view to come to rest."""
    page.evaluate("(v) => { const s = document.querySelector('#galaxymap3d-zoom-slider'); s.value = String(Math.round(v * 1000)); "
                  "s.dispatchEvent(new Event('input', {bubbles: true})); }", share)
    page.wait_for_timeout(int(seconds * 1000))
    return page.evaluate(CAMERA_NOW)


def _slider_distance(share, far=37202.0, floor=1e-12):
    """The camera distance (parsecs) the slider's `share` stands for: one log scale from the whole galaxy to the floor."""
    return math.exp(math.log(far) - share * (math.log(far) - math.log(floor)))


def test_a_fly_to_lands_where_it_was_asked_with_no_snap_back_to_the_stage_range(page, map_site):
    """MAP.333: a region, a sector or a system flown to is reached at the distance asked for, not at the nearest its stage allows,
    and the camera stays there: from the whole galaxy a region of 3 pc is framed from a few parsecs, not from half the galaxy."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, map_site)
    page.wait_for_function("() => { const l = document.querySelector('#galaxymap3d-canvas').__mapRegionList; return l && l.order().length === 2; }", timeout=20000)
    fly = "([c, r]) => document.querySelector('#galaxymap3d-canvas').galaxyGoTo({center: c, radius: r})"
    for centre, radius in (([2.83, 2.83, 0.0], 12.0), ([2.83, 2.83, 0.0], 3.0)):
        page.evaluate(fly, [centre, radius])
        page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=30000)
        _settle(page)
        landed = page.evaluate(CAMERA_NOW)["dist"]
        assert landed < 4.0 * radius, (radius, landed)
        page.wait_for_timeout(1200)
        assert abs(page.evaluate(CAMERA_NOW)["dist"] - landed) < 1e-6 * landed, "the view snapped back after the flight"
    # A sector opened in place, and a system gone in to: a tiny sphere in each is framed closely too.
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_selector(".starmap-sr-list button", state="attached")
    page.locator(".starmap-sr-list button", has_text="Middling Sun").first.evaluate("b => b.click()")
    _settle(page)
    _slider_to(page, 0.6)
    _wait_system(page)
    target = page.evaluate(CAMERA_NOW)["target"]
    page.evaluate(fly, [target, 2e-5])
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=30000)
    _settle(page)
    landed = page.evaluate(CAMERA_NOW)["dist"]
    assert landed < 2e-4, landed
    page.wait_for_timeout(1200)
    assert abs(page.evaluate(CAMERA_NOW)["dist"] - landed) < 1e-6 * landed


def test_the_zoom_slider_is_one_log_scale_and_takes_the_view_through_every_level_and_back(page, map_site):
    """MAP.333: the slider stands for the camera's distance on one scale from the whole galaxy to a moon, the same wherever the view
    is, and a drag carries the view across every level (the windows of 81 to 1, the sector, the system) as the wheel does."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, map_site)
    page.evaluate("([c]) => document.querySelector('#galaxymap3d-canvas').galaxyCentreOn(...c)", [list(STAR_FIELD_CENTER_PC)])
    _settle(page)
    seen = {}
    for share in (0.05, 0.15, 0.2, 0.25):
        camera = _slider_to(page, share)
        asked = _slider_distance(share)
        assert abs(math.log(camera["dist"] / asked)) < 0.05, (share, asked, camera["dist"])
        seen[share] = page.evaluate(STAGE_M)
    # MAP.337: the zoom level is the camera's distance alone; each share down the slider is the same level or a finer one, down to single sectors.
    levels = [seen[share] or 10 ** 6 for share in (0.05, 0.15, 0.2, 0.25)]
    assert levels == sorted(levels, reverse=True) and levels[0] > levels[-1] and seen[0.25] == 1, seen
    # The slider follows the camera it is not holding: the wheel moves it to where the camera's distance stands on the scale.
    shown = page.evaluate("() => Number(document.querySelector('#galaxymap3d-zoom-slider').value) / 1000")
    assert abs(shown - 0.25) < 0.02, shown
    assert abs(_slider_to(page, 0.0)["dist"] / 37202.0 - 1) < 0.2
    assert page.evaluate(STAGE_M) is None, "all the way out is the whole galaxy"
    # Down to a sector and a system opened from a star, and out again.
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_selector(".starmap-sr-list button", state="attached")
    page.locator(".starmap-sr-list button", has_text="Middling Sun").first.evaluate("b => b.click()")
    _settle(page)
    _slider_to(page, 0.6)
    assert (_galaxy_system(page) or {}).get("id"), "the slider carried the view into the system"
    camera = _slider_to(page, 1.0)
    assert camera["dist"] < 1e-9, camera["dist"]
    _slider_to(page, 0.3)
    assert not (_galaxy_system(page) or {}).get("id"), "and back out of it"
    assert _slider_to(page, 0.0)["dist"] > 1e4


def test_the_zoom_range_is_one_absolute_range_in_every_stage_and_the_nearest_thing_sets_the_step(page, map_site):
    """MAP.334: what the camera may do does not depend on the stage: the farthest is the whole galaxy in view in a window, a sector and a
    system alike, the nearest 1e-3 pc outside a system (until the scene has the precision for finer) and 1e-12 pc in one; and the wheel
    reads the nearest thing along the pointer's ray, with the distance to it."""
    page.set_viewport_size(VIEWPORT)
    limits_js = "() => document.querySelector('#galaxymap3d-canvas').galaxyZoomLimits()"
    _open_galaxy(page, map_site)
    page.evaluate("([c]) => document.querySelector('#galaxymap3d-canvas').galaxyCentreOn(...c)", [list(STAR_FIELD_CENTER_PC)])
    _settle(page)
    seen = {}
    for share in (0.15, 0.2, 0.25):
        _slider_to(page, share)
        limits = page.evaluate(limits_js)
        seen[page.evaluate(STAGE_M)] = limits
        assert abs(limits["farthestPc"] / 37202.0 - 1) < 0.01, limits
        assert limits["nearestPc"] <= 1e-3 * (1 + 1e-6), limits
    assert 1 in seen and len(seen) >= 2, seen
    # The nearest thing along a ray: at the middle of the map, with the eye's distance to it.
    box = page.locator(GALAXY_CANVAS).bounding_box()
    thing = page.evaluate("([x, y]) => document.querySelector('#galaxymap3d-canvas').galaxyNearestThing(x, y)",
                          [box["x"] + box["width"] / 2, box["y"] + box["height"] / 2])
    eye = page.evaluate(CAMERA_NOW)["eye"]
    assert thing["kind"] in ("star", "pick", "column", "plane", "body"), thing
    assert abs(thing["distance"] - math.dist(eye, thing["point"])) < 1e-6 * thing["distance"], thing
    # In a system the nearest is a moon framed, the farthest still the whole galaxy.
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_selector(".starmap-sr-list button", state="attached")
    page.locator(".starmap-sr-list button", has_text="Middling Sun").first.evaluate("b => b.click()")
    _settle(page)
    _slider_to(page, 0.6)
    assert (_galaxy_system(page) or {}).get("id"), "the slider carried the view into the system"
    limits = page.evaluate(limits_js)
    assert abs(limits["nearestPc"] / 1e-12 - 1) < 0.01 and abs(limits["farthestPc"] / 37202.0 - 1) < 0.01, limits


def test_a_system_is_drawn_from_the_eyes_own_frame_at_the_nearest_view(page, map_site):
    """MAP.338: at the nearest view (1e-12 pc, 30 km) the camera is held as its target plus an exact offset, not as one galactic position (whose
    step is 1.8e-12 pc at 8 kpc), and the open system is drawn in a second pass from the eye's own frame: the offset is the distance asked for
    to a part in a million, the pass runs, and the canvas has the system's picture in it (it is not blank)."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_selector(".starmap-sr-list button", state="attached")
    page.locator(".starmap-sr-list button", has_text="Middling Sun").first.evaluate("b => b.click()")
    _settle(page)
    _slider_to(page, 0.6)
    _wait_system(page)
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyDeepFrame().active") is False
    camera = _slider_to(page, 1.0, seconds=9.0)
    assert abs(camera["dist"] / 1e-12 - 1) < 0.01, camera["dist"]
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyDeepFrame().passes > 0", timeout=20000)
    deep = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyDeepFrame()")
    assert deep["active"] is True and abs(deep["dist"] / camera["dist"] - 1) < 1e-6, (deep, camera)
    assert 0.1 <= deep["dist"] * deep["scale"] <= 10, deep
    shot = Image.open(io.BytesIO(page.locator(GALAXY_CANVAS).screenshot())).convert("L")
    low, high = shot.getextrema()
    lit = high - low
    assert lit > 20, "the nearest view drew nothing"
    _slider_to(page, 0.6)
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyDeepFrame().active") is False


def test_the_wheel_does_not_stop_at_the_end_of_a_window_whose_sector_was_picked(page, map_site):
    """MAP.333: scrolling out of the window of sectors carries on to the larger windows however it was reached: the wheel
    at the end of the range goes to the parent window (it used to let go of a picked sector first and stop there)."""
    page.set_viewport_size(VIEWPORT)
    _open_galaxy(page, map_site)
    page.evaluate("([c]) => document.querySelector('#galaxymap3d-canvas').galaxyCentreOn(...c)", [list(STAR_FIELD_CENTER_PC)])
    _settle(page)
    _slider_to(page, 0.25)
    assert page.evaluate(STAGE_M) == 1
    box = page.locator(GALAXY_CANVAS).bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    assert _wheel_until(page, 400, lambda: page.evaluate(STAGE_M) != 1, ticks=240), "the wheel stopped at the end of the sector window"



def test_the_wheel_heads_for_the_crowd_in_front_of_the_zoom_and_does_not_fly_through_it(page, dense_site):
    """MAP.301: looking at a point behind the dense region (the camera 40 pc over it, its focus 60 pc under it), a wheel scroll with the
    pointer on the region takes the camera into the region and ends in it, where before it flew through to a star 13 pc beyond; a
    pointer off the region has no crowd to head for."""
    page.set_viewport_size(VIEWPORT)
    x, y, z = DENSE_CENTER_PC
    _open_galaxy(page, dense_site, f"?cam={x:.6f},{y:.6f},{z - 60:.6f},100,0,0,0,1")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyStarFill().count > 0", timeout=20000)
    page.evaluate("() => document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    _settle(page)
    box = page.locator(GALAXY_CANVAS).bounding_box()
    middle = (box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    thing = page.evaluate("([a, b]) => document.querySelector('#galaxymap3d-canvas').galaxyNearestThing(a, b)", list(middle))
    assert thing["kind"] == "crowd" and math.dist(thing["point"], DENSE_CENTER_PC) < 2.0, thing
    aside = page.evaluate("([a, b]) => document.querySelector('#galaxymap3d-canvas').galaxyNearestThing(a, b)", [middle[0] + 150, middle[1]])
    assert aside["kind"] != "crowd", aside

    page.mouse.move(*middle)
    came = float("inf")
    for _ in range(30):
        page.mouse.wheel(0, -250)
        page.wait_for_timeout(120)
        eye = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")["eye"]
        came = min(came, math.dist(eye, DENSE_CENTER_PC))
    # The camera comes into the region and rests in it: before the crowd was looked for it went on through to a star 13 pc beyond the
    # middle. Where it rests within the region depends on the frames the approach took, so the test holds the region, not a point.
    assert came < DENSE_HALO_RADIUS_PC / 2, f"the camera came into the region ({came:.2f} pc from its middle)"
    assert math.dist(eye, DENSE_CENTER_PC) < DENSE_HALO_RADIUS_PC, f"and rests in it, not beyond ({eye})"


SECTOR_KIND_OF_PHENOMENON = {
    "nebula": "nebula", "neutron_star": "neutronStar", "rogue_planet": "roguePlanet", "black_hole": "blackHole",
    "asteroid_field": "asteroidField", "supernova_remnant": "supernovaRemnant", "interstellar_comet": "interstellarComet", "quasar": "quasar",
}
"""The kind the Sector Map gives each type of phenomenon the fixture sector holds (`sectorscene.js`'s `kindOf`)."""

HOVER_AT = """([x, y]) => {
    const canvas = document.querySelector("#galaxymap3d-canvas");
    const tip = document.querySelector("#galaxymap3d-tooltip");
    canvas.dispatchEvent(new PointerEvent("pointermove", {clientX: x, clientY: y, bubbles: true, pointerType: "mouse"}));
    return tip.hidden ? "" : tip.textContent;
}"""


def _picks_by_hover(page, entry):
    """Whether a hover at, or a few pixels round, the screen point where `entry` is drawn names it; returns that point or None."""
    spot = page.evaluate("([x, y, z]) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(x, y, z)", [entry["x"], entry["y"], entry["z"]])
    if not spot:
        return None
    for dx, dy in [(0, 0)] + [(a, b) for a in (-4, 0, 4) for b in (-4, 0, 4) if a or b]:
        point = [spot["x"] + dx, spot["y"] + dy]
        if page.evaluate("([x, y]) => document.elementFromPoint(x, y) === document.querySelector('#galaxymap3d-canvas')", point) and \
                entry["name"] in page.evaluate(HOVER_AT, point):
            return point
    return None


def test_the_opened_sector_draws_and_picks_one_of_every_kind_it_holds(page, map_site):
    """MAP.302: at sector zoom (the sector opened in place) every kind of thing the fixture sector holds is in the scene, drawn on the
    map, named when the pointer is over it and picked by a click: stars (a binary and a brown dwarf among them), a nebula, a supernova
    remnant, an asteroid field, a black hole, a neutron star, a quasar, rogue planets and an interstellar comet. Rogue planets are the one kind
    the Menu leaves off until asked for (they outnumber the stars), so the page is opened with them shown."""
    page.add_init_script(NO_WALL_ICONS)
    page.set_viewport_size(VIEWPORT)
    _open(page, f"{map_site}/galaxy?sector={PRIME_CODE}&open=1&show=roguePlanet", GALAXY_CANVAS)
    page.wait_for_function("() => { const open = document.querySelector('#galaxymap3d-canvas').galaxyOpenSector();"
                           " return !!open && open.share >= 0.999; }", timeout=30000)
    _settle(page)
    entries = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxySectorEntries()")
    held = {"star"} | {SECTOR_KIND_OF_PHENOMENON[p["type"]] for p in PHENOMENA}
    assert held <= {entry["kind"] for entry in entries}, f"the scene lacks {sorted(held - {e['kind'] for e in entries})}"
    assert not [e for e in entries if e["hidden"]], "with the rogue planets shown (they are off until asked for) nothing is left off the map"
    names = {entry["name"] for entry in entries}
    assert {name for _id, name, *_rest in SYSTEMS} | {p["name"] for p in PHENOMENA} <= names

    box = page.locator(GALAXY_CANVAS).bounding_box()
    page.evaluate("() => document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    missed = []
    for kind in sorted(held):
        of_kind = [entry for entry in entries if entry["kind"] == kind and not entry["hidden"]]
        found = None
        for entry in of_kind:
            point = _picks_by_hover(page, entry)
            if point:
                found = (entry, point)
                break
        if not found:
            missed.append(kind)
            continue
        entry, point = found
        page.mouse.click(*point)
        title = _info_title(page)
        assert title == entry["name"], f"a click on the {kind} {entry['name']!r} picked {title!r}"
        page.mouse.click(box["x"] + 3, box["y"] + 3)  # empty space between picks (it keeps the selection; the next pick replaces it)
    assert not missed, f"no hover named a {missed} where it is drawn"
