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

import io
import re
from urllib.parse import parse_qs, urlparse

import pytest

sync_api = pytest.importorskip("playwright.sync_api")
Image = pytest.importorskip("PIL.Image")

from tests.map_site_support import PHENOMENA, SECTOR_ID, SECTORS, SYSTEMS, map_site  # noqa: E402,F401 -- fixture

VIEWPORT = {"width": 1280, "height": 900}
GL_ARGS = ["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"]
SECTOR_CANVAS = "#starmap-canvas"
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
def page(gl_browser):
    context = gl_browser.new_context(viewport=VIEWPORT, device_scale_factor=1, reduced_motion="reduce",
                                     color_scheme="dark")
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    try:
        yield page
    finally:
        context.close()
    assert not errors, f"uncaught exceptions: {errors}"


def _open(page, url, ready):
    response = page.goto(url, wait_until="load")
    assert response.status == 200, f"{url} answered {response.status}"
    page.wait_for_selector(ready, state="attached")
    if not page.evaluate("!!document.createElement('canvas').getContext('webgl2')"):
        pytest.skip("no WebGL in this Chromium")
    page.wait_for_timeout(300)


def _shot(page, selector):
    page.wait_for_timeout(150)
    return page.locator(selector).screenshot()


# --- The Sector Map -------------------------------------------------------------------

def _open_sector(page, base, query=""):
    _open(page, f"{base}/sector/{SECTOR_ID}{query}", SECTOR_CANVAS)
    page.wait_for_selector(".starmap-sr-list button", state="attached")


def _info_title(page):
    heading = page.locator("#starmap-info h3")
    return heading.inner_text() if heading.count() else None


def _bright_spots(page, selector, limit=40, white=False):
    """Screen points (page coordinates) of the brightest small spots on a
    canvas, brightest first: where its points of light are drawn (with
    `white`, only near-white ones: a star's whitened core, not the Galaxy
    Map's lit blocks)."""
    image = Image.open(io.BytesIO(_shot(page, selector))).convert("RGB")
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
    info = page.locator("#starmap-info")
    assert "Star type" in info.inner_text()
    view = info.locator("a", has_text="View system")
    assert view.count() == 1 and view.get_attribute("href").endswith(f"/system/{system_id}")
    nav = [a.inner_text() for a in info.locator("a").all()]
    assert "Nav from here" in nav and "Nav to here" in nav
    assert _shot(page, SECTOR_CANVAS) != before, "the picked star gets a highlight ring"


def test_sector_map_click_on_empty_space_keeps_the_selection(page, map_site):
    _open_sector(page, map_site)
    _x, _y, name = _click_a_star(page)
    box = page.locator(SECTOR_CANVAS).bounding_box()
    page.mouse.click(box["x"] + 3, box["y"] + 3)
    assert _info_title(page) == name


def test_sector_map_drag_turns_and_picks_nothing(page, map_site):
    _open_sector(page, map_site)
    x, y, name = _click_a_star(page)
    page.locator("#starmap-info").evaluate("panel => panel.textContent = ''")
    before = _shot(page, SECTOR_CANVAS)
    page.mouse.move(x, y)
    page.mouse.down()
    page.mouse.move(x + 60, y + 15, steps=6)
    page.mouse.up()
    assert _info_title(page) is None, "a drag that starts on a star doesn't pick it"
    assert _shot(page, SECTOR_CANVAS) != before, "dragging turns the map"


def test_sector_map_arrow_keys_turn_the_view_and_keep_the_scale(page, map_site):
    _open_sector(page, map_site)
    canvas = page.locator(SECTOR_CANVAS)
    scale = page.locator("#starmap-scale-label").inner_text()
    canvas.focus()
    for key in ("ArrowLeft", "ArrowRight", "ArrowUp", "ArrowDown"):
        before = _shot(page, SECTOR_CANVAS)
        page.keyboard.press(key)
        assert _shot(page, SECTOR_CANVAS) != before, f"{key} didn't turn the Sector Map"
    assert page.locator("#starmap-scale-label").inner_text() == scale, "turning doesn't change the scale"
    before = _shot(page, SECTOR_CANVAS)
    page.keyboard.press("a")
    assert _shot(page, SECTOR_CANVAS) == before, "other keys do nothing"


def _scale_line(page):
    width = page.locator("#starmap-scale-bar").evaluate("bar => parseFloat(bar.style.width)")
    return page.locator("#starmap-scale-label").inner_text(), width


def _ly(label):
    number = float(re.match(r"[\d.,]+", label).group(0).replace(",", ""))
    return number / 1000 if "mly" in label else number


def test_sector_map_scale_line_follows_the_zoom(page, map_site):
    _open_sector(page, map_site)
    label, width = _scale_line(page)
    assert re.fullmatch(r"[\d.,]+ (ly|mly)", label), label
    assert 30 <= width <= 160, "the scale line is about 70 px long"
    ly_per_px = _ly(label) / width

    zoom_in = page.locator('#starmap-controls [data-action="zoom-in"]')
    zoom_out = page.locator('#starmap-controls [data-action="zoom-out"]')
    zoom_in.click()
    in_label, in_width = _scale_line(page)
    assert _ly(in_label) / in_width < ly_per_px, "+ shows fewer light-years per pixel"

    box = page.locator(SECTOR_CANVAS).bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    page.mouse.wheel(0, 300)
    page.wait_for_timeout(100)
    wheel_label, wheel_width = _scale_line(page)
    assert _ly(wheel_label) / wheel_width > _ly(in_label) / in_width, "scrolling down zooms out"

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

    page.locator('#starmap-controls [data-action="reset"]').click()
    assert _scale_line(page) == (label, width), "Reset view goes back to the opening zoom"


def test_sector_map_show_on_map_selects_a_rogue_planet(page, map_site):
    _open_sector(page, map_site)
    rogue = next(p for p in PHENOMENA if p["type"] == "rogue_planet")
    button = page.locator(f'[data-map-target="rogue_planet:{rogue["id"]}"]')
    assert button.count() == 1
    # The Contents row lists rogue planets in a folded group.
    button.evaluate("b => { const d = b.closest('details'); if (d) d.open = true; }")
    button.click()
    assert _info_title(page) == rogue["name"]
    assert page.evaluate("document.activeElement.id") == "starmap-canvas"
    view = page.locator("#starmap-info a", has_text="View phenomenon")
    assert view.get_attribute("href").endswith(f"/phenomenon/rogue_planet/{rogue['id']}")


def test_sector_map_rogue_planet_markers_toggle(page, map_site):
    _open_sector(page, map_site)
    toggle = page.locator('#starmap-controls [data-action="toggle-rogue-markers"]')
    background = toggle.evaluate("b => getComputedStyle(b).backgroundColor")
    assert toggle.get_attribute("aria-pressed") == "false", "off by default (MAP.83)"
    before = _shot(page, SECTOR_CANVAS)
    toggle.click()
    assert toggle.get_attribute("aria-pressed") == "true"
    assert toggle.evaluate("b => getComputedStyle(b).backgroundColor") != background, "highlighted while on"
    assert _shot(page, SECTOR_CANVAS) != before
    toggle.click()
    assert toggle.get_attribute("aria-pressed") == "false"
    assert toggle.evaluate("b => getComputedStyle(b).backgroundColor") == background


def test_sector_map_screen_reader_list_selects(page, map_site):
    _open_sector(page, map_site)
    names = page.locator(".starmap-sr-list button").all_inner_texts()
    assert SYSTEM_NAMES <= set(names)
    assert {p["name"] for p in PHENOMENA} <= set(names)
    page.locator(".starmap-sr-list button", has_text="Fixture Pulsar").evaluate("b => b.click()")
    assert _info_title(page) == "Fixture Pulsar"


# --- The Galaxy Map -------------------------------------------------------------------

def _crumbs(page):
    """Every step to here, galaxy first: the Steps menu's list, which
    holds them all however much of the breadcrumb line is folded into
    "…" (MAP.93, MAP.94)."""
    return page.eval_on_selector_all("#galaxymap3d-steps [data-steps-panel] li",
                                     "els => els.map(e => e.textContent.trim())")


def _query(page):
    return urlparse(page.url).query


def _settle(page):
    """Waits out a stage's flight: until the scale bar (its width follows
    the camera) holds still, a few seconds at most on a busy machine."""
    page.wait_for_timeout(250)
    last = None
    for _ in range(40):
        now = page.evaluate("() => (document.querySelector('#galaxymap3d-scale') || {}).innerHTML || ''")
        if now == last:
            return
        last = now
        page.wait_for_timeout(150)


def _open_galaxy(page, base, query=""):
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
            canvas.dispatchEvent(new PointerEvent("pointermove", {clientX: x, clientY: y, bubbles: true, pointerType: "mouse"}));
            const text = tip.hidden ? "" : tip.textContent;
            if (text && (!re || re.test(text))) return [x, y, text];
        }
    }
    return null;
}"""


HOVER_RETRIES = 20


def _hover_choice(page, wanted=None, steps=32):
    """Moves the pointer over the map (as pointer events, for speed) until
    the tooltip names a choice matching `wanted` (a regex, JavaScript
    syntax); returns (x, y, text) or None, the pointer left there. The
    map ignores the pointer while a stage's flight runs, which takes longer
    on a busy machine, so the scan is tried again for a few seconds."""
    found = page.evaluate(HOVER_SCAN, [wanted, steps])
    for _ in range(HOVER_RETRIES):
        if found:
            break
        page.wait_for_timeout(250)
        found = page.evaluate(HOVER_SCAN, [wanted, steps])
    if found:
        page.mouse.move(found[0], found[1])
    return found


def _click_choice(page, wanted=None):
    found = _hover_choice(page, wanted)
    assert found, f"nothing to click on the Galaxy Map ({wanted}); at {_crumbs(page)}"
    x, y, text = found
    before = (page.url, _crumbs(page))
    if re.match(r"Sector .*, generated$", text):
        # A generated sector opens its own page.
        with page.expect_navigation():
            page.mouse.click(x, y)
        return text
    page.mouse.click(x, y)
    _settle(page)
    assert (page.url, _crumbs(page)) != before, f"clicking {text!r} did nothing"
    return text


def _on_galaxy(page):
    return urlparse(page.url).path == "/galaxy"


def _label_of(tooltip):
    """The breadcrumb label a tooltip's choice leads to (see
    test_web_browser_maps.py)."""
    return re.split(r"[:,] ", tooltip, maxsplit=1)[0].strip()


GENERATED_CHOICE = r"(^|, )([1-9][\d,.]*\S* (of .+ )?sectors )?generated$"
"""A tooltip for a choice holding a generated sector (the fixture's
sectors, all near the core)."""


def _walk_down(page, stages=16):
    """Clicks generated choices down from the galaxy, up to `stages`
    clicks or until one opens a generated sector's own page (the last
    entry then); returns the stages seen, [(query, crumbs)]."""
    steps = [(_query(page), _crumbs(page))]
    for _ in range(stages):
        text = _click_choice(page, GENERATED_CHOICE)
        if not _on_galaxy(page):
            steps.append((page.url, None))
            return steps
        steps.append((_query(page), _crumbs(page)))
        assert _crumbs(page)[-1] == _label_of(text), f"clicked {text!r}, breadcrumb now {_crumbs(page)}"
    return steps


def test_galaxy_map_hover_names_a_choice_and_leaving_hides_it(page, map_site):
    _open_galaxy(page, map_site)
    found = _hover_choice(page)
    assert found, "hovering the whole galaxy shows no choice"
    assert re.match(r"(Quarter|Block|Arc)\b", found[2]) or found[2], found[2]
    box = page.locator(GALAXY_CANVAS).bounding_box()
    page.mouse.move(box["x"] + box["width"] + 200, box["y"] + 10)
    page.wait_for_timeout(100)
    assert not page.locator("#galaxymap3d-tooltip").is_visible(), "leaving the map hides the tooltip"


def test_galaxy_map_clicks_walk_down_to_a_generated_sector(page, map_site):
    _open_galaxy(page, map_site)
    assert _crumbs(page) == ["Galaxy"] and _query(page) == ""
    steps = _walk_down(page)
    url, crumbs = steps[-1]
    assert crumbs is None and re.search(r"/sector/\d+$", url), f"never opened a sector: {steps[-2:]}"
    assert int(url.rsplit("/", 1)[1]) in {sector[0] for sector in SECTORS}
    for query, _crumbs_seen in steps[1:-1]:
        assert re.fullmatch(r"(at=\d+\.\d+\.\d+\.-?\d+)?(&?p=[asL0-9.~,-]+)?", query), query
    # MAP.56: arc, slab, segment (entering each block), slab, ... down to
    # a sector; never a 3 x 3 region.
    kinds = " ".join(q for q, _ in steps[:-1])
    for marker in ("p=a", "L", "at=243.", "at=27.", "at=3."):
        assert marker in kinds, f"no {marker} stage on the way down: {[q for q, _ in steps]}"
    assert not re.search(r"[=,]r\d", kinds), kinds


def test_galaxy_map_back_forward_and_url_state(page, map_site):
    _open_galaxy(page, map_site)
    steps = _walk_down(page, stages=3)
    queries = [q for q, _ in steps]
    page.go_back()
    _settle(page)
    assert _query(page) in queries[:-1]
    at = queries.index(_query(page))
    assert _crumbs(page) == steps[at][1]
    page.go_forward()
    _settle(page)
    assert _query(page) in queries[at + 1:]

    # The map's own Back, Forward and Up.
    here, crumbs = _query(page), _crumbs(page)
    page.click('#galaxymap3d-controls [data-action="back"]')
    _settle(page)
    assert _query(page) != here
    page.click('#galaxymap3d-controls [data-action="forward"]')
    _settle(page)
    assert _query(page) == here
    page.click('#galaxymap3d-controls [data-action="up"]')
    _settle(page)
    assert _crumbs(page) == crumbs[:-1]

    # A stage's URL opens it directly; Whole galaxy goes home.
    deep_query, deep_crumbs = steps[-1]
    _open_galaxy(page, map_site, "?" + deep_query)
    assert _crumbs(page) == deep_crumbs
    page.click('#galaxymap3d-controls [data-action="reset"]')
    _settle(page)
    assert _crumbs(page) == ["Galaxy"] and _query(page) == ""


def test_galaxy_map_keys_pick_go_up_and_home(page, map_site):
    _open_galaxy(page, map_site)
    canvas = page.locator(GALAXY_CANVAS)
    canvas.focus()
    tooltip = page.locator("#galaxymap3d-tooltip")
    page.keyboard.press("ArrowRight")
    assert tooltip.is_visible(), "an arrow key marks a choice"
    first = tooltip.inner_text()
    page.keyboard.press("ArrowRight")
    assert tooltip.inner_text() != first, "the next arrow moves to the next choice"
    page.keyboard.press("Enter")
    _settle(page)
    assert len(_crumbs(page)) == 2, "Enter takes the marked choice"
    one_down = _crumbs(page)
    page.keyboard.press("ArrowRight")
    page.keyboard.press("Enter")
    _settle(page)
    if len(_crumbs(page)) == 2:  # a slab choice first: Enter on a layer
        page.keyboard.press("ArrowRight")
        page.keyboard.press("Enter")
        _settle(page)
    assert len(_crumbs(page)) > 2
    page.keyboard.press("Escape")
    _settle(page)
    assert len(_crumbs(page)) < len(one_down) + 2
    page.keyboard.press("Home")
    _settle(page)
    assert _crumbs(page) == ["Galaxy"]


SLAB_TIP = r"^(Slab|Layer)s? -?\d+( to -?\d+)?: "
"""A tooltip naming a whole slab (or a layer of sectors), not one block."""


def test_galaxy_map_hover_while_picking_a_slab_lights_the_whole_slab(page, map_site):
    """MAP.91: whenever the next pick is a slab, at every level and in the
    3 by 3 by 3 view of sectors too, hovering any cube names (and the
    slab's button and line light up for) its whole slab, and a click picks
    that slab."""
    _open_galaxy(page, map_site)
    slab_stages = []
    for _ in range(16):
        if not _on_galaxy(page):
            break
        if page.locator(".galaxy-slab-button").count():
            heading = page.locator("#galaxymap3d-slabs-heading").inner_text()
            found = _hover_choice(page, GENERATED_CHOICE)
            assert found, f"no slab to hover at {_crumbs(page)}"
            text = found[2]
            assert re.match(SLAB_TIP, text), f"hovering a {heading} pick named {text!r}, not a slab"
            lit = page.locator(".galaxy-slab-button.is-lit .galaxy-slab-name")
            assert lit.count() == 1, "one slab button lights"
            assert lit.inner_text() == _label_of(text), f"the buttons mark {lit.inner_text()!r}, the map {text!r}"
            assert page.locator(".galaxy-slab-leader.is-lit").count() == 1
            # Every cube of the slab names that same slab.
            names = set()
            for fx in (0.3, 0.5, 0.7):
                for fy in (0.3, 0.5, 0.7):
                    box = page.locator(GALAXY_CANVAS).bounding_box()
                    page.mouse.move(box["x"] + fx * box["width"], box["y"] + fy * box["height"])
                    tip = page.locator("#galaxymap3d-tooltip")
                    if tip.is_visible():
                        names.add(tip.inner_text())
            assert all(re.match(SLAB_TIP, name) for name in names), names
            page.mouse.move(found[0], found[1])
            slab_stages.append(heading)
            label = _label_of(text)
            page.mouse.click(found[0], found[1])
            _settle(page)
            assert _crumbs(page)[-1] == label, f"clicking {text!r} led to {_crumbs(page)}"
            continue
        _click_choice(page, GENERATED_CHOICE)
    assert "Layers" in slab_stages, f"never met the 3 by 3 by 3 view of sectors: {slab_stages}"
    assert any(h != "Layers" for h in slab_stages), f"never met a slab pick of blocks: {slab_stages}"


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
    _click_choice(page, GENERATED_CHOICE)
    assert not reset_view.is_disabled(), f"no free stage by {_crumbs(page)}"
    stage_scale = _galaxy_scale(page)
    assert stage_scale != whole, "a closer stage has a shorter scale"
    box = page.locator(GALAXY_CANVAS).bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    page.mouse.wheel(0, -400)
    _settle(page)
    assert _galaxy_scale(page) != stage_scale, "the wheel zooms and the scale follows"
    page.click('#galaxymap3d-menu summary')
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
    _click_choice(page)
    stage_query = _query(page)
    star = page.locator("#galaxymap3d-crumbs .galaxy-bookmark")
    assert star.count() == 1
    star.click()
    saved = _bookmarks(page)
    assert len(saved) == 1 and saved[0]["kind"] == "stage" and stage_query in saved[0]["value"]
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
    _click_choice(page)
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
    info = page.locator("#starmap-info")
    links = [a.inner_text() for a in info.locator("a").all()]
    assert links == ["Use as destination"] or (len(links) == 1 and "destination" in links[0].lower()), links
    page.locator(".starmap-sr-list button", has_text="Fixture Pulsar").evaluate("b => b.click()")
    links = [a.inner_text() for a in info.locator("a").all()]
    assert len(links) == 1 and "View" not in links[0], links


def _click_a_galaxy_star(page):
    url = page.url
    for x, y in _bright_spots(page, GALAXY_CANVAS, limit=40, white=True):
        page.mouse.click(x, y)
        heading = page.locator("#galaxymap3d-info h3")
        if heading.count() and heading.inner_text() == "Bright star":
            return
        if page.url != url:  # missed, and took a block instead
            page.go_back()
            _settle(page)
    pytest.fail("no click on the Galaxy Map picked a star")


@pytest.mark.parametrize("pick", [False, True])
def test_galaxy_map_star_links_to_its_system_except_while_picking(page, map_site, pick):
    _open_galaxy(page, map_site, "?at=27.27.0.0" + ("&pick=to" if pick else ""))
    page.wait_for_timeout(1000)
    _click_a_galaxy_star(page)
    links = page.locator("#galaxymap3d-info a", has_text="View system")
    assert links.count() == (0 if pick else 1)


# --- NAV.40: bookmarks keep a course pick ---------------------------------------------

SEED_BOOKMARKS = """(entries) => {
    const db = document.querySelector("[data-bookmark-db]").getAttribute("data-bookmark-db");
    localStorage.setItem("planetgen.bookmarks." + db, JSON.stringify(entries));
}"""

PICK_BOOKMARKS = [
    {"name": "Saved view", "kind": "stage", "value": "/galaxy?at=27.27.0.0", "created": 1},
    {"name": "Middling Sun", "kind": "system", "value": "system:702", "url": "/system/702", "created": 2},
    {"name": "Fixture Second", "kind": "sector", "value": "0·0·2", "url": "/sector/8", "sectorId": 8, "created": 3},
]


def _menu_links(page, menu):
    menu.locator("summary").click()
    return dict(menu.locator("[data-bookmarks-panel] a").evaluate_all(
        "els => els.map(a => [a.querySelector('.bookmarks-name').textContent, a.getAttribute('href')])"))


def _assert_pick_links(links):
    assert links["Middling Sun"] == "/nav?from=system:701&to=system:702", links
    assert links["Fixture Second"] == "/sector/8?pick=to&from=system:701", links
    view = parse_qs(urlparse(links["Saved view"]).query)
    assert view == {"at": ["27.27.0.0"], "pick": ["to"], "from": ["system:701"]}, links


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
    params = parse_qs(_query(page))
    assert params["pick"] == ["to"] and params["from"] == ["system:701"] and params["at"] == ["27.27.0.0"]
    assert page.locator(".pick-banner").count() == 1
    # A step on the map keeps the pick in the URL.
    page.click('#galaxymap3d-controls [data-action="up"]')
    _settle(page)
    assert "at" not in parse_qs(_query(page)) or parse_qs(_query(page))["at"] != ["27.27.0.0"]
    params = parse_qs(_query(page))
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
            more: !!ol.querySelector(".galaxy-crumb-menu"), visible: ol.offsetParent !== null,
            starTop: Math.round(document.querySelector("#galaxymap3d-crumbs .galaxy-bookmark").getBoundingClientRect().top),
            olTop: Math.round(ol.getBoundingClientRect().top), olHeight: ol.getBoundingClientRect().height};
}"""


def _open_deep_galaxy(page, map_site):
    """The Galaxy Map at the sector level of a fixture sector (many steps)."""
    _open_sector(page, map_site)
    href = page.locator(".badges a", has_text="Show on Galaxy Map").get_attribute("href")
    _open_galaxy(page, map_site, href[len("/galaxy"):])
    return _crumbs(page)


def test_galaxy_breadcrumb_stays_on_one_line_folding_its_middle(page, map_site):
    """MAP.93: at any width the breadcrumb is one line: the first step,
    "…" for the steps that don't fit, the last ones and the current one;
    it fits again on resize, and "…" lists the folded steps."""
    page.set_viewport_size({"width": 1280, "height": 900})
    full = _open_deep_galaxy(page, map_site)
    assert len(full) >= 8, full
    for width in (1280, 900, 700, 600):
        page.set_viewport_size({"width": width, "height": 900})
        page.wait_for_timeout(200)
        line = page.evaluate(CRUMB_LINE)
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
    more = page.locator("#galaxymap3d-crumbs .galaxy-crumb-menu")
    more.locator("summary").click()
    folded = more.locator("li").all_inner_texts()
    line_items = page.locator("#galaxymap3d-crumbs > ol > li").all_inner_texts()
    kept = line_items[2:]
    assert line_items[0] == "Galaxy" and folded, line_items
    assert [line_items[0]] + folded + kept == full, (line_items, folded, full)
    more.locator("button").first.click()
    _settle(page)
    assert _crumbs(page) == full[:2], "the first folded step goes there"
    assert not page.evaluate(CRUMB_LINE)["more"] or page.evaluate(CRUMB_LINE)["overflow"] <= 1
    # Wider again: everything shown, no "…".
    page.set_viewport_size({"width": 1280, "height": 900})
    page.wait_for_timeout(200)
    assert not page.evaluate(CRUMB_LINE)["more"] or page.evaluate(CRUMB_LINE)["overflow"] <= 1


def test_galaxy_phone_steps_button_between_back_and_forward(page, map_site):
    """MAP.94: at phone width the breadcrumb line gives way to a round
    Steps button between Back and Forward, listing every step with the
    current one marked; Reset stays beside the arrows."""
    page.set_viewport_size({"width": 390, "height": 844})
    full = _open_deep_galaxy(page, map_site)
    line = page.evaluate(CRUMB_LINE)
    assert not line["visible"], "no breadcrumb line on a phone"
    order = page.locator("#galaxymap3d-controls > *").evaluate_all(
        "els => els.map(e => e.dataset.action || e.id || e.className)")
    assert order.index("back") + 1 == order.index("galaxymap3d-steps") == order.index("forward") - 1, order
    assert page.locator('#galaxymap3d-controls [data-action="reset"]').is_visible()
    steps = page.locator("#galaxymap3d-steps")
    box = steps.locator("summary").bounding_box()
    assert abs(box["width"] - box["height"]) < 2, "round"
    steps.locator("summary").click()
    items = steps.locator("li").all_inner_texts()
    assert items == full, (items, full)
    assert steps.locator("[aria-current]").inner_text() == full[-1]
    steps.locator("button", has_text="Galaxy").first.click()
    _settle(page)
    assert _crumbs(page) == ["Galaxy"]
    assert not steps.evaluate("d => d.open"), "taking a step closes the menu"
    page.set_viewport_size({"width": 1280, "height": 900})
    page.wait_for_timeout(200)
    assert not steps.is_visible() and page.evaluate(CRUMB_LINE)["visible"]


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
            if (!tip.hidden) seen.add(tip.textContent.replace(/ \\(nothing generated here to pick\\)$/, ""));
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
    later stage the choice) under the pointer just as outside pick mode,
    generated or not; only an empty one can't be taken."""
    _open_galaxy(page, map_site)
    _click_choice(page, GENERATED_CHOICE)
    stages = ["", "?" + _query(page)]
    for query in stages:
        _open_galaxy(page, map_site, query)
        browsing = _choices_seen(page)
        _open_galaxy(page, map_site, query + ("&" if query else "?") + "pick=to&from=system:701")
        picking = _choices_seen(page)
        assert len(browsing) > 1 and picking == browsing, (query, len(picking), len(browsing))
    # An empty arc lights but a click there stays put.
    _open_galaxy(page, map_site, "?pick=to&from=system:701")
    found = _hover_choice(page, r"nothing generated here to pick\)$")
    assert found, "no empty arc to hover"
    before = (page.url, _crumbs(page))
    page.mouse.click(found[0], found[1])
    _settle(page)
    assert (page.url, _crumbs(page)) == before


LEADERS = """() => {
    const canvas = document.querySelector("#galaxymap3d-canvas").getBoundingClientRect();
    const svg = document.querySelector(".galaxy-slab-leaders");
    const base = svg.getBoundingClientRect();
    return Array.from(document.querySelectorAll(".galaxy-slab-button")).map(button => {
        const line = svg.querySelector('.galaxy-slab-leader[data-slab="' + button.dataset.slab + '"]');
        const points = line.getAttribute("points").trim().split(" ").map(p => p.split(",").map(Number));
        const at = points.map(p => [p[0] + base.left, p[1] + base.top]);
        const b = button.getBoundingClientRect();
        return {
            slab: button.dataset.slab,
            points: at,
            start: at[0],
            end: at[at.length - 1],
            button: [b.left, b.top, b.right, b.bottom],
            canvas: [canvas.left, canvas.top, canvas.right, canvas.bottom],
            off: line.classList.contains("is-off"),
        };
    });
}"""


def _cross(a, b, c, d):
    """Whether segments ab and cd cross (touching ends don't count)."""
    def side(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    return (side(a, b, c) * side(a, b, d) < 0) and (side(c, d, a) * side(c, d, b) < 0)


def _check_leaders(leaders):
    assert len(leaders) >= 2, leaders
    for line in leaders:
        x0, y0, x1, y1 = line["button"]
        # Each line leaves its own button ...
        assert x0 - 1 <= line["start"][0] <= x1 + 1 and y0 - 1 <= line["start"][1] <= y1 + 1, line
        # ... and ends on the map.
        cx0, cy0, cx1, cy1 = line["canvas"]
        assert cx0 <= line["end"][0] <= cx1 and cy0 <= line["end"][1] <= cy1, line
    # Ordered by their slabs' height on screen, so no two lines cross.
    for n, a in enumerate(leaders):
        for b in leaders[n + 1:]:
            for p, q in zip(a["points"], a["points"][1:]):
                for r, t in zip(b["points"], b["points"][1:]):
                    assert not _cross(p, q, r, t), (a, b)


def test_galaxy_map_slab_buttons_have_lines_that_follow_the_view(page, map_site):
    """MAP.54 and MAP.76: no slab slider; one button per slab, each with a
    line from it to its slab on the map that is redrawn when the view
    turns; hovering a button lights its slab and line, and clicking it
    picks the slab."""
    _open_galaxy(page, map_site)
    for _ in range(4):
        if page.locator(".galaxy-slab-button").count():
            break
        _click_choice(page, GENERATED_CHOICE)
    assert page.locator(".galaxy-slab-button").count() >= 2, _crumbs(page)
    assert not page.locator("input[type=range]").count(), "the slab slider is gone"
    before = page.evaluate(LEADERS)
    _check_leaders(before)
    # Turning the view moves the lines' ends with their slabs.
    page.locator(GALAXY_CANVAS).scroll_into_view_if_needed()
    box = page.locator(GALAXY_CANVAS).bounding_box()
    page.mouse.move(box["x"] + box["width"] * 0.5, box["y"] + 20)
    page.mouse.down()
    page.mouse.move(box["x"] + box["width"] * 0.7, box["y"] + 40, steps=8)
    page.mouse.up()
    page.wait_for_timeout(200)
    after = page.evaluate(LEADERS)
    _check_leaders(after)
    ends = {line["slab"]: line["end"] for line in before}
    assert any(abs(line["end"][0] - ends[line["slab"]][0]) + abs(line["end"][1] - ends[line["slab"]][1]) > 2
               for line in after), f"the lines didn't follow the turn: {before} {after}"
    # A button lights its slab and line, and picks it.
    button = page.locator(".galaxy-slab-button").first
    label = button.locator(".galaxy-slab-name").inner_text()
    button.hover()
    assert "is-lit" in button.get_attribute("class")
    lit = page.locator(".galaxy-slab-leader.is-lit")
    assert lit.count() == 1 and lit.get_attribute("data-slab") == button.get_attribute("data-slab")
    button.click()
    _settle(page)
    assert _crumbs(page)[-1] == label


def test_galaxy_map_slab_lines_on_a_phone_leave_the_buttons_below_the_map(gl_browser, map_site):
    """MAP.76: at phone width the buttons sit in one column below the map
    and their lines still reach their slabs without crossing."""
    context = gl_browser.new_context(viewport={"width": 390, "height": 844}, device_scale_factor=1,
                                     reduced_motion="reduce")
    page = context.new_page()
    try:
        _open_galaxy(page, map_site)
        for _ in range(4):
            if page.locator(".galaxy-slab-button").count():
                break
            _click_choice(page, GENERATED_CHOICE)
        leaders = page.evaluate(LEADERS)
        _check_leaders(leaders)
        assert all(line["button"][1] >= line["canvas"][3] for line in leaders), "buttons below the map"
        lefts = {round(line["button"][0]) for line in leaders}
        assert len(lefts) == 1, "one column"
    finally:
        context.close()


FRAME = "() => document.querySelector('#galaxymap3d-canvas').galaxyFrame()"


def _framed(frame, least=0.75):
    """The stage's blocks all on the map, and filling much of it one way."""
    assert frame, "no frame"
    assert frame["left"] >= -1 and frame["top"] >= -1, frame
    assert frame["right"] <= frame["width"] + 1 and frame["bottom"] <= frame["height"] + 1, frame
    fill = max((frame["right"] - frame["left"]) / frame["width"], (frame["bottom"] - frame["top"]) / frame["height"])
    assert fill >= least, f"the stage fills only {fill:.0%} of the map: {frame}"
    return fill


@pytest.mark.parametrize("query", ["", "?p=a0.180", "?p=a0.180,L0"])
def test_galaxy_map_frames_the_whole_stage_at_any_size_and_turn(page, map_site, query):
    """MAP.53 and MAP.78: every stage (the galaxy, an arc, a slab of it) is
    framed whole on the map, fitted to the map's own size, so a bigger
    window shows it bigger; it refits when the window is resized, and
    stays whole while Shift and the arrow keys turn it about its middle."""
    _open_galaxy(page, map_site, query)
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
