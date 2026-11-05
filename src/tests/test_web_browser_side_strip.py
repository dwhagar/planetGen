"""
MAP.375: the Galaxy Map's linked side strip for height, in a real browser on the dense fixture galaxy of `map_site_support.py`.

- The Menu's Side strip button opens an edge-on strip under the map with the wedge's stars and a readout, and closes it again.
- A star picked on the map is ringed in the strip and the readout gives its exact height, and a star picked in the strip is
  selected on the map (its panel shows), the link in both directions.
- On a phone the strip is a drawer fixed to the bottom of the window, opened from the Menu and put away with its Close button.

Skipped without Playwright, Chromium or WebGL, like the other map browser tests.
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_browser_fixture_maps import (  # noqa: E402,F401 -- fixtures
    GALAXY_CANVAS,
    VIEWPORT,
    _click_star,
    _settle,
    dense_site,
    gl_browser,
    open_dense_region,
    page,
)

STRIP = "#galaxymap3d-strip"
STRIP_CANVAS = "#galaxymap3d-strip canvas"
BUTTON = "#galaxymap3d-controls [data-action=side-strip]"
SNAPSHOT = "() => document.querySelector('#galaxymap3d-canvas').galaxySideStrip()"


def _open_strip(page):
    page.locator("#galaxymap3d-menu > summary").click()
    page.locator(BUTTON).click()
    page.wait_for_selector(STRIP, state="visible")
    page.wait_for_function("() => { const c = document.querySelector('#galaxymap3d-strip canvas'); return c && Number(c.dataset.stars) > 0; }",
                           timeout=20000)


def _star_with_room(stars):
    return next(s for s in stars if all(abs(s["x"] - o["x"]) + abs(s["y"] - o["y"]) > 30 for o in stars if o is not s))


@pytest.mark.names_hidden
def test_the_side_strip_opens_from_the_menu_and_a_pick_in_either_view_selects_in_both(page, dense_site):
    page.set_viewport_size(VIEWPORT)
    open_dense_region(page, dense_site, 14.0)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyPickableStars().length > 20", timeout=20000)
    assert page.locator(STRIP).is_hidden(), "the strip is an option: closed until asked for"
    assert page.evaluate(SNAPSHOT) is None

    _open_strip(page)
    assert page.locator(BUTTON).get_attribute("aria-pressed") == "true"
    snapshot = page.evaluate(SNAPSHOT)
    assert snapshot["on"] and snapshot["wedge"]["widthDeg"] == 60, snapshot["wedge"]
    assert len(snapshot["stars"]) > 20
    assert "Click a star in the strip" in page.locator(f"{STRIP} .galaxy-strip-readout").inner_text()

    # The map to the strip: a star clicked on the map is ringed in the strip, its exact height read out.
    stars = [s for s in page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyPickableStars()") if s["visible"]]
    chosen = _star_with_room([s for s in stars if page.evaluate("(p) => document.elementFromPoint(p[0], p[1]) && document.elementFromPoint(p[0], p[1]).id === 'galaxymap3d-canvas'",
                                                                  [page.locator(GALAXY_CANVAS).bounding_box()["x"] + s["x"], page.locator(GALAXY_CANVAS).bounding_box()["y"] + s["y"]])])
    _click_star(page, chosen)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info').dataset.state === 'selected'")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-strip canvas').dataset.selected === 'ringed'", timeout=10000)
    readout = page.locator(f"{STRIP} .galaxy-strip-readout").inner_text()
    assert "outside this wedge" not in readout, "a selection outside the wedge turns the wedge to face it"
    assert readout.startswith("Height ") and " pc" in readout and "from the center" in readout, readout

    # The strip to the map: with the map's selection put away (a click on bare map), a star clicked in the strip is selected on the map.
    mapbox = page.locator(GALAXY_CANVAS).bounding_box()
    gaps = [(mapbox["width"] * fx / 20, mapbox["height"] * fy / 20) for fx in range(2, 19) for fy in range(2, 19)
            if all(abs(mapbox["width"] * fx / 20 - o["x"]) + abs(mapbox["height"] * fy / 20 - o["y"]) > 40 for o in stars)]
    bare = page.evaluate("""(spots) => spots.filter(([x, y]) => {
        const found = document.elementFromPoint(x, y);
        return found && found.id === 'galaxymap3d-canvas';
    })""", [[mapbox["x"] + x, mapbox["y"] + y] for x, y in gaps])
    assert bare, ("no spot of the map is free of its panels", mapbox, len(gaps), len(stars))
    page.mouse.click(*bare[0])
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info').dataset.state !== 'selected'")
    page.wait_for_function("() => document.querySelector('#galaxymap3d-strip canvas').dataset.selected === 'none'", timeout=10000)
    assert "Click a star in the strip" in page.locator(f"{STRIP} .galaxy-strip-readout").inner_text()

    page.locator(STRIP_CANVAS).scroll_into_view_if_needed()
    box = page.locator(STRIP_CANVAS).bounding_box()
    snapshot = page.evaluate(SNAPSHOT)
    # The dense cluster is a few pixels across in a strip that spans the whole wedge, so click a star with a pixel to itself when
    # there is one (a click picks the nearest star, the brightest of those on the same pixel).
    visible = [s for s in snapshot["stars"] if s["visible"]]
    def room(s):
        return min((abs(s["x"] - o["x"]) + abs(s["y"] - o["y"]) for o in snapshot["stars"] if o is not s), default=99)
    other = max(visible, key=room)
    page.mouse.click(box["x"] + other["x"], box["y"] + other["y"])
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info').dataset.state === 'selected'", timeout=10000)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-strip canvas').dataset.selected === 'ringed'", timeout=10000)
    assert page.locator(f"{STRIP} .galaxy-strip-readout").inner_text().startswith("Height ")

    # Closing: the Close button puts the strip away and the Menu's button follows it.
    page.locator(f"{STRIP} [data-strip-close]").click()
    assert page.locator(STRIP).is_hidden()
    assert page.locator(BUTTON).get_attribute("aria-pressed") == "false"


@pytest.mark.names_hidden
@pytest.mark.viewport(390, 800)
def test_on_a_phone_the_side_strip_is_a_drawer_at_the_bottom_of_the_window(page, dense_site):
    open_dense_region(page, dense_site, 14.0)
    page.locator("#galaxymap3d-menu > summary").click()
    page.locator(BUTTON).click()
    page.wait_for_selector(STRIP, state="visible")
    page.wait_for_function("() => { const c = document.querySelector('#galaxymap3d-strip canvas'); return c && Number(c.dataset.stars) > 0; }",
                           timeout=20000)
    drawer = page.evaluate("""() => {
        const el = document.querySelector('#galaxymap3d-strip');
        const box = el.getBoundingClientRect();
        return { position: getComputedStyle(el).position, bottom: box.bottom, left: box.left, right: box.right,
                 menuOpen: document.querySelector('#galaxymap3d-menu').open, innerHeight: window.innerHeight, innerWidth: window.innerWidth };
    }""")
    assert drawer["position"] == "fixed", drawer
    assert abs(drawer["bottom"] - drawer["innerHeight"]) < 1.5, drawer
    assert drawer["left"] <= 1 and drawer["right"] >= drawer["innerWidth"] - 1, drawer
    assert drawer["menuOpen"] is False, "the Menu folds away so the drawer is not under it"
    # The strip is drawn at the phone's lower star count and still picks.
    snapshot = page.evaluate(SNAPSHOT)
    assert snapshot["on"] and len(snapshot["stars"]) <= 2000
    page.locator(f"{STRIP} [data-strip-close]").click()
    assert page.locator(STRIP).is_hidden()
