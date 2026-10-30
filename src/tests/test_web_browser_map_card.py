"""
MAP.347: the Galaxy Map's info popup is a card with a first view that fits in a quarter of the map's height without scrolling (its name, what
kind of thing it is, three values and its buttons) and, behind an Expand button, the rest of what it says in a region that scrolls on its own,
never the whole card. Checked for a star, a system, a region, a nebula and a black hole, at 1080p and at phone width.

Skipped without Playwright, Chromium or WebGL, like the other map browser tests.
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_browser_fixture_maps import (  # noqa: E402,F401 -- fixtures
    GALAXY_CANVAS, REGION_SPOTS, _info_title, _on_the_canvas, _open_galaxy, _open_habitable_card, _open_menu, _open_sector, _open_window, _poi_layers,
    _select_the_nebula, _settle, gl_browser, map_site, page,
)

SIZES = [(1920, 1080), (390, 844)]

CARD = """() => {
    const info = document.querySelector("#galaxymap3d-info");
    const view = document.querySelector("#galaxymap3d-canvas").closest(".starmap-viewport").getBoundingClientRect();
    const box = info.getBoundingClientRect();
    const rest = info.querySelector(".map-info-rest");
    const inside = (r) => r.left >= box.left - 1 && r.right <= box.right + 1 && r.top >= box.top - 1 && r.bottom <= box.bottom + 1;
    const first = [...info.children].filter((e) => e !== rest && e.getClientRects().length);
    const expand = info.querySelector(".map-info-expand");
    return {
        title: info.querySelector("h3").textContent, kind: (info.querySelector(".map-info-kind") || {}).textContent || null,
        rows: info.querySelectorAll(":scope > dl dt").length,
        height: box.height, view: view.height,
        withinMap: box.left >= view.left - 1 && box.right <= view.right + 1 && box.top >= view.top - 1 && box.bottom <= view.bottom + 1,
        cardScrolls: info.scrollHeight > info.clientHeight + 1 || info.scrollWidth > info.clientWidth + 1,
        firstVisible: first.every((e) => inside(e.getBoundingClientRect())),
        buttons: [...info.querySelectorAll(":scope > .map-info-actions > *")].filter((e) => e.getClientRects().length).map((e) => [e.textContent, inside(e.getBoundingClientRect())]),
        expandable: !!expand && !expand.hidden, expanded: info.hasAttribute("data-expanded"),
        restHidden: !rest || rest.hidden,
        restScrolls: !!rest && !rest.hidden && rest.scrollHeight > rest.clientHeight + 1,
        restOverflow: rest ? getComputedStyle(rest).overflowY : null,
    };
}"""


def _star(page, base):
    _open_sector(page, base)
    page.locator(".starmap-sr-list button").first.dispatch_event("click")


def _system(page, base, width, height):
    _open_habitable_card(page, base, "Fixture Eden", width, height)
    page.wait_for_timeout(500)


def _hole(page, base):
    _open_window(page, base, 27, 27, 0, 0, g=3)
    page.wait_for_timeout(500)
    page.locator(".starmap-sr-markers button", has_text="Fixture Maw").evaluate("b => b.click()")


def _region(page, base):
    _open_galaxy(page, base)
    page.wait_for_function("() => { const g = document.querySelector('#galaxymap3d-canvas').__mapRegions; return g && g.drawn() > 0; }", timeout=20000)
    _poi_layers(page)
    _open_menu(page)
    page.locator("#galaxymap3d-menu > summary").click()
    page.evaluate("() => window.scrollTo(0, 0)")    # the clicks above scrolled the page to the Menu's buttons
    box = page.locator(GALAXY_CANVAS).bounding_box()
    spots = [s for s in page.evaluate(REGION_SPOTS) if s["visible"] and _on_the_canvas(page, (box["x"] + s["x"], box["y"] + s["y"]))]
    assert spots, page.evaluate(REGION_SPOTS)
    where = (box["x"] + spots[0]["x"], box["y"] + spots[0]["y"])
    page.mouse.move(*where)
    page.wait_for_timeout(300)
    page.mouse.click(*where)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")


KINDS = {
    "star": lambda page, base, width, height: _star(page, base),
    "system": _system,
    "region": lambda page, base, width, height: _region(page, base),
    "nebula": lambda page, base, width, height: (_open_sector(page, base), _select_the_nebula(page)),
    "black hole": lambda page, base, width, height: _hole(page, base),
}


@pytest.mark.parametrize("kind", list(KINDS))
@pytest.mark.parametrize("width,height", SIZES)
def test_the_card_shows_a_first_view_in_a_quarter_of_the_map_and_expands_into_a_region_that_scrolls(page, map_site, kind, width, height):
    """MAP.347: for each kind of thing at 1080p and phone width, the first view (name, kind, up to three values, the buttons) is fully visible
    in at most a quarter of the map's height (half at phone width), the card never scrolls, and Expand opens the rest inside the card in a region that scrolls on
    its own."""
    page.set_viewport_size({"width": width, "height": height})
    KINDS[kind](page, map_site, width, height)
    assert _info_title(page)
    page.wait_for_timeout(800)
    first = page.evaluate(CARD)
    assert first["withinMap"] and first["firstVisible"] and not first["cardScrolls"], first
    assert first["kind"], first
    assert first["rows"] <= 3, first
    # A quarter of the map's height at 1080p; a phone's card holds the same lines in a narrower column (its values wrap), so half of it there.
    assert first["height"] <= first["view"] / (4 if width >= 1000 else 2) + 1, first
    assert all(inside for _, inside in first["buttons"]) and first["buttons"], first
    assert first["restHidden"] and not first["expanded"], first
    if not first["expandable"]:
        return
    page.locator("#galaxymap3d-info .map-info-expand").click()
    open_ = page.evaluate(CARD)
    assert open_["expanded"] and not open_["restHidden"], open_
    assert open_["withinMap"] and not open_["cardScrolls"] and open_["firstVisible"], open_
    assert open_["restOverflow"] == "auto", open_
    assert all(inside for _, inside in open_["buttons"]), open_
    page.locator("#galaxymap3d-info .map-info-expand").click()
    assert page.evaluate(CARD)["restHidden"]


LONG_CARD = """() => {
    const panel = document.querySelector("#galaxymap3d-canvas").__mapInfo;
    panel.show({title: "A system with a very long name that goes on and on and on", kind: "Star system",
        fields: Array.from({length: 40}, (_, n) => ["Row " + n, "A value with a good deal of text so that it has to wrap onto a second line " + n]),
        buttons: [{label: "Go to it", onClick: () => {}}]});
}"""


@pytest.mark.parametrize("width,height", [(1280, 720), (800, 560), (390, 640)])
def test_a_card_far_too_long_for_the_window_scrolls_what_is_behind_expand_and_keeps_its_name_and_buttons(page, map_site, width, height):
    """MAP.347: a card of 40 long rows, expanded, stays inside the map's window, scrolls its rest region and not itself, and its name and
    buttons are in view whatever the rest has scrolled to."""
    page.set_viewport_size({"width": width, "height": height})
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_timeout(500)
    page.evaluate(LONG_CARD)
    page.locator("#galaxymap3d-info .map-info-expand").click()
    got = page.evaluate(CARD)
    assert got["expanded"] and got["withinMap"] and not got["cardScrolls"] and got["restScrolls"], got
    page.evaluate("() => { const r = document.querySelector('#galaxymap3d-info .map-info-rest'); r.scrollTop = r.scrollHeight; }")
    after = page.evaluate(CARD)
    assert after["withinMap"] and not after["cardScrolls"] and after["firstVisible"], after
    assert all(inside for _, inside in after["buttons"]), after
