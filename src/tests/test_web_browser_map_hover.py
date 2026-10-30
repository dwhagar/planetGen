"""
MAP.348: hovering a star, a system, a region, a nebula or a black hole on the Galaxy Map shows a small preview: the card's first view (its
name, what kind of thing it is and up to three values), in a corner of the map that is clear of the pointer, of the hovered thing and of the
open card. It takes no pointer event, so a click still reaches the map, and a touch screen shows none.

Skipped without Playwright, Chromium or WebGL, like the other map browser tests.
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_browser_fixture_maps import (  # noqa: E402,F401 -- fixtures
    GALAXY_CANVAS, MARKERS, REGION_SPOTS, _info_title, _markers_until, _on_the_canvas, _open_galaxy, _open_menu, _open_sector, _open_star_field,
    _open_window, _poi_layers, _screen_of, _settle, gl_browser, map_site, page,
)

PREVIEW = ".map-hover-preview"

STATE = """() => {
    const node = document.querySelector(".map-hover-preview");
    if (!node || node.hidden) return null;
    const own = node.getBoundingClientRect();
    const host = node.parentElement.getBoundingClientRect();
    const card = document.querySelector("#galaxymap3d-info");
    const cardBox = card && card.dataset.state === "selected" ? card.getBoundingClientRect() : null;
    const overlap = (a, b) => a.left < b.right && a.right > b.left && a.top < b.bottom && a.bottom > b.top;
    return {
        title: node.querySelector("h3").textContent,
        kind: (node.querySelector(".map-info-kind") || {}).textContent || null,
        rows: node.querySelectorAll("dl dt").length,
        rect: {left: own.left, right: own.right, top: own.top, bottom: own.bottom},
        inHost: own.left >= host.left - 1 && own.right <= host.right + 1 && own.top >= host.top - 1 && own.bottom <= host.bottom + 1,
        overCard: !!cardBox && overlap(own, cardBox),
        pointerEvents: getComputedStyle(node).pointerEvents,
        hit: (() => { const e = document.elementFromPoint((own.left + own.right) / 2, (own.top + own.bottom) / 2); return !!e && node.contains(e); })(),
    };
}"""

OFFSETS = ((0, 0), (0, -6), (6, -6), (-6, -6), (0, -9), (8, 0), (-8, 0), (0, 6), (0, 9), (12, 0), (-12, 0))


def _hover_for_preview(page, spot, title=None):
    """Moves the pointer about `spot` (a name label over a mark takes the pointer first) until the preview is up; returns where the pointer is
    and what the preview says."""
    for dx, dy in OFFSETS:
        at = (spot[0] + dx, spot[1] + dy)
        page.mouse.move(at[0] - 30, at[1] - 30)
        page.mouse.move(*at)
        try:
            page.wait_for_function(
                "t => { const n = document.querySelector('.map-hover-preview'); return !!n && !n.hidden && (!t || n.querySelector('h3').textContent === t); }",
                arg=title, timeout=1500)
        except sync_api.TimeoutError:
            continue
        return at, page.evaluate(STATE)
    raise AssertionError(f"no pixel near {spot} shows the preview of {title!r}")


def _check_preview(page, at, state, title=None):
    assert state, "the preview is up"
    if title:
        assert state["title"] == title, state
    assert state["kind"], state
    assert 1 <= state["rows"] <= 3, state
    assert state["inHost"] and not state["overCard"], state
    assert state["pointerEvents"] == "none" and not state["hit"], state
    rect = state["rect"]
    assert not (rect["left"] <= at[0] <= rect["right"] and rect["top"] <= at[1] <= rect["bottom"]), (at, state)


def _click_through(page, at, title):
    """The preview takes no click: pressing where the pointer is opens the card of what is there, and the preview goes."""
    page.mouse.click(*at)
    page.wait_for_function("t => document.querySelector('#galaxymap3d-info h3') && document.querySelector('#galaxymap3d-info h3').textContent === t",
                           arg=title, timeout=5000)
    assert page.evaluate(STATE) is None, "the press puts the preview away"


def _marker(page, name):
    found = _markers_until(page, lambda entries: any(e["name"] == name and not e["count"] for e in entries))
    entry = next(e for e in found if e["name"] == name and not e["count"])
    spot = _screen_of(page, (entry["x"], entry["y"], entry["z"]))
    assert spot, entry
    return spot


def test_a_system_hovered_shows_its_first_view_and_the_click_goes_through(page, map_site):
    """MAP.348: a star system's mark: the preview names it, says it is a star system and gives up to three values; the card opens on a click."""
    page.set_viewport_size({"width": 1280, "height": 900})
    _open_star_field(page, map_site)
    page.wait_for_timeout(500)
    _poi_layers(page, "marker-habitable")
    _markers_until(page, lambda found: any(e["kind"] == "marker-habitable" for e in found))
    page.locator("#galaxymap3d-menu > summary").click()
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    at, state = _hover_for_preview(page, _marker(page, "Fixture Eden"), "Fixture Eden")
    _check_preview(page, at, state, "Fixture Eden")
    assert state["kind"] == "Star system", state
    _click_through(page, at, "Fixture Eden")


def test_a_region_hovered_shows_its_first_view(page, map_site):
    """MAP.348: a charted region's sphere: the preview names it and says it is a region."""
    page.set_viewport_size({"width": 1280, "height": 900})
    _open_galaxy(page, map_site)
    page.wait_for_function("() => { const g = document.querySelector('#galaxymap3d-canvas').__mapRegions; return g && g.drawn() > 0; }", timeout=20000)
    _poi_layers(page)
    page.locator("#galaxymap3d-menu > summary").click()
    page.evaluate("() => window.scrollTo(0, 0)")
    box = page.locator(GALAXY_CANVAS).bounding_box()
    spots = [s for s in page.evaluate(REGION_SPOTS) if s["visible"] and _on_the_canvas(page, (box["x"] + s["x"], box["y"] + s["y"]))]
    assert spots
    where = (box["x"] + spots[0]["x"], box["y"] + spots[0]["y"])
    at, state = _hover_for_preview(page, where)
    _check_preview(page, at, state)
    assert state["kind"] == "Region", state
    _click_through(page, at, state["title"])


def test_a_black_hole_hovered_shows_its_first_view(page, map_site):
    """MAP.348: a black hole's mark."""
    page.set_viewport_size({"width": 1280, "height": 900})
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_timeout(500)
    _poi_layers(page, "marker-black-hole")
    page.locator("#galaxymap3d-menu > summary").click()
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    at, state = _hover_for_preview(page, _marker(page, "Fixture Maw"), "Fixture Maw")
    _check_preview(page, at, state, "Fixture Maw")
    assert "Black Hole" in state["kind"], state


def _sector_thing(page, kind):
    entries = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxySectorEntries()")
    assert entries, "an opened sector"
    return [e for e in entries if e["kind"] == kind and not e["hidden"]]


@pytest.mark.parametrize("kind", ["star", "nebula"])
def test_a_star_and_a_nebula_of_an_opened_sector_show_their_first_view(page, map_site, kind):
    """MAP.348: the things of an opened sector: a star, and a nebula."""
    page.set_viewport_size({"width": 1280, "height": 900})
    _open_sector(page, map_site)
    _settle(page)
    things = _sector_thing(page, kind)
    assert things, kind
    shown = None
    for thing in things[:12]:
        spot = _screen_of(page, (thing["x"], thing["y"], thing["z"]))
        if not spot or not _on_the_canvas(page, spot):
            continue
        try:
            shown = _hover_for_preview(page, spot)
            break
        except AssertionError:
            continue
    assert shown, f"no {kind} of the sector shows a preview"
    at, state = shown
    _check_preview(page, at, state)


def test_a_touch_shows_no_preview(page, map_site):
    """MAP.348: a touch screen has no hover: a touch pointer over a thing shows nothing extra."""
    page.set_viewport_size({"width": 1280, "height": 900})
    _open_sector(page, map_site)
    _settle(page)
    star = _sector_thing(page, "star")[0]
    spot = _screen_of(page, (star["x"], star["y"], star["z"]))
    page.evaluate("""([x, y]) => {
        const canvas = document.querySelector('#galaxymap3d-canvas');
        canvas.dispatchEvent(new PointerEvent('pointermove', {pointerType: 'touch', clientX: x, clientY: y, bubbles: true}));
    }""", list(spot))
    page.wait_for_timeout(600)
    assert page.evaluate(STATE) is None
