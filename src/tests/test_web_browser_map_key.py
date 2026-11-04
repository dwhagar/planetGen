"""
MAP.234: the Galaxy Map's key. The help dialog has a Key section; the map has an on-map Key that names what is drawn at the view, a line for
each kind of thing while it is drawn, that folds away and stays folded; and every layer the map registers has a line, a legend or a reason it
has none, so a new layer cannot be added without one.

Skipped without Playwright, Chromium or WebGL, like the other map browser tests.
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_browser_fixture_maps import (  # noqa: E402,F401 -- fixtures
    OPEN_PRIME, _open_galaxy, _open_sector, _settle, gl_browser, map_site, page,
)

KEY = "#galaxymap3d-layer-legend"
ROWS = f"{KEY} .map-layer-legend-line"
LAYERS = "document.querySelector('#galaxymap3d-canvas').__mapLayers"


def _drawn_ids(page):
    return page.evaluate(f"() => {LAYERS}.keys().filter((k) => k.drawn).map((k) => k.id)")


def _row_ids(page):
    return page.eval_on_selector_all(ROWS, "els => els.map((e) => e.dataset.layer)")


def test_every_layer_of_the_galaxy_map_has_a_key_line_a_legend_or_a_reason(page, map_site):
    """MAP.234: `unkeyed()` names the layers of the map's registry that the key says nothing about; it is empty on the Galaxy Map and on
    a sector opened in it, so a new layer (or a new kind of thing in a sector) cannot be added without a line."""
    _open_galaxy(page, map_site)
    assert page.evaluate(f"() => {LAYERS}.unkeyed()") == []
    _open_galaxy(page, map_site, OPEN_PRIME)
    assert page.evaluate(f"() => {LAYERS}.unkeyed()") == []


def test_the_key_has_a_line_for_each_kind_drawn_at_the_view_and_none_for_what_is_not(page, map_site):
    """MAP.234: the lines of the on-map Key are the key lines whose thing is drawn at the view (the stars, the blocks of the galaxy, the
    regions), none for a course that is not shown, and an opened sector adds its own kinds."""
    _open_galaxy(page, map_site)
    page.wait_for_function(f"() => {LAYERS}.keys().some((k) => k.drawn)", timeout=20000)
    page.wait_for_timeout(1000)    # the key follows the view every 400 ms
    drawn = _drawn_ids(page)
    assert "key-stars" in drawn, drawn
    assert "key-course" not in drawn and "key-direct" not in drawn, drawn
    rows = _row_ids(page)
    assert set(drawn) <= set(rows), (drawn, rows)
    assert "key-course" not in rows, rows
    _open_galaxy(page, map_site, OPEN_PRIME)
    page.wait_for_function(f"() => {LAYERS}.keys().some((k) => k.drawn && k.id === 'key-stars')", timeout=20000)
    page.wait_for_timeout(1000)
    rows = _row_ids(page)
    assert "key-stars" in rows, rows
    assert set(_drawn_ids(page)) <= set(rows), (_drawn_ids(page), rows)
    text = page.locator(KEY).inner_text()
    assert "stars" in text.lower(), text


def test_the_key_folds_away_and_stays_folded(page, map_site):
    """MAP.234: the Key is a fold the reader can close; its lines are then out of sight, and a fresh visit finds it closed."""
    _open_galaxy(page, map_site)
    page.wait_for_function(f"() => {LAYERS}.keys().some((k) => k.drawn)", timeout=20000)
    page.wait_for_timeout(1000)
    fold = page.locator(f"{KEY} details.map-key")
    assert fold.evaluate("d => d.open") is True
    assert page.locator(ROWS).first.is_visible()
    page.locator(f"{KEY} details.map-key > summary").click()
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-layer-legend details.map-key').open")
    assert not page.locator(ROWS).first.is_visible()
    _open_galaxy(page, map_site)
    page.wait_for_function(f"() => {LAYERS}.keys().some((k) => k.drawn)", timeout=20000)
    assert page.locator(f"{KEY} details.map-key").evaluate("d => d.open") is False
    page.locator(f"{KEY} details.map-key > summary").click()
    page.wait_for_function("() => document.querySelector('#galaxymap3d-layer-legend details.map-key').open")
    assert page.locator(ROWS).first.is_visible()


def test_the_help_dialog_has_a_key_section_naming_what_the_map_draws(page, map_site):
    """MAP.234: the help dialog's Key section names the blocks, stars, regions, superbubbles, nebulae and remnants, the marks, an opened
    sector's kinds and the course lines."""
    _open_galaxy(page, map_site)
    text = page.locator("#galaxymap3d-help").text_content()
    heading = page.locator("#galaxymap3d-help h3", has_text="Key")
    assert heading.count() == 1
    section = text[text.index("Key"):]
    for word in ("Blocks", "Points of light", "Soft spheres", "superbubbles", "nebulae", "supernova remnants", "Marks and symbols", "In an opened sector", "Course"):
        assert word in section, word


def test_the_wheel_over_the_key_zooms_the_map(page, map_site):
    """MAP.234: the Key sits over the map, so a wheel turned over its fold zooms the map as it does anywhere else."""
    _open_galaxy(page, map_site)
    page.wait_for_function(f"() => {LAYERS}.keys().some((k) => k.drawn)", timeout=20000)
    box = page.locator(f"{KEY} details.map-key > summary").bounding_box()
    before = page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera().zoom")
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    page.mouse.wheel(0, -400)
    page.wait_for_timeout(600)
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera().zoom") != before


def test_the_scrolling_key_lines_can_be_reached_from_the_keyboard_and_have_a_name(page, map_site):
    """TEST.519: the lines scroll inside their box, so the box is a tab stop and a named group (axe: scrollable-region-focusable)."""
    _open_galaxy(page, map_site)
    page.wait_for_function(f"() => document.querySelector('{KEY} .map-key-lines')", timeout=20000)
    box = page.locator(f"{KEY} .map-key-lines")
    assert box.get_attribute("tabindex") == "0"
    assert box.get_attribute("role") == "group" and box.get_attribute("aria-label")
