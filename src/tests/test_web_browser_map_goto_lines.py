"""
MAP.349: every line of a Galaxy Map card that names a container or a place has its own Go to button, the map's one Go to: the region and its
family on the card of anything inside a charted region (and on the region's own card), each system of a nearest-systems list, and each member
of a crowd. A browser test clicks each and the camera lands on it.

Skipped without Playwright, Chromium or WebGL, like the other map browser tests.
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_browser_fixture_maps import (  # noqa: E402,F401 -- fixtures
    GALAXY_CANVAS, LABELS, MARKERS, OPEN_PRIME, REGION_SPOTS, _expand_card, _markers_until, _on_the_canvas, _open_galaxy, _open_habitable_card, _open_menu,
    _open_window, _poi_layers, _screen_of, _settle, gl_browser, map_site, page,
)

CAMERA = "() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()"


def _landed(page):
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=30000)
    _settle(page)


def test_a_region_card_has_its_region_and_family_lines_each_with_a_go_to(page, map_site):
    """MAP.349: the card of a charted region (its sphere picked) has the "Region:" and "In:" lines like any card inside it, each with a Go to;
    the family's takes the camera to its 81-block."""
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
    page.mouse.move(where[0] + 1, where[1])
    page.mouse.move(*where)
    page.wait_for_timeout(300)
    page.mouse.click(*where)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    _expand_card(page)
    lines = page.locator("#galaxymap3d-info .map-info-region .map-info-place")
    lines.first.wait_for(state="visible", timeout=20000)
    texts = lines.all_inner_texts()
    assert texts[0].startswith("Region: ") and texts[1].startswith("In: "), texts
    assert all(lines.nth(i).locator(".map-info-go").count() == 1 for i in range(len(texts)))
    before = page.evaluate(CAMERA)["dist"]
    lines.nth(1).locator(".map-info-go").evaluate("b => b.click()")
    _landed(page)
    assert page.evaluate(CAMERA)["dist"] != pytest.approx(before, rel=1e-3), "the camera went to the family"


def test_each_system_of_a_nearest_list_has_a_go_to_that_opens_it(page, map_site):
    """MAP.349: a system card's nearest-systems list gives each system a Go to beside its link, which opens that system's sector and the
    system in place (what the find box does with its printed ID). The fixture's neighbours of a marker have no printed ID, so the list the
    server sends is replaced by one naming a system by a printed ID in the fixture's sector."""
    code = OPEN_PRIME.split("sector=")[1].split("&")[0]
    near = {"nearest": [{"id": f"{code}-A2C-000", "name": "Fixture Neighbour", "distance_ly": 3.2, "url": "/system/9002"}]}
    page.route("**/galaxy/system/*/nearest", lambda route: route.fulfill(json=near))
    _open_habitable_card(page, map_site, "Fixture Eden")
    _expand_card(page)
    items = page.locator("#galaxymap3d-info .map-info-nearest-list li")
    items.first.wait_for()
    assert all(items.nth(i).locator(".map-info-go").count() == 1 for i in range(items.count())), items.all_inner_texts()
    items.first.locator(".map-info-go").evaluate("b => b.click()")
    page.wait_for_function("() => location.search.includes('open=1')", timeout=30000)
    _landed(page)
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyOpenSector()"), "the sector is open in the map"


def test_each_member_of_a_crowd_has_a_go_to(page, map_site):
    """MAP.349: a crowd's card lists its systems, each with a Go to of its own that flies to it."""
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_timeout(1000)
    ids = ("marker-inhabited", "marker-habitable", "marker-uninhabitable", "marker-gas-giant", "marker-asteroid-field")
    _poi_layers(page, *ids)
    found = _markers_until(page, lambda entries: any(e["kind"] in ids and e["count"] for e in entries))
    crowd = max((e for e in found if e["kind"] in ids and e["count"]), key=lambda e: e["count"])
    page.locator("#galaxymap3d-menu > summary").click()
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    spot = _screen_of(page, (crowd["x"], crowd["y"], crowd["z"]))
    assert spot
    page.mouse.move(*spot)
    page.mouse.move(spot[0] + 1, spot[1])
    page.wait_for_timeout(200)
    page.mouse.click(*spot)
    page.wait_for_function("() => document.querySelector('#galaxymap3d-info h3')")
    _expand_card(page)
    lines = page.locator("#galaxymap3d-info .map-info-crowd li")
    assert lines.count() >= 2
    assert all(lines.nth(i).locator(".map-info-go").count() == 1 for i in range(lines.count()))
    before = page.evaluate(CAMERA)["target"]
    lines.first.locator(".map-info-go").evaluate("b => b.click()")
    _landed(page)
    after = page.evaluate(CAMERA)["target"]
    assert max(abs(a - b) for a, b in zip(before, after)) > 0.5, (before, after)
