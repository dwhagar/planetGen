"""
MAP.336 (zoom redesign Z3): one Go to, the find-a-place box and the log ruler on the Galaxy Map, in a real browser on the fixture
galaxy of `map_site_support.py`.

- Every move to a named place is one `goTo`: the box takes coordinates, a sector designation, an object ID, the name of a charted region,
  the names the server finds and the places visited before, and a result is gone to with Enter or a click.
- The ruler is one logarithmic scale for the camera's distance, the whole galaxy at the top and a moon at the bottom; each named stop is a
  Go to for that scale about the middle of the map, and a readout names the distance.
- A Go to a region lands framed on its charted part, or on the cluster of its heaviest sectors when the reader's setting says so.

Skipped without Playwright, Chromium or WebGL, like the other map browser tests.
"""

import json
import math

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.map_site_support import REGION_RADIUS_PC, SECTORS, _sector_center_pc  # noqa: E402
from tests.test_web_browser_fixture_maps import (  # noqa: E402,F401 -- fixtures
    GALAXY_CANVAS,
    OPEN_PRIME,
    PRIME_CODE,
    _camera,
    _crumbs,
    _open_galaxy,
    _raw_query,
    _settle,
    gl_browser,
    map_site,
    page,
)

BOX = "#galaxymap3d-address-input"
RESULTS = "#galaxymap3d-matches"
CANVAS_JS = "document.querySelector('#galaxymap3d-canvas')"


def _results(page):
    return page.eval_on_selector_all(f"{RESULTS} [role=option]", "els => els.map(e => ({text: e.textContent, kind: e.dataset.kind}))")


def _groups(page):
    return page.eval_on_selector_all(f"{RESULTS} .galaxy-find-group h4", "els => els.map(e => e.textContent)")


def _type(page, text):
    page.locator(BOX).fill(text)
    page.wait_for_timeout(450)    # the box waits out typing, and the server's answer


def _go_to(page, text):
    _type(page, text)
    page.locator(BOX).press("Enter")
    _settle(page)
    page.wait_for_function(f"() => !{CANVAS_JS}.galaxyMoving()", timeout=30000)
    return _camera(page)


def test_the_slash_key_reaches_the_find_box_and_coordinates_are_gone_to(page, map_site):
    _open_galaxy(page, map_site)
    page.locator(GALAXY_CANVAS).focus()
    page.keyboard.press("/")
    assert page.evaluate("() => document.activeElement && document.activeElement.id") == BOX[1:]
    _type(page, "3000 0 0")
    assert _groups(page) == ["Go to this place"]
    assert _results(page) == [{"text": "3000, 0, 0 pccoordinates", "kind": "point"}]
    camera = _go_to(page, "3000 0 0")
    assert camera["target"] == pytest.approx([3000.0, 0.0, 0.0], abs=1e-6)
    assert camera["dist"] < 40, camera["dist"]
    assert page.locator(RESULTS).is_hidden(), "the list closes once a place is gone to"
    assert page.locator(BOX).input_value() == ""


def test_a_designation_and_an_object_id_open_their_sector(page, map_site):
    _open_galaxy(page, map_site)
    _go_to(page, PRIME_CODE)
    assert f"sector={PRIME_CODE}" in _raw_query(page)
    assert "open=1" not in _raw_query(page)
    _open_galaxy(page, map_site)
    _type(page, f"moon:{PRIME_CODE}-0000005-003")
    assert [r["kind"] for r in _results(page)][:1] == ["object"]
    page.locator(BOX).press("Enter")
    page.wait_for_function(f"() => new URL(location.href).searchParams.get('open') === '1'", timeout=30000)
    assert f"sector={PRIME_CODE}" in _raw_query(page)


def test_regions_are_found_by_name_and_remembered_as_recent_places(page, map_site):
    _open_galaxy(page, map_site)
    page.wait_for_function(f"() => {CANVAS_JS}.__mapRegions && {CANVAS_JS}.__mapRegions.drawn() > 0", timeout=20000)
    page.locator(BOX).focus()
    page.wait_for_timeout(300)
    assert "Places" in _groups(page) and "Recent" not in _groups(page), "an empty box offers the way to the whole galaxy"
    _type(page, "zarh")
    assert _groups(page) == ["Regions"]
    assert [r["text"] for r in _results(page)][0].startswith("Fixture Zarh Region")
    camera = _go_to(page, "zarh")
    assert math.dist(camera["target"], [28.876, 21.498, 0.0]) < 0.1, camera["target"]
    page.locator(BOX).focus()
    page.wait_for_timeout(300)
    assert _groups(page)[0] == "Recent"
    assert _results(page)[0]["text"].startswith("Fixture Zarh Region")
    # A recent place is gone to again from the whole galaxy.
    page.click('#galaxymap3d-controls [data-action="reset"]')
    _settle(page)
    page.locator(BOX).focus()
    page.wait_for_timeout(300)
    page.locator(f"{RESULTS} [data-kind=recent]").first.click()
    _settle(page)
    assert math.dist(_camera(page)["target"], [28.876, 21.498, 0.0]) < 0.1


def test_names_the_server_finds_are_listed_and_a_system_opens_in_its_sector(page, map_site):
    seen = []

    def answer(route):
        seen.append(route.request.url)
        route.fulfill(content_type="application/json", body=json.dumps({"matches": [{
            "kind": "system", "ref": "system:701", "id": "701", "name": "Middling Sun", "sector_id": SECTORS[0][0],
            "sector_name": "Fixture Prime", "parents": [], "shell_index": 1, "band_index": 0, "wedge_index": 1}]}))

    page.route("**/galaxy/locate*", answer)
    _open_galaxy(page, map_site)
    _type(page, "middling")
    assert _groups(page) == ["Names"]
    assert _results(page) == [{"text": "Middling Sunsystem", "kind": "system"}]
    assert any("q=middling" in url for url in seen)
    page.locator(f"{RESULTS} [role=option]").first.click()
    page.wait_for_function("() => new URL(location.href).searchParams.get('open') === '1'", timeout=30000)
    assert f"sector={PRIME_CODE}" in _raw_query(page)


def test_the_log_ruler_names_its_scales_and_each_stop_is_a_go_to(page, map_site):
    _open_galaxy(page, map_site, OPEN_PRIME)
    stops = page.eval_on_selector_all("#galaxymap3d-ruler .starmap-ruler-stop",
                                      "els => els.map(e => ({key: e.dataset.stop, y: e.getBoundingClientRect().y, title: e.title}))")
    assert [s["key"] for s in stops] == ["galaxy", "arm", "region", "sector", "neighbourhood", "system", "planet", "moon", "surface"]
    ys = [s["y"] for s in stops]
    assert ys == sorted(ys) and ys[0] < ys[-1] - 100, "the whole galaxy is at the top and the surface at the bottom"
    assert all(s["title"].startswith("Go to ") for s in stops)
    readout = page.locator("#galaxymap3d-altitude")
    assert readout.inner_text().startswith("Distance from focus: ")
    before = _camera(page)
    page.locator('#galaxymap3d-ruler [data-stop="sector"]').click()
    _settle(page)
    page.wait_for_function(f"() => !{CANVAS_JS}.galaxyMoving()", timeout=30000)
    at_sector = _camera(page)
    assert at_sector["dist"] == pytest.approx(4.0, rel=0.02), at_sector["dist"]
    assert at_sector["target"] == pytest.approx(before["target"], abs=1e-6), "a stop keeps the middle of the map"
    assert page.locator("#galaxymap3d-zoom-slider").input_value() != "0"
    page.locator('#galaxymap3d-ruler [data-stop="galaxy"]').click()
    _settle(page)
    page.wait_for_function(f"() => !{CANVAS_JS}.galaxyMoving()", timeout=30000)
    far = _camera(page)
    assert far["dist"] > 10000, far["dist"]
    assert page.locator("#galaxymap3d-zoom-slider").input_value() == "0", "the slider is at the top at the whole galaxy"
    page.locator('#galaxymap3d-ruler [data-stop="planet"]').click()
    _settle(page)
    page.wait_for_function(f"() => !{CANVAS_JS}.galaxyMoving()", timeout=30000)
    assert _camera(page)["dist"] == pytest.approx(1e-8, rel=0.02)
    text = readout.inner_text()
    assert text.startswith("Distance from focus: ") and " pc" not in text.split("(")[0], "the readout is in the reader's units, not parsecs only"


def test_the_ruler_readout_names_the_distance_in_parsecs_and_light_years(page, map_site):
    _open_galaxy(page, map_site)
    page.evaluate(f"() => {CANVAS_JS}.galaxyGoTo({{center: [0, 0, 0], dist: 40}})")
    _settle(page)
    page.wait_for_timeout(500)
    text = page.locator("#galaxymap3d-altitude").inner_text()
    assert text.startswith("Distance from focus: 40 pc") and "130 ly" in text, text


def test_the_crumbs_and_up_one_level_are_go_tos_to_what_contains_the_view(page, map_site):
    _open_galaxy(page, map_site, OPEN_PRIME)
    up = page.locator('#galaxymap3d-controls [data-action="up"]')
    assert up.inner_text() == "Up one level"
    crumbs = _crumbs(page)
    assert len(crumbs) >= 3, crumbs
    up.click()
    _settle(page)
    assert len(_crumbs(page)) == len(crumbs) - 1 or "open=1" not in _raw_query(page)
    page.locator("#galaxymap3d-crumbs a, #galaxymap3d-crumbs button").filter(has_text=crumbs[1]).first.click()
    _settle(page)
    assert _crumbs(page)[-1] == crumbs[1]


def _region_target():
    """The Go to target of the fixture's Fixture Beez Region (the 9-block 9.0.0.0): its block and the sphere its extent bounds."""
    from planetgen.galaxy.spherical import sector_position_pc

    centre = sector_position_pc(*SECTORS[0][1:4], 4.0)
    return {"region": {"block": {"m": 9, "shell": 0, "band": 0, "wedge": 0}, "center": list(centre), "radius": REGION_RADIUS_PC}}


def _stable_camera(page):
    """The camera once it has come to rest for good (a landing can go on from the frame to a second flight)."""
    last = None
    for _ in range(30):
        page.wait_for_timeout(300)
        now = page.evaluate(f"() => {CANVAS_JS}.galaxyCamera()")
        if now and not now["flying"] and last and now["dist"] == last["dist"]:
            return now
        last = now
    return last


def test_go_to_a_region_lands_framed_by_default_and_on_the_densest_cluster_by_setting(page, map_site):
    target = _region_target()
    centre = target["region"]["center"]
    _open_galaxy(page, map_site)
    assert page.evaluate(f"(t) => {CANVAS_JS}.galaxyGoTo(t)", target)
    framed = _stable_camera(page)
    assert math.dist(framed["target"], centre) < 0.1, framed["target"]
    assert math.asin(min(1.0, REGION_RADIUS_PC / framed["dist"])) <= math.radians(25) + 1e-6
    assert framed["dist"] > 25, "framed on the whole charted part of the region"

    page.click('#galaxymap3d-controls [data-action="reset"]')
    _settle(page)
    page.evaluate("() => localStorage.setItem('planetgen.pref.map.regionLanding', 'densest')")
    assert page.evaluate(f"(t) => {CANVAS_JS}.galaxyGoTo(t)", target)
    dense = _stable_camera(page)
    assert dense["dist"] < framed["dist"] * 0.6, (dense["dist"], framed["dist"])
    assert math.dist(dense["target"], centre) < REGION_RADIUS_PC, dense["target"]

    # A word that is not one of the two reads as the default.
    page.click('#galaxymap3d-controls [data-action="reset"]')
    _settle(page)
    page.evaluate("() => localStorage.setItem('planetgen.pref.map.regionLanding', 'closest')")
    assert page.evaluate(f"(t) => {CANVAS_JS}.galaxyGoTo(t)", target)
    again = _stable_camera(page)
    assert again["dist"] == pytest.approx(framed["dist"], rel=0.02)


# --- UX.236: the Preferences page sets what the map reads ------------------------------------------------------------------------------

STORED = "(key) => localStorage.getItem('planetgen.pref.' + key)"


def test_a_setting_made_on_the_preferences_page_changes_where_a_go_to_a_region_lands(page, map_site):
    """UX.236: the page lists each setting with its words and its default; choosing Land on the densest part stores map.regionLanding, the
    Galaxy Map then lands a Go to a region on the cluster of its heaviest sectors, and Reset to defaults brings the framed landing back."""
    target = _region_target()
    page.goto(f"{map_site}/preferences", wait_until="load")
    page.wait_for_selector("[data-preference='map.regionLanding'] input[type=radio]", state="attached")
    landing = page.locator("[data-preference='map.regionLanding']")
    assert "Where Go to a region lands" in landing.inner_text() and "Default: Frame the region" in landing.inner_text()
    assert page.locator("[data-preference='map.blackHoleDensity']").inner_text().count("Default: 120 pixels") == 1
    assert page.locator("input[value=frame]").is_checked() and not page.locator("input[value=densest]").is_checked()
    assert page.evaluate(STORED, "map.regionLanding") is None, "a page that only looks stores nothing"

    page.locator("input[value=densest]").check()
    assert page.evaluate(STORED, "map.regionLanding") == "densest"
    assert "Saved" in page.locator("[data-preferences-status]").inner_text()

    _open_galaxy(page, map_site)
    assert page.evaluate(f"(t) => {CANVAS_JS}.galaxyGoTo(t)", target)
    dense = _stable_camera(page)
    page.click('#galaxymap3d-controls [data-action="reset"]')
    _settle(page)

    page.goto(f"{map_site}/preferences", wait_until="load")
    page.wait_for_selector("input[value=densest]", state="attached")
    assert page.locator("input[value=densest]").is_checked(), "the page shows what is stored when it is opened again"
    page.get_by_role("button", name="Reset to defaults").click()
    assert page.evaluate(STORED, "map.regionLanding") is None and page.locator("input[value=frame]").is_checked()

    _open_galaxy(page, map_site)
    assert page.evaluate(f"(t) => {CANVAS_JS}.galaxyGoTo(t)", target)
    framed = _stable_camera(page)
    assert dense["dist"] < framed["dist"] * 0.6, (dense["dist"], framed["dist"])


def test_a_number_on_the_preferences_page_is_stored_when_in_range_and_refused_when_not(page, map_site):
    """UX.236: map.blackHoleDensity takes 30 to 1000; a number outside that is not stored and the page says so, an emptied box is the
    default again, and the other settings stay as they were."""
    page.goto(f"{map_site}/preferences", wait_until="load")
    box = page.locator("input[data-key='map.blackHoleDensity']")
    box.wait_for(state="attached")
    page.locator("input[value=densest]").check()
    box.fill("60")
    box.dispatch_event("change")
    assert page.evaluate(STORED, "map.blackHoleDensity") == "60"
    box.fill("5")
    box.dispatch_event("change")
    assert page.evaluate(STORED, "map.blackHoleDensity") == "60", "an out-of-range number is not stored"
    assert "not between 30 and 1000" in page.locator("[data-preferences-status]").inner_text()
    box.fill("")
    box.dispatch_event("change")
    assert page.evaluate(STORED, "map.blackHoleDensity") is None
    assert page.evaluate(STORED, "map.regionLanding") == "densest", "the other setting stays"
    page.reload(wait_until="load")
    page.wait_for_selector("input[value=densest]", state="attached")
    assert page.locator("input[value=densest]").is_checked() and box.input_value() == ""


def test_the_gear_menu_leads_to_the_preferences_page(page, map_site):
    page.goto(f"{map_site}/", wait_until="load")
    page.locator("header .site-gear sl-icon-button").click()
    link = page.locator("header .site-gear-panel a", has_text="Preferences")
    link.wait_for(state="visible")
    link.click()
    page.wait_for_url("**/preferences")
    assert page.locator("h2", has_text="Preferences").count() == 1
