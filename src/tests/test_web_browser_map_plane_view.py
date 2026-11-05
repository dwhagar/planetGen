"""
MAP.371: the Galaxy Map's plane view. The Menu's Plane view button tilts the camera to a preset and lays a grid on a plane through the focus:
rings at round distances inside a 60 degree wedge of the galaxy, from the galactic center, with a label on each ring; "Plane at galactic
midplane" puts the grid at height 0 instead; a selection moves the rings to it; turning the view off takes the grid away and the stage's own tilt
back. Skipped without Playwright, Chromium or WebGL, like the other map browser tests.
"""

import re

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_browser_fixture_maps import (  # noqa: E402,F401 -- fixtures
    GALAXY_CANVAS, _click_galaxy_stars, _open_galaxy, _open_star_field, _settle, gl_browser, map_site, page)

CAMERA = "() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()"
PLANE = "() => document.querySelector('#galaxymap3d-canvas').galaxyPlane()"
TILT_DEG = 62
LABEL = re.compile(r"^\d+(\.\d+)? (pc|kpc)$")


def _calm(page):
    page.wait_for_function(f"() => {{ const c = ({CAMERA})(); return c && !c.flying; }}")
    page.wait_for_timeout(200)


def _menu_click(page, action):
    if not page.evaluate("() => document.querySelector('#galaxymap3d-menu').open"):
        page.click("#galaxymap3d-menu > summary")
    page.click(f'[data-action="{action}"]')
    _calm(page)


def test_the_plane_view_tilts_the_camera_and_lays_a_ring_grid_inside_a_sixty_degree_wedge(page, map_site):
    _open_galaxy(page, map_site)
    assert page.evaluate(PLANE)["on"] is False and page.evaluate(CAMERA)["tilt"] < 1, "the galaxy opens straight down, no grid"
    _menu_click(page, "plane-view")
    assert page.get_attribute('[data-action="plane-view"]', "aria-pressed") == "true"
    assert abs(page.evaluate(CAMERA)["tilt"] - TILT_DEG) < 1
    page.wait_for_function(f"() => ({PLANE})().rings > 0")
    plane = page.evaluate(PLANE)
    assert plane["on"] and plane["visible"] and plane["wedgeDeg"] == 60 and plane["segments"] > 0
    assert plane["around"] == "galactic center"
    assert plane["maxOffsetDeg"] <= 30.01, "every ring lies inside the 60 degree wedge"
    assert plane["labels"] and all(LABEL.match(label["text"]) for label in plane["labels"]), plane["labels"]
    _menu_click(page, "plane-view")
    assert abs(page.evaluate(CAMERA)["tilt"]) < 1, "off, the camera takes the stage's own tilt again"
    off = page.evaluate(PLANE)
    assert off["on"] is False and off["visible"] is False and off["labels"] == []


def test_the_grid_follows_the_focus_and_the_toggle_between_the_two_planes_works(page, map_site):
    _open_star_field(page, map_site)
    _menu_click(page, "plane-view")
    page.wait_for_function(f"() => ({PLANE})().rings > 0")
    first = page.evaluate(PLANE)
    focus = page.evaluate(CAMERA)["target"]
    assert first["anchor"][:2] == pytest.approx(focus[:2], abs=1e-6 * max(1, abs(focus[0])) + 1e-3 * first["step"])
    assert first["planeZ"] == pytest.approx(focus[2], abs=1e-3 * first["step"] + 1e-6), "the plane passes through the focus"
    # Move the focus on across the view with Alt and an arrow: the grid is rebuilt about the new focus.
    page.locator(GALAXY_CANVAS).focus()
    for _ in range(12):
        page.keyboard.press("Alt+ArrowRight")
    page.wait_for_function(f"() => {{ const a = ({PLANE})().anchor; return Math.hypot(a[0] - {first['anchor'][0]}, a[1] - {first['anchor'][1]}) > 0; }}")
    # The grid is rebuilt about the focus once it has moved a twentieth of the view (a ring step is a quarter of the view).
    page.wait_for_function(f"""() => {{ const p = ({PLANE})(), c = ({CAMERA})();
        return !c.flying && Math.hypot(p.anchor[0] - c.target[0], p.anchor[1] - c.target[1]) <= 0.25 * p.step + 1e-6; }}""")
    moved = page.evaluate(PLANE)
    assert moved["anchor"][:2] != first["anchor"][:2], "the grid moved with the focus"
    assert moved["maxOffsetDeg"] <= 30.01
    # Up the screen of the tilted view is up out of the midplane too: the focus leaves height 0, so the two planes differ.
    page.locator(GALAXY_CANVAS).focus()
    for _ in range(12):
        page.keyboard.press("Alt+ArrowUp")
    page.wait_for_function(f"""() => {{ const p = ({PLANE})(), c = ({CAMERA})();
        return !c.flying && Math.abs(c.target[2]) > 1e-3 && Math.abs(p.planeZ - c.target[2]) <= 0.25 * p.step + 1e-6; }}""")
    assert abs(page.evaluate(CAMERA)["target"][2]) > 1e-3 and page.evaluate(PLANE)["planeZ"] != 0
    _menu_click(page, "plane-midplane")
    assert page.get_attribute('[data-action="plane-midplane"]', "aria-pressed") == "true"
    page.wait_for_function(f"() => ({PLANE})().midplane === true && ({PLANE})().planeZ === 0")
    assert page.evaluate(PLANE)["anchor"][2] == 0
    _menu_click(page, "plane-midplane")
    page.wait_for_function(f"() => ({PLANE})().midplane === false")
    assert page.evaluate(PLANE)["planeZ"] == pytest.approx(page.evaluate(CAMERA)["target"][2], abs=1e-3 * page.evaluate(PLANE)["step"] + 1e-6)


def test_the_midplane_button_is_off_until_the_plane_view_is_on(page, map_site):
    _open_galaxy(page, map_site)
    page.click("#galaxymap3d-menu > summary")
    assert page.is_disabled('[data-action="plane-midplane"]')
    page.click('[data-action="plane-view"]')
    _calm(page)
    assert not page.is_disabled('[data-action="plane-midplane"]')


def test_a_selection_draws_the_rings_around_it_inside_the_wedge(page, map_site):
    _open_star_field(page, map_site)
    _menu_click(page, "plane-view")
    page.wait_for_function(f"() => ({PLANE})().rings > 0")
    assert page.evaluate(PLANE)["around"] == "galactic center"
    page.click("#galaxymap3d-menu > summary")    # close the menu: it sits over the map
    _click_galaxy_stars(page, 1)
    page.wait_for_function(f"() => ({PLANE})().around === 'the selection'")
    selected = page.evaluate(PLANE)
    assert selected["rings"] > 0 and selected["maxOffsetDeg"] <= 30.01
    assert selected["labels"] and all(LABEL.match(label["text"]) for label in selected["labels"])


def test_ring_distances_are_round_and_named_in_parsecs_and_kiloparsecs(page, map_site):
    _open_galaxy(page, map_site)
    script = """async () => {
        const src = [...document.querySelectorAll('script[type=module]')].map(s => s.src).find(s => s.includes('galaxymap3d.js'));
        const mod = await import(src.replace('galaxymap3d.js', 'galaxyplane.js'));
        return [[0.9, 3, 7, 12, 480, 2600, 9000].map(mod.roundStep), [500, 1000, 2500, 20, 0.5].map(mod.distanceLabel)];
    }"""
    steps, labels = page.evaluate(script)
    assert steps == [1, 2, 5, 10, 500, 2000, 10000]
    assert labels == ["500 pc", "1 kpc", "2.5 kpc", "20 pc", "0.5 pc"]


def _stem_ends(plane):
    """The stems as ((x, y, z), (x, y, z)) pairs in galaxy coordinates, from the local ones the plane keeps about its anchor."""
    flat, anchor = plane["stemPositions"], plane["anchor"]
    points = [tuple(flat[i + k] + anchor[k] for k in range(3)) for i in range(0, len(flat), 3)]
    return list(zip(points[0::2], points[1::2]))


def test_stems_drop_from_chosen_stars_to_the_plane_inside_the_wedge_within_the_budget(page, map_site):
    """MAP.372: one stem per chosen object, straight to its shadow on the grid plane, inside the wedge, never one under 2 px, never more
    than the budget; with the plane at the galaxy's midplane the stars off it have stems. The cost of choosing and building them is measured."""
    import math

    _open_star_field(page, map_site)
    _menu_click(page, "plane-view")
    _menu_click(page, "plane-midplane")
    page.wait_for_function(f"() => ({PLANE})().stems > 0")
    plane = page.evaluate(PLANE)
    assert plane["planeZ"] == 0 and 0 < plane["stems"] <= plane["stemBudget"]
    assert plane["stemShortestPx"] >= 2, "a stem shorter than 2 px is skipped"
    ends = _stem_ends(plane)
    assert len(ends) == plane["stems"]
    centre = math.radians(plane["centreDeg"])
    for top, foot in ends:
        assert top[:2] == pytest.approx(foot[:2], abs=1e-3), "a stem is straight up from its shadow"
        assert foot[2] == pytest.approx(plane["planeZ"], abs=1e-3) and top[2] != pytest.approx(foot[2], abs=1e-6)
        offset = abs((math.atan2(top[1], top[0]) - centre + math.pi) % (2 * math.pi) - math.pi)
        assert math.degrees(offset) <= 30.01, "stems live inside the 60 degree wedge"
    kinds = plane["stemKinds"]
    assert kinds["bright"] > 0 and sum(kinds.values()) == plane["stems"]
    # The frame cost at the budget: choosing and building the stems takes well under a frame (a running mean of the last runs, ms).
    page.wait_for_function(f"() => ({PLANE})().stemMs > 0")
    cost = page.evaluate(PLANE)["stemMs"]
    print(f"MAP.372 stem rebuild: {cost:.2f} ms mean at {plane['stems']} stems")
    assert cost < 8, cost
    _menu_click(page, "plane-view")
    assert page.evaluate(PLANE)["stems"] == 0 and page.evaluate(PLANE)["stemPositions"] == []


def test_the_stem_of_a_selected_star_is_the_strongest_and_the_plane_through_it_gives_it_none(page, map_site):
    _open_star_field(page, map_site)
    _menu_click(page, "plane-view")
    _menu_click(page, "plane-midplane")
    page.click("#galaxymap3d-menu > summary")    # close the menu: it sits over the map
    _click_galaxy_stars(page, 1)
    page.wait_for_function(f"() => ({PLANE})().around === 'the selection'")
    page.wait_for_function(f"() => ({PLANE})().stems > 0")
    plane = page.evaluate(PLANE)
    selected_z = plane["stemKinds"]["selected"]
    assert selected_z in (0, 1), "the selection's stem is drawn first when it stands off the midplane (one in the fixture's field is on it)"
    assert plane["stems"] <= plane["stemBudget"]
