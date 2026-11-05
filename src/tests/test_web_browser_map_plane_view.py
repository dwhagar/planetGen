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


TRIPLE = re.compile(r"^(?P<r>[\d.]+) (?P<unit>pc|kpc), (?P<lon>[\d.]+)°, (?P<lat>[+-][\d.]+)°(?P<rest> from the selection)?$")


def _parse_triple(text):
    found = TRIPLE.match(text)
    assert found, text
    scale = 1000 if found["unit"] == "kpc" else 1
    return float(found["r"]) * scale, float(found["lon"]), float(found["lat"]), bool(found["rest"])


def _spherical(point, origin):
    import math

    dx, dy, dz = (point[i] - origin[i] for i in range(3))
    r = math.sqrt(dx * dx + dy * dy + dz * dz)
    return r, math.degrees(math.atan2(dy, dx)) % 360, math.degrees(math.asin(dz / r)) if r else 0.0


def _set_stretch(page, value):
    page.evaluate("""(value) => { const input = document.querySelector('#galaxymap3d-plane-stretch');
        input.value = String(value); input.dispatchEvent(new Event('input', {bubbles: true})); }""", value)


def test_the_stems_switch_and_the_height_stretch_work_and_the_grid_says_how_tall_heights_are_drawn(page, map_site):
    """MAP.373: Stems turns the stems off and on; the stretch (1 to 10) draws each stem on past its star to a bead at that many times the
    height, and the grid's note says so."""
    _open_star_field(page, map_site)
    page.click("#galaxymap3d-menu > summary")
    assert page.is_disabled('[data-action="plane-stems"]') and page.is_disabled("#galaxymap3d-plane-stretch")
    page.click('[data-action="plane-view"]')
    _calm(page)
    assert not page.is_disabled('[data-action="plane-stems"]') and not page.is_disabled("#galaxymap3d-plane-stretch")
    _menu_click(page, "plane-midplane")
    page.wait_for_function(f"() => ({PLANE})().stems > 0")
    plane = page.evaluate(PLANE)
    assert plane["stemsOn"] and plane["stretch"] == 1 and plane["beads"] == 0 and plane["note"] == "Heights at true scale"
    _set_stretch(page, 5)
    page.wait_for_function(f"() => ({PLANE})().stretch === 5 && ({PLANE})().beads > 0")
    stretched = page.evaluate(PLANE)
    assert stretched["beads"] == stretched["stems"] and "5 times taller than true" in stretched["note"]
    assert page.inner_text("#galaxymap3d-plane-stretch-out") == "5 times"
    page.wait_for_function("() => document.querySelector('.map-plane-note').textContent.includes('5 times')")
    _set_stretch(page, 99)
    page.wait_for_function(f"() => ({PLANE})().stretch === 10")
    _set_stretch(page, 1)
    page.wait_for_function(f"() => ({PLANE})().stretch === 1 && ({PLANE})().beads === 0")
    _menu_click(page, "plane-stems")
    assert page.get_attribute('[data-action="plane-stems"]', "aria-pressed") == "false"
    page.wait_for_function(f"() => ({PLANE})().stems === 0")
    assert page.evaluate(PLANE)["stemPositions"] == []
    _menu_click(page, "plane-stems")
    page.wait_for_function(f"() => ({PLANE})().stems > 0")


def test_the_label_of_a_stem_is_the_stars_spherical_triple_from_the_galactic_center_and_matches_its_stored_position(page, map_site):
    from tests.map_site_support import POINT_FIELD, STAR_FIELD, STAR_FIELD_CENTER_PC

    _open_star_field(page, map_site)
    _menu_click(page, "plane-view")
    _menu_click(page, "plane-midplane")
    page.wait_for_function(f"() => ({PLANE})().coords.length > 0")
    plane = page.evaluate(PLANE)
    assert plane["context"] == "galaxy"
    matched = 0
    for label in plane["coords"]:
        r, lon, lat, relative = _parse_triple(label["text"])
        assert not relative
        want = _spherical(label["star"], (0, 0, 0))
        assert r == pytest.approx(want[0], rel=6e-3, abs=0.06) and abs(lon - want[1]) <= 0.06 and abs(lat - want[2]) <= 0.06, (label, want)
        if any(max(abs(label["star"][i] - star[k]) for i, k in enumerate("xyz")) < 0.01 for star in [*STAR_FIELD, *POINT_FIELD, {"x": STAR_FIELD_CENTER_PC[0] + 6.0, "y": -6.0, "z": 2.0}]):    # the last is Fixture Eden, a named world
            matched += 1
    assert matched, ("a label is for an object of the fixture, at its stored position", [label["star"] for label in plane["coords"]])


def test_with_a_star_selected_the_other_labels_are_distance_bearing_and_elevation_from_it(page, map_site):
    import math

    _open_star_field(page, map_site)
    _menu_click(page, "plane-view")
    _menu_click(page, "plane-midplane")
    page.click("#galaxymap3d-menu > summary")    # close the menu: it sits over the map
    _click_galaxy_stars(page, 1)
    page.wait_for_function(f"() => ({PLANE})().around === 'the selection' && ({PLANE})().coords.some((c) => c.selected)")
    plane = page.evaluate(PLANE)
    chosen = next(label for label in plane["coords"] if label["selected"])
    assert not _parse_triple(chosen["text"])[3], "the selected star's own label is its position in the context"
    others = [label for label in plane["coords"] if not label["selected"]]
    for label in others:
        r, bearing, elevation, relative = _parse_triple(label["text"])
        assert relative
        want = _spherical(label["star"], chosen["star"])
        assert r == pytest.approx(want[0], rel=6e-3, abs=0.06) and abs(bearing - want[1]) <= 0.06 and abs(elevation - want[2]) <= 0.06
        assert r == pytest.approx(math.dist(label["star"], chosen["star"]), rel=6e-3, abs=0.06)
    assert plane["coords"], "the selection is labelled"


def test_in_an_opened_sector_the_triple_is_measured_from_the_sector_center(page, map_site):
    from tests.test_web_browser_fixture_maps import OPEN_PRIME

    _open_galaxy(page, map_site, OPEN_PRIME)
    _menu_click(page, "plane-view")
    page.wait_for_function(f"() => ({PLANE})().context === 'sector'")
    assert page.evaluate(PLANE)["context"] == "sector"
    for label in page.evaluate(PLANE)["coords"]:
        assert TRIPLE.match(label["text"]), label
