"""
A sector's own scene fades in as the camera comes up to it (MAP.152), in a real browser on the generated
database `test_web_a11y.py` builds.

In a window of sectors the wheel brings the camera toward the sector at the middle; within 1.5 sector edges of
its centre the sector's scene (its stars and phenomena) is fetched and fades in, whole at 1 edge; it is dropped
past 3 edges. None of that is "opening" the sector: the address and the panel stay the window's, and a click on
the sector then opens it from the scene already there (no second fetch).

Skipped without Playwright, Chromium or MySQL, like the other map browser tests.
"""

import math
import re

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_browser_maps import (  # noqa: E402,F401 -- fixtures
    GL_ARGS,
    VIEWPORT,
    admin_token,
    base_url,
    gl_browser,
    page,
    sample_job,
    sample_job_tree,
    sample_params,
    site_app,
    site_db,
)
from tests.test_web_browser_sector_open import _designation, _ready, _sector  # noqa: E402

CANVAS = "#galaxymap3d-canvas"


def _smoothstep(lo, hi, x):
    t = min(1.0, max(0.0, (x - lo) / (hi - lo)))
    return t * t * (3 - 2 * t)


def _expected_share(distance, edge):
    """sectorreach.js's law: 0 from 1.5 edges out, 1 from 1 edge in, smoothstep between."""
    return 1 - _smoothstep(1.0, 1.5, distance / edge)


def _camera(page):
    return page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera()")


def _ambient(page):
    return page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxyAmbient()")


def _distance(page, centre):
    eye = _camera(page)["eye"]
    return math.dist(eye, centre)


def _wheel(page, spot, delta):
    page.mouse.move(*spot)
    page.mouse.wheel(0, delta)
    page.wait_for_timeout(120)


class _Scene:
    """The sector the test brings the camera up to, and the requests for its scene."""

    def __init__(self, page, base_url, site_app, sample_params):
        self.id = sample_params["sector_id"]
        detail = _sector(site_app, self.id)
        self.centre = (detail["center_x_pc"], detail["center_y_pc"], detail["center_z_pc"])
        scene = page.request.get(f"{base_url}/sector/{self.id}/scene").json()
        self.edge = scene["halfEdgePc"] * 2
        self.code = _designation(detail)
        self.requests = []
        page.on("request", lambda request: self.requests.append(request.url) if re.search(rf"/sector/{self.id}/scene", request.url) else None)

    def fetched(self):
        return len([url for url in self.requests if url.endswith(f"/sector/{self.id}/scene")])


@pytest.fixture()
def reach(page, base_url, site_app, sample_params):
    scene = _Scene(page, base_url, site_app, sample_params)
    page.goto(f"{base_url}/galaxy?sector={scene.code}", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    scene.requests.clear()
    return scene


def _spot(page, reach):
    spot = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(reach.centre))
    assert spot, "the sector's centre is not on the map"
    return (spot["x"], spot["y"])


def _zoom_in_until(page, reach, edges, limit=40):
    """Wheels toward the sector until the camera is within `edges` sector edges of its centre (or no closer is
    possible); returns the distance in edges."""
    for _ in range(limit):
        distance = _distance(page, reach.centre) / reach.edge
        if distance <= edges:
            return distance
        before = _camera(page)["zoom"]
        _wheel(page, _spot(page, reach), -40)
        if _camera(page)["zoom"] >= before - 1e-9:
            break
    return _distance(page, reach.centre) / reach.edge


def test_no_scene_is_loaded_while_the_camera_is_far(page, reach):
    """At the window's fit nothing of the sector's own scene is up and nothing is fetched for it."""
    assert _distance(page, reach.centre) / reach.edge > 1.5
    page.wait_for_timeout(600)
    assert _ambient(page) is None
    assert reach.fetched() == 0


def test_the_scene_fades_in_with_the_distance_and_is_dropped_far_out(page, reach):
    """Coming up to the sector brings its scene in by the law (0 at 1.5 edges, 1 at 1), without opening it; going
    away again drops it past 3 edges."""
    distance = _zoom_in_until(page, reach, 1.4)
    assert distance <= 1.4, f"the window's nearest zoom stops {distance:.2f} edges from the sector's centre"
    # The scene is up once it is fetched; its share is set by the frame after (slow under load).
    page.wait_for_function("() => (document.querySelector('#galaxymap3d-canvas').galaxyAmbient() || {}).share > 0", timeout=30000)
    page.wait_for_timeout(300)
    ambient = _ambient(page)
    assert ambient["id"] == reach.id
    distance = _distance(page, reach.centre) / reach.edge
    assert ambient["share"] == pytest.approx(_expected_share(distance * reach.edge, reach.edge), abs=0.02), (ambient, distance)
    assert "open=1" not in page.url, "the scene was only brought in: the sector is not opened"
    assert reach.fetched() == 1
    # Nearer is more of it, never less, and it is whole by one edge.
    shares = [ambient["share"]]
    for _ in range(40):
        if _distance(page, reach.centre) / reach.edge < 0.9:
            break
        _wheel(page, _spot(page, reach), -40)
        page.wait_for_timeout(100)
        shares.append(_ambient(page)["share"])
    page.wait_for_function("() => (document.querySelector('#galaxymap3d-canvas').galaxyAmbient() || {}).share > 0.999", timeout=30000)
    shares.append(_ambient(page)["share"])
    assert shares == sorted(shares), shares
    # The loop stops at 0.999, the share from which the disc counts as whole; the step that crosses it may land short of 1.
    assert shares[-1] >= 0.999, shares
    assert _camera(page)["zoom"] > 0, "the camera is still in the window"
    # Away again: dropped past 3 edges, and the address never changed.
    for _ in range(80):
        if _distance(page, reach.centre) / reach.edge > 3.2:
            break
        _wheel(page, _spot(page, reach), 40)
    far = _distance(page, reach.centre) / reach.edge
    assert far > 3.2, f"the window's farthest zoom stops {far:.2f} edges from the sector's centre"
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyAmbient() === null", timeout=15000)
    assert "open=1" not in page.url


def test_a_wheel_tick_right_up_to_the_sector_opens_it_where_the_camera_is(page, reach):
    """Zoomed in until the scene is whole, one more wheel tick opens the sector with no flight (the camera does not
    move) and no second fetch; the address becomes the opened sector's."""
    _zoom_in_until(page, reach, 0.6)
    _wheel(page, _spot(page, reach), -600)   # held at the window's nearest zoom, half an edge from the sector
    page.wait_for_function("() => (document.querySelector('#galaxymap3d-canvas').galaxyAmbient() || {}).share > 0.99", timeout=15000)
    before = _camera(page)
    assert before["dist"] / reach.edge < 0.6
    page.wait_for_timeout(500)   # the wheel quiet: a new gesture
    _wheel(page, _spot(page, reach), -40)
    page.wait_for_function("() => location.search.includes('open=1')", timeout=20000)
    _ready(page)
    after = _camera(page)
    assert after["dist"] == pytest.approx(before["dist"], rel=1e-6), (before["dist"], after["dist"])
    assert reach.fetched() == 1
    assert _ambient(page) is None
    # ... and the wheel goes on in from there (the opened sector's own zoom range).
    _wheel(page, _spot(page, reach), -40)
    page.wait_for_timeout(300)
    assert _camera(page)["dist"] < after["dist"]


def _system_requests(page):
    """Collects the requests for any system's scene."""
    seen = []
    page.on("request", lambda request: seen.append(request.url) if re.search(r"/system/[^/]+/scene", request.url) else None)
    return seen


def _star_point(page, base_url, reach):
    """A star of the sector with a system, in the map's frame (pc)."""
    scene = page.request.get(f"{base_url}/sector/{reach.id}/scene").json()
    unit = scene["halfEdgePc"] / (scene["sceneHalfPx"] or 160)
    cx, cy, cz = scene["centerPc"]
    for star in scene["stars"]:
        if star.get("endpoint", "").startswith("system:") and star.get("href"):
            return star, (cx + star["x"] * unit, cy - star["y"] * unit, cz + star["z"] * unit)
    pytest.skip("the sector has no star with a system")


def _system_ambient(page):
    return page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxySystemAmbient()")


def test_a_systems_disc_fades_in_as_the_camera_comes_up_to_its_star(page, base_url, site_app, sample_params):
    """In the opened sector, with the star at the middle of the view, the wheel brings the camera down to where the
    system's disc is whole: its scene is fetched once, fades in by the disc's size, is not "opened" (no system in the
    address), and one more tick opens it from that scene."""
    scene = _Scene(page, base_url, site_app, sample_params)
    page.goto(f"{base_url}/galaxy?sector={scene.code}&open=1", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    star, point = _star_point(page, base_url, scene)
    seen = _system_requests(page)
    page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyCentreOn(...p)", list(point))
    _ready(page)
    assert _system_ambient(page) is None
    position = _still_spot(page, point)
    shares = []
    for _ in range(80):
        _wheel(page, position, -60)
        # Near enough for a disc to show, the camera is given time (a system is asked for after it has stayed in range,
        # and the fetch takes longer under load) before it goes on in: the nearest zoom would open the system.
        near = math.dist(_camera(page)["eye"], point) < 0.05
        page.wait_for_timeout(900 if near else 60)
        ambient = _system_ambient(page)
        if ambient:
            shares.append(ambient["share"])
        if ambient and ambient["share"] >= 0.999:
            break
    assert shares, "the system's scene never came in"
    assert shares == sorted(shares), shares
    assert shares[-1] >= 0.999, shares
    # A small system is first loaded at the 50 AU it is taken for, found too small to stay and dropped; its real size is
    # remembered, so it is loaded once more when its disc is big enough, and not again.
    assert 1 <= len(seen) <= 2, seen
    fetched = len(seen)
    assert "system=" not in page.url
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxySystem()") is None
    # On to the nearest zoom, then one more tick (the wheel quiet first) opens the system from the scene that is up.
    _wheel(page, position, -2000)
    page.wait_for_timeout(600)
    _wheel(page, position, -60)
    page.wait_for_function("() => location.search.includes('system=')", timeout=20000)
    _ready(page)
    assert len(seen) == fetched, "opening the system fetched its scene again"
    assert page.evaluate("() => document.querySelector('#galaxymap3d-canvas').galaxySystem()") is not None
    assert _system_ambient(page) is None


def _anchor_under(page, spot):
    return page.evaluate("(s) => document.querySelector('#galaxymap3d-canvas').galaxyAnchorUnder(...s)", list(spot))


def _still_spot(page, point):
    """Where `point` is on the screen once the camera has stopped moving (the centring is a flight): the same place
    on three reads 100 ms apart. Read at once, a loaded machine gives a spot the star has left by the time the
    anchor is asked for (TEST.231)."""
    last, steady = None, 0
    for _ in range(100):
        spot = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(point))
        here = (spot["x"], spot["y"]) if spot else None
        steady = steady + 1 if here and last and math.dist(here, last) < 0.25 else 0
        if steady >= 3:
            return here
        last = here
        page.wait_for_timeout(100)
    pytest.fail(f"the camera never stopped: {point} was at {last}")


def _shift_nearest_to(page, base_url, reach, point, position):
    """A few pixels' shift from the star's spot (within the 14 px the wheel takes for a star) at which that star is still
    the one nearest to the pointer on the screen, by at least a pixel. The sector's stars are drawn from the database,
    so another star can lie nearer to a fixed shift than the star itself (TEST.231) and the wheel rightly takes it."""
    scene = page.request.get(f"{base_url}/sector/{reach.id}/scene").json()
    unit = scene["halfEdgePc"] / (scene["sceneHalfPx"] or 160)
    cx, cy, cz = scene["centerPc"]
    others = []
    for star in scene["stars"]:
        there = (cx + star["x"] * unit, cy - star["y"] * unit, cz + star["z"] * unit)
        if there != tuple(point):
            spot = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(there))
            if spot:
                others.append((spot["x"], spot["y"]))
    for shift in ((5, -4), (-5, 4), (4, 5), (-4, -5), (6, 0), (0, 6), (-6, 0), (0, -6)):
        click = (position[0] + shift[0], position[1] + shift[1])
        if all(math.dist(click, other) >= math.hypot(*shift) + 1 for other in others):
            return shift
    pytest.skip("another star lies near the target at every shift tried")


def test_the_wheel_anchors_on_the_star_under_the_pointer_not_the_mesh_behind_it(page, base_url, site_app, sample_params):
    """MAP.245: in the opened sector the point the wheel zooms toward is the star under (or a few pixels from) the
    pointer, not the sector's choice mesh behind it; a step in toward it ends nearer the star every time and stops
    short of it, and the eye never goes past the star."""
    scene = _Scene(page, base_url, site_app, sample_params)
    page.goto(f"{base_url}/galaxy?sector={scene.code}&open=1", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    star, point = _star_point(page, base_url, scene)
    page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyCentreOn(...p)", list(point))
    _ready(page)
    position = _still_spot(page, point)
    for shift in ((0, 0), _shift_nearest_to(page, base_url, scene, point, position)):
        anchor = _anchor_under(page, (position[0] + shift[0], position[1] + shift[1]))
        assert anchor == pytest.approx(list(point), abs=1e-9), (shift, anchor, point)
    distances = [math.dist(_camera(page)["eye"], point)]
    for _ in range(12):
        _wheel(page, position, -60)
        distances.append(math.dist(_camera(page)["eye"], point))
        if "system=" in page.url:
            break
    assert all(later < earlier for earlier, later in zip(distances, distances[1:])), distances
    assert distances[-1] > 0
    # Each tick keeps the star where it is on the screen, at the pointer.
    after = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(point))
    assert (after["x"], after["y"]) == pytest.approx(position, abs=2.0)


def test_a_double_click_on_empty_space_flies_half_way_toward_it(page, base_url):
    """MAP.245: on empty space (nothing picked) the flight brings the eye half way to the point under the pointer (the
    focus plane's, here), which stays where it is on the screen."""
    page.goto(f"{base_url}/galaxy", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    rect = page.evaluate("() => { const r = document.querySelector('#galaxymap3d-canvas').getBoundingClientRect(); return [r.left, r.top, r.width, r.height]; }")
    empty = (rect[0] + 6, rect[1] + 6)
    anchor = _anchor_under(page, empty)
    assert anchor, "no point under the corner of the map"
    before = _camera(page)
    page.mouse.dblclick(*empty)
    after = _landed(page, anchor)
    assert after["dist"] == pytest.approx(before["dist"] * 0.5, rel=0.05), (before["dist"], after["dist"])
    spot_after = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(anchor))
    assert (spot_after["x"], spot_after["y"]) == pytest.approx(empty, abs=2.0)


def test_the_wheel_anchor_on_empty_space_is_where_the_column_of_stars_passes_the_threshold(page, base_url, site_app, sample_params):
    """MAP.247: in the opened sector a cursor ray that meets no star has its point where the column of the galaxy's
    stars along it passes half the view's distance in local-average density, not on the focus plane: the fixture
    sector lies near the galaxy's core, where the density is far above the average, so the point is the nearest one
    allowed (a quarter of the way to the focus plane)."""
    scene = _Scene(page, base_url, site_app, sample_params)
    page.goto(f"{base_url}/galaxy?sector={scene.code}&open=1", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    rect = page.evaluate("() => { const r = document.querySelector('#galaxymap3d-canvas').getBoundingClientRect(); return [r.left, r.top, r.width, r.height]; }")
    camera = _camera(page)
    for corner in ((rect[0] + 6, rect[1] + 6), (rect[0] + rect[2] - 6, rect[1] + rect[3] - 6)):
        anchor = _anchor_under(page, corner)
        assert anchor, f"no point under {corner}"
        away = math.dist(camera["eye"], anchor)
        # The focus plane is at least `dist` from the eye along any ray of the view (at most about 1.2 times it).
        assert 0.2 * camera["dist"] <= away <= 0.35 * camera["dist"], (away, camera["dist"])


def _fly_setup(page, base_url, site_app, sample_params):
    """The sector opened, and a star of it with a system: the star's point in the map (pc) and where it is on screen."""
    scene = _Scene(page, base_url, site_app, sample_params)
    page.goto(f"{base_url}/galaxy?sector={scene.code}&open=1", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    star, point = _star_point(page, base_url, scene)
    _still(page)
    spot = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(point))
    assert spot, "the star is not on the map"
    return point, (spot["x"], spot["y"])


def _still(page):
    """Until no flight is under way and the camera is where it was a moment ago: `_ready` goes by the scale label, which
    a camera that is still settling on a loaded machine can leave unchanged for a quarter second, and a click aimed at
    where the star was then misses it (TEST.262)."""
    last = None
    for _ in range(60):
        camera = _camera(page)
        now = (camera["flying"], camera["target"], camera["dist"])
        if not camera["flying"] and now == last:
            return
        last = now
        page.wait_for_timeout(250)


def _landed(page, point):
    """Waits for the flight to end and reads where the camera is."""
    page.wait_for_function("() => document.querySelector('#galaxymap3d-canvas').galaxyCamera().flying === false", timeout=30000)
    page.wait_for_timeout(200)
    return _camera(page)


def test_a_double_click_on_a_star_flies_the_camera_to_it(page, base_url, site_app, sample_params):
    """MAP.146: two clicks on a star take the camera to look at it from where its system's disc is whole, without
    opening the system (the address is the sector's still)."""
    point, spot = _fly_setup(page, base_url, site_app, sample_params)
    before = _camera(page)
    page.mouse.dblclick(*spot)
    camera = _landed(page, point)
    assert math.dist(camera["target"], point) < 1e-6 * max(1.0, before["dist"]), (camera["target"], point)
    assert camera["dist"] < before["dist"] * 0.5, (before["dist"], camera["dist"])
    assert "system=" not in page.url


def test_a_double_click_flies_even_when_the_page_is_slow_to_answer_the_first_click(page, base_url, site_app, sample_params):
    """TEST.252: the two clicks are told apart by the time of their events, not by when the page gets to them, so a
    first click that takes the page over the double click window (a busy machine) is still the first half of a pair."""
    point, spot = _fly_setup(page, base_url, site_app, sample_params)
    before = _camera(page)
    page.evaluate("""() => {
        const canvas = document.querySelector('#galaxymap3d-canvas');
        window.__events = [];
        for (const type of ['pointerdown', 'pointerup', 'dblclick']) {
            canvas.addEventListener(type, (event) => window.__events.push([type, event.timeStamp, performance.now()]), true);
        }
        let first = true;
        canvas.addEventListener('pointerup', () => {
            if (!first) return;
            first = false;
            const start = performance.now();
            while (performance.now() - start < 500);
        });
    }""")
    page.mouse.dblclick(*spot)
    camera = _landed(page, point)
    # The events' times are in the message so a failure under load shows whether the pair was told apart (TEST.262).
    events = page.evaluate("() => window.__events")
    assert math.dist(camera["target"], point) < 1e-6 * max(1.0, before["dist"]), (
        camera["target"], point, spot, before, events)


def test_the_go_to_button_of_a_selected_star_flies_there_too(page, base_url, site_app, sample_params):
    """The info panel of a selected star has a Go to button."""
    point, spot = _fly_setup(page, base_url, site_app, sample_params)
    before = _camera(page)
    page.mouse.click(*spot)
    button = page.get_by_role("button", name="Go to")
    button.wait_for(state="visible", timeout=15000)
    button.click()
    camera = _landed(page, point)
    assert math.dist(camera["target"], point) < 1e-6 * max(1.0, before["dist"]), (camera["target"], point)
    assert camera["dist"] < before["dist"] * 0.5


def test_enter_with_a_star_selected_flies_to_it(page, base_url, site_app, sample_params):
    point, spot = _fly_setup(page, base_url, site_app, sample_params)
    before = _camera(page)
    page.mouse.click(*spot)
    page.get_by_role("button", name="Go to").wait_for(state="visible", timeout=15000)
    page.locator(CANVAS).focus()
    page.keyboard.press("Enter")
    camera = _landed(page, point)
    assert math.dist(camera["target"], point) < 1e-6 * max(1.0, before["dist"]), (camera["target"], point)
    assert camera["dist"] < before["dist"] * 0.5
