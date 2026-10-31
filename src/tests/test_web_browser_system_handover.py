"""
The wheel carries the Galaxy Map into a star system and back out as one move (MAP.326): in the opened sector, scrolling
in on a star brings its system's scene up and then opens it where the camera is (no flight, no page load), and scrolling
out closes it again, for a single star, a pair and a planet with moons.

Skipped without Playwright, Chromium or MySQL, like the other map browser tests.
"""

import math

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
    sample_run,
    sample_params,
    site_app,
    site_db,
)
from tests.test_web_browser_sector_open import _ready, _scene  # noqa: E402
from tests.test_web_browser_sector_reach import CANVAS, _Scene, _camera, _wheel  # noqa: E402

CANVAS_JS = "document.querySelector('#galaxymap3d-canvas')"


def _kinds(client):
    """A generated system of each kind: {kind: system id}."""
    found = {}
    for system in client.get("/api/systems?limit=100").get_json()["items"]:
        scene = client.get(f"/api/systems/{system['id']}/scene").get_json()
        moons = sum(len(p.get("moons", [])) for p in scene["planets"])
        kind = "moons" if moons >= 3 and len(scene["stars"]) == 1 else "binary" if len(scene["stars"]) == 2 and scene["planets"] else "single" if len(scene["stars"]) == 1 and scene["planets"] else None
        if kind and kind not in found:
            found[kind] = system["id"]
    return found


def _star_point_of(reach, system_id):
    """The map-frame position (pc) of the sector star that opens `system_id`."""
    scene = _scene(reach.app, reach.id)
    unit = scene["halfEdgePc"] / (scene["sceneHalfPx"] or 160)
    cx, cy, cz = scene["centerPc"]
    for star in scene["stars"]:
        if star.get("endpoint") == f"system:{system_id}" and star.get("href"):
            return (cx + star["x"] * unit, cy - star["y"] * unit, cz + star["z"] * unit)
    pytest.skip(f"the sector has no star for system {system_id}")


def _state(page):
    return page.evaluate(f"() => ({{system: {CANVAS_JS}.galaxySystem() ? {CANVAS_JS}.galaxySystem().id : null, "
                         f"ambient: ({CANVAS_JS}.galaxySystemAmbient() || {{}}).share || 0, url: location.search}})")


@pytest.mark.parametrize("kind", ["single", "binary", "moons"])
def test_the_wheel_carries_the_map_into_a_system_and_back_out_with_no_page_load(page, base_url, site_app, sample_params, kind):
    systems = _kinds(site_app.test_client())
    if kind not in systems:
        pytest.skip(f"the generated sector has no {kind} system")
    reach = _Scene(page, base_url, site_app, sample_params)
    point = _star_point_of(reach, systems[kind])
    page.goto(f"{base_url}/galaxy?sector={reach.code}&open=1", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    page.evaluate("() => { window.__same_page = true; }")
    page.evaluate(f"(p) => {CANVAS_JS}.galaxyCentreOn(...p)", list(point))
    _ready(page)
    page.wait_for_timeout(1500)
    states = []
    for _ in range(220):
        spot = page.evaluate(f"(p) => {CANVAS_JS}.galaxyScreenOf(...p)", list(point))
        assert spot, "the star left the map"
        _wheel(page, (spot["x"], spot["y"]), -60)
        states.append(_state(page))
        if states[-1]["system"]:
            break
    assert states[-1]["system"], f"never entered the system: {states[-1]}"
    assert any(s["ambient"] > 0 for s in states), "the system's scene never faded in before it opened"
    assert page.evaluate("() => window.__same_page === true"), "the page was loaded again"
    assert "system=" in page.url
    inside = page.evaluate(f"() => {CANVAS_JS}.galaxySystem()")
    assert inside["refs"], "the system has bodies on the map"
    # And out again: the wheel the other way leaves the system for the sector, still the same page.
    for _ in range(260):
        _wheel(page, (page.viewport_size["width"] / 2, page.viewport_size["height"] / 2), 60)
        if not _state(page)["system"]:
            break
    assert not _state(page)["system"], "the wheel out never left the system"
    assert page.evaluate("() => window.__same_page === true")


def test_the_wheel_carries_the_map_from_the_whole_galaxy_into_a_system(page, base_url, site_app, sample_params):
    """MAP.326: from the whole-galaxy view, aiming the wheel at one star all the way in, no pause and no click, the
    stages hand over to each other and the system opens (the camera used to stop 0.5 pc short, the nearest zoom being
    judged from the view's middle and not from the star under the pointer); the wheel out goes back the same way."""
    systems = _kinds(site_app.test_client())
    reach = _Scene(page, base_url, site_app, sample_params)
    point = _star_point_of(reach, systems.get("moons") or next(iter(systems.values())))
    page.goto(f"{base_url}/galaxy", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    page.evaluate("() => { window.__same_page = true; }")
    trail = []
    for n in range(600):
        spot = page.evaluate(f"(p) => {CANVAS_JS}.galaxyScreenOf(...p)", list(point))
        if not spot:
            trail.append(("offscreen", n))
            break
        _wheel(page, (spot["x"], spot["y"]), -60)
        state = _state(page)
        cam = _camera(page)
        trail.append((n, round(math.dist(cam["eye"], point), 6), state["system"] is not None, state["url"][:60]))
        if state["system"]:
            break
    assert _state(page)["system"], trail[-3:]
    assert page.evaluate("() => window.__same_page === true"), "the page was loaded again"
    # The wheel back out carries the map out of the system, the sector and on up the stages with no page load.
    for _ in range(400):
        _wheel(page, (page.viewport_size["width"] / 2, page.viewport_size["height"] / 2), 120)
        if "sector=" not in page.url:
            break
    assert "sector=" not in page.url, page.url
    assert page.evaluate("() => window.__same_page === true")


def test_the_stars_and_bodies_of_a_system_opened_in_the_sector_turn_on_their_axes_and_the_menu_switch_holds_them_still(page, base_url, site_app, sample_params):
    """MAP.328: a system opened in the sector of the Galaxy Map draws its stars and bodies as spheres, and they turn about their
    stored axes by spin.js's law (the same spin time for every body, none faster than one turn in 20 s); a body with no stored
    period stays where it is; the Menu's Spin switch holds them still; the galaxy view itself turns nothing."""
    systems = _kinds(site_app.test_client())
    system_id = systems.get("moons") or systems.get("single") or next(iter(systems.values()))
    reach = _Scene(page, base_url, site_app, sample_params)
    page.emulate_media(reduced_motion="no-preference")
    page.goto(f"{base_url}/galaxy?sector={reach.code}&open=1&system={system_id}", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    page.wait_for_function(f"() => {{ const s = {CANVAS_JS}.galaxySystem(); return s && s.spins().length > 0; }}", timeout=30000)
    sample = f"""async (wait) => {{
        const read = () => {CANVAS_JS}.galaxySystem().spins();
        const first = read();
        await new Promise((resolve) => setTimeout(resolve, wait));
        return {{first: first, second: read()}};
    }}"""
    got = page.evaluate(sample, 900)
    turning = [a for a in got["first"] if a["rate"] > 0]
    assert len(turning) >= 2, "bodies with a stored period are drawn"
    elapsed = {a["ref"]: b["at"] - a["at"] for a, b in zip(got["first"], got["second"])}
    assert all(0.3 < e < 3.0 for e in elapsed.values()), elapsed
    for a, b in zip(got["first"], got["second"]):
        dot = abs(sum(x * y for x, y in zip(a["quaternion"], b["quaternion"])))
        turned = 2 * math.acos(min(1.0, dot))
        if a["locked"]:
            # UX.244: a tidally locked body keeps its face toward what it orbits, so it moves only as its orbit does.
            assert turned < 0.02, (a["ref"], turned)
        else:
            assert turned == pytest.approx(a["rate"] * (b["at"] - a["at"]), abs=2e-3), (a["ref"], turned, a["rate"])
        assert 2 * math.pi / 600 - 1e-9 <= a["rate"] <= 2 * math.pi / 20 + 1e-9 or a["rate"] == 0
    # The Menu's Spin switch holds every body still, and on again sets them turning.
    page.locator("#galaxymap3d-menu > summary").click()
    button = page.locator('#galaxymap3d-kinds button[data-layer="spin"]')
    assert button.count() == 1 and button.get_attribute("aria-pressed") == "true"
    button.click()
    held = page.evaluate(sample, 600)
    assert held["first"][0]["at"] == held["second"][0]["at"], "spin time stopped"
    for a, b in zip(held["first"], held["second"]):
        if a["locked"]:
            assert a["quaternion"] == pytest.approx(b["quaternion"], abs=1e-3)   # UX.244: it follows its orbit, not the spin
        else:
            assert a["quaternion"] == b["quaternion"]
    button.click()
    again = page.evaluate(sample, 600)
    assert again["second"][0]["at"] > again["first"][0]["at"]


def test_a_flight_from_the_whole_galaxy_to_a_planet_switches_no_stage_and_every_scene_follows_the_camera(page, base_url, site_app, sample_params):
    """MAP.335: the camera flown from the whole galaxy to a star lands at the whole-galaxy stage it left; the sector's scene and
    then the system's come up because the eye is near them (the server names the sectors near the eye, `/galaxy/near`, and no data
    window is asked for), the wheel opens the sector and the system where the camera is, and a flight to a planet ends at it."""
    systems = _kinds(site_app.test_client())
    system_id = systems.get("moons") or next(iter(systems.values()))
    reach = _Scene(page, base_url, site_app, sample_params)
    point = _star_point_of(reach, system_id)
    requests = []
    page.on("request", lambda request: requests.append(request.url))
    page.goto(f"{base_url}/galaxy", wait_until="load")
    page.wait_for_selector(CANVAS, state="attached")
    _ready(page)
    page.evaluate("() => { window.__same_page = true; }")
    stage = page.evaluate(f"() => {CANVAS_JS}.galaxyStage()")
    assert stage["win"] is None, stage
    asked = len(requests)
    page.evaluate(f"(p) => {CANVAS_JS}.galaxyGoTo({{center: p, dist: 0.0001}})", list(point))
    page.wait_for_function(f"() => !{CANVAS_JS}.galaxyMoving()", timeout=60000)
    # The sector's scene, then the system's, come up from where the eye is.
    page.wait_for_function(f"() => ({CANVAS_JS}.galaxyAmbient() || {{}}).share > 0.99", timeout=60000)
    page.wait_for_function(f"() => ({CANVAS_JS}.galaxySystemAmbient() || {{}}).share > 0", timeout=60000)
    assert page.evaluate(f"() => {CANVAS_JS}.galaxyStage()") == stage, "the flight switched the stage"
    # In to the star: the sector, then its system, open where the camera is.
    for _ in range(80):
        spot = page.evaluate(f"(p) => {CANVAS_JS}.galaxyScreenOf(...p)", list(point))
        assert spot, "the star left the map"
        _wheel(page, (spot["x"], spot["y"]), -60)
        if _state(page)["system"]:
            break
    assert _state(page)["system"], _state(page)
    assert page.evaluate(f"() => {CANVAS_JS}.galaxyStage()") == stage, "opening the sector and the system switched the stage"
    refs = page.evaluate(f"() => {CANVAS_JS}.galaxySystem().refs")
    planets = [ref for ref in refs if ref.startswith("planet:")]
    assert planets, refs
    assert page.evaluate(f"(r) => {CANVAS_JS}.galaxyGoTo({{body: r}})", planets[0])
    page.wait_for_function(f"() => !{CANVAS_JS}.galaxyMoving()", timeout=60000)
    assert page.evaluate(f"() => {CANVAS_JS}.galaxyStage()") == stage
    assert page.evaluate("() => window.__same_page === true"), "the page was loaded again"
    later = requests[asked:]
    assert any("/galaxy/near" in url for url in later), "no sector was asked for by the eye's place"
    assert not any("/galaxy/window" in url or "/galaxy/stage" in url for url in later), [u for u in later if "/galaxy/" in u]
    assert _camera(page)["dist"] < 1e-3
