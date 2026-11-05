"""
VIEW.33: the sky of a system (`static/skysphere.js`) in a real browser. A packed sky is drawn as one points object round the camera and
each star lands where the light-time arithmetic puts it: a star at rest in the middle of the view, a star moving across the line of sight
shifted by about v / c (the vertex shader agrees with the JavaScript twin `retardedOffset` to a pixel and a half), a nearby star that drifts as the orbit clock
runs and one that does not, the observer's own velocity, a star past the limiting magnitude that is not drawn, and a phone that draws the brightest.

The page is empty and the static files are served to it from disk, so no database or site is needed. Skipped without Playwright or Chromium.
"""

import base64
import mimetypes
import os
from types import SimpleNamespace

import numpy as np
import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from planetgen.web.lib.skypack import KIND_CODES, pack_sky  # noqa: E402

_STATIC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "html", "static")
ORIGIN = "http://sky.test"
SIZE = 400
EPOCH = 1_800_000_000.0
LIGHT_KMS = 299792.458

# [name, offset pc, absolute velocity km/s, absolute magnitude, temperature K]
STARS = [
    ("rest", (10.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0, 6000.0),
    ("fast", (1.0, 0.3, 0.2), (0.0, 29000.0, 0.0), 5.0, 9000.0),
    ("drift", (1.5, -0.5, 0.3), (0.0, -100.0, 0.0), 3.5, 4500.0),
    ("dim", (500.0, 100.0, 0.0), (0.0, 0.0, 0.0), 5.0, 6000.0),
    ("high", (4.0, -1.5, 1.0), (0.0, 0.0, 0.0), 2.0, 12000.0),
]
OBSERVER_KMS = (0.0, 0.0, 15000.0)


@pytest.fixture(scope="module")
def browser():
    try:
        playwright = sync_api.sync_playwright().start()
    except Exception as exc:  # noqa: BLE001 -- any failure here means "no browser"
        pytest.skip(f"Playwright could not start ({exc})")
    try:
        chromium = playwright.chromium.launch(args=["--use-gl=swiftshader", "--enable-unsafe-swiftshader"])
    except Exception as exc:  # noqa: BLE001
        playwright.stop()
        pytest.skip(f"No Chromium for Playwright ({str(exc).splitlines()[0]})")
    try:
        yield chromium
    finally:
        chromium.close()
        playwright.stop()


def _packed(observer=OBSERVER_KMS, stars=STARS):
    count = len(stars)
    sky = SimpleNamespace(
        uid=list(range(count)), offset_pc=np.array([s[1] for s in stars], dtype=float),
        velocity_kms=np.array([s[2] for s in stars], dtype=float), dimmed_absolute_magnitude=np.array([s[3] for s in stars], dtype=float),
        temperature_k=np.array([s[4] for s in stars], dtype=float), tier=[next(iter(KIND_CODES))] * count)
    return base64.b64encode(pack_sky(sky, EPOCH, 7.5, observer_velocity_kms=observer)).decode("ascii")


@pytest.fixture
def sky_page(browser):
    page = browser.new_page(viewport={"width": SIZE, "height": SIZE})

    def serve(route):
        name = route.request.url[len(ORIGIN):].split("?")[0]
        if name == "/":
            return route.fulfill(status=200, content_type="text/html", body="<!doctype html><title>sky</title><body></body>")
        path = os.path.join(_STATIC, name[len("/static/"):]) if name.startswith("/static/") else None
        if path and os.path.isfile(path):
            kind = mimetypes.guess_type(path)[0] or "application/octet-stream"
            return route.fulfill(status=200, content_type="text/javascript" if path.endswith(".js") else kind, path=path)
        return route.fulfill(status=404, body="")

    page.route(f"{ORIGIN}/**", serve)
    page.goto(f"{ORIGIN}/")
    yield page
    page.close()


# Draws the packed sky, once per `years`, and returns for every star the spot the twin says (retardedOffset, the camera's projection) and the
# centroid of the light found there (a window of 15 px), with the light's total; and the pixel count of the whole frame that is lit.
DRAW = """async ([b64, yearsList, options]) => {
  const THREE = await import("/static/vendor/three.module.min.js");
  const pack = await import("/static/skypack.js");
  const sphere = await import("/static/skysphere.js");
  const bytes = Uint8Array.from(atob(b64), (c) => c.charCodeAt(0));
  const sky = pack.unpackSky(bytes);
  const canvas = document.createElement("canvas");
  document.body.appendChild(canvas);
  const renderer = new THREE.WebGLRenderer({ canvas, preserveDrawingBuffer: true, antialias: false });
  renderer.setPixelRatio(1);
  renderer.setSize(%(size)d, %(size)d, false);
  renderer.setClearColor(0x000000, 1);
  const camera = new THREE.PerspectiveCamera(60, 1, 0.1, 1000);
  camera.up.set(0, 0, 1);
  camera.position.set(0, 0, 0);
  camera.lookAt(1, 0, 0);
  camera.updateMatrixWorld();
  const scene = new THREE.Scene();
  const drawn = sphere.createSkySphere(THREE, sky, { pixelRatio: 1, ...options });
  scene.add(drawn.object);
  const gl = renderer.getContext();
  const out = { count: drawn.count, frames: [] };
  for (const years of yearsList) {
    drawn.update(camera, { years });
    renderer.render(scene, camera);
    const pixels = new Uint8Array(%(size)d * %(size)d * 4);
    gl.readPixels(0, 0, %(size)d, %(size)d, gl.RGBA, gl.UNSIGNED_BYTE, pixels);
    const light = (x, y) => { const i = 4 * ((%(size)d - 1 - y) * %(size)d + x); return pixels[i] + pixels[i + 1] + pixels[i + 2]; };
    let lit = 0;
    for (let y = 0; y < %(size)d; y++) for (let x = 0; x < %(size)d; x++) if (light(x, y) > 30) lit++;
    const stars = [];
    for (let i = 0; i < sky.count; i++) {
      const twin = sphere.retardedOffset(
        [sky.offsetPc[3 * i], sky.offsetPc[3 * i + 1], sky.offsetPc[3 * i + 2]],
        [sky.velocityKms[3 * i], sky.velocityKms[3 * i + 1], sky.velocityKms[3 * i + 2]], sky.observerVelocityKms, years);
      const ndc = new THREE.Vector3(...twin.apparent).normalize().multiplyScalar(10).project(camera);
      const expected = [(ndc.x + 1) / 2 * %(size)d, (1 - ndc.y) / 2 * %(size)d];
      let sum = 0, sx = 0, sy = 0;
      for (let y = Math.round(expected[1]) - 7; y <= Math.round(expected[1]) + 7; y++) {
        for (let x = Math.round(expected[0]) - 7; x <= Math.round(expected[0]) + 7; x++) {
          if (x < 0 || y < 0 || x >= %(size)d || y >= %(size)d) continue;
          const v = light(x, y);
          sum += v; sx += v * (x + 0.5); sy += v * (y + 0.5);
        }
      }
      stars.push({ index: sky.index[i], expected, found: sum > 0 ? [sx / sum, sy / sum] : null, light: sum, distance: twin.distance, tau: twin.tau });
    }
    out.frames.push({ years, lit, stars });
  }
  drawn.dispose();
  return out;
}""" % {"size": SIZE}


def _draw(page, years, observer=OBSERVER_KMS, **options):
    return page.evaluate(DRAW, [_packed(observer), years, options])


def _near(found, expected, tolerance=1.5):
    return found is not None and abs(found[0] - expected[0]) <= tolerance and abs(found[1] - expected[1]) <= tolerance


def test_every_star_lands_where_the_light_time_arithmetic_puts_it(sky_page):
    """The vertex shader and the JavaScript twin agree to a pixel and a half for a star at rest, one crossing the line of sight at a tenth of
    the speed of light, a slow nearby one and one above the line of sight, with the observer moving; the first sits in the middle of the view."""
    out = _draw(sky_page, [0.0])
    stars = {s["index"]: s for s in out["frames"][0]["stars"]}
    for index, (name, *_rest) in enumerate(STARS):
        if name == "dim":
            continue
        star = stars[index]
        assert star["light"] > 200, (name, star)
        assert _near(star["found"], star["expected"]), (name, star["found"], star["expected"])
    assert _near(stars[0]["found"], [SIZE / 2, SIZE / 2], 1.0)       # the star at rest straight ahead
    # The fast star is visibly shifted by its light time: its apparent place is not where it stands.
    fast = STARS[1]
    at_rest_spot = np.array(fast[1]) / np.linalg.norm(fast[1])
    assert stars[1]["tau"] > 0.9 / 0.3066 * 0.9 and at_rest_spot[1] > 0
    still = sky_page.evaluate("""async ([b64]) => {
      const sphere = await import("/static/skysphere.js");
      const pack = await import("/static/skypack.js");
      const sky = pack.unpackSky(Uint8Array.from(atob(b64), (c) => c.charCodeAt(0)));
      const moving = sphere.retardedOffset([1, 0.3, 0.2], [0, 29000, 0], [0, 0, 15000], 0).apparent;
      const none = sphere.retardedOffset([1, 0.3, 0.2], [0, 0, 0], [0, 0, 0], 0).apparent;
      return [moving, none];
    }""", [_packed()])
    shift = np.degrees(np.arccos(np.dot(still[0], still[1]) / np.linalg.norm(still[0]) / np.linalg.norm(still[1])))
    assert shift > 4.0, shift          # about v / c = 0.097 rad, 5.5 degrees (a third of that is the observer's own motion)


def test_a_star_past_the_limiting_magnitude_is_not_drawn(sky_page):
    out = _draw(sky_page, [0.0])
    dim = out["frames"][0]["stars"][3]
    assert dim["light"] == 0, dim


def test_the_orbit_clock_drifts_a_near_star_and_leaves_a_distant_one_alone(sky_page):
    """A thousand years on, the slow star at 1.6 pc has moved across the view and the star at 10 pc, at rest, has not; the shader still agrees
    with the twin at the later time (the observer is still here: a moving observer sweeps every star)."""
    out = _draw(sky_page, [0.0, 1000.0], observer=(0.0, 0.0, 0.0))
    before, after = (frame["stars"] for frame in out["frames"])
    assert _near(after[2]["found"], after[2]["expected"]) and _near(before[2]["found"], before[2]["expected"])
    moved = np.hypot(after[2]["expected"][0] - before[2]["expected"][0], after[2]["expected"][1] - before[2]["expected"][1])
    assert moved > 10.0, moved
    assert _near(before[0]["found"], after[0]["found"], 0.01)
    assert _near(after[4]["found"], after[4]["expected"])


def test_a_phone_draws_the_brightest_stars_and_the_limit_is_a_setting(sky_page):
    out = sky_page.evaluate(DRAW, [_packed(), [0.0], {"maxStars": 2}])
    assert out["count"] == 2
    lit = [star for star in out["frames"][0]["stars"] if star["light"] > 200]
    assert len(lit) == 2                       # the two nearest in the apparent magnitude at the epoch (the fast star and the drifter)
    dark = sky_page.evaluate(DRAW, [_packed(), [0.0], {"limit": -1.0}])
    assert all(star["light"] == 0 for star in dark["frames"][0]["stars"])      # a limit under their magnitudes draws none of them
