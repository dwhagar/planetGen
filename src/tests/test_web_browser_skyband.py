"""
VIEW.34: the glow band of the sky a system sees (`static/skyband.js`) in a real browser. A glow image is sampled by direction on a sphere round the
camera: a bright stripe along the galactic equator shows where the map says and nowhere else, the middle of the view reads the display law's level
for the value the image holds there, the grain varies the light from cell to cell by no more than its strength and changes with the seed, and with
no grain the dither leaves the band flat to within one level.

The page is empty and the static files are served to it from disk, so no database or site is needed. Skipped without Playwright or Chromium.
"""

import mimetypes
import os

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

_STATIC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "html", "static")
ORIGIN = "http://sky.test"
SIZE = 400
STRIPE = (2.0, 1.0, 0.5)      # the value of the stripe along the equator, red green blue


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


# Draws a band of `kind` ("stripe": STRIPE within 8 degrees of the equator, 0 elsewhere; "flat": STRIPE everywhere) and returns, per probe
# (a window of 21 x 21 px about an offset (dx, dy) from the middle, y up), the mean and the standard deviation of each channel in 0..255, and the
# twin's answer for the middle: `displayLevel` of the sampled value times 255.
DRAW = """async ([kind, probes, options]) => {
  const THREE = await import("/static/vendor/three.module.min.js");
  const band = await import("/static/skyband.js");
  const stripe = %(stripe)s;
  const image = band.makeImage(256, 128, (d) => kind === "flat" || Math.abs(Math.asin(d[2])) < 8 * Math.PI / 180 ? stripe : [0, 0, 0]);
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
  const drawn = band.createSkyBand(THREE, image, options || {});
  scene.add(drawn.object);
  drawn.update(camera);
  renderer.render(scene, camera);
  const gl = renderer.getContext();
  const pixels = new Uint8Array(%(size)d * %(size)d * 4);
  gl.readPixels(0, 0, %(size)d, %(size)d, gl.RGBA, gl.UNSIGNED_BYTE, pixels);
  const half = 10;
  const out = { probes: probes.map(([dx, dy]) => {
    const sums = [0, 0, 0], squares = [0, 0, 0];
    let n = 0;
    for (let y = -half; y <= half; y++) for (let x = -half; x <= half; x++) {
      const i = 4 * ((%(size)d / 2 - (dy + y)) * %(size)d + (%(size)d / 2 + dx + x));
      for (let k = 0; k < 3; k++) { sums[k] += pixels[i + k]; squares[k] += pixels[i + k] * pixels[i + k]; }
      n++;
    }
    const mean = sums.map((s) => s / n);
    return { mean, sd: squares.map((q, k) => Math.sqrt(Math.max(0, q / n - mean[k] * mean[k]))) };
  }), twin: stripe.map((v) => band.displayLevel(v) * 255) };
  drawn.dispose();
  return out;
}""" % {"size": SIZE, "stripe": list(STRIPE)}


def _draw(page, kind, probes, **options):
    return page.evaluate(DRAW, [kind, probes, options])


def test_the_stripe_shows_on_the_equator_and_nowhere_else(sky_page):
    # 20 degrees above the middle is 74 px up (focal length 346 px), well outside the 8 degree stripe; the equator is dead ahead.
    out = _draw(sky_page, "stripe", [[0, 0], [0, 80], [0, -80]], grain=0.0)
    middle, above, below = out["probes"]
    assert max(above["mean"]) < 1.0 and max(below["mean"]) < 1.0
    for channel in range(3):
        assert abs(middle["mean"][channel] - out["twin"][channel]) <= 1.5, (channel, middle["mean"], out["twin"])
    assert middle["mean"][0] > middle["mean"][1] > middle["mean"][2] > 5          # the stripe's colour, red over green over blue


def test_with_no_grain_the_dither_leaves_the_band_flat_to_one_level(sky_page):
    out = _draw(sky_page, "flat", [[0, 0], [100, 60], [-120, -90]], grain=0.0)
    for probe in out["probes"]:
        for channel in range(3):
            assert abs(probe["mean"][channel] - out["twin"][channel]) <= 1.5, (probe, out["twin"])
            assert probe["sd"][channel] <= 1.0, probe


def test_the_grain_varies_the_light_within_its_strength_and_changes_with_the_seed(sky_page):
    flat = _draw(sky_page, "flat", [[0, 0]], grain=0.0)["probes"][0]
    grainy = _draw(sky_page, "flat", [[0, 0]], grain=0.2, seed=1)["probes"][0]
    other = _draw(sky_page, "flat", [[0, 0]], grain=0.2, seed=2)["probes"][0]
    again = _draw(sky_page, "flat", [[0, 0]], grain=0.2, seed=1)["probes"][0]
    assert grainy["sd"][0] > 2.0 * max(flat["sd"][0], 0.3)            # cells of 0.12 degrees are 2 px across: the patch shows many
    assert grainy["sd"][0] < 0.2 * grainy["mean"][0]                   # but never beyond the strength
    assert abs(grainy["mean"][0] - flat["mean"][0]) < 0.1 * flat["mean"][0]
    assert again == grainy                                             # the same seed draws the same sky
    assert other["mean"] != grainy["mean"]                             # another seed does not
