"""
VIEW.36: the nebulae and remnants of the sky a system sees (`static/skysprites.js`) in a real browser. A nebula under the naked-eye limit is a soft disk
in its family's colour, brightest in the middle (the chord of a ball); a supernova remnant is a ring, dark in the middle; a nebula over the limit, a
remnant not yet seen and one that holds the observer are not drawn; a cloud with no light of its own dims the star behind it.

The page is empty and the static files are served to it from disk, so no database or site is needed. Skipped without Playwright or Chromium.
"""

import mimetypes
import os

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

_STATIC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "html", "static")
ORIGIN = "http://sky.test"
SIZE = 400
RADIUS = 0.2          # radians, about 69 px at the camera below


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


def _record(**over):
    record = {"id": 1, "family": "emission", "direction": [1, 0, 0], "angularRadiusRad": RADIUS, "inside": False, "extinctionMag": 0.0,
              "surfaceBrightness": 18.0, "visible": True, "seen": True}
    record.update(over)
    return record


# Draws the records (and, when `withStar`, a bright point of white light dead center drawn between the two sprite layers) and returns the light
# (r, g, b of the pixel) at the offsets in pixels from the center along +x of the view (to the right).
DRAW = """async ([records, probes, withStar]) => {
  const THREE = await import("/static/vendor/three.module.min.js");
  const sprites = await import("/static/skysprites.js");
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
  const drawn = sprites.createSkySprites(THREE, records, {});
  drawn.objects.forEach((object) => scene.add(object));
  if (withStar) {
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(new Float32Array([40, 0, 0]), 3));
    const material = new THREE.PointsMaterial({ color: 0xffffff, size: 20, sizeAttenuation: false, depthTest: false, depthWrite: false });
    const star = new THREE.Points(geometry, material);
    star.renderOrder = -10;
    star.frustumCulled = false;
    scene.add(star);
  }
  drawn.update(camera);
  renderer.render(scene, camera);
  const gl = renderer.getContext();
  const pixels = new Uint8Array(%(size)d * %(size)d * 4);
  gl.readPixels(0, 0, %(size)d, %(size)d, gl.RGBA, gl.UNSIGNED_BYTE, pixels);
  const at = ([dx, dy]) => { const i = 4 * ((%(size)d / 2 - dy) * %(size)d + (%(size)d / 2 + dx)); return [pixels[i], pixels[i + 1], pixels[i + 2]]; };
  const out = { count: drawn.count, light: probes.map(at) };
  drawn.dispose();
  return out;
}""" % {"size": SIZE}


def _draw(page, records, probes, with_star=False):
    return page.evaluate(DRAW, [records, probes, with_star])


def _sum(rgb):
    return sum(rgb)


def test_a_nebula_is_a_soft_disk_brightest_in_the_middle(sky_page):
    out = _draw(sky_page, [_record()], [[0, 0], [40, 0], [64, 0], [120, 0]])
    center, mid, rim, outside = (_sum(p) for p in out["light"])
    assert out["count"] == {"glow": 1, "dark": 0}
    assert center > 40
    assert center > mid > 0
    assert rim < mid
    assert outside == 0
    # Emission is pink: red leads blue.
    assert out["light"][0][0] > out["light"][0][2] > 0


def test_a_remnant_is_a_ring_dark_in_the_middle(sky_page):
    out = _draw(sky_page, [_record(family="supernova-remnant", surfaceBrightness=18.0)], [[0, 0], [30, 0], [66, 0], [120, 0]])
    center, inner, ring, outside = (_sum(p) for p in out["light"])
    assert ring > 3 * max(center, 1)
    assert ring > inner
    assert outside == 0
    assert out["light"][2][0] > out["light"][2][2]          # orange: red over blue


def test_a_faint_nebula_a_remnant_not_yet_seen_and_one_holding_the_observer_are_not_drawn(sky_page):
    records = [_record(visible=False, surfaceBrightness=26.0), _record(family="supernova-remnant", seen=False), _record(inside=True)]
    out = _draw(sky_page, records, [[0, 0], [60, 0]])
    assert out["count"] == {"glow": 0, "dark": 0}
    assert [_sum(p) for p in out["light"]] == [0, 0]


def test_a_dark_cloud_dims_the_star_behind_it_and_leaves_the_sky_beside_it_alone(sky_page):
    cloud = _record(family="dark", surfaceBrightness=None, visible=False, extinctionMag=3.0)
    bare = _draw(sky_page, [], [[0, 0]], with_star=True)
    dimmed = _draw(sky_page, [cloud], [[0, 0], [120, 0]], with_star=True)
    assert dimmed["count"] == {"glow": 0, "dark": 1}
    assert _sum(bare["light"][0]) > 600
    # A = 3 mag through the middle: 1 - 10^(-1.2) = 94 % of the light is gone.
    assert _sum(dimmed["light"][0]) < 0.25 * _sum(bare["light"][0])
    assert _sum(dimmed["light"][1]) == 0
