"""
VIEW.35: the controller of the System View's sky (`static/systemsky.js`) in a real browser, against stand-in routes. The packed sky and the glow image come
from the real writers (`pack_sky`, half floats), a view that only adds the layers to a scene stands in for `createSystemView`, and the camera looks along the
scene's x axis. Checked: the sky is drawn (a star straight ahead, the glow band along the galaxy's equator laid in the scene's reference plane, so it runs
sideways and not up), the status line says what is drawn, a body in front of a star hides it, the sky goes when it is turned off and comes back when it is
turned on, a sky still filling in is asked for again, and a sky the route refuses leaves nothing attached.

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
WIDTH, HEIGHT = 64, 32


def _band_bytes():
    """A glow image of 2.0 within 8 degrees of the equator and nothing else, as the server sends it (`rgb16f`, row 0 at the north pole)."""
    lat = (0.5 - (np.arange(HEIGHT) + 0.5) / HEIGHT) * 180.0
    rows = np.where(np.abs(lat) < 8.0, 2.0, 0.0).astype(np.float16)
    image = np.repeat(rows[:, None, None], WIDTH, axis=1).repeat(3, axis=2)
    return image.astype("<f2").tobytes()


def _packed(report):
    stars = [((10.0, 0.0, 0.0), 0.0, 6000.0)]          # straight ahead, 10 pc, absolute magnitude 0
    sky = SimpleNamespace(
        uid=[0], offset_pc=np.array([s[0] for s in stars], dtype=float), velocity_kms=np.zeros((1, 3)),
        dimmed_absolute_magnitude=np.array([s[1] for s in stars], dtype=float), temperature_k=np.array([s[2] for s in stars], dtype=float),
        tier=[next(iter(KIND_CODES))])
    meta = {"system": "A", "nebulae": [], "report": report, "band": {"key": "k" * 32, "width": WIDTH, "height": HEIGHT, "format": "rgb16f",
                                                                      "href": "/system/A/sky/band?limit=7.5"}}
    return pack_sky(sky, EPOCH, 7.5, meta=meta)


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


def _page(browser, report, sky_status=200):
    page = browser.new_page(viewport={"width": SIZE, "height": SIZE})
    requests = []

    def serve(route):
        name = route.request.url[len(ORIGIN):]
        path_only = name.split("?")[0]
        requests.append(name)
        if path_only == "/":
            return route.fulfill(status=200, content_type="text/html", body="<!doctype html><title>sky</title><body></body>")
        if path_only == "/system/A/sky":
            if sky_status != 200:
                return route.fulfill(status=sky_status, content_type="application/json", body='{"error": "no"}')
            return route.fulfill(status=200, content_type="application/vnd.planetgen.sky", body=_packed(report))
        if path_only == "/system/A/sky/band":
            return route.fulfill(status=200, content_type="application/octet-stream", body=_band_bytes())
        path = os.path.join(_STATIC, path_only[len("/static/"):]) if path_only.startswith("/static/") else None
        if path and os.path.isfile(path):
            kind = mimetypes.guess_type(path)[0] or "application/octet-stream"
            return route.fulfill(status=200, content_type="text/javascript" if path.endswith(".js") else kind, path=path)
        return route.fulfill(status=404, body="")

    page.route(f"{ORIGIN}/**", serve)
    page.goto(f"{ORIGIN}/")
    page.requests = requests
    return page


# Builds the stand-in view and the controller, loads the sky, renders once with an optional blocker sphere and reads probe pixels (offsets from the middle,
# y up). Leaves `window.sky`, `window.shot(probes, withBlocker)` and `window.view` for the steps after.
SETUP = """async ([options]) => {
  const THREE = await import("/static/vendor/three.module.min.js");
  const mod = await import("/static/systemsky.js");
  const canvas = document.createElement("canvas");
  document.body.appendChild(canvas);
  const renderer = new THREE.WebGLRenderer({ canvas, preserveDrawingBuffer: true, antialias: false });
  renderer.setPixelRatio(1);
  renderer.setSize(%(size)d, %(size)d, false);
  renderer.setClearColor(0x000000, 1);
  const camera = new THREE.PerspectiveCamera(60, 1, 0.1, 100000);
  camera.up.set(0, 1, 0);
  camera.position.set(0, 0, 0);
  camera.lookAt(1, 0, 0);
  camera.updateMatrixWorld();
  const stage = new THREE.Scene();
  const blocker = new THREE.Mesh(new THREE.SphereGeometry(1, 24, 12), new THREE.MeshBasicMaterial({ color: 0x224488 }));
  blocker.position.set(5, 0, 0);
  blocker.visible = false;
  stage.add(blocker);
  let layer = null;
  const view = {
    attachSky(next) { layer = next; next.objects.forEach((o) => stage.add(o)); },
    detachSky() { if (layer) layer.objects.forEach((o) => stage.remove(o)); layer = null; },
    attached: () => layer !== null,
  };
  const notes = [];
  const storage = { data: {}, getItem(k) { return k in this.data ? this.data[k] : null; }, setItem(k, v) { this.data[k] = v; } };
  window.timers = [];
  const sky = mod.createSystemSky(THREE, {
    skyUrl: "/system/A/sky", view, sceneEpochUnix: %(epoch)r, storage, onState: (s, n) => notes.push([s.kind, n]),
    setTimeout: (fn, ms) => { window.timers.push([fn, ms]); return window.timers.length; }, clearTimeout: () => {}, ...options,
  });
  window.sky = sky; window.view = view; window.notes = notes;
  window.shot = (probes, withBlocker) => {
    blocker.visible = !!withBlocker;
    if (layer) layer.update(camera, { years: 0, pixelRatio: 1 });
    renderer.render(stage, camera);
    const gl = renderer.getContext();
    const pixels = new Uint8Array(%(size)d * %(size)d * 4);
    gl.readPixels(0, 0, %(size)d, %(size)d, gl.RGBA, gl.UNSIGNED_BYTE, pixels);
    return probes.map(([dx, dy]) => {
      let total = 0;
      for (let y = -3; y <= 3; y++) for (let x = -3; x <= 3; x++) {
        const i = 4 * ((%(size)d / 2 - (dy + y)) * %(size)d + (%(size)d / 2 + dx + x));
        total += pixels[i] + pixels[i + 1] + pixels[i + 2];
      }
      return total;
    });
  };
  return null;
}""" % {"size": SIZE, "epoch": EPOCH}


def test_the_sky_is_drawn_and_the_galaxys_equator_runs_sideways(browser):
    page = _page(browser, {"filling": False, "built_to_pc": 900.0})
    page.evaluate(SETUP, [{}])
    state = page.evaluate("async () => { const s = await window.sky.reload(); return [s.kind, s.count, s.band]; }")
    assert state == ["ready", 1, True]
    assert page.evaluate("() => window.view.attached()") is True
    middle, right, up = page.evaluate("() => window.shot([[0, 0], [90, 0], [0, 90]], false)")
    assert middle > 150                       # the star ahead and the glow under it
    assert right > 30                         # the glow runs sideways: the galaxy's equator lies in the scene's reference plane
    assert up == 0                            # and not up
    assert page.evaluate("() => window.notes.at(-1)[1]") .startswith("Sky of 1 stars")
    assert any(r.startswith("/system/A/sky/band?limit=7.5&key=") for r in page.requests)


def test_a_body_in_front_of_a_star_hides_it(browser):
    page = _page(browser, {"filling": False})
    page.evaluate(SETUP, [{}])
    page.evaluate("() => window.sky.reload()")
    open_sky = page.evaluate("() => window.shot([[0, 0]], false)")[0]
    blocked = page.evaluate("() => window.shot([[0, 0]], true)")[0]
    assert open_sky > 150
    lit = page.evaluate("() => window.shot([[0, 0]], true)")[0]
    assert blocked == lit and blocked > 0                  # the body's own colour, no star light added over it
    assert blocked != open_sky


def test_turning_the_sky_off_removes_it_and_on_brings_it_back(browser):
    page = _page(browser, {"filling": False})
    page.evaluate(SETUP, [{}])
    page.evaluate("() => window.sky.reload()")
    page.evaluate("() => window.sky.disable()")
    assert page.evaluate("() => window.view.attached()") is False
    assert page.evaluate("() => window.shot([[0, 0]], false)")[0] == 0
    assert page.evaluate("() => window.notes.at(-1)[0]") == "off"
    page.evaluate("async () => { await window.sky.enable(); }")
    assert page.evaluate("() => window.view.attached()") is True
    assert page.evaluate("() => window.shot([[0, 0]], false)")[0] > 150


def test_a_sky_still_filling_in_says_so_and_is_asked_for_again(browser):
    page = _page(browser, {"filling": True, "built_to_pc": 120.0})
    page.evaluate(SETUP, [{}])
    page.evaluate("() => window.sky.reload()")
    assert page.evaluate("() => window.notes.at(-1)") [0] == "filling"
    assert "built to 120 pc, filling in" in page.evaluate("() => window.notes.at(-1)[1]")
    assert page.evaluate("() => window.timers.length") == 1 and page.evaluate("() => window.timers[0][1]") == 30000
    before = len(page.requests)
    page.evaluate("async () => { window.timers[0][0](); await new Promise((r) => setTimeout(r, 400)); }")
    assert len([r for r in page.requests[before:] if r.startswith("/system/A/sky?")]) == 1
    assert page.evaluate("() => window.view.attached()") is True


def test_a_sky_the_route_refuses_leaves_the_plain_background(browser):
    page = _page(browser, {}, sky_status=404)
    page.evaluate(SETUP, [{}])
    page.evaluate("() => window.sky.reload()")
    assert page.evaluate("() => window.view.attached()") is False
    assert "no place in the galaxy" in page.evaluate("() => window.notes.at(-1)[1]")
