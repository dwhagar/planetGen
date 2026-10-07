// tests/js/sectormap.test.mjs -- static/sectormap.js, the Sector Map
// (TEST.58): its zoom buttons, wheel and their clamp, Reset, turning by
// drag and arrow keys (the tilt clamp), and picking: a click on a star,
// the screen-reader list, the "Show on map" buttons, and the click a drag
// swallows. The scene is lib/starmap.py's real JSON (test_js_unit.py's
// fixtures) drawn by a renderer that draws nothing; picking is real
// three.js projection and raycasting.

import { test } from "node:test";
import assert from "node:assert/strict";

import { FakeEvent, h, installDom, STATIC_URL } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

installDom("http://localhost/sector/1");
const THREE = await import(new URL("vendor/three.module.min.js", STATIC_URL).href);
const SM = await import(new URL("sectormap.js", STATIC_URL).href);

const SIZE = 400;
const FOV = 45;
const MAX_ZOOM = 2.5;

// The zoom floor: 0.2, or half the opening zoom when a crowded sector
// opens near it, so - always has somewhere to go (MAP.114).
function minZoom(m) {
  const opening = m.data.defaultZoom > 0 && m.data.defaultZoom <= 1 ? m.data.defaultZoom : 1;
  return Math.min(0.2, opening / 2);
}

function setUp(data) {
  const win = installDom("http://localhost/sector/1");
  data = structuredClone(data || F.sectorMap.data);
  const canvasEl = h("canvas", { id: "starmap-canvas", tabindex: "0" });
  canvasEl.rect = { left: 0, top: 0, width: SIZE, height: SIZE };
  canvasEl.clientWidth = SIZE;
  canvasEl.clientHeight = SIZE;
  const viewport = h("div", { class: "starmap-viewport" }, [canvasEl, h("div", { id: "starmap-scale" }, [
    h("span", { id: "starmap-scale-bar" }), h("span", { id: "starmap-scale-label" }),
  ])]);
  const controls = h("div", { id: "starmap-controls" },
    F.sectorMap.actions.map((action) => h("button", { type: "button", "data-action": action, "aria-pressed": action === "toggle-rogue-markers" ? "false" : null }, action)));
  const info = h("div", { id: "starmap-info" });
  const showButtons = h("ul", {}, [
    h("button", { "data-map-target": "nebula:1", hidden: true }, "Show on map"),
    h("button", { "data-map-target": "nebula:999", hidden: true }, "Show on map"),
  ]);
  document.body.append(viewport, controls, info, showButtons);
  const renderer = {
    calls: 0, camera: null, scene: null, pixelRatio: 1,
    setPixelRatio(r) { this.pixelRatio = r; },
    getPixelRatio() { return this.pixelRatio; },
    setClearColor() {},
    setSize() {},
    render(scene, camera) { this.calls += 1; this.scene = scene; this.camera = camera; },
  };
  SM.initStarmap(canvasEl, data, { renderer });
  const camera = renderer.camera;
  const referenceDistance = data.sceneHalfPx / Math.tan(THREE.MathUtils.degToRad(FOV / 2));
  return {
    win, data, canvasEl, controls, info, showButtons, renderer, camera,
    zoom: () => referenceDistance / camera.position.length(),
    polar: () => Math.acos(camera.position.y / camera.position.length()),
    azimuth: () => Math.atan2(camera.position.x, camera.position.z),
    button: (action) => controls.querySelector(`[data-action="${action}"]`),
  };
}

function near(actual, expected, what) {
  assert.ok(Math.abs(actual - expected) < 1e-9 * Math.max(1, Math.abs(expected)), `${what || ""} ${actual} != ${expected}`);
}

// Where a scene point lands on the canvas.
function screenOf(m, x, y, z) {
  const p = new THREE.Vector3(x, y, z).project(m.camera);
  return { x: ((p.x + 1) / 2) * SIZE, y: ((1 - p.y) / 2) * SIZE };
}

function pointer(m, type, x, y) {
  m.canvasEl.dispatchEvent(new FakeEvent(type, { clientX: x, clientY: y, pointerId: 1, button: 0 }));
}

test("every control button the panel draws is wired", () => {
  const m = setUp();
  assert.deepEqual(F.sectorMap.actions.sort(), ["reset", "toggle-rogue-markers", "zoom-in", "zoom-out"]);
  for (const action of F.sectorMap.actions) {
    if (action === "reset") m.button("zoom-in").click();
    if (action === "zoom-out") m.button("zoom-in").click();
    const before = [m.camera.position.toArray().join(), m.button(action).getAttribute("aria-pressed")].join("|");
    m.button(action).click();
    const after = [m.camera.position.toArray().join(), m.button(action).getAttribute("aria-pressed")].join("|");
    assert.notEqual(after, before, `${action} changed nothing`);
  }
});

test("opens at the panel's default zoom", () => {
  const m = setUp();
  near(m.zoom(), m.data.defaultZoom > 0 && m.data.defaultZoom <= 1 ? m.data.defaultZoom : 1, "zoom");
  near(m.polar(), THREE.MathUtils.degToRad(72), "polar");
});

test("+ and - step the zoom by 0.15 and stop at the floor and 2.5", () => {
  const m = setUp();
  const start = m.zoom();
  m.button("zoom-in").click();
  near(m.zoom(), start + 0.15);
  m.button("zoom-in").click();
  m.button("zoom-out").click();
  near(m.zoom(), start + 0.15);
  for (let n = 0; n < 40; n++) m.button("zoom-in").click();
  near(m.zoom(), MAX_ZOOM);
  for (let n = 0; n < 40; n++) m.button("zoom-out").click();
  near(m.zoom(), minZoom(m));
});

test("- zooms out past the opening view, even for a sector that opens at 0.2", () => {
  const data = structuredClone(F.sectorMap.data);
  data.defaultZoom = 0.2;
  const m = setUp(data);
  m.button("zoom-out").click();
  assert.ok(m.zoom() < 0.2 - 1e-9, `zoom ${m.zoom()} is not below the opening 0.2`);
  m.button("reset").click();
  near(m.zoom(), 0.2);
});

test("the wheel zooms by 0.08 a notch and keeps the page still", () => {
  const m = setUp();
  const start = m.zoom();
  const wheel = new FakeEvent("wheel", { deltaY: -120 });
  m.canvasEl.dispatchEvent(wheel);
  assert.ok(wheel.defaultPrevented);
  near(m.zoom(), start + 0.08);
  for (let n = 0; n < 60; n++) m.canvasEl.dispatchEvent(new FakeEvent("wheel", { deltaY: 120 }));
  near(m.zoom(), minZoom(m));
});

test("Reset view puts back the zoom and the turn", () => {
  const m = setUp();
  const start = m.camera.position.clone();
  m.button("zoom-in").click();
  pointer(m, "pointerdown", 200, 200);
  pointer(m, "pointermove", 260, 230);
  pointer(m, "pointerup", 260, 230);
  assert.ok(m.camera.position.distanceTo(start) > 1);
  m.button("reset").click();
  assert.ok(m.camera.position.distanceTo(start) < 1e-9 * start.length());
});

test("dragging turns the view; the tilt stops 2 degrees short of either pole", () => {
  const m = setUp();
  const azimuth = m.azimuth();
  pointer(m, "pointerdown", 200, 200);
  pointer(m, "pointermove", 250, 200);
  assert.ok(Math.abs(m.azimuth() - azimuth) > 0.2, "turned about the vertical");
  pointer(m, "pointermove", 250, -5000);
  near(m.polar(), THREE.MathUtils.degToRad(178), "dragged up: from below");
  pointer(m, "pointermove", 250, 50000);
  near(m.polar(), THREE.MathUtils.degToRad(2), "dragged down: from above");
  pointer(m, "pointerup", 250, 50000);
});

test("the arrow keys turn the view by 6 degrees and keep the same clamp", () => {
  const m = setUp();
  const polar = m.polar();
  const up = new FakeEvent("keydown", { key: "ArrowUp" });
  m.canvasEl.dispatchEvent(up);
  assert.ok(up.defaultPrevented);
  near(m.polar(), polar - THREE.MathUtils.degToRad(6));
  for (let n = 0; n < 40; n++) m.canvasEl.dispatchEvent(new FakeEvent("keydown", { key: "ArrowUp" }));
  near(m.polar(), THREE.MathUtils.degToRad(2));
  const azimuth = m.azimuth();
  m.canvasEl.dispatchEvent(new FakeEvent("keydown", { key: "ArrowLeft" }));
  assert.ok(Math.abs(m.azimuth() - azimuth) > 0.05);
  const other = new FakeEvent("keydown", { key: "a" });
  m.canvasEl.dispatchEvent(other);
  assert.equal(other.defaultPrevented, false, "other keys are left alone");
});

test("clicking a star shows it; clicking empty space changes nothing", () => {
  const m = setUp();
  const star = m.data.stars[0];
  const at = screenOf(m, star.x, star.y, star.z);
  m.canvasEl.dispatchEvent(new FakeEvent("click", { clientX: at.x, clientY: at.y }));
  assert.equal(m.info.querySelector("h3").textContent, star.name);
  assert.ok(m.info.querySelector(`a[href="${star.href}"]`), "with a link to the system page");
  m.canvasEl.dispatchEvent(new FakeEvent("click", { clientX: 1, clientY: 1 }));
  assert.equal(m.info.querySelector("h3").textContent, star.name, "still showing the star");
});

test("each star stays clickable as the view zooms", () => {
  const m = setUp();
  for (const step of ["zoom-in", "zoom-in", "zoom-out", "zoom-out", "zoom-out"]) {
    m.button(step).click();
    for (const star of m.data.stars) {
      const at = screenOf(m, star.x, star.y, star.z);
      m.canvasEl.dispatchEvent(new FakeEvent("click", { clientX: at.x, clientY: at.y }));
      assert.equal(m.info.querySelector("h3").textContent, star.name, `${star.name} after ${step}`);
    }
  }
});

test("the click that ends a drag picks nothing", () => {
  const m = setUp();
  const star = m.data.stars[0];
  pointer(m, "pointerdown", 10, 10);
  pointer(m, "pointermove", 40, 10);
  pointer(m, "pointerup", 40, 10);
  const at = screenOf(m, star.x, star.y, star.z);
  m.canvasEl.dispatchEvent(new FakeEvent("click", { clientX: at.x, clientY: at.y }));
  assert.equal(m.info.querySelector("h3"), null);
  m.canvasEl.dispatchEvent(new FakeEvent("click", { clientX: at.x, clientY: at.y }));
  assert.equal(m.info.querySelector("h3").textContent, star.name, "the next click picks again");
});

test("the screen-reader list has a button per entry that selects it", () => {
  const m = setUp();
  const list = document.querySelector("ul.starmap-sr-list");
  const buttons = list.querySelectorAll("button");
  const entries = m.data.stars.concat(m.data.clouds).concat(m.data.neighbors);
  assert.equal(buttons.length, entries.length);
  buttons.forEach((button, n) => {
    button.click();
    const heading = m.info.querySelector("h3");
    assert.ok(heading, `entry ${n} shows something`);
  });
  const cloud = m.data.clouds.findIndex((c) => c.kind === "nebula");
  buttons[m.data.stars.length + cloud].click();
  assert.equal(m.info.querySelector("h3").textContent, m.data.clouds[cloud].name);
});

test("Show on map buttons appear only for clouds on the map, and select them", () => {
  const m = setUp();
  const [known, unknown] = m.showButtons.querySelectorAll("button");
  assert.equal(known.hidden, false);
  assert.equal(unknown.hidden, true);
  known.click();
  assert.equal(m.info.querySelector("h3").textContent, m.data.clouds.find((c) => c.key === "nebula:1").name);
  assert.equal(document.activeElement, m.canvasEl, "focus moves to the map");
});

test("the rogue-planet toggle starts off, then marks them: rings, a bigger point, an easier pick", () => {
  const m = setUp();
  const toggle = m.button("toggle-rogue-markers");
  const markers = () => m.renderer.scene.children.filter((o) => o.isGroup && o.children.some((c) => c.isSprite && c.material.sizeAttenuation === false));
  const points = m.renderer.scene.children.find((o) => o.isPoints);
  const rogue = m.data.clouds.find((c) => c.kind === "roguePlanet");
  const index = m.data.stars.length + m.data.clouds.filter((c) => c.light).indexOf(rogue);
  const size = () => points.geometry.getAttribute("pointSize").array[index];
  const glow = () => points.geometry.getAttribute("pointGlow").array[index];
  assert.equal(markers().length, 1, "one group of rogue markers");
  assert.equal(markers()[0].visible, false, "off by default");
  assert.equal(size(), rogue.light.sizePx);
  assert.equal(glow(), 0, "no glow unmarked");

  // A click a few pixels off an unmarked rogue planet misses it.
  const at = screenOf(m, rogue.x, rogue.y, rogue.z);
  m.canvasEl.dispatchEvent(new FakeEvent("click", { clientX: at.x + 7, clientY: at.y }));
  assert.notEqual((m.info.querySelector("h3") || {}).textContent, rogue.name);

  toggle.click();
  assert.equal(toggle.getAttribute("aria-pressed"), "true");
  assert.equal(markers()[0].visible, true);
  assert.equal(size(), rogue.markedLight.sizePx);
  assert.ok(glow() > 0, "glows marked");
  m.canvasEl.dispatchEvent(new FakeEvent("click", { clientX: at.x + 7, clientY: at.y }));
  assert.equal(m.info.querySelector("h3").textContent, rogue.name, "marked, the same click takes it");

  toggle.click();
  assert.equal(toggle.getAttribute("aria-pressed"), "false");
  assert.equal(markers()[0].visible, false);
  assert.equal(size(), rogue.light.sizePx);
});

test("the scale bar follows the zoom", () => {
  const m = setUp();
  const label = document.getElementById("starmap-scale-label");
  const bar = document.getElementById("starmap-scale-bar");
  assert.ok(label.textContent, "a scale label");
  const before = label.textContent + bar.style.width;
  m.button("zoom-in").click();
  m.button("zoom-in").click();
  assert.notEqual(label.textContent + bar.style.width, before);
});
