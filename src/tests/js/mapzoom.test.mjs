// tests/js/mapzoom.test.mjs -- static/mapzoom.js, the viewBox zoom and pan
// under the phenomenon diagram (TEST.58): the buttons, the wheel, the
// zoom clamp (zoomedBox), drag to pan and the click a drag swallows.

import { test } from "node:test";
import assert from "node:assert/strict";

import { FakeEvent, h, importPage, installDom } from "./fakedom.mjs";

const SIZE_PX = 400;

// An <svg> laid out SIZE_PX square at the page's top left, with the
// screen transform a browser would give it (viewBox scaled to fit).
function svg(viewBox, data) {
  const el = h("svg", Object.assign({ viewBox: viewBox }, data || {}));
  el.rect = { left: 0, top: 0, width: SIZE_PX, height: SIZE_PX };
  el.createSVGPoint = () => ({
    x: 0, y: 0,
    matrixTransform(m) { return { x: m.a * this.x + m.e, y: m.d * this.y + m.f }; },
  });
  el.getScreenCTM = () => {
    const [x, y, w, hgt] = el.getAttribute("viewBox").split(/\s+/).map(Number);
    const sx = SIZE_PX / w;
    const sy = SIZE_PX / hgt;
    const m = { a: sx, d: sy, e: -x * sx, f: -y * sy };
    m.inverse = () => ({ a: 1 / sx, d: 1 / sy, e: x, f: y });
    return m;
  };
  document.body.appendChild(el);
  return el;
}

function box(el) {
  const [x, y, w, hgt] = el.getAttribute("viewBox").split(/\s+/).map(Number);
  return { x, y, w, h: hgt };
}

function close(actual, expected, message) {
  for (const key of Object.keys(expected)) {
    assert.ok(Math.abs(actual[key] - expected[key]) <= 1e-9 * Math.max(1, Math.abs(expected[key])),
      `${message || ""} ${key}: ${actual[key]} != ${expected[key]}`);
  }
}

async function setUp(viewBox, data, options) {
  installDom("http://localhost/phenomenon/nebula/1");
  await importPage("mapzoom.js");
  const el = svg(viewBox, data);
  const buttons = { zoomInBtn: h("button"), zoomOutBtn: h("button"), resetBtn: h("button") };
  const changes = [];
  const controller = window.planetgenInitSvgZoomPan(el, Object.assign({}, buttons, {
    onChange: (b) => changes.push(b),
  }, options || {}));
  return { el, controller, changes, ...buttons };
}

function pointer(el, type, x, y, extra) {
  const event = new FakeEvent(type, Object.assign({ clientX: x, clientY: y, pointerId: 1, button: 0 }, extra || {}));
  el.dispatchEvent(event);
  return event;
}

test("the + and - buttons zoom by 1.4 about the middle, Reset goes back", async () => {
  const { el, changes, zoomInBtn, zoomOutBtn, resetBtn } = await setUp("-50 -50 100 100", { "data-min-view-size": "1", "data-max-view-size": "1000" });
  zoomInBtn.click();
  close(box(el), { x: -50 / 1.4, y: -50 / 1.4, w: 100 / 1.4, h: 100 / 1.4 }, "zoomed in");
  zoomOutBtn.click();
  zoomOutBtn.click();
  close(box(el), { x: -50 * 1.4, y: -50 * 1.4, w: 140, h: 140 }, "zoomed out");
  resetBtn.click();
  assert.deepEqual(box(el), { x: -50, y: -50, w: 100, h: 100 });
  assert.equal(changes.length, 4, "onChange runs after every move");
  assert.deepEqual(changes[3], { x: -50, y: -50, w: 100, h: 100 });
});

test("zooming in stops at data-min-view-size and out at data-max-view-size", async () => {
  const { el, zoomInBtn, zoomOutBtn } = await setUp("-5 -5 10 10", { "data-min-view-size": "2", "data-max-view-size": "40" });
  for (let n = 0; n < 20; n++) zoomInBtn.click();
  close(box(el), { x: -1, y: -1, w: 2, h: 2 }, "deepest zoom");
  for (let n = 0; n < 40; n++) zoomOutBtn.click();
  close(box(el), { x: -20, y: -20, w: 40, h: 40 }, "widest zoom");
});

test("without data-max-view-size the starting view is the widest", async () => {
  const { el, zoomInBtn, zoomOutBtn } = await setUp("0 0 100 100", { "data-min-view-size": "10" });
  zoomOutBtn.click();
  assert.deepEqual(box(el), { x: 0, y: 0, w: 100, h: 100 }, "- does nothing at the start");
  zoomInBtn.click();
  zoomOutBtn.click();
  close(box(el), { x: 0, y: 0, w: 100, h: 100 });
});

test("a view that opens at the zoom-out limit can't zoom out (the large nebula case)", async () => {
  // lib/phenomenonmap.py caps the default view at the 1 ly maximum, so a
  // nebula 0.5 ly or more across opens here; TEST.55 checks the page.
  const ly = 63241.077;
  const { el, zoomOutBtn } = await setUp(`${-ly / 2} ${-ly / 2} ${ly} ${ly}`, { "data-min-view-size": "1", "data-max-view-size": String(ly) });
  const before = box(el);
  zoomOutBtn.click();
  close(box(el), before);
});

test("options override the data attributes", async () => {
  const { el, zoomInBtn } = await setUp("0 0 100 100", { "data-min-view-size": "1" }, { minSize: 50 });
  for (let n = 0; n < 5; n++) zoomInBtn.click();
  close(box(el), { x: 25, y: 25, w: 50, h: 50 });
});

test("a non-square view keeps its shape when clamped", async () => {
  const { el, zoomInBtn } = await setUp("0 0 200 100", { "data-min-view-size": "50", "data-max-view-size": "400" });
  for (let n = 0; n < 10; n++) zoomInBtn.click();
  close(box(el), { x: 75, y: 37.5, w: 50, h: 25 });
});

test("the wheel zooms about the point under the pointer", async () => {
  const { el } = await setUp("0 0 100 100", { "data-min-view-size": "1", "data-max-view-size": "1000" });
  // 100 px across 400 px is SVG x = 25, y = 75.
  const event = new FakeEvent("wheel", { clientX: 100, clientY: 300, deltaY: -1 });
  el.dispatchEvent(event);
  assert.ok(event.defaultPrevented, "the page doesn't scroll under the map");
  const b = box(el);
  close(b, { w: 100 / 1.15, h: 100 / 1.15 });
  // The anchor stays at the same place on screen.
  close({ x: b.x + 0.25 * b.w, y: b.y + 0.75 * b.h }, { x: 25, y: 75 }, "anchor");
  el.dispatchEvent(new FakeEvent("wheel", { clientX: 100, clientY: 300, deltaY: 1 }));
  close(box(el), { x: 0, y: 0, w: 100, h: 100 });
});

test("dragging pans the view, and the click that ends a drag is swallowed", async () => {
  const { el, changes } = await setUp("0 0 100 100", { "data-min-view-size": "1", "data-max-view-size": "1000" });
  const marker = h("a", { href: "#x" });
  el.appendChild(marker);
  let clicks = 0;
  marker.addEventListener("click", () => { clicks += 1; });

  pointer(el, "pointerdown", 200, 200);
  pointer(el, "pointermove", 202, 201);
  assert.equal(changes.length, 0, "under 4 px is not a drag");
  pointer(el, "pointermove", 240, 180);
  // 40 px right, 20 px up on a 400 px map of 100 units: 10 left, 5 down.
  close(box(el), { x: -10, y: 5, w: 100, h: 100 });
  pointer(el, "pointerup", 240, 180);
  marker.dispatchEvent(new FakeEvent("click", { button: 0 }));
  assert.equal(clicks, 0, "the drag's click doesn't follow the marker's link");
  marker.dispatchEvent(new FakeEvent("click", { button: 0 }));
  assert.equal(clicks, 1, "only that one click is swallowed");
});

test("a plain click (no drag) reaches the marker and doesn't move the view", async () => {
  const { el, changes } = await setUp("0 0 100 100", { "data-min-view-size": "1" });
  const marker = h("a", { href: "#x" });
  el.appendChild(marker);
  let clicks = 0;
  marker.addEventListener("click", () => { clicks += 1; });
  pointer(marker, "pointerdown", 200, 200);
  pointer(marker, "pointermove", 201, 201);
  pointer(marker, "pointerup", 201, 201);
  marker.dispatchEvent(new FakeEvent("click", { button: 0 }));
  assert.equal(clicks, 1);
  assert.equal(changes.length, 0);
});

test("only the left button drags, and only the pointer that went down", async () => {
  const { el, changes } = await setUp("0 0 100 100", { "data-min-view-size": "1" });
  pointer(el, "pointerdown", 200, 200, { button: 2 });
  pointer(el, "pointermove", 300, 300);
  assert.equal(changes.length, 0);
  pointer(el, "pointerdown", 200, 200);
  pointer(el, "pointermove", 300, 300, { pointerId: 2 });
  assert.equal(changes.length, 0);
});

test("the controller's reset and getBox", async () => {
  const { el, controller, zoomInBtn } = await setUp("-1 -2 30 30", { "data-min-view-size": "1" });
  zoomInBtn.click();
  close(controller.getBox(), box(el));
  controller.reset();
  assert.deepEqual(controller.getBox(), { x: -1, y: -2, w: 30, h: 30 });
});
