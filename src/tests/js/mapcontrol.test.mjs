// tests/js/mapcontrol.test.mjs -- static/mapcontrol.js, the camera and
// input controller the Galaxy Map's drill-down and the Sector Map share
// (MAP.64): zoom policies, the orbit and pan moves, and telling a drag
// from a click (both maps' own tests cover how each uses it).

import { test } from "node:test";
import assert from "node:assert/strict";

import { FakeEvent, h, installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/");
const THREE = await import(new URL("vendor/three.module.min.js", STATIC_URL).href);
const MC = await import(new URL("mapcontrol.js", STATIC_URL).href);

const press = (el, type, x, y, extra) =>
  el.dispatchEvent(new FakeEvent(type, Object.assign({ clientX: x, clientY: y, pointerId: 1, button: 0, pointerType: "mouse" }, extra)));

test("a range policy holds the distance between its multiples of the fit", () => {
  const range = MC.zoomPolicy(MC.ZOOM_RANGE, 0.5, 2);
  assert.equal(MC.clampDistance(range, 10, 100), 50);
  assert.equal(MC.clampDistance(range, 1000, 100), 200);
  assert.equal(MC.clampDistance(range, 120, 100), 120);
  assert.ok(MC.canZoom(range));
});

test("a locked policy never zooms, a free one never clamps", () => {
  const locked = MC.zoomPolicy(MC.ZOOM_LOCKED);
  assert.equal(MC.canZoom(locked), false);
  assert.equal(MC.clampDistance(locked, 10, 100), 100);
  const free = MC.zoomPolicy(MC.ZOOM_FREE);
  assert.equal(MC.clampDistance(free, 1e9, 100), 1e9);
  assert.equal(MC.canZoom(null), false);
});

test("orbitByDrag and orbitByKey turn and tilt within the clamp", () => {
  const clamp = (phi) => Math.max(0.1, Math.min(1, phi));
  const view = { theta: 0, phi: 0.5 };
  MC.orbitByDrag(view, 10, -1000, 0.01, clamp);
  assert.ok(Math.abs(view.theta + 0.1) < 1e-12);
  assert.equal(view.phi, 1);
  assert.ok(MC.orbitByKey(view, "ArrowLeft", 0.2, clamp));
  assert.ok(Math.abs(view.theta - 0.1) < 1e-12);
  assert.ok(MC.orbitByKey(view, "ArrowUp", 5, clamp));
  assert.equal(view.phi, 0.1);
  assert.equal(MC.orbitByKey(view, "a", 0.2, clamp), false);
});

test("panInScreenPlane moves the target against the drag", () => {
  const camera = new THREE.PerspectiveCamera(60, 1, 1, 1000);
  camera.position.set(0, 0, 10);
  camera.lookAt(0, 0, 0);
  camera.updateMatrixWorld();
  const target = MC.panInScreenPlane(THREE, camera, [0, 0, 0], 10, 4, 0.5);
  assert.ok(Math.abs(target[0] + 5) < 1e-9 && Math.abs(target[1] - 2) < 1e-9 && Math.abs(target[2]) < 1e-9);
});

test("wheelPixels turns lines and pages into pixels and caps them", () => {
  assert.equal(MC.wheelPixels({ deltaY: 3, deltaMode: 1 }, 400, 200), 99);
  assert.equal(MC.wheelPixels({ deltaY: 1, deltaMode: 2 }, 400, 200), 200);
  assert.equal(MC.wheelPixels({ deltaY: -50, deltaMode: 0 }, 400, 200), -50);
});

test("a press that hardly moves clicks; past the threshold it drags", () => {
  const canvas = h("canvas");
  const seen = [];
  MC.createPointerControl(canvas, {
    attach: true, dragClickPx: 6,
    onDragStart: () => seen.push("start"),
    onDrag: (dx, dy, pan) => seen.push(["drag", dx, dy, pan]),
    onClick: (event, type) => seen.push(["click", type]),
  });
  press(canvas, "pointerdown", 0, 0);
  press(canvas, "pointermove", 3, 2);
  press(canvas, "pointerup", 3, 2);
  assert.deepEqual(seen, [["click", "mouse"]]);
  seen.length = 0;
  press(canvas, "pointerdown", 0, 0);
  press(canvas, "pointermove", 3, 2);
  press(canvas, "pointermove", 6, 3);
  press(canvas, "pointerup", 6, 3);
  assert.deepEqual(seen, ["start", ["drag", 3, 1, false]]);
});

test("a pan press never clicks, and the right button never clicks", () => {
  const canvas = h("canvas");
  const seen = [];
  MC.createPointerControl(canvas, {
    attach: true, dragClickPx: 6, buttons: [0, 2],
    isPan: (event) => event.button === 2 || event.shiftKey,
    onClick: () => seen.push("click"),
  });
  press(canvas, "pointerdown", 0, 0, { shiftKey: true });
  press(canvas, "pointerup", 0, 0, { shiftKey: true });
  press(canvas, "pointerdown", 0, 0, { button: 2 });
  press(canvas, "pointerup", 0, 0, { button: 2 });
  press(canvas, "pointerdown", 0, 0, { button: 1 });
  press(canvas, "pointerup", 0, 0, { button: 1 });
  assert.deepEqual(seen, []);
});

test("turnAtOnce follows from the first move and eats the click after a drag", () => {
  const canvas = h("canvas");
  const seen = [];
  MC.createPointerControl(canvas, {
    attach: true, dragClickPx: 4, measure: "path", turnAtOnce: true, clickOn: "click",
    onDrag: (dx) => seen.push(dx),
    onClick: () => seen.push("click"),
  });
  press(canvas, "pointerdown", 0, 0);
  press(canvas, "pointermove", 2, 0);
  press(canvas, "pointermove", 0, 0);
  press(canvas, "pointerup", 0, 0);
  canvas.dispatchEvent(new FakeEvent("click", { clientX: 0, clientY: 0 }));
  assert.deepEqual(seen, [2, -2, "click"], "4 px of path is still a click");
  seen.length = 0;
  press(canvas, "pointerdown", 0, 0);
  press(canvas, "pointermove", 3, 0);
  press(canvas, "pointermove", 0, 0);
  press(canvas, "pointerup", 0, 0);
  canvas.dispatchEvent(new FakeEvent("click", { clientX: 0, clientY: 0 }));
  assert.deepEqual(seen, [3, -3], "6 px of path back to the start is a drag");
});

test("two fingers pinch: the ratio is the first spread over the spread now", () => {
  const canvas = h("canvas");
  const ratios = [];
  let started = 0;
  const control = MC.createPointerControl(canvas, {
    dragClickPx: 6,
    pinch: { canStart: () => true, start: () => started++, move: (r) => ratios.push(r) },
    onClick: () => ratios.push("click"),
  });
  const touch = (fn, id, x) => fn(new FakeEvent("p", { pointerId: id, pointerType: "touch", clientX: x, clientY: 0, button: 0 }));
  touch(control.down, 1, 0);
  touch(control.down, 2, 100);
  touch(control.move, 2, 200);
  touch(control.up, 2, 200);
  touch(control.up, 1, 0);
  assert.equal(started, 1);
  assert.deepEqual(ratios, [0.5]);
});

test("zoomAboutAnchor scales the target about the anchor so the anchor holds its place", () => {
  assert.deepEqual(MC.zoomAboutAnchor([10, 0, 0], [0, 0, 0], 0.5), [5, 0, 0]);
  assert.deepEqual(MC.zoomAboutAnchor([10, 4, 2], [10, 4, 2], 0.1), [10, 4, 2], "anchored on the target: nothing moves");
  assert.deepEqual(MC.zoomAboutAnchor([1, 2, 3], [5, 6, 7], 1), [1, 2, 3]);
});

test("pointOnFocusPlane crosses the plane through the target at right angles to the view", () => {
  const point = MC.pointOnFocusPlane([0, 0, 10], [0, 0, -1], [0, 0, 0], [0, 0, -1]);
  assert.deepEqual(point, [0, 0, 0]);
  // A slanted ray reaches the plane further out along the ray.
  const slant = Math.SQRT1_2;
  const off = MC.pointOnFocusPlane([0, 0, 10], [slant, 0, -slant], [0, 0, 0], [0, 0, -1]);
  assert.ok(Math.abs(off[0] - 10) < 1e-9 && Math.abs(off[2]) < 1e-9, String(off));
  assert.equal(MC.pointOnFocusPlane([0, 0, 10], [1, 0, 0], [0, 0, 0], [0, 0, -1]), null, "along the plane");
  assert.equal(MC.pointOnFocusPlane([0, 0, 10], [0, 0, 1], [0, 0, 0], [0, 0, -1]), null, "away from it");
});
