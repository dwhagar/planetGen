// tests/js/mapcore.test.mjs -- static/mapcore.js, the helpers the three
// maps share (MAP.63): the scale bar's nice numbers and screen span, the
// info panel's fields, reading the scene JSON, and picking on screen.

import { test } from "node:test";
import assert from "node:assert/strict";

import { h, installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/");
const THREE = await import(new URL("vendor/three.module.min.js", STATIC_URL).href);
const M = await import(new URL("mapcore.js", STATIC_URL).href);

test("niceScaleValue snaps to 1, 2 or 5 times a power of ten", () => {
  assert.equal(M.niceScaleValue(1.2), 1);
  assert.equal(M.niceScaleValue(2.9), 2);
  assert.equal(M.niceScaleValue(6.283), 5);
  assert.equal(M.niceScaleValue(8), 10);
  assert.equal(M.niceScaleValue(0.034), 0.02);
  for (const bad of [0, -3, NaN, Infinity]) assert.equal(M.niceScaleValue(bad), 0);
});

test("worldUnitsPerPixel is the view's height at that distance over the canvas height", () => {
  const camera = new THREE.PerspectiveCamera(90, 1, 1, 1e6);
  assert.ok(Math.abs(M.worldUnitsPerPixel(camera, 50, 100) - 1) < 1e-12);
  assert.ok(Math.abs(M.worldUnitsPerPixel(camera, 50, 0) - 100) < 1e-9, "a zero-height canvas counts as 1 px");
});

test("addField leaves out an empty value but keeps 0", () => {
  const dl = h("dl");
  M.addField(dl, "Empty", "");
  M.addField(dl, "Missing", null);
  M.addField(dl, "Zero", 0);
  M.addField(dl, "Name", "<b>x</b>");
  const terms = dl.querySelectorAll("dt").map((dt) => dt.textContent);
  assert.deepEqual(terms, ["Zero", "Name"]);
  assert.equal(dl.querySelectorAll("dd")[1].textContent, "<b>x</b>", "plain text, never markup");
});

test("readSceneData parses the JSON block, null when missing or broken", () => {
  assert.deepEqual(M.readSceneData({ textContent: '{"a": 1}' }), { a: 1 });
  assert.equal(M.readSceneData({ textContent: "{" }), null);
  assert.equal(M.readSceneData(null), null);
});

test("formatAddress", () => {
  assert.equal(M.formatAddress(5, -1, 20), "ring 5 layer -1 slot 20");
});

test("nearestOnScreen picks the nearest entry within its reach", () => {
  const camera = new THREE.PerspectiveCamera(90, 1, 1, 1000);
  camera.position.set(0, 0, 100);
  camera.lookAt(0, 0, 0);
  camera.updateMatrixWorld();
  const rect = { left: 0, top: 0, width: 200, height: 200 };
  const a = { x: 0, y: 0, z: 0, name: "a" };
  const b = { x: 10, y: 0, z: 0, name: "b" };
  const behind = { x: 0, y: 0, z: 200, name: "behind" };
  const reach = () => 8;
  // 100 units away with a 90 degree view: 1 unit = 1 px on a 200 px canvas.
  assert.equal(M.nearestOnScreen([a, b], camera, rect, 103, 100, { reach }).entry, a);
  assert.equal(M.nearestOnScreen([a, b], camera, rect, 108, 100, { reach }).entry, b);
  assert.equal(M.nearestOnScreen([a, b], camera, rect, 100, 130, { reach }), null, "out of reach");
  assert.equal(M.nearestOnScreen([a, b], camera, rect, 103, 100, { reach, accept: (e) => e !== a }).entry, b);
  assert.equal(M.nearestOnScreen([behind], camera, rect, 100, 100, { reach }), null, "behind the camera");
  const twin = { x: 0, y: 0, z: 0, name: "twin" };
  assert.equal(M.nearestOnScreen([a, twin], camera, rect, 100, 100, { reach }).entry, a, "a tie goes to the first");
  assert.equal(M.nearestOnScreen([a, twin], camera, rect, 100, 100, { reach, lastWins: true }).entry, twin);
  assert.equal(M.nearestOnScreen([a], camera, { left: 0, top: 0, width: 0, height: 0 }, 0, 0, { reach }), null);
});
