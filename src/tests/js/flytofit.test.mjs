// tests/js/flytofit.test.mjs -- static/flytofit.js, where the camera sits for a target to fill the view (MAP.319).

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/");
const F = await import(new URL("flytofit.js", STATIC_URL).href);

const near = (a, b, eps = 1e-9) => assert.ok(Math.abs(a - b) <= eps * Math.max(1, Math.abs(b)), `${a} is not ${b}`);

test("the distance is the radius over the sine of half the smaller field of view, plus the margin", () => {
  // A square window: both fields are 50 degrees.
  near(F.fitDistance(10, 50, 1), (10 / Math.sin((25 * Math.PI) / 180)) * 1.08);
  assert.equal(F.FIT_MARGIN, 0.08);
});

test("a wide window is limited by its height, a tall one by its width", () => {
  const square = F.fitDistance(10, 50, 1);
  near(F.fitDistance(10, 50, 2), square);
  assert.ok(F.fitDistance(10, 50, 0.5) > square, "a tall window has the narrower field, so the camera is farther");
  const narrow = Math.atan(Math.tan((25 * Math.PI) / 180) * 0.5);
  near(F.fitDistance(10, 50, 0.5), (10 / Math.sin(narrow)) * 1.08);
});

test("the sphere fills the smaller side of the view and no more", () => {
  for (const [fov, aspect] of [[50, 1], [50, 1.6], [50, 0.6], [30, 2.4]]) {
    const radius = 7;
    const distance = F.fitDistance(radius, fov, aspect, { margin: 0 });
    const half = F.halfFieldRad(fov, aspect);
    // The sphere's edge, seen from the camera, is exactly half the smaller field.
    near(Math.asin(radius / distance), half);
  }
});

test("a sphere smaller than the least is framed as the least", () => {
  near(F.fitDistance(0.01, 50, 1, { minRadius: F.MIN_VIEW_RADIUS_PC }), F.fitDistance(F.MIN_VIEW_RADIUS_PC, 50, 1));
  near(F.fitDistance(5, 50, 1, { minRadius: F.MIN_VIEW_RADIUS_PC }), F.fitDistance(5, 50, 1));
  assert.equal(F.MIN_VIEW_RADIUS_PC, 2);
});

test("a radius that is not a number frames the least, never NaN", () => {
  assert.equal(F.fitDistance(NaN, 50, 1), 0);
  near(F.fitDistance(NaN, 50, 1, { minRadius: 2 }), F.fitDistance(2, 50, 1));
});

test("the reach is the farthest edge of the parts from the middle", () => {
  const parts = [{ x: 3, y: 0, z: 0, r: 1 }, { x: 0, y: -2, z: 0, r: 0 }, { x: NaN, y: 0, z: 0, r: 99 }];
  assert.equal(F.reachFrom([0, 0, 0], parts), 4);
  assert.equal(F.reachFrom([3, 0, 0], parts), Math.hypot(3, 2));
  assert.equal(F.reachFrom([0, 0, 0], []), 0);
});

test("an orbit's reach is its farthest vertex", () => {
  assert.equal(F.reachOfPositions([1, 0, 0, 0, 5, 0, 0, 0, 2]), 5);
  assert.equal(F.reachOfPositions([]), 0);
});
