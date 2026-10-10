// tests/js/nearfield.test.mjs -- static/nearfield.js, the fly-through map's
// near field (MAP.149): the depth fade, the see-through focus tube, the
// prominence of a region, and the GLSL twin the shaders use.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const N = await importPage("nearfield.js");

const D = 100;

test("things nearer than a sixth of the focus distance are gone and half of it is whole", () => {
  assert.equal(N.depthFade(0, D), 0);
  assert.equal(N.depthFade(0.15 * D, D), 0);
  assert.equal(N.depthFade(0.5 * D, D), 1);
  assert.equal(N.depthFade(3 * D, D), 1);
  const mid = N.depthFade(0.325 * D, D);
  assert.ok(Math.abs(mid - 0.5) < 1e-12, String(mid));
  let last = 0;
  for (let depth = 0; depth <= D; depth += 1) {
    const share = N.depthFade(depth, D);
    assert.ok(share >= last, "never brightens toward the camera");
    last = share;
  }
});

test("the tube thins what stands between the camera and the focus and nothing else", () => {
  // Inside the tube, well in front of the focus: the minimum.
  assert.ok(Math.abs(N.tubeFade(0, 0, 0.4 * D, D, 0) - N.TUBE_MIN_ALPHA) < 1e-12);
  // Beside it, past 1.6 radii: untouched.
  const rt = N.tubeRadius(D, 0);
  assert.equal(rt, 0.12 * D);
  assert.equal(N.tubeFade(1.6 * rt + 0.01, 0, 0.4 * D, D, 0), 1);
  // At and behind the focus: untouched, however central.
  assert.equal(N.tubeFade(0, 0, D, D, 0), 1);
  assert.equal(N.tubeFade(0, 0, 2 * D, D, 0), 1);
  // Continuous across the focus plane: a point just in front is nearly 1.
  assert.ok(N.tubeFade(0, 0, D - 0.01 * rt, D, 0) > 0.99);
  // A wide focus makes a wide tube.
  assert.equal(N.tubeRadius(D, 30), 30);
  assert.ok(Math.abs(N.tubeFade(25, 0, 0.4 * D, D, 30) - N.TUBE_MIN_ALPHA) < 1e-12);
});

test("the near field is the depth fade times the tube, worked in view space", () => {
  // Camera at the origin looking down -z; the focus is at z = -D.
  assert.equal(N.nearField(0, 0, -0.1 * D, D, 0), 0);
  assert.equal(N.nearField(0, 0, -D, D, 0), 1);
  assert.ok(Math.abs(N.nearField(0, 0, -0.6 * D, D, 0) - N.TUBE_MIN_ALPHA) < 1e-12);
  assert.equal(N.nearField(0.5 * D, 0, -0.6 * D, D, 0), 1);
  // Unknown distance: nothing is faded.
  assert.equal(N.nearField(0, 0, -1, 0, 0), 1);
});

test("a world point goes to view space through the camera's inverse matrix", () => {
  // A camera at (10, 0, 0) looking down -x with +z up: its inverse maps
  // world (0, 0, 0) to z = -10 in front of it.
  const inverse = [
    0, 0, -1, 0, // column 0: world x -> view z (negated)
    0, 1, 0, 0, // column 1: world y -> view y
    1, 0, 0, 0, // column 2: world z -> view x
    0, 0, 10, 1, // translation
  ];
  const v = N.viewSpace(inverse, 0, 0, 0);
  assert.deepEqual(v, [0, 0, 10]);
  const w = N.viewSpace(inverse, 20, 0, 0);
  assert.deepEqual(w, [0, 0, -10]);
  assert.equal(N.nearFieldAtWorld(inverse, 20, 0, 0, 10, 0), N.nearField(0, 0, -10, 10, 0));
});

test("a faint region keeps between the floor and the ceiling and thins with distance", () => {
  assert.equal(N.contextOpacity(0, 10), N.CONTEXT_MAX);
  assert.equal(N.contextOpacity(1000, 10), N.CONTEXT_MIN);
  assert.ok(N.contextOpacity(20, 10) < N.contextOpacity(15, 10));
  assert.equal(N.prominence("focus", 5, 10), 1);
  assert.equal(N.prominence("container", 5, 10), N.CONTAINER_PROMINENCE);
  assert.equal(N.prominence("context", 5, 10), N.contextOpacity(5, 10));
  // Fainter regions show shallower stars: 2.5 magnitudes at prominence 0.
  assert.equal(N.limitingMagnitude(20, 1), 20);
  assert.equal(N.limitingMagnitude(20, 0), 17.5);
});

test("picking needs a share of 0.35", () => {
  assert.equal(N.visibleEnough(0.35), true);
  assert.equal(N.visibleEnough(0.34), false);
  assert.equal(N.visibleEnough(N.TUBE_MIN_ALPHA), false);
});

test("a camera is inside the block whose polar bounds hold it", () => {
  const b = { r0: 100, r1: 104, t0: 0, t1: Math.PI / 2, z0: -2, z1: 2 };
  assert.equal(N.blockContains(b, 102, 1, 0), true);
  assert.equal(N.blockContains(b, 102, 1, 2), false, "the top wall is outside");
  assert.equal(N.blockContains(b, 99, 0, 0), false);
  assert.equal(N.blockContains(b, 104, 0, 0), false);
  assert.equal(N.blockContains(b, 102, -1, 0), false, "the other side of the bearing");
  assert.equal(N.blockContains({ ...b, t0: 3 * Math.PI / 2, t1: 2 * Math.PI }, 102, -1, 0), true, "bearings wrap to 0..2 pi");
  const c = N.blockCenter(b);
  assert.ok(Math.abs(c[0] - 102 * Math.cos(Math.PI / 4)) < 1e-9);
  assert.ok(Math.abs(c[1] - 102 * Math.sin(Math.PI / 4)) < 1e-9);
  assert.equal(c[2], 0);
});

test("the shader's twin carries the same constants", () => {
  const glsl = N.NEAR_FIELD_GLSL;
  assert.match(glsl, /float nearField\(vec3 view, float D, float radius\)/);
  [N.DEPTH_FADE_START, N.DEPTH_FADE_END, N.TUBE_MIN_ALPHA, N.TUBE_RADIUS_FLOOR, N.TUBE_EDGE].forEach((value) => {
    assert.ok(glsl.includes(value.toFixed(4)), `constant ${value} in the GLSL`);
  });
});
