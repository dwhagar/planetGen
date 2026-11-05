// tests/js/skysphere.test.mjs -- static/skysphere.js (VIEW.33): the light-time arithmetic that the sky's vertex shader repeats, the cap that
// keeps the brightest stars, the colours and the shader text.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const S = await importPage("skysphere.js");

const near = (a, b, eps = 1e-9) => Math.abs(a - b) <= eps * Math.max(1, Math.abs(a), Math.abs(b));

function skyOf(stars) {
  const count = stars.length;
  const flat = (key) => Float32Array.from(stars.flatMap((s) => s[key]));
  return {
    count, epochUnix: 1.8e9, limit: 7.5, velocityScale: 1, observerVelocityKms: [0, 0, 0], meta: {},
    offsetPc: flat("offset"), velocityKms: flat("velocity"), magnitude: Float32Array.from(stars.map((s) => s.M)),
    temperatureK: Float32Array.from(stars.map((s) => s.T)), kind: Uint8Array.from(stars.map(() => 0)),
    index: Float64Array.from(stars.map((_, i) => i)),
  };
}

test("the unit constants: a km/s in parsecs a year and the speed of light in parsecs a year", () => {
  assert.ok(near(S.KMS_TO_PC_PER_YEAR, 1.0227e-6, 1e-4));
  assert.ok(near(S.LIGHT_PC_PER_YEAR, 0.306601, 1e-5));
  assert.ok(near(S.LIGHT_PC_PER_YEAR / S.KMS_TO_PC_PER_YEAR, S.LIGHT_SPEED_KMS, 1e-12));
});

test("a star moving faster than 1,000 km/s takes a third light-time pass", () => {
  assert.equal(S.lightPasses(30), 2);
  assert.equal(S.lightPasses(1000), 2);
  assert.equal(S.lightPasses(1000.5), 3);
});

test("a star at rest is seen where it is, a year after the list's epoch or not", () => {
  const seen = S.retardedOffset([3, 4, 0], [0, 0, 0], [0, 0, 0], 12.5);
  assert.deepEqual(seen.apparent, [3, 4, 0]);
  assert.ok(near(seen.distance, 5));
  assert.ok(near(seen.tau, 5 / S.LIGHT_PC_PER_YEAR));
});

test("a star receding at a tenth of the speed of light is seen at (d0 + v t) / (1 + beta), the exact answer of the light-time equation", () => {
  const beta = 0.1;
  const v = beta * S.LIGHT_SPEED_KMS;
  const years = 4;
  const seen = S.retardedOffset([2, 0, 0], [v, 0, 0], [0, 0, 0], years);
  const exact = (2 + v * S.KMS_TO_PC_PER_YEAR * years) / (1 + beta);
  // Two passes leave an error of beta^2 of the light-time shift, three of beta^3: a star this fast takes three.
  assert.ok(Math.abs(seen.apparent[0] - exact) / exact < 5e-4, `${seen.apparent[0]} vs ${exact}`);
  assert.ok(near(seen.apparent[1], 0) && near(seen.apparent[2], 0));
});

test("a star crossing the line of sight is seen behind where it is by v tau, a shift of about v / c radians", () => {
  const beta = 0.01;
  const v = beta * S.LIGHT_SPEED_KMS;
  const seen = S.retardedOffset([10, 0, 0], [0, v, 0], [0, 0, 0], 0);
  const angle = Math.atan2(seen.apparent[1], seen.apparent[0]);
  assert.ok(Math.abs(angle - -beta) < 1e-4, `${angle}`);
});

test("the observer's own velocity enters only through the star's absolute velocity: a star at rest in the galaxy is seen back along the observer's path", () => {
  // The observer moves +y at 1% of c, the star is at rest in the galaxy: its velocity relative to the observer is -y.
  const v = 0.01 * S.LIGHT_SPEED_KMS;
  const seen = S.retardedOffset([10, 0, 0], [0, -v, 0], [0, v, 0], 0);
  assert.ok(near(seen.apparent[1], 0, 1e-12), "a star at rest in the galaxy has no light-time shift of its own");
  const moving = S.retardedOffset([10, 0, 0], [0, -v, 0], [0, 0, 0], 0);
  assert.ok(moving.apparent[1] > 0, "the same relative motion with a still observer is a star moving, and shifts");
});

test("the sky's years are the clock's years moved by the gap between the scene's epoch and the list's", () => {
  const sky = { epochUnix: 1000 };
  assert.equal(S.skyYears(sky, 1000, 3), 3);
  assert.ok(near(S.skyYears(sky, 1000 + S.YEAR_SECONDS, 3), 4));
  assert.ok(near(S.skyYears(sky, 1000 - S.YEAR_SECONDS, 3), 2));
});

test("a cap keeps the brightest stars at the epoch, nearest first among equals, and returns the sky itself when it fits", () => {
  const stars = [
    { offset: [100, 0, 0], velocity: [0, 0, 0], M: 5, T: 5800 },     // m 15
    { offset: [10, 0, 0], velocity: [1, 2, 3], M: 0, T: 9000 },      // m 0
    { offset: [10, 0, 0], velocity: [0, 0, 0], M: 2, T: 4000 },      // m 2
    { offset: [1, 0, 0], velocity: [0, 0, 0], M: 9, T: NaN },        // m 4
  ];
  const sky = skyOf(stars);
  assert.equal(S.keepBrightest(sky, 4), sky);
  const kept = S.keepBrightest(sky, 3);
  assert.equal(kept.count, 3);
  assert.deepEqual(Array.from(kept.index), [1, 2, 3]);
  assert.deepEqual(Array.from(kept.velocityKms.slice(0, 3)), [1, 2, 3]);
  assert.deepEqual(Array.from(kept.magnitude), [0, 2, 9]);
  assert.equal(kept.offsetPc.length, 9);
  assert.ok(Number.isNaN(kept.temperatureK[2]));
});

test("every star gets a colour: its blackbody, white for none", () => {
  const sky = skyOf([
    { offset: [1, 0, 0], velocity: [0, 0, 0], M: 0, T: 3000 },
    { offset: [1, 0, 0], velocity: [0, 0, 0], M: 0, T: 20000 },
    { offset: [1, 0, 0], velocity: [0, 0, 0], M: 0, T: NaN },
  ]);
  const colours = S.starColours(sky);
  assert.equal(colours.length, 9);
  assert.ok(colours[0] > colours[2], "a cool star is red");
  assert.ok(colours[5] > colours[3], "a hot star is blue");
  assert.deepEqual(Array.from(colours.slice(6)), [1, 1, 1]);
});

test("the vertex shader runs the passes, takes its brightness from the sky's own law and puts the star on a sphere round the camera", () => {
  const text = S.VERTEX_SHADER;
  assert.match(text, /tau = length\(apparent\) \/ LIGHT/);
  assert.match(text, /skyOpacity\(m, limit\)/);
  assert.match(text, /skySizePx\(m\) \* pixelRatio/);
  assert.match(text, /direction \* skyRadius/);
  assert.match(text, /float skyPeak\(float m\)/, "the sky's GLSL twins are in the shader");
  assert.match(text, /if \(length\(relativeKms\) > 1000\.0\)/);
  assert.match(S.FRAGMENT_SHADER, /gl_PointCoord/);
});
