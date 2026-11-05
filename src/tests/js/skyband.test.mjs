// tests/js/skyband.test.mjs -- static/skyband.js (VIEW.34): the glow image contract (units, the equirectangular layout), the display law and the shader text.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const B = await importPage("skyband.js");

const near = (a, b, eps = 1e-9) => Math.abs(a - b) <= eps * Math.max(1, Math.abs(a), Math.abs(b));

test("a value of one is 22 mag/arcsec^2 and each 2.5 magnitudes fainter is a tenth", () => {
  assert.equal(B.valueOfSurfaceBrightness(22), 1);
  assert.ok(near(B.valueOfSurfaceBrightness(24.5), 0.1));
  assert.ok(near(B.surfaceBrightnessOfValue(B.valueOfSurfaceBrightness(21.6)), 21.6));
  assert.equal(B.surfaceBrightnessOfValue(0), Infinity);
});

test("the display level is nothing for no light, rises with it and never passes the peak", () => {
  assert.equal(B.displayLevel(0), 0);
  assert.ok(B.displayLevel(0.4) < B.displayLevel(1.5));
  assert.ok(B.displayLevel(1e6) <= B.PEAK);
  assert.ok(near(B.displayLevel(1e6), B.PEAK, 1e-6));
  // The band at the Sun (21.6 to 23.0) is a visible grey but well under half white.
  const level = B.displayLevel(B.valueOfSurfaceBrightness(21.6));
  assert.ok(level > 0.1 && level < 0.3);
});

test("the map is laid out from the north pole down and the left edge on the anti-center", () => {
  const [u, v] = B.equirectUV([1, 0, 0]);
  assert.ok(near(u, 0.5) && near(v, 0.5));          // toward the center: the middle of the map
  assert.ok(near(B.equirectUV([-1, 0, 0])[0], 0, 1e-9) || near(B.equirectUV([-1, 0, 0])[0], 1, 1e-9));
  assert.ok(near(B.equirectUV([0, 1, 0])[0], 0.75));  // longitude +90 is three quarters of the way
  assert.ok(near(B.equirectUV([0, 0, 1])[1], 0));     // the north pole is the top row
  assert.ok(near(B.equirectUV([0, 0, -1])[1], 1));
});

test("a texel's direction maps back to its own texel center", () => {
  const w = 16, h = 8;
  for (const [c, r] of [[0, 0], [3, 2], [8, 4], [15, 7]]) {
    const [u, v] = B.equirectUV(B.texelDirection(c, r, w, h));
    assert.ok(near(u * w - 0.5, c, 1e-9) && near(v * h - 0.5, r, 1e-9), [c, r, u * w - 0.5, v * h - 0.5].join());
  }
});

test("sampling a constant image gives the constant everywhere; between two texels it is their mean", () => {
  const flat = B.makeImage(16, 8, [0.5, 1, 2]);
  for (const d of [[1, 0, 0], [0, 1, 0], [0, 0, 1], [-0.3, 0.4, -0.9], [-1, 0, 0]]) {
    assert.ok(near(B.sample(flat, d, 0), 0.5) && near(B.sample(flat, d, 1), 1) && near(B.sample(flat, d, 2), 2));
  }
  const ramp = B.makeImage(16, 8, (d) => [d[2] > 0 ? 2 : 0, 0, 0]);
  const between = B.sample(ramp, [1, 0, 0], 0);       // the equator is the line between the two middle rows
  assert.ok(near(between, 1, 1e-9));
});

test("the longitude wraps across the anti-center and the poles clamp", () => {
  const w = 16, h = 8;
  const image = B.makeImage(w, h, (d) => [Math.abs(d[1]) < 0.2 && d[0] < 0 ? 1 : 0, 0, 0]);
  const across = [-1, 0.0001, 0];
  assert.ok(B.sample(image, across, 0) > 0.4);
  assert.ok(B.sample(image, [-1, -0.0001, 0], 0) > 0.4);
  const pole = B.makeImage(w, h, [3, 3, 3]);
  assert.ok(near(B.sample(pole, [0, 0, 1], 0), 3));
});

test("the shader carries the display law and the grain and dither", () => {
  assert.match(B.FRAGMENT_SHADER, /PEAK \* \(1\.0 - exp\(-GAIN \* value\)\)/);
  assert.match(B.FRAGMENT_SHADER, /grainStrength/);
  assert.match(B.FRAGMENT_SHADER, /\/ 255\.0/);
  assert.match(B.FRAGMENT_SHADER, new RegExp("PEAK = " + B.PEAK.toFixed(6).replace(".", "\\.")));
});
