// tests/js/thermalglow.test.mjs -- static/thermalglow.js (MAP.275 part 3): the glow a hot body adds to its own picture.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const T = await importPage("thermalglow.js");

test("a look without a thermal entry, or with no strength, does not glow", () => {
  assert.equal(T.thermalOf({}), null);
  assert.equal(T.thermalOf(null), null);
  assert.equal(T.thermalOf({ thermal: { hex: "#ff8800", strength: 0 } }), null);
});

test("a thermal entry gives its colour and strength", () => {
  const got = T.thermalOf({ thermal: { hex: "#ff8800", strength: 0.5 } });
  assert.deepEqual(got.rgb, [255, 136, 0]);
  assert.equal(got.strength, 0.5);
});

test("the glow adds the colour scaled by the pixel's brightness and the strength, up to 255", () => {
  assert.equal(T.glowed(10, 200, 0, 1), 10);
  assert.ok(T.glowed(10, 200, 255, 0.5) > T.glowed(10, 200, 100, 0.5));
  assert.ok(T.glowed(10, 200, 255, 1) > T.glowed(10, 200, 255, 0.5));
  assert.equal(T.glowed(250, 255, 255, 1), 255);
});
