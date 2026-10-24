// tests/js/mapfit.test.mjs -- static/mapfit.js, the size of a map's window on its page (MAP.291).

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/");
const F = await import(new URL("mapfit.js", STATIC_URL).href);

test("the window is what the browser window has left under the top and above the toolbar", () => {
  assert.equal(F.fitHeight(900, 300, 52, { gap: 16 }), 532);
  assert.equal(F.fitHeight(1000, 300, 52), 1000 - 300 - 52 - F.GAP_PX);
});

test("the window is never shorter than the least", () => {
  assert.equal(F.fitHeight(600, 311, 52), F.MIN_HEIGHT_PX);
  assert.equal(F.fitHeight(600, 311, 52, { min: 100 }), 221);
  assert.equal(F.fitHeight(100, 900, 52), F.MIN_HEIGHT_PX);
});

test("a window that is not a number falls back to the least", () => {
  assert.equal(F.fitHeight(NaN, 10, 10), F.MIN_HEIGHT_PX);
  assert.equal(F.fitHeight(undefined, 10, 10), F.MIN_HEIGHT_PX);
});
