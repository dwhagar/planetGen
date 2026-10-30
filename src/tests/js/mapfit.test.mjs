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

test("a map that fills its window is the window less the buffer above and below (MAP.350)", () => {
  assert.equal(F.fillHeight(950), 950 - 2 * F.BUFFER_PX);
  assert.equal(F.fillHeight(1080, { buffer: 10 }), 1060);
  assert.equal(F.fillHeight(300), F.MIN_HEIGHT_PX, "a short window keeps the least");
  assert.equal(F.fillHeight(NaN), F.MIN_HEIGHT_PX);
  assert.equal(F.fillHeight(950, { covered: 56 }), 950 - 56 - 2 * F.BUFFER_PX, "less what a header at the top of the window covers");
});

test("the page is scrolled so the map has the buffer above it, and left alone where it is already there (MAP.350)", () => {
  assert.equal(F.alignScroll(319, 0, 700), 319 - F.BUFFER_PX);
  assert.equal(F.alignScroll(319, 299, 700), null);
  assert.equal(F.alignScroll(319, 300, 700), null, "a pixel or two off is there");
  assert.equal(F.alignScroll(319, 0, 100), 100, "never past the end of the page");
  assert.equal(F.alignScroll(10, 0, 700), null, "never above the top");
  assert.equal(F.alignScroll(NaN, 0, 700), null);
  assert.equal(F.alignScroll(319, 0, 700, { covered: 56 }), 319 - 56 - F.BUFFER_PX, "under a header held at the top of the window");
});
