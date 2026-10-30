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

test("a map that fills its window is 95% of the window (MAP.350)", () => {
  assert.equal(F.fillHeight(1000), 950);
  assert.equal(F.fillHeight(950), Math.round(950 * 0.95));
  assert.equal(F.fillHeight(1080, { share: 1 }), 1080);
  assert.equal(F.fillHeight(300), F.MIN_HEIGHT_PX, "a short window keeps the least");
  assert.equal(F.fillHeight(NaN), F.MIN_HEIGHT_PX);
  assert.equal(F.fillHeight(1056, { covered: 56 }), 950, "95% of what a header at the top of the window leaves");
});

test("the margin above a filling map is half of what it leaves", () => {
  assert.equal(F.fillBuffer(1000, 950, 0), 25);
  assert.equal(F.fillBuffer(1056, 950, 56), 25);
  assert.equal(F.fillBuffer(300, 320, 0), 0, "never negative");
  assert.equal(F.fillBuffer(NaN, 320, 0), 0);
});

test("the page is scrolled so the map has the buffer above it, and left alone where it is already there (MAP.350)", () => {
  assert.equal(F.alignScroll(319, 0, 700, { buffer: 20 }), 319 - 20);
  assert.equal(F.alignScroll(319, 299, 700, { buffer: 20 }), null);
  assert.equal(F.alignScroll(319, 300, 700, { buffer: 20 }), null, "a pixel or two off is there");
  assert.equal(F.alignScroll(319, 0, 100, { buffer: 20 }), 100, "never past the end of the page");
  assert.equal(F.alignScroll(10, 0, 700, { buffer: 20 }), null, "never above the top");
  assert.equal(F.alignScroll(NaN, 0, 700), null);
  assert.equal(F.alignScroll(319, 0, 700, { covered: 56, buffer: 20 }), 319 - 56 - 20, "under a header held at the top of the window");
});
