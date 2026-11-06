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

test("a map with a row above it that stays in view is the share of the window less the row (TEST.634)", () => {
  assert.equal(F.fillHeight(1000, { above: 40 }), 950 - 40, "the row and the map are the block that is 95% of the window");
  assert.equal(F.fillHeight(1056, { covered: 56, above: 40 }), 950 - 40);
  assert.equal(F.fillHeight(1000, { above: -5 }), 950, "a row above that measures less than nothing takes nothing");
  assert.equal(F.fillHeight(1000, { above: NaN }), 950);
  assert.equal(F.fillHeight(400, { above: 100 }), F.MIN_HEIGHT_PX, "a short window keeps the least");
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

test("a box of fixed shape is as wide as keeps it inside what the window leaves (UX.142)", () => {
  assert.equal(F.fitWidth(1000, 1), 950, "a square diagram in a 1000 px window");
  assert.equal(F.fitWidth(1000, 5 / 3), Math.round(950 * 5 / 3), "the nav chart, 5:3");
  assert.equal(F.fitWidth(1056, 1, { covered: 56 }), 950, "under a header held at the top of the window");
  assert.equal(F.fitWidth(300, 1), F.MIN_HEIGHT_PX, "a short window keeps the least height");
  assert.equal(F.fitWidth(NaN, 2), 2 * F.MIN_HEIGHT_PX);
  assert.equal(F.fitWidth(1000, 0), Infinity, "no shape, no cap");
  assert.equal(F.fitWidth(1000, NaN), Infinity);
});

test("watching a box sets its maximum width and follows the window (UX.142)", () => {
  const view = { style: { maxWidth: "" } };
  const watcher = F.watchWindowBox(view, 5 / 3);
  const first = parseInt(view.style.maxWidth, 10);
  assert.ok(first > 0, "a width is set at once");
  assert.equal(first, F.fitWidth(window.innerHeight, 5 / 3, { covered: 0 }));
  window.innerHeight = 500;
  watcher.refit();
  assert.equal(parseInt(view.style.maxWidth, 10), F.fitWidth(500, 5 / 3, { covered: 0 }));
  watcher.stop();
});

test("a filling map scrolls the page to itself once, and not after the reader has pressed a key (TEST.599)", () => {
  F.forgetReader();
  const scrolls = [];
  const original = { scrollTo: window.scrollTo, scrollY: window.scrollY };
  window.scrollY = 0;
  document.documentElement.scrollHeight = 3000;
  window.scrollTo = (x, y) => { scrolls.push(y); window.scrollY = y; };
  const view = { style: { height: "" }, getBoundingClientRect: () => ({ top: 700 - window.scrollY, height: 500, bottom: 1200 - window.scrollY }), closest: () => null };
  const watcher = F.watchWindowFit(view, null, { fill: true });
  try {
    assert.equal(scrolls.length, 1, "the page is scrolled to the map when it is watched");
    window.scrollY = 0;                                  // the page is back at the top, as a reader who has not touched it leaves it
    window.dispatchEvent(new Event("keydown"));          // a Tab to the skip link at the top of the page
    watcher.refit();
    assert.equal(scrolls.length, 1, "the page is left where the reader's key left it");
    const later = F.watchWindowFit(view, null, { fill: true });          // a map built after the key (the 3D view, when its data has come)
    assert.equal(scrolls.length, 1, "a map built after the reader's key does not scroll the page either");
    later.stop();
  } finally {
    watcher.stop();
    window.scrollTo = original.scrollTo;
    window.scrollY = original.scrollY;
  }
});

test("a pointer press and the wheel also leave the page where it is", () => {
  for (const type of ["pointerdown", "wheel", "touchstart"]) {
    F.forgetReader();
    const scrolls = [];
    window.scrollY = 0;
    document.documentElement.scrollHeight = 3000;
    window.scrollTo = (x, y) => { scrolls.push(y); window.scrollY = y; };
    const view = { style: { height: "" }, getBoundingClientRect: () => ({ top: 700 - window.scrollY, height: 500, bottom: 1200 - window.scrollY }), closest: () => null };
    const watcher = F.watchWindowFit(view, null, { fill: true });
    scrolls.length = 0;
    window.scrollY = 0;
    window.dispatchEvent(new Event(type));
    watcher.refit();
    assert.equal(scrolls.length, 0, type);
    watcher.stop();
  }
});

// A page whose map has a row above it (the Galaxy Map's Go to row, TEST.634): the row sits 40 px above the map, the map's top is 700 px down the page.
function rowAndMap(rowHeight) {
  const state = { rowHeight: rowHeight };
  const row = { getBoundingClientRect: () => ({ top: 660 - window.scrollY, height: state.rowHeight, bottom: 660 + state.rowHeight - window.scrollY }) };
  const view = {
    style: { height: "" },
    getBoundingClientRect: () => {
      const height = parseInt(view.style.height, 10) || 500;
      return { top: 700 - window.scrollY, height: height, bottom: 700 + height - window.scrollY };
    },
    closest: () => null,
  };
  return { row, view, state };
}

test("the row and the map are the block that fills the window, and the page rests with the row at the buffer (TEST.634)", () => {
  F.forgetReader();
  const original = { scrollTo: window.scrollTo, scrollY: window.scrollY, innerHeight: window.innerHeight };
  window.scrollY = 0;
  window.innerHeight = 1000;
  document.documentElement.scrollHeight = 3000;
  window.scrollTo = (x, y) => { window.scrollY = y; };
  const { row, view } = rowAndMap(40);
  const watcher = F.watchWindowFit(view, null, { fill: true, block: row });
  try {
    assert.equal(parseInt(view.style.height, 10), 950 - 40, "the map is the row's height shorter");
    const blockHeight = 40 + 910;
    const buffer = Math.round((1000 - blockHeight) / 2);
    assert.equal(window.scrollY, 660 - buffer, "the row's top, not the map's, is at the buffer");
    assert.ok(row.getBoundingClientRect().top >= buffer - 1, "the row is in view");
  } finally {
    watcher.stop();
    Object.assign(window, original);
  }
});

test("a row that is not shown leaves the map to fill the window alone (TEST.634)", () => {
  F.forgetReader();
  const original = { scrollTo: window.scrollTo, scrollY: window.scrollY, innerHeight: window.innerHeight };
  window.scrollY = 0;
  window.innerHeight = 1000;
  document.documentElement.scrollHeight = 3000;
  window.scrollTo = (x, y) => { window.scrollY = y; };
  const { row, view, state } = rowAndMap(0);
  const watcher = F.watchWindowFit(view, null, { fill: true, block: () => row });
  try {
    assert.equal(parseInt(view.style.height, 10), 950, "a hidden row has no height");
    assert.equal(window.scrollY, 700 - 25, "the map's own top is at the buffer");
    state.rowHeight = 40;                                  // the row is shown
    watcher.refit();
    assert.equal(parseInt(view.style.height, 10), 910, "the map gives the row its room");
  } finally {
    watcher.stop();
    Object.assign(window, original);
  }
});
