// tests/js/systemreach.test.mjs -- static/systemreach.js (MAP.152): which star's system the camera wants up in an
// opened sector, and how much of it shows. The spec is docs/design/fly-through-view-distance.md 4.4 and 4.7: a system's
// disc appears at about 20 pixels, is whole at 40 and is dropped below 10.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const R = await importPage("systemreach.js");
const VIEW = { heightPx: 1000, fovDeg: 50 };
const AU_PC = 1 / 206264.806;

test("the numbers are the design's", () => {
  assert.equal(R.LOAD_PX, 20);
  assert.equal(R.UNLOAD_PX, 10);
  assert.ok(R.FULL_PX > R.LOAD_PX);
  assert.ok(Math.abs(R.NOMINAL_RADIUS_PC - 50 * AU_PC) < 1e-12);
});

test("a 50 AU system on a 1000 px map at 50 degrees is 20 px across a few thousand AU away (design 4.4: about 3,000)", () => {
  const d = R.distanceForPixels(20, R.NOMINAL_RADIUS_PC, VIEW.heightPx, VIEW.fovDeg);
  assert.ok(d / AU_PC > 2500 && d / AU_PC < 6000, `${d / AU_PC} AU`);
  assert.ok(Math.abs(R.discPixels(d, R.NOMINAL_RADIUS_PC, VIEW.heightPx, VIEW.fovDeg) - 20) < 1e-9);
});

test("the pixel size falls with the distance and rises with the radius", () => {
  const near = R.discPixels(0.01, R.NOMINAL_RADIUS_PC, 1000, 50);
  const far = R.discPixels(0.02, R.NOMINAL_RADIUS_PC, 1000, 50);
  assert.ok(Math.abs(near / far - 2) < 1e-9);
  assert.ok(R.discPixels(0.01, 2 * R.NOMINAL_RADIUS_PC, 1000, 50) > near);
  assert.equal(R.discPixels(0, 1, 1000, 50), Infinity);
});

test("the share is 0 at 20 px, 1 at 40 px and smooth between", () => {
  assert.equal(R.systemShare(0), 0);
  assert.equal(R.systemShare(20), 0);
  assert.equal(R.systemShare(40), 1);
  assert.equal(R.systemShare(500), 1);
  let last = 0;
  for (let px = 20; px <= 40; px += 1) {
    const share = R.systemShare(px);
    assert.ok(share >= last - 1e-12);
    last = share;
  }
  assert.ok(Math.abs(R.systemShare(30) - 0.5) < 1e-9);
});

const at = (px, key = "a") => {
  // A star `px`-pixels-across away from an eye at the origin.
  const distance = R.distanceForPixels(px, R.NOMINAL_RADIUS_PC, VIEW.heightPx, VIEW.fovDeg);
  return { key, center: [distance, 0, 0] };
};

test("nothing is wanted until a disc is 20 px across", () => {
  assert.equal(R.chooseSystem([0, 0, 0], [at(19)], VIEW, null), null);
  const found = R.chooseSystem([0, 0, 0], [at(21)], VIEW, null);
  assert.equal(found.key, "a");
  assert.ok(found.share > 0 && found.share < 0.05);
});

test("a system up stays between 10 and 20 px with no share, and goes below 10", () => {
  const kept = R.chooseSystem([0, 0, 0], [at(15)], VIEW, "a");
  assert.equal(kept.key, "a");
  assert.equal(kept.share, 0);
  assert.equal(R.chooseSystem([0, 0, 0], [at(9)], VIEW, "a"), null);
  assert.equal(R.chooseSystem([0, 0, 0], [at(15)], VIEW, null), null, "and it does not come in from there");
});

test("a camera hovering near the size does not flip the scene on and off", () => {
  let current = null;
  const states = [];
  for (const px of [21, 19, 21, 12, 19, 11, 21, 9, 19]) {
    const found = R.chooseSystem([0, 0, 0], [at(px)], VIEW, current);
    current = found ? found.key : null;
    states.push(current);
  }
  assert.deepEqual(states, ["a", "a", "a", "a", "a", "a", "a", null, null]);
});

test("a known radius is used in place of the nominal one", () => {
  const c = at(15);
  assert.equal(R.chooseSystem([0, 0, 0], [c], VIEW, null), null);
  const bigger = R.chooseSystem([0, 0, 0], [Object.assign({}, c, { radiusPc: 3 * R.NOMINAL_RADIUS_PC })], VIEW, null);
  assert.equal(bigger.key, "a");
  assert.equal(bigger.share, 1);
});

test("one scene is up at a time: a system takes over only when clearly larger", () => {
  const eye = [0, 0, 0];
  const a = at(30, "a");
  const b = Object.assign(at(30, "b"), { center: [0, R.distanceForPixels(30, R.NOMINAL_RADIUS_PC, 1000, 50), 0] });
  assert.equal(R.chooseSystem(eye, [a, b], VIEW, "a").key, "a");
  assert.equal(R.chooseSystem(eye, [a, b], VIEW, "b").key, "b");
  const bBig = Object.assign({}, b, { center: [0, R.distanceForPixels(40, R.NOMINAL_RADIUS_PC, 1000, 50), 0] });
  assert.equal(R.chooseSystem(eye, [a, bBig], VIEW, "a").key, "b");
});
