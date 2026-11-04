// tests/js/ringpaint.test.mjs -- static/ringpaint.js (MAP.275): the strip of a planet's ring from the ring entry of its look record: the same ring and seed
// give the same pixels, a gap is clearer than the ring beside it, the edges fade, a dusty ring is a smooth haze and a narrow one has a few ringlets.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const R = await importPage("ringpaint.js");

const BROAD = { kind: "broad", in: 1.3, out: 2.3, tau: 0.5, ice: 0.99, alpha: 0.4, rgb: [225, 214, 140], gaps: [[0.6, 0.08, 0.9]] };
const alphaAt = (strip, share) => strip[4 * Math.floor(share * (strip.length / 4)) + 3];

test("the same ring and seed give the same strip, another seed another one", () => {
  assert.equal(R.ringFingerprint(BROAD, 7, 256), R.ringFingerprint(BROAD, 7, 256));
  assert.notEqual(R.ringFingerprint(BROAD, 7, 256), R.ringFingerprint(BROAD, 8, 256));
  assert.equal(R.paintRing(BROAD, 7, 256).length, 256 * 4);
});

test("a resonance gap is clearer than the ring on either side of it", () => {
  const strip = R.paintRing(BROAD, 3, 512);
  const inside = alphaAt(strip, 0.6);
  assert.ok(inside < 0.25 * alphaAt(strip, 0.45) && inside < 0.25 * alphaAt(strip, 0.75), [inside, alphaAt(strip, 0.45), alphaAt(strip, 0.75)]);
  const without = R.paintRing({ ...BROAD, gaps: [] }, 3, 512);
  assert.ok(alphaAt(without, 0.6) > 3 * inside);
});

test("the edges fade and the colour is the ring's, brightened or dimmed a little", () => {
  const strip = R.paintRing({ ...BROAD, gaps: [] }, 5, 512);
  assert.ok(strip[3] < 0.3 * 255 * BROAD.alpha && strip[strip.length - 1] < 0.3 * 255 * BROAD.alpha);
  for (let p = 40; p < 470; p += 40) {
    for (let c = 0; c < 3; c++) assert.ok(Math.abs(strip[4 * p + c] - BROAD.rgb[c]) <= 0.2 * BROAD.rgb[c] + 1);
  }
});

test("a dusty ring is a haze without ringlets and its opacity never passes the peak", () => {
  const dusty = { kind: "dusty", in: 1.2, out: 2.2, tau: 1e-8, ice: 0, alpha: 0.06, rgb: [40, 30, 25], gaps: [] };
  const strip = R.paintRing(dusty, 11, 256);
  let max = 0;
  for (let p = 0; p < 256; p++) max = Math.max(max, strip[4 * p + 3]);
  assert.ok(max <= Math.round(0.06 * 255) && max >= 0.6 * 0.06 * 255, max);
  assert.ok(alphaAt(strip, 0.5) > alphaAt(strip, 0.1));
});

test("a narrow ring is a few ringlets: most of its width is far fainter than its cores", () => {
  const narrow = { kind: "narrow", in: 1.6, out: 2.0, tau: 0.01, ice: 0.1, alpha: 0.06, rgb: [60, 55, 50], gaps: [] };
  const strip = R.paintRing(narrow, 21, 512);
  const alphas = [];
  for (let p = 0; p < 512; p++) alphas.push(strip[4 * p + 3]);
  const peak = Math.max(...alphas);
  assert.ok(peak > 4 * 0.06 * 255 * 0.9, peak);
  assert.ok(alphas.filter((a) => a < 0.3 * peak).length > 0.4 * 512);
});

test("the mesh's texture coordinate is the share of the way from the inner edge to the outer one", () => {
  const attributes = { position: { count: 3, xy: [[1.5, 0], [0, 2], [-2.5, 0]], getX(i) { return this.xy[i][0]; }, getY(i) { return this.xy[i][1]; } },
    uv: { out: [], setXY(i, x, y) { this.out[i] = [x, y]; }, needsUpdate: false } };
  const THREE = { RingGeometry: function () { this.attributes = attributes; } };
  R.ringGeometry(THREE, 1.5, 2.5);
  assert.deepEqual(attributes.uv.out.map((c) => c[0]), [0, 0.5, 1]);
});
